!^CFG COPYRIGHT UM
!========================================================================
module ModUser
  ! This is the default user module which contains empty methods defined
  ! in ModUserEmpty.f90

  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_update_states,              &
       IMPLEMENTED2 => user_calc_sources,               &
       IMPLEMENTED4 => user_read_inputs,                &
       IMPLEMENTED5 => user_set_plot_var,               &
       IMPLEMENTED6 => user_init_session,               &
       IMPLEMENTED7 => user_set_ics

  use ModMain, ONLY: iTest, jTest, kTest, BlkTest, ProcTest, VarTest, &
       UseUserInitSession, UseUserIcs, UseUserSource, UseUserUpdateStates

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.2
  character (len=*), parameter :: &
       NameUserModule = 'HYDRO + IONIZATION EQUILIBRIUM + LEVEL SETS'

  ! Wall parameters
  real :: rInnerTube =  287.5    ! inner radius [micron]
  real :: rOuterTube =  312.5    ! outer radius [micron]
  real :: RhoDimTube = 1430.0    ! density      [kg/m3]
  real :: RhoDimOutside = 6.5    ! density  of Xe outside tube [kg/m3]
  real :: pDimOutside   = 1.1e5  ! pressure of Xe outside tube [Pa]

  ! Description of gold washer around the tube
  logical :: UseGold = .false.
  real :: WidthGold  = 50.0      ! width   [micron]
  real :: RhoDimGold = 20000.0   ! density [kg/m3]

  ! True if the plastic tube extends all the way to X=0, and Be is inside it
  logical :: IsFullTube = .false.

  ! Allow for 2D with cylindrical symmetry around the X axis
  logical :: IsCylindrical = .false.

  ! Treat cells near material interface as a mixture
  logical :: UseMixedCell = .false.
  
  ! Mixed material cell is assumed if the ratio of dominant to total
  ! atomic concentration is below MixLimit
  real :: MixLimit = 0.97

  ! Variables for Hyades file
  logical           :: UseHyadesFile = .false. ! read Hyades file?
  character(len=100):: NameHyadesFile          ! name of hyades file
  integer           :: nCellHyades             ! number of cells
  real              :: xBeHyades = -1.0        ! position of Be-Xe interface
  real, allocatable :: xHyades_C(:)            ! cell center coordinate
  real, allocatable :: StateHyades_VC(:,:)     ! cell centered Hyades state
  integer           :: iRhoHyades = 2          ! index of density
  integer           :: iUxHyades = 1           ! index of velocity
  integer           :: iPHyades = 3            ! index of pressure
  integer           :: iZHyades = 6            ! index of ionization level

contains

  !============================================================================
  subroutine user_read_inputs

    use ModReadParam
    character (len=100) :: NameCommand
    !------------------------------------------------------------------------

    UseUserUpdateStates = .true. ! for internal energy and cylindrical symm.
    UseUserInitSession  = .true. ! to set units for level set variables
    UseUserIcs          = .true. ! to read in Hyades file
    !                              and initialize the level set variables

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#HYADES")
          call read_var('UseHyadesFile', UseHyadesFile)
          call read_var('NameHyadesFile',NameHyadesFile)
       case("#TUBE")
          call read_var('IsFullTube', IsFullTube)
          call read_var('rInnerTube', rInnerTube)
          call read_var('rOuterTube', rOuterTube)
          call read_var('RhoDimTube', RhoDimTube)
          call read_var('RhoDimOutside', RhoDimOutside)
          call read_var('pDimOutside',   pDimOutside)
       case("#GOLD")
          call read_var('UseGold',    UseGold)
          call read_var('WidthGold',  WidthGold)
          call read_var('RhoDimGold', RhoDimGold)
       case("#MIXEDCELL")
          call read_var('UseMixedCell', UseMixedCell)
          if(UseMixedCell)call read_var('MixLimit', MixLimit)
       case("#CYLINDRICAL")
          call read_var('IsCylindrical', IsCylindrical)
       case('#USERINPUTEND')
          EXIT
       case default
          call stop_mpi('ERROR in ModUserCrash: unknown command='//NameCommand)
       end select
    end do

  end subroutine user_read_inputs

  !============================================================================
  subroutine user_set_ics

    use ModProcMH,  ONLY: iProc
    use ModMain,    ONLY: GlobalBlk, nI, nJ, nK, nBlock, UnusedBlk
    use ModPhysics, ONLY: ShockPosition, ShockSlope, &
         Io2No_V, No2Si_V, Si2No_V, UnitRho_, UnitP_, UnitEnergyDens_, &
         inv_gm1
    use ModAdvance, ONLY: State_VGB, Rho_, RhoUx_, RhoUy_, RhoUz_, p_, &
         ExtraEint_, LevelBe_, LevelXe_, LevelPl_
    use ModGeometry,ONLY: x_BLK, y_BLK, z_BLK, y2
    use ModEnergy,  ONLY: calc_energy_ghost
    use ModEos,     ONLY: pressure_to_eint, Be_, Xe_
    use ModPolyimide, ONLY: cAtomicMass_I, cAPolyimide

    real    :: x, y, xBe, DxBe, xTube, DyPl, pSi, RhoSi, EinternalSi
    real    :: DxyGold = -1.0
    logical :: IsError

    integer :: iBlock, i, j, k, iCell, iMaterial
    real    :: Weight1, Weight2

    character(len=*), parameter :: NameSub = "user_set_ics"
    !------------------------------------------------------------------------

    iBlock = GlobalBlk

    if(UseHyadesFile)then
       ! Read in and interpolate Hyades output
       if(.not.allocated(StateHyades_VC)) call read_hyades_file

       do i=-1,nI+2
          ! Find the Hyades points around this position
          x = x_Blk(i,1,1,iBlock)

          do iCell=1, nCellHyades
             if(xHyades_C(iCell) >= x) EXIT
          end do
          if (iCell == 1) call stop_mpi(NameSub // &
               " Hyades solution does not cover the left boundary")

          if(iCell > nCellHyades)then
             ! Cell is beyond the last point of Hyades output: use last cell
             iCell   = nCellHyades
             Weight1 = 0.0
             Weight2 = 1.0
          else
             ! Assign weights for linear interpolation between iCell-1, iCell
             Weight1 = (xHyades_C(iCell) - x) &
                  /    (xHyades_C(iCell) - xHyades_C(iCell-1))
             Weight2 = 1.0 - Weight1
          end if

          do k=-1,nk+2; do j=-1,nJ+2
             ! Interpolate density, momentum and pressure

             State_VGB(Rho_,i,j,k,iBlock) = &
                  ( Weight1*StateHyades_VC(iRhoHyades, iCell-1) &
                  + Weight2*StateHyades_VC(iRhoHyades, iCell) )

             State_VGB(RhoUx_,i,j,k,iBlock) =  State_VGB(Rho_,i,j,k,iBlock) * &
                  ( Weight1*StateHyades_VC(iUxHyades, iCell-1) &
                  + Weight2*StateHyades_VC(iUxHyades, iCell) )

             State_VGB(p_,i,j,k,iBlock) = &
                  ( Weight1*StateHyades_VC(iPHyades, iCell-1) &
                  + Weight2*StateHyades_VC(iPHyades, iCell) )

             ! Set transverse momentum to zero
             State_VGB(RhoUy_:RhoUz_,i,j,k,iBlock) = 0.0

          end do; end do
       end do
    end if

    ! Set level set functions, internal energy, and other values
    do k=-1, nK+2; do j=-1, nJ+2; do i=-1, nI+2 

       x = x_BLK(i,j,k,iBlock)
       y = y_BLK(i,j,k,iBlock)

       if(UseHyadesFile)then
          ! Be - Xe interface is given by Hyades file
          xBe = xBeHyades
       else
          ! Be - Xe interface is at the shock defined by #SHOCKPOSITION
          xBe = ShockPosition - ShockSlope*y
       end if

       ! Set the x coordinate for the left edge of the plastic tube
       if(IsFullTube)then
          xTube = 0.0
       else
          xTube = xBe
       end if

       ! Distance from Be disk: positive for x < xBe
       DxBe = xBe - x

       ! Distance from plastic wall: 
       !     positive for rInnerTube < |y| < rOuterTube only
       DyPl = min(abs(y) - rInnerTube, rOuterTube - abs(y))

       ! Distance from gold washer: positive for xTube < x < xTube + WidthGold

       if(UseGold) DxyGold = &
            min( x - xTube, xTube + WidthGold - x, abs(y) - rOuterTube )

       ! Set plastic tube state
       if(x >= xTube .and. DyPl > 0.0)then
          ! Use the density and pressure given by the #TUBE command
          State_VGB(Rho_,i,j,k,iBlock) = RhoDimTube*Io2No_V(UnitRho_)
          State_VGB(p_  ,i,j,k,iBlock) = pDimOutside*Io2No_V(UnitP_)
          ! Assume that plastic wall is at rest
          State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
       end if

       ! Set pressure and speed outside the tube
       if(abs(y) > rOuterTube) then
          State_VGB(p_,i,j,k,iBlock) = pDimOutside*Io2No_V(UnitP_)
          State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
       end if

       ! Set density outside the tube
       if(x >= xTube .and. abs(y) > rOuterTube) &
            State_VGB(Rho_,i,j,k,iBlock) = RhoDimOutside*Io2No_V(UnitRho_)
       
       ! Set density of gold washer (if present)
       if(DxyGold > 0.0) &
            State_VGB(Rho_,i,j,k,iBlock) = RhoDimGold*Io2No_V(UnitRho_)

       if(IsFullTube)then
          ! Plastic is between rInnerTube and rOuterTube
          State_VGB(LevelPl_,i,j,k,iBlock) = DyPl
          ! Berylium is inside plastic tube and left of Be/Xe interface
          State_VGB(LevelBe_,i,j,k,iBlock) = min(DxBe, rInnerTube - abs(y))
          ! Xenon is right of xBe, inside rInnerTube and outside rOuterTube
          State_VGB(LevelXe_,i,j,k,iBlock) = &
               max(abs(y) - rOuterTube, min( -DxBe, rInnerTube - abs(y)))
       else
          ! Berylium is everywhere left to the end of the plastic tube
          State_VGB(LevelBe_,i,j,k,iBlock) = DxBe
          ! Plastic is right of xBe and between rInnerTube and rOuterTube
          State_VGB(LevelPl_,i,j,k,iBlock) = min( -DxBe, DyPl)
          ! Xenon is right of xBe, inside rInnerTube and outside rOuterTube
          State_VGB(LevelXe_,i,j,k,iBlock) = min( -DxBe, -DyPl)
       end if

       ! Use plastic to represent gold
       !if(UseGold) then
       !   State_VGB(LevelPl_,i,j,k,iBlock) = &
       !        max(State_VGB(LevelPl_,i,j,k,iBlock), DxyGold)
       !   State_VGB(LevelXe_,i,k,k,iBlock) = &
       !        min(State_VGB(LevelXe_,i,j,k,iBlock), -DxyGold)
       !end if

       if(UseMixedCell)then
          ! Use atomic concentrations instead of smooth level set functions

          State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock) = &
               max(0.0, State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock))

          if( State_VGB(LevelXe_,i,j,k,iBlock) > 0.0) &
               State_VGB(LevelXe_,i,j,k,iBlock) = 1.0 / cAtomicMass_I(54)

          if( State_VGB(LevelBe_,i,j,k,iBlock) > 0.0) &
               State_VGB(LevelBe_,i,j,k,iBlock) = 1.0 / cAtomicMass_I(4)

          if( State_VGB(LevelPl_,i,j,k,iBlock) > 0.0) &
               State_VGB(LevelPl_,i,j,k,iBlock) = 1.0 / cAPolyimide
       end if

       ! Multiply level set functions with density unless the 
       ! non-conservative approach is used
       if(.not.UseUserSource) &
            State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock) = &
            State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock) &
            *State_VGB(Rho_,i,j,k,iBlock)

       ! Calculate internal energy for Xe/Be from pressure and density
       if(  State_VGB(Rho_,i,j,k,iBlock) < 1500.0*Io2No_V(UnitRho_) .and. &
           (State_VGB(LevelXe_,i,j,k,iBlock) > 0.0 .or. &
            State_VGB(LevelBe_,i,j,k,iBlock) > 0.0)) then 

          RhoSi = State_VGB(Rho_,i,j,k,iBlock)*No2Si_V(UnitRho_)
          pSi   = State_VGB(p_,i,j,k,iBlock)*No2Si_V(UnitP_)
          if(State_VGB(LevelBe_,i,j,k,iBlock) > 0.0)then
             iMaterial = Be_
          else
             iMaterial = Xe_
          end if

          call pressure_to_eint(pSi, RhoSi, iMaterial, &
               uDensityTotalOut=EinternalSi, IsError=IsError)

          if(IsError)then
             write(*,*) NameSub,' i,j,k,iBlock,iProc=',i,j,k,iBlock,iProc
             write(*,*) NameSub,' x,y,z=',x_BLK(i,j,k,iBlock), &
                  y_BLK(i,j,k,iBlock), z_BLK(i,j,k,iBlock)
             write(*,*) NameSub,' pSI, RhoSi, Material=', pSi, RhoSi, Be_
             call CON_stop(NameSub//': returned error from eos function')
          end if

          State_VGB(ExtraEInt_,i,j,k,iBlock) = &
               EInternalSi*Si2No_V(UnitEnergyDens_) &
               - inv_gm1*State_VGB(P_,i,j,k,iBlock)
       else
          State_VGB(ExtraEInt_,i,j,k,iBlock) = 0.0
       end if

    end do; end do; end do

  end subroutine user_set_ics

  !============================================================================

  subroutine read_hyades_file

    use ModIoUnit, ONLY: UnitTmp_
    use ModPhysics, ONLY: Si2No_V, UnitX_, UnitRho_, UnitU_, UnitP_

    integer :: iError, iCell  
    integer :: nStepHyades, nDimHyades, nEqparHyades, nVarHyades
    real    :: TimeHyades
    real, allocatable:: EqparHyades_I(:), Hyades2No_V(:)
    character(len=100):: StringHeadHyades, NameVarHyades

    character(len=*), parameter :: NameSub = "ModUser::read_hyades_file"
    !-------------------------------------------------------------------------
    open(UnitTmp_, FILE=NameHyadesFile, STATUS="old", IOSTAT=iError)

    if(iError /= 0)call stop_mpi(NameSub // &
         " could not open Hyades file="//NameHyadesFile)

    read(UnitTmp_, "(a)") StringHeadHyades
    read(UnitTmp_, *) &
         nStepHyades, TimeHyades, nDimHyades, nEqparHyades, nVarHyades
    if(nDimHyades /= 1)call stop_mpi(NameSub // &
         " only 1D Hyades file can be read now")

    read(UnitTmp_,*) nCellHyades

    allocate(EqparHyades_I(nEqparHyades))
    read(UnitTmp_,*) EqparHyades_I

    read(UnitTmp_, "(a)") NameVarHyades

    ! Set conversion from Hyades units to normalized units
    allocate(Hyades2No_V(nVarHyades))
    Hyades2No_V = 1.0
    Hyades2No_V(iRhoHyades) = 1000.0 * Si2No_V(UnitRho_) ! g/cm3 -> kg/m3
    Hyades2No_V(iUxHyades)  = 0.01   * Si2No_V(UnitU_)   ! cm/s  -> m/s
    Hyades2No_V(iPHyades)   = 0.1    * Si2No_V(UnitP_)   ! dyne  -> Pa

    allocate(xHyades_C(nCellHyades), StateHyades_VC(nVarHyades, nCellHyades))
    do iCell = 1, nCellHyades
       read(UnitTmp_, *) xHyades_C(iCell), StateHyades_VC(:, iCell)

       ! Convert from CGS to normalized units
       xHyades_C(iCell) = xHyades_C(iCell) * 0.01 * Si2No_V(UnitX_)
       StateHyades_VC(:, iCell) = StateHyades_VC(:, iCell) * Hyades2No_V


       ! Locate the Be-Xe interface based on the ionization level going 
       ! through 5 (Be ionization level can be at most 4)
       if( xBeHyades < 0.0 .and. iCell > 1 )then
          if(  StateHyades_VC(iZHyades, iCell-1) < 5.0 .and.  &
               StateHyades_VC(iZHyades, iCell)   > 5.0 ) &
               xBeHyades = 0.5*(xHyades_C(iCell-1) + xHyades_C(iCell))
       end if
    end do

    close(UnitTmp_)
    deallocate(EqparHyades_I)

    if(xBeHyades < 0.0)call stop_mpi(NameSub // &
         ' could not find Be-Xe interface based on ionization levels')

  end subroutine read_hyades_file

  !============================================================================
  subroutine user_update_states(iStage,iBlock)

    use ModProcMH,  ONLY: iProc
    use ModVarIndexes
    use ModSize
    use ModAdvance, ONLY: State_VGB, Rho_, RhoUy_, p_, ExtraEInt_, &
         LevelXe_, LevelPl_, Flux_VX, Flux_VY, Flux_VZ, Source_VC
    use ModGeometry,ONLY: x_BLK, y_BLK, z_BLK, vInv_CB
    use ModNodes,   ONLY: NodeY_NB
    use ModMain,    ONLY: nStage, Cfl
    use ModPhysics
    use ModEnergy,  ONLY: calc_energy_cell
    use ModEos,     ONLY: eos, eos_mixed_cell, UsePreviousTe

    implicit none

    integer, intent(in):: iStage,iBlock

    integer:: i, j, k, iMaterial, iMaterial_I(1)
    real   :: vInv_C(nI,nJ,nK)
    real   :: PressureSi, Einternal, EinternalSi, RhoSi, RhoToARatioSI_I(0:2)
    logical:: IsError = .false.

    character(len=*), parameter :: NameSub = 'user_update_states'
    !------------------------------------------------------------------------

    UsePreviousTe = .false.

    if(IsCylindrical)then
       ! Multiply fluxes with radius = abs(Y) at the X,Y and Z faces
       do k=1,nK; do j=1, nJ; do i=1, nI+1
          Flux_VX(:,i,j,k)=Flux_VX(:,i,j,k)*abs(y_BLK(i,j,k,iBlock))
       end do; end do; end do
       do k=1,nK; do j=1, nJ+1; do i=1, nI
          Flux_VY(:,i,j,k)=Flux_VY(:,i,j,k)*abs(NodeY_NB(i,j,k,iBlock))
       end do; end do; end do
       do k=1,nK+1; do j=1, nJ; do i=1, nI
          Flux_VZ(:,i,j,k)=Flux_VZ(:,i,j,k)*abs(y_BLK(i,j,k,iBlock))
       end do; end do; end do

       ! Add "geometrical source term" p/r to the radial momentum equation
       ! The "radial" direction is along the Y axis. There is no velocity
       ! in the azimuthal (=Z) direction, so there are no more terms.
       ! NOTE: here we have to use signed radial distance!

       do k=1,nK; do j=1, nJ; do i=1, nI
          Source_VC(RhoUy_,i,j,k) = Source_VC(RhoUy_,i,j,k) &
               + State_VGB(P_,i,j,k,iBlock) / y_BLK(i,j,k,iBlock)
       end do; end do; end do

       ! Multiply volume with radius (=Y) at cell center -> divide inverse vol
       vInv_C = vInv_CB(:,:,:,iBlock)
       do k=1,nK; do j=1, nJ; do i=1, nI
          vInv_CB(i,j,k,iBlock)=vInv_C(i,j,k)/abs(y_BLK(i,j,k,iBlock))
       end do; end do; end do

    end if

    call update_states_MHD(iStage,iBlock)
    
    ! Undo change of volume (fluxes and sources are not used any more)
    if(IsCylindrical) vInv_CB(:,:,:,iBlock) = vInv_C

    ! update of pressure and relaxation energy::

    do k=1,nK; do j=1,nJ; do i=1,nI
       ! Total internal energy ExtraEInt + P/(\gamma -1) transformed to SI
       EInternalSI = No2Si_V(UnitEnergyDens_)*&
            (inv_gm1*State_VGB(P_,i,j,k,iBlock) + &
            State_VGB(ExtraEInt_,i,j,k,iBlock))

       ! Density, transformed to SI
       RhoSI = No2Si_V(UnitRho_)*State_VGB(Rho_,i,j,k,iBlock)

       ! Find maximum level set value. 
       ! Note that for now plastic is now taken as carbon
       iMaterial_I = maxloc(State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock))
       iMaterial = iMaterial_I(1) - 1

       if( UseMixedCell .and. &
            maxval(State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock)) < &
            MixLimit * sum(State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock)) ) then
          ! The cell is mixed if none of the material is dominant
          RhoToARatioSI_I = &
               State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock) * No2Si_V(UnitRho_)
          call eos_mixed_cell(&
               UDensityTotal=EInternalSI, RhoToARatio_I=RhoToARatioSI_I,& 
               PTotalOut=PressureSI) !!! , IsError=IsError)
          if(IsError)write(*,*) NameSub,' EintSI,RhoToARatioSI_I,Material=',&
               EInternalSI, RhoToARatioSI_I, iMaterial
       else
          ! Get pressure from EOS
          call eos(UDensityTotal=EInternalSI, Rho=RhoSI, iMaterial=iMaterial, &
               pTotalOut=PressureSI) !!! , IsError=IsError)
          if(IsError)write(*,*) NameSub,' EintSi, RhoSi, iMaterial=',&
               EInternalSI, RhoSi, iMaterial
       end if
       if(IsError)then
          write(*,*) NameSub,' i,j,k,iBlock,iProc=',i,j,k,iBlock,iProc
          write(*,*) NameSub,' x,y,z=',x_BLK(i,j,k,iBlock), &
               y_BLK(i,j,k,iBlock), z_BLK(i,j,k,iBlock)
          call CON_stop(NameSub//': returned error from eos function')
       end if

       ! Set pressure and ExtraEInt = Total internal energy - P/(gamma -1)
       State_VGB(P_,i,j,k,iBlock) = PressureSI*Si2No_V(UnitP_)
       State_VGB(ExtraEInt_,i,j,k,iBlock) = Si2No_V(UnitEnergyDens_)*&
            (EInternalSI - PressureSI*inv_gm1)

    end do; end do; end do

    call calc_energy_cell(iBlock)

  end subroutine user_update_states

  !===========================================================================

  subroutine user_calc_sources

    use ModMain,     ONLY: nI, nJ, nK, GlobalBlk
    use ModAdvance,  ONLY: State_VGB, LevelXe_, LevelPl_, &
         Source_VC, uDotArea_XI, uDotArea_YI, uDotArea_ZI
    use ModGeometry, ONLY: vInv_CB

    character (len=*), parameter :: NameSub = 'user_calc_sources'
    integer :: i, j, k, iBlock
    !-------------------------------------------------------------------

    iBlock = globalBlk

    ! Add Level*div(u) as a source term so level sets beome advected scalars
    ! Note that all levels use the velocity of the first (and only) fluid
    
    do k=1,nK; do j=1,nJ; do i=1,nI
       Source_VC(LevelXe_:LevelPl_,i,j,k) =                 &
            Source_VC(LevelXe_:LevelPl_,i,j,k)              &
            + State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock)     &
            * vInv_CB(i,j,k,iBlock)*                        &
            ( uDotArea_XI(i+1,j,k,1) - uDotArea_XI(i,j,k,1) &
            + uDotArea_YI(i,j+1,k,1) - uDotArea_YI(i,j,k,1) &
            + uDotArea_ZI(i,j,k+1,1) - uDotArea_ZI(i,j,k,1))
    end do; end do; end do


  end subroutine user_calc_sources

  !===========================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModSize,    ONLY: nI, nJ, nK
    use ModAdvance, ONLY: State_VGB, LevelXe_, LevelPl_

    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(-1:nI+2, -1:nJ+2, -1:nK+2)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    character (len=*), parameter :: Name='user_set_plot_var'

    integer :: i, j, k, iMaterial_I(1)
    !------------------------------------------------------------------------  
    IsFound = NameVar == 'level'
    if(.not. IsFound) RETURN

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0
    
    do k=-1, nK+1; do j=-1, nJ+1; do i=-1,nI+2
       iMaterial_I = maxloc(State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock))
       PlotVar_G(i,j,k) = iMaterial_I(1)
    end do; end do; end do

  end subroutine user_set_plot_var

  !===========================================================================

  subroutine user_init_session

    use ModVarIndexes, ONLY: LevelXe_, LevelPl_, Rho_, UnitUser_V
    character (len=*), parameter :: NameSub = 'user_init_session'
    !-------------------------------------------------------------------

    if(UseUserSource)then
       UnitUser_V(LevelXe_:LevelPl_) = 1.e-6 ! = No2Io_V(UnitX_) = micron
    else if(UseMixedCell) then
       UnitUser_V(LevelXe_:LevelPl_) = UnitUser_V(Rho_)
    else
       UnitUser_V(LevelXe_:LevelPl_) = UnitUser_V(Rho_)*1.e-6
    end if

  end subroutine user_init_session

end module ModUser
