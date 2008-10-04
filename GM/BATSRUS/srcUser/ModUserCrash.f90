!^CFG COPYRIGHT UM
!========================================================================
module ModUser
  ! This is the default user module which contains empty methods defined
  ! in ModUserEmpty.f90

  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_update_states,              &
       IMPLEMENTED2 => user_calc_sources,               &
       IMPLEMENTED3 => user_initial_perturbation,       &
       IMPLEMENTED4 => user_read_inputs,                &
       IMPLEMENTED5 => user_set_plot_var,               &
       IMPLEMENTED6 => user_init_session

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.2
  character (len=*), parameter :: &
       NameUserModule = 'HYDRO + IONIZATION EQUILIBRIUM + LEVEL SETS'

  ! Wall parameters
  real :: rXe = 287.5, rPl = 312.5, RhoPlDim = 1430.0
  ! True if the plastic tube extends all the way to X=0, and Be is inside it
  logical :: IsFullTube = .false.

  ! Allow for 2D with cylindrical symmetry around the X axis
  logical :: IsCylindrical = .false.

  ! Treat cells near material interface as a mixture
  logical :: UseMixedCell = .false.
  
  ! Mixed material cell is assumed if the ratio of dominant to total
  ! atomic concentration is below MixLimit
  real :: MixLimit = 0.97

contains

  subroutine user_read_inputs

    use ModReadParam
    character (len=100) :: NameCommand
    !------------------------------------------------------------------------
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#TUBE")
          call read_var('IsFullTube', IsFullTube)
          call read_var('rInnerTube', rXe)
          call read_var('rOuterTube', rPl)
          call read_var('RhoTubeDim', RhoPlDim)
       case('#MIXEDCELL')
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
  subroutine user_update_states(iStage,iBlock)

    use ModProcMH,  ONLY: iProc
    use ModVarIndexes
    use ModSize
    use ModAdvance, ONLY: State_VGB, Rho_, RhoUy_, p_, ExtraEInt_, &
         LevelXe_, LevelPl_, Flux_VX, Flux_VY, Flux_VZ, time_BLK
    use ModGeometry,ONLY: x_BLK, y_BLK, z_BLK, vInv_CB
    use ModNodes,   ONLY: NodeY_NB
    use ModMain,    ONLY: nStage, Cfl
    use ModPhysics
    use ModEnergy,  ONLY: calc_energy_cell
    use ModEos,     ONLY: eos, eos_mixed_cell, UsePreviousTe

    implicit none

    integer, intent(in):: iStage,iBlock

    integer:: i, j, k, iMaterial, Loc_I(1)
    real   :: PressureSi, Einternal, EinternalSi, RhoSi, RhoToARatioSI_I(0:2)
    real   :: DtFactor
    logical:: IsError = .false.

    character(len=*), parameter :: NameSub = 'user_update_states'
    !------------------------------------------------------------------------

    UsePreviousTe = .false.

    if(IsCylindrical)then
       ! Multiply fluxes with radius (=Y) at face
       do k=1,nK; do j=1, nJ; do i=1, nI+1
          Flux_VX(:,i,j,k)=Flux_VX(:,i,j,k)*y_BLK(i,j,k,iBlock)
       end do; end do; end do
       do k=1,nK; do j=1, nJ+1; do i=1, nI
          Flux_VY(:,i,j,k)=Flux_VY(:,i,j,k)*NodeY_NB(i,j,k,iBlock)
       end do; end do; end do
       do k=1,nK+1; do j=1, nJ; do i=1, nI
          Flux_VZ(:,i,j,k)=Flux_VZ(:,i,j,k)*y_BLK(i,j,k,iBlock)
       end do; end do; end do

       ! Multiply volume with radius (=Y) at cell center
       do k=1,nK; do j=1, nJ; do i=1, nI
          vInv_CB(i,j,k,iBlock)=vInv_CB(i,j,k,iBlock)*y_BLK(i,j,k,iBlock)
       end do; end do; end do
    end if

    call update_states_MHD(iStage,iBlock)

    if(IsCylindrical)then
       ! Undo change of volume
       do k=1,nK; do j=1, nJ; do i=1, nI
          vInv_CB(i,j,k,iBlock)=vInv_CB(i,j,k,iBlock)/y_BLK(i,j,k,iBlock)
       end do; end do; end do

       ! Store DtFactor for adding geometrical source term
       DtFactor = iStage*(Cfl/nStage)
    end if

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
       Loc_I = maxloc(State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock))
       iMaterial = Loc_I(1) - 1

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

       ! Add "geometrical source term" p/r to the radial momentum equation
       ! The "radial" direction is along the Y axis. There is no velocity
       ! in the azimuthal (=Z) direction, so there are no more terms.
       if(IsCylindrical) &
            State_VGB(RhoUy_,i,j,k,iBlock) = State_VGB(RhoUy_,i,j,k,iBlock) &
            + DtFactor * time_BLK(i,j,k,iBlock) &
            * State_VGB(P_,i,j,k,iBlock) / y_BLK(i,j,k,iBlock)

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

  subroutine user_initial_perturbation

    use ModProcMH,  ONLY: iProc
    use ModMain,    ONLY: nI, nJ, nK, nBlock, UnusedBlk, UseUserSource
    use ModPhysics, ONLY: ShockPosition, ShockSlope, ShockRightState_V, &
         Io2No_V, No2Si_V, Si2No_V, UnitRho_, UnitP_, UnitEnergyDens_, &
         inv_gm1
    use ModAdvance, ONLY: State_VGB, Rho_, RhoUx_, RhoUz_, p_, ExtraEint_, &
         LevelBe_, LevelXe_, LevelPl_
    use ModGeometry,ONLY: x_BLK, y_BLK, z_BLK, y2
    use ModEnergy,  ONLY: calc_energy_ghost
    use ModEos,     ONLY: pressure_to_eint, Be_
    use ModPolyimide, ONLY: cAtomicMass_I, cAPolyimide

    integer :: i, j, k, iBlock
    real    :: x, y, xBeSlope, DxBe, DyPl, pSi, RhoSi, EinternalSi
    logical :: IsError

    character (len=*), parameter :: NameSub = 'user_initial_perturbation'
    !-------------------------------------------------------------------
  
    do iBlock = 1, nBlock

       if (unusedBLK(iBlock)) CYCLE

       do k=-1, nK+2; do j=-1, nJ+2; do i=-1, nI+2 

          x = x_BLK(i,j,k,iBlock)
          y = y_BLK(i,j,k,iBlock)

          ! Be - Xe interface is at the shock defined by #SHOCKPOSITION
          xBeSlope = ShockPosition - ShockSlope*y

          ! Distance from Be disk: positive for x < xBeSlope
          DxBe = xBeSlope - x

          ! Distance from plastic wall: positive for rXe < |y| < rPl only
          DyPl = min(abs(y) - rXe, rPl - abs(y))

          ! Set plastic wall state
          if((IsFullTube .or. DxBe < 0.0) .and. DyPl > 0.0)then
             ! Use the density and pressure given by the #MATERIAL command
             State_VGB(Rho_,i,j,k,iBlock) = RhoPlDim*Io2No_V(UnitRho_)
             State_VGB(p_  ,i,j,k,iBlock) = &
                  ShockRightState_V(p_)*Io2No_V(UnitP_)
             ! Assume that plastic wall is at rest
             State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
          end if

          ! Set Xe/Be pressure and speed outside rPl
          if(abs(y) > rPl) then
             State_VGB(p_,i,j,k,iBlock) = ShockRightState_V(p_)*Io2No_V(UnitP_)
             State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
          end if

          if(IsFullTube)then
             ! Berylium is inside plastic tube and left of "shock"
             State_VGB(LevelBe_,i,j,k,iBlock) = min(DxBe, -DyPl)
             ! Plastic is between rInnerPl and rOuterPl
             State_VGB(LevelPl_,i,j,k,iBlock) = DyPl
          else
             ! Berylium is everywhere left to the end of the plastic tube
             State_VGB(LevelBe_,i,j,k,iBlock) = DxBe
             ! Plastic is right of xBe and between rInnerTube and rOuterTube
             State_VGB(LevelPl_,i,j,k,iBlock) = min( -DxBe, DyPl)
          end if

          ! Xenon is right of xBe, inside rInnerTube and outside rOuterTube
          State_VGB(LevelXe_,i,j,k,iBlock) = min( -DxBe, -DyPl)

          if(UseMixedCell)then
             ! Use atomic concentrations instead of level set functions

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

          ! Calculate internal energy from pressure and density
          !!!
          if(.false. .and. State_VGB(LevelBe_,i,j,k,iBlock) > 0.0)then !!!
             RhoSi = State_VGB(Rho_,i,j,k,iBlock)*No2Si_V(UnitRho_)
             pSi   = State_VGB(p_,i,j,k,iBlock)*No2Si_V(UnitP_)
             call pressure_to_eint(pSi, RhoSi, Be_, &
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

       call calc_energy_ghost(iBlock)

    end do

  end subroutine user_initial_perturbation

  !====================================================================

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

    integer :: i, j, k, Loc_I(1)
    !------------------------------------------------------------------------  
    IsFound = NameVar == 'level'
    if(.not. IsFound) RETURN

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0
    
    do k=-1, nK+1; do j=-1, nJ+1; do i=-1,nI+2
       Loc_I = maxloc(State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock))
       PlotVar_G(i,j,k) = Loc_I(1)
    end do; end do; end do

  end subroutine user_set_plot_var

  !===========================================================================

  subroutine user_init_session

    use ModVarIndexes, ONLY: LevelXe_, LevelPl_, Rho_, UnitUser_V
    use ModMain, ONLY: UseUserSource
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
