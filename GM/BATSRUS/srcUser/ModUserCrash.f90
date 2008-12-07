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
       IMPLEMENTED7 => user_set_ics,                    &
       IMPLEMENTED8 => user_material_properties

  use ModMain, ONLY: iTest, jTest, kTest, BlkTest, ProcTest, VarTest, &
       UseUserInitSession, UseUserIcs, UseUserSource, UseUserUpdateStates

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.2
  character (len=*), parameter :: &
       NameUserModule = 'HYDRO + IONIZATION EQUILIBRIUM + LEVEL SETS'

  ! Wall parameters
  real :: xEndTube   =   40.0    ! x coordinate of tube ending
  real :: rInnerTube =  287.5    ! inner radius [micron]
  real :: rOuterTube =  312.5    ! outer radius [micron]
  real :: RhoDimTube = 1430.0    ! density      [kg/m3]
  real :: RhoDimOutside = 6.5    ! density  of Xe outside tube [kg/m3]
  real :: pDimOutside   = 1.1e5  ! pressure of Xe outside tube [Pa]

  ! Allow overwriting the Xe state inside the tube for x > xUniformXe > 0
  real :: xUniformXe = -1.0

  ! Description of gold washer around the tube
  logical :: UseGold = .false.
  real :: WidthGold  = 50.0      ! width   [micron]
  real :: RhoDimGold = 20000.0   ! density [kg/m3]

  ! Allow for 2D with cylindrical symmetry around the X axis
  logical :: IsCylindrical = .false.

  ! Treat cells near material interface as a mixture
  logical :: UseMixedCell = .false.
  
  ! Mixed material cell is assumed if the ratio of dominant to total
  ! atomic concentration is below MixLimit
  real :: MixLimit = 0.97

  ! Variables for Hyades file
  logical           :: UseHyadesFile   = .false. ! read Hyades file?
  character(len=100):: NameHyadesFile            ! name of hyades file
  integer           :: nDimHyades      = -1      ! number of dimensions 
  integer           :: nVarHyades      = -1      ! number of variables
  integer           :: nCellHyades     = -1      ! number of cells
  integer           :: iCellLastHyades = -1      ! cell with maximum X and r=0
  real              :: xBeHyades       = -1.0    ! position of Be-Xe interface
  real, allocatable :: DataHyades_VC(:,:)        ! cell centered Hyades data
  real, allocatable :: LevelHyades_VC(:,:)       ! level set functions
  integer           :: iXHyades        = -1      ! index of x coordinate
  integer           :: iYHyades        = -1      ! index of y coordinate
  integer           :: iRhoHyades      = -1      ! index of density
  integer           :: iUxHyades       = -1      ! index of x velocity
  integer           :: iUyHyades       = -1      ! index of y velocity
  integer           :: iPHyades        = -1      ! index of pressure
  integer           :: iZHyades        = -1      ! index of ionization level
  integer           :: iTeHyades       = -1      ! index of electron temper.
  integer           :: iTrHyades       = -1      ! index of rad. temperature
  integer           :: iMaterialHyades = -1      ! index of material type

  ! Variables related to radiation
  character(len=20) :: TypeOpacity="constant"
  real :: RosselandOpacity(0:2) = 1.0
  real :: PlanckOpacity(0:2) = 10.0

  ! Indexes for lookup tables
  integer:: iTableRhoE = -1, iTableRhoP = -1

contains

  !============================================================================
  subroutine user_read_inputs

    use ModReadParam
    use ModFermiGas,ONLY:read_fermi_gas_param
    character (len=100) :: NameCommand
    character(len=*), parameter :: NameSub = 'user_read_inputs'
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
          call read_var('xEndTube',   xEndTube)
          call read_var('rInnerTube', rInnerTube)
          call read_var('rOuterTube', rOuterTube)
          call read_var('RhoDimTube', RhoDimTube)
          call read_var('RhoDimOutside', RhoDimOutside)
          call read_var('pDimOutside',   pDimOutside)
          call read_var('xUniformXe',    xUniformXe)
       case("#GOLD")
          call read_var('UseGold',    UseGold)
          call read_var('WidthGold',  WidthGold)
          call read_var('RhoDimGold', RhoDimGold)
       case("#MIXEDCELL")
          call read_var('UseMixedCell', UseMixedCell)
          if(UseMixedCell)call read_var('MixLimit', MixLimit)
       case("#CYLINDRICAL")
          call read_var('IsCylindrical', IsCylindrical)
       case("#FERMIGASEFFECT")
          call read_fermi_gas_param
       case("#OPACITY")
          call read_var('TypeOpacity', TypeOpacity)
          select case(TypeOpacity)
          case("constant")
             call read_var('PlanckOpacityXe', PlanckOpacity(0))
             call read_var('PlanckOpacityBe', PlanckOpacity(1))
             call read_var('PlanckOpacityPl', PlanckOpacity(2))
             call read_var('RosselandOpacityXe', RosselandOpacity(0))
             call read_var('RosselandOpacityBe', RosselandOpacity(1))
             call read_var('RosselandOpacityPl', RosselandOpacity(2))
          case default
             call stop_mpi(NameSub//"Wrong TypeOpacity ="//trim(TypeOpacity))
          end select
       case('#USERINPUTEND')
          EXIT
       case default
          call stop_mpi('ERROR in ModUserCrash: unknown command='//NameCommand)
       end select
    end do

  end subroutine user_read_inputs

  !============================================================================
  subroutine user_set_ics

    use ModProcMH,    ONLY: iProc
    use ModMain,      ONLY: GlobalBlk, nI, nJ, nK, UseGrayDiffusion
    use ModPhysics,   ONLY: inv_gm1, ShockPosition, ShockSlope, &
         Io2No_V, No2Si_V, Si2No_V, UnitRho_, UnitP_, UnitEnergyDens_
    use ModAdvance,   ONLY: State_VGB, Rho_, RhoUx_, RhoUz_, p_, &
         ExtraEint_, LevelBe_, LevelXe_, LevelPl_, Eradiation_
    use ModGeometry,  ONLY: x_BLK, y_BLK, z_BLK
    use ModEos,       ONLY: eos, Be_, Xe_, Plastic_
    use ModPolyimide, ONLY: cAtomicMass_I, cAPolyimide
    use ModConst,     ONLY: cRadiation
    use ModLookupTable, ONLY: interpolate_lookup_table

    real    :: x, y, xBe, DxBe, DxyPl, pSi, RhoSi, EinternalSi, TeSi
    real    :: e_I(Xe_:Plastic_)
    real    :: DxyGold = -1.0
    logical :: IsError

    integer :: iBlock, i, j, k, iMaterial, iMaterial_I(1)

    character(len=*), parameter :: NameSub = "user_set_ics"
    !------------------------------------------------------------------------

    iBlock = GlobalBlk

    if(UseHyadesFile)then
       ! Read in and interpolate Hyades output
       if(.not.allocated(DataHyades_VC)) call read_hyades_file

       if(nDimHyades == 1)then
          call interpolate_hyades1d(iBlock)
       else
          call interpolate_hyades2d(iBlock)
       end if
    end if

    ! Set level set functions, internal energy, and other values
    do k=-1, nK+2; do j=-1, nJ+2; do i=-1, nI+2 

       x = x_BLK(i,j,k,iBlock)
       y = y_BLK(i,j,k,iBlock)

       ! Distance from plastic wall: 
       !     positive for rInnerTube < |y| < rOuterTube and x > xEndTube only
       DxyPl = min(abs(y) - rInnerTube, rOuterTube - abs(y), x - xEndTube)

       ! Set plastic tube state
       if(DxyPl > 0.0)then
          ! Use the density and pressure given by the #TUBE command
          State_VGB(Rho_,i,j,k,iBlock) = RhoDimTube*Io2No_V(UnitRho_)
          State_VGB(p_  ,i,j,k,iBlock) = pDimOutside*Io2No_V(UnitP_)
          ! Assume that plastic wall is at rest
          State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
       end if

       ! Set pressure and speed outside the tube. 
       ! For 2D Hyades input do not overwrite values left of xEndTube
       if(abs(y) > rOuterTube &
            .and. (nDimHyades == 1 .or. x > xEndTube) ) then
          State_VGB(p_,i,j,k,iBlock) = pDimOutside*Io2No_V(UnitP_)
          State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
       end if

       ! Set the Xenon state inside the tube for x > xUniformXe if it is set
       if(xUniformXe > 0.0 .and. x > xUniformXe .and. abs(y) < rInnerTube)then
          State_VGB(Rho_,i,j,k,iBlock) = RhoDimOutside*Io2No_V(UnitP_)
          State_VGB(p_,i,j,k,iBlock)   = pDimOutside*Io2No_V(UnitP_)
          State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
       end if

       ! Set density outside the tube
       if(x >= xEndTube .and. abs(y) > rOuterTube) &
            State_VGB(Rho_,i,j,k,iBlock) = RhoDimOutside*Io2No_V(UnitRho_)
       
       ! Distance from gold washer xEndTube < x < xEndTube + WidthGold
       if(UseGold) DxyGold = &
            min(x - xEndTube, xEndTube + WidthGold - x, abs(y) - rOuterTube)

       ! Set density of gold washer (if present)
       if(DxyGold > 0.0) &
            State_VGB(Rho_,i,j,k,iBlock) = RhoDimGold*Io2No_V(UnitRho_)

       if(nDimHyades == 2)then
          State_VGB(LevelPl_,i,j,k,iBlock) = &
               min(-DxyGold, max(State_VGB(LevelPl_,i,j,k,iBlock), DxyPl))
          State_VGB(LevelXe_,i,j,k,iBlock) = &
               max(DxyGold, min(State_VGB(LevelXe_,i,j,k,iBlock), -DxyPl))
          State_VGB(LevelBe_,i,j,k,iBlock) = &
               min(-DxyGold, -DxyPl, State_VGB(LevelBe_,i,j,k,iBlock))
       else

          if(UseHyadesFile)then
             ! Be - Xe interface is given by Hyades file
             xBe = xBeHyades
          else
             ! Be - Xe interface is at the shock defined by #SHOCKPOSITION
             xBe = ShockPosition - ShockSlope*y
          end if

          ! Distance from Be disk: positive for x < xBe
          DxBe = xBe - x

          ! Plastic is between rInnerTube and rOuterTube and right of xEndTube
          State_VGB(LevelPl_,i,j,k,iBlock) = DxyPl

          ! Berylium is left of xBe inside rInnerTube 
          ! and it is left of xEndTube outside
          State_VGB(LevelBe_,i,j,k,iBlock) = &
               max(xEndTube - x, min(DxBe, rInnerTube - abs(y)))

          ! Xenon is right of xBe inside rInnerTube and 
          ! right of xEndTube outside rOuterTube
          State_VGB(LevelXe_,i,j,k,iBlock) = max( &
               min( x - xEndTube, abs(y) - rOuterTube), &
               min( -DxBe, rInnerTube - abs(y)) )

       end if

       if(UseMixedCell)then
          ! Use atomic concentrations instead of smooth level set functions

          if(maxval( State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock) ) <= 0.0)then
             State_VGB(LevelXe_,i,j,k,iBlock) = 1.0 / (3*cAtomicMass_I(54))
             State_VGB(LevelBe_,i,j,k,iBlock) = 1.0 / (3*cAtomicMass_I(4))
             State_VGB(LevelPl_,i,j,k,iBlock) = 1.0 / (3*cAPolyimide)
          else
             State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock) = &
                  max(0.0, State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock))

             if( State_VGB(LevelXe_,i,j,k,iBlock) > 0.0) &
                  State_VGB(LevelXe_,i,j,k,iBlock) = 1.0 / cAtomicMass_I(54)

             if( State_VGB(LevelBe_,i,j,k,iBlock) > 0.0) &
                  State_VGB(LevelBe_,i,j,k,iBlock) = 1.0 / cAtomicMass_I(4)

             if( State_VGB(LevelPl_,i,j,k,iBlock) > 0.0) &
                  State_VGB(LevelPl_,i,j,k,iBlock) = 1.0 / cAPolyimide
          end if

       end if

       ! Multiply level set functions with density unless the 
       ! non-conservative approach is used
       if(.not.UseUserSource) &
            State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock) = &
            State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock) &
            *State_VGB(Rho_,i,j,k,iBlock)

       ! Calculate internal energy from pressure and density
       iMaterial_I = maxloc(State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock))
       iMaterial = iMaterial_I(1) - 1

       RhoSi = State_VGB(Rho_,i,j,k,iBlock)*No2Si_V(UnitRho_)
       pSi   = State_VGB(p_,i,j,k,iBlock)*No2Si_V(UnitP_)

       ! The IsError flag avoids stopping for Fermi degenerated state
       if(iTableRhoP > 0)then
          call interpolate_lookup_table(iTableRhoP, RhoSi, pSi, E_I)
          EinternalSi = e_I(iMaterial)
       else
          call eos(iMaterial,RhoSi,pTotalIn=pSi, &
               ETotalOut=EinternalSi, IsError=IsError)
       end if

       State_VGB(ExtraEInt_,i,j,k,iBlock) = &
            EInternalSi*Si2No_V(UnitEnergyDens_) &
            - inv_gm1*State_VGB(P_,i,j,k,iBlock)

!!!       if(UseGrayDiffusion)then
!!!          ! The IsError flag avoids stopping for Fermi degenerated state
!!!          call eos(iMaterial, RhoSi, pTotalIn=pSi, TeOut=TeSi, IsError=IsError)
!!!
!!!          State_VGB(Eradiation_,i,j,k,iBlock) = cRadiation*TeSi**4 &
!!!               *Si2No_V(UnitEnergyDens_)
!!!       end if

    end do; end do; end do

  end subroutine user_set_ics

  !============================================================================

  subroutine read_hyades_file

    use ModIoUnit,    ONLY: UnitTmp_
    use ModPhysics,   ONLY: Si2No_V, Io2No_V, UnitX_, UnitRho_, UnitU_, &
         UnitP_, UnitTemperature_
    use ModUtilities, ONLY: split_string
    use ModEos,       ONLY: Xe_, Be_, Plastic_
    use ModConst,     ONLY: cKevToK
    use ModMain,      ONLY: UseGrayDiffusion

    integer             :: nStepHyades, nEqparHyades
    integer, allocatable:: nCellHyades_D(:)
    real                :: TimeHyades
    real, allocatable   :: EqparHyades_I(:), Hyades2No_V(:)
    character(len=100)  :: StringHeadHyades, NameVarHyades

    ! Variables for reading in variable names
    integer, parameter:: MaxString = 20
    character(len=10) :: String_I(MaxString)
    integer           :: nString

    ! Variables for setting level set functions
    integer :: iError, i, iCell, iMaterial, jMaterial
    real    :: x, y
    integer, allocatable:: iMaterial_C(:)
    real,    allocatable:: Distance2_C(:)

    character(len=*), parameter :: NameSub = "ModUser::read_hyades_file"
    !-------------------------------------------------------------------------
    open(UnitTmp_, FILE=NameHyadesFile, STATUS="old", IOSTAT=iError)

    if(iError /= 0)call stop_mpi(NameSub // &
         " could not open Hyades file="//NameHyadesFile)

    read(UnitTmp_, "(a)") StringHeadHyades
    read(UnitTmp_, *) &
         nStepHyades, TimeHyades, nDimHyades, nEqparHyades, nVarHyades

    ! Ignor negative value (signaling distorted grid)
    nDimHyades = abs(nDimHyades)
    
    ! Read grid size
    allocate(nCellHyades_D(nDimHyades))
    read(UnitTmp_,*) nCellHyades_D
    nCellHyades = product(nCellHyades_D)

    ! Read equation parameters
    allocate(EqparHyades_I(nEqparHyades))
    read(UnitTmp_,*) EqparHyades_I

    ! Read coordinate, variable and eqpar names
    read(UnitTmp_, "(a)") NameVarHyades
    call split_string(NameVarHyades, MaxString, String_I, nString)

    ! Find the columns for the coordinates and variables
    do i = 1, nDimHyades + nVarHyades
       ! The first nDimHyades strings are for the coordinates
       select case(String_I(i))
       case('x')
          iXHyades   = i
       case('y', 'r')
          iYHyades   = i
       case('rho')
          iRhoHyades = i
       case('ux')
          iUxHyades  = i
       case('uy', 'ur')
          iUyHyades  = i
       case('p')
          iPHyades   = i
       case('te')
          iTeHyades  = i
       case('tr')
          iTrHyades  = i
       case('z')
          iZHyades   = i
       case('material')
          iMaterialHyades = i
       end select

    end do
    ! Check if every coordinate/variable has been found
    if(iRhoHyades < 0)call stop_mpi(NameSub// &
         ' could not find rho in '//trim(NameVarHyades))

    if(iPHyades < 0)call stop_mpi(NameSub// &
         ' could not find p in '//trim(NameVarHyades))
       
    if(iUxHyades < 0)call stop_mpi(NameSub// &
         ' could not find ux in '//trim(NameVarHyades))

    if(iZHyades < 0 .and. iMaterialHyades < 0) call stop_mpi(NameSub// &
         ' could not find neither z nor material in '//trim(NameVarHyades))

    if(nDimHyades > 1)then
       ! y, uy and material are needed in 2D
       if(iYHyades < 0) call stop_mpi(NameSub// &
            ' could not find y in '//trim(NameVarHyades))
       if(iUyHyades < 0) call stop_mpi(NameSub// &
            ' could not find uy in '//trim(NameVarHyades))
       if(iMaterialHyades < 0) call stop_mpi(NameSub// &
            ' could not find material in '//trim(NameVarHyades))
    end if

    ! Set conversion from Hyades units to normalized units
    allocate(Hyades2No_V(nDimHyades + nVarHyades))
    Hyades2No_V = 1.0
    Hyades2No_V(iXHyades)   = 0.01   * Si2No_V(UnitX_)   ! cm    -> m
    Hyades2No_V(iRhoHyades) = 1000.0 * Si2No_V(UnitRho_) ! g/cm3 -> kg/m3
    Hyades2No_V(iUxHyades)  = 0.01   * Si2No_V(UnitU_)   ! cm/s  -> m/s
    Hyades2No_V(iPHyades)   = 0.1    * Si2No_V(UnitP_)   ! dyne  -> Pa

    if(UseGrayDiffusion)then
       if(iTrHyades < 0) call stop_mpi(NameSub// &
            ' could not find radiation temperature in '//trim(NameVarHyades))
       if(iTeHyades < 0) call stop_mpi(NameSub// &
            ' could not find electron temperature in '//trim(NameVarHyades))

       Hyades2No_V(iTeHyades)= cKevToK* Si2No_V(UnitTemperature_) ! KeV   -> K
       Hyades2No_V(iTrHyades)= cKevToK* Si2No_V(UnitTemperature_) ! KeV   -> K
    end if
    if(nDimHyades > 1)then
       Hyades2No_V(iYHyades)  = 0.01 * Si2No_V(UnitX_)   ! cm    -> m
       Hyades2No_V(iUyHyades) = 0.01 * Si2No_V(UnitU_)   ! cm/s  -> m/s
    end if

    ! Read in the data
    allocate(DataHyades_VC(nDimHyades + nVarHyades, nCellHyades))
    do iCell = 1, nCellHyades
       read(UnitTmp_, *) DataHyades_VC(:, iCell)
       ! Convert from CGS to normalized units
       DataHyades_VC(:, iCell) = DataHyades_VC(:, iCell) * Hyades2No_V
    end do
    close(UnitTmp_)

    if(nDimHyades == 1)then

       ! Locate the Be-Xe interface in 1D 
       do iCell = 2, nCellHyades
          if(iMaterialHyades > 0)then
             ! Check if material changes from Be to Xe
             if(  nint(DataHyades_VC(iMaterialHyades, iCell-1)) == Be_ .and. &
                  nint(DataHyades_VC(iMaterialHyades, iCell  )) == Xe_) EXIT
          else
             ! Check if ionization level jumps through 5
             if(  DataHyades_VC(iZHyades, iCell-1) < 5.0 .and.  &
                  DataHyades_VC(iZHyades, iCell)   > 5.0 ) EXIT
          end if
       end do
       if(iCell > nCellHyades)call stop_mpi(NameSub // &
            ' could not find Be/Xe interface')

       xBeHyades = 0.5* &
            ( DataHyades_VC(iXHyades, iCell-1) &
            + DataHyades_VC(iXHyades, iCell))

    else

       ! Fix the pressure in the plastic cells
       where(DataHyades_VC(iMaterialHyades, :) == Plastic_ &
            .and. DataHyades_VC(iPHyades, :) < 1.0) &
               DataHyades_VC(iPHyades, :) = pDimOutside*Io2No_V(UnitP_)

       ! Find cell with maximum X coordinate along the symmetry axis
       iCellLastHyades = nCellHyades_D(1)
       write(*,*) 'Coordinates of rightmost cell along axis:', &
            DataHyades_VC(iXHyades,iCellLastHyades), &
            DataHyades_VC(iYHyades,iCellLastHyades)

       ! Calculate level set functions on the Hyades grid using 
       ! the minimum distance between cells of different materials
       allocate(LevelHyades_VC(Xe_:Plastic_, nCellHyades))
       allocate(Distance2_C(nCellHyades), iMaterial_C(nCellHyades))

       do iCell = 1, nCellHyades
          x         = DataHyades_VC(iXHyades, iCell)
          y         = DataHyades_VC(iYHyades, iCell)
          iMaterial = DataHyades_VC(iMaterialHyades, iCell)

          ! Distance squared from all other points
          Distance2_C = (x - DataHyades_VC(iXHyades,:))**2       &
               +        (y - DataHyades_VC(iYHyades,:))**2

          ! Integer value of material in Hyades grid
          iMaterial_C = DataHyades_VC(iMaterialHyades,:)

          ! For each cell set 3 level set functions
          do jMaterial = Xe_, Plastic_
             if(iMaterial == jMaterial)then
                ! Level is the smallest distance to a different material
                LevelHyades_VC(jMaterial, iCell) =  sqrt(minval &
                     ( Distance2_C, MASK=iMaterial_C /= jMaterial))
             else
                ! Level is -1 times the smallest distance to same material
                LevelHyades_VC(jMaterial, iCell) = - sqrt(minval &
                     ( Distance2_C, MASK=iMaterial_C == jMaterial))
             end if
          end do
       end do
       deallocate(Distance2_C, iMaterial_C)
    end if

    deallocate(EqparHyades_I)

  end subroutine read_hyades_file

  !============================================================================

  subroutine interpolate_hyades1d(iBlock)

    use ModSize,     ONLY: nI, nJ, nK
    use ModAdvance,  ONLY: State_VGB, Rho_, RhoUx_, RhoUy_, RhoUz_, p_, &
         Eradiation_
    use ModGeometry, ONLY: x_BLK
    use ModPhysics,  ONLY: Si2No_V, UnitEnergyDens_, UnitTemperature_, &
         cRadiationNo
    use ModMain,     ONLY: UseGrayDiffusion

    integer, intent(in) :: iBlock

    integer :: i, j, k, iCell
    real :: x, Weight1, Weight2
    real :: Tr
    character(len=*), parameter :: NameSub='interpolate_hyades1d'
    !-------------------------------------------------------------------------
    do i = -1, nI+2
       ! Find the Hyades points around this position
       x = x_Blk(i,1,1,iBlock)

       do iCell=1, nCellHyades
          if(DataHyades_VC(iXHyades, iCell) >= x) EXIT
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
          Weight1 = (DataHyades_VC(iXHyades, iCell) - x) &
               /    (DataHyades_VC(iXHyades, iCell) &
               -     DataHyades_VC(iXHyades, iCell-1))
          Weight2 = 1.0 - Weight1
       end if

       do k = -1,nk+2; do j = -1,nJ+2
          ! Interpolate density, momentum and pressure

          State_VGB(Rho_,i,j,k,iBlock) = &
               ( Weight1*DataHyades_VC(iRhoHyades, iCell-1) &
               + Weight2*DataHyades_VC(iRhoHyades, iCell) )

          State_VGB(RhoUx_,i,j,k,iBlock) =  State_VGB(Rho_,i,j,k,iBlock) * &
               ( Weight1*DataHyades_VC(iUxHyades, iCell-1) &
               + Weight2*DataHyades_VC(iUxHyades, iCell) )

          State_VGB(p_,i,j,k,iBlock) = &
               ( Weight1*DataHyades_VC(iPHyades, iCell-1) &
               + Weight2*DataHyades_VC(iPHyades, iCell) )

          ! Set transverse momentum to zero
          State_VGB(RhoUy_:RhoUz_,i,j,k,iBlock) = 0.0

          ! Radiation energy = cRadiation*Trad**4
          if(UseGrayDiffusion)then
             Tr = ( Weight1*DataHyades_VC(iTrHyades, iCell-1) &
                  + Weight2*DataHyades_VC(iTrHyades, iCell) )

             State_VGB(Eradiation_,i,j,k,iBlock) = cRadiationNo*Tr**4
          end if

       end do; end do
    end do

  end subroutine interpolate_hyades1d

  !============================================================================

  subroutine interpolate_hyades2d(iBlock)

    ! Use Delaunay triangulation to interpolate Hyades grid onto CRASH grid

    use ModSize,     ONLY: nI, nJ, nK
    use ModAdvance,  ONLY: State_VGB, Rho_, RhoUx_, RhoUy_, RhoUz_, p_, &
         LevelXe_, LevelPl_, Eradiation_
    use ModGeometry,    ONLY: x_BLK, y_BLK
    use ModTriangulate, ONLY: calc_triangulation, find_triangle
    use ModMain,        ONLY: UseGrayDiffusion
    use ModPhysics,     ONLY: cRadiationNo

    integer, intent(in) :: iBlock

    integer, save              :: nTriangle
    integer, allocatable, save :: iNodeTriangle_II(:,:)
    real, allocatable,    save :: DataHyades_V(:)

    integer :: i, j, k, iNode1, iNode2, iNode3
    real    :: x, y, Xy_D(2), Weight1, Weight2, Weight3

    character(len=*), parameter :: NameSub='interpolate_hyades2d'
    !-------------------------------------------------------------------------
    if(.not.allocated(iNodeTriangle_II))then
       ! allocate variables and do triangulation
       allocate(iNodeTriangle_II(3,2*nCellHyades))
       allocate(DataHyades_V(nDimHyades + nVarHyades))
       call calc_triangulation( &
            nCellHyades, DataHyades_VC( (/iXHyades, iYHyades/), :), &
            iNodeTriangle_II, nTriangle)
    end if

    ! Interpolate points 
    do j = -1, nJ+2; do i = -1, nI+2

       x = x_Blk(i,j,1,iBlock)
       y = y_Blk(i,j,1,iBlock)

       ! abs(y) = radial distance in Hyades
       Xy_D = (/ x, abs(y) /)

       ! Check if we are at the end of the Hyades grid
       if(x >= DataHyades_VC(iXHyades, iCellLastHyades))then
          iNode1 = iCellLastHyades;  Weight1 = 1.0
          iNode2 = 1;                Weight2 = 0.0
          iNode3 = 1;                Weight3 = 0.0
       else
          ! Find the Hyades triangle around this position
          call find_triangle(&
               nCellHyades, nTriangle, &
               Xy_D, DataHyades_VC( (/iXHyades, iYHyades/),:), &
               iNodeTriangle_II(:,1:nTriangle), &
               iNode1, iNode2, iNode3, Weight1, Weight2, Weight3)
       end if

       DataHyades_V = &
            Weight1*DataHyades_VC(:, iNode1) + &
            Weight2*DataHyades_VC(:, iNode2) + &
            Weight3*DataHyades_VC(:, iNode3)

       do k = -1, nk+2
          ! Interpolate density, momentum and pressure

          State_VGB(Rho_,i,j,k,iBlock)  = DataHyades_V(iRhoHyades)

          State_VGB(p_,i,j,k,iBlock)    = DataHyades_V(iPHyades)

          State_VGB(RhoUx_,i,j,k,iBlock) = &
               DataHyades_V(iRhoHyades) * DataHyades_V(iUxHyades)

          State_VGB(RhoUy_,i,j,k,iBlock) = sign(1.0, y) * &
               DataHyades_V(iRhoHyades) * DataHyades_V(iUyHyades)

          ! Set azimuthal momentum to zero
          State_VGB(RhoUz_,i,j,k,iBlock) = 0.0

          ! Interpolate level set functions
          State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock) = &
               Weight1*LevelHyades_VC(:, iNode1) + &
               Weight2*LevelHyades_VC(:, iNode2) + &
               Weight3*LevelHyades_VC(:, iNode3)


          ! Radiation energy = cRadiation*Trad**4
          if(UseGrayDiffusion) State_VGB(Eradiation_,i,j,k,iBlock) = &
               cRadiationNo * DataHyades_V(iTrHyades)**4

       end do

    end do; end do

  end subroutine interpolate_hyades2d

  !============================================================================

  subroutine user_update_states(iStage,iBlock)

    use ModProcMH,  ONLY: iProc
    use ModVarIndexes
    use ModSize
    use ModAdvance, ONLY: State_VGB, Rho_, RhoUy_, p_, ExtraEInt_, &
         LevelXe_, LevelPl_, Flux_VX, Flux_VY, Flux_VZ, Source_VC, &
         VdtFace_Y, VdtFace_Z
    use ModGeometry,ONLY: x_BLK, y_BLK, z_BLK, vInv_CB
    use ModNodes,   ONLY: NodeY_NB
    use ModPhysics
    use ModEnergy,  ONLY: calc_energy_cell
    use ModEos,     ONLY: eos, Xe_, Plastic_
    use ModLookupTable, ONLY: interpolate_lookup_table

    implicit none

    integer, intent(in):: iStage,iBlock

    integer:: i, j, k, iMaterial, iMaterial_I(1)
    real   :: vInv_C(nI,nJ,nK)
    real   :: PressureSi, EinternalSi, RhoSi
    real   :: RhoToARatioSI_I(Xe_:Plastic_), p_I(Xe_:Plastic_)

    character(len=*), parameter :: NameSub = 'user_update_states'
    !------------------------------------------------------------------------

    if(IsCylindrical)then
       ! Multiply fluxes with radius = abs(Y) at the X and Y faces
       do k=1,nK; do j=1, nJ; do i=1, nI+1
          Flux_VX(:,i,j,k)=Flux_VX(:,i,j,k)*abs(y_BLK(i,j,k,iBlock))
       end do; end do; end do
       do k=1,nK; do j=1, nJ+1; do i=1, nI
          Flux_VY(:,i,j,k)=Flux_VY(:,i,j,k)*abs(NodeY_NB(i,j,k,iBlock))
          ! Upper estimate on the time step restriction takes the smaller r
          VdtFace_Y(i,j,k) = VdtFace_Y(i,j,k)* abs(NodeY_NB(i,j,k,iBlock))/&
                  min(abs(y_BLK(i,j,k,iBlock)),abs(y_BLK(i,j-1,k,iBlock)))
       end do; end do; end do
       ! There are no fluxes and CFL restrictions in the azimuthal direction
       do k=1,nK+1; do j=1, nJ; do i=1, nI
          Flux_VZ(:,i,j,k) = 0.0
          VdtFace_Z(i,j,k) = 0.0
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
          call eos(&
               RhoToARatioSI_I,ETotalIn=EInternalSI,& 
               PTotalOut=PressureSI) 
       else
          ! Get pressure from EOS
          if(iTableRhoE > 0)then
             call interpolate_lookup_table(iTableRhoE, RhoSi, EinternalSi, p_I)
             PressureSi = p_I(iMaterial)
          else
             call eos(iMaterial, Rho=RhoSI,ETotalIn=EInternalSI, &
                  pTotalOut=PressureSI)
          end if
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

    use ModConst,   ONLY: cKtoKev
    use ModSize,    ONLY: nI, nJ, nK
    use ModAdvance, ONLY: State_VGB, Rho_, p_, LevelXe_, LevelPl_, &
         Eradiation_
    use ModPhysics, ONLY: No2Si_V, UnitRho_, UnitP_, UnitTemperature_, &
         cRadiationNo
    use ModEos,     ONLY: eos, Plastic_

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

    real    :: p, Rho, pSi, RhoSi, TeSi
    integer :: i, j, k, iMaterial, iMaterial_I(1)
    logical :: IsError
    !------------------------------------------------------------------------  
    IsFound = .true.
    select case(NameVar)
    case('level', 'material')
       do k=-1, nK+1; do j=-1, nJ+1; do i=-1,nI+2
          iMaterial_I = maxloc(State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock))
          PlotVar_G(i,j,k) = iMaterial_I(1)
       end do; end do; end do
    case('tekev', 'TeKev')
       do k=-1, nK+1; do j=-1, nJ+1; do i=-1,nI+2
          Rho = State_VGB(Rho_,i,j,k,iBlock)
          p   = State_VGB(p_,i,j,k,iBlock)
          PlotVar_G(i,j,k) = p/Rho*No2Si_V(UnitTemperature_) * cKToKev
          iMaterial_I = maxloc(State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock))
          iMaterial   = iMaterial_I(1) - 1

          RhoSi = Rho*No2Si_V(UnitRho_)
          pSi   = p*No2Si_V(UnitP_)
          ! The IsError flag avoids stopping for Fermi degenerated state
          call eos(iMaterial, RhoSi, pTotalIn=pSi,TeOut=TeSi, &
               IsError=IsError)
          PlotVar_G(i,j,k) = TeSi * cKToKev
       end do; end do; end do
    case('tradkev','trkev')
       ! multiply by sign of Erad for debugging purpose
       PlotVar_G = sign(1.0,State_VGB(Eradiation_,:,:,:,iBlock)) &
            *sqrt(sqrt(abs(State_VGB(Eradiation_,:,:,:,iBlock))/cRadiationNo))&
            * No2Si_V(UnitTemperature_) * cKToKev
    case default
       IsFound = .false.
    end select

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0
    
  end subroutine user_set_plot_var

  !===========================================================================

  subroutine user_init_session

    use ModProcMH,      ONLY: iProc, iComm
    use ModVarIndexes,  ONLY: LevelXe_, LevelPl_, Rho_, UnitUser_V
    use ModLookupTable, ONLY: i_lookup_table, make_lookup_table

    character (len=*), parameter :: NameSub = 'user_init_session'
    !-------------------------------------------------------------------

    if(UseUserSource)then
       UnitUser_V(LevelXe_:LevelPl_) = 1.e-6 ! = No2Io_V(UnitX_) = micron
    else if(UseMixedCell) then
       UnitUser_V(LevelXe_:LevelPl_) = UnitUser_V(Rho_)
    else
       UnitUser_V(LevelXe_:LevelPl_) = UnitUser_V(Rho_)*1.e-6
    end if

    iTableRhoE = i_lookup_table('p(rho,e)')
    iTableRhoP = i_lookup_table('e(rho,p)')

    if(iProc==0) &
         write(*,*)NameSub,' iTableRhoE, iTableRhoP = ', iTableRhoE, iTableRhoP

    if(iTableRhoE > 0) &
         call make_lookup_table(iTableRhoE, calc_table_rho_e, iComm)
    if(iTableRhoP > 0) &
         call make_lookup_table(iTableRhoP, calc_table_rho_p, iComm)

  end subroutine user_init_session

  !===========================================================================
  subroutine calc_table_rho_e(Rho, e, p_I)

    use ModEos, ONLY: eos, Xe_, Plastic_

    ! Calculate pressures for Xe_, Be_ and Plastic_ for given Rho and e

    real, intent(in):: Rho, e
    real, intent(out):: p_I(:)

    integer:: iMaterial
    !-----------------------------------------------------------------------
    ! Material index starts from 0 :-( hence the +1
    do iMaterial = Xe_, Plastic_
       call eos(iMaterial, Rho, EtotalIn=e, pTotalOut=p_I(iMaterial+1))
    end do
    
  end subroutine calc_table_rho_e
  !===========================================================================
  subroutine calc_table_rho_p(Rho, p, e_I)

    use ModEos, ONLY: eos, Xe_, Plastic_

    ! Calculate energies for Xe_, Be_ and Plastic_ for given Rho and p

    real, intent(in):: Rho, p
    real, intent(out):: e_I(:)

    integer:: iMaterial
    !-----------------------------------------------------------------------
    ! Material index starts from 0 :-( hence the +1
    do iMaterial = Xe_, Plastic_
       call eos(iMaterial, Rho, PtotalIn=p, eTotalOut=e_I(iMaterial+1))
    end do
    
  end subroutine calc_table_rho_p
  !===========================================================================

  subroutine user_material_properties(State_V, &
       TeSi, AbsorptionOpacitySi, RosselandMeanOpacitySi)

    use ModEos,        ONLY: eos
    use ModPhysics,    ONLY: No2Si_V, UnitRho_, UnitP_
    use ModVarIndexes, ONLY: nVar, Rho_, LevelXe_, LevelPl_, p_

    real, intent(in) :: State_V(1:nVar)
    real, optional, intent(out) :: TeSi                   ! [K]
    real, optional, intent(out) :: AbsorptionOpacitySi    ! [1/m]
    real, optional, intent(out) :: RosselandMeanOpacitySi ! [1/m]

    character (len=*), parameter :: NameSub = 'user_material_properties'

    real    :: pSi, RhoSi, TemperatureSi
    integer :: iMaterial, iMaterial_I(1)
    logical :: IsError
    !-------------------------------------------------------------------------
    iMaterial_I = maxloc(State_V(LevelXe_:LevelPl_))
    iMaterial   = iMaterial_I(1) - 1

    RhoSi = State_V(Rho_)*No2Si_V(UnitRho_)
    pSi   = State_V(p_)*No2Si_V(UnitP_)
    ! The IsError flag avoids stopping for Fermi degenerated state
    call eos(iMaterial, RhoSi, pTotalIn=pSi, TeOut=TemperatureSi, &
         IsError=IsError)

    if(present(TeSi)) TeSi = TemperatureSi

    if(present(AbsorptionOpacitySi)) &
         AbsorptionOpacitySi = PlanckOpacity(iMaterial)*RhoSi

    if(present(RosselandMeanOpacitySi)) &
         RosselandMeanOpacitySi = RosselandOpacity(iMaterial)*RhoSi

  end subroutine user_material_properties

end module ModUser
