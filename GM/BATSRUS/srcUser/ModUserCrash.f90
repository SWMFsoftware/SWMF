!^CFG COPYRIGHT UM
!========================================================================
module ModUser
  ! This is the default user module which contains empty methods defined
  ! in ModUserEmpty.f90

  use ModUserEmpty,                                      &
       IMPLEMENTED1  => user_update_states,              &
       IMPLEMENTED2  => user_calc_sources,               &
       IMPLEMENTED3  => user_read_inputs,                &
       IMPLEMENTED4  => user_set_plot_var,               &
       IMPLEMENTED5  => user_init_session,               &
       IMPLEMENTED6  => user_set_ics,                    &
       IMPLEMENTED7  => user_material_properties,        &
       IMPLEMENTED8  => user_amr_criteria,               &
       IMPLEMENTED9  => user_get_log_var,                &
       IMPLEMENTED10 => user_set_resistivity,            &
       IMPLEMENTED11 => user_action

  use ModMain, ONLY: &
       UseUserInitSession, UseUserIcs, UseUserSource, UseUserUpdateStates, &
       UseUserLogFiles
  use ModSize, ONLY: x_, y_, z_, nI, nJ, nK

  use ModVarIndexes, ONLY: nMaterial, MaterialFirst_, MaterialLast_, Pe_
  use CRASH_ModEos, ONLY: MassMaterial_I => cAtomicMassCRASH_I, &
       Be_, Plastic_, Au_, Ay_
  use BATL_amr, ONLY: BetaProlong
  use BATL_lib, ONLY: MaxLevel
  use CRASH_ModInterfaceNLTE, ONLY: UseNLTE

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.2
  character (len=*), parameter :: &
       NameUserModule = 'HYDRO + IONIZATION EQUILIBRIUM + LEVEL SETS'

  ! There are at most 5 materials: Xe, Be, Plastic, Gold, Acrylic
  integer, parameter :: MaxMaterial = 5

  ! Do we use the plastic, gold and acrylic levels ?
  ! Definitions in ModEos: Xe_=0, Be_=1, Plastic_=2, Au_=3, Ay_=4
  logical, parameter :: UsePl = nMaterial > Plastic_
  logical, parameter :: UseAu = nMaterial > Au_
  logical, parameter :: UseAy = nMaterial > Ay_

  ! Material level values in the order of Xe, Be, Plastic, Gold, Acrylic
  ! The levels that are not in use default to MaterialLast_
  integer, parameter :: LevelXe_ = MaterialFirst_
  integer, parameter :: LevelBe_ = min(MaterialFirst_+Be_,MaterialLast_)
  integer, parameter :: LevelPl_ = min(MaterialFirst_+Plastic_,MaterialLast_)
  integer, parameter :: LevelAu_ = min(MaterialFirst_+Au_,MaterialLast_)
  integer, parameter :: LevelAy_ = min(MaterialFirst_+Ay_,MaterialLast_)

  ! The Maximum Level set index that is used
  integer, parameter :: LevelMax = MaterialLast_

  ! Nozzle that shrinks the circular cross section of the tube 
  ! into an ellipse or smaller circle
  logical :: UseNozzle = .false.
  real    :: xStartNozzle = 0.0
  real    :: xEndNozzle   = 0.0
  real    :: yRatioNozzle = 0.0
  real    :: zRatioNozzle = 0.0

  ! Wall parameters
  logical:: UseTube = .false.
  real :: xEndTube   =   40.0    ! x coordinate of tube ending
  real :: rInnerTube =  287.5    ! inner radius [micron]
  real :: rOuterTube =  312.5    ! outer radius [micron]
  real :: RhoDimTube = 1430.0    ! density      [kg/m3]
  real :: RhoDimOutside = 6.5    ! density  of Xe outside tube [kg/m3]
  real :: pDimOutside   = 1.1e5  ! pressure of Xe outside tube [Pa]
  logical:: UseTubeRadius=.false.! true if tube radius is explicitly set

  !Gas parameters:

  real :: RhoDimInside = 6.5 ! density  of Xe inside tube [kg/m3]


  ! Allow overwriting the Xe state inside the tube for x > xUniformXe > 0
  real :: xUniformXe = -1.0

  ! Use conservative or non-conservative levelsets
  logical :: UseNonConsLevelSet = .true.

  ! Use volume fraction method at the material interface
  logical :: UseVolumeFraction = .false.

  ! Treat cells near material interface as a mixture
  logical :: UseMixedCell = .false.

  ! Mixed material cell is assumed if the ratio of dominant to total
  ! atomic concentration is below MixLimit
  real :: MixLimit = 0.97

  ! Variables for triangulation of 2D initialization file
  integer, allocatable:: iNodeTriangle_II(:,:)
  integer, allocatable:: nTriangle_C(:,:)
  integer, pointer    :: iTriangle_IC(:,:,:)

  ! Variables for reading 2D CRASH file
  logical           :: UseCrash2DFile  = .false. ! read CRAS2D file?
  character(len=100):: NameCrash2DFile           ! name of CRASH2D file
  character(len=5)  :: TypeCrash2DFile           ! type of CRASH2D file
  real              :: TimeCrash2D    = -1.0     ! time of CRASH2D snapshot

  integer           :: nCellRead      = -1       ! number of cells in 2D file
  real, allocatable :: CoordRead_DC(:,:)         ! coords of input file
  real, allocatable :: StateRead_VC(:,:)         ! state  of input file

  ! Variables for Hyades file
  logical           :: UseDelaunay     = .false. ! use Delaunay triangulation?
  logical           :: UseHyadesFile   = .false. ! read Hyades file?
  character(len=100):: NameHyadesFile            ! name of hyades file
  integer           :: nDimHyades      = -1      ! number of dimensions 
  real              :: xBeHyades       = -1.0    ! position of Be-Xe interface

  ! mesh coordinates: ux and ur are always defined on these coordinates
  integer           :: nCellHyadesMesh      = -1 ! number of mesh points
  integer           :: nCellHyadesMesh_D(3) = -1 ! no. mesh cells per dimension
  integer           :: nVarHyadesMesh       = -1 ! number of variables
  real, allocatable :: CoordHyadesMesh_DC(:,:)   ! mesh centered Hyades coords
  real, allocatable :: DataHyadesMesh_VC(:,:)    ! mesh centered Hyades data
  integer           :: iXHyadesMesh         = -1 ! index of x coordinate
  integer           :: iRHyadesMesh         = -1 ! index of r coordinate
  integer           :: iUxHyades            = -1 ! index of x velocity
  integer           :: iUrHyades            = -1 ! index of r velocity
  integer           :: iCellLastHyadesMesh  = -1 ! cell with maximum X and r=0

  ! center of mass coordinates of the zone:
  ! if these coordinates are present all variables except velocities are on
  ! these coordinates, otherwise they are interpolated to the mesh points
  integer           :: nCellHyades          = -1 ! number of mesh points
  integer           :: nCellHyades_D(3)     = -1 ! no. mesh cells per dimension
  integer           :: nVarHyades           = -1 ! number of variables
  real, allocatable :: CoordHyades_DC(:,:)       ! cell centered Hyades coords
  real, allocatable :: DataHyades_VC(:,:)        ! cell centered Hyades data
  integer           :: iXHyades             = -1 ! index of x coordinate
  integer           :: iRHyades             = -1 ! index of r coordinate
  integer           :: iRhoHyades           = -1 ! index of density
  integer           :: iPHyades             = -1 ! index of pressure
  integer           :: iZHyades             = -1 ! index of ionization level
  integer           :: iTeHyades            = -1 ! index of electron temper.
  integer           :: iTiHyades            = -1 ! index of ion temperature
  integer           :: iTrHyades            = -1 ! index of rad. temperature
  integer           :: iMaterialHyades      = -1 ! index of material type
  real, allocatable :: LevelHyades_VC(:,:)       ! level set functions
  integer           :: iCellLastHyades      = -1 ! cell with maximum X and r=0

  ! if we detect xcm as a HYADES variable, all variables except for the
  ! velocities are defined on the zone centers
  ! Default is that these variables are interpolated to the mesh points
  logical :: UseZoneCenter = .false.

  ! Variables for Hyades multi-group file
  logical           :: UseHyadesGroupFile = .false.! read Hyades multi-group ?
  character(len=100):: NameHyadesGroupFile         ! name of multi-group file
  real, allocatable :: EradHyades_VC(:,:)          ! Hyades group energies

  ! Opacity scale factor for sensitivity studies on opacities (UQ only !)
  real :: RosselandScaleFactor_I(0:MaxMaterial-1) = 1.0
  real :: PlanckScaleFactor_I(0:MaxMaterial-1) = 1.0

  ! Gamma law per material (UQ only !)
  logical :: UseGammaLaw = .false.
  real :: Gamma_I(0:MaxMaterial-1) = 5.0/3.0

  ! Factor to suppress the A.S.S.
  real :: AssFactor = -1.0

  ! Indexes for lookup tables
  integer:: iTablePPerRho = -1, iTableThermo = -1
  ! If UseElectronPressure is true (Pe_>1) then nThermo=6, otherwise it is 5
  integer, parameter:: nThermo=4+min(2,Pe_)
  ! Named indices for thermo lookup table
  ! If UseElectronPressure is false, default TeTi_ to Cond_ to keep
  ! the compiler happy
  integer, parameter:: Te_=1, Cv_=2, Gamma_=3, Zavg_=4, Cond_=5, TeTi_=nThermo

  integer:: iTableOpacity = -1
  integer:: iTableOpacity_I(0:MaxMaterial-1) = -1

  ! For comparison with SESAME
  integer:: iTableSesame = -1

  ! Variables for some tests
  logical :: UseWave    = .false.
  real    :: xStartWave = -100.0
  real    :: xEndWave   = +100.0
  real    :: DpWave     =  100.0

  ! electron and ion temperatures read from Hyades input
  real :: Te_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK), Ti_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

  ! Temperature limit for cold plastic (30,000K is a good limit)
  real:: TeMaxColdPlSi = -1.0

  ! Assume equal electron/ion temperature in plastic and gold at Hyades handoff
  logical :: UseEqualTemperatureHyades = .false.

  ! AMR parameters
  real:: RhoMinAmrDim = 20.0   ! kg/m3
  real:: xMaxAmr      = 2500.0 ! microns

  integer:: nLevelShock = MaxLevel
  integer:: nLevelAuInterface = MaxLevel
  integer:: nLevelBeryllium = MaxLevel

  ! Option to start with reduced radiation in the initial condition.
  real :: RadiationScaleFactor = 1.0

  ! Maximum allowed value of the inverse magnetic Reynolds number
  real :: MaxNormalizedResistivity = 1e5

  !Non LTE array to store E/B ratios
  real, allocatable :: EOverB_VGB(:,:,:,:,:) 

  ! Logical for using ModInitialState
  logical:: UseInitialStateDefinition = .false.

contains

  !============================================================================
  subroutine user_read_inputs

    use ModReadParam
    use CRASH_ModEos,        ONLY: read_eos_parameters, NameMaterial_I
    use CRASH_ModMultiGroup, ONLY: read_opacity_parameters
    use ModAdvance,          ONLY: UseElectronPressure
    use ModWaves,            ONLY: FreqMinSI, FreqMaxSI
    use ModConst,            ONLY: cHPlanckEV
    use CRASH_ModInterfaceNLTE, ONLY: read_nlte
    use ModVarIndexes,       ONLY: nVar, NameVar_V
    use ModInitialState,     ONLY: init_initial_state, read_initial_state_param
    use ModUtilities,        ONLY: join_string

    integer :: iMaterial
    real :: EnergyPhotonMin, EnergyPhotonMax
    character (len=100) :: NameCommand
    character (len=500) :: StringVar

    character(len=*), parameter :: NameSub = 'user_read_inputs'
    !------------------------------------------------------------------------

    UseUserLogFiles     = .true. ! to allow integration of rhoxe, rhobe...
    UseUserUpdateStates = .true. ! for internal energy
    UseUserSource       = .true. ! for non-cons levelset and pressure equation
    UseUserInitSession  = .true. ! to set units for level set variables
    UseUserIcs          = .true. ! to read in Hyades file
    !                              and initialize the level set variables

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#STATEDEFINITION")
          UseInitialStateDefinition = .true.
          call join_string(nVar, NameVar_V(1:nVar), StringVar)
          call init_initial_state(StringVar, MaterialFirst_, nMaterial)
          call read_initial_state_param(NameCommand)
       case("#STATEINTERFACE")
          call read_initial_state_param(NameCommand)
       case("#CRASH2D")
          call read_var('UseCrash2DFile',  UseCrash2DFile)
          call read_var('NameCrash2DFile', NameCrash2DFile)
          call read_var('TypeCrash2DFile', TypeCrash2DFile)
          call read_var('TimeCrash2D'    , TimeCrash2D)

       case("#HYADES")
          call read_var('UseHyadesFile', UseHyadesFile)
          call read_var('NameHyadesFile',NameHyadesFile)

       case("#HYADESGROUP")
          call read_var('UseHyadesGroupFile', UseHyadesGroupFile)
          call read_var('NameHyadesGroupFile',NameHyadesGroupFile)

       case("#OPACITY")
          call read_opacity_parameters

       case("#GROUPRANGE")
          call read_var('EnergyPhotonMin', EnergyPhotonMin)  ! in eV
          call read_var('EnergyPhotonMax', EnergyPhotonMax)  ! in eV
          FreqMinSi = EnergyPhotonMin/cHPlanckEV
          FreqMaxSi = EnergyPhotonMax/cHPlanckEV

       case("#TUBE")
          UseTube = .true.
          call read_var('xEndTube',   xEndTube)
          call read_var('rInnerTube', rInnerTube)
          call read_var('rOuterTube', rOuterTube)
          call read_var('RhoDimTube', RhoDimTube)
          call read_var('RhoDimOutside', RhoDimOutside)
          call read_var('pDimOutside',   pDimOutside)
          call read_var('xUniformXe',    xUniformXe)
          UseTubeRadius = .true.

       case("#TUBEDIAMETER")
          call read_var('rInnerTube', rInnerTube)
          call read_var('rOuterTube', rOuterTube)
          UseTubeRadius = .true.

       case('#RHOINSIDE')
          call read_var('RhoDimInside',RhoDimInside)

       case("#LEVELSET")
          call read_var('UseNonConsLevelSet', UseNonConsLevelSet)

       case("#VOLUMEFRACTION")
          call read_var('UseVolumeFraction', UseVolumeFraction)

       case("#MIXEDCELL")
          call read_var('UseMixedCell', UseMixedCell)
          if(UseMixedCell)then
             call read_var('MixLimit', MixLimit)

             ! The mixed cell implementation is limited to 3 materials only
             if(nMaterial>3) call stop_mpi(NameSub // " Gold and Acrylic "// &
                  "are not yet supported in the mixed cell approach")
             if(nMaterial<3) call stop_mpi(NameSub // " the mixed cell "// &
                  "implementation needs 3 materials, including plastic")
          end if

       case("#EOS", "#MATERIAL", "#EOSTABLE", "#OPACITYTABLE", &
            "#USEEOSTABLE", "#USEOPACTABLE")
          call read_eos_parameters(nMaterial, NameCommand)

       case("#OPACITYSCALEFACTOR") ! UQ only
          do iMaterial = 0, nMaterial - 1
             call read_var('PlanckScaleFactor'//NameMaterial_I(iMaterial), &
                  PlanckScaleFactor_I(iMaterial))
          end do
          do iMaterial = 0, nMaterial - 1
             call read_var('RosselandScaleFactor'//NameMaterial_I(iMaterial), &
                  RosselandScaleFactor_I(iMaterial))
          end do

       case("#GAMMALAW") ! UQ only
          call read_var('UseGammaLaw', UseGammaLaw)
          if(UseGammaLaw)then
             if(.not.UseElectronPressure) call stop_mpi( &
                  NameSub//" can not use gamma-law without electron pressure")
             if(UseNLTE) call stop_mpi( &
                  NameSub//" can not use gamma-law with NLTE EOS")

             do iMaterial = 0, nMaterial - 1
                call read_var('Gamma'//NameMaterial_I(iMaterial), &
                     Gamma_I(iMaterial))
             end do
          end if

       case("#THREEDIM")
          ! Not needed anymore
       case("#NOZZLE")
          call read_var('UseNozzle', UseNozzle)
          if(UseNozzle)then
             call read_var('xStartNozzle', xStartNozzle)
             call read_var('xEndNozzle',   xEndNozzle)
             call read_var('yRatioNozzle', yRatioNozzle)
             call read_var('zRatioNozzle', zRatioNozzle)
          end if
       case("#WAVE")
          call read_var('UseWave',    UseWave)
          if(UseWave)then
             call read_var('xStartWave', xStartWave)
             call read_var('xEndWave',   xEndWave)
             call read_var('DpWave',     DpWave)
          end if

       case("#USERAMR")
          call read_var('RhoMinAmr',   RhoMinAmrDim)
          call read_var('xMaxAmr',     xMaxAmr)
          call read_var('BetaProlong', BetaProlong)

       case("#AMRLEVEL")
          call read_var('nLevelShock', nLevelShock)
          call read_var('nLevelAuInterface', nLevelAuInterface)
          call read_var('nLevelBeryllium', nLevelBeryllium)

       case("#PLASTIC")
          call read_var('TeMaxColdPlSi',  TeMaxColdPlSi)

       case("#TEMPERATUREHYADES")
          call read_var('UseEqualTemperatureHyades', UseEqualTemperatureHyades)

       case("#RADIATIONSCALEFACTOR")
          ! Multiply the initial radiation by RadiationScaleFactor
          call read_var('RadiationScaleFactor', RadiationScaleFactor)

       case("#ASS")
          ! Multiply plastic pressure and negative Uy/Ur by AssFactor
          call read_var('AssFactor', AssFactor)

       case("#MAXRESISTIVITY")
          ! Maximum allowed value of the inverse magnetic Reynolds number
          call read_var('MaxNormalizedResistivity', MaxNormalizedResistivity)

       case("#NLTE")
          call read_nlte

       case('#USERINPUTEND')
          EXIT
       case default
          call stop_mpi('ERROR in ModUserCrash: unknown command='//NameCommand)
       end select

    end do

    if(UseMixedCell) UseNonConsLevelSet = .false.

  end subroutine user_read_inputs

  !============================================================================
  subroutine user_set_ics(iBlock)

    use ModMain,        ONLY: nI, nJ, nK, UseRadDiffusion
    use ModPhysics,     ONLY: inv_gm1, ShockPosition, ShockSlope, &
         Io2No_V, No2Si_V, Si2No_V, UnitRho_, UnitP_, UnitEnergyDens_, &
         UnitTemperature_, UnitN_, PeMin, ExtraEintMin
    use ModAdvance,     ONLY: State_VGB, UseElectronPressure
    use ModVarIndexes,  ONLY: Rho_, RhoUx_, RhoUz_, p_, ExtraEint_, &
         Pe_, Erad_, WaveFirst_, WaveLast_, Te0_, nWave
    use ModGeometry,    ONLY: Xyz_DGB
    use ModConst,       ONLY: cPi, cAtomicMass
    use CRASH_ModEos,   ONLY: eos, Xe_, Plastic_
    use CRASH_ModMultiGroup, ONLY:get_planck_g_from_temperature
    use ModInitialState,ONLY: get_initial_state

    integer, intent(in) :: iBlock

    real    :: x, y, z, r, xBe, DxBe, DxyPl, EinternalSi, Te0SI, PlanckSI
    real    :: TeSi, PeSi, Natomic, NatomicSi, RhoSi, pSi, p, Te

    integer :: i, j, k, iMaterial, iP, iGroup

    character(len=*), parameter :: NameSub = "user_set_ics"
    !------------------------------------------------------------------------

    if(UseCrash2dFile)then

       call timing_start('interpolate_crash2d')
       call interpolate_crash2d(iBlock)
       call timing_stop('interpolate_crash2d')

       RETURN
    end if

    call timing_start(NameSub)

    if(UseHyadesFile)then

       call timing_start('interpolate_hyades')
       ! interpolate Hyades output
       if(nDimHyades == 1)then
          call interpolate_hyades1d(iBlock)
       else
          call interpolate_hyades2d(iBlock)
       end if
       call timing_stop('interpolate_hyades')

    end if

    if(UseElectronPressure .and. UseTube)then
       call stop_mpi(NameSub //" electron energy does not yet work " &
            //"with plastic tube or gold washer")
    end if

    ! Set level set functions, internal energy, and other values
    do k=1, nK; do j=1, nJ; do i=1, nI 

       x = Xyz_DGB(x_,i,j,k,iBlock)
       y = Xyz_DGB(y_,i,j,k,iBlock)
       z = Xyz_DGB(z_,i,j,k,iBlock)

       call set_yzr(x, y, z, r)

       if(UseTube)then
          ! Distance from plastic wall: 
          ! positive for rInnerTube < |y| < rOuterTube and x > xEndTube only
          DxyPl = &
               min(r - rInnerTube, rOuterTube - r, x - xEndTube)

          ! Set plastic tube state
          if(DxyPl > 0.0)then

             ! Use the density and pressure given by the #TUBE command
             State_VGB(Rho_,i,j,k,iBlock) = RhoDimTube*Io2No_V(UnitRho_)
             State_VGB(p_  ,i,j,k,iBlock) = pDimOutside*Io2No_V(UnitP_)
             State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0

             call set_small_radiation_energy
             if(nDimHyades == 2)then
                State_VGB(LevelPl_,i,j,k,iBlock) =  DxyPl
                State_VGB(LevelXe_,i,j,k,iBlock) =  &
                     max(State_VGB(LevelBe_,i,j,k,iBlock), r - rOuterTube)
             end if
          end if

          ! Set pressure and speed outside the tube. 
          ! For 1D Hyades input do not overwrite values left of xEndTube
          if(r > rOuterTube .and. (nDimHyades == 1 .or. x > xEndTube) ) then
             State_VGB(Rho_,i,j,k,iBlock) = RhoDimOutside*Io2No_V(UnitRho_)
             State_VGB(p_  ,i,j,k,iBlock) = pDimOutside*Io2No_V(UnitP_)
             State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
             call set_small_radiation_energy
             if(nDimHyades == 2)then
                State_VGB(LevelXe_,i,j,k,iBlock) =  &
                     min(r - rOuterTube,x - xEndTube) 
                State_VGB(LevelPl_,i,j,k,iBlock) =  rOuterTube - r
             end if
          end if

          ! Set the Xe state inside the tube for x > xUniformXe if it is set
          if(xUniformXe > 0.0 .and. x > xUniformXe .and. r < rInnerTube)then
             State_VGB(Rho_,i,j,k,iBlock) = RhoDimInside*Io2No_V(UnitRho_)
             State_VGB(p_  ,i,j,k,iBlock) = pDimOutside*Io2No_V(UnitP_)
             State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
             call set_small_radiation_energy
          end if
       end if ! UseTube

       ! Create sound wave by making a pressure hump (for testing)
       if(UseWave .and. x > xStartWave .and. x < xEndWave) &
            State_VGB(p_,i,j,k,iBlock) = State_VGB(p_,i,j,k,iBlock) + &
            Io2No_V(UnitP_)*DpWave &
            *sin( cPi*(x - xStartWave)/(xEndWave - xStartWave) )**2

       if(nDimHyades /= 2)then

          if(UseHyadesFile)then
             ! Be - Xe interface is given by Hyades file
             xBe = xBeHyades
          else
             ! Be - Xe interface is at the shock defined by #SHOCKPOSITION
             xBe = ShockPosition - ShockSlope*Xyz_DGB(y_,i,j,k,iBlock)
          end if

          ! Distance from Be disk: positive for x < xBe
          DxBe = xBe - x

          ! Add a plastic tube if required
          if(UseTube)then
             ! Distance from plastic wall: 
             ! positive for rInnerTube < |y| < rOuterTube and x > xEndTube only
             DxyPl = min(r - rInnerTube, rOuterTube - r, x - xEndTube)

             ! Berylium is left of xBe inside rInnerTube 
             ! and it is left of xEndTube outside
             State_VGB(LevelBe_,i,j,k,iBlock) = &
                  max(xEndTube - x, min(DxBe, rInnerTube - r))

             ! Xenon is right of xBe inside rInnerTube and 
             ! right of xEndTube outside rOuterTube
             State_VGB(LevelXe_,i,j,k,iBlock) = max( &
                  min( x - xEndTube, r - rOuterTube), &
                  min( -DxBe, rInnerTube - r) )

             ! Plastic 
             State_VGB(LevelPl_,i,j,k,iBlock) = DxyPl
          else
             ! If there is no plastic tube, things are easy
             State_VGB(LevelBe_,i,j,k,iBlock) =  DxBe
             State_VGB(LevelXe_,i,j,k,iBlock) = -DxBe
             do iMaterial = Plastic_, nMaterial - 1
                State_VGB(LevelXe_+iMaterial,i,j,k,iBlock) = -1e30
             end do
          end if

       end if ! nDimHyades /= 2
       !When the hyades file is used for Xe, and the actual material is
       !different, rescale the density of xenon:
       !if(UseHyadesFile.and.&
       !maxloc(State_VGB(LevelXe_:LevelMax,i,j,k,iBlock), 1)==1)&
       !State_VGB(Rho_,i,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock)*&
       !RhoDimInside/RhoHe
       if(UseMixedCell)then
          if(nMaterial>3) call stop_mpi(NameSub // " Gold and Acrylic "// &
               "are not yet supported in the mixed cell approach")
          if(nMaterial<3) call stop_mpi(NameSub //"mixed cell needs plastic")

          ! Use atomic concentrations instead of smooth level set functions

          ! Used materials: Xe, Be, Pl, and optionally Au
          if(maxval( State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock) ) <= 0.0)then
             ! Ignore Au
             State_VGB(LevelXe_:LevelMax,i,j,k,iBlock) = &
                  1.0/(3*MassMaterial_I(Xe_:nMaterial-1))
          else
             State_VGB(LevelXe_:LevelMax,i,j,k,iBlock) = &
                  max(0.0, State_VGB(LevelXe_:LevelMax,i,j,k,iBlock))

             where( State_VGB(LevelXe_:LevelMax,i,j,k,iBlock) > 0.0) &
                  State_VGB(LevelXe_:LevelMax,i,j,k,iBlock) = &
                  1./MassMaterial_I(0:nMaterial-1)
          end if

       end if

       if(UseInitialStateDefinition)then
          call get_initial_state( (/x, y/), State_VGB(:,i,j,k,iBlock) )
          call set_initial_temperature
       end if

       ! Multiply level set functions with density unless the 
       ! non-conservative approach is used
       if(.not.UseNonConsLevelSet) &
            State_VGB(LevelXe_:LevelMax,i,j,k,iBlock) = &
            State_VGB(LevelXe_:LevelMax,i,j,k,iBlock) &
            *State_VGB(Rho_,i,j,k,iBlock)

    end do; end do; end do

    iP = p_
    if(UseElectronPressure) iP = Pe_

    ! Set the remaining State_VGB quantities that involve
    ! user_material_properties
    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(UseElectronPressure)then

          TeSi = Te_G(i,j,k)*No2Si_V(UnitTemperature_)
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, TeIn=TeSi, &
               PressureOut=PeSi, NatomicOut=NatomicSi)

          State_VGB(Pe_,i,j,k,iBlock) = max(PeMin, PeSi*Si2No_V(UnitP_))

          if(UsePl .and. UseEqualTemperatureHyades .and. .not.UseMixedCell &
               .and. any(State_VGB(LevelPl_:LevelMax,i,j,k,iBlock) > 0.0))then
             ! equal electron/ion temperature for plastic and gold
             iMaterial = maxloc(State_VGB(LevelXe_:LevelMax,i,j,k,iBlock), 1)-1
             RhoSi = State_VGB(Rho_,i,j,k,iBlock)*No2Si_V(UnitRho_)
             pSi   = State_VGB(p_,i,j,k,iBlock)*No2Si_V(UnitP_)
             call eos(iMaterial, RhoSi, PtotalIn=pSi, TeOut=TeSi)
             Te = TeSi*Si2No_V(UnitTemperature_)
             Natomic = RhoSi/(cAtomicMass*MassMaterial_I(iMaterial)) &
                  *Si2No_V(UnitN_)
             p = Natomic*Te
             State_VGB(Pe_,i,j,k,iBlock) = &
                  max(PeMin,State_VGB(p_,i,j,k,iBlock) - p)
             State_VGB(p_,i,j,k,iBlock) = p
          else if(UsePl .and. State_VGB(LevelPl_,i,j,k,iBlock) > 0.0 &
               .and. TeSi < TeMaxColdPlSi)then
             ! Subtract electron pressure from the total pressure
             State_VGB(p_,i,j,k,iBlock)  = max(PeMin, &
                  State_VGB(p_,i,j,k,iBlock) - State_VGB(Pe_,i,j,k,iBlock))
          else
             Natomic = NatomicSi*Si2No_V(UnitN_)
             State_VGB(p_,i,j,k,iBlock)  = Natomic*Ti_G(i,j,k)
          end if
       end if

       ! Calculate internal energy
       call user_material_properties(State_VGB(:,i,j,k,iBlock), &
            i, j, k, iBlock, EinternalOut=EinternalSi)

       State_VGB(ExtraEint_,i,j,k,iBlock) = max(ExtraEintMin, &
            EinternalSi*Si2No_V(UnitEnergyDens_) &
            - inv_gm1*State_VGB(iP,i,j,k,iBlock))

       if(UseRadDiffusion) &
            State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock) = &
            State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock)*RadiationScaleFactor
       if(UseNlte)then
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i,j,k,iBlock, TeOut=Te0SI)
          State_VGB(Te0_,i,j,k,iBlock) = Te0SI * Si2No_V(UnitTemperature_)
       
          !Calculate "E_rad"
          EOverB_VGB(:,i,j,k,iBlock) = No2Si_V(UnitEnergyDens_)* &
               State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock)
             
          !Calculate "over B"
          Te0Si = State_VGB(Te0_,i,j,k,iBlock) * No2Si_V(UnitTemperature_)
          do iGroup = 1, nWave
             call get_planck_g_from_temperature(iGroup, Te0Si, EgSI=PlanckSi)
             if(EOverB_VGB(iGroup,i,j,k,iBlock)>0.0)then
                EOverB_VGB(iGroup,i,j,k,iBlock) = EOverB_VGB(iGroup,i,j,k,iBlock)/&
                     max(PlanckSi, 0.01*EOverB_VGB(iGroup,i,j,k,iBlock))
             else
                EOverB_VGB(iGroup,i,j,k,iBlock)=0.0
             end if
          end do
          

          EInternalSI = (State_VGB(iP,i,j,k,iBlock)*inv_gm1+&
                 State_VGB(ExtraEInt_,i,j,k,iBlock))*No2Si_V(UnitEnergyDens_)
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, EInternalIn=EInternalSI, &
               PressureOut=PeSi)

          State_VGB(iP,i,j,k,iBlock) = max(PeMin, PeSi*Si2No_V(UnitP_))
    

       ! Calculate internal energy

       State_VGB(ExtraEint_,i,j,k,iBlock) = max(ExtraEintMin, &
            EinternalSi*Si2No_V(UnitEnergyDens_) &
            - inv_gm1*State_VGB(iP,i,j,k,iBlock))
       end if

    end do; end do; end do

    call timing_stop(NameSub)
  contains
    !==========================================================================
    subroutine set_small_radiation_energy

      use ModMain,ONLY: UseRadDiffusion
      use ModPhysics,ONLY: cRadiationNo, Si2No_V, UnitTemperature_
      !----------------------------------------------------------------------
      if(.not.UseRadDiffusion)RETURN

      State_VGB(Erad_,i,j,k,iBlock) = cRadiationNo * &
           (500.0 * Si2No_V(UnitTemperature_))**4

    end subroutine set_small_radiation_energy
    !========================================================================
    subroutine set_initial_temperature

      use CRASH_ModMultiGroup, ONLY: get_energy_g_from_temperature
      use ModVarIndexes,  ONLY: nWave, WaveFirst_
      use ModPhysics,     ONLY: cRadiationNo

      real    :: Tr = 0.0259                      ! [eV?]
      real    :: TrSi, EgSi
      integer :: iWave
      !---------------------------------------------------------------------
      TeSi = 2.1e6 *(0.025/180.) 
      Te_G(i,j,k) = TeSi*Si2No_V(UnitTemperature_)
      Ti_G(i,j,k) = TeSi*Si2No_V(UnitTemperature_)

      if (nWave == 1) then
         ! Set small initial radiation-energy density everywhere.
         State_VGB(Erad_,i,j,k,iBlock) = cRadiationNo*Tr**4
      else
         TrSi = 11600.

         do iWave = 1, nWave
            call get_energy_g_from_temperature(iWave, TrSi,EgSI=EgSi)
            State_VGB(WaveFirst_+iWave-1,i,j,k,iBlock) = &
                 EgSi*Si2No_V(UnitEnergyDens_)
         enddo
      endif

    end subroutine set_initial_temperature

  end subroutine user_set_ics


  !============================================================================

  subroutine read_crash2d_file

    use ModProcMH,   ONLY: iProc
    use ModAdvance,  ONLY: nVar
    use ModMain,     ONLY: Time_Simulation
    use ModPlotFile, ONLY: read_plot_file
    use ModIoUnit,   ONLY: UnitTmp_

    integer:: nDimRead, nVarRead, nCellRead_D(2), iSnapshot
    character(len=500):: NameVarRead

    character(len=*), parameter :: NameSub = "ModUser::read_crash2d_file"
    !-------------------------------------------------------------------------
    do iSnapShot = 1, 100000
       
       call read_plot_file(NameCrash2DFile, UnitTmp_, TypeCrash2DFile, &
            TimeOut = Time_Simulation, nDimOut = nDimRead, nVarOut = nVarRead,&
            nOut_D = nCellRead_D, NameVarOut = NameVarRead)

       if(nDimRead /= 2 .and. iProc==0)then
          write(*,*)'ERROR: nDimRead=', nDimRead,' instead of 2'
          call stop_mpi(NameSub//': cannot use file '//trim(NameCrash2DFile))
       end if

       if(nVarRead /= nVar)then
          write(*,*)'ERROR: nVarRead=',nVarRead,' is different from nVar=',nVar
          call stop_mpi(NameSub// &
               ': incorrect number of variables in '//trim(NameCrash2DFile))
       end if

       ! total number of mesh cells
       nCellRead = product(nCellRead_D)

       if(allocated(CoordRead_DC)) deallocate(CoordRead_DC, StateRead_VC)
       allocate(CoordRead_DC(2,nCellRead), StateRead_VC(nVar,nCellRead))

       ! Read in the data
       call read_plot_file(NameCrash2DFile, UnitTmp_, TypeCrash2DFile, &
            CoordOut_DI = CoordRead_DC, VarOut_VI = StateRead_VC)

       if(time_simulation >= TimeCrash2D) EXIT

    end do
    close(UnitTmp_)

    if(iProc==0)then
       write(*,*) NameSub, ': read 2D file     = ', trim(NameCrash2DFile)
       if(iSnapshot > 1) &
            write(*,*) NameSub, ': snapshot number  =',iSnapshot
       write(*,*) NameSub, ': simulation time  =', time_simulation
       write(*,*) NameSub, ': number of cells  =', nCellRead
       write(*,*) NameSub, ': number of vars   =', nVarRead
       write(*,*) NameSub, ': coord array size =', size(CoordRead_DC)
       write(*,*) NameSub, ': state array size =', size(StateRead_VC)
       write(*,*) NameSub, ': coord/var/par names = ', trim(NameVarRead)
    end if

  end subroutine read_crash2d_file

  !============================================================================

  subroutine interpolate_crash2d(iBlock)

    ! Use Delaunay triangulation to interpolate 2D CRASH grid onto 
    ! a different (3D?) CRASH grid

    use ModProcMH,      ONLY: iProc
    use ModSize,        ONLY: nI, nJ, nK
    use ModAdvance,     ONLY: State_VGB, RhoUy_, RhoUz_, By_, Bz_, UseB
    use ModGeometry,    ONLY: Xyz_DGB
    use ModTriangulate, ONLY: calc_triangulation, find_triangle

    integer, intent(in) :: iBlock

    integer, save              :: nTriangle
    real, save :: xMin, xMax, rMax

    integer :: i, j, k, iNode1, iNode2, iNode3
    real    :: x, y, z, r, Weight1, Weight2, Weight3
    real    :: RhoUr, Bphi

    character(len=*), parameter :: NameSub='interpolate_crash2d'
    !-------------------------------------------------------------------------
    if(.not.allocated(iNodeTriangle_II))then

       if(.not.allocated(CoordRead_DC) .and. iProc==0) &
            call CON_stop(NameSub// &
            ': no 2D CRASH file was read. Add #CRASH2D command')

       ! allocate variables for triangulation, nTriangle_C for fast search.
       allocate(iNodeTriangle_II(3,2*nCellRead), nTriangle_C(200,200))
       nTriangle_C(1,1) = -1

       ! Triangulate the 2D CRASH mesh
       call calc_triangulation( &
            nCellRead, CoordRead_DC, iNodeTriangle_II, nTriangle)

       ! Determine (approximate) coordinate limits
       xMin = minval(CoordRead_DC(1,:))
       xMax = maxval(CoordRead_DC(1,:))
       rMax = maxval(CoordRead_DC(2,:))
    end if

    ! Interpolate points 
    do j = 1, nJ; do i = 1, nI; do k = 1, nk

       x = Xyz_DGB(x_,i,j,k,iBlock)
       y = Xyz_DGB(y_,i,j,k,iBlock)
       z = Xyz_DGB(z_,i,j,k,iBlock)

       call set_yzr(x, y, z, r, rMax)

       ! Stay within the size of the 2D file
       x = max(xMin, min(xMax,x))

       ! Find the triangle around this position
       call find_triangle(&
            nCellRead, nTriangle, &
            (/x, r/), CoordRead_DC, &
            iNodeTriangle_II, &
            iNode1, iNode2, iNode3, Weight1, Weight2, Weight3, &
            nTriangle_C=nTriangle_C, iTriangle_IC=iTriangle_IC)

       ! Interpolate state (the CRASH input is already smooth)
       State_VGB(:,i,j,k,iBlock) = &
            Weight1*StateRead_VC(:,iNode1) + &
            Weight2*StateRead_VC(:,iNode2) + &
            Weight3*StateRead_VC(:,iNode3)

       ! Done if running in 2D
       if(nK == 1) CYCLE

       ! Correct the y,z velocity components to point in radial direction
       RhoUr = State_VGB(RhoUy_,i,j,k,iBlock)
       State_VGB(RhoUy_:RhoUz_,i,j,k,iBlock) = (/y, z/)/r * RhoUr

       ! Done with hydro
       if(.not. UseB) CYCLE

       ! Correct the magnetic field to point in the azimuthal direction
       Bphi = State_VGB(Bz_,i,j,k,iBlock)
       State_VGB(By_:Bz_,i,j,k,iBlock) = (/-z, y/)/r * Bphi

    end do; end do; end do

  end subroutine interpolate_crash2d

  !============================================================================

  subroutine user_action(NameAction)

    use ModProcMH, ONLY: iProc
    use ModInitialState, ONLY: clean_initial_state

    character(len=*), intent(in):: NameAction

    character(len=*), parameter:: NameSub = 'user_action'
    !-------------------------------------------------------------------------
    if(iProc==0)write(*,*) NameSub,' called with action ',NameAction

    select case(NameAction)
    case('initial condition done')
       if(allocated(CoordRead_DC)) deallocate( &
            CoordRead_DC, StateRead_VC)
       if(allocated(iNodeTriangle_II)) deallocate( &
            iNodeTriangle_II, nTriangle_C, iTriangle_IC)
       if(allocated(DataHyadesMesh_VC)) deallocate( &
            DataHyadesMesh_VC, DataHyades_VC, &
            CoordHyades_DC, CoordHyadesMesh_DC)

       if(UseInitialStateDefinition)then
          call clean_initial_state
          UseInitialStateDefinition = .false.
       end if
    end select

  end subroutine user_action

  !============================================================================

  subroutine read_hyades_file

    use ModAdvance,    ONLY: UseElectronPressure
    use ModPhysics,    ONLY: Si2No_V, Io2No_V, UnitX_, UnitRho_, UnitU_, &
         UnitP_, UnitTemperature_, UnitEnergyDens_
    use ModPlotFile,   ONLY: read_plot_file
    use ModUtilities,  ONLY: split_string
    use CRASH_ModEos,  ONLY: Xe_, Be_, Plastic_, Au_, Ay_
    use ModConst,      ONLY: cKevToK, cHPlanckEV
    use ModMain,       ONLY: UseRadDiffusion, Time_Simulation
    use ModVarIndexes, ONLY: nWave
    use ModWaves,      ONLY: FreqMinSi, FreqMaxSi

    integer :: nVarHyadesAll
    integer :: iVarZone, iVarMesh, iVar, iCellMesh, j
    integer, allocatable :: iVarZone_I(:), iVarMesh_I(:)

    real                :: TimeHyades
    real, allocatable   :: Hyades2No_V(:)
    character(len=100)  :: NameVarHyades

    ! Variables for variable names
    integer, parameter:: MaxString = 20
    character(len=10) :: String_I(MaxString)
    integer           :: nString

    ! Variables for reading in coordinates and variables
    real, allocatable:: Var_VI(:,:)

    ! Variables for setting level set functions
    integer :: i, iCell, iMaterial, jMaterial
    real    :: x, r
    integer, allocatable:: iMaterial_C(:)
    real,    allocatable:: Distance2_C(:)

    ! variables for HYADES multi-group
    integer, parameter :: nGroupMax = 100
    integer :: nGroupHyades, nGroup
    integer :: iGroup, iGroupFirst, iGroupLast
    real :: EnergyGroupHyades_I(0:nGroupMax) ! Photon energy (unit of keV)
    real :: EnergyGroupMin, EnergyGroupMax   ! unit of keV
    real :: DeltaLogEnergy

    character(len=*), parameter :: NameSub = "ModUser::read_hyades_file"
    !-------------------------------------------------------------------------

    nCellHyades_D = 1
    nCellHyadesMesh_D = 1
    call read_plot_file(NameHyadesFile, &
         TimeOut = TimeHyades, nDimOut = nDimHyades, nVarOut = nVarHyadesAll, &
         nOut_D = nCellHyadesMesh_D, NameVarOut = NameVarHyades)

    ! reset simulation time to HYADES time
    Time_Simulation = TimeHyades

    ! extract coordinate, variable and eqpar names
    call split_string(NameVarHyades, MaxString, String_I, nString)

    ! total number of mesh cells
    nCellHyadesMesh = product(nCellHyadesMesh_D)

    ! Do not use zone centers unless xcm is in the variable list
    UseZoneCenter = index(NameVarHyades,'xcm') > 0

    if(UseZoneCenter)then
       nCellHyades_D(1:nDimHyades) = nCellHyadesMesh_D(1:nDimHyades) - 1
       nCellHyades   = product(nCellHyades_D)
    else
       ! The zone centered HYADES variables are interpolated to the mesh center
       nCellHyades_D = nCellHyadesMesh_D
       nCellHyades   = nCellHyadesMesh
    end if

    iVarMesh = 0
    iVarZone = 0

    allocate( &
         iVarMesh_I(nDimHyades + nVarHyadesAll), &
         iVarZone_I(nDimHyades + nVarHyadesAll) )

    ! Find the columns for the coordinates and variables
    do i = 1, nDimHyades + nVarHyadesAll
       ! The first nDimHyades strings are for the coordinates
       select case(String_I(i))
       case('x', 'y', 'r', 'ux', 'uy', 'ur')
          ! Staggered mesh variables
          iVarMesh = iVarMesh + 1
          iVarMesh_I(iVarMesh) = i
          select case(String_I(i))
          case('x')
             iXHyadesMesh = iVarMesh
             ! if zone centers are not used, this means that the variables are
             ! interpolated to the mesh centers
             if(.not.UseZoneCenter)then
                iVarZone = iVarZone + 1
                iVarZone_I(iVarZone) = i
                iXHyades = iVarZone
             end if
          case('y', 'r')
             iRHyadesMesh = iVarMesh
             if(.not.UseZoneCenter)then
                iVarZone = iVarZone + 1
                iVarZone_I(iVarZone) = i
                iRHyades = iVarZone
             end if
          case('ux')
             iUxHyades  = iVarMesh
          case('uy', 'ur')
             iUrHyades  = iVarMesh
          end select

       case default
          ! zone centered mesh variables
          iVarZone = iVarZone + 1
          iVarZone_I(iVarZone) = i
          select case(String_I(i))
          case('xcm')
             iXHyades   = iVarZone
          case('rcm')
             iRHyades   = iVarZone
          case('rho')
             iRhoHyades = iVarZone
          case('p')
             iPHyades   = iVarZone
          case('te')
             iTeHyades  = iVarZone
          case('ti')
             iTiHyades  = iVarZone
          case('tr')
             iTrHyades  = iVarZone
          case('z')
             iZHyades   = iVarZone
          case('material')
             iMaterialHyades = iVarZone
          end select
       end select
    end do

    nVarHyades = iVarZone - nDimHyades
    nVarHyadesMesh = iVarMesh - nDimHyades

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
       if(iRHyadesMesh < 0) call stop_mpi(NameSub// &
            ' could not find y/r in '//trim(NameVarHyades))
       if(UseZoneCenter .and. iRHyades < 0) call stop_mpi(NameSub// &
            ' could not find rcm in '//trim(NameVarHyades))
       if(iUrHyades < 0) call stop_mpi(NameSub// &
            ' could not find uy/ur in '//trim(NameVarHyades))
       if(iMaterialHyades < 0) call stop_mpi(NameSub// &
            ' could not find material in '//trim(NameVarHyades))
    end if

    ! Set conversion from Hyades units to normalized units
    allocate(Hyades2No_V(nDimHyades + nVarHyadesAll))
    Hyades2No_V = 1.0
    Hyades2No_V(iVarZone_I(iXHyades))  = 0.01*Si2No_V(UnitX_)   !cm    -> m
    Hyades2No_V(iVarZone_I(iRhoHyades))= 1e3 *Si2No_V(UnitRho_) !g/cm3 -> kg/m3
    Hyades2No_V(iVarZone_I(iPHyades))  = 0.1 *Si2No_V(UnitP_)   !dyne  -> Pa

    Hyades2No_V(iVarMesh_I(iXHyadesMesh)) = 0.01*Si2No_V(UnitX_) ! cm    -> m
    Hyades2No_V(iVarMesh_I(iUxHyades))    = 0.01*Si2No_V(UnitU_) ! cm/s  -> m/s

    if(UseRadDiffusion .or. UseElectronPressure)then
       if(iTeHyades < 0) call stop_mpi(NameSub// &
            ' could not find electron temperature in '//trim(NameVarHyades))

       Hyades2No_V(iVarZone_I(iTeHyades)) = &
            cKevToK*Si2No_V(UnitTemperature_) ! KeV   -> K
    end if

    if(UseRadDiffusion)then
       if(iTrHyades < 0) call stop_mpi(NameSub// &
            ' could not find radiation temperature in '//trim(NameVarHyades))

       Hyades2No_V(iVarZone_I(iTrHyades)) = &
            cKevToK*Si2No_V(UnitTemperature_) ! KeV   -> K
    end if

    if(UseElectronPressure)then
       if(iTiHyades < 0) call stop_mpi(NameSub// &
            ' could not find ion temperature in '//trim(NameVarHyades))

       Hyades2No_V(iVarZone_I(iTiHyades)) = &
            cKevToK*Si2No_V(UnitTemperature_) ! KeV   -> K
    end if

    if(nDimHyades > 1)then
       Hyades2No_V(iVarZone_I(iRHyades)) = 0.01*Si2No_V(UnitX_)  ! cm    -> m

       Hyades2No_V(iVarMesh_I(iRHyadesMesh)) = &
            0.01*Si2No_V(UnitX_)  ! cm    -> m
       Hyades2No_V(iVarMesh_I(iUrHyades)) = &
            0.01*Si2No_V(UnitU_)  ! cm/s  -> m/s
    end if

    ! Read in the data
    allocate( &
         DataHyadesMesh_VC(nDimHyades + nVarHyadesMesh, nCellHyadesMesh), &
         DataHyades_VC(nDimHyades + nVarHyades, nCellHyades), &
         CoordHyades_DC(nDimHyades,nCellHyades), &
         CoordHyadesMesh_DC(nDimHyades,nCellHyadesMesh), &
         Var_VI(nVarHyadesAll,nCellHyadesMesh) )

    call read_plot_file(NameHyadesFile, &
         CoordOut_DI = CoordHyadesMesh_DC, VarOut_VI = Var_VI)

    iCell = 0
    iCellMesh = 0
    do j = 1, nCellHyadesMesh_D(2); do i = 1, nCellHyadesMesh_D(1)
       iCellMesh = iCellMesh + 1

       ! Convert from CGS to normalized units and store in DataHyadesMesh_VC
       DataHyadesMesh_VC(:nDimHyades,iCellMesh) = &
            CoordHyadesMesh_DC(:,iCellMesh)*Hyades2No_V(:nDimHyades)
       ! the staggered velocity components
       do iVar = nDimHyades + 1, nDimHyades + nVarHyadesMesh
          DataHyadesMesh_VC(iVar,iCellMesh) = &
               Var_VI(iVarMesh_I(iVar)-nDimHyades,iCellMesh) &
               *Hyades2No_V(iVarMesh_I(iVar))
       end do

       ! The first rows in HYADES only contains mesh data
       if(UseZoneCenter .and. (i == 1 .or. j == 1)) CYCLE

       iCell = iCell + 1

       if(UseZoneCenter)then
          do iVar = 1, nDimHyades
             DataHyades_VC(iVar,iCell) = &
                  Var_VI(iVarZone_I(iVar)-nDimHyades,iCellMesh) &
                  *Hyades2No_V(iVarZone_I(iVar))
          end do
       else
          DataHyades_VC(:nDimHyades,iCell) = &
               DataHyadesMesh_VC(:nDimHyades,iCell)
       end if
       ! zone centered variables
       do iVar = nDimHyades + 1, nDimHyades + nVarHyades
          DataHyades_VC(iVar,iCell) = &
               Var_VI(iVarZone_I(iVar)-nDimHyades,iCellMesh) &
               *Hyades2No_V(iVarZone_I(iVar))
       end do
    end do; end do

    deallocate(Var_VI)

    if(iMaterialHyades > 0)then
       ! Vacuum (5) --> Polyimid
       where(nint(DataHyades_VC(iMaterialHyades, :)) == 5) &
            DataHyades_VC(iMaterialHyades, :) = Plastic_
       if(.not.UseAy)then
          ! Acrylic (4) --> Polyimid
          where(nint(DataHyades_VC(iMaterialHyades, :)) == Ay_) &
               DataHyades_VC(iMaterialHyades, :) = Plastic_
       end if
       if(.not.UseAu)then
          ! Gold (3) --> Polyimid
          where(nint(DataHyades_VC(iMaterialHyades, :)) == Au_) &
               DataHyades_VC(iMaterialHyades, :) = Plastic_
       end if
    end if

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

       ! Fix the pressure where it is set to some very small value
       where(DataHyades_VC(iPHyades, :) < 1e-10) &
            DataHyades_VC(iPHyades, :) = pDimOutside*Io2No_V(UnitP_)

       ! Find cell with maximum X coordinate along the symmetry axis
       iCellLastHyades     = nCellHyades_D(1)
       iCellLastHyadesMesh = nCellHyadesMesh_D(1)

       ! Calculate level set functions on the Hyades grid using 
       ! the minimum distance between cells of different materials
       allocate(LevelHyades_VC(0:nMaterial-1, nCellHyades))

       if(UseMixedCell)then
          ! Simply set 1.0 the levelset function corresponding to the material
          LevelHyades_VC = -1.0
          do iCell = 1, nCellHyades
             LevelHyades_VC(nint(DataHyades_VC(iMaterialHyades,iCell)),iCell) &
                  = 1.0
          end do
       else
          ! Determine distance functions
          allocate(Distance2_C(nCellHyades), iMaterial_C(nCellHyades))
          do iCell = 1, nCellHyades
             x         = DataHyades_VC(iXHyades, iCell)
             r         = DataHyades_VC(iRHyades, iCell)
             iMaterial = DataHyades_VC(iMaterialHyades, iCell)

             ! Distance squared from all other points
             Distance2_C = (x - DataHyades_VC(iXHyades,:))**2       &
                  +        (r - DataHyades_VC(iRHyades,:))**2

             ! Integer value of material in Hyades grid
             iMaterial_C = DataHyades_VC(iMaterialHyades,:)

             ! For each cell set 3 (or 4) level set functions
             do jMaterial = 0, nMaterial-1
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
    end if


    if(.not.UseHyadesGroupFile) RETURN

    ! read HYADES multi-group file
    call read_plot_file(NameHyadesGroupFile, nVarOut = nGroupHyades, &
         ParamOut_I = EnergyGroupHyades_I)

    ! Read in the data
    ! The HYADES file contains monochromatic radiation energy density
    allocate(Var_VI(nGroupHyades,nCellHyades))

    call read_plot_file(NameHyadesGroupFile, VarOut_VI = Var_VI)

    ! Convert the monochromatic radiation energy density to
    ! radiation energy density and convert from CGS to normalized units
    do iCell = 1, nCellHyades
       Var_VI(:,iCell) = Var_VI(:,iCell) &
            *(EnergyGroupHyades_I(1:nGroupHyades) &
            - EnergyGroupHyades_I(0:nGroupHyades-1)) &
            *0.1*Si2No_V(UnitEnergyDens_)   ! erg/cm^3 -> J/m^3
    end do

    ! Based on FreqMinSi and FreqMaxSi determine which Hyades groups
    ! are to be used.
    ! Initial minimum and maximum photon energy is in keV
    EnergyGroupMin = FreqMinSi*cHPlanckEV*1e-3
    EnergyGroupMax = FreqMaxSi*cHPlanckEV*1e-3

    ! Reset the minimum group energy (minimum group energy of HYADES does
    ! not follow a logarithmic scale)
    DeltaLogEnergy = &
         (log(EnergyGroupHyades_I(nGroupHyades))-log(EnergyGroupHyades_I(1))) &
         /(nGroupHyades - 1)
    EnergyGroupHyades_I(0) = &
         exp(log(EnergyGroupHyades_I(nGroupHyades)) &
         -   DeltaLogEnergy*nGroupHyades)

    ! Truncate the number of groups supplied by HYADES depending
    ! on the user supplied FreqMaxSi.
    ! Find maximum group index for which EnergyGroupHyades_I(nGroup)
    ! < EnergyGroupMax (units of keV)
    do iGroup = 1, nGroupHyades
       if(EnergyGroupHyades_I(iGroup) > EnergyGroupMax)then
          iGroupLast = iGroup - 1
          EXIT
       end if
       iGroupLast = iGroup
    end do

    ! Truncate from below according to FreqMinSi
    do iGroup = iGroupLast, 1, -1
       if(EnergyGroupHyades_I(iGroup-1) < EnergyGroupMin)then
          iGroupFirst = iGroup + 1
          EXIT
       end if
       iGroupFirst = iGroup
    end do

    nGroup = iGroupLast - iGroupFirst + 1

    if(nGroup /= nWave)then
       write(*,*)NameSub, 'nWave should be reset to ', nGroup
       call stop_mpi(NameSub//' reconfigure and recompile !')
    end if

    allocate(EradHyades_VC(nGroup,nCellHyades))

    do iCell = 1, nCellHyades
       EradHyades_VC(:,iCell) = Var_VI(iGroupFirst:iGroupLast,iCell)
    end do

    deallocate(Var_VI)

    ! convert minimum and maximum photon energy in keV to frequencies in Herz
    FreqMinSi = EnergyGroupHyades_I(iGroupFirst-1)*1e3/cHPlanckEV
    FreqMaxSi = EnergyGroupHyades_I(iGroupLast)*1e3/cHPlanckEV

  end subroutine read_hyades_file

  !============================================================================

  subroutine interpolate_hyades1d(iBlock)

    use BATL_size,           ONLY: nJ, nK, MinI, MaxI
    use CRASH_ModMultiGroup, ONLY: get_energy_g_from_temperature
    use ModAdvance,    ONLY: State_VGB, Rho_, RhoUx_, RhoUy_, RhoUz_, p_, &
         Erad_, UseElectronPressure
    use ModGeometry,   ONLY: Xyz_DGB
    use ModPhysics,    ONLY: Si2No_V, No2Si_V, UnitTemperature_, &
         cRadiationNo, UnitEnergyDens_
    use ModMain,       ONLY: UseRadDiffusion
    use ModVarIndexes, ONLY: nWave, WaveFirst_, WaveLast_

    integer, intent(in) :: iBlock

    integer :: i, j, k, iCell, iWave
    real :: x, Weight1, Weight2
    real :: Tr
    real :: TrSi, EgSi
    character(len=*), parameter :: NameSub='interpolate_hyades1d'
    !-------------------------------------------------------------------------
    do i = MinI, MaxI
       ! Find the Hyades points around this position
       x = Xyz_DGB(x_,i,1,1,iBlock)

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

       do k = 1, nK; do j = 1, nJ
          ! Interpolate density, momentum and pressure

          State_VGB(Rho_,i,j,k,iBlock) = &
               ( Weight1*DataHyades_VC(iRhoHyades, iCell-1) &
               + Weight2*DataHyades_VC(iRhoHyades, iCell) )

          State_VGB(RhoUx_,i,j,k,iBlock) =  State_VGB(Rho_,i,j,k,iBlock) * &
               ( Weight1*DataHyadesMesh_VC(iUxHyades, iCell-1) &
               + Weight2*DataHyadesMesh_VC(iUxHyades, iCell) )

          if(UseElectronPressure)then
             Te_G(i,j,k) = ( Weight1*DataHyades_VC(iTeHyades, iCell-1) &
                  +          Weight2*DataHyades_VC(iTeHyades, iCell) )
             Ti_G(i,j,k) = ( Weight1*DataHyades_VC(iTiHyades, iCell-1) &
                  +          Weight2*DataHyades_VC(iTiHyades, iCell) )
          end if

          State_VGB(p_,i,j,k,iBlock) = &
               ( Weight1*DataHyades_VC(iPHyades, iCell-1) &
               + Weight2*DataHyades_VC(iPHyades, iCell) )

          ! Set transverse momentum to zero
          State_VGB(RhoUy_:RhoUz_,i,j,k,iBlock) = 0.0

          if(UseRadDiffusion)then
             if(UseHyadesGroupFile)then
                State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock) = &
                     ( Weight1*EradHyades_VC(:,iCell-1) &
                     + Weight2*EradHyades_VC(:,iCell) )
             else
                ! Start from hyades radiation temperature
                ! Total radiation energy = cRadiation*Trad**4
                Tr = ( Weight1*DataHyades_VC(iTrHyades, iCell-1) &
                     + Weight2*DataHyades_VC(iTrHyades, iCell) )

                if(nWave ==1)then
                   State_VGB(Erad_,i,j,k,iBlock) = cRadiationNo*Tr**4
                else
                   TrSi = Tr*No2Si_V(UnitTemperature_)
                   do iWave = 1, nWave
                      call get_energy_g_from_temperature(iWave, TrSi,EgSI=EgSi)
                      State_VGB(WaveFirst_+iWave-1,i,j,k,iBlock) = &
                           EgSi*Si2No_V(UnitEnergyDens_)
                   end do
                end if
             end if
          end if

       end do; end do
    end do

  end subroutine interpolate_hyades1d

  !============================================================================

  subroutine interpolate_hyades2d(iBlock)

    ! Use triangulation to interpolate Hyades grid onto CRASH grid

    use CRASH_ModMultiGroup, ONLY: get_energy_g_from_temperature
    use ModSize,        ONLY: nI, nJ, nK
    use ModAdvance,     ONLY: State_VGB, Rho_, RhoUx_, RhoUy_, RhoUz_, p_, &
         Erad_, UseElectronPressure
    use ModGeometry,    ONLY: Xyz_DGB, y2
    use ModTriangulate, ONLY: calc_triangulation, mesh_triangulation, &
         find_triangle
    use ModMain,        ONLY: UseRadDiffusion
    use ModPhysics,     ONLY: cRadiationNo, No2Si_V, Si2No_V, &
         UnitTemperature_, UnitEnergyDens_
    use ModVarIndexes,  ONLY: nWave, WaveFirst_, WaveLast_

    integer, intent(in) :: iBlock

    integer, save              :: nTriangle, nTriangleMesh
    integer, allocatable, save :: &
         iNodeTriangle_II(:,:), iNodeTriangleMesh_II(:,:), &
         nTriangle_C(:,:), nTriangleMesh_C(:,:)
    integer, pointer, save     :: iTriangle_IC(:,:,:), iTriangleMesh_IC(:,:,:)
    real, allocatable,    save :: DataHyades_V(:), DataHyadesMesh_V(:)
    real                       :: LevelHyades_V(0:nMaterial-1)
    real                       :: EradHyades_V(nWave)

    integer :: i, j, k, iNode1, iNode2, iNode3, iWave
    real    :: x, y, z, r, Weight1, Weight2, Weight3
    real    :: WeightNode_I(3), WeightMaterial_I(0:nMaterial-1), Weight
    real    :: TrSi, EgSi

    integer :: iMaterial, iMaterialNode_I(3)

    character(len=*), parameter :: NameSub='interpolate_hyades2d'
    !-------------------------------------------------------------------------
    if(.not.allocated(iNodeTriangle_II))then
       ! allocate variables and do triangulation
       allocate(iNodeTriangle_II(3,2*nCellHyades))
       allocate(DataHyades_V(nDimHyades + nVarHyades))

       CoordHyades_DC = DataHyades_VC( (/iXHyades, iRHyades/),:)
       ! Make CoordHyadesMesh_DC dimensionless
       if(UseZoneCenter) CoordHyadesMesh_DC = DataHyadesMesh_VC(:nDimHyades,:)

       if(UseDelaunay)then
          call calc_triangulation( &
               nCellHyades, CoordHyades_DC, iNodeTriangle_II, nTriangle)
       else
          call mesh_triangulation(nCellHyades_D(1), nCellHyades_D(2), &
               CoordHyades_DC, iNodeTriangle_II, nTriangle)
       end if
       ! Use a uniform mesh for fast search algorithm
       allocate(nTriangle_C(nCellHyades_D(1),nCellHyades_D(2)))
       nTriangle_C(1,1) = -1

       allocate(DataHyadesMesh_V(nDimHyades + nVarHyadesMesh))
       if(UseZoneCenter)then
          allocate(iNodeTriangleMesh_II(3,2*nCellHyadesMesh))
          if(UseDelaunay)then
             call calc_triangulation( nCellHyadesMesh, CoordHyadesMesh_DC, &
                  iNodeTriangleMesh_II, nTriangleMesh)
          else
             call mesh_triangulation( &
                  nCellHyadesMesh_D(1), nCellHyadesMesh_D(2), &
                  CoordHyadesMesh_DC, iNodeTriangleMesh_II, nTriangleMesh)
          end if
          allocate( &
               nTriangleMesh_C(nCellHyadesMesh_D(1),nCellHyadesMesh_D(2)))
          nTriangleMesh_C(1,1) = -1
       end if
    end if

    ! Interpolate points 
    do j = 1, nJ; do i = 1, nI; do k = 1, nk

       x = Xyz_DGB(x_,i,j,k,iBlock)
       y = Xyz_DGB(y_,i,j,k,iBlock)
       z = Xyz_DGB(z_,i,j,k,iBlock)

       call set_yzr(x, y, z, r, rMax=y2)

       ! Check if we are at the end of the Hyades grid
       if(x >= DataHyades_VC(iXHyades,iCellLastHyades))then
          iNode1 = iCellLastHyades;  Weight1 = 1.0
          iNode2 = 1;                Weight2 = 0.0
          iNode3 = 1;                Weight3 = 0.0
       else
          ! Find the Hyades triangle around this position
          call find_triangle(&
               nCellHyades, nTriangle, &
               (/x, r/), CoordHyades_DC, &
               iNodeTriangle_II, &
               iNode1, iNode2, iNode3, Weight1, Weight2, Weight3, &
               nTriangle_C=nTriangle_C, iTriangle_IC=iTriangle_IC)
       end if

       ! Check if the 3 points consist of the same material or not
       ! If the materials are different use the points with the
       ! material that has the largest total weight

       ! Weight and material of the 3 nodes of the triangle
       WeightNode_I    = (/ Weight1, Weight2, Weight3 /)
       iMaterialNode_I = DataHyades_VC(iMaterialHyades, &
            (/iNode1, iNode2, iNode3/) )

       if(maxval(iMaterialNode_I) /= minval(iMaterialNode_I))then

          ! Add up the weights for all materials
          do iMaterial = 0, nMaterial - 1
             WeightMaterial_I(iMaterial) = sum(WeightNode_I, &
                  MASK = (iMaterialNode_I == iMaterial) )
          end do

          ! Find the dominant material
          iMaterial   = maxloc(WeightMaterial_I, 1) - 1
          Weight      = WeightMaterial_I(iMaterial)

          where(iMaterialNode_I == iMaterial) 
             ! Reset weights so they add up to 1 for the dominant material
             WeightNode_I = WeightNode_I / Weight
          elsewhere
             ! Other materials get zero weight
             WeightNode_I = 0.0
          end where

       end if

       DataHyades_V = &
            WeightNode_I(1)*DataHyades_VC(:, iNode1) + &
            WeightNode_I(2)*DataHyades_VC(:, iNode2) + &
            WeightNode_I(3)*DataHyades_VC(:, iNode3)

       LevelHyades_V = &
            WeightNode_I(1)*LevelHyades_VC(:, iNode1) + &
            WeightNode_I(2)*LevelHyades_VC(:, iNode2) + &
            WeightNode_I(3)*LevelHyades_VC(:, iNode3)

       if(UseHyadesGroupFile)then
          EradHyades_V = &
               WeightNode_I(1)*EradHyades_VC(:,iNode1) + &
               WeightNode_I(2)*EradHyades_VC(:,iNode2) + &
               WeightNode_I(3)*EradHyades_VC(:,iNode3)
       end if

       if(UseZoneCenter)then
          ! Check if we are at the end of the Hyades grid
          if(x >= DataHyadesMesh_VC(iXHyadesMesh, iCellLastHyadesMesh))then
             iNode1 = iCellLastHyadesMesh;  Weight1 = 1.0
             iNode2 = 1;                    Weight2 = 0.0
             iNode3 = 1;                    Weight3 = 0.0
          else
             ! Find the Hyades triangle around this position
             call find_triangle(&
                  nCellHyadesMesh, nTriangleMesh, (/x, r/), &
                  CoordHyadesMesh_DC, &
                  iNodeTriangleMesh_II, &
                  iNode1, iNode2, iNode3, Weight1, Weight2, Weight3, &
                  nTriangle_C=nTriangleMesh_C, &
                  iTriangle_IC=iTriangleMesh_IC)
          end if
          ! Weight of the 3 nodes of the triangle
          WeightNode_I    = (/ Weight1, Weight2, Weight3 /)
       end if
       DataHyadesMesh_V = &
            WeightNode_I(1)*DataHyadesMesh_VC(:, iNode1) + &
            WeightNode_I(2)*DataHyadesMesh_VC(:, iNode2) + &
            WeightNode_I(3)*DataHyadesMesh_VC(:, iNode3)

       ! Interpolate density, momentum and pressure

       State_VGB(Rho_,i,j,k,iBlock)  = DataHyades_V(iRhoHyades)

       State_VGB(RhoUx_,i,j,k,iBlock) = &
            DataHyades_V(iRhoHyades) * DataHyadesMesh_V(iUxHyades)

       State_VGB(RhoUy_:RhoUz_,i,j,k,iBlock) = (/y, z/)/r * &
            DataHyades_V(iRhoHyades) * DataHyadesMesh_V(iUrHyades)

       ! Interpolate level set functions
       State_VGB(LevelXe_:LevelMax,i,j,k,iBlock) = LevelHyades_V

       if(UseElectronPressure)then
          Te_G(i,j,k) = DataHyades_V(iTeHyades)
          Ti_G(i,j,k) = DataHyades_V(iTiHyades)
       end if
       State_VGB(p_,i,j,k,iBlock)  = DataHyades_V(iPHyades)

       if(UsePl .and. AssFactor > 0.0)then
          ! Reduce plastic pressure
          if(State_VGB(LevelPl_,i,j,k,iBlock)>0. .or. &
               UseAy .and. State_VGB(LevelAy_,i,j,k,iBlock)>0.) &
               State_VGB(p_,i,j,k,iBlock) &
               = AssFactor*State_VGB(p_,i,j,k,iBlock)

          ! Reduce inward velocity
          if(State_VGB(RhoUy_,i,j,k,iBlock)<0.0) &
               State_VGB(RhoUy_,i,j,k,iBlock) &
               = AssFactor*State_VGB(RhoUy_,i,j,k,iBlock)
       end if

       if(UseRadDiffusion)then
          if(UseHyadesGroupFile)then
             State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock) = EradHyades_V
          else
             ! Start from hyades radiation temperature
             ! Total radiation energy = cRadiation*Trad**4
             if(nWave == 1)then
                State_VGB(Erad_,i,j,k,iBlock) = &
                     cRadiationNo * DataHyades_V(iTrHyades)**4
             else
                TrSi = DataHyades_V(iTrHyades)*No2Si_V(UnitTemperature_)
                do iWave = 1, nWave
                   call get_energy_g_from_temperature(iWave, TrSi, EgSI = EgSi)
                   State_VGB(WaveFirst_+iWave-1,i,j,k,iBlock) = &
                        EgSi*Si2No_V(UnitEnergyDens_)
                end do
             end if
          end if
       end if

    end do; end do; end do

  end subroutine interpolate_hyades2d

  !============================================================================

  subroutine user_update_states(iStage,iBlock)

    use ModSize,     ONLY: nI, nJ, nK
    use ModAdvance,  ONLY: State_VGB, p_, ExtraEint_, &
         UseNonConservative, IsConserv_CB, UseElectronPressure
    use ModPhysics,  ONLY: inv_gm1, Si2No_V, No2Si_V, &
         UnitP_, UnitEnergyDens_, ExtraEintMin, UnitTemperature_
    use ModEnergy,   ONLY: calc_energy_cell
    use ModVarIndexes, ONLY: Pe_, Te0_, nWave, WaveFirst_, WaveLast_
    use ModVarIndexes, ONLY: DefaultState_V
    use ModSize, ONLY:MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use CRASH_ModMultiGroup, ONLY: get_planck_g_from_temperature

    integer, intent(in):: iStage,iBlock

    integer:: i, j, k, iP, iGroup
    real   :: PressureSi, EinternalSi, Te0Si, PlanckSi
    logical:: IsConserv

    character(len=*), parameter :: NameSub = 'user_update_states'
    !------------------------------------------------------------------------
    if(UseNLTE.and. iStage==1) then
       !\
       !Reset value of E_rad/B
       !/
       !Ghost cells are needed to calculate temperature or electron density 
       !gradients
       do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI

          if(State_VGB(Te0_,i,j,k,iBlock)==DefaultState_V(Te0_))then
             !Not filled-in ghost cell or initial conditions
             EOverB_VGB(:,i,j,k,iBlock) = 1.0
             CYCLE
          end if
          !Calculate "E_rad"
          EOverB_VGB(:,i,j,k,iBlock) = No2Si_V(UnitEnergyDens_)* &
               State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock)
             
          !Calculate "over B"
          Te0Si = State_VGB(Te0_,i,j,k,iBlock) * No2Si_V(UnitTemperature_)
          do iGroup = 1, nWave
             call get_planck_g_from_temperature(iGroup, Te0Si, EgSI=PlanckSi)
             if(EOverB_VGB(iGroup,i,j,k,iBlock)>0.0)then
                EOverB_VGB(iGroup,i,j,k,iBlock) = EOverB_VGB(iGroup,i,j,k,iBlock)/&
                     max(PlanckSi, 0.01*EOverB_VGB(iGroup,i,j,k,iBlock))
             else
                EOverB_VGB(iGroup,i,j,k,iBlock)=0.0
             end if
          end do
       end do; end do; end do
    end if


    call update_states_MHD(iStage,iBlock)

    iP = p_
    if(UseElectronPressure) iP = Pe_

    ! update of pressure, ionization and total energies
    do k=1,nK; do j=1,nJ; do i=1,nI

       if(allocated(IsConserv_CB))then
          IsConserv = IsConserv_CB(i,j,k,iBlock)
       else
          IsConserv = .not. UseNonConservative
       end if

       if(IsConserv)then
          ! Single temperature:
          !   From update_states_MHD, we obtained p^*, e^* with ideal gamma
          !   and ExtraEInt^* with pure advection.
          !   Total energy density E^n+1  = e^* + ExtraEInt^* is conserved.
          !   Total internal energy Eint^n+1 = p^*/(g-1) + ExtraEInt^*
          ! Two temperature:
          !   From update_states_MHD, we obtained Pe^*, e^* with ideal gamma
          !   and ExtraEInt^* with pure advection.
          !   Total energy density E^n+1  = e^* + ExtraEInt^* is conserved.
          !   Electron internal energy Eint^n+1 = Pe^*/(g-1) + ExtraEInt^*
          EinternalSi = No2Si_V(UnitEnergyDens_)*&
               (inv_gm1*State_VGB(iP,i,j,k,iBlock) + &
               State_VGB(ExtraEint_,i,j,k,iBlock))

          ! Single temperature: determine p^n+1 = EOS( rho^n+1, Eint^n+1)
          ! Two temperature:   determine Pe^n+1 = EOS( rho^n+1, Eint^n+1)
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, &
               EinternalIn=EinternalSi, PressureOut=PressureSi)

          ! Set normalized pressure (electron pressure for two temperature)
          State_VGB(iP,i,j,k,iBlock) = PressureSi*Si2No_V(UnitP_)
       else
          ! Single temperature:
          !   Calculate total internal energy Eint^n+1 from pressure p^n+1
          ! Two temperature:
          !   Calculate electron internal energy Eint^n+1 from
          !   electron pressure Pe^n+1
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, EinternalOut=EinternalSi)
       end if

       ! Set ExtraEint^n+1 = Eint^n+1 - p^n+1/(g -1)
       State_VGB(ExtraEint_,i,j,k,iBlock) = max(ExtraEintMin, &
            Si2No_V(UnitEnergyDens_)*EinternalSi &
            - inv_gm1*State_VGB(iP,i,j,k,iBlock))

    end do; end do; end do

    ! Calculate e^n+1 using updated P^n+1 (Pe^n+1) for single (two-)temperature
    ! Note that this is identical to e^n+1 = e^* + ExtraEInt^* - ExtraEint^n+1
    call calc_energy_cell(iBlock)

  end subroutine user_update_states

  !============================================================================

  subroutine user_calc_sources(iBlock)

    use ModMain,       ONLY: nI, nJ, nK
    use ModAdvance,    ONLY: State_VGB, &
         Source_VC, uDotArea_XI, uDotArea_YI, uDotArea_ZI, &
         IsConserv_CB, UseNonConservative, UseElectronPressure
    use BATL_lib,      ONLY: CellVolume_GB
    use ModPhysics,    ONLY: g
    use ModVarIndexes, ONLY: p_, Pe_

    integer, intent(in) :: iBlock

    integer :: i, j, k, iP
    real :: DivU, GammaEos
    logical :: IsConserv

    character (len=*), parameter :: NameSub = 'user_calc_sources'
    !-------------------------------------------------------------------
    if(.not.(UseNonConsLevelSet .or. UseNonConservative)) RETURN

    iP = p_
    if(UseElectronPressure) iP = Pe_

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       ! Note that the velocity of the first (and only) fluid is used
       DivU            =        uDotArea_XI(i+1,j,k,1) - uDotArea_XI(i,j,k,1)
       if(nJ > 1) DivU = DivU + uDotArea_YI(i,j+1,k,1) - uDotArea_YI(i,j,k,1)
       if(nK > 1) DivU = DivU + uDotArea_ZI(i,j,k+1,1) - uDotArea_ZI(i,j,k,1)
       DivU = DivU/CellVolume_GB(i,j,k,iBlock)


       ! Add Level*div(u) as a source term so level sets beome advected scalars
       if(UseNonConsLevelSet) &
            Source_VC(LevelXe_:LevelMax,i,j,k) = &
            Source_VC(LevelXe_:LevelMax,i,j,k) &
            + State_VGB(LevelXe_:LevelMax,i,j,k,iBlock)*DivU

       if(allocated(IsConserv_CB))then
          IsConserv = IsConserv_CB(i,j,k,iBlock)
       else
          IsConserv = .not. UseNonConservative
       end if

       ! Fix adiabatic compression source for pressure
       if(.not.IsConserv)then
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, GammaOut=GammaEos)

          Source_VC(iP,i,j,k) = Source_VC(iP,i,j,k) &
               -(GammaEos-g)*State_VGB(iP,i,j,k,iBlock)*DivU
       end if
    end do; end do; end do

  end subroutine user_calc_sources

  !===========================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModConst,   ONLY: cKtoKev, cBoltzmann
    use ModAdvance, ONLY: State_VGB, UseElectronPressure
    use ModPhysics, ONLY: No2Si_V, No2Io_V, UnitRho_, UnitP_, &
         UnitTemperature_, cRadiationNo, No2Si_V
    use ModGeometry, ONLY: r_BLK, Xyz_DGB
    use ModVarIndexes, ONLY: Rho_, p_, nWave, WaveFirst_, WaveLast_
    use CRASH_ModEos, ONLY: Xe_, Be_, Plastic_, Au_, Ay_
    use BATL_size,    ONLY: nI, nJ, nK, MinI, MaxI
    use BATL_lib,     ONLY: IsRzGeometry

    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    real    :: WaveEnergy
    real    :: PiSi, TiSi, NatomicSi
    real    :: OpacityPlanckSi_W(nWave)
    real    :: OpacityRosselandSi_W(nWave)
    integer :: i, j, k, iMaterial, iLevel, iWave

    ! Do not use MinJ,MinK,MaxJ,MaxK here to avoid pgf90 compilation error...
    integer, parameter:: jMin = 1 - 2*min(1,nJ-1), jMax = nJ + 2*min(1,nJ-1)
    integer, parameter:: kMin = 1 - 2*min(1,nK-1), kMax = nK + 2*min(1,nK-1)

    ! Optical depth for Xenon, Berylium, plastic and Gold for radiography
    ! The optical depth of Acrylic is taken to be the same as polyimide
    real, parameter:: RadioDepth_I(0:MaxMaterial-1) = &
         (/ 79.4, 0.36, 2.24, 1e10, 2.24 /)

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    !------------------------------------------------------------------------  
    IsFound = .true.
    select case(NameVar)
    case('level', 'material')
       do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
          PlotVar_G(i,j,k) = &
               maxloc(State_VGB(LevelXe_:LevelMax,i,j,k,iBlock), 1) - 1
       end do; end do; end do
    case('tekev', 'TeKev')
       NameIdlUnit = 'KeV'
       !if(Te0_>1)then
       !   do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
       !      PlotVar_G(i,j,k) = PlotVar_G(i,j,k) * &
       !           No2Si_V(UnitTemperature_) * cKToKev
       !   end do; end do; end do
       !else
          do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, TeOut = PlotVar_G(i,j,k))
             PlotVar_G(i,j,k) = PlotVar_G(i,j,k) * cKToKev
          end do; end do; end do
       !end if
    case('tikev', 'TiKev')
       NameIdlUnit = 'KeV'
       if(UseElectronPressure)then
          do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, NatomicOut=NatomicSi)
             PiSi = State_VGB(p_,i,j,k,iBlock)*No2Si_V(UnitP_)
             TiSi = PiSi/(cBoltzmann*NatomicSi)
             PlotVar_G(i,j,k) = TiSi*cKToKev
          end do; end do; end do
       else
          ! Te = Ti at all times, use Te
          do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  i, j, k, iBlock, TeOut = PlotVar_G(i,j,k))
             PlotVar_G(i,j,k) = PlotVar_G(i,j,k) * cKToKev
          end do; end do; end do
       end if
    case('tradkev','trkev')
       ! radiation temperature is physically meaningless, but only
       ! used as a measure of the total radiation energy !!!
       ! multiply by sign of Erad for debugging purpose
       NameIdlUnit = 'KeV'
       do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
          WaveEnergy = 0.0
          do iWave = WaveFirst_, WaveLast_
             WaveEnergy = WaveEnergy + State_VGB(iWave,i,j,k,iBlock)
          end do
          PlotVar_G(i,j,k) = sign(1.0,WaveEnergy) &
               *sqrt(sqrt(abs(WaveEnergy)/cRadiationNo))&
               * No2Si_V(UnitTemperature_) * cKToKev
       end do; end do; end do
    case('planck')
       do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, &
               OpacityPlanckOut_W = OpacityPlanckSi_W)
          PlotVar_G(i,j,k) = OpacityPlanckSi_W(1)
       end do; end do; end do
    case('ross')
       do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, &
               OpacityRosselandOut_W = OpacityRosselandSi_W)
          PlotVar_G(i,j,k) = OpacityRosselandSi_W(1)
       end do; end do; end do
    case('zavg')
       do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, AverageIonChargeOut = PlotVar_G(i,j,k))
       end do; end do; end do
    case('cond')
       do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, HeatCondOut = PlotVar_G(i,j,k))
       end do; end do; end do
    case('teti')
       do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, TeTiRelaxOut = PlotVar_G(i,j,k))
       end do; end do; end do
    case('gamma')
       do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, GammaOut = PlotVar_G(i,j,k))
       end do; end do; end do
    case('usersphere10')
       ! Test function for LOS images: sphere with "density" 
       !    100 - r^2 inside r=10, and 0 outside.
       if(IsRzGeometry)then
          ! In R-Z geometry the "sphere" is a circle in the X-Y plane
          PlotVar_G = max(0.0, &
               100 - Xyz_DGB(x_,:,:,:,iBlock)**2 - Xyz_DGB(y_,:,:,:,iBlock)**2)
       else
          ! In Cartesian geometry it is real sphere
          PlotVar_G = max(0.0, 100 - r_BLK(:,:,:,iBlock)**2)
       end if
    case('radiograph')
       if(UseMixedCell)then
          do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
             PlotVar_G(i,j,k) = -sum(RadioDepth_I(0:nMaterial-1) &
                  *MassMaterial_I(0:nMaterial-1) &
                  *State_VGB(LevelXe_:LevelMax,i,j,k,iBlock))
          end do; end do; end do
       else
          do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
             iMaterial = &
                  maxloc(State_VGB(LevelXe_:LevelMax,i,j,k,iBlock), 1) - 1
             PlotVar_G(i,j,k) = -RadioDepth_I(iMaterial) &
                  *State_VGB(Rho_,i,j,k,iBlock)
          end do; end do; end do
       end if
       if(IsDimensional) PlotVar_G(MinI:MaxI,jMin:jMax,kMin:kMax) = &
            No2Io_V(UnitRho_)*PlotVar_G(MinI:MaxI,jMin:jMax,kMin:kMax)
    case('xe', 'be', 'pl', 'au', 'ay')
       select case(NameVar)
       case('xe')
          iLevel = LevelXe_; iMaterial = Xe_
       case('be')
          iLevel = LevelBe_; iMaterial = Be_
       case('pl')
          iLevel = LevelPl_; iMaterial = Plastic_
       case('au')
          iLevel = LevelAu_; iMaterial = Au_
       case('ay')
          iLevel = LevelAy_; iMaterial = Ay_
       end select
       PlotVar_G(MinI:MaxI,jMin:jMax,kMin:kMax) = &
            State_VGB(iLevel,MinI:MaxI,jMin:jMax,kMin:kMax,iBlock)

    case('rhoxe', 'rhobe', 'rhopl', 'rhoau', 'rhoay')
       select case(NameVar)
       case('rhoxe')
          iLevel = LevelXe_; iMaterial = Xe_
       case('rhobe')
          iLevel = LevelBe_; iMaterial = Be_
       case('rhopl')
          iLevel = LevelPl_; iMaterial = Plastic_
       case('rhoau')
          iLevel = LevelAu_; iMaterial = Au_
       case('rhoay')
          iLevel = LevelAy_; iMaterial = Ay_
       end select
       if(UseMixedCell)then
          PlotVar_G(MinI:MaxI,jMin:jMax,kMin:kMax) = MassMaterial_I(iMaterial)&
               *State_VGB(iLevel,MinI:MaxI,jMin:jMax,kMin:kMax,iBlock)
       else
          do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
             if(iMaterial == &
                  maxloc(State_VGB(LevelXe_:LevelMax,i,j,k,iBlock),1) - 1) then
                PlotVar_G(i,j,k) = State_VGB(Rho_,i,j,k,iBlock)
             else
                PlotVar_G(i,j,k) = 0.0
             end if
          end do; end do; end do
       end if
       if(IsDimensional) PlotVar_G(MinI:MaxI,jMin:jMax,kMin:kMax) = &
            No2Io_V(UnitRho_)*PlotVar_G(MinI:MaxI,jMin:jMax,kMin:kMax)
    case default
       IsFound = .false.
    end select

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

  end subroutine user_set_plot_var

  !============================================================================
  subroutine user_get_log_var(VarValue, TypeVar, Radius)

    use ModAdvance,    ONLY: State_VGB, tmp1_BLK
    use ModVarIndexes, ONLY: Rho_
    use CRASH_ModEos,  ONLY: Xe_, Be_, Plastic_, Au_, Ay_
    use ModMain,       ONLY: nI, nJ, nK, nBlock, Unused_B
    use ModGeometry,   ONLY: DomainVolume

    real, external :: integrate_BLK

    real, intent(out)            :: VarValue
    character (len=*), intent(in):: TypeVar
    real, intent(in), optional :: Radius

    integer:: iMaterial, iLevel, iBlock, i, j, k

    character (len=*), parameter :: NameSub = 'user_get_log_var'
    !-------------------------------------------------------------------------
    select case(TypeVar)
    case('rhoxe', 'rhobe', 'rhopl', 'rhoau', 'rhoay')
       select case(TypeVar)
       case('rhoxe')
          iLevel = LevelXe_; iMaterial = Xe_
       case('rhobe')
          iLevel = LevelBe_; iMaterial = Be_
       case('rhopl')
          iLevel = LevelPl_; iMaterial = Plastic_
       case('rhoau')
          iLevel = LevelAu_; iMaterial = Au_
       case('rhoay')
          iLevel = LevelAy_; iMaterial = Ay_
       end select
       tmp1_BLK(:,:,:,1:nBlock) = 0.0
       if(UseMixedCell)then
          do iBlock = 1, nBlock
             if(Unused_B(iBlock)) CYCLE
             tmp1_BLK(1:nI,1:nJ,1:nK,iBlock) = MassMaterial_I(iMaterial) &
                  *State_VGB(iLevel,1:nI,1:nJ,1:nK,iBlock)
          end do
       else
          do iBlock = 1, nBlock
             if(Unused_B(iBlock)) CYCLE
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                if(iMaterial == &
                     maxloc(State_VGB(LevelXe_:LevelMax,i,j,k,iBlock),1) - 1) &
                     tmp1_BLK(i,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock)
             end do; end do; end do
          end do
       end if
       VarValue = integrate_BLK(1,tmp1_BLK)/DomainVolume
    case default
       VarValue = -7777.0
    end select

  end subroutine user_get_log_var

  !===========================================================================

  subroutine user_init_session

    use ModProcMH,      ONLY: iProc, iComm
    use ModVarIndexes,  ONLY: Rho_, UnitUser_V, Te0_
    use ModLookupTable, ONLY: i_lookup_table, make_lookup_table
    use ModPhysics,     ONLY: &
         No2Io_V, UnitX_
    use ModConst,       ONLY: cHPlanckEV
    use ModWaves,       ONLY: nWave, FreqMinSI, FreqMaxSI, DoAdvectWaves
    use ModIo,          ONLY: restart
    use CRASH_ModMultiGroup, ONLY: nGroup, IsLogMultiGroupGrid
    use CRASH_ModEos,        ONLY: nMaterialEos, NameMaterial_I 
    use CRASH_ModEosTable,   ONLY: check_eos_table, check_opac_table
    use ModSize,             ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, MaxBlock
    use ModVarIndexes,       ONLY: nWave
    use CRASH_ModInterfaceNLTE, ONLY: check_nlte,radiom_version=>printversion

    integer:: iMaterial
    logical:: IsFirstTime = .true.
    logical,save:: UseNLTESaved = .false.
    character (len=*), parameter :: NameSub = 'user_init_session'
    !-------------------------------------------------------------------

    if(UseGammaLaw .and. (UseMixedCell .or. UseVolumeFraction)) &
         call stop_mpi(NameSub// &
         " can not use gamma-law with mixed cell or volume fraction")
    
    ! Pass number of materials and waves=groups to the CRASH EOS library
    nMaterialEos = nMaterial
    nGroup       = nWave

    ! The units always have to be reset, because set_physics sets them
    if(UseMixedCell) then
       UnitUser_V(LevelXe_:LevelMax) = UnitUser_V(Rho_)
    elseif(UseNonConsLevelSet)then
       UnitUser_V(LevelXe_:LevelMax) = No2Io_V(UnitX_)
    else
       UnitUser_V(LevelXe_:LevelMax) = UnitUser_V(Rho_)*No2Io_V(UnitX_)
    end if

    if( .not.IsFirstTime .and. UseNLTE .and. .not. UseNLTESaved) &
         call CON_stop(NameSub// &
         ': NLTE can be turned on only in the first session')
 
    ! The rest of the initialization should be done once
    if(.not.IsFirstTime) RETURN
    IsFirstTime = .false.

    ! Create separate EOS table for each material
    ! This should be done with UseNLTE = .false., because
    ! to create the EOS tables with UseNLTE = .true. would
    ! create NLTE tables for E_rad/B = 0
    UseNLTESaved = UseNLTE    
    UseNLTE      = .false.

    call check_eos_table(iComm = iComm)

    UseNLTE = UseNLTESaved

    !Set the photon energy range
    !First, check if the values of FreqMinSI and FreqSI are set:

    !If the frequency range IN HERZ has been alredy set, then skip
    if(FreqMinSI <= 0) then
       !Reset the minimum photon energy to be 0.1 eV
       FreqMinSI = 0.1 /cHPlanckEV
    end if

    if(FreqMaxSI <= 0) then
       !Reset the maximum photon energy to be 10 keV
       FreqMaxSI = 10000.0 /cHPlanckEV
    end if

    ! Read in Hyades output
    if(UseHyadesFile .and. .not. restart) call read_hyades_file

    ! Read in 2D CRASH output
    if(UseCrash2DFile .and. .not. restart) call read_crash2d_file

    ! create separate opacity lookup tables for each material
    call check_opac_table(FreqMinSI, FreqMaxSI, iComm = iComm)

    if(UseNLTE)then
       call check_nlte
       if(iProc==0)call radiom_version
       if(Te0_==1)call CON_stop('Use state vector with Te0_ component with NLTE!')
       if(.not.allocated(EOverB_VGB))then
          allocate(EOverB_VGB(nWave,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
          EOverB_VGB = 1.0
       end if
    end if

    ! The above call sets IsLogMultigroupGrid variable
    ! If the energy grid is not logarithmic, do not advect waves
    if(.not.IsLogMultigroupGrid) DoAdvectWaves = .false.

    ! these EOS tables (if used) combine multiple materials
    iTablePPerRho   = i_lookup_table('pPerRho(e/rho,rho)')
    iTableThermo    = i_lookup_table('Thermo(rho,p/rho)')
    iTableOpacity   = i_lookup_table('Opacity(rho,T)')

    ! another form of opacity tables
    do iMaterial = 0, nMaterial - 1
       iTableOpacity_I(iMaterial) = &
            i_lookup_table('Opacity'//NameMaterial_I(iMaterial)//'(rho,T)')
    end do

    ! SESAME gray opacity table
    iTableSesame = i_lookup_table('Sesame(rho,T)')

    if(iProc==0) write(*,*) NameSub, &
         ' iTablePPerRho, Thermo, Opacity, Opacity_I, iTableSesame = ', &
         iTablePPerRho, iTableThermo, iTableOpacity, &
         iTableOpacity_I, iTableSesame

    if(iTablePPerRho > 0) &
         call make_lookup_table(iTablePPerRho, calc_table_value, iComm)
    if(iTableThermo > 0) &
         call make_lookup_table(iTableThermo, calc_table_value, iComm)
    if(iTableOpacity > 0) &
         call make_lookup_table(iTableOpacity, calc_table_value, iComm)
    do iMaterial = 0, nMaterial-1
       if(iTableOpacity_I(iMaterial) > 0) &
            call make_lookup_table(iTableOpacity_I(iMaterial), &
            calc_table_value, iComm)
    end do

    if(iTableSesame > 0) &
         call make_lookup_table(iTableSesame, calc_table_value, iComm)

  end subroutine user_init_session

  !===========================================================================
  subroutine calc_table_value(iTable, Arg1, Arg2, Value_V)

    use ModAdvance,     ONLY: UseElectronPressure
    use ModProcMH,      ONLY: iProc
    use CRASH_ModEos,   ONLY: eos, Ay_, Plastic_
    use ModConst,       ONLY: cProtonMass, cBoltzmann
    use ModLookupTable, ONLY: interpolate_lookup_table
    use ModVarIndexes,  ONLY: nWave

    integer, intent(in):: iTable
    real, intent(in)   :: Arg1, Arg2
    real, intent(out)  :: Value_V(:)

    real:: Rho, p, e, Te, Cv, Gamma, Zavg, HeatCond, TeTiRelax
    real:: OpacityPlanck_W(nWave), OpacityRosseland_W(nWave)
    integer:: iMaterial
    character(len=*), parameter:: NameSub = 'ModUser::calc_table_value'
    !-----------------------------------------------------------------------

    if(iTable == iTablePPerRho)then
       ! Calculate p/e for Xe_, Be_ and Plastic_ for given Rho and e/Rho
       ! Au_ and Ay_ are optional
       Rho = Arg2
       e   = Arg1*Rho
       do iMaterial = 0, nMaterial-1
          if(UseElectronPressure)then
             call eos(iMaterial, Rho, eElectronIn=e, pElectronOut=p)
          else
             call eos(iMaterial, Rho, EtotalIn=e, pTotalOut=p)
          end if

          ! Material index starts from 0 :-( hence the +1
          Value_V(iMaterial+1) = p/Rho
       end do
    elseif(iTable == iTableThermo)then
       ! Calculate cV, gamma, HeatCond and Te for Xe_, Be_ and Plastic_ 
       ! for given Rho and p/Rho
       ! Au_ and Ay_ are optional
       Rho = Arg1
       p   = Arg2*Rho
       do iMaterial = 0, nMaterial-1
          if(UseElectronPressure)then
             call eos(iMaterial, Rho, pElectronIn=p, &
                  TeOut=Te, CvElectronOut=Cv, GammaEOut=Gamma, zAverageOut=Zavg, &
                  HeatCond=HeatCond, TeTiRelax=TeTiRelax)

             Value_V(Te_   +iMaterial*nThermo) = Te
             Value_V(Cv_   +iMaterial*nThermo) = Cv
             Value_V(Gamma_+iMaterial*nThermo) = Gamma
             Value_V(Zavg_ +iMaterial*nThermo) = Zavg
             Value_v(Cond_ +iMaterial*nThermo) = HeatCond
             Value_V(TeTi_ +iMaterial*nThermo) = TeTiRelax
          else

             call eos(iMaterial, Rho, PtotalIn=p, &
                  TeOut=Te, CvTotalOut=Cv, GammaOut=Gamma, zAverageOut=Zavg, &
                  HeatCond=HeatCond)

             ! Note that material index starts from 0
             if(Te > 0.0)then
                Value_V(Te_   +iMaterial*nThermo) = Te
                Value_V(Cv_   +iMaterial*nThermo) = Cv
                Value_V(Gamma_+iMaterial*nThermo) = Gamma
                Value_V(Zavg_ +iMaterial*nThermo) = Zavg
                Value_V(Cond_ +iMaterial*nThermo) = HeatCond
             else
                ! The eos() function returned impossible values, take ideal gas
                Value_V(Te_   +iMaterial*nThermo) = &
                     p/Rho*cProtonMass/cBoltzmann
                Value_V(Cv_   +iMaterial*nThermo) = 1.5*Rho
                Value_V(Gamma_+iMaterial*nThermo) = 5./3.
                Value_V(Zavg_ +iMaterial*nThermo) = 0.0
                Value_V(Cond_ +iMaterial*nThermo) = 0.0
             end if
          end if
       end do
    elseif(iTable == iTableOpacity)then
       ! Calculate gray specific opacities for all materials
       ! for given Rho and Te
       ! Au_ and Ay_ are optional
       Rho = Arg1
       Te  = Arg2
       do iMaterial = 0, nMaterial-1
          call eos(iMaterial, Rho, TeIn=Te, &
               OpacityPlanckOut_I=OpacityPlanck_W, &
               OpacityRosselandOut_I=OpacityRosseland_W)
          Value_V(1 + iMaterial*2) = OpacityPlanck_W(1)/Rho
          Value_V(2 + iMaterial*2) = OpacityRosseland_W(1)/Rho
       end do
    elseif(any(iTable == iTableOpacity_I(0:nMaterial-1)))then

       if(iProc == 0 .and. size(Value_V) /= 2*nWave)then
          write(*,*) 'ERROR ',NameSub, &
               ' number of table elements=', size(Value_V),&
               ' does not agree with 2*nWave=', 2*nWave
          call stop_mpi(NameSub//': Config.pl -nWave=..; make CRASH')
       end if
       ! Calculate multigroup specific opacities for one material
       ! for given Rho and Te
       Rho = Arg1
       Te  = Arg2
       do iMaterial = 0, nMaterial-1
          if(iTable == iTableOpacity_I(iMaterial)) EXIT
       end do

       if(UseAy .and. iMaterial == Ay_)then
          ! copy the opacity for acrylic from polyimide
          call interpolate_lookup_table(iTableOpacity_I(Plastic_), &
               Rho, Te, Value_V, DoExtrapolate = .false.)
       else
          call eos(iMaterial, Rho, TeIn=Te, &
               OpacityPlanckOut_I=OpacityPlanck_W, &
               OpacityRosselandOut_I=OpacityRosseland_W)
          Value_V(1:nWave)         = OpacityPlanck_W/Rho
          Value_V(nWave+1:2*nWave) = OpacityRosseland_W/Rho
       end if
    elseif(iTable == iTableSesame)then
       ! Calculate pressure and internal energy for all materials
       Rho = Arg1
       Te  = Arg2
       do iMaterial = 0, nMaterial-1
          call eos(iMaterial, Rho, TeIn=Te, pTotalOut=p, eTotalOut=e)
          Value_V(1 + iMaterial*2) = p
          Value_V(2 + iMaterial*2) = e/Rho
       end do
    else
       write(*,*)NameSub,' iTable=', iTable
       call stop_mpi(NameSub//' invalid value for iTable')
    endif

  end subroutine calc_table_value
  !===========================================================================

  subroutine user_material_properties(State_V, i, j, k, iBlock, iDir, &
       EinternalIn, TeIn, NatomicOut, AverageIonChargeOut, &
       EinternalOut, TeOut, PressureOut, &
       CvOut, GammaOut, HeatCondOut, IonHeatCondOut, TeTiRelaxOut, &
       OpacityPlanckOut_W, OpacityRosselandOut_W, PlanckOut_W)

    ! The State_V vector is in normalized units, all other physical
    ! quantities are in SI.
    !
    ! If the electron energy is used, then EinternalIn, EinternalOut,
    ! PressureOut, CvOut refer to the electron internal energies,
    ! electron pressure, and electron specific heat, respectively.
    ! Otherwise they refer to the total (electron + ion) internal energies,
    ! total (electron + ion) pressure, and the total specific heat.

    use CRASH_ModEos,  ONLY: eos, Xe_, Be_, Plastic_
    use CRASH_ModMultiGroup, ONLY: get_planck_g_from_temperature
    use ModMain,       ONLY: nI
    use ModAdvance,    ONLY: State_VGB, UseElectronPressure
    use ModPhysics,    ONLY: No2Si_V, UnitRho_, UnitP_, &
         inv_gm1
    use ModVarIndexes, ONLY: nVar, Rho_, p_, nWave, &
         Pe_, ExtraEint_
    use ModLookupTable,ONLY: interpolate_lookup_table
    use ModConst,      ONLY: cAtomicMass
    use CRASH_ModInterfaceNLTE, ONLY: NLTE_EOS

    real, intent(in) :: State_V(nVar)
    integer, optional, intent(in):: i, j, k, iBlock, iDir  ! cell/face index
    real, optional, intent(in)  :: EinternalIn             ! [J/m^3]
    real, optional, intent(in)  :: TeIn                    ! [K]
    real, optional, intent(out) :: NatomicOut              ! [1/m^3]
    real, optional, intent(out) :: AverageIonChargeOut     ! dimensionless
    real, optional, intent(out) :: EinternalOut            ! [J/m^3]
    real, optional, intent(out) :: TeOut                   ! [K]
    real, optional, intent(out) :: PressureOut             ! [Pa]
    real, optional, intent(out) :: CvOut                   ! [J/(K*m^3)]
    real, optional, intent(out) :: GammaOut                ! dimensionless
    real, optional, intent(out) :: HeatCondOut             ! [J/(m*K*s)]
    real, optional, intent(out) :: IonHeatCondOut          ! [J/(m*K*s)]
    real, optional, intent(out) :: TeTiRelaxOut            ! [1/s]
    real, optional, intent(out) :: &
         OpacityPlanckOut_W(nWave)                         ! [1/m]
    real, optional, intent(out) :: &
         OpacityRosselandOut_W(nWave)                      ! [1/m]

    ! Multi-group specific interface. The variables are respectively:
    !  Group Planckian spectral energy density
    real, optional, intent(out) :: PlanckOut_W(nWave)      ! [J/m^3]

    logical :: IsMix
    integer :: iMaterial, jMaterial
    real    :: pSi, RhoSi, TeSi, EinternalSi, LevelSum
    real    :: Value_V(nMaterial*nThermo)
    ! Our gray opacity table are for three materials
    real    :: Opacity_V(2*max(3,nMaterial))
    real    :: GroupOpacity_W(2*nWave)
    real, dimension(0:nMaterial-1) :: pPerRho_I, Weight_I
    real :: ePerRho
    real :: RhoToARatioSi_I(0:Plastic_) = 0.0
    real :: Level_I(3), LevelLeft, LevelRight

    ! multi-group variables
    integer :: iWave
    real :: PlanckSi

    character (len=*), parameter :: NameSub = 'user_material_properties'
    !-------------------------------------------------------------------------
    ! Density, transformed to SI
    RhoSi = No2Si_V(UnitRho_)*State_V(Rho_)

    if(present(EinternalIn)) EinternalSi = max(1e-30, EinternalIn)
    if(useNLTE.and..not.present(i))&
            call CON_stop('i,j,k,iBlock should be passed to user_material...for NLTE=T')

    ! The electron temperature may be needed for the opacities
    ! Initialize to negative value to see if it gets set
    TeSi = -7.70

    ! Find maximum level set value. 
    iMaterial   = maxloc(State_V(LevelXe_:LevelMax), 1) - 1

    ! By default use weight 1 for the material with maximum level
    IsMix               = .false.
    Weight_I            = 0.0
    Weight_I(iMaterial) = 1.0

    ! Calculate the weights
    if(UseVolumeFraction)then
       if(present(i) .and. .not. present(iDir))then
          ! This implementation is for 1D only !!!
          if(i>=0.and.i<=nI+1)then
             ! Divide by density to get actual levelset function
             Level_I = State_VGB(LevelXe_,i-1:i+1,j,k,iBlock) &
                  /    State_VGB(Rho_    ,i-1:i+1,j,k,iBlock)
             ! Calculate face values for the level set function
             LevelLeft  = 0.5*(Level_I(1)+Level_I(2))
             LevelRight = 0.5*(Level_I(2)+Level_I(3))
             ! Cell is mixed if face values change signs
             IsMix = LevelLeft*LevelRight < 0
             if(IsMix)then
                ! Make weight proportional to the cell volume fraction
                Weight_I(Xe_)      = max(LevelLeft,  LevelRight) &
                     /               abs(LevelLeft - LevelRight)
                Weight_I(Be_)      = 1 - Weight_I(Xe_)
                if(UsePl) Weight_I(Plastic_:nMaterial-1) = 0.0
             end if
          end if
       end if
    elseif(UseMixedCell)then
       ! Shall we use mixed material cells?
       LevelSum = sum(State_V(LevelXe_:LevelMax))
       IsMix = maxval(State_V(LevelXe_:LevelMax)) < MixLimit*LevelSum 

       if(IsMix)then
          ! Use number densities for eos() or weights in look up tables.
          if(UsePl) RhoToARatioSi_I = &
               State_V(LevelXe_:LevelXe_+Plastic_)*No2Si_V(UnitRho_)
          Weight_I = State_V(LevelXe_:LevelMax)/LevelSum
       end if
    end if

    ! get the atomic concentration
    if(present(NatomicOut))then
       if(IsMix)then
          NatomicOut = sum(RhoToARatioSi_I)/cAtomicMass
       else
          NatomicOut = RhoSi/(cAtomicMass*MassMaterial_I(iMaterial))
       end if
    end if

    if(UseElectronPressure)then
       call get_electron_thermo
    else
       call get_thermo
    end if

    if(present(TeOut)) TeOut = TeSi

    if(present(OpacityPlanckOut_W) .or. present(OpacityRosselandOut_W))then

       if(iTableOpacity > 0 .and. nWave == 1)then
          if(RhoSi <= 0 .or. TeSi <= 0) call lookup_error(&
               'Gray opacity(Rho,Te)', RhoSi, TeSi)

          call interpolate_lookup_table(iTableOpacity, RhoSi, TeSi, &
               Opacity_V, DoExtrapolate = .false.)

          Opacity_V(1:2*nMaterial:2) = Opacity_V(1:2*nMaterial:2) &
               *PlanckScaleFactor_I(0:nMaterial-1)
          Opacity_V(2:2*nMaterial:2) = Opacity_V(2:2*nMaterial:2) &
               *RosselandScaleFactor_I(0:nMaterial-1)
          if(UseVolumeFraction)then
             if(present(OpacityPlanckOut_W)) OpacityPlanckOut_W &
                  = sum(Weight_I*Opacity_V(1:2*nMaterial:2)) * RhoSi
             if(present(OpacityRosselandOut_W)) OpacityRosselandOut_W &
                  = sum(Weight_I*Opacity_V(2:2*nMaterial:2)) * RhoSi
          else
             if(present(OpacityPlanckOut_W)) OpacityPlanckOut_W &
                  = Opacity_V(2*iMaterial + 1) * RhoSi
             if(present(OpacityRosselandOut_W)) OpacityRosselandOut_W &
                  = Opacity_V(2*iMaterial + 2) * RhoSi
          end if

       elseif(all(iTableOpacity_I(0:nMaterial-1) > 0))then
          if(UseVolumeFraction)then
             if(present(OpacityPlanckOut_W)) OpacityPlanckOut_W = 0
             if(present(OpacityRosselandOut_W)) OpacityRosselandOut_W = 0
             do jMaterial = 0, nMaterial-1

                if(RhoSi <= 0 .or. TeSi <= 0) call lookup_error( &
                     'Group opacity(Rho,Te,jMaterial)', RhoSi, TeSi, jMaterial)

                call interpolate_lookup_table(iTableOpacity_I(jMaterial), &
                     RhoSi, TeSi, GroupOpacity_W, DoExtrapolate = .false.)

                GroupOpacity_W(1:nWave) = GroupOpacity_W(1:nWave) &
                     *PlanckScaleFactor_I(jMaterial)
                GroupOpacity_W(nWave+1:) = GroupOpacity_W(nWave+1:) &
                     *RosselandScaleFactor_I(jMaterial)

                if(present(OpacityPlanckOut_W)) &
                     OpacityPlanckOut_W = OpacityPlanckOut_W &
                     + Weight_I(jMaterial)*GroupOpacity_W(1:nWave) * RhoSi
                if(present(OpacityRosselandOut_W)) &
                     OpacityRosselandOut_W =  OpacityRosselandOut_W &
                     + Weight_I(jMaterial)*GroupOpacity_W(nWave+1:)*RhoSi
             end do
          else

             if(RhoSi <= 0 .or. TeSi <= 0) call lookup_error( &
                  'Group opacity(Rho,Te,iMaterial)', RhoSi, TeSi, iMaterial)

             call interpolate_lookup_table(iTableOpacity_I(iMaterial), &
                  RhoSi, TeSi, GroupOpacity_W, DoExtrapolate = .false.)

             GroupOpacity_W(1:nWave) = GroupOpacity_W(1:nWave) &
                  *PlanckScaleFactor_I(iMaterial)
             GroupOpacity_W(nWave+1:) = GroupOpacity_W(nWave+1:) &
                  *RosselandScaleFactor_I(iMaterial)

             if(present(OpacityPlanckOut_W)) OpacityPlanckOut_W &
                  = GroupOpacity_W(1:nWave)*RhoSi
             if(present(OpacityRosselandOut_W)) OpacityRosselandOut_W &
                  = GroupOpacity_W(nWave+1:)*RhoSi
          end if
       else
          ! inline opacities
          if(IsMix)then
             call eos(RhoToARatioSi_I, TeIn=TeSi, &
                  OpacityPlanckOut_I=OpacityPlanckOut_W, &
                  OpacityRosselandOut_I=OpacityRosselandOut_W)
          else
             if(.not.UseNLTE)then
                call eos(iMaterial, RhoSi, TeIn=TeSi, &
                  OpacityPlanckOut_I=OpacityPlanckOut_W, &
                  OpacityRosselandOut_I=OpacityRosselandOut_W)
             else
                call NLTE_EOS(iMaterial, RhoSi, TeIn=TeSi, &
                   EoBIn_I = EOverB_VGB(:,i,j,k,iBlock),   &
                   OpacityPlanckOut_I=OpacityPlanckOut_W, & 
                   OpacityRosselandOut_I=OpacityRosselandOut_W)
             end if

             if(present(OpacityPlanckOut_W)) OpacityPlanckOut_W &
                  = OpacityPlanckOut_W*PlanckScaleFactor_I(iMaterial)

             if(present(OpacityRosselandOut_W)) OpacityRosselandOut_W &
                  = OpacityRosselandOut_W*RosselandScaleFactor_I(iMaterial)

          end if

       end if
    end if

    if(present(PlanckOut_W))then
       do iWave = 1, nWave
          call get_planck_g_from_temperature(iWave, TeSi, EgSI=PlanckSi)
          PlanckOut_W(iWave) = PlanckSi
       end do
    end if

  contains

    !========================================================================

    subroutine get_thermo

      !----------------------------------------------------------------------

      ! Obtain the pressure from EinternalIn or TeIn or State_V
      ! Do this for various cases: mixed cell or not, lookup tables or not
      if(present(EinternalIn))then
         ! Obtain the pressure from EinternalIn
         if(iTablePPerRho > 0)then
            ! Use lookup table
            if(RhoSi <= 0 .or. EinternalSi <= 0) call lookup_error( &
                 'pPerRho(Rho,Einternal)', RhoSi, EinternalSi)

            call interpolate_lookup_table(iTablePPerRho, &
                 EinternalSi/RhoSi, RhoSi, pPerRho_I, DoExtrapolate = .false.)
            ! Use a number density weighted average
            pSi = RhoSi*sum(Weight_I*pPerRho_I)
         else
            ! Use EOS function
            if(IsMix)then
               call eos(RhoToARatioSi_I, eTotalIn=EinternalIn, &
                    pTotalOut=pSi, TeOut=TeSi, CvTotalOut=CvOut, &
                    GammaOut=GammaOut, zAverageOut=AverageIonChargeOut, &
                    HeatCond=HeatCondOut)
            else
               if(.not.UseNLTE)then
                  call eos(iMaterial, Rho=RhoSi, eTotalIn=EinternalIn, &
                    pTotalOut=pSi, TeOut=TeSi, CvTotalOut=CvOut, &
                    GammaOut=GammaOut, zAverageOut=AverageIonChargeOut, &
                    HeatCond=HeatCondOut)
               else
                   call NLTE_EOS(iMaterial, Rho=RhoSi, eTotalIn=EinternalIn, &
                    EoBIn_I = EOverB_VGB(:,i,j,k,iBlock),        &
                    pTotalOut=pSi, TeOut=TeSi, CvTotalOut=CvOut, &
                    GammaOut=GammaOut, zAverageOut=AverageIonChargeOut, &
                    HeatCond=HeatCondOut)
                end if
            end if
         end if
      elseif(present(TeIn))then
         ! Calculate pressure from electron temperature
         TeSi = TeIn
         if( IsMix ) then
            call eos(RhoToARatioSi_I, TeIn=TeIn, &
                 eTotalOut=EinternalOut, pTotalOut=pSi, &
                 CvTotalOut=CvOut, GammaOut=GammaOut, &
                 zAverageOut=AverageIonChargeOut, HeatCond=HeatCondOut)
         else
            if(.not.UseNLTE)then
               call eos(iMaterial, Rho=RhoSi, TeIn=TeIn, &
                 eTotalOut=EinternalOut, pTotalOut=pSi, &
                 CvTotalOut=CvOut, GammaOut=GammaOut, &
                 zAverageOut=AverageIonChargeOut, HeatCond=HeatCondOut)
            else
               call NLTE_EOS(iMaterial, Rho=RhoSi, TeIn=TeIn, &
                 EoBIn_I = EOverB_VGB(:,i,j,k,iBlock),        &
                 eTotalOut=EinternalOut, pTotalOut=pSi, &
                 CvTotalOut=CvOut, GammaOut=GammaOut, &
                 zAverageOut=AverageIonChargeOut, HeatCond=HeatCondOut)
            end if
         end if
      else
         ! Pressure is simply part of State_V
         pSi = State_V(p_)*No2Si_V(UnitP_)
         if(present(EinternalOut))then
            ! Obtain the internal energy from pressure
            ! Use the reversed lookup procedure
            if(iTablePPerRho > 0)then

               if(RhoSi <= 0 .or. pSi <= 0) call lookup_error( &
                    'PPerRho(p/Rho,Rho,E/Rho_out)', pSi, RhoSi)

               EinternalOut = 0.0
               do jMaterial = 0, nMaterial-1
                  if(Weight_I(jMaterial) < 1e-6) CYCLE
                  call interpolate_lookup_table(iTablePPerRho, &
                       iVal=jMaterial+1, &
                       ValIn=pSi/RhoSi, Arg2In=RhoSi, &
                       Value_V=pPerRho_I, Arg1Out=ePerRho, &
                       DoExtrapolate = .false.)
                  ! Use a number density weighted average
                  EinternalOut = EinternalOut &
                       + Weight_I(jMaterial)*ePerRho*RhoSi
               end do
            else
               if(IsMix)then
                  call eos(RhoToARatioSi_I, pTotalIn=pSi, &
                       EtotalOut=EinternalOut, TeOut=TeSi, &
                       CvTotalOut=CvOut, GammaOut=GammaOut, &
                       zAverageOut=AverageIonChargeOut, HeatCond=HeatCondOut)
               else
                  if(UseNLTE)then
                       ! call CON_stop('NLTE energy cannot be found from pressure')
                     EinternalOut = pSi*inv_gm1 + State_VGB(ExtraEint_,i,j,k,iBlock)*&
                          No2Si_V(UnitP_)
                     EInternalOut = max(EInternalOut,1e-30)
                     call NLTE_EOS(iMaterial, RhoSi, eTotalIn=EinternalOut, &
                       EoBIn_I = EOverB_VGB(:,i,j,k,iBlock),        &
                       TeOut=TeSi,        &
                       CvTotalOut=CvOut, GammaOut=GammaOut, &
                       zAverageOut=AverageIonChargeOut, HeatCond=HeatCondOut)
                  else
                     call eos(iMaterial,RhoSi,pTotalIn=pSi, &
                       EtotalOut=EinternalOut, TeOut=TeSi, &
                       CvTotalOut=CvOut, GammaOut=GammaOut, &
                       zAverageOut=AverageIonChargeOut, HeatCond=HeatCondOut)
                  end if
               end if
            end if
         end if
      end if

      if(present(PressureOut)) PressureOut = pSi

      if(present(TeOut) .or. present(CvOut) .or. present(GammaOut) .or. &
           present(AverageIonChargeOut) .or. present(HeatCondOut) .or. &
           present(OpacityPlanckOut_W) .or. &
           present(OpacityRosselandOut_W) .or. &
           present(PlanckOut_W))then

         if(iTableThermo > 0)then

            if(RhoSi <= 0 .or. pSi <= 0) call lookup_error( &
                 'thermo(Rho,p)', RhoSi, pSi)

            call interpolate_lookup_table(iTableThermo, RhoSi, pSi/RhoSi, &
                 Value_V, DoExtrapolate = .false.)

            if(UseVolumeFraction)then
               TeSi = sum(Weight_I*Value_V(Te_   :nMaterial*nThermo:nThermo))
               if(present(CvOut))  CvOut  &
                    = sum(Weight_I*Value_V(Cv_   :nMaterial*nThermo:nThermo))
               if(present(GammaOut)) GammaOut &
                    = sum(Weight_I*Value_V(Gamma_:nMaterial*nThermo:nThermo))
               if(present(AverageIonChargeOut)) AverageIonChargeOut &
                    = sum(Weight_I*Value_V(Zavg_ :nMaterial*nThermo:nThermo))
               if(present(HeatCondOut)) HeatCondOut &
                    = sum(Weight_I*Value_V(Cond_ :nMaterial*nThermo:nThermo))
            else
               TeSi = Value_V(Te_   +iMaterial*nThermo)
               if(present(CvOut))  CvOut &
                    = Value_V(Cv_   +iMaterial*nThermo)
               if(present(GammaOut)) GammaOut &
                    = Value_V(Gamma_+iMaterial*nThermo)
               if(present(AverageIonChargeOut)) AverageIonChargeOut &
                    = Value_V(Zavg_ +iMaterial*nThermo)
               if(present(HeatCondOut)) HeatCondOut &
                    = Value_V(Cond_ +iMaterial*nThermo)
            end if

         elseif(TeSi < 0.0) then
            ! If TeSi is not set yet then we need to calculate things here
            if(IsMix) then
               call eos(RhoToARatioSi_I, pTotalIn=pSi, &
                    eTotalOut = EinternalOut, TeOut=TeSi, &
                    CvTotalOut=CvOut, GammaOut=GammaOut, &
                    zAverageOut=AverageIonChargeOut, HeatCond=HeatCondOut)
            else
               if(.not.UseNLTE)then
                  call eos(iMaterial, RhoSi, pTotalIn=pSi, &
                    eTotalOut = EinternalOut, TeOut=TeSi, &
                    CvTotalOut=CvOut, GammaOut=GammaOut, &
                    zAverageOut=AverageIonChargeOut, HeatCond=HeatCondOut)
               else
                  EinternalSi = pSi*inv_gm1+&
                       State_VGB(ExtraEint_,i,j,k,iBlock)*No2Si_V(UnitP_)
                  EInternalSi = max(EInternalSi,1e-30)
                  call NLTE_EOS(iMaterial, RhoSi, eTotalIn=EinternalSi,&
                       EoBIn_I = EOverB_VGB(:,i,j,k,iBlock),        &
                       TeOut=TeSi, &
                       CvTotalOut=CvOut, GammaOut=GammaOut, &
                       zAverageOut=AverageIonChargeOut, HeatCond=HeatCondOut)  
               end if
            end if
         end if
      end if

    end subroutine get_thermo

    !========================================================================

    subroutine get_electron_thermo

      !----------------------------------------------------------------------

      ! Obtain the pressure from EinternalIn or TeIn or State_V
      ! Do this for mixed cell or not, lookup tables or not
      if(present(EinternalIn))then
         ! Obtain electron pressure from the true electron internal energy
         if(iTablePPerRho > 0)then
            ! Use lookup table
            if(RhoSi <= 0 .or. EinternalSi <= 0) call lookup_error( &
                 'pPerRho_e(Rho,Eint)', RhoSi, EinternalSi)

            call interpolate_lookup_table(iTablePPerRho, &
                 EinternalSi/RhoSi, RhoSi, pPerRho_I, DoExtrapolate = .false.)

            ! Use a number density weighted average
            pSi = RhoSi*sum(Weight_I*pPerRho_I)
         else
            ! Use inline electron EOS
            if(IsMix)then
               call eos(RhoToARatioSi_I, eElectronIn=EinternalIn, &
                    pElectronOut=pSi, TeOut=TeSi, CvElectronOut=CvOut, &
                    GammaEOut=GammaOut, zAverageOut=AverageIonChargeOut, &
                    HeatCond=HeatCondOut, TeTiRelax=TeTiRelaxOut)
            else
               if(.not.UseNLTE)then
                  call eos(iMaterial, Rho=RhoSi, eElectronIn=EinternalIn, &
                    pElectronOut=pSi, TeOut=TeSi, CvElectronOut=CvOut, &
                    GammaEOut=GammaOut, zAverageOut=AverageIonChargeOut, &
                    HeatCond=HeatCondOut, TeTiRelax=TeTiRelaxOut)
               else
                  call NLTE_EOS(iMaterial, Rho=RhoSi, eElectronIn=EinternalIn, &
                    EoBIn_I = EOverB_VGB(:,i,j,k,iBlock),         &
                    pElectronOut=pSi, TeOut=TeSi, CvElectronOut=CvOut, &
                    GammaEOut=GammaOut, zAverageOut=AverageIonChargeOut, &
                    HeatCond=HeatCondOut, TeTiRelax=TeTiRelaxOut)
               end if
            end if
         end if
      elseif(present(TeIn))then
         ! Calculate electron pressure from electron temperature
         TeSi = TeIn
         if(IsMix) then
            call eos(RhoToARatioSi_I, TeIn=TeIn, &
                 eElectronOut=EinternalOut, &
                 pElectronOut=pSi, CvElectronOut=CvOut, &
                 GammaEOut=GammaOut, zAverageOut=AverageIonChargeOut, &
                 HeatCond=HeatCondOut, TeTiRelax=TeTiRelaxOut)
         else
            if(.not.UseNLTE)then
               call eos(iMaterial, Rho=RhoSi, TeIn=TeIn, &
                 eElectronOut=EinternalOut, &
                 pElectronOut=pSi, CvElectronOut=CvOut, &
                 GammaEOut=GammaOut, zAverageOut=AverageIonChargeOut, &
                 HeatCond=HeatCondOut, TeTiRelax=TeTiRelaxOut)
            else
              call NLTE_EOS(iMaterial, Rho=RhoSi, TeIn=TeIn, &
                 EoBIn_I = EOverB_VGB(:,i,j,k,iBlock),         &
                 eElectronOut=EinternalOut, &
                 pElectronOut=pSi, CvElectronOut=CvOut, &
                 GammaEOut=GammaOut, zAverageOut=AverageIonChargeOut, &
                 HeatCond=HeatCondOut, TeTiRelax=TeTiRelaxOut)
            end if
         end if
      else
         ! electron pressure is State_V(Pe_)
         ! Use this pressure to calculate the true electron internal energy
         pSi = State_V(Pe_)*No2Si_V(UnitP_)
         if(present(EinternalOut))then
            if(iTablePPerRho > 0)then

               if(RhoSi <= 0 .or. pSi <= 0) call lookup_error( &
                    'PPerRho_e(p/Rho,Rho,E/Rho_out)', pSi, RhoSi)

               EinternalOut = 0.0
               do jMaterial = 0, nMaterial-1
                  if(Weight_I(jMaterial) < 1e-6) CYCLE
                  call interpolate_lookup_table(iTablePPerRho, &
                       iVal=jMaterial+1, &
                       ValIn=pSi/RhoSi, Arg2In=RhoSi, &
                       Value_V=pPerRho_I, Arg1Out=ePerRho, &
                       DoExtrapolate = .false.)
                  ! Use a number density weighted average
                  EinternalOut = EinternalOut &
                       + Weight_I(jMaterial)*ePerRho*RhoSi

               end do
            else
               if(IsMix)then
                  call eos(RhoToARatioSi_I, pElectronIn=pSi, &
                       eElectronOut=EinternalOut, TeOut=TeSi, &
                       CvElectronOut=CvOut, GammaEOut=GammaOut, &
                       zAverageOut=AverageIonChargeOut, HeatCond=HeatCondOut, &
                       TeTiRelax=TeTiRelaxOut)
               else
                  if(UseNLTE)then
                     !call CON_stop('NLTE energy cannot be found from pressure')
                     EInternalOut = pSi*inv_gm1+&
                          State_VGB(ExtraEint_,i,j,k,iBlock)*No2Si_V(UnitP_)
                     EInternalOut = max(EInternalOut,1e-30)
                     call NLTE_EOS(iMaterial, RhoSi, eElectronIn=EinternalOut,&
                          EoBIn_I = EOverB_VGB(:,i,j,k,iBlock),         &
                          TeOut=TeSi, &
                          CvElectronOut=CvOut, GammaEOut=GammaOut, &
                          zAverageOut=AverageIonChargeOut, HeatCond=HeatCondOut, &
                          TeTiRelax=TeTiRelaxOut)
                  else
                     call eos(iMaterial, RhoSi, pElectronIn=pSi, &
                       eElectronOut=EinternalOut, TeOut=TeSi, &
                       CvElectronOut=CvOut, GammaEOut=GammaOut, &
                       zAverageOut=AverageIonChargeOut, HeatCond=HeatCondOut, &
                       TeTiRelax=TeTiRelaxOut)
                  end if
               end if
            end if
         end if
      end if

      if(present(PressureOut)) PressureOut = pSi

      if(present(TeOut) .or. present(CvOut) .or. present(GammaOut) .or. &
           present(AverageIonChargeOut) .or. &
           present(HeatCondOut) .or. present(TeTiRelaxOut) .or. &
           present(OpacityPlanckOut_W) .or. &
           present(OpacityRosselandOut_W) .or. &
           present(PlanckOut_W))then

         if(iTableThermo > 0)then
            if(RhoSi <= 0 .or. pSi <= 0) call lookup_error( &
                 'thermo_e(Rho,p)', RhoSi, pSi)

            call interpolate_lookup_table(iTableThermo, RhoSi, pSi/RhoSi, &
                 Value_V, DoExtrapolate = .false.)

            if(UseVolumeFraction)then
               TeSi = sum(Weight_I*Value_V(Te_   :nMaterial*nThermo:nThermo))
               if(present(CvOut))  CvOut  &
                    = sum(Weight_I*Value_V(Cv_   :nMaterial*nThermo:nThermo))
               if(present(GammaOut)) GammaOut &
                    = sum(Weight_I*Value_V(Gamma_:nMaterial*nThermo:nThermo))
               if(present(AverageIonChargeOut)) AverageIonChargeOut &
                    = sum(Weight_I*Value_V(Zavg_ :nMaterial*nThermo:nThermo))
               if(present(HeatCondOut)) HeatCondOut &
                    = sum(Weight_I*Value_V(Cond_ :nMaterial*nThermo:nThermo))
               if(present(TeTiRelaxOut)) TeTiRelaxOut &
                    = sum(Weight_I*Value_V(TeTi_ :nMaterial*nThermo:nThermo))
            else
               TeSi = Value_V(Te_   +iMaterial*nThermo)
               if(present(CvOut)) CvOut &
                    = Value_V(Cv_   +iMaterial*nThermo)
               if(present(GammaOut)) GammaOut &
                    = Value_V(Gamma_+iMaterial*nThermo)
               if(present(AverageIonChargeOut)) AverageIonChargeOut &
                    = Value_V(Zavg_ +iMaterial*nThermo)
               if(present(HeatCondOut)) HeatCondOut &
                    = Value_V(Cond_ +iMaterial*nThermo)
               if(present(TeTiRelaxOut)) TeTiRelaxOut &
                    = Value_V(TeTi_ +iMaterial*nThermo)
            end if

         elseif(TeSi < 0.0) then
            ! If TeSi is not set yet then we need to calculate things here
            if(IsMix) then
               call eos(RhoToARatioSi_I, pElectronIn=pSi, &
                    eElectronOut=EinternalOut, TeOut=TeSi, &
                    CvElectronOut=CvOut, GammaEOut=GammaOut, &
                    zAverageOut=AverageIonChargeOut, &
                    HeatCond=HeatCondOut, TeTiRelax=TeTiRelaxOut)
            else
               if(.not.UseNLTE)then
                  call eos(iMaterial, RhoSi, pElectronIn=pSi, &
                    eElectronOut=EinternalOut, TeOut=TeSi, &
                    CvElectronOut=CvOut, GammaEOut=GammaOut, &
                    zAverageOut=AverageIonChargeOut, HeatCond=HeatCondOut, &
                    TeTiRelax=TeTiRelaxOut)
               else
                  EinternalSi = pSi*inv_gm1+&
                       State_VGB(ExtraEint_,i,j,k,iBlock)*No2Si_V(UnitP_)
                  EInternalSi = max(EInternalSi,1e-30)
                  call NLTE_EOS(iMaterial, RhoSi, eElectronIn=EinternalSi,&
                       EoBIn_I = EOverB_VGB(:,i,j,k,iBlock),         &
                       TeOut=TeSi, &
                       CvElectronOut=CvOut, GammaEOut=GammaOut, &
                       zAverageOut=AverageIonChargeOut, HeatCond=HeatCondOut, &
                       TeTiRelax=TeTiRelaxOut)
               end if
            end if
         end if
      end if

      if(UseGammaLaw .and. Gamma_I(iMaterial) > 1.0)then
         if(present(GammaOut)) GammaOut = Gamma_I(iMaterial)
         if(present(PressureOut) .and. present(EinternalIn))then
            pSi = EinternalIn*(Gamma_I(iMaterial) - 1)
            PressureOut = pSi
         elseif(present(EinternalOut))then
            EinternalOut = pSi/(Gamma_I(iMaterial) - 1)
         end if
         if(present(CvOut)) CvOut = pSi/TeSi/(Gamma_I(iMaterial)-1)
      end if

    end subroutine get_electron_thermo

    !========================================================================

    subroutine lookup_error(String, Arg1, Arg2, iArg)

      use ModProcMH, ONLY: iProc
      use ModGeometry, ONLY: Xyz_DGB
      use ModVarIndexes, ONLY: ExtraEint_

      character(len=*),  intent(in) :: String
      real,              intent(in) :: Arg1, Arg2
      integer, optional, intent(in) :: iArg

      !---------------------------------------------------------------------
      write(*,*) 'ERROR for lookup arguments of '//String//': ', Arg1, Arg2
      if(present(iArg)) write(*,*) 'iArg =', iArg

      write(*,*)'ERROR at i,j,k,iBlock,iProc=', i, j, k, iBlock, iProc
      write(*,*)'ERROR at x,y,z=', &
           Xyz_DGB(x_,i,j,k,iBlock), Xyz_DGB(y_,i,j,k,iBlock), Xyz_DGB(z_,i,j,k,iBlock)
      write(*,*)'ERROR pressure, ExtraEint=', State_V(p_)*No2Si_V(UnitP_), &
           State_V(ExtraEint_)*No2Si_V(UnitP_)
      write(*,*)'ERROR State_V=', State_V
      call stop_mpi('lookup_error')

    end subroutine lookup_error

  end subroutine user_material_properties

  !===========================================================================
  subroutine user_amr_criteria(iBlock, UserCriteria, TypeCriteria, IsFound)

    ! Set UserCriteria = 1.0 for refinement, 0.0 for coarsening.

    use BATL_lib,    ONLY: iNode_B, iTree_IA, Level_
    use ModSize,     ONLY: nI, nJ, nK
    use ModAdvance,  ONLY: State_VGB, Rho_
    use ModMain,     ONLY: UseLaserHeating
    use ModPhysics,  ONLY: Io2No_V, UnitRho_
    use ModGeometry, ONLY: Xyz_DGB, CellSize_DB, MinDxValue

    ! Variables required by this user subroutine
    integer, intent(in)          :: iBlock
    real, intent(out)            :: UserCriteria
    character (len=*),intent(in) :: TypeCriteria
    logical ,intent(inout)       :: IsFound

    logical:: IsXe_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    logical:: IsAu_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real   :: RhoMin
    integer:: i, j, k, iMin, iMax, jMin, jMax, kMin, kMax, nLevel, iNode
    !------------------------------------------------------------------
    IsFound = .true.

    UserCriteria = 0.0

    ! If block is beyond xMaxAmr, do not refine
    if(Xyz_DGB(x_,1,1,1,iBlock) >= xMaxAmr) RETURN

    iNode = iNode_B(iBlock)
    nLevel = iTree_IA(Level_,iNode)

    if( (CellSize_DB(x_,iBlock) - MinDxValue) > 1e-6)then
       iMin = 0; iMax = nI+1; jMin = 0; jMax = nJ+1; kMin = 0; kMax = nK+1
    else
       iMin = -1; iMax = nI+2; jMin = -1; jMax = nJ+2; kMin = -1; kMax = nK+2
    end if
    if(nJ == 1)then
       jMin = 1; jMax = 1
    endif
    if(nK == 1)then
       kMin = 1; kMax = 1
    end if

    ! If there is a Xe interface anywhere in the block, refine
    do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
       IsXe_G(i,j,k) = State_VGB(LevelXe_,i,j,k,iBlock) &
            >   maxval(State_VGB(LevelBe_:LevelMax,i,j,k,iBlock))
    end do; end do; end do

    UserCriteria = 1.0

    if(UseLaserHeating)then
       ! Refine all beryllium to the right of x = -5 micron
       if(any(maxloc(State_VGB(LevelXe_:LevelMax, &
            iMin:iMax,jMin:jMax,kMin:kMax,iBlock),1)-1 == Be_) &
            .and. Xyz_DGB(x_,nI,nJ,nK,iBlock)+0.5*CellSize_DB(x_,iBlock) > -5.0)then
          if(nLevel >= nLevelBeryllium) UserCriteria = 0.5
          RETURN
       end if
    end if

    if(any(IsXe_G(iMin:iMax,jMin:jMax,kMin:kMax)) .and. &
         .not. all(IsXe_G(iMin:iMax,jMin:jMax,kMin:kMax)))then
       if(nLevel >= nLevelShock) UserCriteria = 0.5
       RETURN
    end if

    if(UseAu)then
       ! If there is a Au interface anywhere in the block, refine
       do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
          if(UseAy)then
             IsAu_G(i,j,k) = State_VGB(LevelAu_,i,j,k,iBlock) &
                  > max(maxval(State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock)), &
                  State_VGB(LevelAy_,i,j,k,iBlock))
          else
             IsAu_G(i,j,k) = State_VGB(LevelAu_,i,j,k,iBlock) &
                  >   maxval(State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock))
          end if
       end do; end do; end do

       if(any(IsAu_G(iMin:iMax,jMin:jMax,kMin:kMax)) .and. &
            .not. all(IsAu_G(iMin:iMax,jMin:jMax,kMin:kMax)))then
          if(nLevel >= nLevelAuInterface) UserCriteria = 0.5
          RETURN
       end if
    end if

    ! If Xe density exceeds RhoMin, refine
    RhoMin = RhoMinAmrDim*Io2No_V(UnitRho_)
    do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
       if(IsXe_G(i,j,k) .and. State_VGB(Rho_,i,j,k,iBlock) > RhoMin)then
          if(nLevel >= nLevelShock) UserCriteria = 0.5
          RETURN
       end if
    end do; end do; end do

    ! No need to refine
    UserCriteria = 0.0

  end subroutine user_amr_criteria

  !===========================================================================

  subroutine set_yzr(x, y, z, r, rMax)

    real, intent(in)   :: x
    real, intent(inout):: y, z
    real, intent(out)  :: r

    real, intent(in), optional:: rMax

    ! For nozzle geometry divide the Y and Z coordinates by a shrink factor
    !
    !  = 1                            (for x <= xStartNozzle)
    !  = yRatioNozzle or zRatioNozzle (for x >= xEndNozzle)
    !
    ! The factor is linearly varying in the [xStartNozzle,xEndNozzle] interval.
    !
    ! Then calculate the 2D (x-r) coordinate pair from input coordinates
    !
    ! Then limit the radial distance if necessary

    real :: Factor
    real :: y0, z0, r0, rInner
    !-------------------------------------------------------------------------
    if(UseNozzle)then
       Factor = max(0.0, min(1.0, &
            (x - xStartNozzle)/(xEndNozzle - xStartNozzle)))

       ! Stretch ellipse into circle
       y0 = y
       z0 = z
       y  = y0/(1 + Factor*(yRatioNozzle-1))
       z  = z0/(1 + Factor*(zRatioNozzle-1))

       if(nK > 1 .and. zRatioNozzle > 0.99 .and. UseTubeRadius &
            .and. Factor > 0.0)then
          ! Keep tube width 
          ! Radial distance of stretched position
          r = sqrt(y**2 + z**2)
          if(r > rInnerTube)then
             r0     = sqrt(y0**2 + z0**2) ! original radial distance
             rInner = rInnerTube*r0/r     ! rInnerTube shrinks to rInner
             r = rInnerTube + r0 - rInner ! radius shifted instead of stretched
             y = r/r0 * y0                ! corresponding Y coordinate
             z = r/r0 * z0                ! corresponding Z coordinate
          end if
       end if
    end if

    if(nK > 1)then
       ! in 3D
       r = sqrt(y**2 + z**2)
    else
       ! in 2D
       r = abs(y)
       z = 0.0
    end if
    
    if(present(rMax))then

       if(r > rMax)then
          ! Shrink coordinates in the radial direction to rMax
          y = y*rMax/r
          z = z*rMax/r
          r = rMax
       end if

    end if
  end subroutine set_yzr

  !============================================================================

  subroutine user_set_resistivity(iBlock, Eta_G)
    ! This subrountine set the eta for every block
    use ModAdvance,    ONLY: State_VGB
    use ModConst,      ONLY: cBoltzmann, cElectronCharge, cMu
    use ModPhysics,    ONLY: Si2No_V, UnitX_, UnitT_
    use ModSize,       ONLY: nI, nJ, nK

    integer, intent(in) :: iBlock
    real,    intent(out):: Eta_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK) 

    integer :: i, j, k
    real :: TeSi, EtaCoef, HeatCondSi

    character (len=*), parameter :: NameSub = 'user_set_resistivity'
    !--------------------------------------------------------------------------

    EtaCoef = (cBoltzmann/cElectronCharge)**2 &
         *Si2No_V(UnitX_)**2/(cMu*Si2No_V(UnitT_))

    do k = MinK,MaxK; do j = MinJ,MaxJ; do i = MinI,MaxI
       call user_material_properties(State_VGB(:,i,j,k,iBlock), TeOut=TeSi, &
            HeatCondOut=HeatCondSi)

       if(HeatCondSi==0.0)then
          Eta_G(i,j,k) = MaxNormalizedResistivity
       else
          Eta_G(i,j,k) = min(EtaCoef*TeSi/HeatCondSi, MaxNormalizedResistivity)
       end if

    end do; end do; end do

  end subroutine user_set_resistivity

end module ModUser
