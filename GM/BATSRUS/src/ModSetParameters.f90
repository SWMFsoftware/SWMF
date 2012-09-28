!^CFG COPYRIGHT UM
subroutine MH_set_parameters(TypeAction)

  ! Set input parameters for the Global Magnetosphere (GM) module

  use ModProcMH
  use ModMain
  use ModAdvance
  use ModB0, ONLY: UseB0Source, UseCurlB0, DoUpdateB0, DtUpdateB0, &
       read_b0_param, init_mod_b0
  use ModGeometry, ONLY : init_mod_geometry, TypeGeometry, nMirror_D, &
       x1,x2,y1,y2,z1,z2,XyzMin_D,XyzMax_D,RadiusMin,RadiusMax,&
       MinBoundary,MaxBoundary,&
       read_gen_radial_grid, set_gen_radial_grid, NameGridFile
  use ModNodes, ONLY : init_mod_nodes
  use ModImplicit                                       !^CFG IF IMPLICIT
  use ModImplHypre, ONLY: hypre_read_param              !^CFG IF IMPLICIT
  use ModPhysics
  use ModProject                                        !^CFG IF PROJECTION
  use ModCT, ONLY : init_mod_ct, DoInitConstrainB       !^CFG IF CONSTRAINB
  use ModBlockData, ONLY: clean_block_data
  use BATL_lib, ONLY: read_amr_criteria, read_amr_geometry, &
       DoCritAmr, DoAutoAmr, DoStrictAmr, BetaProlong,&
       init_mpi, IsCartesianGrid, IsCartesian, IsRzGeometry, &
       IsCylindrical, IsRLonLat, IsLogRadius, IsGenRadius
  use ModAMR
  use ModParallel, ONLY : proc_dims
  use ModRaytrace                                       !^CFG IF RAYTRACE
  use ModIO
  use CON_planet,       ONLY: read_planet_var, check_planet_var, NamePlanet
  use ModPlanetConst
  use CON_axes,         ONLY: init_axes,get_axes
  use ModUtilities,     ONLY: check_dir, fix_dir_name, DoFlush, &
       split_string, join_string
  use CON_planet,       ONLY: get_planet
  use ModTimeConvert,   ONLY: time_int_to_real, time_real_to_int
  use ModReadParam
  use ModMessagePass,   ONLY: DoOneCoarserLayer
  use ModFaceValue,     ONLY: &
       UseTvdResChange, UseAccurateResChange, &
       UseVolumeIntegral4, UseFaceIntegral4, UseLimiter4, &
       DoLimitMomentum, BetaLimiter, TypeLimiter, read_face_value_param
  use ModPartSteady,    ONLY: UsePartSteady, MinCheckVar, MaxCheckVar, &
       RelativeEps_V, AbsoluteEps_V
  use ModUser,          ONLY: user_read_inputs, user_init_session, &
       NameUserModule, VersionUserModule
  use ModBoundaryCells, ONLY: init_mod_boundary_cells
  use ModPointImplicit, ONLY: read_point_implicit_param, UsePointImplicit
  use ModRestartFile,   ONLY: read_restart_parameters, init_mod_restart_file
  use ModHallResist,    ONLY: &
       UseHallResist, HallFactor, HallCmaxFactor, &
       PoleAngleHall, dPoleAngleHall, rInnerHall, DrInnerHall, &
       NameHallRegion, x0Hall, y0Hall, z0Hall, rSphereHall, DrSphereHall, &
       xSizeBoxHall, DxSizeBoxHall, &
       ySizeBoxHall, DySizeBoxHall, &
       zSizeBoxHall, DzSizeBoxHall, &
       UseBiermannBattery
  use ModHeatConduction, ONLY: read_heatconduction_param !^CFG IF IMPLICIT
  use ModRadDiffusion,   ONLY: read_rad_diffusion_param  !^CFG IF IMPLICIT
  use ModResistivity, ONLY: UseResistivity, &            !^CFG IF DISSFLUX
       read_resistivity_param, init_mod_resistivity      !^CFG IF DISSFLUX
  use ModMultiFluid, ONLY: MassIon_I,ChargeIon_I,nIonFluid, iFluid, &
       DoConserveNeutrals, DoOhNeutralBc, &
       uBcFactor_I, RhoBcFactor_I, RhoNeutralsISW_dim, &
       UxNeutralsISW_dim, UyNeutralsISW_dim, UzNeutralsISW_dim, &
       TNeutralsISW_dim, mProtonMass

  use ModMultiIon, ONLY: multi_ion_set_parameters
  use ModSolarwind, ONLY: UseSolarwindFile, NameSolarwindFile, &
       read_solar_wind_file, normalize_solar_wind_data
  use ModSatelliteFile, ONLY: nSatellite, &
       read_satellite_parameters, read_satellite_input_files
  use ModGroundMagPerturb
  use ModGmGeoindices, ONLY: &
       DoWriteIndices, DoCalcKp, nKpMins, dtWriteIndices, init_mod_geoindices
  use ModFaceFlux, ONLY: face_flux_set_parameters, TypeFluxNeutral, &
       UseClimit, UsePoleDiffusion
  use ModLookupTable, ONLY: read_lookup_table_param
  use ModIonoVelocity,ONLY: read_iono_velocity_param
  use ModTimeStepControl, ONLY: read_time_step_control_param
  use ModLaserHeating,    ONLY: read_laser_heating_param
  use ModIoUnit, ONLY: io_unit_new
  use ModNumConst, ONLY: cDegToRad

  !CORONA SPECIFIC PARAMETERS
  use EEE_ModMain, ONLY: EEE_set_parameters
  use ModMagnetogram, ONLY: set_parameters_magnetogram, &
       read_magnetogram_file, read_potential_field
  use ModExpansionFactors,ONLY: NameModelSW, CoronalT0Dim, read_wsa_coeff
  use ModCoronalHeating,  ONLY: read_corona_heating, &
       init_coronal_heating, UseCoronalHeating, DoOpenClosedHeat
  use ModRadiativeCooling,ONLY: UseRadCooling,&
       read_modified_cooling, check_cooling_param, read_chromosphere
  use ModWaves, ONLY: read_waves_param, check_waves
  use ModLdem, ONLY: UseLdem, NameLdemFile, iRadiusLdem, read_ldem
  use ModViscosity, ONLY: Viscosity_set_parameters, UseViscosity,&
                          Viscosity_init

  implicit none

  character (len=17) :: NameSub='MH_set_parameters'

  ! Arguments

  ! TypeAction determines if we read or check parameters
  character (len=*), intent(in) :: TypeAction

  ! Local variables
  integer :: iFile, i, j
  logical :: IsUninitialized      = .true.
  real :: local_root_dx

  !  logical :: HdfUninitialized      = .true.
  logical :: DoReadSolarwindFile  = .false.
  logical :: DoReadSatelliteFiles = .false.
  logical :: DoReadMagnetometerFile=.false.
  logical :: IsMirrorX,  IsMirrorY,  IsMirrorZ

  ! The name of the command
  character (len=lStringLine) :: NameCommand, StringLine

  ! Temporary variables
  logical :: DoEcho=.false.
  integer :: nVarRead=0, nGMin
  character (len=lStringLine) :: NameEquationRead="?"

  character (len=50) :: plot_string,log_string
  character (len=3)  :: plot_area, plot_var
  character (len=2)  :: NameCompRead="??"
  integer :: qtotal, nIJKRead_D(3)

  integer            :: TimingDepth=-1
  character (len=10) :: TimingStyle='cumu'

  ! Variables for checking/reading #STARTTIME command
  real (Real8_)         :: StartTimeCheck = -1.0_Real8_

  ! Variable for #UNIFORMAXIS
  logical:: UseUniformAxis = .true.

  ! Variables for the #GRIDRESOLUTION and #GRIDLEVEL commands
  real    :: InitialResolution = -1.0

  ! Variables for checking the user module
  character (len=lStringLine) :: NameUserModuleRead='?'
  real                        :: VersionUserModuleRead=0.0
  integer :: iSession, iPlotFile, iVar

  character(len=10) :: NamePrimitive_V(nVar)

  !-------------------------------------------------------------------------
  NameSub(1:2) = NameThisComp

  iSession = i_session_read()

  ! Initialize BATL
  call init_mpi(iComm)
 
  if(IsUninitialized)then
     call set_namevar
     call set_defaults
     IsUninitialized=.false.
  end if

  if(iSession>1)then
     restart=.false.           ! restart in session 1 only
  end if

  if(DoReadSatelliteFiles)then
     call read_satellite_input_files
     DoReadSatelliteFiles = .false.
  end if

  if(DoReadMagnetometerFile)then
     call read_mag_input_file
     DoReadMagnetometerFile = .false.
  end if

  select case(TypeAction)
  case('CHECK')
     if(iProc==0)write(*,*) NameSub,': CHECK iSession =',iSession

     if(StartTimeCheck > 0.0 .and. abs(StartTime - StartTimeCheck) > 0.001)then
        write(*,*)NameSub//' WARNING: '//NameThisComp//'::StartTimeCheck=', &
             StartTimeCheck,' differs from CON::StartTime=', &
             StartTime,' !!!'
        if(UseStrict)then
           call stop_mpi('Fix #STARTTIME/#SETREALTIME commands in PARAM.in')
        else
           ! Fix iStartTime_I array
           call time_real_to_int(StartTime, iStartTime_I)
        end if
     end if

     ! Planet NONE in GM means that we do not use a body
     if (NameThisComp=='GM' .and. NamePlanet == 'NONE' .and. iSession == 1)then
        body1 = .false.
        ! Change the default conservative criteria when there is no planet
        ! and the #CONSERVATIVECRITERIA command did not occur
        if(i_line_command("#CONSERVATIVECRITERIA") < 0)then
           nConservCrit = 0
           deallocate(TypeConservCrit_I)
           if(i_line_command("#NONCONSERVATIVE") < 0) &
                UseNonConservative = .false.
        end if
     end if

     ! In standalone mode set and obtain GM specific parameters 
     ! in CON_planet and CON_axes

     if(NameThisComp == 'GM') then
        ! Initialize axes
        call init_axes(StartTime)
        call get_axes(Time_Simulation, MagAxisTiltGsmOut = ThetaTilt)
        call get_planet(DipoleStrengthOut = DipoleStrengthSi)
     end if

     if(IsStandAlone .and. NameThisComp=='GM') then
        ! Check and set some planet variables (e.g. DoUpdateB0)
        call check_planet_var(iProc==0, time_accurate)

        if(body1)then
           call get_planet(UseRotationOut = UseRotatingBc)
        else
           UseRotatingBc = .false.
        end if

        ! Obtain some planet parameters
        if(DipoleStrengthSi == 0.0 .and. .not. UseUserB0 )then
           UseB0       = .false.
           UseB0Source = .false.
           UseCurlB0   = .false.
           DoUpdateB0  = .false.
           DtUpdateB0  = -1.0
        else
           call get_planet( &
                DoUpdateB0Out = DoUpdateB0, DtUpdateB0Out = DtUpdateB0)
        end if

     end if

     call correct_parameters

     ! initialize module variables
     call init_mod_advance
     call init_mod_geometry
     call init_mod_boundary_cells
     call init_mod_nodes
     if(UseB0)         call init_mod_b0
     if(UseRaytrace)   call init_mod_raytrace  !^CFG IF RAYTRACE
     if(UseConstrainB) call init_mod_ct        !^CFG IF CONSTRAINB
     if(UseImplicit.or.UseSemiImplicit) &      !^CFG IF IMPLICIT
          call init_mod_implicit               !^CFG IF IMPLICIT
     if (DoWriteIndices) call init_mod_geoindices

     ! clean dynamic storage
     call clean_block_data

     ! set physics uses dimensional solar wind data
     if (DoReadSolarwindFile) call read_solar_wind_file

     call set_physics_constants

     ! Normalization of solar wind data requires normalization in set_physics
     if (DoReadSolarwindFile) call normalize_solar_wind_data

     call set_extra_parameters

     ! initialize ModEqution (e.g. variable units)
     call init_mod_equation

     if(UseResistivity)call init_mod_resistivity !^CFG IF DISSFLUX
     if(UseViscosity) call Viscosity_init

     if(UseMagnetogram)then
        if(i_line_command("#MAGNETOGRAM") > 0)then
           call read_magnetogram_file(NamePlotDir)
        elseif(i_line_command("#READPOTENTIALFIELD") > 0)then
           call read_potential_field(NamePlotDir)
        end if
     end if


     if(UseEmpiricalSW .and. i_line_command("#EMPIRICALSW") > 0)&
          call set_empirical_model(NameModelSW, BodyTDim_I(IonFirst_))

     if(UseCoronalHeating)call init_coronal_heating
     call check_cooling_param

     ! Initialize user module and allow user to modify things
     if(UseUserInitSession)call user_init_session

     ! if using open closed heating initialize auxilary WSA grid
     if(DoOpenClosedHeat .and. i_line_command("#OPENCLOSEDHEAT")>0)&
          call set_empirical_model('WSA', CoronalT0Dim)

     call check_waves

     if((iProc==0 .or. UseTimingAll) .and. IsStandAlone)then
        call timing_active(UseTiming)
        if(iSession==1)then
           call timing_step(0)
           if(UseTimingAll)then
              iUnitTiming = io_unit_new()
              write(NameTimingFile,'(a,i6.6)') 'timing.pe', iProc
              open(iUnitTiming, FILE=NameTimingFile, STATUS='replace')
              call timing_iounit(iUnitTiming)
              call timing_comp_proc("  ",iProc)
           end if
        end if
        call timing_depth(TimingDepth)
        call timing_report_style(TimingStyle)
     end if

     RETURN
  case('read','Read','READ')
     if(iProc==0)then
        write(*,*) NameSub,': READ iSession =',iSession,&
             ' iLine=',i_line_read(),' nLine =',n_line_read()
     end if
  case default
     call stop_mpi(NameSub//': TypeAction='//TypeAction// &
          ' must be "CHECK" or "READ"!')
  end select

  ! Read solarwindfile if #SOLARWINDFILE command is in this session
  DoReadSolarwindFile = .false.

  !\
  ! Read parameters from the text
  !/
  READPARAM: do
     if(.not.read_line(StringLine) )then
        IsLastRead = .true.
        EXIT READPARAM
     end if

     if(.not.read_command(NameCommand)) CYCLE READPARAM
     !        write (*,*) 'SAVEPLOTNAME iFile=', iFile

     select case(NameCommand)

     case("#BATLPROLONG")
        call read_var('BetaProlong',BetaProlong )

     case("#COMPONENT")
        call read_var('NameComp', NameCompRead)
        if(NameThisComp /= NameCompRead)then
           NameThisComp = NameCompRead
           call set_defaults
        end if

     case("#DESCRIPTION")
        call check_stand_alone

     case("#BEGIN_COMP","#END_COMP")
        call check_stand_alone
        i = len_trim(NameCommand)
        if(StringLine(i+2:i+3) /= NameThisComp)&
             call stop_mpi(NameSub//' ERROR: the component is not '// &
             NameThisComp//'in '//trim(StringLine))

     case("#END")
        call check_stand_alone
        IslastRead=.true.
        EXIT READPARAM

     case("#RUN")
        call check_stand_alone
        IslastRead=.false.
        EXIT READPARAM

     case("#STOP")
        call check_stand_alone
        call read_var('MaxIteration',nIter)
        call read_var('tSimulationMax',t_max)

     case("#CPUTIMEMAX")
        call check_stand_alone
        call read_var('CpuTimeMax',cputime_max)

     case("#CHECKSTOPFILE")
        call check_stand_alone
        call read_var('DoCheckStopfile',check_stopfile)

     case("#PROGRESS")
        call check_stand_alone
        call read_var('DnProgressShort',dn_progress1)
        call read_var('DnProgressLong',dn_progress2)

     case("#TIMEACCURATE")
        call check_stand_alone
        call read_var('DoTimeAccurate',time_accurate)

     case("#ECHO")
        call check_stand_alone
        call read_var('DoEcho',DoEcho)
        if(iProc==0)call read_echo_set(DoEcho)

     case("#FLUSH")
        call read_var('DoFlush',DoFlush)

     case("#TEST")
        call read_var('StringTest',test_string)

     case("#TESTXYZ")
        coord_test=.true.
        UseTestCell=.true.
        call read_var('xTest',Xtest)
        call read_var('yTest',Ytest)
        call read_var('zTest',Ztest)

     case("#TESTIJK")
        coord_test=.false.
        UseTestCell=.true.
        call read_var('iTest',Itest)
        call read_var('jTest',Jtest)
        call read_var('kTest',Ktest)
        call read_var('iBlockTest',BLKtest)
        call read_var('iProcTest',PROCtest)

     case("#TESTVAR")
        call read_var('iVarTest ',VARtest)

     case("#TESTDIM")
        call read_var('iDimTest ',DIMtest)

     case("#TESTTIME")
        call read_var('nIterTest',ITERtest)
        call read_var('TimeTest',Ttest)

     case("#STRICT")
        call read_var('UseStrict',UseStrict)

     case("#VERBOSE")
        call read_var('lVerbose',lVerbose)

     case("#DEBUG")
        call read_var('DoDebug',okdebug)
        call read_var('DoDebugGhost',ShowGhostCells)

     case("#TIMING")
        call read_var('UseTiming',UseTiming)
        if(UseTiming)then
           call read_var('DnTiming',dn_timing)
           call read_var('nDepthTiming',TimingDepth)
           call read_var('TypeTimingReport',TimingStyle)
           UseTimingAll = index(TimingStyle,'all') > 0
           TimingStyle  = TimingStyle(1:4)
        end if

     case("#OPTIMIZEMPI")
        call read_var('UseOptimizeMpi', UseOptimizeMpi)

     case("#OUTERBOUNDARY")
        call read_var('TypeBcEast'  ,TypeBc_I(1))  
        call read_var('TypeBcWest'  ,TypeBc_I(2))
        call read_var('TypeBcSouth' ,TypeBc_I(3))
        call read_var('TypeBcNorth' ,TypeBc_I(4))
        call read_var('TypeBcBottom',TypeBc_I(5))
        call read_var('TypeBcTop'   ,TypeBc_I(6))    

     case("#INNERBOUNDARY")
        call read_var('TypeBcInner',TypeBc_I(body1_))
        ! Backward compatibility 
        if (TypeBc_I(body1_) == 'coronatoih') TypeBc_I(body1_)= 'buffergrid'
        if(UseBody2) call read_var('TypeBcBody2',TypeBc_I(body2_)) 

     case("#TIMESTEPPING", "#RUNGEKUTTA", "#RK")
        call read_var('nStage',  nStage)
        call read_var('CflExpl', Cfl)
        CflOrig = Cfl
        ExplCfl = Cfl                                   !^CFG IF IMPLICIT
        UseHalfStep = NameCommand == "#TIMESTEPPING" .and. nStage <= 2

     case("#FIXEDTIMESTEP")
        call read_var('UseDtFixed',UseDtFixed)
        if(UseDtFixed)call read_var('DtFixedDim', DtFixedDim)

     case("#PARTSTEADY")
        call read_var('UsePartSteady',UsePartSteady)

     case("#PARTSTEADYCRITERIA","#STEADYCRITERIA")
        call read_var('MinCheckVar',MinCheckVar)
        call read_var('MaxCheckVar',MaxCheckVar)
        do iVar=MinCheckVar, MaxCheckVar
           call read_var('RelativeEps',RelativeEps_V(iVar))
           call read_var('AbsoluteEps',AbsoluteEps_V(iVar))
        end do

     case("#POINTIMPLICIT")
        call read_point_implicit_param

     case("#IMPLICIT", &                                !^CFG IF IMPLICIT BEGIN
          "#IMPLCRITERIA", "#IMPLICITCRITERIA", "#STEPPINGCRITERIA", &
          "#PARTIMPL", "#PARTIMPLICIT",     &
          "#SEMIIMPL", "#SEMIIMPLICIT",     &
          "#IMPLSCHEME", "#IMPLICITSCHEME", &
          "#IMPLSTEP", "#IMPLICITSTEP",     &
          "#IMPLCHECK", "#IMPLICITCHECK",   &
          "#IMPLENERGY", "#IMPLICITENERGY", &
          "#NEWTON", "#JACOBIAN", "#PRECONDITIONER", &
          "#KRYLOV", "#KRYLOVSIZE")
        call read_implicit_param(NameCommand)           

     case("#HYPRE")
        call hypre_read_param                           !^CFG END IMPLICIT

     case("#RESISTIVITY", "#RESISTIVITYOPTIONS", &      !^CFG IF DISSFLUX BEGIN
          "#RESISTIVITYREGION", "#RESISTIVEREGION")
        call read_resistivity_param(NameCommand)
        !                                               ^CFG END DISSFLUX

     case("#VISCOSITY", "#VISCOSITYREGION")
        call Viscosity_set_parameters(NameCommand)

     case("#HALLRESISTIVITY")
        call read_var('UseHallResist',  UseHallResist)
        if(UseHallResist)then
           call read_var('HallFactor',  HallFactor)
           call read_var('HallCmaxFactor', HallCmaxFactor)
        end if

     case("#HALLREGION")
        call read_var('NameHallRegion', NameHallRegion)

        i = index(NameHallRegion, '0')
        if(i < 1 .and. NameHallRegion /= 'all')then
           call read_var("x0Hall", x0Hall)
           call read_var("y0Hall", y0Hall)
           call read_var("z0Hall", z0Hall)
        else
           x0Hall = 0.0; y0Hall = 0.0; z0Hall = 0.0
           if(i>1)NameHallRegion = &
                NameHallRegion(1:i-1)//NameHallRegion(i+1:len(NameHallRegion))
        end if

        select case(NameHallRegion)
        case("all", "user")
        case("sphere")
           call read_var("rSphereHall",rSphereHall)
           call read_var("DrSphereHall",DrSphereHall)
        case("box")
           call read_var("xSizeBoxHall ",xSizeBoxHall)
           call read_var("DxSizeBoxHall",DxSizeBoxHall)
           call read_var("ySizeBoxHall ",ySizeBoxHall)
           call read_var("DySizeBoxHall",DySizeBoxHall)
           call read_var("zSizeBoxHall ",zSizeBoxHall)
           call read_var("DzSizeBoxHall",DzSizeBoxHall)
        case default
           call stop_mpi(NameSub//': unknown NameHallRegion='&
                //NameHallRegion)
        end select

     case("#HALLPOLEREGION")
        call read_var("PoleAngleHall ", PoleAngleHall)
        call read_var("dPoleAngleHall", dPoleAngleHall)
        PoleAngleHall  =  PoleAngleHall*cDegToRad
        dPoleAngleHall = dPoleAngleHall*cDegToRad

     case("#HALLINNERREGION")
        call read_var("rInnerHall ", rInnerHall)
        call read_var("DrInnerHall", DrInnerHall)

     case("#BIERMANNBATTERY")
        call read_var("UseBiermannBattery", UseBiermannBattery)

     case("#MINIMUMPRESSURE")
        do iFluid = 1, nFluid
           call read_var('pMinDim', pMinDim_I(iFluid))
        end do

     case("#ELECTRONPRESSURE")
        call read_var('PeMinSi', PeMinSi)

     case("#ANISOTROPICPRESSURE")
        call read_var('UseConstantTau', UseConstantTau)
        if(UseConstantTau) &
             call read_var('TauInstabilitySi', TauInstabilitySi)
        call read_var('TauGlobalSi', TauGlobalSi)

     case("#EXTRAINTERNALENERGY")
        call read_var('ExtraEintMinSi', ExtraEintMinSi)

     case("#RADIATION", "#HEATFLUXLIMITER", "#ACCURATERADIATION")
        !^CFG IF IMPLICIT BEGIN 
        call read_rad_diffusion_param(NameCommand)

     case("#HEATCONDUCTION", "#WEAKFIELDCONDUCTION", &
          "#IONHEATCONDUCTION", "#HEATFLUXREGION")
        call read_heatconduction_param(NameCommand) !^CFG END IMPLICIT

     case("#SAVELOGFILE")
        call read_var('DoSaveLogfile',save_logfile)
        if(save_logfile)then
           if(iProc==0)call check_dir(NamePlotDir)
           nfile=max(nfile,logfile_)
           call read_var('StringLog',log_string)
           call read_var('DnSaveLogfile',dn_output(logfile_))
           call read_var('DtSaveLogfile',dt_output(logfile_))

           ! Log variables
           if(index(log_string,'VAR')>0 .or. index(log_string,'var')>0)then
              plot_dimensional(logfile_) = index(log_string,'VAR')>0
              log_time='step time'
              call read_var('log_vars',log_vars)
           elseif(index(log_string,'RAW')>0 .or. index(log_string,'raw')>0)then
              plot_dimensional(logfile_) = index(log_string,'RAW')>0
              log_time='step time'
              log_vars='dt '//NameConservativeVar//' Pmin Pmax'
           elseif(index(log_string,'MHD')>0 .or. index(log_string,'mhd')>0)then
              plot_dimensional(logfile_) = index(log_string,'MHD')>0
              log_time='step date time'
              log_vars=NameConservativeVar//' Pmin Pmax'
           elseif(index(log_string,'FLX')>0 .or. index(log_string,'flx')>0)then
              plot_dimensional(logfile_) = index(log_string,'FLX')>0
              log_time='step date time'
              log_vars='rho pmin pmax rhoflx pvecflx e2dflx'
           else
              call stop_mpi('Log variables (mhd,MHD,var,VAR,flx) missing'&
                   //' from log_string='//log_string)
           end if

           ! Determine the time output format to use in the logfile.
           ! This is loaded by default above, but can be input in the 
           ! log_string line.
           if(index(log_string,'none')>0) then
              log_time = 'none'
           elseif((index(log_string,'step')>0) .or. &
                (index(log_string,'date')>0) .or. &
                (index(log_string,'time')>0)) then
              log_time = ''
              if(index(log_string,'step')>0) log_time = 'step'
              if(index(log_string,'date')>0) &
                   write(log_time,'(a)') log_time(1:len_trim(log_time))&
                   //' date'
              if(index(log_string,'time')>0) &
                   write(log_time,'(a)') log_time(1:len_trim(log_time))&
                   //' time'
           end if

           ! If any flux variables are used - input a list of radii
           ! at which to calculate the flux
           if (index(log_vars,'flx')>0) &
                call read_var('log_R_str',log_R_str)
        end if

     case("#SAVEINITIAL")
        call read_var('DoSaveInitial',DoSaveInitial)

     case("#SAVEPLOT")
        call read_var('nPlotFile', nPlotFile)

        if(nPlotFile>0 .and. iProc==0)call check_dir(NamePlotDir)

        nfile = max(nfile, plot_ + nPlotFile)
        if (nFile > MaxFile .or. nPlotFile > MaxPlotFile) call stop_mpi(&
             'The number of ouput files is too large in #SAVEPLOT:'&
             //' nPlotFile > MaxPlotFile .or. nFile > MaxFile')

        do iFile = Plot_ + 1, Plot_ + nPlotFile
           call read_var('StringPlot',plot_string)

           ! Plotting frequency
           call read_var('DnSavePlot',dn_output(ifile))
           call read_var('DtSavePlot',dt_output(ifile))

           ! Plotting area
           if(index(plot_string,'cut')>0)then
              plot_area='cut'
              call read_var('xMinCut',plot_range(1,ifile))
              call read_var('xMaxCut',plot_range(2,ifile))
              call read_var('yMinCut',plot_range(3,ifile))
              call read_var('yMaxCut',plot_range(4,ifile))
              call read_var('zMinCut',plot_range(5,ifile))
              call read_var('zMaxCut',plot_range(6,ifile))
           elseif(index(plot_string,'slc')>0)then
              plot_area='slc'
              call read_var('xMinCut',plot_range(1,ifile))
              call read_var('xMaxCut',plot_range(2,ifile))
              call read_var('yMinCut',plot_range(3,ifile))
              call read_var('yMaxCut',plot_range(4,ifile))
              call read_var('zMinCut',plot_range(5,ifile))
              call read_var('zMaxCut',plot_range(6,ifile))
              call read_var('xPoint',plot_point(1,ifile))
              call read_var('yPoint',plot_point(2,ifile))
              call read_var('zPoint',plot_point(3,ifile))
              call read_var('xNormal',plot_normal(1,ifile))
              call read_var('yNormal',plot_normal(2,ifile))
              call read_var('zNormal',plot_normal(3,ifile))
           elseif(index(plot_string,'dpl')>0)then
              plot_area='dpl'
              call read_var('xMinCut',plot_range(1,ifile))
              call read_var('xMaxCut',plot_range(2,ifile))
              call read_var('yMinCut',plot_range(3,ifile))
              call read_var('yMaxCut',plot_range(4,ifile))
              call read_var('zMinCut',plot_range(5,ifile))
              call read_var('zMaxCut',plot_range(6,ifile))
           elseif (index(plot_string,'blk')>0) then
              plot_area='blk'
              call read_var('xPoint',plot_point(1,ifile))
              call read_var('yPoint',plot_point(2,ifile))
              call read_var('zPoint',plot_point(3,ifile))
           elseif(index(plot_string,'lin')>0)then     !^CFG IF RAYTRACE BEGIN
              iPlotFile = iFile - Plot_
              plot_area='lin'
              call read_var('NameLine', NameLine_I(iPlotFile), &
                   IsUpperCase=.true.)
              call read_var('IsSingleLine',IsSingleLine_I(iPlotFile))
              call read_var('nLine',nLine_I(iPlotFile))
              if(nLine_I(iPlotFile)==1)IsSingleLine_I(iPlotFile)=.true.
              if(nLine_I(iPlotFile) > MaxLine)then
                 if(iProc==0)then
                    write(*,*)NameSub,' WARNING: nLine=',nLine_I(iPlotFile),&
                         ' exceeds MaxLine=',MaxLine
                    write(*,*)NameSub,' WARNING reducing nLine to MaxLine'
                 end if
                 nLine_I(iPlotFile) = MaxLine
              end if
              do i = 1, nLine_I(iPlotFile)
                 call read_var('xStartLine',XyzStartLine_DII(1,i,iPlotFile))
                 call read_var('yStartLine',XyzStartLine_DII(2,i,iPlotFile))
                 call read_var('zStartLine',XyzStartLine_DII(3,i,iPlotFile))
                 call read_var('IsParallel',IsParallelLine_II(i,iPlotFile))
              end do                                  !^CFG END RAYTRACE
           elseif (index(plot_string,'eqr')>0)then
              plot_area='eqr'
              call read_var('nRadius',   plot_range(1,ifile))
              call read_var('nLon',      plot_range(2,ifile))
              call read_var('RadiusMin', plot_range(3,ifile))
              call read_var('RadiusMax', plot_range(4,ifile))
           elseif (index(plot_string,'ieb')>0)then
              plot_area='ieb'
           elseif (index(plot_string,'lcb')>0)then
              plot_area='lcb'
              call read_var('Radius', plot_range(1,ifile))
              call read_var('nLon',   plot_range(2,ifile))
           elseif (index(plot_string,'sph')>0)then
   	      plot_area='sph'
	      call read_var('Radius',plot_range(1,ifile))
           elseif (index(plot_string,'los')>0) then
              plot_area='los'
              ! Line of sight vector
              ! Satellite position
              call read_var('ObsPosX',ObsPos_DI(1,ifile))
              call read_var('ObsPosY',ObsPos_DI(2,ifile))
              call read_var('ObsPosZ',ObsPos_DI(3,ifile))
              ! Offset angle
              call read_var('OffsetAngle',offset_angle(ifile))
              offset_angle(ifile) = offset_angle(ifile)*cDegToRad
              ! read max dimensions of the 2d image plane
              call read_var('rSizeImage',r_size_image(ifile))
              ! read the position of image origin relative to grid origin
              call read_var('xOffset',xoffset(ifile))
              call read_var('yOffset',yoffset(ifile))
              ! read the occulting radius
              call read_var('rOccult',radius_occult(ifile))
              ! read the limb darkening parameter
              call read_var('MuLimbDarkening',mu_los)
              ! read the number of pixels
              call read_var('nPix',n_pix_r(ifile))
              ! if it is an EUV plot using a long table then read in the name 
              ! of the specific lookup table (will be matched to the name read
              ! in by the lookuptable command).
              if (index(plot_string,'TBL')>0&
                   .or.index(plot_string,'tbl')>0) &
                   call read_var('NameLosTable',NameLosTable(ifile))            
           elseif (index(plot_string,'rfr')>0) then
              ! Refractive radiowave image 
              plot_area='rfr'
              ! Observer position
              call read_var('ObsPosX', ObsPos_DI(1,ifile))
              call read_var('ObsPosY', ObsPos_DI(2,ifile))
              call read_var('ObsPosZ', ObsPos_DI(3,ifile))
              ! read number of radiowave frequencies, i.e. # of plots
              !call read_var('nRadioFrequency', nRadioFrequency)
              call read_var('StringRadioFrequency', &
                   StringRadioFrequency_I(iFile))
              call read_var('X_Size_Image', X_Size_Image(ifile))
              call read_var('Y_Size_Image', Y_Size_Image(ifile))
              ! read the number of pixels
              call read_var('n_Pix_X', n_Pix_X(ifile))            
              call read_var('n_Pix_Y', n_Pix_Y(ifile))            
           elseif (index(plot_string,'ion')>0) then
              plot_area='ion'
           elseif(index(plot_string,'1d')>0)then
              plot_area='1d_'
           elseif(index(plot_string,'2d')>0)then
              plot_area='2d_'
           elseif(index(plot_string,'3d')>0)then
              plot_area='3d_'
           elseif(index(plot_string,'x=0') > 0)then
              plot_area = 'x=0'
           elseif(index(plot_string,'y=0') > 0)then
              plot_area = 'y=0'
           elseif(index(plot_string,'z=0') > 0)then
              plot_area = 'z=0'
           else
              call stop_mpi('Area (1d,2d,3d,x=0,y=0,z=0,cut,sph...) missing'&
                   //' from plot_string='//plot_string)
           end if

           ! Plot file format
           if(index(plot_string,'idl') >0 )then
              plot_form(ifile)='idl'
              if ((plot_area /= 'ion')&
                   .and. plot_area /= 'sph' &
                   .and. plot_area /= 'los' &
                   .and. plot_area /= 'rfr' &
                   .and. plot_area /= 'lin' &        !^CFG IF RAYTRACE
                   .and. plot_area /= 'eqr' &        !^CFG IF RAYTRACE
                   ) call read_var('DxSavePlot',plot_dx(1,ifile))
              if(.not.IsCartesianGrid)plot_dx(1,ifile)=-1.0 

              ! Extract the type of idl plot file: default is real4
              TypeIdlFile_I(iFile) = 'real4' 
              if(index(plot_string,'idl_real8') > 0) &
                   TypeIdlFile_I(iFile) = 'real8'
              if(index(plot_string,'idl_ascii') > 0) &
                   TypeIdlFile_I(iFile) = 'ascii'
           elseif(index(plot_string, 'hdf') > 0) then
              ! With these values VisIt recognises the files as timesteps
              ! with the general defaults it does not. 
              IsPlotName_n = .true.
              IsPlotName_t = .false.
              IsPlotName_e = .false.
              plot_form(iFile)='hdf'
           elseif(index(plot_string,'tec')>0)then 
              plot_form(ifile)='tec'
              plot_dx(1,ifile)=0.
           else
              call stop_mpi('Format (idl,tec) missing from plot_string='&
                   //plot_string)
           end if

           ! Plot variables
           if(index(plot_string,'VAR')>0 .or. index(plot_string,'var')>0 )then
              plot_var='var'
              plot_dimensional(ifile) = index(plot_string,'VAR')>0
              call read_var('NameVars',plot_vars(ifile))
              call read_var('NamePars',plot_pars(ifile))
              !                                         ^CFG  IF RAYTRACE BEGIN
           elseif(index(plot_string,'RAY')>0.or.index(plot_string,'ray')>0)then
              plot_var='ray'
              plot_dimensional(ifile) = index(plot_string,'RAY')>0
              plot_vars(ifile)='bx by bz theta1 phi1 theta2 phi2 status blk'
              plot_pars(ifile)='rbody'
              !                                         ^CFG END RAYTRACE
           elseif(index(plot_string,'RAW')>0.or.index(plot_string,'raw')>0)then
              plot_var='raw'
              plot_dimensional(ifile)=index(plot_string,'RAW')>0
              plot_vars(ifile)=NameConservativeVar//' p b1x b1y b1z absdivB'
              plot_pars(ifile)='g rbody'
           elseif(index(plot_string,'MHD')>0.or.index(plot_string,'mhd')>0)then
              plot_var='mhd'
              plot_dimensional(ifile) = index(plot_string,'MHD')>0
              plot_vars(ifile) = NamePrimitiveVar//' jx jy jz'
              plot_pars(ifile)='g rbody'
           elseif(index(plot_string,'HD')>0.or.index(plot_string,'hd')>0)then
              plot_var='hd'
              plot_dimensional(ifile) = index(plot_string,'HD')>0
              plot_vars(ifile) = NamePrimitiveVar
              plot_pars(ifile)='g'
              if(rBody>0.0)plot_pars(ifile)='g rbody'
           elseif(index(plot_string,'ALL')>0.or.index(plot_string,'all')>0)then
              ! This is intended for restart with a different dimensionality
              plot_var='all'
              plot_dimensional(ifile) = .false.
              call join_string(nVar, NameVar_V(1:nVar), plot_vars(ifile))
              plot_pars(ifile)='g'
           elseif(index(plot_string,'FUL')>0.or.index(plot_string,'ful')>0)then
              plot_var='ful'
              plot_dimensional(ifile) = index(plot_string,'FUL')>0
              plot_vars(ifile) = NamePrimitiveVar//' b1x b1y b1z e jx jy jz'
              plot_pars(ifile)='g rbody'
           elseif(index(plot_string,'FLX')>0.or.index(plot_string,'flx')>0)then
              plot_var='flx'
              plot_dimensional(ifile) = index(plot_string,'FLX')>0
              plot_vars(ifile)='rho mr br p jr pvecr'
              plot_pars(ifile)='g rbody'
           elseif(index(plot_string,'SOL')>0.or.index(plot_string,'sol')>0)then
              plot_var='sol'
              plot_dimensional(ifile) = index(plot_string,'SOL')>0
              plot_vars(ifile)='wl pb' ! white light
              plot_pars(ifile)='mu'
           elseif(index(plot_string,'EUV')>0.or.index(plot_string,'euv')>0)then
              plot_var='euv'
              plot_dimensional(ifile) = index(plot_string,'EUV')>0
              plot_vars(ifile)='euv171 euv195 euv284' ! main euv bands
              plot_pars(ifile)='mu'
           elseif(index(plot_string,'SXR')>0.or.index(plot_string,'sxr')>0)then
              plot_var='sxr'
              plot_dimensional(ifile) = index(plot_string,'SXR')>0
              plot_vars(ifile)='sxr' ! soft x-ray band
              plot_pars(ifile)='mu'
           elseif(index(plot_string,'TBL')>0.or.index(plot_string,'tbl')>0)then
              plot_var='tbl'
              plot_dimensional(ifile) = index(plot_string,'TBL')>0
              plot_vars(ifile)='tbl' ! will read a table in write_plot_los 
              plot_pars(ifile)='mu'
           elseif(index(plot_string,'RWI')>0.or.index(plot_string,'rwi')>0)then 
              plot_var='rwi'
              plot_dimensional(ifile) = .false.
              plot_vars(ifile)='' ! Intensity
              plot_pars(ifile)=''
           elseif(index(plot_string,'pos')>0.or.index(plot_string,'POS')>0)then
              plot_var='pos'
              plot_dimensional(ifile) = index(plot_string,'POS')>0
              if(plot_area /= 'lin')call stop_mpi(&
                   'Variable "pos" can only be used with area "lin" !')
           elseif(index(plot_string,'eqr')>0)then
              plot_var ='eqr'
              plot_dimensional(ifile) = .true.
           elseif(index(plot_string,'NUL')>0.or.index(plot_string,'nul')>0)then
              plot_var ='nul'
              plot_dimensional(ifile) = .true.
              plot_vars(ifile)=''
              plot_pars(ifile)=''
           elseif(index(plot_string,'INT')>0.or.index(plot_string,'int')>0)then
              plot_var ='int'
              plot_dimensional(ifile) = index(plot_string,'INT')>0
              plot_vars(ifile)=''
              plot_pars(ifile)=''
           elseif(index(plot_string,'BBK')>0.or.index(plot_string,'bbk')>0)then
              plot_var='blk'
              plot_dimensional(ifile) = index(plot_string,'BBK')>0
              plot_vars(ifile)='dx pe blk blkall'
           else
              call stop_mpi('Variable definition missing from plot_string=' &
                   //plot_string)
           end if

           ! Set equation parameters for 3D unstructured IDL files to describe
           ! block structure and the dipole. Needed by CCMC.
           if(plot_area == '3d_' .and. plot_form(iFile) == 'idl' &
                .and. plot_dx(1, iFile) < 0.0) &
                plot_pars(iFile) = 'g c th p1 p2 p3 NX NY NZ R'

           plot_type(iFile) = plot_area//'_'//plot_var
        end do

     case("#SAVEPLOTNAME")
        call read_var('IsPlotName_n',IsPlotName_n)
        call read_var('IsPlotName_t',IsPlotName_t)
        call read_var('IsPlotName_e',IsPlotName_e)
        ! Will set only _n true when not time accurate automatically.
        ! Set _n true if both _t and _e are false.
        if(.not.IsPlotName_t .and. .not.IsPlotName_e)&
             IsPlotName_n=.true.

     case("#PLOTFILENAME")
        call read_var('NameMaxTimeUnit', NameMaxTimeUnit)

     case("#SAVELOGNAME")
        call read_var('IsLogName_n',IsLogName_n)
        call read_var('IsLogName_e',IsLogName_e)
        ! Set _n true if _e is false.
        if(.not.IsLogName_e) IsLogName_n=.true.

     case("#SAVEPLOTSAMR")
        call read_var('DoSavePlotsAmr',save_plots_amr)

     case("#SAVEBINARY")
        call read_var('DoSaveBinary',save_binary)

     case("#GRIDRESOLUTION","#GRIDLEVEL","#AREARESOLUTION","#AREALEVEL",&
          "#AMRREGION")

        call read_amr_geometry(NameCommand,UseStrictIn=UseStrict, &
             InitLevelInOut=initial_refine_levels, &
             InitResInOut=InitialResolution)

     case("#AMRLEVELS")
        call read_var('MinBlockLevel',min_block_level)
        call read_var('MaxBlockLevel',max_block_level)
        if(iSession==1)then
           DoSetLevels=.true.
        else
           call set_levels
        end if

     case("#AMRRESOLUTION")
        call read_var('DxCellMin',min_cell_dx)
        call read_var('DxCellMax',max_cell_dx)
        ! See also end of correct_parameters

     case("#DOAMRPROFILE")
        call read_var('DoAmrPofile',DoProfileAmr)

     case("#AMR")
        call read_var('DnRefine',DnAmr)
        DoAmr = DnAmr > 0
        DtAmr = -1.0
        if (DoAmr)then
           call read_var('DoAutoRefine',automatic_refinement)
           if (automatic_refinement) call read_amr_criteria("#AMR")
        end if

     case("#AMRLIMIT", "#AMRTYPE")
        call read_amr_criteria(NameCommand)

     case("#DOAMR")
        call read_var('DoAmr',DoAmr)
        if(DoAmr) then
           call read_var('DnAmr',DnAmr)
           call read_var('DtAmr',DtAmr)
           call read_var('IsStrictAmr'  ,DoStrictAmr)
        end if

     case("#AMRCRITERIA")

        DoCritAmr = .true.
        DoAutoAmr = .true.
        automatic_refinement = DoAutoAmr ! for now
        call read_var('nCriteria',nAmrCriteria)
        call init_mod_amr(nAmrCriteria)
        call read_amr_criteria(NameCommand, &
             nCritInOut=nAmrCriteria, NameCritOut_I=RefineCrit,&
             NameStatVarIn_V= NameVar_V,&
             nStateVarIn = nVar,ReadExtraOut=UseSunEarth)
        if (UseSunEarth) then
           call read_var('xEarth'  ,xEarth)
           call read_var('yEarth'  ,yEarth)
           call read_var('zEarth'  ,zEarth)
           call read_var('InvD2Ray',InvD2Ray)
        end if
        
     case("#AMRCRITERIALEVEL","#AMRCRITERIARESOLUTION")

        DoCritAmr = .true.
        DoAutoAmr = .true.
        automatic_refinement = DoAutoAmr ! for now
        call read_var('nCriteria',nAmrCriteria)
        call init_mod_amr(nAmrCriteria)
        call read_amr_criteria(NameCommand, &
             nCritInOut=nAmrCriteria, NameCritOut_I=RefineCrit,&
             NameStatVarIn_V= NameVar_V,&
             nStateVarIn = nVar,ReadExtraOut=UseSunEarth)
        if(nAmrCriteria<0)call stop_mpi(NameSub// &
             ' ERROR: nAmrCriteria must be positiv.')
        if (UseSunEarth) then
           call read_var('xEarth'  ,xEarth)
           call read_var('yEarth'  ,yEarth)
           call read_var('zEarth'  ,zEarth)
           call read_var('InvD2Ray',InvD2Ray)
        end if

     case("#SCHEME")
        call read_var('nOrder'  ,nOrder)
        ! Set default value for nStage. Can be overwritten if desired.
        nStage = nOrder

        call read_var('TypeFlux',FluxType, IsUpperCase=.true.)
        BetaLimiter = 1.0
        if(nOrder > 1 .and. FluxType /= "SIMPLE")then
           call read_var('TypeLimiter', TypeLimiter)
           if(TypeLimiter /= 'minmod') &
                call read_var('LimiterBeta', BetaLimiter)
        else
           TypeLimiter = "no"
        end if

        nGMin = 2 ! Minimum number of ghost cell layers
        if(nOrder==4)then
           call read_var('UseVolumeIntegral4', UseVolumeIntegral4)
           if(UseVolumeIntegral4) nGMin = nGMin + 1
           call read_var('UseFaceIntegral4', UseFaceIntegral4)
           if(nDim == 1) UseFaceIntegral4 = .false.
           if(TypeLimiter /= 'no')then
              nGMin = nGMin + 1
              call read_var('UseLimiter4', UseLimiter4)
              if(UseLimiter4) nGMin = nGMin + 1
           end if
        end if
        if(nGMin /= nG .and. iProc==0)then
           write(*,*)'The code is configured with nG=',nG,' ghost cell layers.'
           write(*,*)'The selected scheme requires nGMin=',nGMin,' layers!'
           if(nGMin > nG)then
              write(*,*)'Either change settings or reconfigure and recompile!'
              call CON_stop(NameSub//': insufficient number of ghost cells')
           end if
        end if

     case('#LIMITER', '#RESCHANGE', '#RESOLUTIONCHANGE', '#TVDRESCHANGE', &
          '#LIMITPTOTAL')
        call read_face_value_param(NameCommand)

     case("#NONCONSERVATIVE")
        call read_var('UseNonConservative',UseNonConservative)

     case("#CONSERVATIVECRITERIA")
        call read_var('nConservCrit',nConservCrit)
        if(nConservCrit > 0) then
           if(allocated(TypeConservCrit_I)) deallocate(TypeConservCrit_I)
           allocate( TypeConservCrit_I(nConservCrit) )
           do i=1,nConservCrit
              call read_var('TypeConservCrit',TypeConservCrit_I(i),&
                   IsLowerCase=.true.)
              select case(TypeConservCrit_I(i))
                 !\
                 ! Geometry based criteria: 
                 !/
              case('r','radius')
                 !    non-conservative scheme is used for r < rConserv
                 TypeConservCrit_I(i) = 'r'
                 call read_var('rConserv',rConserv)
              case('parabola','paraboloid')
                 !    non-conservative scheme is used for 
                 !    x < xParabolaConserv - (y**2+z**2)/yParabolaConserv
                 TypeConservCrit_I(i) = 'parabola'
                 call read_var('xParabolaConserv',xParabolaConserv)
                 call read_var('yParabolaConserv',yParabolaConserv)
                 !\
                 ! Physics based criteria
                 !/
              case('p')
                 ! Balsara/Ryu switch 1
                 call read_var('pCoeffConserv',pCoeffConserv)
              case('gradp','jumpp')
                 ! Balsara/Ryu switch 2
                 call read_var('GradPCoeffConserv',GradPCoeffConserv)
              case default
                 if(UseStrict)then
                    call stop_mpi(NameSub//&
                         ' ERROR: unknown TypeConservCrit_I=' &
                         //TypeConservCrit_I(i))
                 else
                    if(iProc==0)write(*,'(a)') NameSub // &
                         ' WARNING: ignoring unknown TypeConservCrit_I=',&
                         trim(TypeConservCrit_I(i))//' !!!'
                 end if
              end select
              ! Check if the same criterion has been used before
              if(any(TypeConservCrit_I(1:i-1)==TypeConservCrit_I(i)))then
                 if(iProc==0)write(*,'(a)')NameSub // &
                      ' WARNING: multiple use of criterion ',&
                      trim(TypeConservCrit_I(i))
              end if
           end do
        end if

     case("#TIMESTEPCONTROL", "#CONTROLTIMESTEP", &
          "#CONTROLDECREASE", "#CONTROLINCREASE", &
          "#CONTROLFACTOR", "#CONTROLVAR", "#CONTROLINIT")
        call read_time_step_control_param(NameCommand)

     case("#UPDATECHECK")
        call read_var("UseUpdateCheck",UseUpdateCheck)          
        if(UseUpdateCheck)then
           call read_var("RhoMinPercent", percent_max_rho(1))        
           call read_var("RhoMaxPercent", percent_max_rho(2))
           call read_var("pMinPercent",   percent_max_p(1))
           call read_var("pMaxPercent",   percent_max_p(2))
        end if

     case("#PROLONGATION")
        call read_var('nOrderProlong',prolong_order)
        call read_var('TypeProlong' ,prolong_type)

     case("#MESSAGEPASS","#OPTIMIZE")               
        call read_var('TypeMessagePass', optimize_message_pass)

     case('#CLIMIT')
        call face_flux_set_parameters(NameCommand)

     case("#BORIS")                                  !^CFG IF BORISCORR BEGIN
        if(.not.UseB)CYCLE READPARAM
        call read_var('UseBorisCorrection', boris_correction)   
        if(boris_correction) then
           call read_var('BorisClightFactor', boris_cLight_factor)
           if(IsMhd)then
              UseBorisSimple = .false.
           else
              ! For non-MHD equations only simplified Boris correction is possible
              UseBorisSimple   = .true.
              boris_correction = .false.
           end if
        else
           boris_cLIGHT_factor = 1.0
        end if                                       !^CFG END BORISCORR

     case("#BORISSIMPLE")                            !^CFG IF SIMPLEBORIS BEGIN
        call read_var('UseBorisSimple',UseBorisSimple)
        if(UseBorisSimple) then
           call read_var('BorisClightFactor',boris_cLIGHT_factor)
           boris_correction=.false.                     !^CFG IF BORISCORR
        else
           boris_cLIGHT_factor = 1.0
        end if                                       !^CFG END SIMPLEBORIS

     case("#DIVB")
        if(.not.UseB)CYCLE READPARAM
        call read_var('UseDivbSource'   ,UseDivbSource)   
        call read_var('UseDivbDiffusion',UseDivbDiffusion)!^CFG IF DIVBDIFFUSE
        call read_var('UseProjection'   ,UseProjection)  !^CFG IF PROJECTION
        call read_var('UseConstrainB'   ,UseConstrainB)  !^CFG IF CONSTRAINB

        !^CFG IF CONSTRAINB BEGIN
        if (UseProjection.and.UseConstrainB.and.iProc==0) &
             call stop_mpi('Do not use projection and constrain B together!')
        !^CFG END CONSTRAINB
        !^CFG IF DIVBDIFFUSE BEGIN
        !^CFG IF PROJECTION BEGIN
        if (UseProjection.and.UseDivbSource.and.iProc==0) then
           write(*,'(a)')NameSub // &
                ' WARNING: using divbsource and projection together !!!'
           if (UseStrict) call stop_mpi('Correct ')
        end if
        !^CFG END PROJECTION
        !^CFG END DIVBDIFFUSE
        !^CFG IF CONSTRAINB BEGIN
        if (UseConstrainB.and.UseDivbSource.and.iProc==0) then
           write(*,'(a)')NameSub // &
                ' WARNING: using divbsource and constrain B together !!!'
           if (UseStrict) call stop_mpi('Correct ')
        end if
        !^CFG END CONSTRAINB

        if (iProc==0 &
             .and..not.UseHyperbolicDivb &
             .and..not.UseDivbSource &
             .and..not.UseProjection &      !^CFG IF PROJECTION
             .and..not.UseConstrainB &      !^CFG IF CONSTRAINB
             .and..not.UseDivbDiffusion &   !^CFG IF DIVBDIFFUSE
             ) then
           write(*,*) NameSub // &
                'WARNING: you should use some div B control method !!!'
           if (UseStrict) call stop_mpi('Correct PARAM.in!')
        end if
        ! Make sure that divbmax will be calculated      !^CFG IF PROJECTION
        DivbMax = -1.0                                   !^CFG IF PROJECTION
        ! reinitialize constrained transport if needed   !^CFG IF CONSTRAINB
        DoInitConstrainB = .true.                        !^CFG IF CONSTRAINB

     case("#USEB0", "#DIVBSOURCE", "#USECURLB0")
        if(.not.is_first_session())CYCLE READPARAM
        call read_b0_param(NameCommand)

     case("#HYPERBOLICDIVB")
        if(.not.UseB)CYCLE READPARAM
        if(Hyp_ == 1)then
           if(iProc==0)then
              write(*,*) NameSub // 'WARNING: ',&
                   'there is no hyperbolic scalar variable in the equation module!'
              if (UseStrict) call stop_mpi('Correct PARAM.in or change equation!')
           end if
        else
           call read_var('UseHyperbolicDivB',UseHyperbolicDivB)
           if(UseHyperbolicDivB) then
              call read_var('SpeedHypDim',SpeedHypDim)
              call read_var('HypDecay'   ,HypDecay)
           end if
        endif

     case("#PROJECTION")                              !^CFG IF PROJECTION BEGIN
        if(.not.UseB)CYCLE READPARAM
        call read_var('TypeProjectIter' ,proj_method)
        call read_var('TypeProjectStop' ,proj_typestop)
        call read_var('RelativeLimit'   ,proj_divbcoeff)
        call read_var('AbsoluteLimit'   ,proj_divbconst)
        call read_var('MaxMatvec'       ,proj_matvecmax)
        ! Make sure that DivbMax is recalculated
        DivbMax = -1.0                                !^CFG END PROJECTION

     case("#CORRECTP")                                !^CFG IF PROJECTION BEGIN
        call read_var('pRatioLow',Pratio_lo)
        call read_var('pRatioHigh',Pratio_hi)
        if(Pratio_lo>=Pratio_hi)&
             call stop_mpi(NameSub//' ERROR: Pratio_lo>=Pratio_hi')
        !                                              ^CFG END PROJECTION

     case("#IOUNITS")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('TypeIoUnit',TypeIoUnit,IsUpperCase=.true.)

     case("#NORMALIZATION")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('TypeNormalization',TypeNormalization,IsUpperCase=.true.)
        select case(TypeNormalization)
        case('NONE')
           No2Si_V = 1.0
        case('READ')
           call read_var('No2SiUnitX',   No2Si_V(UnitX_))
           call read_var('No2SiUnitU',   No2Si_V(UnitU_))
           call read_var('No2SiUnitRho', No2Si_V(UnitRho_))
        case('PLANETARY', 'SOLARWIND')
           ! Depends on other commands, defined in set_physics
        case('USER')
           ! Call user_normalization later in set_units (see set_physics.f90)
           ! to set the normalization units
        case default
           call stop_mpi(NameSub//' ERROR: unknown TypeNormalization=' &
                //TypeNormalization)
        end select

     case("#SHOCKTUBE")
        UseShockTube = .true.
        call split_string(NamePrimitiveVar,nVar,NamePrimitive_V,nVarRead)
        do i=1,nVar
           call read_var(NamePrimitive_V(i)//' left',ShockLeftState_V(i))
        end do
        do i=1,nVar
           call read_var(NamePrimitive_V(i)//' right',ShockRightState_V(i))
        end do

     case("#SHOCKPOSITION")
        call read_var('ShockPosition',ShockPosition)
        call read_var('ShockSlope',ShockSlope)

     case("#SOLARWINDFILE", "#UPSTREAM_INPUT_FILE")
        call read_var('UseSolarWindFile',UseSolarwindFile)
        DoReadSolarwindFile = UseSolarwindFile
        if (UseSolarwindFile) &
             call read_var('NameSolarWindFile', NameSolarWindFile)
        !                                               ^CFG IF RAYTRACE BEGIN

     case("#RAYTRACE")
        call read_var('UseRaytrace',UseRaytrace)
        if(UseRaytrace)then
           call read_var('UseAccurateTrace', UseAccurateTrace)
           call read_var('DtExchangeRay',    DtExchangeRay)
           call read_var('DnRaytrace',       DnRaytrace)
        end if
     case("#RAYTRACELIMIT")
        call read_var('RayLengthMax', RayLengthMax)
     case("#IE")
        call read_var('DoTraceIE', DoTraceIE)
        !                                              ^CFG END RAYTRACE

     case("#IECOUPLING")
        call read_iono_velocity_param

     case("#IMCOUPLING","#IM")                        !^CFG IF RCM BEGIN
        call read_var('TauCoupleIm',TauCoupleIm)
        if(TauCoupleIm < 1.0)then
           TauCoupleIM = 1.0/TauCoupleIM
           if(iProc==0)then
              write(*,'(a)')NameSub//' WARNING: TauCoupleIm should be >= 1'
              if(UseStrict)call stop_mpi('Correct PARAM.in!')
              write(*,*)NameSub//' using the inverse:',TauCoupleIm
           end if
        end if
        call read_var('DoImSatTrace',DoImSatTrace)
        if(NameCommand == "#IMCOUPLING")then
           call read_var('DoCoupleImPressure', DoCoupleImPressure)
           call read_var('DoCoupleImDensity',  DoCoupleImDensity)
           call read_var('DoFixPolarRegion',   DoFixPolarRegion)
           if(DoFixPolarRegion)then
              call read_var('rFixPolarRegion',    rFixPolarRegion)
              do iFluid = IonFirst_, nFluid
                 call read_var('PolarNDim',  PolarNDim_I(iFluid))
                 call read_var('PolarTDim',  PolarTDim_I(iFluid))
              end do
           end if
        end if
     case("#IMCOUPLINGSMOOTH")
        call read_var('dLatSmoothIm', dLatSmoothIm)

     case("#MULTIFLUIDIM")
        if(.not.is_first_session())CYCLE READPARAM
        ! couple GM and IM in multi-fluid (all, Hp, Op) mode                
        call read_var('DoMultiFluidIMCoupling', DoMultiFluidIMCoupling)
        !^CFG END RCM

     case('#ANISOPRESSUREIM')
        if(.not.is_first_session())CYCLE READPARAM
        ! couple GM and IM in anisotropic pressure mode
        call read_Var('DoAnisoPressureIMCoupling', DoAnisoPressureIMCoupling)

     case("#RBSATCOUPLING")
        call read_var('DoRbSatTrace',DoRbSatTrace)
     case("#USERFLAGS", "#USER_FLAGS")
        call read_var('UseUserInnerBcs'         ,UseUserInnerBcs)
        call read_var('UseUserSource'           ,UseUserSource)
        call read_var('UseUserPerturbation'     ,UseUserPerturbation)
        call read_var('UseUserOuterBcs'         ,UseUserOuterBcs)
        call read_var('UseUserICs'              ,UseUserICs)
        call read_var('UseUserSpecifyRefinement',UseUserSpecifyRefinement)
        call read_var('UseUserLogFiles'         ,UseUserLogFiles)
        call read_var('UseUserWritePlot'        ,UseUserWritePlot)
        call read_var('UseUserAMR'              ,UseUserAMR)
        call read_var('UseUserEchoInput'        ,UseUserEchoInput)
        call read_var('UseUserB0'               ,UseUserB0) 
        call read_var('UseUserInitSession'      ,UseUserInitSession)
        call read_var('UseUserUpdateStates'     ,UseUserUpdateStates)

     case("#USERINPUTBEGIN")
        call user_read_inputs
        ! Make sure that MassIon_I is consistent with MassFluid_I
        MassIon_I = MassFluid_I(IonFirst_:IonLast_)

     case("#CODEVERSION")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('CodeVersion',CodeVersionRead)
        if(abs(CodeVersionRead-CodeVersion)>0.005.and.iProc==0)&
             write(*,'(a,f6.3,a,f6.3,a)')NameSub//&
             ' WARNING: CodeVersion in file=',CodeVersionRead,&
             ' but '//NameThisComp//' version is ',CodeVersion,' !!!'

     case("#EQUATION")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('NameEquation',NameEquationRead)
        call read_var('nVar',        nVarRead)
        if(NameEquationRead /= NameEquation .and. iProc==0)then
           write(*,'(a)')'BATSRUS was compiled with equation '// &
                NameEquation//' which is different from '// &
                NameEquationRead
           call stop_mpi(NameSub//' ERROR: incompatible equation names')
        end if
        if(nVarRead /= nVar .and. iProc==0)then
           write(*,'(a,i2,a,i2)')&
                'BATSRUS was compiled with nVar=',nVar, &
                ' which is different from nVarRead=',nVarRead
           call stop_mpi(NameSub//' ERROR: Incompatible number of variables')
        end if

     case("#NEWRESTART","#RESTARTINDIR","#RESTARTINFILE",&
          "#PRECISION","#BLOCKLEVELSRELOADED")
        if(.not.is_first_session())CYCLE READPARAM
        call read_restart_parameters(NameCommand)

     case("#SAVERESTART", "#RESTARTOUTDIR","#RESTARTOUTFILE")
        call read_restart_parameters(NameCommand)

     case("#PLOTDIR")
        call read_var("NamePlotDir",NamePlotDir)
        call fix_dir_name(NamePlotDir)
        if (iProc==0) call check_dir(NamePlotDir)

     case("#SATELLITE")
        if(.not.is_first_session())CYCLE READPARAM
        call read_satellite_parameters(NameCommand)
        DoReadSatelliteFiles = nSatellite > 0

     case('#STEADYSTATESATELLITE')
        call read_satellite_parameters(NameCommand)

     case('#GEOMAGINDICES')
        ! See ModGmGeoindices.f90 for more information.
        DoWriteIndices = .true. ! Activiate geoindices output file.
        DoCalcKp = .true.       ! Kp calculated (no others available.)
        call read_var('nKpWindow', nKpMins)
        call read_var('DtOutput' , dtWriteIndices)
        dt_output(indexfile_) = dtWriteIndices
        dn_output(indexfile_) = -1  ! See ModGmGeoindices.f90.
        nFile = max(nFile, indexfile_)

     case("#MAGNETOMETER")
        DoReadMagnetometerFile = .true.
        save_magnetometer_data = .true.
        call read_var('MagInputFile', MagInputFile)

        if (iProc==0) call check_dir(NamePlotDir)

        call read_var('DnOutput', dn_output(magfile_))
        call read_var('DtOutput', dt_output(magfile_)) 
        nFile = max(nFile, magfile_) 

     case("#GRIDGEOMETRY")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('TypeGeometry', TypeGeometry, IsLowerCase=.true.)
        ! need to read in the general grid file      
        if(TypeGeometry == 'spherical_genr') then
           call read_var('NameGridFile',NameGridFile)
           call read_gen_radial_grid(NameGridFile)
        end if

     case("#GRIDSYMMETRY")
        nMirror_D = 1
        call read_var('IsMirrorX', IsMirrorX)
        if(IsMirrorX) nMirror_D(1) = 2
        call read_var('IsMirrorY', IsMirrorY)
        if(IsMirrorY) nMirror_D(2) = 2
        call read_var('IsMirrorZ', IsMirrorZ)
        if(IsMirrorZ) nMirror_D(3) = 2

     case("#LIMITRADIUS")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('RadiusMin', RadiusMin)
        call read_var('RadiusMax', RadiusMax)

     case("#UNIFORMAXIS")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('UseUniformAxis', UseUniformAxis)

     case("#FIXAXIS")
        call read_var('UsePoleDiffusion', UsePoleDiffusion)
        call read_var('DoFixAxis',DoFixAxis)
        call read_var('rFixAxis',rFixAxis)
        call read_var('r2FixAxis',r2FixAxis)

     case("#GRID")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('nRootBlockX',proc_dims(1)) 
        call read_var('nRootBlockY',proc_dims(2))
        call read_var('nRootBlockZ',proc_dims(3))

        call read_var('xMin',x1)
        call read_var('xMax',x2)
        call read_var('yMin',y1)
        call read_var('yMax',y2)
        call read_var('zMin',z1)
        call read_var('zMax',z2)

     case("#USERMODULE")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('NameUserModule',NameUserModuleRead)
        call read_var('VersionUserModule',VersionUserModuleRead)
        if(NameUserModuleRead /= NameUserModule .or. &
             abs(VersionUserModule-VersionUserModuleRead)>0.001) then
           if(iProc==0)write(*,'(4a,f5.2)') NameSub, &
                ' WARNING: code is compiled with user module ',NameUserModule,&
                ' version',VersionUserModule
           if(UseStrict)call stop_mpi('Select the correct user module!')
        end if

     case("#CHECKGRIDSIZE")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('nI',nIJKRead_D(1))
        call read_var('nJ',nIJKRead_D(2))
        call read_var('nK',nIJKRead_D(3))
        if(any(nIJK_D/=nIJKRead_D).and.iProc==0)then
           write(*,*)'Code is compiled with nI,nJ,nK=',nIJK_D
           call stop_mpi('Change nI,nJ,nK in ModSize.f90 and recompile!')
        end if
        call read_var('MinBlockALL',qtotal)
        if(qtotal>nBLK*nProc .and. iProc==0)then
           write(*,*)'nBLK*nProc=',nBLK*nProc
           call stop_mpi('Use more processors'//&
                ' or increase nBLK in ModSize and recompile!')
        end if

     case("#AMRINITPHYSICS")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('nRefineLevelIC',nRefineLevelIC)

     case("#GAMMA")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('Gamma',g)
        !\
        ! Compute gamma related values.
        !/
        gm1     = g - 1.0
        gm2     = g - 2.0
        gp1     = g + 1.0
        inv_g   = 1.0 / g
        inv_gm1 = 1.0 /gm1
        g_half  = 0.5*g

     case("#LOOKUPTABLE")
        call read_lookup_table_param

     case("#PLASMA")
        do iFluid = IonFirst_, nFluid
           call read_var('MassFluid', MassFluid_I(iFluid))
        end do
        MassIon_I = MassFluid_I(IonFirst_:IonLast_)
        do iFluid = 1, nIonFluid
           call read_var('ChargeIon', ChargeIon_I(iFluid))
        end do

        call read_var('ElectronTemperatureRatio', ElectronTemperatureRatio)
        ElectronPressureRatio = ElectronTemperatureRatio

        !averageioncharge is only useful when there is only one ion
        if(nIonFluid==1)then  
           AverageIonCharge = ChargeIon_I(1)
           ElectronPressureRatio = ElectronTemperatureRatio*AverageIonCharge
        end if

        PePerPtotal = ElectronPressureRatio/(1 + ElectronPressureRatio)

     case("#MULTISPECIES")
        call read_var('DoReplaceDensity', DoReplaceDensity)
        call read_var('SpeciesPercentCheck',SpeciesPercentCheck)

     case("#UNIFORMIONVELOCITY")
        call read_var('UseUniformIonVelocity', UseUniformIonVelocity)

     case("#NEUTRALFLUID")
        call read_var('DoConserveNeutrals', DoConserveNeutrals)
        call read_var('TypeFluxNeutral',    TypeFluxNeutral)

     case("#MULTIION", "#MHDIONS", "#COLLISION")
        call multi_ion_set_parameters(NameCommand)

     case('#USERBOUNDARY', '#EXTRABOUNDARY')
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('UseExtraBoundary', UseExtraBoundary)
        if(UseExtraBoundary)then
           call read_var('TypeBc_I(ExtraBc_)', TypeBc_I(ExtraBc_))      
           ! Backward compatibility
           if (TypeBc_I(ExtraBc_) == 'coronatoih') &
                TypeBc_I(ExtraBc_)= 'buffergrid'
        end if

     case("#FACEBOUNDARY")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('MinBoundary',MinBoundary)
        call read_var('MaxBoundary',MaxBoundary)

     case("#SOLARWIND")
        !if(.not.is_first_session())CYCLE READPARAM
        call read_var('SwNDim',  SW_n_dim)
        call read_var('SwTDim'  ,SW_T_dim)
        call read_var('SwUxDim' ,SW_Ux_dim)
        call read_var('SwUyDim' ,SW_Uy_dim)
        call read_var('SwUzDim' ,SW_Uz_dim)
        call read_var('SwBxDim' ,SW_Bx_dim)
        call read_var('SwByDim' ,SW_By_dim)
        call read_var('SwBzDim' ,SW_Bz_dim)

     case("#OUTFLOWPRESSURE")
        call read_var('UseOutflowPressure', UseOutflowPressure)
        if(UseOutflowPressure) call read_var('pOutflowSi', pOutflowSi)

     case("#MAGNETOSPHERE","#BODY")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('UseBody',body1)
        if(body1)then
           call read_var('rBody', rBody)
           if(NameThisComp=='GM')&
                call read_var('rCurrents' ,Rcurrents)
           do iFluid = IonFirst_, nFluid
              call read_var('BodyNDim', BodyNDim_I(iFluid))
              call read_var('BodyTDim', BodyTDim_I(iFluid))
           end do
        end if

     case("#POLARBOUNDARY")
        do iFluid = IonFirst_, nFluid
           call read_var('PolarNDim',  PolarNDim_I(iFluid))
           call read_var('PolarTDim',  PolarTDim_I(iFluid))
           call read_var('PolarUDim',  PolarUDim_I(iFluid))
        end do
        call read_var('PolarLatitude', PolarLatitude)
        PolarTheta = (90-PolarLatitude)*cDegToRad

     case("#CPCPBOUNDARY")
        call read_var('UseCpcpBc', UseCpcpBc)
        if(UseCpcpBc)then
           if(nFluid > 1)call stop_mpi("CPCBOUNDARY works for single fluid only")
           call read_var('Rho0Cpcp',   Rho0Cpcp)
           call read_var('RhoPerCpcp', RhoPerCpcp)
        end if
     case("#GRAVITY")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('UseGravity',UseGravity)
        if(UseGravity)then
           call read_var('iDirGravity',GravityDir)
           if(GravityDir /= 0) call read_var('GravitySi', GravitySi)
        end if

     case("#SECONDBODY")                        !^CFG IF SECONDBODY BEGIN
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('UseBody2',UseBody2)
        if(UseBody2)then
           call read_var('rBody2',rBody2)
           call read_var('xBody2',xBody2)
           call read_var('yBody2',yBody2)
           call read_var('zBody2',zBody2)
           call read_var('rCurrentsBody2',rCurrentsBody2)
           call read_var('RhoDimBody2',RhoDimBody2)
           call read_var('tDimBody2'  ,tDimBody2)
           call read_var('UseOrbit'  ,UseOrbit)
           if(UseOrbit)then
              call read_var('OrbitPeriod [days]', OrbitPeriod)
              ! Convert orbit period from days to seconds
              OrbitPeriod = OrbitPeriod*cSecondPerDay 
           end if
        end if

     case("#DIPOLEBODY2")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('BdpDimBody2x',BdpDimBody2_D(1))
        call read_var('BdpDimBody2y',BdpDimBody2_D(2))
        call read_var('BdpDimBody2z',BdpDimBody2_D(3))
        !                                           ^CFG END SECONDBODY

     case('#PLANET','#MOON','#COMET','#IDEALAXES','#ROTATIONAXIS',&
          '#MAGNETICAXIS','#MAGNETICCENTER','#ROTATION','#DIPOLE', &
          '#NONDIPOLE','#UPDATEB0')

        call check_stand_alone
        if(.not.is_first_session())CYCLE READPARAM

        call read_planet_var(NameCommand)

     case("#COORDSYSTEM","#COORDINATESYSTEM")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('TypeCoordSystem',TypeCoordSystem,IsUpperCase=.true.)
        select case(NameThisComp)
        case('GM')
           if(TypeCoordSystem /= 'GSM' .and. TypeCoordSystem /= 'GSE') &
                call stop_mpi(NameSub// &
                ' ERROR: cannot handle coordinate system '//TypeCoordSystem)
        case('IH','OH')
           select case(TypeCoordSystem)
           case('HGI')
              ! Note: transformation from HGR to HGI does not work properly
              !       unless HGR happens to be the same as HGC 
              !       (ie. aligned with HGI at the initial time)
              DoTransformToHgi = UseRotatingFrame
              UseRotatingFrame = .false.
           case('HGC','HGR')
              UseRotatingFrame = .true.
           case default
              call stop_mpi(NameSub// &
                   ' ERROR: cannot handle coordinate system '&
                   //TypeCoordSystem)
           end select
        case('SC','LC')
           select case(TypeCoordSystem)
           case('HGR','HGC')
              UseRotatingFrame = .true.
              if(iProc==0)then
                 write(*,*)NameSub//' setting .UseRotatingFrame = T'
              end if
           case('HGI')
              if(iProc==0) write(*,*) NameSub,&
                   ' WARNING: inertial SC/LC is less accurate'
              UseRotatingFrame = .false.
           case default
              call stop_mpi(NameSub// &
                   ' ERROR: cannot handle coordinate system '&
                   //TypeCoordSystem)
           end select
        end select

     case("#NSTEP")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('nStep',n_step)

     case("#NPREVIOUS")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('nPrev',n_prev)             !^CFG IF IMPLICIT
        call read_var('DtPrev',dt_prev)           !^CFG IF IMPLICIT

     case("#STARTTIME", "#SETREALTIME")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('iYear'  ,iStartTime_I(1))
        call read_var('iMonth' ,iStartTime_I(2))
        call read_var('iDay'   ,iStartTime_I(3))
        call read_var('iHour'  ,iStartTime_I(4))
        call read_var('iMinute',iStartTime_I(5))
        call read_var('iSecond',iStartTime_I(6))
        iStartTime_I(7) = 0
        if(IsStandAlone)then
           call time_int_to_real(iStartTime_I, StartTime)
        else
           ! Check if things work out or not
           call time_int_to_real(iStartTime_I, StartTimeCheck)
        end if

     case("#TIMESIMULATION")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('tSimulation',time_simulation)

     case("#HELIOUPDATEB0")
        if(.not.UseB0)CYCLE READPARAM
        call read_var('DtUpdateB0', DtUpdateB0)
        DoUpdateB0 = DtUpdateB0 > 0.0

     case("#HELIODIPOLE")
        if(.not.is_first_session())CYCLE READPARAM
        if(.not.UseB0)CYCLE READPARAM
        call read_var('HelioDipoleStrengthSi',DipoleStrengthSi)
        call read_var('HelioDipoleTilt'      ,ThetaTilt)
        ThetaTilt = ThetaTilt * cDegToRad

     case("#HELIOROTATION", "#INERTIAL")
        if(iProc==0)write(*,*) NameSub, ' WARNING: ',&
             ' #HELIOROTATION / #INERTIAL command is obsolete and ignored'

     case("#HELIOBUFFERGRID")
        if(.not.is_first_session())CYCLE READPARAM
        if(NameThisComp == 'SC'.or. NameThisComp == 'LC') &
             call stop_mpi(NameSub//' ERROR:'// &
             ' #HELIOBUFFERGRID command can be used in IH,OH components only')
        call read_var('rBuffMin',  rBuffMin)
        call read_var('rBuffMax',  rBuffMax)
        call read_var('nThetaBuff',nThetaBuff)
        call read_var('nPhiBuff',  nPhiBuff)

     case("#BUFFERGRID")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('nRBuff',      nRBuff)
        call read_var('nPhiBuff',    nPhiBuff)
        call read_var('nThetaBuff',  nThetaBuff)
        call read_var('rBuffMin',    BufferMin_D(BuffR_))
        call read_var('rBuffMax',    BufferMax_D(BuffR_))
        call read_var('PhiBuffMin',  BufferMin_D(BuffPhi_))
        call read_var('PhiBuffMax',  BufferMax_D(BuffPhi_))
        call read_var('LatBuffMin',  BufferMin_D(BuffTheta_))
        call read_var('LatBuffMax',  BufferMax_D(BuffTheta_))
        ! Update rBuffMin, rBuffMax 
        ! in case #HELIOBUFFERGRID was not used
        rBuffMin = BufferMin_D(BuffR_)
        rBuffMax = BufferMax_D(BuffR_)

        ! Convert degrees to radians, latitude to co-latitude
        BufferMin_D(BuffPhi_)   = BufferMin_D(BuffPhi_) * cDegToRad
        BufferMax_D(BuffPhi_)   = BufferMax_D(BuffPhi_) * cDegToRad
        BufferMin_D(BuffTheta_) = cHalfPi - BufferMin_D(BuffTheta_) * cDegToRad
        BufferMax_D(BuffTheta_) = cHalfPi - BufferMax_D(BuffTheta_) * cDegToRad

     case("#THINCURRENTSHEET")
        call read_var('DoThinCurrentSheet', DoThinCurrentSheet)

        ! OUTERHELIOSPHERE SPECIFIC COMMANDS

     case("#OHNEUTRALS")
        call read_var('RhoNeutralsISW_dim' ,RhoNeutralsISW_dim)
        call read_var('TNeutralsISW_dim' ,TNeutralsISW_dim)
        call read_var('UxNeutralsISW_dim' ,UxNeutralsISW_dim)
        call read_var('UyNeutralsISW_dim' ,UyNeutralsISW_dim)
        call read_var('UzNeutralsISW_dim' ,UzNeutralsISW_dim)
        call read_var('mProtonMass',mProtonMass)

     case("#OHBOUNDARY")
        call read_var('DoOhNeutralBc',DoOhNeutralBc)
        do iFluid = 2,nFluid
           call read_var('RhoBcFactor', RhoBcFactor_I(iFluid))
           call read_var('uBcFactor'  , uBcFactor_I(iFluid))
        end do

        !CORONA SPECIFIC COMMANDS

     case("#MAGNETOGRAM", "#READPOTENTIALFIELD")
        call read_var('UseMagnetogram', UseMagnetogram)
        if(UseMagnetogram)&
             call set_parameters_magnetogram(NameCommand)

     case("#SAVEPOTENTIALFIELD")
        call set_parameters_magnetogram(NameCommand)

     case('#LDEM')
        call read_var('UseLdem', UseLdem)
        if(UseLdem) then
           call read_var('NameLdemFile', NameLdemFile)
           call read_var('iRadiusLdem', iRadiusLdem)
           call read_ldem(NamePlotDir)
        end if

     case("#EMPIRICALSW")
        call read_var('NameModelSW', NameModelSW)
        UseEmpiricalSW = NameModelSW /= 'none'

     case("#WSACOEFF")
        call read_wsa_coeff 

     case("#CORONALHEATING", "#LONGSCALEHEATING", "#ACTIVEREGIONHEATING", &
          "#ELECTRONHEATING")
        call read_corona_heating(NameCommand)

     case("#OPENCLOSEDHEAT")
        call read_var('DoOpenClosedHeat', DoOpenClosedHeat)

     case("#RADIATIVECOOLING")
        call read_var('UseRadCooling', UseRadCooling)

     case("#CHROMOSPHERE")
        call read_chromosphere

     case("#TRANSITIONREGION")
        call read_modified_cooling

     case("#ADVECTWAVES", "#ALFVENWAVES", "#WAVEPRESSURE", &
          "#FREQUENCY", "#SPECTRUM", "#WAVEREFLECTION")
        call read_waves_param(NameCommand)

     case("#LASERPULSE", "#LASERBEAM", "#LASERBEAMS", "#LASERBEAMPROFILE", &
          "#LASERRAYTEST")
        call read_laser_heating_param(NameCommand)

     case("#CME", "#ARCH", "#TD99FLUXROPE", "#GL98FLUXROPE", "#SHEARFLOW")
        call EEE_set_parameters(NameCommand)

     case("#STAR")
        UseStar=.true.
        call read_var('RadiusStar',         RadiusStar)
        call read_var('MassStar',           MassStar)
        call read_var('RotationPeriodStar', RotationPeriodStar)        

     case default
        if(iProc==0) then
           write(*,*) NameSub // ' WARNING: unknown #COMMAND ' // &
                trim(NameCommand),' !!!'
           if(UseStrict)call stop_mpi('Correct PARAM.in!')
        end if
     end select
  end do READPARAM

  ! end reading parameters

contains

  !===========================================================================
  subroutine check_stand_alone

    if(IsStandAlone) RETURN
    if(iProc==0) write(*,*) NameSub,' WARNING: command '//trim(NameCommand)//&
         ' is allowed in stand alone mode only !!!'
    if(UseStrict)call stop_mpi(NameSub//' Correct PARAM.in')

  end subroutine check_stand_alone

  !===========================================================================

  logical function is_first_session()

    is_first_session = iSession == 1

    if(iSession/=1 .and. iProc==0)then
       write(*,*)NameSub//' WARNING: command '//trim(NameCommand)// &
            ' can be used in the first session only !!!'
       if(UseStrict)call stop_mpi('Correct PARAM.in')
    end if

  end function is_first_session

  !===========================================================================

  subroutine set_namevar

    use ModUtilities, ONLY: join_string

    integer :: iWave
    character(len=3):: NameWave
    integer :: iMaterial
    character(len=2):: NameMaterial
    !-------------------------------------------------------------------------

    ! Fix the NameVar_V string for waves
    if(WaveLast_ > 1)then
       do iWave = 1, nWave
          write(NameWave,'(a,i2.2)') 'I',iWave
          NameVar_V(WaveFirst_+iWave-1) = NameWave
       end do
    end if

    ! Fix the NameVar_V string for material levels
    if(MaterialLast_ > 1)then
       do iMaterial = 1, nMaterial
          write(NameMaterial,'(a,i1.1)') 'M',iMaterial
          NameVar_V(MaterialFirst_+iMaterial-1) = NameMaterial
       end do
    end if

    ! space separated NameVar string containing all variable names
    call join_string(nVar, NameVar_V, NameVarCouple)

  end subroutine set_namevar

  !===========================================================================

  subroutine set_defaults

    !\
    ! Default plot and restart directories depend on NameThisComp
    !/
    NamePlotDir(1:2) = NameThisComp

    ! Set defaults for restart files
    call init_mod_restart_file

    CodeVersionRead = -1.

    ! Default coordinate systems
    select case(NameThisComp)
    case('IH','OH')
       ! Do not start line with Type... to avoid an Emacs indentation bug
       UseRotatingFrame  = .false.
       UseRotatingBc     = .false.; TypeCoordSystem   = 'HGI'
    case('SC','LC')
       UseRotatingFrame  = .true.
       UseRotatingBc     = .false.; TypeCoordSystem   = 'HGR'
    case('GM')
       UseRotatingFrame  = .false.
       UseRotatingBc     = .true.;  TypeCoordSystem   = 'GSM'
    case('EE')
       UseRotatingFrame  = .false.
       UseRotatingBc     = .false.; TypeCoordSystem   = 'GSM'
    end select

    !Do not set B0 field in IH and OH
    if(NameThisComp/='IH'.and.NameThisComp/='OH')then
       UseB0=UseB
    else
       UseB0=.false.
    end if

    ! Do not update B0 for LC, SC or IH by default
    if(NameThisComp /= 'GM')then
       DoUpdateB0 = .false.
       DtUpdateB0 = -1.0
    end if

    ! Initialize StartTime to the default values
    ! For SWMF it is set during 'CHECK'
    if(IsStandAlone)call time_int_to_real(iStartTime_I, StartTime)

    iterTEST        =-1
    tTEST           =1.0E30
    Itest           =1
    Jtest           =1
    Ktest           =1
    BLKtest         =1
    PROCtest        =0
    VARtest         =1
    DIMtest         =1

    nStage        = 2
    Cfl           = 0.8
    CflOrig       = 0.8
    UseDtFixed    = .false.
    dt            = 0.0
    dt_BLK        = 0.0

    initial_refine_levels = 0
    nRefineLevelIC        = 0

    min_block_level =  0
    max_block_level = 99
    min_cell_dx =     0.
    max_cell_dx = 99999.

    maxTotalBlocks = nBLK*nProc

    DnAmr=-1
    DoAmr=.false.
    automatic_refinement = .false.

    nOrder = 2
    FluxType = 'RUSANOV'               !^CFG IF RUSANOVFLUX
    !FluxType = 'SOKOLOV'              !^CFG UNCOMMENT IF NOT RUSANOVFLUX

    ! Default implicit parameters      !^CFG IF IMPLICIT BEGIN
    UseImplicit      = .false.
    ImplCritType     = 'dt'

    if(nByteReal>7)then
       JacobianEps   = 1.E-12
    else
       JacobianEps   = 1.E-6
    end if                            !^CFG END IMPLICIT

    UseDivbSource   =  UseB .and. nDim > 1
    UseDivbDiffusion= .false.         !^CFG IF DIVBDIFFUSE
    UseProjection   = .false.         !^CFG IF PROJECTION
    UseConstrainB   = .false.         !^CFG IF CONSTRAINB

    UseB0Source     = UseB0
    UseHyperbolicDivB = Hyp_ > 1

    UseUpdateCheck  = .true.
    ! The use of (/../) is correct F90, but it is replaced
    ! with setting the elements to avoid a compiler bug in
    ! the Portland Group F90 compiler PGF90/any Linux/x86 3.3-2
    !    percent_max_rho = (/40., 400./)
    !    percent_max_p   = (/40., 400./)
    percent_max_rho(1) = 40.
    percent_max_rho(2) = 400.
    percent_max_p(1)   = 40.
    percent_max_p(2)   = 400.

    optimize_message_pass = 'opt'

    plot_dimensional      = .true.

    restart           = .false.

    !\
    ! Give some "reasonable" default values
    !/

    DipoleStrengthSi = 0.0

    UseBody2 = .false.                                !^CFG IF SECONDBODY BEGIN
    RBody2 =-1.0
    xBody2 = 0.0
    yBody2 = 0.0			
    zBody2 = 0.0			
    BdpBody2_D  = 0.0
    rCurrentsBody2 = 0.0
    RhoDimBody2 = 1.0    ! n/cc
    TDimBody2   = 10000.0! K                          !^CFG END SECONDBODY

    MassIon_I = MassFluid_I(IonFirst_:IonLast_) ! Ion masses

    ! Default limits for  buffergrid boundary conditions
    if (any(TypeBc_I(body1_:ExtraBc_)=='buffergrid')) then
       BufferMin_D(BuffR_)     = 19.
       BufferMin_D(BuffPhi_)   = 0.
       BufferMin_D(BuffTheta_) = 0.
       BufferMax_D(BuffR_)     = 21.
       BufferMax_D(BuffPhi_)   = 360.
       BufferMax_D(BuffTheta_) = 180.
    end if

    nAmrCriteria = 0
    call init_mod_amr(nAmrCriteria)

    !\
    ! Set component dependent defaults
    !/

    GravityDir=0
    if(allocated(TypeConservCrit_I)) deallocate(TypeConservCrit_I)

    select case(NameThisComp)
    case('LC','SC','IH','OH')
       ! Body parameters
       UseGravity=.true.
       body1      =.true.
       Rbody      = 1.00
       Rcurrents  =-1.00

       ! Non Conservative Parameters
       UseNonConservative   = .false.
       nConservCrit         = 0
       rConserv             = -1.

       ! Boundary Conditions
       TypeBc_I(1:6)  = 'float'
       TypeBc_I(body1_)      = 'unknown'
       BodyTDim_I            = 2.85E06    ! K
       BodyNDim_I(IonFirst_) = 1.50E8     ! /cc  protons
       BodyNDim_I(IonFirst_+1:nFluid) = BodyNDim_I(IonFirst_)*cTiny

       ! Normalization and I/O units
       TypeNormalization     = "SOLARWIND"
       TypeIoUnit            = "HELIOSPHERIC"

    case('GM')
       ! Body Parameters
       UseGravity=.false.
       body1      =.true.
       Rbody      = 3.00
       Rcurrents  = 4.00

       ! Non Conservative Parameters
       UseNonConservative   = .true.
       nConservCrit         = 1
       allocate( TypeConservCrit_I(nConservCrit) )
       TypeConservCrit_I(1) = 'r'
       rConserv             = 2*rBody

       ! Boundary Conditions and Normalization
       TypeBc_I(1)        ='outflow'
       TypeBc_I(2)        ='inflow'
       TypeBc_I(3:6)  ='fixed'
       TypeBc_I(body1_)='ionosphere'
       BodyTDim_I    = 25000.0          ! K
       BodyNDim_I    = 5.0              ! /cc

       ! Normalization and I/O units
       TypeNormalization = "PLANETARY"
       TypeIoUnit        = "PLANETARY"

    case('EE')
       ! Body Parameters
       UseGravity = .true.
       GravityDir = 3
       GravitySi  = -273.06667
       body1      = .false.
       Rbody      = 1.00
       Rcurrents  =-1.00

       ! Non Conservative Parameters
       UseNonConservative = .false.
       nConservCrit       = 0
       rConserv           = -1.

       ! Boundary Conditions and Normalization
       TypeBc_I(1:4) = 'periodic'
       TypeBc_I(5)         = 'fixvalue'
       TypeBc_I(6)         = 'no_inflow'
       TypeBc_I(body1_)       = 'unknown'
       BodyTDim_I    = 25000.0          ! K
       BodyNDim_I    = 5.0              ! /cc

       ! Normalization and I/O units
       TypeNormalization = "READ"
       No2Si_V(UnitX_)   = 1.5e5
       No2Si_V(UnitU_)   = 6.4e3
       No2Si_V(UnitRho_) = 2.7e-4
       TypeIoUnit        = "HELIOSPHERIC"

    end select

  end subroutine set_defaults

  !=========================================================================
  subroutine correct_parameters

    use ModMultiFluid, ONLY: UseMultiIon
    use ModWaves, ONLY: UseAlfvenWaves, UseWavePressure
    use ModImplHypre, ONLY: IsHypreAvailable

    ! option and module parameters
    character (len=40) :: Name
    real               :: Version
    logical            :: IsOn

    real    :: BetaProlongOrig = 0.0
    logical :: IsFirstCheck = .true.
    !---------------------------------------------------------------------

    !\
    ! Check for some combinations of things that cannot be accepted as input
    !/
    if (iProc==0) write (*,*) ' '

    if(IsFirstCheck)then
       call correct_grid_geometry
       if(UseConstrainB .or. UseFaceIntegral4) then
          ! Extend face index range in the orthogonal direction
          ! The CT scheme needs 1 extra layers
          ! the 4th order face integrals need 1 and 2 ghost layers
          iMinFace = 0; iMaxFace = nI+1
          jMinFace = 1 - min(1,nJ-1); jMaxFace = nJ + min(1,nJ-1)
          kMinFace = 1 - min(1,nK-1); kMaxFace = nK + min(1,nK-1)
          if(UseFaceIntegral4)then
             iMinFace2 = -1; iMaxFace2 = nI+2
             jMinFace2 = 1 - 2*min(1,nJ-1); jMaxFace2 = nJ + 2*min(1,nJ-1)
             kMinFace2 = 1 - 2*min(1,nK-1); kMaxFace2 = nK + 2*min(1,nK-1)
          else
             iMinFace2 = iMinFace; iMaxFace2 = iMaxFace
             jMinFace2 = jMinFace; jMaxFace2 = jMaxFace
             kMinFace2 = kMinFace; kMaxFace2 = kMaxFace
          end if
       end if
    end if

    ! This depends on the grid geometry set above
    call correct_plot_range

    ! Fix resolutions (this depends on domain size set above)
    if(InitialResolution > 0.0) initial_refine_levels = nint( &
         alog(((XyzMax_D(x_)-XyzMin_D(x_)) / (proc_dims(x_) * nI))  &
         / InitialResolution) / alog(2.0) )

    if(UseTiming)then
       call timing_version(IsOn,Name,Version)
       if(.not.IsOn)then
          if(iProc==0)then
             write(*,'(a)')NameSub//' WARNING: TIMING module is OFF !!!'
             if(UseStrict)call stop_mpi( &
                  'Correct PARAM.in or switch TIMING ON')
             write(*,*)'setting UseTiming=.false.'
          end if
          UseTiming=.false.
       end if
    end if

    ! Check flux type selection
    select case(FluxType)
    case('SIMPLE','Simple')
       FluxType='Simple'
    case('ROE','Roe')                                !^CFG IF ROEFLUX BEGIN
       FluxType='Roe'
       if(UseAlfvenWaves .or. UseWavePressure)then
          if(iProc==0) write(*,'(a)')NameSub // &
               'Wave transport and wave pressure do not work with ' // &
               trim(FluxType)
          if(UseStrict)&                             !^CFG IF AWFLUX
               call stop_mpi('Correct PARAM.in')    
          FluxType='Sokolov'                         !^CFG IF AWFLUX
       end if
    case('ROEOLD','RoeOld')             
       FluxType='RoeOld'                             !^CFG END ROEFLUX
    case('RUSANOV','TVDLF','Rusanov')                !^CFG IF RUSANOVFLUX
       FluxType='Rusanov'                            !^CFG IF RUSANOVFLUX
    case('LINDE','HLLEL','Linde')                    !^CFG IF LINDEFLUX
       FluxType='Linde'                              !^CFG IF LINDEFLUX
    case('SOKOLOV','AW','Sokolov')                   !^CFG IF AWFLUX
       FluxType='Sokolov'                            !^CFG IF AWFLUX
    case('HLLD')                                     !^CFG IF HLLDFLUX
    case('GODUNOV','Godunov')
       FluxType='Godunov'
    case default
       if(iProc==0)then
          write(*,'(a)')NameSub // &
               'WARNING: unknown value for FluxType=' // trim(FluxType)
          if(UseStrict) &                            !^CFG IF RUSANOVFLUX
               call stop_mpi('Correct PARAM.in!')
          write(*,*)'setting FluxType=Rusanov'       !^CFG IF RUSANOVFLUX
       end if
       FluxType='Rusanov'                            !^CFG IF RUSANOVFLUX
    end select

    ! Set flux type for neutral fluids
    select case(TypeFluxNeutral)
    case('default')
       select case(FluxType)
       case('Rusanov', 'Linde')
          TypeFluxNeutral = FluxType
       case default
          TypeFluxNeutral = 'Linde'
       end select
    case('RUSANOV','TVDLF','Rusanov')
       TypeFluxNeutral = 'Rusanov'
    case('LINDE','HLLEL','Linde')
       TypeFluxNeutral = 'Linde'
    case default
       if(iProc==0)then
          write(*,'(a)')NameSub// &
               ' WARNING: Unknown value for TypeFluxNeutral='// &
               trim(TypeFluxNeutral)//' !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' setting TypeFluxNeutral=Linde'
       end if
       TypeFluxNeutral = 'Linde'
    end select

    ! Check flux type selection for implicit   !^CFG IF IMPLICIT BEGIN
    select case(FluxTypeImpl)
    case('default')
       FluxTypeImpl = FluxType
    case('ROE','Roe')                                !^CFG IF ROEFLUX BEGIN
       FluxTypeImpl='Roe'
    case('ROEOLD','RoeOld')
       FluxTypeImpl='RoeOld'                         !^CFG END ROEFLUX
    case('RUSANOV','TVDLF','Rusanov')                !^CFG IF RUSANOVFLUX
       FluxTypeImpl='Rusanov'                        !^CFG IF RUSANOVFLUX
    case('LINDE','HLLEL','Linde')                    !^CFG IF LINDEFLUX
       FluxTypeImpl='Linde'                          !^CFG IF LINDEFLUX
    case('SOKOLOV','AW','Sokolov')                   !^CFG IF AWFLUX
       FluxTypeImpl='Sokolov'                        !^CFG IF AWFLUX
    case('HLLD')                                     !^CFG IF HLLDFLUX
    case('GODUNOV','Godunov')
       FluxTypeImpl='Godunov'
    case default
       if(iProc==0)then
          write(*,'(a)')NameSub// &
               ' WARNING: Unknown value for FluxTypeImpl='// &
               trim(FluxTypeImpl)//' !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' setting FluxTypeImpl=',trim(FluxType)
       end if
       FluxTypeImpl=FluxType
    end select                               !^CFG END IMPLICIT

    if(UseSemiImplicit .or. nStage == 1)then !^CFG IF IMPLICIT BEGIN
       UseBDF2 = .false.
    elseif (time_accurate .and. nStage == 2)then
       UseBDF2 = .true.
    end if                                   !^CFG END IMPLICIT

    ! Make sure periodic boundary conditions are symmetric
    if(any(TypeBc_I(1:2)=='periodic')) TypeBc_I(1:2)='periodic'
    if(any(TypeBc_I(3:4)=='periodic')) TypeBc_I(3:4)='periodic'
    if(any(TypeBc_I(5:6)=='periodic')) TypeBc_I(5:6)='periodic'

    if(UseConstrainB .and. .not.time_accurate)then  !^CFG IF CONSTRAINB BEGIN
       if(iProc==0)then
          write(*,'(a)')NameSub//&
               ' WARNING: constrain_B works for time accurate run only !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' setting UseConstrainB=F UseDivbSource=T'
       end if
       UseConstrainB=.false.
       UseDivbSource=.true.
    end if
    if(UseConstrainB .and. (UseTvdReschange .or. UseAccurateResChange))then
       if(iProc==0)write(*,'(a)')NameSub//&
            ' WARNING: cannot use TVD or accurate schemes at res. change' // &
            ' with ConstrainB'
       !if(UseStrict)call stop_mpi('Correct PARAM.in!')
       if(iProc==0)write(*,*)NameSub// &
            ' setting UseTvdReschange=F UseAccurateResChange=F'
       UseTvdReschange      = .false.
       UseAccurateResChange = .false.
    end if
    if (UseConstrainB) optimize_message_pass = 'all'
    if (UseConstrainB .and. DoAmr)then
       if(iProc==0)write(*,'(a)')NameSub//&
            ' WARNING: cannot use AMR with constrained transport'
       if(UseStrict)call stop_mpi('Correct PARAM.in!')
       if(iProc==0)write(*,*)NameSub//' setting DoAmr=F'
       DoAmr = .false.
    end if                                            !^CFG END CONSTRAINB

    if(UseViscosity) optimize_message_pass = 'all'

    if ( UseHallResist &
         .or. UseResistivity &                       !^CFG IF DISSFLUX
         ) optimize_message_pass = 'all'

    !^CFG IF IMPLICIT BEGIN
    if ( UseRadDiffusion .and. (UseFullImplicit .or. UseSemiImplicit)) &
         optimize_message_pass = 'all'
    !^CFG END IMPLICIT

    !Check for magnetogram

    if(UseEmpiricalSW.and..not.UseMagnetogram)&
         call stop_mpi(&
         'Empirical Solar Wind model requires magnetogram')

    if(DoOpenClosedHeat.and.(.not.UseMagnetogram.and.&
         (iSession==1.and.i_line_command('#PFSSM')<0)))&
         call stop_mpi(&
         'The heating in the closed field region requires magnetogram')

    ! Check algorithm at resolution changes
    if(prolong_order/=1 .and. optimize_message_pass(1:3)=='all')&
         call stop_mpi(NameSub// &
         'The prolongation order=2 requires message_pass_dir')

    if(nOrder == 1)then
       BetaProlongOrig = BetaProlong    
       BetaProlong = 0.0
    else
       BetaProlong = max(BetaProlong, BetaProlongOrig)
    end if

    if(nK == 1 .and. UseTvdResChange) then
       ! in 1D and 2D only accurate reschange is implemented
       UseTvdResChange      = .false.
       UseAccurateResChange = .true.
    end if

    if (UseAccurateResChange) optimize_message_pass = 'all'

    ! Check test processor
    if(ProcTest > nProc)then
       if(iProc==0) write(*,'(a)')NameSub//&
            ' WARNING: procTEST > nProc, setting procTEST=0 !!!'
       procTEST=0
    end if

    if(IsRzGeometry .and. UseB)then
       if(UseMultiIon) &
            call stop_mpi('RZ geometry is not implemented for multi-ion')
       if(UseHallResist) &
            call stop_mpi('RZ geometry is not implemented for Hall-MHD')
    end if

    if(.not.IsCartesian)then
       if(UseProjection)call stop_mpi(&                   !^CFG IF PROJECTION
            'Only Cartesian works with projection')       !^CFG IF PROJECTION
       if(UseConstrainB)call stop_mpi(&                   !^CFG IF CONSTRAINB
            'Only Cartesian works with constrain B')      !^CFG IF CONSTRAINB
       if(UseDivBDiffusion)call stop_mpi(&                !^CFG IF DIVBDIFFUSE
            'Only Cartesian works with divB diffusion')   !^CFG IF DIVBDIFFUSE
    end if

    if(UseHyperbolicDivb)then
       if(.false.&
            .or.UseDivbDiffusion & !^CFG IF DIVBDIFFUSE
            .or. UseConstrainB &   !^CFG IF CONSTRAINB
            .or. UseProjection &   !^CFG IF PROJECTION
            )then
          if(iProc==0)then
             write(*,'(a)')NameSub// &
                  ' WARNING: Hyperbolic cleaning can only be combined', &
                  ' with 8-wave scheme, not any other div B scheme'
             if (UseStrict) call stop_mpi('Correct PARAM.in!')
             write(*,*)NameSub//' Setting UseHyperbolicDivb = .false.'
          end if
          UseHyperbolicDivb = .false.
       end if
    end if

    ! Check CFL number
    if(.not.time_accurate .and. iProc==0)then
       if(boris_correction)then                     !^CFG IF BORISCORR BEGIN
          if(cfl>0.65) then
             write(*,'(a)')NameSub// &
                  ' WARNING: CFL number above 0.65 may be unstable'// &
                  ' for local timestepping with Boris correction !!!'
             if (UseStrict) call stop_mpi('Correct PARAM.in!')
          end if
       else                                         !^CFG END BORISCORR
          if(cfl>0.85)then
             write(*,'(a)')NameSub// &
                  ' WARNING: CFL number above 0.85 may be unstable'// &
                  ' for local timestepping !!!'
             if (UseStrict) call stop_mpi('Correct PARAM.in!')
          end if
       end if                                       !^CFG IF BORISCORR
    end if

    !Boris correction checks                        !^CFG IF BORISCORR BEGIN
    if((FluxType=='Roe' &                              !^CFG IF ROEFLUX BEGIN
         .or. FluxTypeImpl=='Roe' &                       !^CFG IF IMPLICIT
         ) .and. boris_correction)then
       if (iProc == 0) then
          write(*,'(a)')NameSub//&
               ' WARNING: Boris correction not available for Roe flux !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' Setting boris_correction = .false.'
       end if
       boris_correction = .false.
       boris_cLIGHT_factor = 1.00
    end if                                             !^CFG END ROEFLUX
    ! Finish checks for boris                       !^CFG END BORISCORR

    ! Check parameters for implicit                 !^CFG IF IMPLICIT BEGIN

    if(UseClimit .and. .not. UseImplicit)then
       if(iProc==0)then
          write(*,'(a)')'UseClimit=T requires (part) implicit scheme'
          call stop_mpi('Correct PARAM.in')
       end if
    end if

    if(UsePartImplicit .and. ImplCritType=='dt' .and.&
         (.not.time_accurate .or. .not.UseDtFixed))then
       if(iProc==0)then
          write(*,'(a)')'Part implicit scheme with ImplCritType=dt'
          write(*,'(a)')'requires time accurate run with fixed time step'
          call stop_mpi('Correct PARAM.in')
       end if
    end if

    if(PrecondType == "HYPRE")then
       if(.not.IsHypreAvailable)call  stop_mpi(NameSub// &
            ' empty HYPRE module! Use Config.pl -hypre')

       if(.not.UseSemiImplicit)call stop_mpi(NameSub// &
            ' HYPRE preconditioner only works with semi-implicit scheme')

       if(KrylovType == "CG")call stop_mpi(NameSub// &
            ' HYPRE preconditioner does not work with CG. Use GMRES in #KRYLOV command')

    endif

    if(nKrylovVector > KrylovMatvecMax)then
       if(iProc==0)then
          write(*,'(a)')'nKrylovVector>KrylovMatvecMax is useless!'
          write(*,*)'reducing nKrylovVector to KrylovMatvecMax'
       endif
       nKrylovVector = KrylovMatvecMax
    endif

    if(PrecondSide /= 'left' .and. PrecondSide/='symmetric' .and. &
         PrecondSide/='right')then
       if(iProc==0)then
          write(*,'(3a)')NameSub//' WARNING: PrecondSide=',PrecondSide,&
               ' is not one of left, symmetric, or right !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' setting PrecondSide=symmetric'
       end if
       PrecondSide='symmetric'
    end if

    select case(KrylovInitType)
    case('explicit','scaled','nul')
    case default
       if(iProc==0)then
          write(*,'(2a)')NameSub//' WARNING: KrylovInitType=',&
               KrylovInitType, &
               ' is not one of explicit, scaled, or nul !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' setting KrylovInitType=nul'
       end if
       KrylovInitType='nul'
    end select

    if(KrylovInitType/='nul'.and.PrecondSide=='right')then
       if(iProc==0)then
          write(*,'(2a)')NameSub//&
               ' WARNING: PrecondSide=right only works with',&
               ' KrylovInitType=nul !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' setting KrylovInitType=nul'
       end if
       KrylovInitType='nul'
    endif

    if(.not.time_accurate.and.UseBDF2)then
       if(iProc==0)then
          write(*,'(a)') NameSub//&
               ' WARNING: BDF2 is only available for time accurate run !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' setting UseBDF2=.false.'
       end if
       UseBDF2=.false.
    end if

    if(UseImplicit .and. UsePartSteady)then
       if(iProc==0)then
          write(*,'(a)') NameSub//&
               ' WARNING: Implicit and part steady schemes do not'//&
               ' work together !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*) NameSub//' setting UsePartSteady=F'
       endif
       UsePartSteady = .false.
    end if

    if(UsePointImplicit .and. UseFullImplicit)then
       if(iProc==0)then
          write(*,'(a)') NameSub//&
               ' WARNING: Full and point implicit schemes do not'//&
               ' work together !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*) NameSub//' setting UsePointImplicit=F'
       endif
       UsePointImplicit = .false.
    end if

    !Finish checks for implicit                     !^CFG END IMPLICIT

    !Set min_block_level and max_block_level for #AMRRESOLUTION in session 1
    if(i_line_command("#AMRRESOLUTION", iSessionIn = 1) > 0)then
       local_root_dx = (XyzMax_D(x_)-XyzMin_D(x_))/real(proc_dims(1)*nI)
       if    (max_cell_dx < -1.E-6) then
          min_block_level = -1
       elseif(max_cell_dx <  1.E-6) then
          min_block_level = 99
       else
          do j=1,99
             min_block_level = j-1
             if ( local_root_dx/(2**j) < max_cell_dx) EXIT
          end do
       end if
       if    (min_cell_dx < -1.E-6) then
          max_block_level = -1
       elseif(min_cell_dx <  1.E-6) then
          max_block_level = 99
       else
          do j=1,99
             max_block_level = j-1
             if ( local_root_dx/(2**j) < min_cell_dx) EXIT
          end do
       end if
    end if

    if(TypeGeometry == 'cartesian')then                               
       if(UsePoleDiffusion .or. DoFixAxis)then
          UsePoleDiffusion = .false.
          DoFixAxis = .false.
          if(iProc == 0) write(*,*) NameSub// &
               ' setting UsePoleDiffusion and DoFixAxis to FALSE'
       end if
    end if

    if(DoFixAxis .and. .not. UseUniformAxis) &
         call stop_mpi("DoFixAxis=T only works with UseUniformAxis=T!")

    IsFirstCheck = .false.

  end subroutine correct_parameters

  !===========================================================================
  subroutine correct_grid_geometry

    use ModGeometry, ONLY: LogRGen_I
    use BATL_lib, ONLY: init_batl, CoordMin_D, CoordMax_D
    use ModUser,  ONLY : user_specify_refinement
    character(len=20):: TypeGeometryBatl
    !-----------------------------------------------------------------------

    if(i_line_command("#GRID", iSessionIn = 1) < 0) &
         call stop_mpi(NameSub // &
         ' #GRID command must be specified in the first session!')

    if(product(proc_dims) > nBLK*nProc .and. iProc==0)then
       write(*,*)'Root blocks will not fit on all processors, check nBLK'
       write(*,*)'product(proc_dims)   =',product(proc_dims)
       write(*,*)'nBLK,nProc,nBLK*nProc=',nBLK,nProc,nBLK*nProc
       call stop_mpi('product(proc_dims) > nBLK*nProc!')
    end if

    ! Set XyzMin_D, XyzMax_D based on 
    ! #GRID, #GRIDGEOMETRY, and #LIMITRADIUS
    select case(TypeGeometry)
    case('cartesian', 'round')
       XyzMin_D = (/x1, y1, z1/)
       XyzMax_D = (/x2, y2, z2/)
    case('rz')
       z1 = -0.5
       z2 = +0.5
       XyzMin_D = (/x1, y1, z1/)
       XyzMax_D = (/x2, y2, z2/)
    case('spherical', 'spherical_lnr', 'spherical_genr')
       !             R,   Phi, Latitude
       XyzMin_D = (/ 0.0, 0.0, -cHalfPi/)
       XyzMax_D = (/ &
            sqrt(max(x1**2,x2**2) + max(y1**2,y2**2) + max(z1**2,z2**2)), &
            cTwoPi, cHalfPi /)
    case('cylindrical', 'cylindrical_lnr', 'cylindrical_genr')
       !            R,   Phi, Z
       XyzMin_D = (/0.0, 0.0, z1/) 
       XyzMax_D = (/sqrt(max(x1**2,x2**2) + max(y1**2,y2**2)), cTwoPi, 0.0/)
    end select

    !if(index(TypeGeometry,'_lnr')>0) XyzMax_D(R_)=log(XyzMax_D(R_))
    !if(index(TypeGeometry,'_genr')>0) call radius_to_gen(XyzMax_D(R_))

    if(i_line_command("#LIMITRADIUS", iSessionIn = 1) > 0) then 
       XyzMin_D(1) = RadiusMin
       XyzMax_D(1) = RadiusMax
    else
       if(Body1 .and. rBody > 0.0)then
          ! Set inner boundary to match rBody for spherical coordinates
          if(TypeGeometry(1:3)=='sph' .or. TypeGeometry(1:3)=='cyl') &
               XyzMin_D(1) = rBody
       end if
       RadiusMin = XyzMin_D(1)
       RadiusMax = XyzMax_D(1)
    end if

    ! Set default MaxBoundary if it was not set by #FACEBOUNDARY command
    if(i_line_command("#FACEBOUNDARY", iSessionIn = 1) < 0)then
       ! Use all face based BCs by default for spherical geometry
       if(TypeGeometry(1:9) == 'spherical')   MaxBoundary = 6

       ! Use face based boundaries by default for cylindrical geometry 
       ! except for top and bottom 
       if(TypeGeometry == 'cylindrical') MaxBoundary = 4
    end if

    ! Make sure MinBoundary and MaxBoundary cover face only boundaries
    if(UseBody2) MinBoundary = min(Body2_, MinBoundary)   !^CFG IF SECONDBODY
    if(body1) then   
       MinBoundary = min(Body1_, MinBoundary)
       MaxBoundary = max(Body1_, MaxBoundary)
    end if
    if(UseExtraBoundary) then                        
       MinBoundary = min(ExtraBc_, MinBoundary)
       MaxBoundary = max(ExtraBc_, MaxBoundary)
    end if
    MaxBoundary = min(MaxBoundary, 6)
    MinBoundary = max(MinBoundary, body2_)

    if(index(TypeGeometry,'_genr') < 1) call set_gen_radial_grid

    ! Initialize BATL
    if(TypeGeometry(1:9)=='spherical') then
       TypeGeometryBatl = 'rlonlat'//TypeGeometry(10:20)
    else
       TypeGeometryBatl = TypeGeometry
    end if

    call init_batl(XyzMin_D(1:nDim), XyzMax_D(1:nDim), MaxBlock, &
         TypeGeometryBatl, TypeBc_I(1:2*nDim-1:2) == 'periodic', &
         proc_dims(1:nDim), UseRadiusIn=.true., UseDegreeIn=.false.,&
         RgenIn_I = exp(LogRGen_I), UseUniformAxisIn=UseUniformAxis,&
         user_amr_geometry=user_specify_refinement)

    if(IsLogRadius .or. IsGenRadius)then
       ! Overwrite radial coordinates if necessary
       XyzMin_D(1) = CoordMin_D(1)
       XyzMax_D(1) = CoordMax_D(1)
    end if

    ! Fix grid size in ignored directions
    if(nDim == 1)then
       y1 = -0.5; XyzMin_D(2) = -0.5
       y2 = +0.5; XyzMax_D(2) = +0.5
    end if
    if(nDim < 3)then
       z1 = -0.5; XyzMin_D(3) = -0.5
       z2 = +0.5; XyzMax_D(3) = +0.5
    end if

  end subroutine correct_grid_geometry

  !============================================================================

  subroutine correct_plot_range

    use ModGeometry, ONLY: XyzMin_D, XyzMax_D, nIJK_D
    use ModParallel, ONLY: proc_dims
    use BATL_lib,    ONLY: radius_to_gen
    use ModIO

    implicit none

    integer :: iFile, iDim, iMin, iMax
    real    :: Ratio, PlotRes_D(3)

    real    :: CellSizeMax_D(3)
    real    :: SmallSize_D(3)   ! Used for 2D plot areas

    character(len=*), parameter:: NameSub = 'correct_plot_range'
    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    if(DoTestMe)write(*,*) NameSub,' XyzMin_D, XyzMax_D=', XyzMin_D, XyzMax_D

    ! Largest cell size and a much smaller distance for 2D cuts
    CellSizeMax_D = (XyzMax_D - XyzMin_D)/(nIJK_D*proc_dims)
    SmallSize_D   = cTiny*CellSizeMax_D

    if(DoTestMe)write(*,*)NameSub,' CellSizeMax_D=',CellSizeMax_D

    PLOTFILELOOP: do iFile = Plot_+1, Plot_ + nPlotFile

       plot_area = plot_type(iFile)(1:3)

       if(DoTestMe)write(*,*)'iFile, plot_area=',iFile, plot_area

       !Don't do anything to plot_type for specific plots.
       if(plot_area == 'lcb')CYCLE

       ! Fix plot range for sph, x=0, y=0, z=0 areas
       select case(plot_area)
       case('sph')
          if(IsCartesianGrid)then
             plot_dx(1,ifile) = 1.0    ! set to match write_plot_sph
             plot_dx(2:3,ifile) = 1.0  ! angular resolution in degrees
             plot_range(2,ifile)= plot_range(1,ifile) + 1.e-4 !so that R/=0
             plot_range(3,ifile)= 0.   - 0.5*plot_dx(2,ifile)
             plot_range(4,ifile)= 90.0 + 0.5*plot_dx(2,ifile)
             plot_range(5,ifile)= 0.   - 0.5*plot_dx(3,ifile)
             plot_range(6,ifile)= 360.0- 0.5*plot_dx(3,ifile)
          elseif(IsRLonLat)then
             plot_dx(1,ifile) = -1.0
             if(IsLogRadius) plot_range(1,ifile) = log(plot_range(1,ifile))
             if(IsGenRadius) call radius_to_gen(plot_range(1,ifile))
             plot_range(2,ifile)= plot_range(1,ifile) + 1.e-4 !so that R/=0
             do i=Phi_,Theta_
                plot_range(2*i-1,ifile) = XyzMin_D(i)
                plot_range(2*i,ifile)   = XyzMax_D(i)  
             end do
             plot_area='r=r' ! to disable the write_plot_sph routine
          else
             call stop_mpi(NameSub// &
                  ' Sph-plot is not implemented for geometry= '&
                  //TypeGeometry)
          end if

          ! There is nothing else to do for sph area
          CYCLE PLOTFILELOOP

       case('x=0')
          plot_range(1:5:2, iFile) = XyzMin_D
          plot_range(2:6:2, iFile) = XyzMax_D
          if( plot_form(iFile)=='tec' .or. IsCartesianGrid )then
             ! Limit plot range along x direction to be very small
             plot_range(1, iFile) = -SmallSize_D(x_)
             plot_range(2, iFile) = +SmallSize_D(x_)
             if(IsRlonLat) then
                plot_range(1, iFile) = XyzMin_D(x_)
                plot_range(2, iFile) = XyzMin_D(x_)+SmallSize_D(x_)
             end if
          else
             ! Limit Phi direction around cHalfPi
             plot_range(3, iFile) = cHalfPi - SmallSize_D(Phi_)
             plot_range(4, iFile) = cHalfPi + SmallSize_D(Phi_)
          end if

       case('y=0')
          plot_range(1:5:2, iFile) = XyzMin_D
          plot_range(2:6:2, iFile) = XyzMax_D
          if(plot_form(iFile) == 'idl' .and. (IsRLonLat.or.IsCylindrical)) then
             ! Limit plot range in Phi direction to be small around 180 degs
             plot_range(3, iFile) = cPi - SmallSize_D(y_)
             plot_range(4, iFile) = cPi + SmallSize_D(y_)
          else
             ! Limit plot range along y (or Phi) direction to be very small
             plot_range(3, iFile) = -SmallSize_D(y_)
             plot_range(4, iFile) = +SmallSize_D(y_)
          end if

       case('z=0', '2d_')
          ! Limit plot range along z direction to be very small
          plot_range(1:5:2, iFile) = XyzMin_D
          plot_range(2:6:2, iFile) = XyzMax_D
          plot_range(5, iFile) = -SmallSize_D(z_)
          plot_range(6, iFile) = +SmallSize_D(z_)

       case('1d_')
          ! Limit plot range along 2nd and 3rd directions to be very small
          plot_range(1, iFile) = XyzMin_D(1)
          plot_range(2, iFile) = XyzMax_D(1)
          plot_range(3:5:2, iFile) = -SmallSize_D(y_:z_)
          plot_range(4:6:2, iFile) = +SmallSize_D(y_:z_)
          
       case('3d_')
          plot_range(1:5:2, iFile) = XyzMin_D
          plot_range(2:6:2, iFile) = XyzMax_D

       end select

       ! Reduce plot range in ignored dimensions
       if(nJ == 1) plot_range(3,iFile) = -SmallSize_D(y_)
       if(nJ == 1) plot_range(4,iFile) = +SmallSize_D(y_)
       if(nK == 1) plot_range(5,iFile) = -SmallSize_D(z_)
       if(nK == 1) plot_range(6,iFile) = +SmallSize_D(z_)

       if(DoTestMe)write(*,*)'For file ',ifile-plot_,&
            ' original range   =',plot_range(:,ifile)

       plot_range(1:5:2, iFile) = max(plot_range(1:5:2, iFile), XyzMin_D)
       plot_range(2:6:2, iFile) = min(plot_range(2:6:2, iFile), XyzMax_D)

       ! Regular grid is not (yet) working in generalized coordinates
       ! because multiple pieces are used in the domain for x=0 and y=0 area
       if(plot_form(iFile) == 'idl' .and. .not.IsCartesianGrid) &
            plot_dx(1, iFile) = -1.0

       ! For plot_dx = 0.0 or -1.0 there is no need to adjust cut range
       if(plot_dx(1, iFile) <= cTiny)then

          plot_dx(:, iFile) = plot_dx(1, iFile) ! Define y and z

          CYCLE PLOTFILELOOP

       end if

       ! Make sure that plot resolution is a power of 2 fraction of cell size
       Ratio     = CellSizeMax_D(x_)/plot_dx(1, iFile)
       Ratio     = 2.0**nint(alog(Ratio)/alog(2.0))
       PlotRes_D = CellSizeMax_D / Ratio
       plot_dx(1:3, iFile) = PlotRes_D

       ! Make sure that plotting range is placed at an integer multiple of dx

       do iDim = 1, 3
          iMin = 2*iDim - 1; iMax = iMin+1

          ! Skip ignored dimensions of 2D and 1D cuts
          if(plot_range(iMax, iFile) - plot_range(iMin, iFile) &
               <= 1.5*PlotRes_D(iDim)) CYCLE

          ! Shift plot range slightly outward
          plot_range(iMin,iFile) = plot_range(iMin,iFile) - SmallSize_D(iDim)
          plot_range(iMax,iFile) = plot_range(iMax,iFile) + SmallSize_D(iDim)

          ! Round plot range to multiple of plot resolution
          plot_range(iMin:iMax, iFile) = XyzMin_D(iDim) + PlotRes_D(iDim)* &
               nint( (plot_range(iMin:iMax,iFile) &
               -      XyzMin_D(iDim))/PlotRes_D(iDim) )
       end do

       if(DoTestMe)write(*,*)'For file ',ifile-plot_,&
            ' adjusted range   =',plot_range(:,ifile)

    end do PLOTFILELOOP

  end subroutine correct_plot_range

  !===========================================================================
  subroutine set_extra_parameters

    use ModMultiFluid, ONLY: UseMultiIon
    use ModTimeStepControl, ONLY: UseTimeStepControl,TimeStepControlInit
    !--------------------------------------------------------------------------
    ! We need normalization for dt
    if(UseDtFixed)then
       if(.not.time_accurate)then
          if(iProc==0)then
             write(*,'(a)')NameSub//&
                  ' WARNING: UseDtFixed=T and not time accurate run',&
                  ' cannot be used together !!!'
             if (UseStrict) call stop_mpi('Correct PARAM.in')
             write(*,*)NameSub//' setting UseDtFixed to false...'
          end if
          UseDtFixed=.false.
       else ! fixed time step for time accurate
          DtFixed = DtFixedDim * Io2No_V(UnitT_)
          DtFixedOrig = DtFixed                   ! Store the initial setting
          Dt = DtFixed
          Cfl=1.0
       end if
    end if

    if(UseTimeStepControl)then
       ! Reduce initial time step / Cfl number. 
       ! The original values are stored in DtFixedOrig and CflOrig
       if(UseDtFixed)then
          DtFixed = TimeStepControlInit*DtFixed
       else
          Cfl     = TimeStepControlInit*Cfl
       end if
       TimeStepControlInit = 1.0
    else
       if(UseDtFixed)then
          DtFixed = DtFixedOrig
       else
          Cfl     = CflOrig
       end if
    end if

    DoOneCoarserLayer = .not. (nOrder>1 .and. &
         (UseTvdResChange .or. UseAccurateResChange))
    DoLimitMomentum = boris_correction .and. DoOneCoarserLayer

!!! momentum limiting fails for multiion: to be debugged
    if(UseMultiIon)DoLimitMomentum = .false.
!!!

  end subroutine set_extra_parameters

end subroutine MH_set_parameters
!=======================================================================
subroutine set_levels
  use ModAMR, ONLY: min_block_level, max_block_level
  use BATL_lib, ONLY: iTree_IA, MinLevel_, MaxLevel_

  iTree_IA(MinLevel_,:) = min_block_level
  iTree_IA(MaxLevel_,:) = max_block_level

end subroutine set_levels
