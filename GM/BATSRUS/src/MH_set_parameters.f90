!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine MH_set_parameters(TypeAction)

  ! Set input parameters for BATS-R-US

  use ModProcMH
  use ModMain
  use ModAdvance
  use ModB0, ONLY: UseB0Source, UseCurlB0, DoUpdateB0, DtUpdateB0, &
       read_b0_param, init_mod_b0
  use ModGeometry, ONLY: init_mod_geometry, TypeGeometry, nMirror_D, &
       x1,x2,y1,y2,z1,z2,XyzMin_D,XyzMax_D,RadiusMin,RadiusMax,&
       MinFaceBoundary, MaxFaceBoundary, CoordDimMin_D, CoordDimMax_D, &
       read_gen_radial_grid, set_gen_radial_grid, NameGridFile
  use ModNodes, ONLY: init_mod_nodes
  use ModImplicit
  use ModSemiImplicit, ONLY: read_semi_impl_param, init_mod_semi_impl
  use ModPartImplicit, ONLY: read_part_impl_param, init_mod_part_impl
  use ModImplHypre, ONLY: hypre_read_param
  use ModPhysics
  use ModProject
  use ModCT, ONLY : init_mod_ct, DoInitConstrainB
  use ModBlockData, ONLY: clean_block_data
  use BATL_lib, ONLY: read_amr_criteria, read_region_param, &
       DoCritAmr, DoAutoAmr, DoStrictAmr, BetaProlong,&
       init_mpi, IsCartesianGrid, IsCartesian, &
       IsRzGeometry, IsCylindrical, IsRLonLat, IsLogRadius, IsGenRadius
  use ModAMR
  use ModRaytrace
  use ModIO
  use CON_planet,       ONLY: read_planet_var, check_planet_var, NamePlanet
  use ModPlanetConst
  use CON_axes,         ONLY: init_axes, get_axes, &
       dLongitudeHgr, dLongitudeHgrDeg, dLongitudeHgi, dLongitudeHgiDeg
  use ModUtilities,     ONLY: fix_dir_name, check_dir, make_dir, DoFlush, &
       split_string, join_string
  use CON_planet,       ONLY: get_planet
  use ModTimeConvert,   ONLY: time_int_to_real, time_real_to_int
  use ModReadParam
  use ModMessagePass,   ONLY: DoOneCoarserLayer
  use ModFaceValue,     ONLY: &
       UseTvdResChange, UseAccurateResChange, &
       UseVolumeIntegral4, UseFaceIntegral4, UseLimiter4, nGUsed, &
       DoLimitMomentum, BetaLimiter, TypeLimiter, read_face_value_param, &
       TypeLimiter5, UseCweno,&
       iVarSmooth_V, iVarSmoothIndex_I, &
       StringLowOrderRegion, iRegionLowOrder_I
  use ModPartSteady,    ONLY: UsePartSteady, MinCheckVar, MaxCheckVar, &
       RelativeEps_V, AbsoluteEps_V
  use ModBoundaryCells, ONLY: init_mod_boundary_cells
  use ModPointImplicit, ONLY: read_point_implicit_param, UsePointImplicit
  use ModRestartFile,   ONLY: read_restart_parameters, init_mod_restart_file, &
       DoChangeRestartVariables, nVarRestart, UseRestartWithFullB, &
       NameRestartInDir, NameRestartOutDir
  use ModHallResist,    ONLY: &
       UseHallResist, read_hall_param
  use ModParticleFieldLine, ONLY: read_particle_line_param
  use ModHeatConduction, ONLY: read_heatconduction_param
  use ModHeatFluxCollisionless, ONLY: read_heatflux_param
  use ModRadDiffusion,   ONLY: read_rad_diffusion_param
  use ModResistivity, ONLY: UseResistivity, &
       read_resistivity_param, init_mod_resistivity
  use ModMultiFluid, ONLY: MassIon_I,ChargeIon_I,nIonFluid, iFluid, &
       DoConserveNeutrals, DoOhNeutralBc, &
       uBcFactor_I, RhoBcFactor_I, RhoNeutralsISW_dim, &
       UxNeutralsISW_dim, UyNeutralsISW_dim, UzNeutralsISW_dim, &
       TNeutralsISW_dim, mProtonMass

  use ModMultiIon, ONLY: multi_ion_set_parameters
  use ModSolarwind, ONLY: UseSolarwindFile, read_solar_wind_param, &
       read_solar_wind_file, normalize_solar_wind_data
  use ModSatelliteFile, ONLY: nSatellite, &
       read_satellite_parameters, read_satellite_input_files
  use ModGroundMagPerturb, ONLY: read_magperturb_param, init_mod_magperturb
  use ModFaceFlux, ONLY: face_flux_set_parameters, TypeFluxNeutral, &
       UseClimit, UsePoleDiffusion, DoBurgers
  use ModLookupTable,     ONLY: read_lookup_table_param
  use ModIeCoupling,      ONLY: read_ie_velocity_param
  use ModTimeStepControl, ONLY: read_time_step_control_param
  use ModLaserHeating,    ONLY: read_laser_heating_param
  use ModLocalTimeStep,   ONLY: read_localstep_param
  use ModIoUnit, ONLY: io_unit_new
  use ModNumConst, ONLY: cDegToRad
  use ModSort, ONLY: sort_quick

  use ModViscosity, ONLY: UseViscosity, viscosity_read_param, viscosity_init
  use ModPIC, ONLY: pic_read_param, pic_init_region, UsePic
  use ModFaceBoundary, ONLY: read_face_boundary_param
  !CORONA SPECIFIC PARAMETERS
  use EEE_ModMain, ONLY: EEE_set_parameters
  use ModMagnetogram, ONLY: set_parameters_magnetogram, &
       read_magnetogram_file, read_potential_field,     &
       read_new_magnetogram_file, read_new_potential_field
  use ModExpansionFactors,ONLY: NameModelSW, CoronalT0Dim, read_wsa_coeff
  use ModCoronalHeating,  ONLY: read_corona_heating, &
       init_coronal_heating, UseCoronalHeating, DoOpenClosedHeat
  use ModFieldLineThread, ONLY: read_threads 
  use ModThreadedLC,      ONLY: init_threaded_lc, read_threaded_bc
  use ModRadiativeCooling,ONLY: UseRadCooling,&
       read_modified_cooling, check_cooling_param, read_chromosphere
  use ModCoarseAxis, ONLY: read_coarse_axis_param
  use ModWaves, ONLY: read_waves_param, check_waves
  use ModLdem, ONLY: UseLdem, NameLdemFile, iRadiusLdem, read_ldem

  use ModUser, ONLY: NameUserModule, VersionUserModule
  use ModUserInterface ! user_read_inputs, user_init_session
  use ModConserveFlux, ONLY: DoConserveFlux
  use ModVarIndexes, ONLY: UseMultiSpecies, MassSpecies_V, SpeciesFirst_, &
       SpeciesLast_
  use BATL_lib, ONLY: Dim2_, Dim3_, &
       create_grid, set_high_geometry, get_region_indexes, &
       rRound0, rRound1
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
  logical :: IsMirrorX,  IsMirrorY,  IsMirrorZ

  ! The name of the command
  character (len=lStringLine) :: NameCommand, StringLine

  ! Temporary variables
  logical :: DoEcho=.false.
  integer :: nVarRead=0, nVarEquationRead = 0
  character (len=lStringLine) :: NameEquationRead="?"
  logical :: IsReadNameVarRestart = .false.
  character(len=lStringLine) :: NameVarRestartRead=''

  character (len=50) :: plot_string,log_string
  character (len=3)  :: plot_area, plot_var
  character (len=2)  :: NameCompRead="??"
  integer :: MinBlockAll, nIJKRead_D(3), nRootRead_D(3)=1

  integer            :: TimingDepth=-1
  character (len=10) :: TimingStyle='cumu'

  ! Variables for checking/reading #STARTTIME command
  real(Real8_)       :: StartTimeCheck   = -1.0_Real8_
  real               :: tSimulationCheck = -1.0

  ! Variable for #UNIFORMAXIS
  logical:: UseUniformAxis = .true.

  ! Variables for checking the user module
  character (len=lStringLine) :: NameUserModuleRead='?'
  real                        :: VersionUserModuleRead=0.0
  integer :: iSession, iPlotFile, iVar, iDim

  integer :: iSpecies

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

  if(iSession > 1)then
     restart=.false.           ! restart in session 1 only
  end if

  if(DoReadSatelliteFiles)then
     call read_satellite_input_files
     DoReadSatelliteFiles = .false.
  end if

  select case(TypeAction)
  case('CHECK')
     if(iProc==0)write(*,*) NameSub,': CHECK iSession =',iSession

     ! Make output and check input directories
     if(iProc==0) call make_dir(NamePlotDir)
     if(iProc==0 .and. save_restart_file) call make_dir(NameRestartOutDir)
     if(iProc==0 .and. restart) call check_dir(NameRestartInDir)

     if(StartTimeCheck > 0.0 .and. tSimulationCheck > 0.0)then
        if(abs(StartTime+time_simulation - StartTimeCheck-tSimulationCheck)&
             > 0.001)then
           write(*,*)NameSub//' WARNING: '// &
                NameThisComp//'::StartTimeCheck+tSimulationCheck=', &
                StartTimeCheck + tSimulationCheck, &
                ' differs from CON::StartTime+tSimulation=', &
                StartTime + time_simulation,' !!!'
           if(UseStrict)then
              call stop_mpi('Fix #STARTTIME/#SETREALTIME commands in PARAM.in')
           else
              ! Fix iStartTime_I array
              call time_real_to_int(StartTime, iStartTime_I)
           end if
        end if
        StartTimeCheck   = -1.0
        tSimulationCheck = -1.0
     end if
     if(UseEndTime)then
        t_max = EndTime - StartTime
        nIter = -1
        if(IsStandAlone)then
           if(.not.time_accurate)call CON_stop( &
                '#ENDTIME command cannot be used in steady-state mode')
           if(.not.IsLastRead) call CON_stop(&
                '#ENDTIME command can be used in the last session only')
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

     ! Initialize axes (coordinate transformation matrices)
     call init_axes(StartTime)
     if(NameThisComp == 'GM') then
        ! Set and obtain GM specific parameters from CON_planet and CON_axes
        call get_axes(Time_Simulation, MagAxisTiltGsmOut = ThetaTilt)
        call get_planet(DipoleStrengthOut = DipoleStrengthSi)
     end if

     if(MonopoleStrengthSi /= 0.0)then
        UseB0       = .true.
        UseB0Source = .false.
        UseCurlB0   = .false.
        DoUpdateB0  = .false.
        DtUpdateB0  = -1.0
     elseif(IsStandAlone .and. NameThisComp=='GM') then
        ! Check and set some planet variables (e.g. DoUpdateB0)
        call check_planet_var(iProc==0, time_accurate)

        if(body1)then
           call get_planet(UseRotationOut = UseRotatingBc)
        else
           UseRotatingBc = .false.
        end if

        ! Obtain some planet parameters
        if(.not.UseB0 .or. &
             (DipoleStrengthSi == 0.0 .and. .not. UseUserB0))then
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

     ! Set number and indexes of vector variables
     call init_vector_variables

     call correct_parameters

     ! initialize module variables
     call init_mod_advance
     call init_mod_geometry
     call init_mod_boundary_cells
     call init_mod_nodes
     if(UseB0)           call init_mod_b0
     if(UseRaytrace)     call init_mod_raytrace
     if(UseConstrainB)   call init_mod_ct
     if(UseImplicit)     call init_mod_part_impl
     if(UseSemiImplicit) call init_mod_semi_impl
     call init_mod_magperturb

     call get_region_indexes(StringLowOrderRegion, iRegionLowOrder_I)

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

     if(UseResistivity)call init_mod_resistivity
     if(UseViscosity) call viscosity_init

     if(UseMagnetogram)then
        if(UseNewMagnetogram)then
           if(i_line_command("#NEWMAGNETOGRAM") > 0)then
              call read_new_magnetogram_file(NamePlotDir, iProc, nProc, iComm)
           elseif(i_line_command("#READNEWPOTENTIALFIELD") > 0)then
              call read_new_potential_field(NamePlotDir, iProc, nProc, iComm)
           end if
           tMagnetogram = Time_Simulation
        end if
        if(i_line_command("#MAGNETOGRAM") > 0)then
           call read_magnetogram_file(NamePlotDir, iProc, nProc, iComm)
        elseif(i_line_command("#READPOTENTIALFIELD") > 0)then
           call read_potential_field(NamePlotDir, iProc, nProc, iComm)
        end if
     end if

     if(UsePic)  call pic_init_region

     if(UseEmpiricalSW .and. i_line_command("#EMPIRICALSW") > 0)&
          call set_empirical_model(NameModelSW, BodyTDim_I(IonFirst_))

     if(UseCoronalHeating)call init_coronal_heating
     call check_cooling_param

     if(UseFieldLineThreads.and.iSession==1)call init_threaded_lc
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

     select case(NameCommand)

     case("#BATLPROLONG")
        call read_var('BetaProlong', BetaProlong )

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
        call read_var('DoEcho', DoEcho)
        if(iProc==0)call read_echo_set(DoEcho)

     case("#FLUSH")
        call read_var('DoFlush', DoFlush)

     case("#TESTINFO")
        call read_var('DoWriteCallSequence', DoWriteCallSequence)

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
        do i = 1, 2*nDim
           call read_var('TypeBc', TypeBc_I(i))  
        end do

     case("#TIMESTEPPING", "#RUNGEKUTTA", "#RK")
        call read_var('nStage',  nStage)
        call read_var('CflExpl', Cfl)
        CflOrig = Cfl
        ExplCfl = Cfl
        UseHalfStep = NameCommand == "#TIMESTEPPING" .and. nStage <= 2

     case('#LOCALTIMESTEP')
        call read_localstep_param(NameCommand, iSession)

     case("#FIXEDTIMESTEP", "#TIMESTEPLIMIT")
        call               read_var('UseDtLimit', UseDtFixed)
        if(UseDtFixed)call read_var('DtLimitDim', DtFixedDim)

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

     case("#IMPLICIT", &
          "#IMPLCRITERIA", "#IMPLICITCRITERIA", "#STEPPINGCRITERIA", &
          "#PARTIMPL", "#PARTIMPLICIT",     &
          "#IMPLSCHEME", "#IMPLICITSCHEME", &
          "#IMPLSTEP", "#IMPLICITSTEP",     &
          "#IMPLCHECK", "#IMPLICITCHECK",   &
          "#IMPLENERGY", "#IMPLICITENERGY", &
          "#NEWTON", "#JACOBIAN", "#PRECONDITIONER", &
          "#KRYLOV", "#KRYLOVSIZE")
        call read_part_impl_param(NameCommand)           

     case("#SEMIIMPL", "#SEMIIMPLICIT", &
          "#SEMICOEFF", "#SEMIIMPLCOEFF", "#SEMIIMPLICITCOEFF", &
          "#SEMIPRECOND", "#SEMIPRECONDITIONER",&
          "#SEMIKRYLOV", "#SEMIKRYLOVSIZE","#SEMIIMPLICITSTABLE")
        call read_semi_impl_param(NameCommand)

     case("#HYPRE")
        call hypre_read_param

     case("#PIC", "#PICREGION", "#PICUNIT", "#PICCOUPLE", "#PICBALANCE")
        call pic_read_param(NameCommand)

     case("#VISCOSITY", "#VISCOSITYREGION","#ARTIFICIALVISCOSITY")
        call viscosity_read_param(NameCommand)

     case("#RESISTIVITY", "#RESISTIVITYOPTIONS", &
          "#RESISTIVITYREGION", "#RESISTIVEREGION",&
          '#MESSAGEPASSRESISTIVITY')
        call read_resistivity_param(NameCommand)

     case("#HALLRESISTIVITY", "#HALLREGION", "#BIERMANNBATTERY")
        call read_hall_param(NameCommand)

     case("#MINIMUMDENSITY")
        do iFluid = 1, nFluid
           call read_var('RhoMinDim', RhoMinDim_I(iFluid))
        end do

     case("#MINIMUMPRESSURE")
        do iFluid = 1, nFluid
           call read_var('pMinDim', pMinDim_I(iFluid))
        end do
        if(UseElectronPressure) call read_var('PeMinDim', PeMinDim)

     case("#MINIMUMTEMPERATURE")
        do iFluid = 1, nFluid
           call read_var('TMinDim', TMinDim_I(iFluid))
        end do
        if(UseElectronPressure) call read_var('TeMinDim', TeMinDim)

     case("#ELECTRONPRESSURE")
        call read_var('PeMinSi', PeMinSi)

     case("#ELECTRONENTROPY")
        call read_var('UseElectronEntropy', UseElectronEntropy)

     case("#ANISOTROPICPRESSURE")
        call read_var('UseConstantTau', UseConstantTau)
        if(UseConstantTau) &
             call read_var('TauInstabilitySi', TauInstabilitySi)
        call read_var('TauGlobalSi', TauGlobalSi)

     case("#EXTRAINTERNALENERGY")
        call read_var('ExtraEintMinSi', ExtraEintMinSi)

     case("#RADIATION", "#HEATFLUXLIMITER", "#ACCURATERADIATION")
        call read_rad_diffusion_param(NameCommand)

     case("#HEATCONDUCTION", "#WEAKFIELDCONDUCTION", "#IONHEATCONDUCTION")
        call read_heatconduction_param(NameCommand)

     case("#HEATFLUXREGION", "#HEATFLUXCOLLISIONLESS")
        call read_heatflux_param(NameCommand)

     case("#PARTICLELINE")
        call read_particle_line_param(NameCommand)

     case("#SAVELOGFILE")
        call read_var('DoSaveLogfile',save_logfile)
        if(save_logfile)then
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

           ! Recognize coordinate system names
           if (index(log_string,'GEO') > 0) TypeCoordPlot_I(logfile_) = 'GEO'
           if (index(log_string,'GSE') > 0) TypeCoordPlot_I(logfile_) = 'GSE'
           if (index(log_string,'GSM') > 0) TypeCoordPlot_I(logfile_) = 'GSM'
           if (index(log_string,'MAG') > 0) TypeCoordPlot_I(logfile_) = 'MAG'
           if (index(log_string,'SMG') > 0) TypeCoordPlot_I(logfile_) = 'SMG'
           if (index(log_string,'HGR') > 0) TypeCoordPlot_I(logfile_) = 'HGR'
           if (index(log_string,'HGI') > 0) TypeCoordPlot_I(logfile_) = 'HGI'
           if (index(log_string,'HGC') > 0) TypeCoordPlot_I(logfile_) = 'HGC'

           ! If any flux variables are used - input a list of radii
           ! at which to calculate the flux
           if (index(log_vars,'flx')>0) &
                call read_var('log_R_str',log_R_str)

        end if

     case("#SAVEINITIAL")
        call read_var('DoSaveInitial',DoSaveInitial)

     case("#SAVEPLOT")
        call read_var('nPlotFile', nPlotFile)
        nfile = max(nfile, plot_ + nPlotFile)
        if (nFile > MaxFile .or. nPlotFile > MaxPlotFile) call stop_mpi(&
             'The number of ouput files is too large in #SAVEPLOT:'&
             //' nPlotFile > MaxPlotFile .or. nFile > MaxFile')

        do iFile = Plot_ + 1, Plot_ + nPlotFile
           call read_var('StringPlot',plot_string)

           ! Plotting frequency
           call read_var('DnSavePlot',dn_output(ifile))
           call read_var('DtSavePlot',dt_output(ifile))

           ! Default resolution (original AMR grid)
           plot_dx(:,ifile) = -1.0

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
           elseif (index(plot_string,'pnt')>0) then
              plot_area='pnt'
           elseif(index(plot_string,'lin')>0)then
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
              end do
           elseif (index(plot_string,'eqr')>0)then
              plot_area='eqr'
              call read_var('nRadius',   plot_range(1,ifile))
              call read_var('nLon',      plot_range(2,ifile))
              call read_var('RadiusMin', plot_range(3,ifile))
              call read_var('RadiusMax', plot_range(4,ifile))
              plot_range(5,ifile) =   0.0
              plot_range(6,ifile) = 360.0
           elseif (index(plot_string,'eqb')>0)then
              plot_area='eqb'
              call read_var('nRadius',   plot_range(1,ifile))
              call read_var('nLon',      plot_range(2,ifile))
              call read_var('RadiusMin', plot_range(3,ifile))
              call read_var('RadiusMax', plot_range(4,ifile))
              call read_var('LongitudeMin', plot_range(5,ifile))
              call read_var('LongitudeMax', plot_range(6,ifile))
           elseif (index(plot_string,'ieb')>0)then
              plot_area='ieb'
           elseif (index(plot_string,'lcb')>0)then
              plot_area='lcb'
              call read_var('Radius', plot_range(1,ifile))
              call read_var('nLon',   plot_range(2,ifile))
           elseif (index(plot_string,'sph')>0)then
   	      plot_area='sph'
	      call read_var('Radius',plot_range(1,ifile))
           elseif (index(plot_string, 'shl')>0)then
              plot_area = 'shl'
              call read_var('TypeCoord', TypeCoordPlot_I(iFile))
              call read_var('rMin',   plot_range(1,iFile))
              call read_var('rMax',   plot_range(2,iFile))
              if (plot_range(1, iFile) /= plot_range(2,iFile)) &
                   call read_var('dR',   plot_dx(1,iFile))
              call read_var('LonMin', plot_range(3,iFile))
              call read_var('LonMax', plot_range(4,iFile))
              if (plot_range(3, iFile) /= plot_range(4,iFile)) &
                   call read_var('dLon', plot_dx(2,iFile))
              call read_var('LatMin', plot_range(5,iFile))
              call read_var('LatMax', plot_range(6,iFile))
              if (plot_range(5, iFile) /= plot_range(6,iFile)) &
                   call read_var('dLat', plot_dx(3,iFile))

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
              call read_var('xSizeImage', X_Size_Image(ifile))
              call read_var('ySizeImage', Y_Size_Image(ifile))
              ! read the number of pixels
              call read_var('nPixX', n_Pix_X(ifile))            
              call read_var('nPixY', n_Pix_Y(ifile))
           elseif(index(plot_string,'buf')>0)then
              plot_area='buf'
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
           if(index(plot_string,'idl') > 0)then
              plot_form(ifile)='idl'
              if (       plot_area /= 'sph' &
                   .and. plot_area /= 'shl' &
                   .and. plot_area /= 'los' &
                   .and. plot_area /= 'rfr' &
                   .and. plot_area /= 'lin' &
                   .and. plot_area /= 'eqr' &
                   .and. plot_area /= 'eqb' &
                   .and. plot_area /= 'buf' &
                   ) call read_var('DxSavePlot',plot_dx(1,ifile))

              ! Extract the type of idl plot file: default is real4
              TypeFile_I(iFile) = 'real4' 
              if(index(plot_string,'idl_real8') > 0) &
                   TypeFile_I(iFile) = 'real8'
              if(index(plot_string,'idl_ascii') > 0) &
                   TypeFile_I(iFile) = 'ascii'
              if(index(plot_string,'idl_tec') > 0) &
                   TypeFile_I(iFile) = 'tec'
           elseif(index(plot_string, 'hdf') > 0) then
              ! With these values VisIt recognises the files as timesteps
              ! with the general defaults it does not. 
              IsPlotName_n = .true.
              IsPlotName_t = .false.
              IsPlotName_e = .false.
              plot_form(iFile)='hdf'
              TypeFile_I(iFile) = 'hdf5'
           elseif(index(plot_string,'tec')>0)then 
              plot_form(ifile)='tec'
              TypeFile_I(iFile) = 'tec'
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
           elseif(index(plot_string,'RAY')>0.or.index(plot_string,'ray')>0)then
              plot_var='ray'
              plot_dimensional(ifile) = index(plot_string,'RAY')>0
              if(DoMapEquatorRay)then
                 plot_vars(ifile)='bx by bz req1 phi1 req2 phi2 status blk'
              else
                 plot_vars(ifile)='bx by bz theta1 phi1 theta2 phi2 status blk'
              end if
              plot_pars(ifile)='rbody'
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
           elseif(index(plot_string,'eqb')>0)then
              plot_var ='eqb'
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

     case("#PLOTFILENAME")
        call read_var('NameMaxTimeUnit', NameMaxTimeUnit)

     case("#LOSPLOT")
        call read_var('UseLosSimple', UseLosSimple)

     case("#SAVELOGNAME")
        call read_var('IsLogName_n',IsLogName_n)
        call read_var('IsLogName_e',IsLogName_e)
        ! Set _n true if _e is false.
        if(.not.IsLogName_e) IsLogName_n=.true.

     case("#SAVEPLOTSAMR")
        call read_var('DoSavePlotsAmr',save_plots_amr)

     case("#SAVEBINARY")
        call read_var('DoSaveBinary',save_binary)

     case("#GRIDRESOLUTION","#GRIDLEVEL","#REGION","#AMRREGION")
        call read_region_param(NameCommand, UseStrictIn=UseStrict)

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

     case('#CONSERVEFLUX')
        call read_var('DoConserveFlux', DoConserveFlux)

     case("#SCHEME")
        if(iSession>1) nOrderOld = nOrder
        call read_var('nOrder'  ,nOrder)
        ! Set default value for nStage. Can be overwritten if desired.
        nStage = nOrder
        ! Use RK3 for MP5 scheme
        if(nOrder > 4) nStage = 3
        UseHalfStep = nStage <= 2

        call read_var('TypeFlux',FluxType, IsUpperCase=.true.)
        BetaLimiter = 1.0
        if(nOrder > 1 .and. FluxType /= "SIMPLE")then
           call read_var('TypeLimiter', TypeLimiter)
           if(TypeLimiter /= 'minmod') &
                call read_var('LimiterBeta', BetaLimiter)
        else
           TypeLimiter = "no"
        end if

        if(nOrder==5) then
           ! The commands below can be reset in #SCHEME5
           UseFDFaceFlux    = .true. 
           UseCweno         = .false.
           !HighResChange does not work for 1D.
           if(nJ > 1) UseHighResChange = .true. 
           UseHighOrderAMR  = .true. 

           DoConserveFlux = .false. 
           UseTvdReschange = .false. 
           UseAccurateResChange = .false.            
        endif

     case("#SCHEME4")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('UseVolumeIntegral4', UseVolumeIntegral4)
        call read_var('UseFaceIntegral4',   UseFaceIntegral4)
        call read_var('UseLimiter4',        UseLimiter4)
        ! There is no face integral in 1D
        if(nDim == 1) UseFaceIntegral4 = .false.

     case("#SCHEME5")
        ! DoInterpolateFlux may not be used anymore. 
        ! call read_var('DoInterpolateFlux', DoInterpolateFlux)

        ! If UseFDFaceFlux is true. Use ECHO scheme, which based on:
        ! L. Del Zanna, O. Zanotti, N. Bucciantini, P. Londrillo,&
        ! Astronomy and Astrophysics, 473 (2007), pp.11-30.
        call read_var('UseFDFaceFlux', UseFDFaceFlux)
        call read_Var('TypeLimiter5', TypeLimiter5, IsLowerCase=.true.)
        call read_var('UseHighResChange', UseHighResChange)
        call read_var('UseHighOrderAMR',UseHighOrderAMR)

        ! If it is not 'cweno', mp5 scheme will be used. 
        UseCweno = TypeLimiter5 == 'cweno'

        ! The following lines are related to cweno scheme, and it needs
        ! more tests. 
        !if(UseCweno) call read_var('UsePerVarLimiter', UsePerVarLimiter)
        ! if(UseCweno .and. .not. DoInterpolateFlux) then
        !    ! Density and velocity use density as smooth indicator. 
        !    ! Other variables use themselves.
        !    iVarSmooth_V(1:Uz_) = Rho_
        !    do iVar = Uz_+1, nVar
        !       iVarSmooth_V(iVar) = iVar
        !    enddo
        !    call sort_smooth_indicator
        ! endif

        if(UseFDFaceFlux) then
           DoConserveFlux   = .false. 
        endif

        if(UseHighResChange) then
           UseTvdReschange = .false. 
           UseAccurateResChange = .false. 
        endif

     case('#BURGERSEQUATION')
        call read_var('DoBurgers', DoBurgers)

     case('#LIMITER', '#RESCHANGE', '#RESOLUTIONCHANGE', '#TVDRESCHANGE', &
          '#LIMITPTOTAL', '#FLATTENING', '#LOWORDERREGION')
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
        call read_var('nOrderProlong', nOrderProlong)
        UseTvdResChange = .false.
        UseAccurateResChange = .false.

     case("#MESSAGEPASS","#OPTIMIZE")               
        call read_var('TypeMessagePass', optimize_message_pass)

     case('#CLIMIT')
        call face_flux_set_parameters(NameCommand)

     case("#BORIS")
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
        end if

     case("#BORISSIMPLE")
        call read_var('UseBorisSimple',UseBorisSimple)
        if(UseBorisSimple) then
           call read_var('BorisClightFactor',boris_cLIGHT_factor)
           boris_correction=.false.
        else
           boris_cLIGHT_factor = 1.0
        end if

     case("#DIVB")
        if(.not.UseB)CYCLE READPARAM
        call read_var('UseDivbSource'   ,UseDivbSource)   
        call read_var('UseDivbDiffusion',UseDivbDiffusion)
        call read_var('UseProjection'   ,UseProjection)
        call read_var('UseConstrainB'   ,UseConstrainB)

        if (UseProjection.and.UseConstrainB.and.iProc==0) &
             call stop_mpi('Do not use projection and constrain B together!')
        if (UseProjection.and.UseDivbSource.and.iProc==0) then
           write(*,'(a)')NameSub // &
                ' WARNING: using divbsource and projection together !!!'
           if (UseStrict) call stop_mpi('Correct PARAM.in')
        end if
        if (UseConstrainB.and.UseDivbSource.and.iProc==0) then
           write(*,'(a)')NameSub // &
                ' WARNING: using divbsource and constrain B together !!!'
           if (UseStrict) call stop_mpi('Correct PARAM.in')
        end if

        if (iProc==0 .and. nDim > 1 &
             .and..not.UseHyperbolicDivb &
             .and..not.UseDivbSource &
             .and..not.UseProjection &
             .and..not.UseConstrainB &
             .and..not.UseDivbDiffusion) then
           write(*,'(a)') NameSub // &
                ' WARNING: you should use some div B control method !!!'
           if (UseStrict) call stop_mpi('Correct PARAM.in!')
        end if
        ! Make sure that divbmax will be calculated
        DivbMax = -1.0
        ! reinitialize constrained transport if needed
        DoInitConstrainB = .true.

     case("#USEB0", "#DIVBSOURCE", "#USECURLB0", "#MONOPOLEB0", "#B0FACTOR")
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

     case("#PROJECTION")
        if(.not.UseB)CYCLE READPARAM
        call read_var('TypeProjectIter' ,proj_method)
        call read_var('TypeProjectStop' ,proj_typestop)
        call read_var('RelativeLimit'   ,proj_divbcoeff)
        call read_var('AbsoluteLimit'   ,proj_divbconst)
        call read_var('MaxMatvec'       ,proj_matvecmax)
        ! Make sure that DivbMax is recalculated
        DivbMax = -1.0

     case("#CORRECTP")
        call read_var('pRatioLow',Pratio_lo)
        call read_var('pRatioHigh',Pratio_hi)
        if(Pratio_lo>=Pratio_hi)&
             call stop_mpi(NameSub//' ERROR: Pratio_lo>=Pratio_hi')

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

     case("#UNIFORMSTATE")
        UseShockTube = .true.
        call split_string(NamePrimitiveVar, nVar, NamePrimitive_V, nVarRead, &
             UseArraySyntaxIn=.true.)
        do i=1,nVar
           call read_var(NamePrimitive_V(i), ShockLeftState_V(i))
        end do
        ShockRightState_V = ShockLeftState_V

     case("#SHOCKTUBE")
        UseShockTube = .true.
        call split_string(NamePrimitiveVar, nVar, NamePrimitive_V, nVarRead, &
             UseArraySyntaxIn=.true.)
        do i=1,nVar
           call read_var(NamePrimitive_V(i)//' left', ShockLeftState_V(i))
        end do
        do i=1,nVar
           call read_var(NamePrimitive_V(i)//' right', ShockRightState_V(i))
        end do

     case("#SHOCKPOSITION")
        call read_var('ShockPosition',ShockPosition)
        call read_var('ShockSlope',ShockSlope)

     case("#SOLARWINDFILE", "#UPSTREAM_INPUT_FILE", "#REFRESHSOLARWINDFILE")
        call read_solar_wind_param(NameCommand)
        DoReadSolarwindFile = UseSolarwindFile

     case("#RAYTRACE")
        call read_var('UseRaytrace',UseRaytrace)
        if(UseRaytrace)then
           call read_var('UseAccurateTrace', UseAccurateTrace)
           call read_var('DtExchangeRay',    DtExchangeRay)
           call read_var('DnRaytrace',       DnRaytrace)
        end if
     case("#RAYTRACELIMIT")
        call read_var('RayLengthMax', RayLengthMax)
     case("#RAYTRACEEQUATOR")
        call read_var('DoMapEquatorRay', DoMapEquatorRay)
     case("#IE")
        call read_var('DoTraceIE', DoTraceIE)
     case("#IECOUPLING")
        call read_ie_velocity_param
     case("#IMCOUPLING","#IM")
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

     case("#CHANGEVARIABLES")
        call read_var('DoChangeRestartVariables',DoChangeRestartVariables)

     case("#EQUATION")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('NameEquation',NameEquationRead)
        call read_var('nVar',        nVarEquationRead)
        if(NameEquationRead /= NameEquation .and. iProc==0 &
             .and. .not. DoChangeRestartVariables)then   
           write(*,'(a)')'BATSRUS was compiled with equation '// &
                NameEquation//' which is different from '// &
                NameEquationRead
           call stop_mpi(NameSub//' ERROR: incompatible equation names')

        end if
        if(nVarEquationRead /= nVar .and. iProc==0 .and. &
             .not. DoChangeRestartVariables)then
           write(*,'(a,i2,a,i2)')&
                'BATSRUS was compiled with nVar=',nVar, &
                ' which is different from nVarEquationRead=',nVarEquationRead
           call stop_mpi(NameSub//' ERROR: Incompatible number of variables')
        end if

     case("#NEWRESTART","#RESTARTINDIR","#RESTARTINFILE",&
          "#PRECISION","#BLOCKLEVELSRELOADED")
        if(.not.is_first_session())CYCLE READPARAM
        call read_restart_parameters(NameCommand)

     case("#SAVERESTART", "#RESTARTOUTDIR","#RESTARTOUTFILE")
        call read_restart_parameters(NameCommand)

     case('#RESTARTBLOCKDATA')
        call read_restart_parameters(NameCommand)

     case("#RESTARTVARIABLES")
        ! This reads the names of the variables saved in the input
        ! restart file. 
        call read_var('NameVarRestartRead', NameVarRestartRead)
        IsReadNameVarRestart = .true.

     case("#RESTARTWITHFULLB")
        UseRestartWithFullB = .true.

     case("#PLOTDIR")
        call read_var("NamePlotDir",NamePlotDir)
        call fix_dir_name(NamePlotDir)

     case("#SATELLITE")
        if(.not.is_first_session())CYCLE READPARAM
        call read_satellite_parameters(NameCommand)
        DoReadSatelliteFiles = nSatellite > 0

     case('#STEADYSTATESATELLITE')
        call read_satellite_parameters(NameCommand)

     case('#GEOMAGINDICES')
        call read_magperturb_param(NameCommand)
        nFile = max(nFile, indexfile_)

     case("#MAGNETOMETER")
        call read_magperturb_param(NameCommand)
        nFile = max(nFile, magfile_) 

     case("#MAGNETOMETERGRID")
        call read_magperturb_param(NameCommand)
        nFile = max(nFile, maggridfile_) 

     case("#GRIDGEOMETRY", "#GRIDGEOMETRYLIMIT")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('TypeGeometry', TypeGeometry, IsLowerCase=.true.)
        ! need to read in the general grid file      
        if(TypeGeometry == 'spherical_genr') then
           call read_var('NameGridFile',NameGridFile)
           call read_gen_radial_grid(NameGridFile)
        end if
        if(TypeGeometry == 'roundcube') then
           call read_var('rRound0',rRound0)
           call read_var('rRound1',rRound1)
        end if        
        if(NameCommand == '#GRIDGEOMETRYLIMIT')then
           do iDim = 1, nDim
              call read_var('CoordDimMin_D', CoordDimMin_D(iDim))
              call read_var('CoordDimMax_D', CoordDimMax_D(iDim))
           end do
           XyzMin_D = CoordDimMin_D
           XyzMax_D = CoordDimMax_D
           if(TypeGeometry(1:9) == 'spherical')then
              RadiusMin = XyzMin_D(1)
              RadiusMax = XyzMax_D(1)
              XyzMin_D(Dim2_) = XyzMin_D(Dim2_)*cDegToRad
              XyzMax_D(Dim2_) = XyzMax_D(Dim2_)*cDegToRad
              XyzMin_D(Dim3_) = XyzMin_D(Dim3_)*cDegToRad
              XyzMax_D(Dim3_) = XyzMax_D(Dim3_)*cDegToRad
           elseif(TypeGeometry(1:11) == 'cylindrical')then
              RadiusMin = XyzMin_D(1)
              RadiusMax = XyzMax_D(1)
              XyzMin_D(Dim2_) = XyzMin_D(Dim2_)*cDegToRad
              XyzMax_D(Dim2_) = XyzMax_D(Dim2_)*cDegToRad
           end if
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

     case('#COARSEAXIS')
        call read_coarse_axis_param

     case("#GRID")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('nRootBlockX', nRootRead_D(1)) 
        call read_var('nRootBlockY', nRootRead_D(2))
        call read_var('nRootBlockZ', nRootRead_D(3))

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
             abs(VersionUserModule-VersionUserModuleRead)>0.001 .and. &
             .not. DoChangeRestartVariables) then
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
           call stop_mpi('Change nI,nJ,nK with Config.pl -g and recompile!')
        end if
        call read_var('MinBlockALL',MinBlockAll)
        if(MinBlockAll > MaxBlock*nProc .and. iProc==0)then
           write(*,*)'MaxBlock*nProc=', MaxBlock*nProc
           call stop_mpi('Use more processors'//&
                ' or increase MaxBlock with Config.pl -g and recompile!')
        end if

     case("#AMRINITPHYSICS")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('nRefineLevelIC',nRefineLevelIC)

     case("#GAMMA")
        if(.not.is_first_session())CYCLE READPARAM
        do iFluid = 1, nFluid
           call read_var('Gamma_I',Gamma_I(iFluid))
        end do
        ! Derived values for fluids
        GammaMinus1_I    = Gamma_I - 1.0
        where(GammaMinus1_I /= 0.0)
           InvGammaMinus1_I = 1.0 / GammaMinus1_I
        elsewhere
           ! This should not be used (isothermal case)
           InvGammaMinus1_I = 1.5
        end where

        ! Isothermal case (for ions?)
        if(any(Gamma_I == 1.0))then
           UseNonConservative = .true.
           nConservCrit = 0
           if(allocated(TypeConservCrit_I)) deallocate(TypeConservCrit_I)
           if(iProc==0) &
                write(*,*) NameSub,': for gamma=1 UseNonConservative is set to TRUE'
        endif

        ! Scalar values for the first fluid for simpler code
        Gamma          = Gamma_I(1)
        GammaMinus1    = GammaMinus1_I(1)
        InvGammaMinus1 = InvGammaMinus1_I(1)

        if(UseElectronPressure)then
           call read_var('GammaElectron', GammaElectron)
           ! Derived values for electron
           GammaElectronMinus1    = GammaElectron - 1.0
           InvGammaElectronMinus1 = 1.0 / GammaElectronMinus1
        else
           ! Default values for electrons are the same as first fluid
           ! so ideal MHD works as expected
           GammaElectron          = Gamma
           GammaElectronMinus1    = GammaMinus1
           InvGammaElectronMinus1 = InvGammaMinus1
        end if

     case("#LOOKUPTABLE")
        call read_lookup_table_param

     case("#PLASMA")
        if(UseMultiSpecies) then
           do iSpecies = SpeciesFirst_, SpeciesLast_
              call read_var('MassSpecies', MassSpecies_V(iSpecies))
           enddo           
           do iSpecies = SpeciesFirst_, SpeciesLast_
              call read_var('ChargeSpecies', ChargeSpecies_I(iSpecies))
           enddo           
        else
           do iFluid = IonFirst_, nFluid
              call read_var('MassFluid', MassFluid_I(iFluid))
           end do
           MassIon_I = MassFluid_I(IonFirst_:IonLast_)
           do iFluid = 1, nIonFluid
              call read_var('ChargeIon', ChargeIon_I(iFluid))
           end do
        endif
        
        call read_var('ElectronTemperatureRatio', ElectronTemperatureRatio)
        ElectronPressureRatio = ElectronTemperatureRatio

        !averageioncharge is only useful when there is only one ion
        if(nIonFluid==1 .and. .not. UseMultiSpecies)then  
           AverageIonCharge = ChargeIon_I(1)
           ElectronPressureRatio = ElectronTemperatureRatio*AverageIonCharge
        end if

        PePerPtotal = ElectronPressureRatio/(1 + ElectronPressureRatio)

     case("#MULTISPECIES")
        call read_var('DoReplaceDensity', DoReplaceDensity)
        call read_var('SpeciesPercentCheck',SpeciesPercentCheck)

     case("#NEUTRALFLUID")
        call read_var('DoConserveNeutrals', DoConserveNeutrals)
        call read_var('TypeFluxNeutral',    TypeFluxNeutral)

     case("#MULTIION", "#MHDIONS", "#COLLISION", "#MULTIIONSTATE")
        call multi_ion_set_parameters(NameCommand)

     case('#USERBOUNDARY', '#EXTRABOUNDARY')
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('UseExtraBoundary', UseExtraBoundary)
        if(UseExtraBoundary)then
           call read_var('TypeBc_I(ExtraBc_)', TypeBc_I(ExtraBc_))      
        end if

     case("#FACEBOUNDARY")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('MinFaceBoundary', MinFaceBoundary)
        call read_var('MaxFaceBoundary', MaxFaceBoundary)

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

     case("#INNERBOUNDARY", "#POLARBOUNDARY", "#CPCPBOUNDARY", &
          "#MAGNETICINNERBOUNDARY")
        call read_face_boundary_param(NameCommand)

     case("#GRAVITY")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('UseGravity',UseGravity)
        if(UseGravity)then
           call read_var('iDirGravity',GravityDir)
           if(GravityDir /= 0) call read_var('GravitySi', GravitySi)
        end if

     case("#SECONDBODY")
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
           call read_var('UseBody2Orbit', UseBody2Orbit)
           if(UseBody2Orbit)then
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

     case('#PLANET','#MOON','#COMET','#IDEALAXES','#ROTATIONAXIS',&
          '#MAGNETICAXIS','#MAGNETICCENTER','#ROTATION','#DIPOLE', &
          '#NONDIPOLE','#UPDATEB0')

        call check_stand_alone
        if(.not.is_first_session())CYCLE READPARAM

        call read_planet_var(NameCommand)

     case("#ROTATEHGR")
        call check_stand_alone
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('dLongitudeHgr', dLongitudeHgrDeg)
        dLongitudeHgr = dLongitudeHgrDeg * cDegToRad

     case("#ROTATEHGI")
        call check_stand_alone
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('dLongitudeHgi', dLongitudeHgiDeg)
        dLongitudeHgi = dLongitudeHgiDeg * cDegToRad        

     case("#COORDSYSTEM","#COORDINATESYSTEM")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('TypeCoordSystem',TypeCoordSystem,IsUpperCase=.true.)
        select case(NameThisComp)
        case('GM')
           if (TypeCoordSystem == 'GEO') UseRotatingFrame = .true.
           if(TypeCoordSystem /= 'GSM' .and. TypeCoordSystem /= 'GSE' &
                .and. TypeCoordSystem /= 'GEO' ) &
                call stop_mpi(NameSub// &
                ' ERROR: cannot handle coordinate system '//TypeCoordSystem)
        case('IH','OH')
           select case(TypeCoordSystem)
           case('HGI')
              ! If rotating frame was on in the previous session then
              ! we need to transform from HGR/HGC to HGI system.
              ! Note: This only works if the twoo coordinate systems are aligned 
              ! at the initial time (i.e. HGR = HGC).
              if(UseRotatingFrame) iSignRotationIC = +1
              UseRotatingFrame = .false.
           case('HGC','HGR')
              UseRotatingFrame = .true.
           case default
              call stop_mpi(NameSub// &
                   ' ERROR: cannot handle coordinate system '&
                   //TypeCoordSystem)
           end select
        case('SC')
           select case(TypeCoordSystem)
           case('HGR','HGC')
              UseRotatingFrame = .true.
              if(iProc==0)then
                 write(*,*)NameSub//' setting .UseRotatingFrame = T'
              end if
           case('HGI')
              if(iProc==0) write(*,*) NameSub,&
                   ' WARNING: inertial SC is less accurate'
              UseRotatingFrame = .false.
           case default
              call stop_mpi(NameSub// &
                   ' ERROR: cannot handle coordinate system '&
                   //TypeCoordSystem)
           end select
        case('EE')
           select case(TypeCoordSystem)
           case('HGR','HGC')
              UseRotatingFrame = .true.
           case('HGI','GSM')
              UseRotatingFrame = .false.
           case default
              call stop_mpi(NameSub// &
                   ' ERROR: cannot handle coordinate system '&
                   //TypeCoordSystem)
           end select
        end select

     case("#ROTATINGINITIALCONDITION")
        call read_var('iSignRotationIC', iSignRotationIC)

     case("#NSTEP")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('nStep',n_step)

     case("#NPREVIOUS")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('nPrev',n_prev)
        call read_var('DtPrev',dt_prev)

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

     case("#TIMEEND", "#ENDTIME")
        UseEndTime = .true.
        call read_var('iYear'  ,iEndTime_I(1))
        call read_var('iMonth' ,iEndTime_I(2))
        call read_var('iDay'   ,iEndTime_I(3))
        call read_var('iHour'  ,iEndTime_I(4))
        call read_var('iMinute',iEndTime_I(5))
        call read_var('iSecond',iEndTime_I(6))
        iEndTime_I(7) = 0
        call time_int_to_real(iEndTime_I, EndTime)


     case("#TIMESIMULATION")
        if(.not.is_first_session())CYCLE READPARAM
        if(IsStandAlone)then
           call read_var('tSimulation',time_simulation)
        else
           call read_var('tSimulation',tSimulationCheck)
        end if

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
        if(NameThisComp == 'SC') &
             call stop_mpi(NameSub//' ERROR:'// &
             ' #HELIOBUFFERGRID command can be used in IH,OH components only')
        call read_var('rBuffMin',   BufferMin_D(BuffR_))
        call read_var('rBuffMax',   BufferMax_D(BuffR_))
        call read_var('nThetaBuff', nThetaBuff)
        call read_var('nPhiBuff',   nPhiBuff)

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

     case("#NEWMAGNETOGRAM", "#READNEWPOTENTIALFIELD")
        call read_var('UseNewMagnetogram', UseNewMagnetogram)
        if(UseNewMagnetogram)&
             call set_parameters_magnetogram(NameCommand)

     case("#SAVEPOTENTIALFIELD", "#B0GRID")
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
          "#HEATPARTITIONING", "#POYNTINGFLUX", "#KOLMOGOROV")
        call read_corona_heating(NameCommand)

     case("#OPENCLOSEDHEAT")
        call read_var('DoOpenClosedHeat', DoOpenClosedHeat)

     case("#RADIATIVECOOLING")
        call read_var('UseRadCooling', UseRadCooling)

     case("#CHROMOSPHERE")
        call read_chromosphere

     case("#TRANSITIONREGION")
        call read_modified_cooling

     case("#FIELDLINETHREAD")
        call read_threads(iSession)

     case("#THREADEDBC")
        call read_threaded_bc   

     case("#ADVECTWAVES", "#ALFVENWAVES", "#WAVEPRESSURE", &
          "#FREQUENCY", "#SPECTRUM", "#WAVEREFLECTION")
        call read_waves_param(NameCommand)

     case("#LASERPULSE", "#LASERBEAM", "#LASERBEAMS", "#LASERBEAMPROFILE", &
          "#LASERRAYTEST")
        call read_laser_heating_param(NameCommand)

     case("#CME", "#ARCH", "#TD99FLUXROPE", "#GL98FLUXROPE", "#SHEARFLOW", &
          "#CMS")
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
       UseRotatingFrame  = .false.
       UseRotatingBc     = .false.
       TypeCoordSystem   = 'HGI'
    case('SC','EE')
       UseRotatingFrame  = .true.
       UseRotatingBc     = .false.
       TypeCoordSystem   = 'HGR'
    case('GM')
       UseRotatingFrame  = .false.
       UseRotatingBc     = .true.
       TypeCoordSystem   = 'GSM'
    end select

    !Do not set B0 field in IH and OH
    if(NameThisComp/='IH'.and.NameThisComp/='OH')then
       UseB0=UseB
    else
       UseB0=.false.
    end if

    ! Do not update B0 except in GM
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

    nRefineLevelIC        = 0

    min_block_level =  0
    max_block_level = 99
    min_cell_dx =     0.
    max_cell_dx = 99999.

    DnAmr=-1
    DoAmr=.false.
    automatic_refinement = .false.

    nOrder = 2
    FluxType = 'RUSANOV'

    ! Default implicit parameters
    UseImplicit      = .false.
    ImplCritType     = 'dt'

    UseDivbSource   =  UseB .and. nDim > 1
    UseDivbDiffusion= .false.
    UseProjection   = .false.
    UseConstrainB   = .false.

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

    optimize_message_pass = 'all'

    plot_dimensional      = .true.

    restart           = .false.

    !\
    ! Give some "reasonable" default values
    !/

    DipoleStrengthSi = 0.0

    UseBody2 = .false.
    RBody2 =-1.0
    xBody2 = 0.0
    yBody2 = 0.0			
    zBody2 = 0.0			
    BdpBody2_D  = 0.0
    rCurrentsBody2 = 0.0
    RhoDimBody2 = 1.0    ! n/cc
    TDimBody2   = 10000.0! K

    MassIon_I = MassFluid_I(IonFirst_:IonLast_) ! Ion masses

    nAmrCriteria = 0
    call init_mod_amr(nAmrCriteria)

    !\
    ! Set component dependent defaults
    !/

    GravityDir=0
    if(allocated(TypeConservCrit_I)) deallocate(TypeConservCrit_I)

    select case(NameThisComp)
    case('SC','IH','OH','EE')
       ! Body parameters
       UseGravity=.true.
       body1      =.true.
       if(NameThisComp == 'EE') body1 = .false.
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

    end select

  end subroutine set_defaults

  !=========================================================================
  subroutine correct_parameters

    use ModMultiFluid, ONLY: UseMultiIon
    use ModWaves, ONLY: UseAlfvenWaves, UseWavePressure
    use ModRestartFile, ONLY: NameVarRestart_V

    ! option and module parameters
    character (len=40) :: Name
    real               :: Version
    logical            :: IsOn

    real    :: BetaProlongOrig = 0.0
    logical :: IsFirstCheck = .true.
    character(len=4) :: NameVarTemp_V(100) = ''
    !---------------------------------------------------------------------

    !\
    ! Check for some combinations of things that cannot be accepted as input
    !/
    if (iProc==0) write (*,*) ' '

    if(IsFirstCheck)then
       call correct_grid_geometry

       if( (.not.IsCartesian .or. MinFaceBoundary <= MaxFaceBoundary) &
            .and. (UseVolumeIntegral4 .or. UseFaceIntegral4))then
          if(iProc==0)then
             if(.not. IsCartesian) write(*,'(a)')NameSub//&
                  ': UseVolumeIntegral4/UseFaceIntegral4 are implemented ', &
                  'for Cartesian grid only!'

             if(MinFaceBoundary <= MaxFaceBoundary) write(*,'(a)')NameSub//&
                  ': UseVolumeIntegral4/UseFaceIntegral4 are implemented ', &
                  'for cell based boundaries only!'

             if (UseStrict) call stop_mpi('Correct PARAM.in!')
             write(*,*)NameSub//' Setting Use*Integral4 = .false.'
          end if
          UseVolumeIntegral4 = .false.
          UseFaceIntegral4   = .false.
       end if

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
    else if(nOrderOld /=5 .and. nOrder==5 .and. UseFDFaceFlux) then
       ! calculate high-order geometry coefficients
       call set_high_geometry(UseFDFaceFluxIn=UseFDFaceFlux)
       call create_grid
       nOrderOld = 5
    endif ! IsFirstCheck

    ! Get the number of used ghost cell layers
    select case(nOrder)
    case(1, 2)
       nGUsed = nOrder
    case(4)
       ! 4th order interpolation formula needs 2 ghost cells
       nGUsed = 2 
       ! Volume integral needs an extra ghost cell layer
       if(UseVolumeIntegral4) nGUsed = nGUsed + 1
       if(TypeLimiter /= 'no')then
          ! PPM limiter needs another 1 or 2 ghost cell layers
          nGUsed = nGUsed + 1
          if(UseLimiter4) nGUsed = nGUsed + 1
       end if
    case(5)
       ! MP5 scheme needs 3 ghost cell layers
       nGUsed = 3
    case default
       if(iProc==0)then
          write(*,*)'ERROR: nOrder=', nOrder
          call stop_mpi(NameSub//': Invalid value for nOrder')
       end if
    end select

    ! Check if there are enough ghost cells
    ! We could dynamically (re)allocate State_VGB etc. with nGUsed ?!
    if(nGUsed /= nG .and. iProc==0)then
       write(*,*)'The code is configured with nG=',nG,' ghost cell layers.'
       write(*,*)'The selected scheme requires nGUsed=',nGUsed,' layers!'
       if(nGUsed > nG)then
          write(*,*)'Either change settings or reconfigure and recompile!'
          call stop_mpi(NameSub//': insufficient number of ghost cells')
       end if
    end if

    ! This depends on the grid geometry set above
    call correct_plot_range

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
    case('ROE','Roe')
       FluxType='Roe'
       if(UseAlfvenWaves .or. UseWavePressure)then
          if(iProc==0) write(*,'(a)')NameSub // &
               'Wave transport and wave pressure do not work with ' // &
               trim(FluxType)
          if(UseStrict)&
               call stop_mpi('Correct PARAM.in')
          FluxType='Sokolov'
       end if
    case('ROEOLD','RoeOld')             
       FluxType='RoeOld'
    case('RUSANOV','TVDLF','Rusanov')
       FluxType='Rusanov'
    case('LINDE','HLLEL','Linde')
       FluxType='Linde'
    case('SOKOLOV','AW','Sokolov')
       FluxType='Sokolov'
    case('HLLD')
    case('GODUNOV','Godunov')
       FluxType='Godunov'
    case default
       if(iProc==0)then
          write(*,'(a)')NameSub // &
               'WARNING: unknown value for FluxType=' // trim(FluxType)
          if(UseStrict) &
               call stop_mpi('Correct PARAM.in!')
          write(*,*)'setting FluxType=Rusanov'
       end if
       FluxType='Rusanov'
    end select

    ! Set flux type for neutral fluids
    select case(TypeFluxNeutral)
    case('default')
       select case(FluxType)
       case('Rusanov', 'Linde', 'Sokolov', 'Godunov')
          TypeFluxNeutral = FluxType
       case default
          TypeFluxNeutral = 'Linde'
       end select
    case('RUSANOV','TVDLF','Rusanov')
       TypeFluxNeutral = 'Rusanov'
    case('LINDE','HLLE','Linde')
       TypeFluxNeutral = 'Linde'
    case('SOKOLOV','AW','Sokolov')
       TypeFluxNeutral = 'Sokolov'
    case('GODUNOV','Godunov')
       TypeFluxNeutral = 'Godunov'
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

    ! Check flux type selection for implicit
    select case(FluxTypeImpl)
    case('default')
       FluxTypeImpl = FluxType
    case('ROE','Roe')
       FluxTypeImpl='Roe'
    case('ROEOLD','RoeOld')
       FluxTypeImpl='RoeOld'
    case('RUSANOV','TVDLF','Rusanov')
       FluxTypeImpl='Rusanov'
    case('LINDE','HLLEL','Linde')
       FluxTypeImpl='Linde'
    case('SOKOLOV','AW','Sokolov')
       FluxTypeImpl='Sokolov'
    case('HLLD')
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
    end select

    ! Check flux types
    if( (FluxType(1:3)=='Roe' .or. FluxTypeImpl(1:3)=='Roe' .or. &
         FluxType=='HLLD' .or.  FluxTypeImpl=='HLLD') &
         .and. .not.UseB)then
       if (iProc == 0) then
          write(*,'(a)')NameSub//&
               ' WARNING: HLLD/Roe(Old) flux is only implemented for MHD !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' Setting TypeFlux(Impl) = Linde'
       end if
       if(FluxType(1:3)=='Roe' .or. FluxType=='HLLD') &
            FluxType     = 'Linde'
       if(FluxTypeImpl(1:3)=='Roe' .or. FluxTypeImpl=='HLLD') &
            FluxTypeImpl = 'Linde'
    end if

    if((FluxType=='Godunov' .or. FluxTypeImpl=='Godunov') &
         .and. UseB)then
       if (iProc == 0) then
          write(*,'(a)')NameSub//&
               ' WARNING: Godunov flux is only implemented for HD !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' Setting TypeFlux(Impl) = Linde'
       end if
       if(FluxType=='Godunov')     FluxType     = 'Linde'
       if(FluxTypeImpl=='Godunov') FluxTypeImpl = 'Linde'
    end if

    UseBdf2 = nStage > 1 .and. time_accurate

    ! Make sure periodic boundary conditions are symmetric
    if(any(TypeBc_I(1:2)=='periodic')) TypeBc_I(1:2)='periodic'
    if(any(TypeBc_I(3:4)=='periodic')) TypeBc_I(3:4)='periodic'
    if(any(TypeBc_I(5:6)=='periodic')) TypeBc_I(5:6)='periodic'

    if(UseConstrainB .and. .not.time_accurate)then
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
    end if

    if (UseHallResist .or. UseResistivity .or. UseViscosity) &
         optimize_message_pass = 'all'

    if (UseRadDiffusion .and. (UseFullImplicit .or. UseSemiImplicit)) &
         optimize_message_pass = 'all'

    !Check for magnetogram

    if(UseEmpiricalSW.and..not.UseMagnetogram)&
         call stop_mpi(&
         'Empirical Solar Wind model requires magnetogram')

    if(DoOpenClosedHeat.and.(.not.UseMagnetogram.and.&
         (iSession==1.and.i_line_command('#PFSSM')<0)))&
         call stop_mpi(&
         'The heating in the closed field region requires magnetogram')

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

    ! Accurate res change algorithm and 4th order finite volume scheme
    ! both need corners and edges
    if (UseAccurateResChange .or. nOrder == 4) &
         optimize_message_pass = 'all'

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
       if(UseProjection)call stop_mpi(&
            'Only Cartesian works with projection')
       if(UseConstrainB)call stop_mpi(&
            'Only Cartesian works with constrain B')
       if(UseDivBDiffusion)call stop_mpi(&
            'Only Cartesian works with divB diffusion')
    end if

    if(UseHyperbolicDivb)then
       if(UseDivbDiffusion .or. UseConstrainB .or. UseProjection )then
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
       if(boris_correction)then
          if(Cfl > 0.65) then
             write(*,'(a)')NameSub// &
                  ' WARNING: CFL number above 0.65 may be unstable'// &
                  ' for local timestepping with Boris correction !!!'
             if (UseStrict) call stop_mpi('Correct PARAM.in!')
          end if
       else
          if(cfl > 0.85)then
             write(*,'(a)')NameSub// &
                  ' WARNING: CFL number above 0.85 may be unstable'// &
                  ' for local timestepping !!!'
             if (UseStrict) call stop_mpi('Correct PARAM.in!')
          end if
       end if
    end if

    ! Boris correction checks
    if((FluxType(1:3)=='Roe' .or. FluxTypeImpl(1:3)=='Roe') &
         .and. boris_correction)then
       if (iProc == 0) then
          write(*,'(a)')NameSub//&
               ' WARNING: Boris correction not available for Roe flux !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' Setting boris_correction = .false.'
       end if
       boris_correction = .false.
       boris_cLIGHT_factor = 1.0
    end if

    if(boris_correction .and. UseHallResist .and. &
         TypeSemiImplicit /= 'resisthall' .and. TypeSemiImplicit /= 'resistivity')then
       if(iProc==0)then
          write(*,'(a)') NameSub//&
               ' WARNING: Boris correction only works with semi-implicit Hall MHD!!!'
          if (UseStrict) call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' Setting boris_correction = .false.'
       end if
       boris_correction = .false.
       boris_cLIGHT_factor = 1.0
    end if

    ! Check parameters for implicit
    if(UseImplicit .and. UseSemiImplicit .and. nStage > 1)then
       if(iProc==0)write(*,'(a)') NameSub//&
            ' WARNING: Full/part implicit + semi implicit schemes with'//&
            ' nStage=2 (BDF2) scheme is inconsistent if there is overlap!!!'
       if(UseStrict)then
          if(iProc == 0)write(*,*)NameSub//' setting nStage=1'
          nStage = 1
       endif
    end if

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

    !Finish checks for implicit

    !Set min_block_level and max_block_level for #AMRRESOLUTION in session 1
    if(i_line_command("#AMRRESOLUTION", iSessionIn = 1) > 0)then
       local_root_dx = (XyzMax_D(x_)-XyzMin_D(x_))/real(nRootRead_D(1)*nI)
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

    if(NameThisComp == 'SC' .and. TypeCoordSystem == 'HGI')then
       if(iProc == 0)then
          write(*,'(a)') NameSub//&
               ' WARNING: SC only works with rotating frame!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*) NameSub//' setting TypeCoordSystem = HGC'
       end if
       TypeCoordSystem  = 'HGC'
       iSignRotationIC = 0
       UseRotatingFrame = .true.
    end if

    if(TypeCoordSystem == 'HGI' .and. NameThisComp /= 'OH' &
         .and. .not.time_accurate)then
       if(iProc == 0)then
          write(*,'(a)') NameSub//&
               ' WARNING: there is no steady state solution in HGI system!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*) NameSub//' setting TypeCoordSystem = HGC'
       end if
       TypeCoordSystem  = 'HGC'
       iSignRotationIC = 0
       UseRotatingFrame = .true.
    end if

    ! Update check does not work with Runge-Kutta schemes
    ! because the final update is a linear combination of all stages.
    if(.not.UseHalfStep) UseUpdateCheck = .false.

    ! Use first order prolongation for the first stage of high 
    ! resolution change.
    if(UseHighResChange) nOrderProlong = 1

    ! If the state variables are changed when restarting, the names of the
    ! state variables saved in the restart file must be specified, so they
    ! can be properly copied into State_VGB.
    if(DoChangeRestartVariables .and. .not. IsReadNameVarRestart &
         .and. NameEquationRead /= NameEquation) then
       if(iProc == 0)then
          write(*,'(a)') NameSub//&
               ' : Could not find  #RESTARTVARIABLES command, needed '// &
               'when changing state variables at restart.'
       end if
       call stop_mpi('Correct PARAM.in')
    end if

    ! Check that the number of variables listed in #RESTARTVARIABLES
    ! matches the number appearing in the #EQUATION command
    ! (this check is useful if the #RESTARTVARIABLES command is added
    ! manually to older restart header files, but might become redundant 
    ! in the future).
    if(IsReadNameVarRestart) then
       call split_string(NameVarRestartRead,100, NameVarTemp_V, &
            nVarRestart)
       if (nVarRestart /= nVarEquationRead) then
          write(*,*)'Number of variables in #EQUATION command is different'//&
               ' than the number of variables listed in #RESTARTVARIABLES.'
          call stop_mpi(NameSub//': Correct PARAM.in file!')
       else
          ! Save array of variable names
          if (allocated(NameVarRestart_V)) deallocate(NameVarRestart_V)
          allocate(NameVarRestart_V(nVarRestart))
          NameVarRestart_V = NameVarTemp_V
       end if
    end if

    UseHighOrderAMR = UseHighOrderAMR .and. DoAmr
    IsFirstCheck = .false.

  end subroutine correct_parameters

  !===========================================================================
  subroutine correct_grid_geometry

    use ModGeometry, ONLY: LogRGen_I
    use BATL_lib, ONLY: init_batl, CoordMin_D, CoordMax_D, &
         IsRotatedCartesian

    character(len=20):: TypeGeometryBatl

    character(len=*), parameter:: NameSub='correct_grid_geometry'
    !-----------------------------------------------------------------------

    if(i_line_command("#GRID", iSessionIn = 1) < 0) &
         call stop_mpi(NameSub // &
         ' #GRID command must be specified in the first session!')

    if(product(nRootRead_D) > MaxBlock*nProc .and. iProc==0)then
       write(*,*)'Not enough grid blocks allocated for root blocks'
       write(*,*)'Number of root blocks=',product(nRootRead_D)
       write(*,*)'MaxBlock, nProc, MaxBlock*nProc=', &
            MaxBlock, nProc, MaxBlock*nProc
       call stop_mpi(NameSub//': insufficient number of grid blocks')
    end if

    ! Set XyzMin_D, XyzMax_D based on 
    ! #GRID, #GRIDGEOMETRY, and #LIMITRADIUS
    ! #GRIDGEOMETRYLIMIT already sets XyzMin_D, XyzMax_D so that it does not
    ! have to be reset here
    if(.not.i_line_command("#GRIDGEOMETRYLIMIT", iSessionIn = 1) > 0) then
       select case(TypeGeometry)
       case('cartesian' ,'rotatedcartesian')
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
          XyzMax_D = (/sqrt(max(x1**2,x2**2) + max(y1**2,y2**2)), cTwoPi, z2/)
       case('roundcube')
          if(rRound0 > rRound1)then
             ! Cartesian outside, so use x1..z2
             XyzMin_D = (/x1, y1, z1/)
             XyzMax_D = (/x2, y2, z2/)
          else
             ! Round outside, so fit this inside x1..z2
             if(nDim==2) XyzMax_D = &
                  min(abs(x1), abs(x2), abs(y1), abs(y2)) &
                  /sqrt(2.0)
             if(nDim==3) XyzMax_D = &
                  min(abs(x1), abs(x2), abs(y1), abs(y2), abs(z1), abs(z2)) &
                  /sqrt(3.0)
             XyzMin_D = -XyzMax_D
          end if

       case default
          call stop_mpi(NameSub//': unknown TypeGeometry='//TypeGeometry)
       end select

       if(i_line_command("#LIMITRADIUS", iSessionIn = 1) > 0) then 
          XyzMin_D(1) = RadiusMin
          XyzMax_D(1) = RadiusMax
          if(TypeGeometry == 'roundcube' .and. rRound1 > rRound0) &
               XyzMax_D = RadiusMax/sqrt(real(nDim))
       else
          if(Body1 .and. rBody > 0.0)then
             ! Set inner boundary to match rBody for spherical coordinates
             if(TypeGeometry(1:3)=='sph' .or. TypeGeometry(1:3)=='cyl') &
                  XyzMin_D(1) = rBody
          end if
          RadiusMin = XyzMin_D(1)
          RadiusMax = XyzMax_D(1)
       end if
    end if

    ! Set defaults for MinFaceBoundary and MaxFaceBoundary 
    ! if they were not set by #FACEBOUNDARY command
    ! Using face BC is necessary if there is brick cut out of the
    ! spherical/cylindrical grid. 
    if(i_line_command("#FACEBOUNDARY", iSessionIn = 1) < 0)then
       ! Use all face based BCs by default for spherical geometry
       if(TypeGeometry(1:9) == 'spherical')then
          MinFaceBoundary = 1
          MaxFaceBoundary = 6
       end if

       ! Use face based boundaries by default for cylindrical geometry 
       ! except for top and bottom 
       if(TypeGeometry == 'cylindrical')then
          MinFaceBoundary = 1
          MaxFaceBoundary = 4
       end if
    end if

    ! Make sure MinFaceBoundary and MaxFaceBoundary cover face only boundaries
    if(UseBody2) MinFaceBoundary = min(Body2_, MinFaceBoundary)
    if(body1) then   
       MinFaceBoundary = min(Body1_, MinFaceBoundary)
       MaxFaceBoundary = max(Body1_, MaxFaceBoundary)
    end if
    if(UseExtraBoundary) then                        
       MinFaceBoundary = min(ExtraBc_, MinFaceBoundary)
       MaxFaceBoundary = max(ExtraBc_, MaxFaceBoundary)
    end if
    MaxFaceBoundary = min(MaxFaceBoundary, 6)
    MinFaceBoundary = max(MinFaceBoundary, body2_)

    if(index(TypeGeometry,'_genr') < 1) call set_gen_radial_grid

    ! Initialize BATL
    if(TypeGeometry(1:9)=='spherical') then
       TypeGeometryBatl = 'rlonlat'//TypeGeometry(10:20)
    else
       TypeGeometryBatl = TypeGeometry
    end if

    call init_batl(XyzMin_D(1:nDim), XyzMax_D(1:nDim), MaxBlock, &
         TypeGeometryBatl, TypeBc_I(1:2*nDim-1:2) == 'periodic', &
         nRootRead_D(1:nDim), UseRadiusIn=.true., UseDegreeIn=.false.,&
         RgenIn_I = exp(LogRGen_I), UseUniformAxisIn=UseUniformAxis,&
         UseFDFaceFluxIn=UseFDFaceFlux, iVectorVarIn_I=iVectorVar_I)

    if(IsRotatedCartesian)then
       ! Fix x1, x2 .. z2 to include the full rotated domain
       x2 = sum(abs(CoordMin_D)) + sum(abs(CoordMax_D)); y2 = x2; z2 = x2
       x1 = -x2; y1 = x1; z1 = x1
    end if

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

    use BATL_lib,    ONLY: radius_to_gen, Phi_, Theta_, nRoot_D, &
         CoordMin_D, CoordMax_D, nIJK_D
    use ModKind,     ONLY: nByteReal
    use ModIO

    implicit none

    integer :: iFile
    real    :: Ratio, PlotRes_D(3)

    real    :: CellSizeMax_D(3), Cut
    real    :: SmallSize_D(3)   ! Used for 2D plot areas

    character(len=*), parameter:: NameSub = 'correct_plot_range'
    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    if(DoTestMe)write(*,*) NameSub,' CoordMin_D, CoordMax_D=', &
         CoordMin_D, CoordMax_D

    ! Largest cell size and a much smaller distance for 2D cuts
    CellSizeMax_D = (CoordMax_D - CoordMin_D)/(nIJK_D*nRoot_D)

    if(nByteReal == 8)then
       SmallSize_D   = 1e-9*CellSizeMax_D
    else
       SmallSize_D   = 1e-6*CellSizeMax_D
    end if

    if(DoTestMe)write(*,*)NameSub,' CellSizeMax_D=',CellSizeMax_D

    PLOTFILELOOP: do iFile = Plot_+1, Plot_ + nPlotFile

       plot_area = plot_type(iFile)(1:3)

       if(DoTestMe)write(*,*)'iFile, plot_area=',iFile, plot_area

       ! Fix plot range for various plot areas
       select case(plot_area)
       case('shl', 'eqb', 'eqr', 'lcb')
          ! These plot areas read all ranges from PARAM.in
          CYCLE PLOTFILELOOP
       case('cut')
          if(IsLogRadius) plot_range(1:2,iFile) = log(plot_range(1:2,iFile))
          if(IsGenRadius) then
             call radius_to_gen(plot_range(1,iFile))
             call radius_to_gen(plot_range(2,iFile))
          end if
          if(Phi_ > 0) plot_range(2*Phi_-1:2*Phi_,iFile) = &
               cDegToRad*plot_range(2*Phi_-1:2*Phi_,iFile) 
          if(Theta_ > 0) plot_range(2*Theta_-1:2*Theta_,iFile) = &
               cDegToRad*plot_range(2*Theta_-1:2*Theta_,iFile)
          do iDim = 1, nDim
             if(plot_range(2*iDim-1,iFile) < plot_range(2*iDim,iFile)) CYCLE
             Cut = 0.5*(plot_range(2*iDim-1,iFile) + plot_range(2*iDim-1,iFile))
             plot_range(2*iDim-1,iFile) = Cut - SmallSize_D(iDim)
             plot_range(2*iDim,iFile)   = Cut + SmallSize_D(iDim)
          end do
          do iDim = nDim+1, MaxDim
             plot_range(2*iDim-1,iFile) = 0.0
             plot_range(2*iDim,iFile)   = 0.0
          end do
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
                plot_range(2*i-1,ifile) = CoordMin_D(i)
                plot_range(2*i,ifile)   = CoordMax_D(i)  
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
          plot_range(1:5:2, iFile) = CoordMin_D
          plot_range(2:6:2, iFile) = CoordMax_D
          if( plot_form(iFile)=='tec' .or. IsCartesianGrid )then
             ! Limit plot range along x direction to be very small
             plot_range(1, iFile) = -SmallSize_D(x_)
             plot_range(2, iFile) = +SmallSize_D(x_)
             if(IsRlonLat) then
                plot_range(1, iFile) = CoordMin_D(x_)
                plot_range(2, iFile) = CoordMin_D(x_)+SmallSize_D(x_)
             end if
          else
             ! Limit Phi direction around cHalfPi
             plot_range(3, iFile) = cHalfPi - SmallSize_D(Phi_)
             plot_range(4, iFile) = cHalfPi + SmallSize_D(Phi_)
          end if

       case('y=0')
          plot_range(1:5:2, iFile) = CoordMin_D
          plot_range(2:6:2, iFile) = CoordMax_D
          if(plot_form(iFile) == 'idl' .and. (IsRLonLat .or. IsCylindrical) &
               .and. CoordMin_D(2) < cPi .and. CoordMax_D(2) > cPi) then
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
          plot_range(1:5:2, iFile) = CoordMin_D
          plot_range(2:6:2, iFile) = CoordMax_D
          plot_range(5, iFile) = -SmallSize_D(z_)
          plot_range(6, iFile) = +SmallSize_D(z_)

       case('1d_')
          ! Limit plot range along 2nd and 3rd directions to be very small
          plot_range(1, iFile) = CoordMin_D(1)
          plot_range(2, iFile) = CoordMax_D(1)
          plot_range(3:5:2, iFile) = -SmallSize_D(y_:z_)
          plot_range(4:6:2, iFile) = +SmallSize_D(y_:z_)

       case('3d_')
          plot_range(1:5:2, iFile) = CoordMin_D
          plot_range(2:6:2, iFile) = CoordMax_D

       end select

       ! Reduce plot range in ignored dimensions
       if(nJ == 1) plot_range(3,iFile) = -SmallSize_D(y_)
       if(nJ == 1) plot_range(4,iFile) = +SmallSize_D(y_)
       if(nK == 1) plot_range(5,iFile) = -SmallSize_D(z_)
       if(nK == 1) plot_range(6,iFile) = +SmallSize_D(z_)

       if(DoTestMe)write(*,*)'For file ',ifile-plot_,&
            ' original range   =',plot_range(:,ifile)

       plot_range(1:5:2, iFile) = max(plot_range(1:5:2, iFile), CoordMin_D)
       plot_range(2:6:2, iFile) = min(plot_range(2:6:2, iFile), CoordMax_D)

       if(DoTestMe)write(*,*)'For file ',ifile-plot_,&
            ' limited range   =',plot_range(:,ifile)

       ! For plot_dx = 0.0 or -1.0 there is no need to adjust cut range
       if(plot_dx(1, iFile) <= cTiny)then

          plot_dx(:, iFile) = plot_dx(1, iFile) ! Define y and z

          CYCLE PLOTFILELOOP

       end if

       ! Make sure that plot resolution is a power of 2 fraction of cell size
       Ratio     = CellSizeMax_D(x_)/plot_dx(1, iFile)
       Ratio     = 2.0**nint(log(Ratio)/log(2.0))
       PlotRes_D = CellSizeMax_D / Ratio
       plot_dx(1:3, iFile) = PlotRes_D

       ! Make sure that plotting range is placed at an integer multiple of dx
       call adjust_plot_range(PlotRes_D(1), plot_range(:,iFile))
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
       DtFixed = DtFixedDim * Io2No_V(UnitT_)
       DtFixedOrig = DtFixed                   ! Store the initial setting
       Dt = DtFixed
       if(time_accurate) Cfl=1.0
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
    if(UseHighResChange) DoOneCoarserLayer = .false.

    DoLimitMomentum = boris_correction .and. DoOneCoarserLayer

!!! momentum limiting fails for multiion: to be debugged
    if(UseMultiIon)DoLimitMomentum = .false.
!!!

  end subroutine set_extra_parameters
  !============================================================================
  subroutine sort_smooth_indicator

    ! The variables using the same smooth indicator should be 
    ! calculated one by one. The smooth indicator 
    ! itself is calculated first. 
    ! iVarSmoothIndex_I is the calculation order.

    real:: iVarSmoothReal_V(nVar)
    integer:: iVar

    do iVar = 1, nVar
       if(iVarSmooth_V(iVar) == iVar) then
          iVarSmoothReal_V(iVar) = real(iVar) - 0.5
       else
          iVarSmoothReal_V(iVar) = real(iVarSmooth_V(iVar))
       endif
    enddo

    call sort_quick(nVar,iVarSmoothReal_V,iVarSmoothIndex_I)

  end subroutine sort_smooth_indicator
  !============================================================================

end subroutine MH_set_parameters
!=======================================================================
subroutine set_levels
  use ModAMR, ONLY: min_block_level, max_block_level
  use BATL_lib, ONLY: iTree_IA, MinLevel_, MaxLevel_

  iTree_IA(MinLevel_,:) = min_block_level
  iTree_IA(MaxLevel_,:) = max_block_level

end subroutine set_levels
