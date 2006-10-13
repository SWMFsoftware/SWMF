!^CFG COPYRIGHT UM
subroutine MH_set_parameters(TypeAction)

  ! Set input parameters for the Global Magnetosphere (GM) module

  use ModProcMH
  use ModMain
  use ModAdvance
  use ModGeometry, ONLY : init_mod_geometry, &
       TypeGeometry,UseCovariant,UseVertexBasedGrid,is_axial_geometry,  & !^CFG IF COVARIANT
       allocate_face_area_vectors,allocate_old_levels,rTorusLarge,rTorusSmall,& !^CFG IF COVARIANT
       x1,x2,y1,y2,z1,z2,XyzMin_D,XyzMax_D,MinBoundary,MaxBoundary
  use ModNodes, ONLY : init_mod_nodes
  use ModImplicit                                       !^CFG IF IMPLICIT
  use ModPhysics
  use ModProject                                        !^CFG IF PROJECTION
  use ModCT, ONLY : init_mod_ct, DoInitConstrainB       !^CFG IF CONSTRAINB
  use ModBlockData, ONLY: clean_block_data
  use ModAMR
  use ModParallel, ONLY : proc_dims
  use ModRaytrace                                       !^CFG IF RAYTRACE
  use ModIO
  use CON_planet,       ONLY: read_planet_var, check_planet_var
  use CON_axes,         ONLY: init_axes
  use ModUtilities,     ONLY: check_dir, fix_dir_name, upper_case, lower_case,&
       DoFlush

  use CON_planet,       ONLY: get_planet
  use ModTimeConvert,   ONLY: time_int_to_real, time_real_to_int
  use ModCoordTransform,ONLY: rot_matrix_x, rot_matrix_y, rot_matrix_z
  use ModReadParam
  use ModMPCells,       ONLY: iCFExchangeType,DoOneCoarserLayer
  use ModFaceValue,     ONLY: UseTvdResChange, UseAccurateResChange, &
       DoLimitMomentum, &                              !^CFG IF BORISCORR
       BetaLimiter, TypeLimiter
  use ModPartSteady,    ONLY: UsePartSteady, MinCheckVar, MaxCheckVar, &
       RelativeEps_V, AbsoluteEps_V
  use ModUser,          ONLY: user_read_inputs, user_init_session, &
       NameUserModule, VersionUserModule
  use ModBoundaryCells, ONLY: SaveBoundaryCells,allocate_boundary_cells
  use ModPointImplicit, ONLY: UsePointImplicit, BetaPointImpl
  use ModrestartFile,   ONLY: read_restart_parameters
  use ModHallResist,    ONLY: UseHallResist, HallFactor, HallCmaxFactor,&
       NameHallRegion, R1Hall, R2Hall, HallWidth

  implicit none

  character (len=17) :: NameSub='MH_set_parameters'

  ! Arguments

  ! TypeAction determines if we read or check parameters
  character (len=*), intent(in) :: TypeAction

  ! Local variables
  integer :: ifile, i,j, iError
  real :: local_root_dx

  logical :: IsUninitialized=.true.
  logical :: read_new_upstream=.false.
  logical :: DoReadSatelliteFiles=.false.

  ! The name of the command
  character (len=lStringLine) :: NameCommand, StringLine

  ! Variables to remember for multiple calls
  character (len=lStringLine) :: UpstreamFileName='???'

  ! Temporary variables
  logical :: DoEcho
  integer :: nVarRead
  character (len=lStringLine) :: NameEquationRead

  character (len=50) :: plot_string,log_string,satellite_string
  character (len=3)  :: plot_area, plot_var, satellite_var, satellite_form
  character (len=2)  :: NameCompRead
  integer :: problem_type_r, qtotal, nIJKRead_D(3)

  integer            :: TimingDepth=-1
  character (len=10) :: TimingStyle='cumm'

  ! Variables for checking/reading #STARTTIME command
  real (Real8_)         :: StartTimeCheck = -1.0_Real8_

  ! Variables for the #GRIDRESOLUTION and #GRIDLEVEL commands
  character (len=lStringLine) :: NameArea
  integer                     :: nLevelArea
  real                        :: AreaResolution, RadiusArea
  logical                     :: DoReadAreaCenter
  real, dimension(3)          :: XyzStartArea_D, XyzEndArea_D
  real                        :: xRotateArea, yRotateArea, zRotateArea

  ! Variables for checking the user module
  character (len=lStringLine) :: NameUserModuleRead
  real                        :: VersionUserModuleRead

  integer :: iSession, iPlotFile, iVar
  !-------------------------------------------------------------------------
  NameSub(1:2) = NameThisComp

  iSession = i_session_read()

  if(IsUninitialized)then
     call set_defaults
     IsUninitialized=.false.
  end if

  if(iSession>1)then
     restart=.false.           ! restart in session 1 only
     read_new_upstream=.false. ! upstream file reading in session 1 only
  end if

  if(DoReadSatelliteFiles)then
     call read_satellite_input_files
     DoReadSatelliteFiles = .false.
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

     ! In standalone mode set and obtain GM specific parameters 
     ! in CON_planet and CON_axes
     if(IsStandAlone .and. NameThisComp=='GM') then
        ! Check and set some planet variables (e.g. DoUpdateB0)
        call check_planet_var(iProc==0, time_accurate)

        ! Initialize axes
        call init_axes(StartTime)

        if(body1)then
           call get_planet(UseRotationOut = UseRotatingBc)
        else
           UseRotatingBc = .false.
        end if

        ! Obtain some planet parameters
        if(Bdp_dim == 0.0)then
           DoUpdateB0 = .false.
           Dt_UpdateB0 = -1.0
        else
           call get_planet( &
                DoUpdateB0Out = DoUpdateB0, DtUpdateB0Out = Dt_UpdateB0)
        end if

     end if

     if(problem_type==-1) call stop_mpi(&
          NameSub//': ERROR problem type must be set in PARAM.in')

     call correct_parameters

     call check_parameters

     call check_plot_range

     call check_options

     ! initialize module variables
     call init_mod_advance
     call init_mod_geometry
     call init_mod_nodes
     call init_mod_raytrace                    !^CFG IF RAYTRACE
     if(UseConstrainB) call init_mod_ct        !^CFG IF CONSTRAINB
     if(UseImplicit)   call init_mod_implicit  !^CFG IF IMPLICIT

     ! clean dynamic storage
     call clean_block_data

     if (read_new_upstream) &
          call read_upstream_input_file(UpstreamFileName)

     call set_physics_constants

     call set_extra_parameters

     ! initialize ModEqution (e.g. variable units)
     call init_mod_equation

     ! Initialize user module and allow user to modify things
     if(UseUserInitSession)call user_init_session

     if(iProc==0 .and. IsStandAlone)then
        call timing_active(UseTiming)
        if(iSession==1)call timing_step(0)
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
     case("#COMPONENT")
        call read_var('NameComp',NameCompRead)
        if(NameCompRead /= NameThisComp)&
             call stop_mpi(NameSub//' ERROR: BATSRUS is running as component '&
             //NameThisComp//' and not as '//NameCompRead)
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
        end if
     case("#OUTERBOUNDARY")
        call read_var('TypeBcEast'  ,TypeBc_I(east_))  
        call read_var('TypeBcWest'  ,TypeBc_I(west_))
        call read_var('TypeBcSouth' ,TypeBc_I(south_))
        call read_var('TypeBcNorth' ,TypeBc_I(north_))
        call read_var('TypeBcBottom',TypeBc_I(bot_))
        call read_var('TypeBcTop'   ,TypeBc_I(top_))    
     case("#INNERBOUNDARY")
        call read_var('TypeBcInner',TypeBc_I(body1_))
        !                                              ^CFG IF SECONDBODY BEGIN
        if(UseBody2) &                                      
             call read_var('TypeBcBody2',TypeBc_I(body2_)) 
        !                                              ^CFG END SECONDBODY
     case("#TIMESTEPPING")
        call read_var('nStage',nSTAGE)
        call read_var('CflExpl',Cfl)
        ExplCfl = Cfl                                   !^CFG IF IMPLICIT
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
        call read_var('UsePointImplicit', UsePointImplicit)
        if(UsePointImplicit)then
           call read_var('BetaPointImplicit',BetaPointImpl)
           UsePartImplicit = .false.                   !^CFG IF IMPLICIT
           UseFullImplicit = .false.                   !^CFG IF IMPLICIT
        end if
     case("#PARTLOCAL")                                !^CFG IF IMPLICIT BEGIN 
        call read_var('UsePartLocal',UsePartLocal)

     case("#IMPLICIT")                                 
        call read_var('UsePointImplicit', UsePointImplicit)
        call read_var('UsePartImplicit', UsePartImplicit)
        call read_var('UseFullImplicit', UseFullImplicit)

        UseImplicit = UseFullImplicit .or. UsePartImplicit

        ! For implicit scheme it is better to use unsplit dB0/dt evaluation
        DoSplitDb0Dt = .not.UseImplicit

        if(UsePartImplicit  .and. UseFullImplicit) call stop_mpi(&
             'Only one of UsePartImplicit and UseFullImplicit can be true')

        if(UsePointImplicit .and. UseImplicit) call stop_mpi(&
             'Only one of Point-/Full-/Partimmplicit can be true')

        if(UseImplicit)call read_var('ImplCFL',ImplCFL)

     case("#IMPLCRITERIA", "#IMPLICITCRITERIA", "#STEPPINGCRITERIA")
        call read_var('TypeImplCrit',ImplCritType)
        select case(ImplCritType)
        case('R','r')
           call read_var('rImplicit'   ,Rimplicit)
        case('test','dt')
        case default
           if(iProc==0)then
              write(*,'(a)')NameSub// &
                   ' WARNING: invalid ImplCritType='//trim(ImplCritType)// &
                   ' !!!'
              if(UseStrict)call stop_mpi('Correct PARAM.in!')
              write(*,*)NameSub//' Setting ImplCritType=R'
           end if
           ImplCritType='R'
        end select
     case("#PARTIMPL", "#PARTIMPLICIT")
        call read_var('UsePartImplicit2',UsePartImplicit2)
     case("#IMPLSCHEME", "#IMPLICITSCHEME")
        call read_var('nOrderImpl',nORDER_impl)
        call read_var('TypeFluxImpl',FluxTypeImpl)
     case("#IMPLSTEP", "#IMPLICITSTEP")
        call read_var('ImplCoeff ',ImplCoeff0)
        call read_var('UseBDF2   ',UseBDF2)
        call read_var('ImplSource',implsource)
     case("#IMPLCHECK", "#IMPLICITCHECK")
        call read_var('RejectStepLevel' ,   RejectStepLevel)
        call read_var('RejectStepFactor',   RejectStepFactor)
        call read_var('ReduceStepLevel' ,   ReduceStepLevel)
        call read_var('ReduceStepFactor',   ReduceStepFactor)
        call read_var('IncreaseStepLevel' , IncreaseStepLevel)
        call read_var('IncreaseStepFactor', IncreaseStepFactor)
     case("#NEWTON")
        call read_var('UseConservativeImplicit',UseConservativeImplicit)
        call read_var('UseNewton',UseNewton)
        if(UseNewton)then
           call read_var('UseNewMatrix ',NewMatrix)
           call read_var('MaxIterNewton',NewtonIterMax)
        endif
     case("#JACOBIAN")
        call read_var('TypeJacobian',JacobianType)
        call read_var('JacobianEps', JacobianEps)
     case("#PRECONDITIONER")
        call read_var('TypePrecondSide',PrecondSide)
        call read_var('TypePrecond'    ,PrecondType)
        call read_var('GustafssonPar'  ,GustafssonPar)
        if(GustafssonPar<=0.)GustafssonPar=0.
        if(GustafssonPar>1.)GustafssonPar=1.
     case("#KRYLOV")
        call read_var('TypeKrylov'     ,KrylovType)
        call read_var('TypeInitKrylov' ,KrylovInitType)
        call read_var('ErrorMaxKrylov' ,KrylovErrorMax)
        call read_var('MaxMatvecKrylov',KrylovMatvecMax)
        nKrylovVector = KrylovMatvecMax
     case("#KRYLOVSIZE")
        call read_var('nKrylovVector',nKrylovVector)    !^CFG END IMPLICIT

     case("#HEATFLUX")                                  !^CFG IF DISSFLUX BEGIN
        call read_var('UseHeatFlux'   ,UseHeatFlux)
        call read_var('UseSpitzerForm',UseSpitzerForm)
        if (.not.UseSpitzerForm) then
           call read_var('Kappa0Heat'  ,Kappa0Heat)
           call read_var('ExponentHeat',ExponentHeat)
        end if
     case("#RESISTIVEFLUX")
        call read_var('UseResistFlux' ,UseResistFlux)
        if(UseResistFlux)then
           call read_var('UseSpitzerForm',UseSpitzerForm)
           if (.not.UseSpitzerForm) then
              call read_var('TypeResist',TypeResist)
              call read_var('Eta0Resist',Eta0Resist)
              if (TypeResist=='Localized'.or. &
                   TypeResist=='localized') then
                 call read_var('Alpha0Resist',Alpha0Resist)
                 call read_var('yShiftResist',yShiftResist)
                 call read_var('TimeInitRise',TimeInitRise)
                 call read_var('TimeConstLev',TimeConstLev)
              end if
           end if
           call read_var('UseAnomResist',UseAnomResist)
           if (UseAnomResist) then
              call read_var('Eta0AnomResist',   Eta0AnomResist)
              call read_var('EtaAnomMaxResist', EtaAnomMaxResist)
              call read_var('jCritResist',      jCritResist)
           end if
        end if
        !                                               ^CFG END DISSFLUX
     case("#HALLRESISTIVITY")
        call read_var('UseHallResist',  UseHallResist)
        if(UseHallResist)then
           call read_var('HallFactor',  HallFactor)
           call read_var('HallCmaxFactor', HallCmaxFactor)
        end if
     case("#HALLREGION")
        call read_var('NameHallRegion', NameHallRegion)

        select case(NameHallRegion)
        case("all", "user")
        case("shell")
           call read_var("R1Hall",R1Hall)
           call read_var("R2Hall",R2Hall)
           call read_var("HallWidth",HallWidth)
        case default
           call stop_mpi(NameSub//': unknown NameHallRegion='&
                //NameHallRegion)
        end select

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
     case("#SAVEPLOT")
        call read_var('nPlotFile',nplotfile)

        if(nPlotFile>0 .and. iProc==0)call check_dir(NamePlotDir)

        nfile=max(nfile,plot_+nplotfile)
        if (nfile > maxfile .or. nplotfile > maxplotfile)call stop_mpi(&
             'The number of ouput files is too large in #SAVEPLOT:'&
             //' nplotfile>maxplotfile .or. nfile>maxfile')
        do iFile=plot_+1,plot_+nplotfile

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
              call read_var('NameLine',NameLine_I(iPlotFile))
              call upper_case(NameLine_I(iPlotFile))
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
              do i=1,nLine_I(iPlotFile)
                 call read_var('xStartLine',XyzStartLine_DII(1,i,iPlotFile))
                 call read_var('yStartLine',XyzStartLine_DII(2,i,iPlotFile))
                 call read_var('zStartLine',XyzStartLine_DII(3,i,iPlotFile))
                 call read_var('IsParallel',IsParallelLine_II(i,iPlotFile))
              end do                                  !^CFG END RAYTRACE
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
           elseif (index(plot_string,'ion')>0) then
              plot_area='ion'
           else
              do i=1,3
                 plot_range(2*i-1,ifile)=XyzMin_D(i)
                 plot_range(2*i  ,ifile)=XyzMax_D(i)
              end do
              if(index(plot_string,'x=0')>0)then
                 plot_area='x=0'
                 if(.false.&
                      .or.(.not.is_axial_geometry()) & !^CFG IF COVARIANT
                      .or.index(plot_string,'tec')>0)then     
                    plot_range(1,ifile)=-cTiny*(XyzMax_D(x_)-XyzMin_D(x_)) &
                         /(nCells(1)*proc_dims(1))
                    plot_range(2,ifile)=+cTiny*(XyzMax_D(x_)-XyzMin_D(x_)) &
                         /(nCells(1)*proc_dims(1))
                 else                              !^CFG IF COVARIANT BEGIN
                    plot_range(3,ifile)=cHalfPi                 &
                         -cTiny*(XyzMax_D(Phi_)-XyzMin_D(Phi_)) &
                         /(nCells(Phi_)*proc_dims(Phi_))            
                    plot_range(4,ifile)=cHalfPi                 &
                         +cTiny*(XyzMax_D(Phi_)-XyzMin_D(Phi_)) &
                         /(nCells(Phi_)*proc_dims(Phi_))  !^CFG END COVARIANT 
                 end if                            

              elseif(index(plot_string,'y=0')>0)then
                 plot_area='y=0'
                 plot_range(3,ifile)=-cTiny*(XyzMax_D(2)-XyzMin_D(2))&
                      /(nCells(2)*proc_dims(2)) 
                 plot_range(4,ifile)=+cTiny*(XyzMax_D(2)-XyzMin_D(2))&
                      /(nCells(2)*proc_dims(2)) 
              elseif(index(plot_string,'z=0')>0)then
                 plot_area='z=0'
                 plot_range(5,ifile)=-cTiny*(XyzMax_D(3)-XyzMin_D(3))&
                      /(nCells(3)*proc_dims(3)) 
                 plot_range(6,ifile)=+cTiny*(XyzMax_D(3)-XyzMin_D(3))&
                      /(nCells(3)*proc_dims(3)) 
              elseif(index(plot_string,'3d')>0)then
                 plot_area='3d_'
              else
                 call stop_mpi('Area (x=0,y=0,z=0,3d,cut,sph,ion) missing'&
                      //' from plot_string='//plot_string)
              endif
           end if

           ! Plot file format
           if(index(plot_string,'idl')>0)then
              plot_form(ifile)='idl'
              if ((plot_area /= 'ion')&
                   .and. plot_area /= 'sph' &
                   .and. plot_area /= 'los' &
                   .and. plot_area /= 'lin' &        !^CFG IF RAYTRACE
                   ) call read_var('DxSavePlot',plot_dx(1,ifile))
              if(is_axial_geometry())plot_dx(1,ifile)=-1.0 !^CFG IF COVARIANT 
   
           elseif(index(plot_string,'tec')>0)then 
              plot_form(ifile)='tec'
              plot_dx(1,ifile)=0.
           else
              call stop_mpi('Format (idl,tec) missing from plot_string='&
                   //plot_string)
           end if
           if (plot_area == 'sph') then
              select case(TypeGeometry)                !^CFG IF COVARIANT
              case('cartesian')                        !^CFG IF COVARIANT
                 plot_dx(1,ifile) = 1.0    ! set to match value in write_plot_sph
                 plot_dx(2:3,ifile) = 1.0  ! set to degrees desired in angular resolution
                 plot_range(2,ifile)= plot_range(1,ifile) + 1.e-4   ! so that R/=0
                 plot_range(3,ifile)= 0.   - 0.5*plot_dx(2,ifile)
                 plot_range(4,ifile)= 90.0 + 0.5*plot_dx(2,ifile)
                 plot_range(5,ifile)= 0.   - 0.5*plot_dx(3,ifile)
                 plot_range(6,ifile)= 360.0- 0.5*plot_dx(3,ifile)
              case('spherical')                    !^CFG IF COVARIANT BEGIN
                 plot_dx(1,ifile) = -1.0   
                 plot_range(2,ifile)= plot_range(1,ifile) + 1.e-4   ! so that R/=0 
                 do i=Phi_,Theta_
                    plot_range(2*i-1,ifile)=XyzMin_D(i)
                    plot_range(2*i,ifile)=XyzMax_D(i)  
                 end do
                 plot_area='R=0'           ! to disable the write_plot_sph routine
              case('spherical_lnr')           
                 plot_dx(1,ifile) = -1.0  
                 plot_range(1,ifile)=alog(max(plot_range(1,ifile),cTiny)) 
                 plot_range(2,ifile)= plot_range(1,ifile) + 1.e-4   ! so that R/=0 
                 do i=Phi_,Theta_
                    plot_range(2*i-1,ifile)=XyzMin_D(i)
                    plot_range(2*i,ifile)=XyzMax_D(i)  
                 end do
                 plot_area='R=0'           ! to disable the write_plot_sph routine
              case default
                 call stop_mpi(NameSub//' Sph-plot is not implemented for geometry= '&
                      //TypeGeometry)
              end select                           !^CFG END COVARIANT 
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
              plot_pars(ifile)='R_ray'
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
           elseif(index(plot_string,'pos')>0.or.index(plot_string,'POS')>0)then
              plot_var='pos'
              plot_dimensional(ifile) = index(plot_string,'POS')>0
              if(plot_area /= 'lin')call stop_mpi(&
                   'Variable "pos" can only be used with area "lin" !')
           else
              call stop_mpi('Variable definition missing from plot_string=' &
                   //plot_string)
           end if

           plot_type(ifile)=plot_area//'_'//plot_var
        end do
     case("#SAVEPLOTSAMR")
        call read_var('DoSavePlotsAmr',save_plots_amr)
     case("#SAVEBINARY")
        call read_var('DoSaveBinary',save_binary)

     case("#GRIDRESOLUTION","#GRIDLEVEL","#AREARESOLUTION","#AREALEVEL")
        if(index(NameCommand,"RESOLUTION")>0)then
           call read_var('AreaResolution',AreaResolution)
        else
           call read_var('nLevelArea',nLevelArea)
        end if
        call read_var('NameArea',NameArea)
        call lower_case(NameArea)

        NameArea = adjustl(NameArea)

        if(NameArea(1:4) == 'init')then
           ! 'init' or 'initial' means that the initial resolution is set,
           ! and no area is created. Replaces the #AMRINIT 'none' nlevel.
           ! Convert resolution to level if necessary
           if(index(NameCommand,"RESOLUTION")>0) &
                nLevelArea = nint( &
                alog(((XyzMax_D(x_)-XyzMin_D(x_)) / (proc_dims(x_) * nI))  &
                / AreaResolution) / alog(2.0) )

           ! Set the initial levels and the refinement type to 'none'
           initial_refine_levels = nLevelArea
           InitialRefineType     = 'none'

           ! No area is created, continue reading the parameters
           CYCLE READPARAM
        end if

        ! Convert level to resolution if necessary
        if(index(NameCommand,"LEVEL")>0) &
             AreaResolution = (XyzMax_D(x_)-XyzMin_D(x_)) &
             / (proc_dims(x_) * nI * 2.0**nLevelArea)

        nArea = nArea + 1
        if(nArea > MaxArea)then
           if(UseStrict)call stop_mpi(NameSub// &
                ' ERROR: Too many grid areas were defined')
           if(iProc == 0)then
              write(*,*)NameSub," nArea = ",nArea
              write(*,*)NameSub," WARNING: Too many grid areas were defined"
              write(*,*)NameSub," ignoring command ",NameCommand
           end if
           nArea = MaxArea
        end if

        Area_I(nArea) % Resolution = AreaResolution

        ! Set the default center to be the origin, 
        ! the size to be 1 and no rotatoin
        Area_I(nArea)%Center_D  = 0.0
        Area_I(nArea)%Size_D    = 1.

        ! Check for the word rotated in the name
        i = index(NameArea,'rotated')
        Area_I(nArea)%DoRotate = i > 0
        if(i>0) NameArea = NameArea(1:i-1)//NameArea(i+7:len(NameArea))

        ! Check for the character '0' in the name
        i = index(NameArea,'0')
        DoReadAreaCenter = i < 1
        if(i>0) NameArea = NameArea(1:i-1)//NameArea(i+1:len(NameArea))

        ! Remove leading spaces
        NameArea = adjustl(NameArea)
        Area_I(nArea)%Name = NameArea

        ! Read center of area if needed
        if(DoReadAreaCenter .and. &
             NameArea /= 'all' .and. NameArea /= 'box') then
           call read_var("xCenter",Area_I(nArea)%Center_D(1))
           call read_var("yCenter",Area_I(nArea)%Center_D(2))
           call read_var("zCenter",Area_I(nArea)%Center_D(3))
        endif

        select case(NameArea)
        case("all")
           ! No geometry info is needed for uniform refinement

        case("box")
           call read_var("xMinBox",XyzStartArea_D(1))
           call read_var("yMinBox",XyzStartArea_D(2))
           call read_var("zMinBox",XyzStartArea_D(3))
           call read_var("xMaxBox",XyzEndArea_D(1))
           call read_var("yMaxBox",XyzEndArea_D(2))
           call read_var("zMaxBox",XyzEndArea_D(3))
           ! Convert to center and size information
           Area_I(nArea)%Center_D = 0.5*   (XyzStartArea_D + XyzEndArea_D)
           Area_I(nArea)%Size_D   = 0.5*abs(XyzEndArea_D - XyzStartArea_D)

           ! Overwrite name with brick
           Area_I(nArea)%Name     = "brick"

        case("brick")
           call read_var("xSize", Area_I(nArea)%Size_D(1))
           call read_var("ySize", Area_I(nArea)%Size_D(2))
           call read_var("zSize", Area_I(nArea)%Size_D(3))

           ! Area size is measured from the center: half of brick size
           Area_I(nArea)%Size_D   = 0.5*Area_I(nArea)%Size_D

        case("shell")
           call read_area_radii
           Area_I(nArea)%Size_D = RadiusArea

        case("sphere")
           call read_var("Radius", RadiusArea)
           Area_I(nArea)%Size_D   = RadiusArea

        case("ringx")
           call read_var("Height", Area_I(nArea)%Size_D(1))
           Area_I(nArea)%Size_D(1)     = 0.5*Area_I(nArea)%Size_D(1)
           call read_area_radii
           Area_I(nArea)%Size_D(2:3)   = RadiusArea

        case("ringy")
           call read_var("Height", Area_I(nArea)%Size_D(2))
           Area_I(nArea)%Size_D(2)     = 0.5*Area_I(nArea)%Size_D(2)
           call read_area_radii
           Area_I(nArea)%Size_D(1:3:2) = RadiusArea

        case("ringz")
           call read_var("Height", Area_I(nArea)%Size_D(3))
           Area_I(nArea)%Size_D(3)     = 0.5*Area_I(nArea)%Size_D(3)
           call read_area_radii
           Area_I(nArea)%Size_D(1:2)   = RadiusArea

        case("cylinderx")
           call read_var("Length", Area_I(nArea)%Size_D(1))
           Area_I(nArea)%Size_D(1)     = 0.5 * Area_I(nArea)%Size_D(1)
           call read_var("Radius", RadiusArea)
           Area_I(nArea)%Size_D(2:3)   = RadiusArea

        case("cylindery")
           call read_var("Length", Area_I(nArea)%Size_D(2))
           Area_I(nArea)%Size_D(2)     = 0.5 * Area_I(nArea)%Size_D(2)
           call read_var("Radius", RadiusArea)
           Area_I(nArea)%Size_D(1:3:2) = RadiusArea

        case("cylinderz")
           call read_var("Length", Area_I(nArea)%Size_D(3))
           Area_I(nArea)%Size_D(3)     = 0.5 * Area_I(nArea)%Size_D(3)
           call read_var("Radius", RadiusArea)
           Area_I(nArea)%Size_D(1:2)   = RadiusArea

        case default
           if(UseStrict) call stop_mpi(NameSub//&
                   ' ERROR: unknown NameArea='//trim(Area_I(nArea)%Name))

           if(iProc == 0) &
                write(*,*) NameSub//' WARNING: unknown NameArea=' // &
                trim(Area_I(nArea)%Name) // ', ignoring command ' // &
                trim(NameCommand)

           nArea = nArea - 1
           CYCLE READPARAM
        end select

        if(Area_I(nArea) % DoRotate)then

           ! Read 3 angles for the rotation matrix in degrees
           call read_var('xRotate',xRotateArea)
           call read_var('yRotate',yRotateArea)
           call read_var('zRotate',zRotateArea)

           ! Rotation matrix rotates around X, Y and Z axes in this order
           Area_I(nArea)%Rotate_DD = matmul( matmul( &
                rot_matrix_z(-zRotateArea*cDegToRad), &
                rot_matrix_y(-yRotateArea*cDegToRad)),&
                rot_matrix_x(-xRotateArea*cDegToRad))

           !write(*,*)'Matrix=',Area_I(nArea)%Rotate_DD

        end if
     case("#AMRLEVELS")
        call read_var('MinBlockLevel',min_block_level)
        call read_var('MaxBlockLevel',max_block_level)
        call read_var('DoFixBodyLevel',fix_body_level)
        if(iSession==1)then
           DoSetLevels=.true.
        else
           call set_levels
        end if
     case("#AMRRESOLUTION")
        call read_var('DxCellMin',min_cell_dx)
        call read_var('DxCellMax',max_cell_dx)
        call read_var('DoFixBodyLevel',fix_body_level)
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
        if(iSession==1)then
           DoSetLevels=.true.
        else
           call set_levels
        end if
     case("#AMR")
        call read_var('DnRefine',dn_refine)
        if (dn_refine > 0)then
           call read_var('DoAutoRefine',automatic_refinement)
           if (automatic_refinement) then
              call read_var('PercentCoarsen',percentCoarsen)
              call read_var('PercentRefine' ,percentRefine)
              call read_var('MaxTotalBlocks',MaxTotalBlocks)
           end if
        end if
     case("#AMRINIT")
        call read_var('TypeRefineInit'  ,InitialRefineType)
        call read_var('nRefineLevelInit',initial_refine_levels)
     case("#AMRCRITERIA")
        call read_var('nRefineCrit',nRefineCrit)
        if(nRefineCrit<0 .or. nRefineCrit>3)&
             call stop_mpi(NameSub//' ERROR: nRefineCrit must be 0, 1, 2 or 3')
        do i=1,nRefineCrit
           call read_var('TypeRefine',RefineCrit(i))
           if(RefineCrit(i)=='Transient'.or.RefineCrit(i)=='transient') then
              call read_var('TypeTransient_I(i)',TypeTransient_I(i))
              call read_var('UseSunEarth'       ,UseSunEarth)
           end if
           if (UseSunEarth) then
              call read_var('xEarth'  ,xEarth)
              call read_var('yEarth'  ,yEarth)
              call read_var('zEarth'  ,zEarth)
              call read_var('InvD2Ray',InvD2Ray)
           end if
        end do
     case("#SCHEME")
        call read_var('nOrder'  ,nOrder)
        call read_var('TypeFlux',FluxType)
        if(nOrder>1)&                                                
             call read_var('TypeLimiter', TypeLimiter)
        if(TypeLimiter == 'minmod') then
           BetaLimiter = 1.0
        else
           call read_var('LimiterBeta', BetaLimiter)
        end if
     case("#NONCONSERVATIVE")
        call read_var('UseNonConservative',UseNonConservative)
     case("#CONSERVATIVECRITERIA")
        call read_var('nConservCrit',nConservCrit)
        if(nConservCrit > 0) then
           if(allocated(TypeConservCrit_I)) deallocate(TypeConservCrit_I)
           allocate( TypeConservCrit_I(nConservCrit) )
           do i=1,nConservCrit
              call read_var('TypeConservCrit',TypeConservCrit_I(i) )
              call lower_case(TypeConservCrit_I(i))
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
              case('gradp')
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
        if(.not.is_axial_geometry()) then               !^CFG IF COVARIANT
           call read_var('TypeMessagePass',optimize_message_pass)
           if(optimize_message_pass=='allold' .or.&
                optimize_message_pass=='oldopt')then
              if(iProc==0)write(*,'(a)')NameSub// &
                   ' WARNING: message_pass mode='// &
                   trim(optimize_message_pass)// &
                   ' is not available any longer, allopt is set !!!'
              optimize_message_pass='allopt'
           end if
        end if                                           !^CFG IF COVARIANT
     case('#RESCHANGE','#RESOLUTIONCHANGE')
        call read_var('UseAccurateResChange',UseAccurateResChange)
        if(UseAccurateResChange) UseTvdResChange=.false.
     case('#TVDRESCHANGE')
        call read_var('UseTvdResChange',UseTvdResChange)
        if(UseTvdResChange) UseAccurateResChange=.false.
     case("#BORIS")
        !                                              ^CFG IF BORISCORR BEGIN
        call read_var('UseBorisCorrection',boris_correction)   
        if(boris_correction) then
           call read_var('BorisClightFactor',boris_cLIGHT_factor)
           UseBorisSimple=.false.     !^CFG IF SIMPLEBORIS
        else
           boris_cLIGHT_factor = 1.0
        end if
        !                                              ^CFG END BORISCORR
     case("#BORISSIMPLE")                            !^CFG IF SIMPLEBORIS BEGIN
        call read_var('UseBorisSimple',UseBorisSimple)
        if(UseBorisSimple) then
           call read_var('BorisClightFactor',boris_cLIGHT_factor)
           boris_correction=.false.                     !^CFG IF BORISCORR
        else
           boris_cLIGHT_factor = 1.0
        end if                                       !^CFG END SIMPLEBORIS
     case("#DIVB")
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
     case("#DIVBSOURCE")
	call read_var('UseB0Source'   ,UseB0Source)
     case("#PROJECTION")                              !^CFG IF PROJECTION BEGIN
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
     case("#SHOCKTUBE")
        do i=1,nVar
           call read_var('LeftState',shock_Lstate(i))
        end do
        do i=1,nVar
           call read_var('RightState',shock_Rstate(i))
        end do
        call read_var('ShockSlope',ShockSlope)
     case("#SOLARWINDFILE", "#UPSTREAM_INPUT_FILE")
        call read_var('UseSolarWindFile',UseUpstreamInputFile)
        if (UseUpstreamInputFile) then
           read_new_upstream = .true.
           call read_var('NameSolarWindFile', UpstreamFileName)
        end if
        !                                               ^CFG IF RAYTRACE BEGIN
     case("#RAYTRACE")
        call read_var('UseAccurateIntegral',UseAccurateIntegral)
        call read_var('UseAccurateTrace'   ,UseAccurateTrace)
        if(UseAccurateTrace .and. .not. UseAccurateIntegral)then
           if(iProc==0)then
              write(*,'(a)')NameSub//' WARNING: '// &
                   'UseAccurateTrace=T requires UseAccurateIntegral=T'
              if(UseStrict)call stop_mpi('Correct PARAM.in!')
              write(*,*)NameSub//' setting UseAccurateIntegral=T'
           end if
           UseAccurateIntegral = .true.
        end if
        call read_var('DtExchangeRay',DtExchangeRay)
        call read_var('DnRaytrace',   DnRaytrace)
        !                                              ^CFG END RAYTRACE
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
        if(NameCommand == "#IMCOUPLING")then
           call read_var('DoCoupleImPressure',DoCoupleImPressure)
           call read_var('DoCoupleImDensity',DoCoupleImDensity)
        end if                                        !^CFG END RCM
     case("#MASSLOADING")
        call read_var('UseMassLoading',UseMassLoading)
        call read_var('DoAccelerateMassLoading',AccelerateMassLoading) 
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
     case("#PROBLEMTYPE")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('iProblemType',problem_type_r) 
        if(problem_type<0)then
           problem_type=problem_type_r
           if (problem_type==problem_dissipation) &
                call read_var('TypeProblemDiss',TypeProblemDiss)
           call set_problem_defaults
        else if(problem_type/=problem_type_r) then
           if(iProc==0)write(*,*)'problem type already set to',&
                problem_type,' /=',problem_type_r
           call stop_mpi('PROBLEM TYPE CANNOT BE OVERWRITTEN !')
        else
           if(iProc==0)write(*,'(a)')NameSub // &
                ' WARNING: problem type set twice !!!'
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
        call read_var('nSatellite',nsatellite)
        if(nSatellite <= 0) CYCLE READPARAM
        save_satellite_data=.true.
        DoReadSatelliteFiles=.true.
        if(iProc==0) call check_dir(NamePlotDir)
        nFile = max(nFile, satellite_ + nSatellite)
        if (nFile > MaxFile .or. nSatellite > MaxSatelliteFile)&
             call stop_mpi(&
             'The number of output files is too large in #SATELLITE:'&
             //' nFile > MaxFile .or. nSatellite > MaxSatelliteFile')

        do iFile = satellite_+1, satellite_ + nSatellite
           call read_var('StringSatellite',satellite_string)

           ! Satellite output frequency
           ! Note that we broke with tradition here so that the
           ! dt_output will always we read!  This may be changed
           ! in later distributions
           call read_var('DnOutput',dn_output(ifile))
           call read_var('DtOutput',dt_output(ifile))

           ! Satellite inputfile name or the satellite name
           call read_var('NameTrajectoryFile',&
                Satellite_name(ifile-satellite_))
           if(index(satellite_string,'eqn')>0 &
                .or. index(satellite_string,'Eqn')>0 .or. &
                index(satellite_string,'EQN')>0 ) then
              UseSatelliteFile(ifile-satellite_) = .false.
           else
              UseSatelliteFile(ifile-satellite_) = .true.
           end if

           ! Satellite variables
           if(index(satellite_string,'VAR')>0 .or. &
                index(satellite_string,'var')>0 )then
              satellite_var='var'
              plot_dimensional(ifile)= index(satellite_string,'VAR')>0
              sat_time(ifile) = 'step date'
              call read_var('NameSatelliteVars',satellite_vars(ifile))
           elseif(index(satellite_string,'MHD')>0 .or. &
                index(satellite_string,'mhd')>0)then
              satellite_var='mhd'
              plot_dimensional(ifile)= index(satellite_string,'MHD')>0
              sat_time(ifile) = 'step date'
              satellite_vars(ifile)='rho ux uy uz bx by bz p jx jy jz'
           elseif(index(satellite_string,'FUL')>0 .or. &
                index(satellite_string,'ful')>0)then
              satellite_var='ful'
              plot_dimensional(ifile)= index(satellite_string,'FUL')>0
              sat_time(ifile) = 'step date'
              satellite_vars(ifile)=&
                   'rho ux uy uz bx by bz b1x b1y b1z p jx jy jz'
           else
              call stop_mpi(&
                   'Variable definition (mhd,ful,var) missing' &
                   //' from satellite_string='//satellite_string)
           end if
           plot_type(ifile) = "satellite"

           ! Determine the time output format to use in the 
           ! satellite files.  This is loaded by default above, 
           ! but can be input in the log_string line.
           if(index(satellite_string,'none')>0) then
              sat_time(ifile) = 'none'
           elseif((index(satellite_string,'step')>0) .or. &
                (index(satellite_string,'date')>0) .or. &
                (index(satellite_string,'time')>0)) then
              sat_time(ifile) = ''
              if(index(satellite_string,'step')>0) &
                   sat_time(ifile) = 'step'
              if(index(satellite_string,'date')>0) &
                   write(sat_time(ifile),'(a)') &
                   sat_time(ifile)(1:len_trim(sat_time(ifile)))&
                   //' date'
              if(index(satellite_string,'time')>0) &
                   write(sat_time(ifile),'(a)') &
                   sat_time(ifile)(1:len_trim(sat_time(ifile)))&
                   //' time'
           end if

        end do
     case('#RESCHANGEBOUNDARY')
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('SaveBoundaryCells',SaveBoundaryCells)
     case('#VERTEXBASEDGRID')          !^CFG IF COVARIANT BEGIN
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('UseVertexBasedGrid',UseVertexBasedGrid)
     case("#COVARIANTGEOMETRY")          
        if(.not.is_first_session())CYCLE READPARAM
        UseCovariant=.true.
        
        call read_var('TypeGeometry',TypeGeometry)      
        if(is_axial_geometry().and.mod(proc_dims(2),2)==1)&
             proc_dims(2)=2*proc_dims(2)
        if(index(TypeGeometry,'spherical')>0)then      
           automatic_refinement=.false.                  
           MaxBoundary = Top_
           DoFixExtraBoundary=.true.
           SaveBoundaryCells=.true.
        end if
        if(index(TypeGeometry,'cylindrical')>0)then
           automatic_refinement=.false.
           MaxBoundary = max(MaxBoundary,North_)
           DoFixExtraBoundary=.true.
           SaveBoundaryCells=.true.
        end if
    case("#LIMITGENCOORD1")                    
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('XyzMin_D(1)',XyzMin_D(1))
        call read_var('XyzMax_D(1)',XyzMax_D(1))
     case('#TORUSSIZE')
        call read_var('rTorusLarge',rTorusLarge)
        call read_var('rTorusSmall',rTorusSmall)
        !^CFG END COVARIANT 
     case("#GRID")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('nRootBlockX',proc_dims(1)) 
        call read_var('nRootBlockY',proc_dims(2))
        call read_var('nRootBlockZ',proc_dims(3))
        !^CFG IF COVARIANT BEGIN
        if( is_axial_geometry()&
             .and.mod(proc_dims(Phi_),2)==1) then
           proc_dims(Phi_)=2*proc_dims(3)
           if(iProc==0)then
              write(*,*)&
                   'Original number of blocks for the polar angle Phi must be even'
              write(*,*)&
                   'Proc_dim(2) is set to be equal to 2*Proc_dim(3) = ',&
                   proc_dims(2)
           end if
        end if
        !^CFG END COVARIANT
        if(product(proc_dims)>nBLK.and.iProc==0)then
           write(*,*)'Root blocks will not fit on 1 processor, check nBLK'
           call stop_mpi('product(proc_dims) > nBLK!')
        end if
        call read_var('xMin',x1)
        call read_var('xMax',x2)
        call read_var('yMin',y1)
        call read_var('yMax',y2)
        call read_var('zMin',z1)
        call read_var('zMax',z2)
        
        call set_xyzminmax  
        

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
        if(any(nCells/=nIJKRead_D).and.iProc==0)then
           write(*,*)'Code is compiled with nI,nJ,nK=',nCells
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
        gm1     = g-cOne
        gm2     = g-cTwo
        gp1     = g+cOne
        inv_g   = cOne/g
        inv_gm1 = cOne/gm1
        g_half  = cHalf*g
     case('#USERBOUNDARY', '#EXTRABOUNDARY')
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('UseExtraBoundary',UseExtraBoundary)
        if(UseExtraBoundary) call read_var('TypeBc_I(ExtraBc_)',&
             TypeBc_I(ExtraBc_))      
        
        if(.not.is_axial_geometry())&             !^CFG IF COVARIANT
             call read_var('DoFixExtraBoundary',&  
             DoFixExtraBoundary)  

     case('#FACEOUTERBC')                      
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('MaxBoundary',MaxBoundary)
        !^CFG IF COVARIANT BEGIN
        select case(TypeGeometry)
        case('spherical','spherical_lnr')
           MaxBoundary = Top_
        case('cylindrical')
           MaxBoundary = max(MaxBoundary,North_)
        end select
        !^CFG END COVARIANT
        if(MaxBoundary>=East_)&
             call read_var('DoFixOuterBoundary',DoFixOuterBoundary) 
    
     case("#SOLARWIND")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('SwRhoDim',SW_rho_dim)
        call read_var('SwTDim'  ,SW_T_dim)
        call read_var('SwUxDim' ,SW_Ux_dim)
        call read_var('SwUyDim' ,SW_Uy_dim)
        call read_var('SwUzDim' ,SW_Uz_dim)
        call read_var('SwBxDim' ,SW_Bx_dim)
        call read_var('SwByDim' ,SW_By_dim)
        call read_var('SwBzDim' ,SW_Bz_dim)
     case("#MAGNETOSPHERE","#BODY")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('UseBody',body1)
        if(body1)then
           call read_var('rBody'     ,Rbody)
           if(NameThisComp=='GM')then
              call read_var('rCurrents' ,Rcurrents)
              call read_var('BodyRhoDim',Body_Rho_Dim)
              call read_var('BodyTDim'  ,Body_T_dim)
           end if
        end if
     case("#GRAVITY")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('UseGravity',UseGravity)
        if(UseGravity)call read_var('iDirGravity',GravityDir)
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
        end if
     case("#DIPOLEBODY2")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('BdpDimBody2x',BdpDimBody2_D(1))
        call read_var('BdpDimBody2y',BdpDimBody2_D(2))
        call read_var('BdpDimBody2z',BdpDimBody2_D(3))
        !                                           ^CFG END SECONDBODY

     case('#PLANET','#IDEALAXES','#ROTATIONAXIS','#MAGNETICAXIS',&
          '#ROTATION','#DIPOLE','#NONDIPOLE','#UPDATEB0')

        call check_stand_alone
        if(.not.is_first_session())CYCLE READPARAM

        call read_planet_var(NameCommand)

     case("#COORDSYSTEM","#COORDINATESYSTEM")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('TypeCoordSystem',TypeCoordSystem)
        call upper_case(TypeCoordSystem)
        select case(NameThisComp)
        case('GM')
           if(TypeCoordSystem /= 'GSM')call stop_mpi(NameSub// &
                ' ERROR: cannot handle coordinate system '&
                //TypeCoordSystem)
        case('IH')
           select case(TypeCoordSystem)
           case('HGI')
              DoTransformToHgi = UseRotatingFrame
              UseRotatingFrame = .false.
           case('HGC')
              UseRotatingFrame = .true.
           case('HGR')
              if(iProc==0)then
                 write(*,*) NameSub, &
                      ' WARNING: corotating IH does not fully work'
                 if(UseStrict)call stop_mpi('Correct PARAM.in!')
                 write(*,*)NameSub//' setting .UseRotatingFrame = T'
              end if
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
     case("#HELIOSPHERE")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('BodyTDim' ,Body_T_dim)
        call read_var('BodyRhoDim',Body_rho_dim)
        call read_var('qSun' ,Qsun)
        call read_var('tHeat',Theat)
        call read_var('rHeat',Rheat)
        call read_var('SigmaHeat',SIGMAheat)
        call read_var('UseFluxRope',UseFluxRope)
        if(UseFluxRope)then
           call read_var('CmeA' ,cme_a)
           call read_var('CmeR1',cme_r1)
           call read_var('CmeR0',cme_r0)
           call read_var('CmeA1',cme_a1)
           call read_var('CmeAlpha',cme_alpha)
           call read_var('CmeRho1',cme_rho1)
           call read_var('CmeRho2',cme_rho2)
           call read_var('ModulationRho',ModulationRho)
           call read_var('ModulationP',ModulationP)
           call read_var('OrientationGl98',OrientationGL98)
           call read_var('LatitudeGl98   ',LatitudeGL98)
           call read_var('LongitudeGl98  ',LongitudeGL98)
        end if
     case("#HELIOUPDATEB0")
        call read_var('DtUpdateB0',dt_updateb0)
        DoUpdateB0 = dt_updateb0 > 0.0
     case("#HELIODIPOLE")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('HelioDipoleStrength',Bdp_dim)
        call read_var('HelioDipoleTilt'    ,ThetaTilt)
        ThetaTilt = ThetaTilt * cDegToRad
     case("#HELIOROTATION", "#INERTIAL")
        if(iProc==0)write(*,*) NameSub, ' WARNING: ',&
             ' #HELIOROTATION / #INERTIAL command is obsolete and ignored'
     case("#HELIOTEST")
        call read_var('DoSendMHD',DoSendMHD)
     case("#HELIOBUFFERGRID")
        if(.not.is_first_session())CYCLE READPARAM
        if(NameThisComp /= "IH")call stop_mpi(NameSub//' ERROR:'// &
             ' #HELIOBUFFERGRID command can be used in IH component only')
        call read_var('rBuffMin',  rBuffMin)
        call read_var('rBuffMax',  rBuffMax)
        call read_var('nThetaBuff',nThetaBuff)
        call read_var('nPhiBuff',  nPhiBuff)
     case("#ARCADE")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('TArcDim'  ,TArcDim)
        call read_var('RhoArcDim',RhoArcDim)
        call read_var('BArcDim'  ,BArcDim)
        call read_var('ByArcDim' ,ByArcDim)
        call read_var('UzArcDim' ,UzArcDim)
        call read_var('phi0Arc'  ,phi0Arc)
        call read_var('muArc'    ,muArc)
        call read_var('expArc'   ,expArc)
        call read_var('widthArc' ,widthArc)
     case("#CME")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('Cmetype',cme_type)
        call read_var('Cmea' ,cme_a)
        call read_var('Cmer1',cme_r1)
        call read_var('Cmer0',cme_r0)
        call read_var('Cmea1',cme_a1)
        call read_var('Cmealpha',cme_alpha)
        call read_var('Cmerho1',cme_rho1)
        call read_var('Cmerho2',cme_rho2)
        call read_var('CmeB1_dim',cme_B1_dim)
        call read_var('Cmev_erupt',cme_v_erupt)
     case("#TESTDISSMHD")                         !^CFG IF DISSFLUX BEGIN
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('UseDefaultUnits',UseDefaultUnits)
        call read_var('Grav0Diss'      ,Grav0Diss)
        call read_var('Beta0Diss'      ,Beta0Diss)
        call read_var('Length0Diss'    ,Length0Diss)
        call read_var('Time0Diss'      ,Time0Diss)
        call read_var('Rho0Diss'       ,Rho0Diss)
        call read_var('Tem0Diss'       ,Tem0Diss)
        call read_var('ThetaDiss'      ,ThetaDiss)
        call read_var('DeltaDiss'      ,DeltaDiss)
        call read_var('EpsilonDiss'    ,EpsilonDiss)
        call read_var('RhoDifDiss'     ,RhoDifDiss)
        call read_var('yShiftDiss'     ,yShiftDiss)
        call read_var('scaleHeightDiss',scaleHeightDiss)
        call read_var('scaleFactorDiss',scaleFactorDiss)
        call read_var('BZ0Diss'        ,BZ0Diss)
        !                                             ^CFG END DISSFLUX
     case("#COMET")
        if(.not.is_first_session())CYCLE READPARAM
        call read_var('Qprod' ,Qprod)
        call read_var('Unr_in' ,Unr_in)
        call read_var('mbar',mbar)
        call read_var('ionization_rate',ionization_rate)
        call read_var('kin_in',kin_in)
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

  subroutine set_defaults

    CodeVersionRead = -1.

    problem_type = -1

    ! Default coordinate systems
    select case(NameThisComp)
    case('IH')
       TypeCoordSystem = 'HGI'
       UseRotatingFrame     = .false.
    case('SC')
       TypeCoordSystem = 'HGR'
       UseRotatingFrame     = .true.
    case('GM')
       TypeCoordSystem = 'GSM'
       UseRotatingFrame     = .false.
    end select

    ! Do not update B0 for SC or IH by default
    if(NameThisComp /= 'GM')then
       DoUpdateB0 = .false.
       Dt_UpdateB0 = -1.0
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

    UsePartLocal  = .false.       !^CFG IF IMPLICIT
    nSTAGE        = 2
    cfl           = 0.80
    UseDtFixed    = .false.
    dt            = 0.0
    dt_BLK        = 0.0

    initial_refine_levels = 4   
    nRefineLevelIC        = 0

    min_block_level =  0
    max_block_level = 99
    min_cell_dx =     0.
    max_cell_dx = 99999.
    fix_body_level = .false.

    percentCoarsen = 0.
    percentRefine  = 0.
    maxTotalBlocks = nBLK*nProc

    dn_refine=-1
    automatic_refinement = .false.

    nOrder = 2
    FluxType = 'Rusanov'               !^CFG IF RUSANOVFLUX
    !FluxType = 'Sokolov'              !^CFG UNCOMMENT IF NOT RUSANOVFLUX
    !UseCovariant = .true.             !^CFG UNCOMMENT IF COVARIANT
  
    ! Default implicit parameters      !^CFG IF IMPLICIT BEGIN
    UsePointImplicit = .false.
    UsePartImplicit  = .false.
    UseFullImplicit  = .false.
    UseImplicit      = .false.
    ImplCritType     = 'dt'

    nOrder_impl   = 1
    FluxTypeImpl  = 'default'

    ImplCoeff0    = 1.0
    UseBDF2       = .false.
    ImplSource    = .false.

    UseConservativeImplicit = .false.
    UseNewton     = .false.
    NewMatrix     = .true.
    NewtonIterMax = 1

    JacobianType  = 'prec'
    if(nByteReal>7)then
       JacobianEps   = 1.E-12
    else
       JacobianEps   = 1.E-6
    end if

    PrecondSide   = 'symmetric'
    PrecondType   = 'MBILU'
    GustafssonPar = 0.5

    KrylovType      = 'gmres' !do not rename the gmres string
    KrylovInitType  = 'nul'
    KrylovErrorMax  = 0.001
    KrylovMatvecMax = 100
    nKrylovVector   = KrylovMatvecMax !^CFG END IMPLICIT

    UseDivbSource   = .true.
    UseDivbDiffusion= .false.         !^CFG IF DIVBDIFFUSE
    UseProjection   = .false.         !^CFG IF PROJECTION
    UseConstrainB   = .false.         !^CFG IF CONSTRAINB

    UseB0Source     = .true.

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

    UseBorisSimple = .false.          !^CFG IF SIMPLEBORIS
    boris_correction = .false.        !^CFG IF BORISCORR
    boris_cLIGHT_factor = 1.0         

    proc_dims(1) = 1
    proc_dims(2) = 1  
    proc_dims(3) = 1

    x1 = -10.
    x2 =  10.
    y1 = -10.
    y2 =  10.
    z1 = -10.
    z2 =  10.

    call set_xyzminmax

    optimize_message_pass = 'allopt'

    UseUpstreamInputFile = .false.

    UseMassLoading        = .false.
    AccelerateMassLoading = .false.

    plot_dimensional      = .true.
    save_satellite_data   = .false.

    restart           = .false.
    read_new_upstream = .false.

  end subroutine set_defaults

  !===========================================================================
  subroutine set_problem_defaults


    !\
    ! Give some "reasonable" default values
    !/

    UseGravity=.false.
    GravityDir=0
    select case(problem_type)
    case( problem_cme, problem_saturn, problem_jupiter, &
         problem_venus, problem_mars)
       UseGravity=.true.
    case(problem_heliosphere)
       UseGravity=.true.    
    case(problem_arcade)
       UseGravity=.true.
       GravityDir=z_
    end select

    select case(problem_type)
    case(problem_earth)
       body1      =.true.
       Rbody      = 3.00
       Rcurrents  = 4.00

    case(problem_saturn, problem_jupiter)
       body1      =.true.
       Rbody      = 3.00
       Rcurrents  = 4.00

    case(problem_rotation)
       body1      =.false.
       Rbody      = 3.00
       Rcurrents  = 0.00

    case(problem_comet)
       body1      =.false.
       Rbody      = -1.00
       Rcurrents  = -1.00

    case(problem_heliosphere, problem_venus, &
         problem_cylinder, problem_sphere)
       body1      =.true.
       Rbody      = 1.00
       Rcurrents  =-1.00

    case(problem_mars)
       body1      =.true.
       Rbody      = 1.041225
       Rcurrents  =-1.00

    case(problem_globalhelio)
       body1= .true.
       Rbody = 30.00
       Rcurrents= -1.0

    case(problem_cme)
       body1      =.true.
       UseBody2      =.false.           !^CFG IF SECONDBODY
       Rbody      = 1.0
       Rcurrents  = -1.0

    case(problem_shocktube)
       shock_Lstate(rho_)   =1.
       shock_Lstate(ux_:uz_)=0.
       shock_Lstate(Bx_)    =0.75
       shock_Lstate(By_)    =1.
       shock_Lstate(Bz_)    =0.
       shock_Lstate(P_)     =1.
       shock_Rstate(rho_)   =0.125
       shock_Rstate(ux_:uz_)=0.
       shock_Rstate(Bx_)    =0.75
       shock_Rstate(By_)    =-1.
       shock_Rstate(Bz_)    =0.
       shock_Rstate(P_)     =0.1
       ShockSlope           =0.0
    case default
       body1      =.false.
       Rbody      = -1.0
       Rcurrents  = -1.0

    end select

    select case(problem_type)                          !^CFG IF NOT COVARIANT BEGIN
    case(problem_earth,problem_saturn,problem_jupiter,problem_rotation)
       UseNonConservative   = .true.
       nConservCrit         = 1
       if(allocated(TypeConservCrit_I)) deallocate(TypeConservCrit_I)
       allocate( TypeConservCrit_I(nConservCrit) )
       TypeConservCrit_I(1) = 'r'
       rConserv             = 2*rBody
    case default                                       !^CFG END COVARIANT
       UseNonConservative   = .false.
       if(allocated(TypeConservCrit_I)) deallocate(TypeConservCrit_I)
       nConservCrit         = 0
       rConserv             = -1.
    end select                                         !^CFG IF NOT COVARIANT

    if(IsStandAlone)then
       select case(problem_type)
       case(problem_earth)
          Bdp_dim = -31100.0
       case(problem_saturn)
          Bdp_dim = 20800.0
       case(problem_jupiter)
          Bdp_dim = 426000.0
       case default
          Bdp_dim = 0.0
       end select
    end if

    select case(problem_type)
    case(problem_earth, problem_saturn, problem_jupiter, problem_venus, &
         problem_mars, problem_sphere,problem_comet)
       TypeBc_I(east_)        ='outflow'
       TypeBc_I(west_)        ='inflow'
!!$     TypeBc_I(south_:top_)  ='float'
       TypeBc_I(south_:top_)  ='fixed'
    case(problem_rotation)
       TypeBc_I(east_:top_)   ='fixedB1'
    case(problem_diffusion)
       TypeBc_I(east_:top_)   ='float'
    case(problem_cylinder)
       TypeBc_I(east_)        ='float'
       TypeBc_I(west_)        ='inflow'
       TypeBc_I(south_:north_)='fixed'
       TypeBc_I(bot_:top_)    ='float'
    case(problem_shocktube)
       TypeBc_I(east_:north_) ='shear'
       TypeBc_I(bot_:top_)    ='float'
    case default
       TypeBc_I(east_:top_)   ='float'
    end select

    select case(problem_type)
    case(problem_earth, problem_saturn, problem_jupiter, problem_rotation)
       TypeBc_I(body1_)='ionosphere'
    case default
       TypeBc_I(body1_)='unknown'
    end select

    ! These are used for problem_heliosphere only
    Qsun  = 25.00
    Theat = 1.75
    Rheat = 1.00
    SIGMAheat = 4.50

    ! These are used by problem_cme only
    cme_type    = 'Low'
    cme_a       = 0.7
    cme_r1      = 1.2
    cme_r0      = 1.0
    cme_a1      = 0.23
    cme_alpha   = 0.0
    cme_rho1    = 2.5E-12
    cme_rho2    = 2.0E-13
    cme_B1_dim  = 1.0 
    cme_v_erupt = 4.0E5


    ! These are used for problem_arcade
    TArcDim  = 1.0E6
    RhoArcDim= 1.0E-12
    BArcDim  = 0.71814
    ByArcDim = 0.0
    UzArcDim = 5.0E3
    phi0Arc  = 0.5
    muArc    = 1.3
    expArc   = 3
    widthArc = 0.5

    ! These are used for problem_comet
    Qprod = 1.0E28           ! 1/s 
    Unr   = 1.0              ! km/s
    mbar  = 17.0             ! AMU
    ionization_rate = 1.0E-6 ! 1/s
    kin   = 1.7E-9           ! cm^3/s

    SW_rho_dim =     -1.0    ! n/cc
    SW_T_dim   =     -1.0    ! K
    SW_Ux_dim  =      0.0    ! km/s
    SW_Uy_dim  =      0.0    ! km/s
    SW_Uz_dim  =      0.0    ! km/s
    SW_Bx_dim  =      0.0    ! nT
    SW_By_dim  =      0.0    ! nT
    SW_Bz_dim  =      0.0    ! nT

    select case (problem_type)
    case (problem_heliosphere, problem_cme)
       Body_rho_dim = 1.50E8
       Body_T_dim   = 2.85E06
    case (problem_earth)
       Body_rho_dim = 5.0    ! n/cc
       Body_T_dim   = 25000.0! K
    case (problem_saturn)
       Body_rho_dim = 0.1    ! n/cc
       Body_T_dim   = 35000.0! K
    case (problem_jupiter)
       Body_rho_dim = 0.4    ! n/cc
       Body_T_dim   = 35000.0! K
    case (problem_rotation)
       Body_rho_dim = 25.0    ! n/cc
       Body_T_dim   = 181712.175   ! K
       !^CFG IF GLOBALHELIOSPHERE BEGIN
    case(problem_globalhelio)
       Body_rho_dim = SW_rho*unitUSER_rho*1.E6
       Body_T_dim = SW_T_dim
       !^CFG END GLOBALHELIOSPHERE
    case default
       Body_rho_dim = 1.0    ! n/cc
       Body_T_dim   = 10000.0! K
    end select

    UseBody2 = .false.                                !^CFG IF SECONDBODY BEGIN
    RBody2 =-1.0
    xBody2 = 0.0
    yBody2 = 0.0			
    zBody2 = 0.0			
    BdpBody2_D  = 0.0
    rCurrentsBody2 = 0.0
    RhoDimBody2 = 5.0    ! n/cc
    TDimBody2   = 25000.0! K                          !^CFG END SECONDBODY

    call set_InitialRefineType

    !
    ! Available criteria values to choose from: (SELECT UP TO 3)
    !    dTempmax drhomax dPmax dP2max dEmax curlUmax curlBmax
    !    divUmax divBmax Vamax betamax fluxmax Zmax
    !
    select case (problem_type)
    case (problem_heliosphere)
       nRefineCrit    = 3
       RefineCrit(1)  = 'geometry'
       RefineCrit(2)  = 'Va'
       RefineCrit(3)  = 'flux'
       !^CFG IF GLOBALHELIOSPHERE BEGIN
    case (problem_globalhelio)
       nRefineCrit    = 3
       RefineCrit(1)  = 'geometry'
       RefineCrit(2)  = 'Va'
       RefineCrit(3)  = 'flux'
       !^CFG END GLOBALHELIOSPHERE
    case (problem_earth,problem_saturn,problem_jupiter)
       nRefineCrit    = 3
       RefineCrit(1)  = 'gradlogP'
       RefineCrit(2)  = 'curlB'
       RefineCrit(3)  = 'Rcurrents'
    case (problem_comet)
       nRefineCrit    = 3
       RefineCrit(1)  = 'gradlogP'
       RefineCrit(2)  = 'divU'
       RefineCrit(3)  = 'curlB'
    case default
       nRefineCrit    = 3
       RefineCrit(1)  = 'gradlogP'
       RefineCrit(2)  = 'divB'
       RefineCrit(3)  = 'curlB'
    end select

  end subroutine set_problem_defaults

  !==========================================================================
  subroutine set_InitialRefineType

    select case (problem_type)
    case (problem_uniform)
       InitialRefineType = 'none'
    case (problem_shocktube)
       InitialRefineType = 'all'
    case (problem_rotation)
       InitialRefineType = '3Dbodyfocus'
    case (problem_heliosphere)
       InitialRefineType = 'all_then_focus'
    case (problem_cme)
       InitialRefineType = 'cme'
    case (problem_earth)
       InitialRefineType = 'magnetosphere'
    case (problem_saturn)
       InitialRefineType = 'magnetosaturn'
    case (problem_jupiter)
       InitialRefineType = 'magnetojupiter'
    case (problem_venus)
       InitialRefineType = 'magnetosphere'
    case (problem_mars)
       InitialRefineType = 'magnetosphere'
    case (problem_cylinder)
       InitialRefineType = 'all'
    case (problem_sphere)
       InitialRefineType = 'none'
    case (problem_comet)
       InitialRefineType = 'comet'
    case default
       InitialRefineType = 'none'
    end select

  end subroutine set_InitialRefineType

  !=========================================================================
  subroutine check_options

    ! option and module parameters
    character (len=40) :: Name
    real               :: Version
    logical            :: IsOn

    if(UseCovariant)then   !^CFG IF COVARIANT BEGIN
       call allocate_face_area_vectors
       if(UseVertexBasedGrid) call allocate_old_levels
       if(UseImplicit)call stop_mpi(&                  !^CFG IF IMPLICIT 
            'Do not use covariant with implicit')      !^CFG IF IMPLICIT
       if(UseProjection)call stop_mpi(&                !^CFG IF PROJECTION
            'Do not use covariant with projection')    !^CFG IF PROJECTION
       if(UseConstrainB)call stop_mpi(&                !^CFG IF CONSTRAINB
            'Do not use covariant with constrain B')   !^CFG IF CONSTRAINB
       if(UseRaytrace) call stop_mpi(&                 !^CFG IF RAYTRACE
            'Do not use covariant with ray tracing')   !^CFG IF RAYTRACE
       if(UseDivBDiffusion)call stop_mpi(&             !^CFG IF DIVBDIFFUSE
            'Do not use covariant with divB diffusion')!^CFG IF DIVBDIFFUSE
       if(boris_correction)call stop_mpi(&             !^CFG IF BORISCORR BEGIN
            'Do not use covariant with the Boris correction')
                                                       !^CFG END BORISCORR
    end if                                             !^CFG END COVARIANT
    if(SaveBoundaryCells)then
       if(index(optimize_message_pass,'all')>0)then
          call allocate_boundary_cells
       else
          if(iProc==0)&
               write(*,'(a)')NameSub//&
               'SaveBoundaryCells does not work '&
               //'with message_pass option='//optimize_message_pass//&
               ', set SaveBoundaryCells to F'
          SaveBoundaryCells=.false.
       end if
    end if
    if(UseImplicit)then                      !^CFG IF IMPLICIT BEGIN
       call OPTION_IMPLICIT(IsOn,Name)
       if(.not.IsOn)then
          if(iProc==0)then
             write(*,'(a)')NameSub//' WARNING: IMPLICIT module is OFF !!!'
             if(UseStrict) &
                  call stop_mpi('Correct PARAM.in or switch IMPLICIT on!')
             write(*,*)NameSub//' setting UseImplicit=.false.'
          end if
          UseImplicit=.false.
       end if
    end if                                   !^CFG END IMPLICIT

    if(UseProjection)then                    !^CFG IF PROJECTION BEGIN
       call OPTION_PROJECTION(IsOn,Name)
       if(.not.IsOn)then
          if(iProc==0)then
             write(*,'(a)')NameSub//' WARNING: PROJECTION module is OFF !!!'
             if(UseStrict)&
                  call stop_mpi('Correct PARAM.in or switch PROJECTION on!')
             write(*,*)NameSub//&
                  ' setting UseProjection=.false. and UseDivbSource=.true.'
          end if
          UseProjection=.false.
          UseDivbSource=.true.
       end if
    end if                                   !^CFG END PROJECTION

    if(UseConstrainB)then                    !^CFG IF CONSTRAINB BEGIN
       call OPTION_CONSTRAIN_B(IsOn,Name)
       if(.not.IsOn)then
          if(iProc==0)then
             write(*,'(a)')NameSub//' WARNING: CONSTRAINB module is OFF !!!'
             if(UseStrict)&
                  call stop_mpi('Correct PARAM.in or switch CONSTRAINB on!')
             write(*,*)NameSub//&
                  ' setting UseConstrainB=.false. and UseDivbSource=.true.'
          end if
          UseConstrainB=.false.
          UseDivbSource=.true.
       end if
    end if                                   !^CFG END CONSTRAINB

    !^CFG  IF RAYTRACE BEGIN
    if(UseRaytrace .or. any(index(plot_type,'ray')>0))then
       call OPTION_RAYTRACING(IsOn,Name)
       if(.not.IsOn)then
          if(iProc==0)write(*,'(a)')NameSub// &
               ' WARNING: RAYTRACING module is OFF !!!'
          if(UseStrict)&
               call stop_mpi('Correct PARAM.in or switch RAYTRACING on!')
          if(UseRaytrace)then
             if(iProc==0)write(*,*)NameSub//' setting UseRaytrace=.false.'
             UseRaytrace=.false.
          end if
          if(any(index(plot_type,'ray')>0))then
             if(iProc==0)write(*,*)'setting plot_type=nul'
             where(index(plot_type,'ray')>0)plot_type='nul'
          end if
       end if
    end if
    !^CFG END RAYTRACE

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

  end subroutine check_options

  !=========================================================================
  subroutine correct_parameters

    !\
    ! Check for some combinations of things that cannot be accepted as input
    !/
    if (iProc==0) write (*,*) ' '

    ! Check flux type selection
    select case(FluxType)
    case('1','roe','Roe','ROE')                      !^CFG IF ROEFLUX
       FluxType='Roe'                                !^CFG IF ROEFLUX
    case('2','rusanov','Rusanov','RUSANOV','TVDLF')  !^CFG IF RUSANOVFLUX
       FluxType='Rusanov'                            !^CFG IF RUSANOVFLUX
    case('3','linde','Linde','LINDE','HLLEL')        !^CFG IF LINDEFLUX
       FluxType='Linde'                              !^CFG IF LINDEFLUX
    case('4','sokolov','Sokolov','SOKOLOV','AW')     !^CFG IF AWFLUX
       FluxType='Sokolov'                            !^CFG IF AWFLUX
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

    ! Check flux type selection for implicit   !^CFG IF IMPLICIT BEGIN
    select case(FluxTypeImpl)
    case('default')
       FluxTypeImpl = FluxType
    case('1','roe','Roe','ROE')                      !^CFG IF ROEFLUX
       FluxTypeImpl='Roe'                            !^CFG IF ROEFLUX
    case('2','rusanov','Rusanov','RUSANOV','TVDLF')  !^CFG IF RUSANOVFLUX
       FluxTypeImpl='Rusanov'                        !^CFG IF RUSANOVFLUX
    case('3','linde','Linde','LINDE','HLLEL')        !^CFG IF LINDEFLUX
       FluxTypeImpl='Linde'                          !^CFG IF LINDEFLUX
    case('4','sokolov','Sokolov','SOKOLOV','AW')     !^CFG IF AWFLUX
       FluxTypeImpl='Sokolov'                        !^CFG IF AWFLUX
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

    if (time_accurate .and. nStage==2) UseBDF2=.true. !^CFG IF IMPLICIT

    ! Make sure periodic boundary conditions are symmetric
    if(any(TypeBc_I(1:2)=='periodic')) TypeBc_I(1:2)='periodic'
    if(any(TypeBc_I(3:4)=='periodic')) TypeBc_I(3:4)='periodic'
    if(any(TypeBc_I(5:6)=='periodic')) TypeBc_I(5:6)='periodic'

    ! Reset initial refine type if it is set to 'default'
    if(InitialRefineType=='default') call set_InitialRefineType
 
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
       if(UseStrict)call stop_mpi('Correct PARAM.in!')
       if(iProc==0)write(*,*)NameSub// &
            ' setting UseTvdReschange=F UseAccurateResChange=F'
       UseTvdReschange      = .false.
       UseAccurateResChange = .false.
    end if
    if (UseConstrainB .and. index(optimize_message_pass,'opt') > 0) then
       if(iProc==0 .and. optimize_message_pass /= 'allopt') then
          write(*,'(a)')NameSub//&
               ' WARNING: constrain_B does not work for'// &
               ' optimize_message_pass='//trim(optimize_message_pass)//' !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' setting optimize_message_pass = all'
       end if
       optimize_message_pass = 'all'
    endif                                            !^CFG END CONSTRAINB

    if (UseHallResist .and. index(optimize_message_pass,'opt') > 0) then
       if(iProc==0 .and. optimize_message_pass /= 'allopt') then
          write(*,'(a)')NameSub//&
               ' WARNING: Hall resistivity does not work for'// &
               ' optimize_message_pass='//trim(optimize_message_pass)//' !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' setting optimize_message_pass = all'
       end if
       optimize_message_pass = 'all'
    endif                                     

    if(prolong_order/=1 .and. optimize_message_pass(1:3)=='all')&
         call stop_mpi(NameSub// &
         'The prolongation order=2 requires message_pass_dir')

    if(UseTvdResChange .and. &
         optimize_message_pass(1:3)=='all' .and. iCFExchangeType/=2)then
       if(iProc==0) write(*,'(a)') NameSub// &
            ' WARNING: TVD limiter at the resolution change' //&
            ' requires iCFExchangeType=2 in message_pass_cells'
       call stop_mpi('Correct PARAM.in or message_pass_cells')
    end if

    ! Check test processor
    if(ProcTest > nProc)then
       if(iProc==0) write(*,'(a)')NameSub//&
            ' WARNING: procTEST > nProc, setting procTEST=0 !!!'
       procTEST=0
    end if

  end subroutine correct_parameters

  !===========================================================================
  subroutine check_parameters

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

    if(UsePartImplicit .and. ImplCritType=='dt' .and.&
         (.not.time_accurate .or. .not.UseDtFixed))then
       if(iProc==0)then
          write(*,'(a)')'Part implicit scheme with ImplCritType=dt'
          write(*,'(a)')'requires time accurate run with fixed time step'
          call stop_mpi('Correct PARAM.in')
       end if
    end if

    if(nKrylovVector>KrylovMatvecMax)then
       if(iProc==0)then
          write(*,'(a)')'nKrylovVector>KrylovMatvecMax is useless!'
          write(*,*)'reducing nKrylovVector to KrylovMatvecMax'
       endif
       nKrylovVector = KrylovMatvecMax
    endif

    if(PrecondSide/='left'.and.PrecondSide/='symmetric'.and.&
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

    if(UsePartImplicit)then
       select case(optimize_message_pass)
       case('oldopt','allold')
          if(iProc==0)then
             write(*,'(2a)') NameSub//&
                  ' WARNING: PartImplicit scheme does not work with',&
                  ' TypeMessagePass=oldopt or allold !!!'
             if(UseStrict)call stop_mpi('Correct PARAM.in!')
             write(*,*) NameSub//' setting optimize_message_pass=allopt'
          end if
          optimize_message_pass='allopt'
       end select
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

    !Finish checks for implicit                     !^CFG END IMPLICIT

  end subroutine check_parameters

  !===========================================================================
  subroutine set_extra_parameters

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
          DtFixed = DtFixedDim/unitSI_t
          DtFixedOrig = DtFixed ! Store the initial setting
          dt = DtFixed
          cfl=1.0
       end if
    end if

    ! You have to have the normalizations first
    if (read_new_upstream)  call normalize_upstream_data
    
    ! Set MinBoundary and MaxBoundary
    if(UseBody2) MinBoundary=min(Body2_,MinBoundary)   !^CFG IF SECONDBODY
    if(body1) then   
       MinBoundary=min(Body1_,MinBoundary)
       MaxBoundary=max(Body1_,MaxBoundary)
    end if
    if(UseExtraBoundary) then                        
       MinBoundary=min(ExtraBc_,MinBoundary)
       MaxBoundary=max(ExtraBc_,MaxBoundary)
    end if
    MaxBoundary=min(MaxBoundary,Top_)
    MinBoundary=max(MinBoundary,body2_)

    DoOneCoarserLayer = .not. (UseTvdResChange .or. UseAccurateResChange)
    DoLimitMomentum = &                                !^CFG IF BORISCORR
         boris_correction .and. DoOneCoarserLayer      !^CFG IF BORISCORR

    if(UseConstrainB) then          !^CFG IF CONSTRAINB BEGIN
       jMinFaceX=0; jMaxFaceX=nJ+1
       kMinFaceX=0; kMaxFaceX=nK+1
       iMinFaceY=0; iMaxFaceY=nI+1
       kMinFaceY=0; kMaxFaceY=nK+1
       iMinFaceZ=0; iMaxFaceZ=nI+1
       jMinFaceZ=0; jMaxFaceZ=nJ+1
    end if                          !^CFG END CONSTRAINB
 
  end subroutine set_extra_parameters

  !==========================================================================
  subroutine read_area_radii

    real :: Radius1, Radius2

    ! Read inner and outer radii for areas shell and ring
    call read_var("Radius1", Radius1)
    call read_var("Radius2", Radius2)

    ! Set outer size of the area
    RadiusArea = max(Radius1, Radius2)
    ! Normalize inner radius to the outer size
    Area_I(nArea)%Radius1  = min(Radius1, Radius2) / RadiusArea

  end subroutine read_area_radii

end subroutine MH_set_parameters
