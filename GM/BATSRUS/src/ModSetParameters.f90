!^CFG COPYRIGHT UM
subroutine MH_set_parameters(TypeAction)

  ! Set input parameters for the Global Magnetosphere (GM) module

  use ModProcMH
  use ModMain
  use ModAdvance
  use ModGeometry, ONLY : init_mod_geometry, &
       TypeGeometry,iVolumeCounterBLK,iVolumeCounterI,& !^CFG IF NOT CARTESIAN
       x1,x2,y1,y2,z1,z2,XyzMin_D,XyzMax_D,MinBoundary,MaxBoundary
  use ModNodes, ONLY : init_mod_nodes
  use ModImplicit                                       !^CFG IF IMPLICIT
  use ModPhysics
  use ModProject                                        !^CFG IF PROJECTION
  use ModCT, ONLY : init_mod_ct, DoInitConstrainB       !^CFG IF CONSTRAINB
  use ModAMR
  use ModParallel, ONLY : UseCorners,proc_dims
  use ModRaytrace                                       !^CFG IF RAYTRACE
  use ModIO
  use ModCompatibility, ONLY: read_compatible_command, SetDipoleTilt
  use CON_planet,       ONLY: read_planet_var, check_planet_var
  use CON_axes,         ONLY: init_axes
  use ModUtilities,     ONLY: check_dir, fix_dir_name, upper_case

  use CON_planet,       ONLY: get_planet
  use ModTimeConvert,   ONLY: time_int_to_real, time_real_to_int
  use ModReadParam
  use ModMpi

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

  ! The name of the command
  character (len=lStringLine) :: NameCommand, StringLine
  logical                     :: DoEcho

  ! Variables to remember for multiple calls
  character (len=100) :: UpstreamFileName='???'

  ! Temporary variables
  integer :: nByteRealRead, nVarRead
  character (len=100) :: NameEquationRead

  character (len=50) :: plot_string,log_string,satellite_string
  character (len=3)  :: plot_area, plot_var, satellite_var, satellite_form
  integer :: problem_type_r, qtotal, nIJKRead_D(3)

  integer            :: TimingDepth=-1
  character (len=10) :: TimingStyle='cumm'

  ! Variables for checking/reading #STARTTIME command
  real (Real8_)         :: StartTimeCheck = -1.0_Real8_

  logical            :: UseSimple=.true.

  integer :: iSession, iPlotFile
  !-------------------------------------------------------------------------
  NameSub(1:2) = NameThisComp

  UseSimple = .false.              !^CFG IF NOT SIMPLE

  iSession = i_session_read()

  if(IsUninitialized)then
     call set_defaults
     IsUninitialized=.false.
  end if

  if(iSession>1)then
     restart=.false.           ! restart in session 1 only
     read_new_upstream=.false. ! upstream file reading in session 1 only
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
        if(UseNewAxes)then
           ! Check and set some planet variables (e.g. DoUpdateB0)
           call check_planet_var(iProc==0, time_accurate)

           ! Initialize axes
           call init_axes(StartTime)

           ! Obtain some planet parameters
           call get_planet(UseRotationOut = UseCorotation, &
                DoUpdateB0Out = DoUpdateB0, DtUpdateB0Out = Dt_UpdateB0)
        else
           DoUpdateB0 = time_accurate .and. Dt_updateB0 > 0.0 .and. &
                UseCorotation .and. .not.SetDipoleTilt
        endif
     end if

     !^CFG IF NOT SIMPLE BEGIN
     if(problem_type==-1) call stop_mpi(&
          NameSub//': ERROR problem type must be set in PARAM.in')
     !^CFG END SIMPLE

     call correct_parameters

     call check_parameters

     call check_plot_range

     call check_options

     ! initialize module variables
     call init_mod_advance
     call init_mod_geometry
     call init_mod_nodes
     if(UseConstrainB) call init_mod_ct        !^CFG IF CONSTRAINB
     if(UseImplicit)   call init_mod_implicit  !^CFG IF IMPLICIT
     if(UseIM)         call init_mod_raytrace  !^CFG IF RCM

     if(UseConstrainB) then          !^CFG IF CONSTRAINB BEGIN
        jMinFaceX=0; jMaxFaceX=nJ+1
        kMinFaceX=0; kMaxFaceX=nK+1
        iMinFaceY=0; iMaxFaceY=nI+1
        kMinFaceY=0; kMaxFaceY=nK+1
        iMinFaceZ=0; iMaxFaceZ=nI+1
        jMinFaceZ=0; jMaxFaceZ=nJ+1
     end if                          !^CFG END CONSTRAINB

     if (read_new_upstream) &
          call read_upstream_input_file(UpstreamFileName)

     call set_physics_constants

     call set_physics_parameters

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
  do
     if(.not.read_line(StringLine) )then
        IsLastRead = .true.
        EXIT
     end if

     if(.not.read_command(NameCommand)) CYCLE

     select case(NameCommand)
     case("#NEWPARAM")
        call check_stand_alone
        call read_var('UseNewParam',UseNewParam)
        call read_var('UseNewAxes', UseNewAxes)
        call read_var('DoTimeAccurate',time_accurate)
        call read_var('UseCorotation',UseCorotation)
     case("#DESCRIPTION")
        call check_stand_alone
     case("#BEGIN_COMP","#END_COMP")
        call check_stand_alone
        i = len_trim(NameCommand)
        if(StringLine(i+2:i+3) /= 'GM')&
             call stop_mpi(NameSub//' ERROR: the component is not GM in '// &
             trim(StringLine))
     case("#END")
        call check_stand_alone
        IslastRead=.true.
        EXIT
     case("#RUN")
        call check_stand_alone
        IslastRead=.false.
        EXIT
     case("#STOP")
        call check_stand_alone
        call read_var('nIter',nIter)
        if(UseNewParam .or. time_accurate) &
             call read_var('t_max',t_max)
     case("#CPUTIMEMAX")
        call check_stand_alone
        call read_var('cputime_max',cputime_max)
     case("#CHECKSTOPFILE")
        call check_stand_alone
        call read_var('check_stopfile',check_stopfile)
     case("#PROGRESS")
        call check_stand_alone
        call read_var('dn_progress1',dn_progress1)
        call read_var('dn_progress2',dn_progress2)
     case("#TIMEACCURATE")
        call check_stand_alone
        call read_var('DoTimeAccurate',time_accurate)
     case("#ECHO")
        call check_stand_alone
        call read_var('DoEcho',DoEcho)
        if(iProc==0)call read_echo_set(DoEcho)
     case("#TEST")
        call read_var('test_string',test_string)
        if(.not.UseNewParam .and. index(test_string,'read_inputs')>0)then
           DoEcho = .true.
           if(iProc==0)call read_echo_set(DoEcho)
        end if
     case("#TESTXYZ")
        coord_test=.true.
        UseTestCell=.true.
        call read_var('Xtest',Xtest)
        call read_var('Ytest',Ytest)
        call read_var('Ztest',Ztest)
     case("#TESTIJK")
        coord_test=.false.
        UseTestCell=.true.
        call read_var('Itest',Itest)
        call read_var('Jtest',Jtest)
        call read_var('Ktest',Ktest)
        call read_var('BLKtest',BLKtest)
        call read_var('PROCtest',PROCtest)
     case("#TESTVAR")
        call read_var('VARtest ',VARtest)
     case("#TESTDIM")
        call read_var('DIMtest ',DIMtest)
     case("#TESTTIME")
        call read_var('Itertest',ITERtest)
        if(UseNewParam .or. time_accurate) &
             call read_var('Ttest',Ttest)
     case("#STRICT")
        call read_var('UseStrict',UseStrict)
     case("#VERBOSE")
        call read_var('lVerbose',lVerbose)
     case("#DEBUG")
        call read_var('okdebug',okdebug)
        call read_var('ShowGhostCells',ShowGhostCells)
     case("#TIMING")
        call read_var('UseTiming',UseTiming)
        if(UseTiming)then
           call read_var('dnTiming',dn_timing)
           call read_var('TimingDepth',TimingDepth)
           call read_var('TimingStyle',TimingStyle)
        end if
     case("#OUTERBOUNDARY")
        !                                              ^CFG IF NOT SIMPLE BEGIN
        call read_var('TypeBc_I(east_)',TypeBc_I(east_))  
        call read_var('TypeBc_I(west_)',TypeBc_I(west_))
        call read_var('TypeBc_I(south_)',TypeBc_I(south_))
        call read_var('TypeBc_I(north_)',TypeBc_I(north_))
        call read_var('TypeBc_I(bot_)',TypeBc_I(bot_))
        call read_var('TypeBc_I(top_)',TypeBc_I(top_))    
        !                                              ^CFG END SIMPLE
     case("#INNERBOUNDARY")
        !                                              ^CFG IF NOT SIMPLE BEGIN
        call read_var('TypeBc_I(body1_)',TypeBc_I(body1_))
        !                                              ^CFG IF SECONDBODY BEGIN
        if(UseBody2) &                                      
             call read_var('TypeBc_I(body2_)',TypeBc_I(body2_)) 
        !                                              ^CFG END SECONDBODY
        !                                              ^CFG END SIMPLE
     case("#TIMESTEPPING")
        if(.not.UseNewParam)call read_var('DoTimeAccurate',time_accurate)
        call read_var('nStage',nSTAGE)
        call read_var('CflExpl',Cfl)
        ExplCfl = Cfl                                   !^CFG IF IMPLICIT
     case("#FIXEDTIMESTEP")
        call read_var('UseDtFixed',UseDtFixed)
        if(UseDtFixed)call read_var('DtFixedDim', DtFixedDim)
        !                                                ^CFG IF IMPLICIT BEGIN

     case("#PARTLOCAL")                                 !^CFG IF NOT SIMPLE
        call read_var('UsePartLocal',UsePartLocal)      !^CFG IF NOT SIMPLE
     case("#IMPLICIT")
        call read_var('UsePointImplicit', UsePointImplicit) !^CFG IF POINTIMPLICIT
        call read_var('UsePartImplicit',  UsePartImplicit)
        call read_var('UseFullImplicit',  UseFullImplicit)

        UseImplicit = UseFullImplicit .or. UsePartImplicit

        ! For implicit scheme it is better to use unsplit dB0/dt evaluation
        DoSplitDb0Dt = .not.UseImplicit

        if(UsePartImplicit  .and. UseFullImplicit) call stop_mpi(&
             'Only one of UsePartImplicit and UseFullImplicit can be true')

        !                                           ^CFG IF POINTIMPLICIT BEGIN
        if(UsePointImplicit .and. UseImplicit) call stop_mpi(&
             'Only one of Point-/Full-/Partimmplicit can be true')
        !                                           ^CFG END POINTIMPLICIT

        if(UseImplicit)call read_var('ImplCFL',ImplCFL)

        !                                             ^CFG IF NOT SIMPLE BEGIN
     case("#IMPLICITCRITERIA", "#STEPPINGCRITERIA")
        call read_var('ImplCritType',ImplCritType)
        select case(ImplCritType)
        case('R','r')
           call read_var('Rimplicit'   ,Rimplicit)
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
     case("#IMPLSCHEME")
        call read_var('nOrderImpl',nORDER_impl)
        call read_var('FluxTypeImpl',FluxTypeImpl)
     case("#IMPLSTEP")
        call read_var('ImplCoeff ',ImplCoeff0)
        call read_var('UseBDF2   ',UseBDF2)
        call read_var('ImplSource',implsource)
     case("#NEWTON")
        call read_var('UseConservativeImplicit',UseConservativeImplicit)
        call read_var('UseNewton',UseNewton)
        if(UseNewton)then
           call read_var('NewMatrix    ',NewMatrix)
           call read_var('NewtonIterMax',NewtonIterMax)
        endif
     case("#JACOBIAN")
        call read_var('JacobianType',JacobianType)
        call read_var('JacobianEps', JacobianEps)
     case("#PRECONDITIONER")
        call read_var('PrecondSide'  ,PrecondSide)
        call read_var('PrecondType'  ,PrecondType)
        call read_var('GustafssonPar',GustafssonPar)
        if(GustafssonPar<=0.)GustafssonPar=0.
        if(GustafssonPar>1.)GustafssonPar=1.
     case("#KRYLOV")
        call read_var('KrylovType'     ,KrylovType)
        call read_var('KrylovInitType' ,KrylovInitType)
        call read_var('KrylovErrorMax' ,KrylovErrorMax)
        call read_var('KrylovMatvecMax',KrylovMatvecMax)
        nKrylovVector = KrylovMatvecMax
     case("#KRYLOVSIZE")
        call read_var('nKrylovVector',nKrylovVector)
        !                                               ^CFG END SIMPLE
        !                                               ^CFG END IMPLICIT
        !                                               ^CFG IF DISSFLUX BEGIN
     case("#HEATFLUX")
        call read_var('UseHeatFlux'   ,UseHeatFlux)
        call read_var('UseSpitzerForm',UseSpitzerForm)
        if (.not.UseSpitzerForm) then
           call read_var('Kappa0Heat'  ,Kappa0Heat)
           call read_var('ExponentHeat',ExponentHeat)
        end if
     case("#RESISTIVEFLUX")
        call read_var('UseResistFlux' ,UseResistFlux)
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
           call read_var('Eta0AnomResist'       ,Eta0AnomResist)
           call read_var('EtaAnomMaxResist'     ,EtaAnomMaxResist)
           call read_var('ThresholdFactorResist',ThresholdFactorResist)
        end if
        !                                               ^CFG END DISSFLUX
     case("#SAVERESTART")
        call read_var('SaveRestartFile',save_restart_file)
        if(save_restart_file)then
           if(iProc==0)call check_dir(NameRestartOutDir)
           call read_var('dnOutput(restart_)',dn_output(restart_))
           if(UseNewParam .or. time_accurate) &
                call read_var('dt_output(restart_)',dt_output(restart_))
           nfile=max(nfile,restart_)
        end if
     case("#SAVELOGFILE")
        call read_var('SaveLogfile',save_logfile)
        if(save_logfile)then
           if(iProc==0)call check_dir(NamePlotDir)
           nfile=max(nfile,logfile_)
           call read_var('log_string',log_string)
           call read_var('dn_output(logfile_)',dn_output(logfile_))
           if(UseNewParam .or. time_accurate) &
                call read_var('dt_output(logfile_)',dt_output(logfile_))

           ! Log variables
           if((index(log_string,'VAR')>0) .or. (index(log_string,'var')>0))then
              plot_dimensional(logfile_) = index(log_string,'VAR')>0
              log_time='step time'
              call read_var('log_vars',log_vars)
           elseif((index(log_string,'RAW')>0) .or.(index(log_string,'raw')>0) )then
              plot_dimensional(logfile_) = index(log_string,'RAW')>0
              log_time='step time'
              log_vars='dt rho mx my mz bx by bz e Pmin Pmax'
           elseif(index(log_string,'MHD')>0 .or. index(log_string,'mhd')>0 )then
              plot_dimensional(logfile_) = index(log_string,'MHD')>0
              log_time='step date time'
              log_vars='rho mx my mz bx by bz e Pmin Pmax'
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
        call read_var('nplotfile',nplotfile)

        if(nPlotFile>0 .and. iProc==0)call check_dir(NamePlotDir)

        nfile=max(nfile,plot_+nplotfile)
        if (nfile > maxfile .or. nplotfile > maxplotfile)call stop_mpi(&
             'The number of ouput files is too large in #SAVEPLOT:'&
             //' nplotfile>maxplotfile .or. nfile>maxfile')
        do iFile=plot_+1,plot_+nplotfile

           call read_var('plot_string',plot_string)

           ! Plotting frequency
           call read_var('dn_output',dn_output(ifile))
           if(UseNewParam .or. time_accurate) &
                call read_var('dt_output',dt_output(ifile))

           ! Plotting area
           if(index(plot_string,'cut')>0)then
              plot_area='cut'
              call read_var('plot_range(x1)',plot_range(1,ifile))
              call read_var('plot_range(x2)',plot_range(2,ifile))
              call read_var('plot_range(y1)',plot_range(3,ifile))
              call read_var('plot_range(y2)',plot_range(4,ifile))
              call read_var('plot_range(z1)',plot_range(5,ifile))
              call read_var('plot_range(z2)',plot_range(6,ifile))
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
           elseif (index(plot_string,'sph')>0)then    !^CFG IF NOT SIMPLE BEGIN
   	      plot_area='sph'
	      call read_var('R_plot',plot_range(1,ifile))
           elseif (index(plot_string,'los')>0) then
              plot_area='los'
              ! Line of sight vector
              call read_var('los_vector(1)',los_vector(1,ifile))
              call read_var('los_vector(2)',los_vector(2,ifile))
              call read_var('los_vector(3)',los_vector(3,ifile))
              ! read max dimensions of the 2d image plane
              call read_var('x_size_image',x_size_image)
              call read_var('y_size_image',y_size_image)
              ! read the position of image origin relative to grid origin
              call read_var('xoffset',xoffset)
              call read_var('yoffset',yoffset)
              ! read the occulting radius
              call read_var('radius_occult',radius_occult)
              ! read the limb darkening parameter
              call read_var('mu_los',mu_los)
              ! read the number of pixels
              call read_var('n_pix_X',n_pix_X)
              call read_var('n_pix_Y',n_pix_Y)        !^CFG END SIMPLE
           elseif (index(plot_string,'ion')>0) then
              plot_area='ion'
           else
              do i=1,3
                 plot_range(2*i-1,ifile)=XyzMin_D(i)
                 plot_range(2*i  ,ifile)=XyzMax_D(i)
              end do
              if(index(plot_string,'x=0')>0)then
                 plot_area='x=0'
                 select case(TypeGeometry)               !^CFG IF NOT CARTESIAN
                 case('cartesian')                       !^CFG IF NOT CARTESIAN
                    plot_range(1,ifile)=-cTiny*(XyzMax_D(1)-XyzMin_D(1))&
                         /(nCells(1)*proc_dims(1))
                    plot_range(2,ifile)=+cTiny*(XyzMax_D(1)-XyzMin_D(1))&
                         /(nCells(1)*proc_dims(1))
                    !                               ^CFG IF NOT CARTESIAN BEGIN
                 case('spherical','spherical_lnr','cylindrical')  
                    plot_range(3,ifile)=cHalfPi&                       
                         -cTiny*(XyzMax_D(2)-XyzMin_D(2))&
                         /(nCells(2)*proc_dims(2))            
                    plot_range(4,ifile)=cHalfPi&                        
                         +cTiny*(XyzMax_D(2)-XyzMin_D(2))&
                         /(nCells(2)*proc_dims(2)) 
                 case default
                    call stop_mpi(NameSub//' ERROR: unknown geometry type = '&
                         //TypeGeometry)
                 end select
                 !                                  ^CFG END CARTESIAN
              elseif(index(plot_string,'y=0')>0)then
                 plot_area='y=0'
                 plot_range(3,ifile)=-cTiny*(XyzMax_D(2)-XyzMin_D(2))&
                      /(nCells(2)*proc_dims(2)) 
                 plot_range(4,ifile)=+cTiny*(XyzMax_D(2)-XyzMin_D(2))&
                      /(nCells(2)*proc_dims(2)) 
              elseif(index(plot_string,'z=0')>0)then
                 plot_area='z=0'
                 select case(TypeGeometry)               !^CFG IF NOT CARTESIAN
                 case('cartesian','cylindrical')         !^CFG IF NOT CARTESIAN
                    plot_range(5,ifile)=-cTiny*(XyzMax_D(3)-XyzMin_D(3))&
                         /(nCells(3)*proc_dims(3)) 
                    plot_range(6,ifile)=+cTiny*(XyzMax_D(3)-XyzMin_D(3))&
                         /(nCells(3)*proc_dims(3)) 
                    !^CFG IF NOT CARTESIAN BEGIN
                 case('spherical','spherical_lnr')
                    plot_range(5,ifile)=cHalfPi&                        
                         -cTiny*(XyzMax_D(3)-XyzMin_D(3))&
                         /(nCells(3)*proc_dims(3))               
                    plot_range(6,ifile)=cHalfPi&                        
                         +cTiny*(XyzMax_D(3)-XyzMin_D(3))&
                         /(nCells(3)*proc_dims(3)) 
                 case default
                    call stop_mpi(NameSub//' ERROR: unknown geometry type = '&
                         //TypeGeometry)
                 end select
                 !^CFG END CARTESIAN
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
                   .and. plot_area /= 'sph' &        !^CFG IF NOT SIMPLE
                   .and. plot_area /= 'los' &        !^CFG IF NOT SIMPLE
                   .and. plot_area /= 'lin' &        !^CFG IF RAYTRACE
                   ) call read_var('plot_dx',plot_dx(1,ifile))
              !                                     ^CFG IF NOT CARTESIAN BEGIN
              if(TypeGeometry=='spherical'.or.TypeGeometry=='spherical_lnr')&
                   plot_dx(1,ifile)=-1.0 
              !                                     ^CFG END CARTESIAN
           elseif(index(plot_string,'tec')>0)then 
              plot_form(ifile)='tec'
              plot_dx(1,ifile)=0.
           else
              call stop_mpi('Format (idl,tec) missing from plot_string='&
                   //plot_string)
           end if
           if (plot_area == 'sph') then            !^CFG IF NOT SIMPLE BEGIN
              select case(TypeGeometry)                !^CFG IF NOT CARTESIAN
              case('cartesian')                        !^CFG IF NOT CARTESIAN
                 plot_dx(1,ifile) = 1.0    ! set to match value in write_plot_sph
                 plot_dx(2:3,ifile) = 1.0  ! set to degrees desired in angular resolution
                 plot_range(2,ifile)= plot_range(1,ifile) + 1.e-4   ! so that R/=0
                 plot_range(3,ifile)= 0.   - 0.5*plot_dx(2,ifile)
                 plot_range(4,ifile)= 90.0 + 0.5*plot_dx(2,ifile)
                 plot_range(5,ifile)= 0.   - 0.5*plot_dx(3,ifile)
                 plot_range(6,ifile)= 360.0- 0.5*plot_dx(3,ifile)
              case('spherical')                    !^CFG IF NOT CARTESIAN BEGIN
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
                 call stop_mpi(NameSub//' ERROR: unknown geometry type = '&
                      //TypeGeometry)
              end select                           !^CFG END CARTESIAN 
           end if                                  !^CFG END SIMPLE

           ! Plot variables
           if(index(plot_string,'VAR')>0 .or. index(plot_string,'var')>0 )then
              plot_var='var'
              plot_dimensional(ifile) = index(plot_string,'VAR')>0
              call read_var('plot_vars(ifile)',plot_vars(ifile))
              call read_var('plot_pars(ifile)',plot_pars(ifile))
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
              plot_vars(ifile)='rho mx my mz bx by bz e p b1x b1y b1z absdivB'
              plot_pars(ifile)='g c'
           elseif(index(plot_string,'MHD')>0.or.index(plot_string,'mhd')>0)then
              plot_var='mhd'
              plot_dimensional(ifile) = index(plot_string,'MHD')>0
              plot_vars(ifile)='rho ux uy uz bx by bz p jx jy jz'
              plot_pars(ifile)='g c'
           elseif(index(plot_string,'FUL')>0.or.index(plot_string,'ful')>0)then
              plot_var='ful'
              plot_dimensional(ifile) = index(plot_string,'FUL')>0
              plot_vars(ifile)=&
                   'rho ux uy uz bx by bz b1x b1y b1z p e jx jy jz'
              plot_pars(ifile)='g c'
           elseif(index(plot_string,'FLX')>0.or.index(plot_string,'flx')>0)then
              plot_var='flx'
              plot_dimensional(ifile) = index(plot_string,'FLX')>0
              plot_vars(ifile)='rho mr br p jr pvecr'
              plot_pars(ifile)='g c'
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
        call read_var('save_plots_amr',save_plots_amr)
     case("#SAVEBINARY")
        call read_var('save_binary',save_binary)
     case("#AMRLEVELS")
        call read_var('min_block_level',min_block_level)
        call read_var('max_block_level',max_block_level)
        call read_var( 'fix_body_level', fix_body_level)    !^CFG IF NOT SIMPLE
        if(iSession==1)then
           DoSetLevels=.true.
        else
           call set_levels
        end if
     case("#AMRRESOLUTION")
        call read_var('min_cell_dx',min_cell_dx)
        call read_var('max_cell_dx',max_cell_dx)
        call read_var('fix_body_level',fix_body_level)      !^CFG IF NOT SIMPLE
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
        call read_var('dnRefine',dn_refine)
        if (dn_refine > 0)then
           call read_var('DoAutomaticRefinement',automatic_refinement)
           if (automatic_refinement) then
              call read_var('percentCoarsen',percentCoarsen)
              call read_var('percentRefine' ,percentRefine)
              call read_var('MaxTotalBlocks',MaxTotalBlocks)
           end if
        end if
     case("#AMRINIT")
        !                                              ^CFG IF NOT SIMPLE BEGIN
        call read_var('InitialRefineType',InitialRefineType)
        call read_var('InitialRefineLevel',initial_refine_levels)
        !                                              ^CFG END SIMPLE
     case("#AMRCRITERIA")
        call read_var('nRefineCrit',nRefineCrit)
        if(nRefineCrit<1 .or. nRefineCrit>3)&
             call stop_mpi('nRefineCrit must be 1,2,or 3')
        do i=1,nRefineCrit
           call read_var('RefineCrit(i)',RefineCrit(i))
           !                                           ^CFG IF NOT SIMPLE BEGIN
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
           !                                           ^CFG END SIMPLE
        end do
     case("#SCHEME")
        call read_var('nOrder'  ,nOrder)
        !                                              ^CFG IF NOT SIMPLE BEGIN
        call read_var('FluxType',FluxType)
        if(nOrder>1)&                                                
             call read_var('LimiterType',limiter_type)
        if(limiter_type=='beta') call read_var('v_limiter_beta_param',&  
             v_limiter_beta_param)
        !                                              ^CFG END SIMPLE
     case("#NONCONSERVATIVE")
        call read_var('UseNonConservative',UseNonConservative) !^CFG IF NOT SIMPLE
     case("#CONSERVATIVECRITERIA")
        !                                              ^CFG IF NOT SIMPLE BEGIN
        call read_var('nConservCrit',nConservCrit)
        if(nConservCrit > 0) then
           if(allocated(TypeConservCrit_I)) deallocate(TypeConservCrit_I)
           allocate( TypeConservCrit_I(nConservCrit) )
           do i=1,nConservCrit
              call read_var('TypeConservCrit_I',TypeConservCrit_I(i) )
              select case(TypeConservCrit_I(i))
                 !\
                 ! Geometry based criteria: 
                 !    non-conservative scheme is used for r < rConserv
                 !/
              case('R','r','radius','Radius')
                 TypeConservCrit_I(i) = 'r'
                 call read_var('rConserv',rConserv)
                 !\
                 ! Physics based criteria
                 !/
              case('p','P')
                 ! Balsara/Ryu switch 1
                 TypeConservCrit_I(i) = 'p'
                 call read_var('pCoeffConserv',pCoeffConserv)
              case('GradP','gradp')
                 ! Balsara/Ryu switch 2
                 TypeConservCrit_I(i) = 'GradP'
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
        !                                              ^CFG END SIMPLE
     case("#UPDATECHECK")
        !                                              ^CFG IF NOT SIMPLE BEGIN
        call read_var("UseUpdateCheck",UseUpdateCheck)          
        if(UseUpdateCheck)then
           call read_var("rhomin[%]", percent_max_rho(1))        
           call read_var("rhomax[%]", percent_max_rho(2))
           call read_var("Pmin[%]",   percent_max_p(1))
           call read_var("Pmax[%]",   percent_max_p(2))
        end if
        !                                              ^CFG END SIMPLE
     case("#PROLONGATION")
        call read_var('prolong_order',prolong_order)  !^CFG IF NOT SIMPLE
        call read_var('prolong_type' ,prolong_type)   !^CFG IF NOT SIMPLE
     case("#MESSAGEPASS","#OPTIMIZE")
        !                                              ^CFG IF NOT SIMPLE BEGIN
        if(TypeGeometry=='cartesian') then               !^CFG IF NOT CARTESIAN
           call read_var('OptimizeMessagePass'    ,optimize_message_pass)
           if(optimize_message_pass=='allold' .or.&
                optimize_message_pass=='oldopt')then
              if(iProc==0)write(*,'(a)')NameSub// &
                   ' WARNING: message_pass mode='// &
                   trim(optimize_message_pass)// &
                   ' is not available any longer, allopt is set !!!'
              optimize_message_pass='allopt'
           end if
        end if                                           !^CFG IF NOT CARTESIAN
        !                                              ^CFG END SIMPLE
     case("#BORIS")
        !                                              ^CFG IF BORISCORR BEGIN
        call read_var('boris_correction',boris_correction)   
        if(boris_correction) then
           call read_var('boris_cLIGHT_factor',boris_cLIGHT_factor)
           UseBorisSimple=.false.     !^CFG IF SIMPLEBORIS
        else
           boris_cLIGHT_factor = 1.0
        end if
        !                                              ^CFG END BORISCORR
     case("#BORISSIMPLE")                            !^CFG IF SIMPLEBORIS BEGIN
        call read_var('UseBorisSimple',UseBorisSimple)
        if(UseBorisSimple) then
           call read_var('boris_cLIGHT_factor',boris_cLIGHT_factor)
           boris_correction=.false.                     !^CFG IF BORISCORR
        else
           boris_cLIGHT_factor = 1.0
        end if                                       !^CFG END SIMPLEBORIS
     case("#DIVB")
        !                                              ^CFG IF NOT SIMPLE BEGIN
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
        !                                              ^CFG END SIMPLE
     case("#DIVBSOURCE")
	call read_var('UseB0Source'   ,UseB0Source)   !^CFG IF NOT SIMPLE
     case("#PROJECTION")                              !^CFG IF PROJECTION BEGIN
        call read_var('proj_method'   ,proj_method)
        call read_var('proj_typestop ',proj_typestop)
        call read_var('proj_divbcoeff',proj_divbcoeff)
        call read_var('proj_divbconst',proj_divbconst)
        call read_var('proj_matvecmax',proj_matvecmax)
        ! Make sure that DivbMax is recalculated
        DivbMax = -1.0                                !^CFG END PROJECTION
     case("#DIVBDIFFUSION")                           !^CFG IF DIVBDIFFUSE
        call read_var('divb_diffcoeff',divb_diffcoeff)!^CFG IF DIVBDIFFUSE
     case("#CORRECTP")                                !^CFG IF PROJECTION BEGIN
        call read_var('Pratio_lo',Pratio_lo)
        call read_var('Pratio_hi',Pratio_hi)
        if(Pratio_lo>=Pratio_hi)&
             call stop_mpi(NameSub//' ERROR: Pratio_lo>=Pratio_hi')
        !                                              ^CFG END PROJECTION
     case("#SHOCKTUBE")                               !^CFG IF NOT SIMPLE BEGIN
        do i=1,nVar
           call read_var('shock_Lstate',shock_Lstate(i))
        end do
        do i=1,nVar
           call read_var('shock_Rstate',shock_Rstate(i))
        end do
        call read_var('ShockSlope',ShockSlope)
        !                                              ^CFG END SIMPLE
     case("#UPSTREAM_INPUT_FILE")
        call read_var('UseUpstreamInputFile',UseUpstreamInputFile)
        if (UseUpstreamInputFile) then
           read_new_upstream = .true.
           call read_var('UpstreamFileName',UpstreamFileName)
           call read_var('Satellite_Y_Pos',Satellite_Y_Pos)
           call read_var('Satellite_Z_Pos',Satellite_Z_Pos)
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
     case("#IM")                                      !^CFG IF RCM
        call read_var('TauCoupleIm',TauCoupleIm)      !^CFG IF RCM
     case("#MASSLOADING")                             !^CFG IF NOT SIMPLE BEGIN
        call read_var('UseMassLoading',UseMassLoading)
        call read_var('AccelerateMassLoading',AccelerateMassLoading) 
     case("#USER_FLAGS")             !^CFG IF USERFILES BEGIN
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
        call read_var('UseUserSetPhysConst'     ,UseUserSetPhysConst)
        call read_var('UseUserUpdateStates'     ,UseUserUpdateStates)
     case("#USERINPUTBEGIN")        
        call user_read_inputs        !^CFG END USERFILES
        !                                              ^CFG END SIMPLE
     case("#CODEVERSION")
        if(.not.is_first_session())CYCLE
        call read_var('CodeVersion',CodeVersionRead)
        if(abs(CodeVersionRead-CodeVersion)>0.005.and.iProc==0)&
             write(*,'(a,f6.3,a,f6.3,a)')NameSub//&
             ' WARNING: CodeVersion in file=',CodeVersionRead,&
             ' but '//NameThisComp//' version is ',CodeVersion,' !!!'
     case("#PRECISION")
        if(.not.is_first_session())CYCLE
        call read_var('nByteReal',nByteRealRead)
        if(nByteReal/=nByteRealRead.and.iProc==0)then
           write(*,'(a,i1,a)')'BATSRUS was compiled with ',nByteReal,&
                ' byte reals'
           call stop_mpi(NameSub//' ERROR: incorrect precision for reals')
        end if
     case("#EQUATION")
        if(.not.is_first_session())CYCLE
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
        if(.not.is_first_session())CYCLE
        !                                        ^CFG IF NOT SIMPLE BEGIN
        call read_var('problem_type',problem_type_r) 
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
                ' WARNING: probem type set twice !!!'
        end if                                  !^CFG END SIMPLE
     case("#RESTARTINDIR")
        if(.not.is_first_session())CYCLE
        call read_var("NameRestartInDir",NameRestartInDir)
        call fix_dir_name(NameRestartInDir)
        if (iProc==0) call check_dir(NameRestartInDir)
     case("#RESTARTOUTDIR")
        call read_var("NameRestartOutDir",NameRestartOutDir)
        call fix_dir_name(NameRestartOutDir)
        if (iProc==0) call check_dir(NameRestartOutDir)
     case("#PLOTDIR")
        call read_var("NamePlotDir",NamePlotDir)
        call fix_dir_name(NamePlotDir)
        if (iProc==0) call check_dir(NamePlotDir)
     case("#NEWRESTART")
        if(.not.is_first_session())CYCLE
        restart=.true.
        restart_reals=.true.
        restart_ghost=.false.
        call read_var('restart_Bface',restart_Bface) !^CFG IF CONSTRAINB
     case("#BLOCKLEVELSRELOADED")
        if(.not.is_first_session())CYCLE
        ! Sets logical for upgrade of restart files 
        ! to include LEVmin and LEVmax
        RestartBlockLevels=.true.
     case("#SATELLITE")
        if(.not.is_first_session())CYCLE
        call read_var('nsatellite',nsatellite)
        if(nsatellite>0) save_satellite_data=.true.
        if(save_satellite_data)then
           if(iProc==0) call check_dir(NamePlotDir)
           nfile=max(nfile,satellite_+nsatellite)
           if (nfile > maxfile .or. nsatellite > maxsatellitefile)&
                call stop_mpi(&
                'The number of output files is too large in #SATELLITE:'&
                //' nfile>maxfile .or. nsatellite>maxsatellitefile')

           do ifile=satellite_+1,satellite_+nsatellite
              call read_var('satellite_string',satellite_string)


              ! Satellite output frequency
              ! Note that we broke with tradition here so that the
              ! dt_output will always we read!  This may be changed
              ! in later distributions
              call read_var('dn_output',dn_output(ifile))
              call read_var('dt_output',dt_output(ifile))

              ! Satellite outputfile name or the satellite name
              call read_var('Satellite_name',&
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
                 call read_var('satellite_vars',satellite_vars(ifile))
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
           call read_satellite_input_files
        end if
     case("#COVARIANTGEOMETRY")
        !                                         ^CFG IF CARTESIAN BEGIN
        call stop_mpi&
             ('Configure BATSRUS with CARTESIAN=OFF to use the command'&
             //NameCommand)
        !                                         ^CFG END CARTESIAN 
        !^CFG IF NOT CARTESIAN BEGIN
        if(.not.is_first_session())CYCLE
        call read_var('TypeGeometry',TypeGeometry)      
        select case(TypeGeometry)                                   
        case('cartesian') 
           !^CFG IF FACEOUTERBC BEGIN
        case('spherical','spherical_lnr')      
           automatic_refinement=.false.                  
           MaxBoundary = Top_
           if(mod(proc_dims(2),2)==1)proc_dims(2)=2*proc_dims(3)
           DoFixExtraBoundary=.true.
           DoFixOuterBoundary=.true.
        case('cylindrical')
           automatic_refinement=.false.
           MaxBoundary = max(MaxBoundary,North_)
           if(mod(proc_dims(2),2)==1)proc_dims(2)=2*proc_dims(3)
           DoFixExtraBoundary=.true.
           DoFixOuterBoundary=.true.         
           !^CFG END FACEOUTERBC
        case default
           call stop_mpi(NameSub//' ERROR: unknown geometry type = '&
                //TypeGeometry)
        end select
        !^CFG END CARTESIAN 
     case("#GRID")
        !                                        ^CFG IF NOT SIMPLE BEGIN
        if(.not.is_first_session())CYCLE
        call read_var('proc_dims(1)',proc_dims(1)) 
        call read_var('proc_dims(2)',proc_dims(2))
        call read_var('proc_dims(3)',proc_dims(3))
        !^CFG IF NOT CARTESIAN BEGIN
        if( (TypeGeometry=='spherical'.or.TypeGeometry=='cylindrical'&
             .or.TypeGeometry=='spherical_lnr')&
             .and.mod(proc_dims(2),2)==1) then
           proc_dims(2)=2*proc_dims(3)
           if(iProc==0)then
              write(*,*)&
                   'Original number of blocks for the polar angle Phi must be even'
              write(*,*)&
                   'Proc_dim(2) is set to be equal to 2*Proc_dim(3) = ',&
                   proc_dims(2)
           end if
        end if
        !^CFG END CARTESIAN
        if(product(proc_dims)>nBLK.and.iProc==0)then
           write(*,*)'Root blocks will not fit on 1 processor, check nBLK'
           call stop_mpi('product(proc_dims) > nBLK!')
        end if
        call read_var('x1',x1)
        call read_var('x2',x2)
        call read_var('y1',y1)
        call read_var('y2',y2)
        call read_var('z1',z1)
        call read_var('z2',z2)      !^CFG END SIMPLE
        select case(TypeGeometry)   !^CFG IF NOT CARTESIAN
        case('cartesian')           !^CFG IF NOT CARTESIAN
           call set_xyzminmax_cart  
        case('spherical')           !^CFG IF NOT CARTESIAN BEGIN
           call set_xyzminmax_sph  
        case('spherical_lnr')         
           call set_xyzminmax_sph2      
        case('cylindrical')         
           call set_xyzminmax_cyl
        case default
           call stop_mpi(NameSub//' ERROR: unknown geometry type = '&
                //TypeGeometry)
        end select                  !^CFG END CARTESIAN

        !case("#LIMITGENCOORD1")                   !^CFG IF NOT CARTESIAN
        ! call read_var('XyzMin_D(1)',XyzMin_D(1)) !^CFG IF NOT CARTESIAN
        ! call read_var('XyzMax_D(1)',XyzMax_D(1)) !^CFG IF NOT CARTESIAN
     case("#CHECKGRIDSIZE")
        if(.not.is_first_session())CYCLE
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
        if(.not.is_first_session())CYCLE
        call read_var('nRefineLevelIC',nRefineLevelIC)
     case("#GAMMA")
        if(.not.is_first_session())CYCLE
        call read_var('gamma',g)
        !\
        ! Compute gamma related values.
        !/
        gm1     = g-cOne
        gm2     = g-cTwo
        gp1     = g+cOne
        inv_g   = cOne/g
        inv_gm1 = cOne/gm1
        g_half  = cHalf*g
     case('#EXTRABOUNDARY')                   !^CFG IF USERFILES BEGIN
        if(.not.is_first_session())CYCLE
        call read_var('UseExtraBoundary',UseExtraBoundary)
        if(UseExtraBoundary) call read_var('TypeBc_I(ExtraBc_)',&
             TypeBc_I(ExtraBc_))      
        !                                      ^CFG IF FACEOUTERBC BEGIN
        if(TypeGeometry=='cartesian')&             !^CFG IF NOT CARTESIAN
             call read_var('DoFixExtraBoundary',&  
             DoFixExtraBoundary)  
        !                                      ^CFG END FACEOUTERBC
        !                                      ^CFG END USERFILES
     case('#FACEOUTERBC')                      !^CFG IF FACEOUTERBC BEGIN
        if(.not.is_first_session())CYCLE
        call read_var('MaxBoundary',MaxBoundary)
        !^CFG IF NOT CARTESIAN BEGIN
        select case(TypeGeometry)
        case('spherical','spherical_lnr')
           MaxBoundary = Top_
        case('cylindrical')
           MaxBoundary = max(MaxBoundary,North_)
        end select
        !^CFG END CARTESIAN
        if(MaxBoundary>=East_)&
             call read_var('DoFixOuterBoundary',DoFixOuterBoundary) 
        !                                       ^CFG END FACEOUTERBC 
     case("#SOLARWIND")
        if(.not.is_first_session())CYCLE
        call read_var('SW_rho_dim',SW_rho_dim)
        call read_var('SW_T_dim'  ,SW_T_dim)
        call read_var('SW_Ux_dim' ,SW_Ux_dim)
        call read_var('SW_Uy_dim' ,SW_Uy_dim)
        call read_var('SW_Uz_dim' ,SW_Uz_dim)
        call read_var('SW_Bx_dim' ,SW_Bx_dim)
        call read_var('SW_By_dim' ,SW_By_dim)
        call read_var('SW_Bz_dim' ,SW_Bz_dim)
     case("#MAGNETOSPHERE","#BODY")
        if(.not.is_first_session())CYCLE
        call read_var('body1',body1)            !^CFG IF NOT SIMPLE BEGIN
        if(body1)then
           call read_var('Rbody'     ,Rbody)
           call read_var('Rcurrents' ,Rcurrents)
           call read_var('BodyRhoDim',Body_Rho_Dim)
           call read_var('BodyTDim'  ,Body_T_dim)
        end if                                  !^CFG END SIMPLE
     case("#GRAVITY")
        if(.not.is_first_session())CYCLE
        call read_var('UseGravity',UseGravity)        !^CFG IF NOT SIMPLE
        if(UseGravity)call read_var('GravityDir',GravityDir) !^CFG IF NOT SIMPLE
     case("#SECONDBODY")                        !^CFG IF SECONDBODY BEGIN
        if(.not.is_first_session())CYCLE
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
        if(.not.is_first_session())CYCLE
        call read_var('BdpDimBody2x',BdpDimBody2_D(1))
        call read_var('BdpDimBody2y',BdpDimBody2_D(2))
        call read_var('BdpDimBody2z',BdpDimBody2_D(3))
        !                                           ^CFG END SECONDBODY

     case('#PLANET','#IDEALAXES','#ROTATIONAXIS','#MAGNETICAXIS',&
          '#ROTATION','#NONDIPOLE','#UPDATEB0')

        call check_stand_alone
        if(.not.is_first_session())CYCLE

        if(.not.UseNewAxes .and. iProc==0) &
             call stop_mpi(NameSub//': ERROR command '// &
             trim(NameCommand)//' cannot be used with UseNewAxes=F')

        call read_planet_var(NameCommand)

     case("#COROTATION","#AXES")

        call check_stand_alone
        if(.not.is_first_session())CYCLE

        if(UseNewAxes .and. iProc==0) &
             call stop_mpi(NameSub//': ERROR command '// &
             trim(NameCommand)//' cannot be used with UseNewAxes=T')

        call read_compatible_command(NameCommand)

     case("#DIPOLE")

        call check_stand_alone
        if(.not.is_first_session())CYCLE

        if(UseNewAxes)then
           call read_planet_var(NameCommand)
        else
           call read_compatible_command(NameCommand)
        end if

     case("#COORDSYSTEM")
        if(.not.is_first_session())CYCLE
        call read_var('TypeCoordSystem',TypeCoordSystem)
        call upper_case(TypeCoordSystem)
        select case(NameThisComp)
        case('GM')
           if(TypeCoordSystem /= 'GSM')call stop_mpi(NameSub// &
                ' GM_ERROR: cannot handle coordinate system '&
                //TypeCoordSystem)
        case('IH','SC')
           if(TypeCoordSystem /= 'HGI')call stop_mpi(NameSub// &
                NameThisComp//'_ERROR: cannot handle coordinate system '&
                //TypeCoordSystem)
        end select
     case("#NSTEP")
        if(.not.is_first_session())CYCLE
        call read_var('n_step',n_step)
     case("#NPREVIOUS")
        if(.not.is_first_session())CYCLE
        call read_var('nPrev',n_prev)             !^CFG IF IMPLICIT
        call read_var('DtPrev',dt_prev)           !^CFG IF IMPLICIT
     case("#STARTTIME", "#SETREALTIME")
        if(.not.is_first_session())CYCLE
        call read_var('year'  ,iStartTime_I(1))
        call read_var('month' ,iStartTime_I(2))
        call read_var('day'   ,iStartTime_I(3))
        call read_var('hour'  ,iStartTime_I(4))
        call read_var('minute',iStartTime_I(5))
        call read_var('second',iStartTime_I(6))
        iStartTime_I(7) = 0
        if(IsStandAlone)then
           call time_int_to_real(iStartTime_I, StartTime)
        else
           ! Check if things work out or not
           call time_int_to_real(iStartTime_I, StartTimeCheck)
        end if
     case("#TIMESIMULATION")
        if(.not.is_first_session())CYCLE
        if(IsStandAlone) &
             call read_var('tSimulation',time_simulation)
        !                                    ^CFG IF NOT SIMPLE BEGIN
     case("#HELIOSPHERE")
        if(.not.is_first_session())CYCLE
        call read_var('Tsun' ,Body_T_dim)
        call read_var('NDsun',Body_rho_dim)
        call read_var('Qsun' ,Qsun)
        call read_var('Theat',Theat)
        call read_var('Rheat',Rheat)
        call read_var('SIGMAheat',SIGMAheat)
        call read_var('UseFluxRope',UseFluxRope)
        if(UseFluxRope)then
           call read_var('cme_a' ,cme_a)
           call read_var('cme_r1',cme_r1)
           call read_var('cme_r0',cme_r0)
           call read_var('cme_a1',cme_a1)
           call read_var('cme_alpha',cme_alpha)
           call read_var('cme_rho1',cme_rho1)
           call read_var('cme_rho2',cme_rho2)
           call read_var('ModulationRho',ModulationRho)
           call read_var('ModulationP',ModulationP)
           call read_var('cRot_x_GL98'  ,cRot_x_GL98)
           call read_var('cRot_y_GL98'  ,cRot_y_GL98)
           call read_var('cRot_z_GL98'  ,cRot_z_GL98)
        end if
     case("#HELIODIPOLE")
        if(.not.is_first_session())CYCLE
        call read_var('HelioDipoleStrength',Bdp_dim)
        call read_var('HelioDipoleTilt'    ,ThetaTilt)
        ThetaTilt = ThetaTilt * cDegToRad
     case("#HELIOROTATION", "#INERTIAL")
        if(.not.is_first_session())CYCLE
        call read_var('UseInertialFrame',UseInertial)
        UseRotatingFrame = .not.UseInertial
        if(UseInertial)then
           call read_var('UseRotatingBC',UseCorotation)
        else
           ! In the corotating frame the inner BC does not rotate
           UseCorotation=.false.
        end if
     case("#HELIOTEST")
        call read_var('DoSendMHD',DoSendMHD)
     case("#ARCADE")
        if(.not.is_first_session())CYCLE
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
        if(.not.is_first_session())CYCLE
        call read_var('cme_type',cme_type)
        call read_var('cme_a' ,cme_a)
        call read_var('cme_r1',cme_r1)
        call read_var('cme_r0',cme_r0)
        call read_var('cme_a1',cme_a1)
        call read_var('cme_alpha',cme_alpha)
        call read_var('cme_rho1',cme_rho1)
        call read_var('cme_rho2',cme_rho2)
        call read_var('cme_B1_dim',cme_B1_dim)
        call read_var('cme_v_erupt',cme_v_erupt)
     case("#TESTDISSMHD")                         !^CFG IF DISSFLUX BEGIN
        if(.not.is_first_session())CYCLE
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
        if(.not.is_first_session())CYCLE
        call read_var('Qprod' ,Qprod)
        call read_var('Unr_in' ,Unr_in)
        call read_var('mbar',mbar)
        call read_var('ionization_rate',ionization_rate)
        call read_var('kin_in',kin_in)                              
        !                                            ^CFG END SIMPLE
     case default
        if(iProc==0) then
           write(*,*) NameSub // ' WARNING: unknown #COMMAND ' // &
                trim(NameCommand),' !!!'
           if(UseStrict)call stop_mpi('Correct PARAM.in!')
        end if
     end select
  end do

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

    problem_type = -1                                      !^CFG IF NOT SIMPLE

    if(NameThisComp=='IH'.or.NameThisComp=='SC')TypeCoordSystem = 'HGI'

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
    !FluxType = 'Sokolov'              !^CFG UNCOMMENT IF RUSANOVFLUX

    ! Default implicit parameters      !^CFG IF IMPLICIT BEGIN
    UsePointImplicit = .false.              !^CFG IF POINTIMPLICIT
    UsePartImplicit  = .false.
    UseFullImplicit  = .false.
    UseImplicit      = .false.
    ImplCritType     = 'dt'

    nOrder_impl   = 1
    FluxTypeImpl  = 'Rusanov'               !^CFG IF RUSANOVFLUX
    !FluxTypeImpl  = FluxType               !^CFG UNCOMMENT IF RUSANOVFLUX

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

    KrylovType      = 'gmres'
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

    UsePointImplicit = .false.        !^CFG IF POINTIMPLICIT
    UseBorisSimple = .false.          !^CFG IF SIMPLEBORIS
    boris_correction = .false.        !^CFG IF BORISCORR
    boris_cLIGHT_factor = 1.0         

    proc_dims(1) = 1
    proc_dims(2) = 1  
    proc_dims(3) = 1

    x1 = -10.                                         !^CFG IF NOT SIMPLE BEGIN
    x2 =  10.
    y1 = -10.
    y2 =  10.
    z1 = -10.
    z2 =  10.

    call set_xyzminmax_cart                           !^CFG END SIMPLE

    optimize_message_pass = 'allopt'
    UseCorners            = .false.

    UseUpstreamInputFile = .false.

    UseMassLoading        = .false.        !^CFG IF USERFILES
    AccelerateMassLoading = .false.        !^CFG IF USERFILES

    plot_dimensional      = .true.
    save_satellite_data   = .false.

    restart           = .false.
    read_new_upstream = .false.

    if(UseSimple)then                               !^CFG IF SIMPLE BEGIN
       proc_dims(1) = 2
       proc_dims(2) = 1
       proc_dims(3) = 1
       if(nProc<2)then
          write(*,*)'Code is running on ',nProc,' processor'
          call stop_mpi('Not less than 2 processors should be used!!')
       end if



       if (nI == 4 .and. nJ == 4 .and. nK == 4) then 
          x1 = -224.                                           
          x2 =   32.
          y1 =  -64.
          y2 =   64.
          z1 =  -64.
          z2 =   64.
          initial_refine_levels = 7
          InitialRefineType = 'magneto12'
       elseif (nI == 6 .and. nJ == 6 .and. nK == 6) then
          x1 = -288.0
          x2 =   96.0
          y1 =  -96.0
          y2 =   96.0
          z1 =  -96.0
          z2 =   96.0
          initial_refine_levels = 7
          InitialRefineType = 'magneto12'
       elseif (nI == 8 .and. nJ == 8 .and. nK == 8) then
          x1 = -224.0
          x2 =   32.0
          y1 =  -64.0
          y2 =   64.0
          z1 =  -64.0
          z2 =   64.0
          initial_refine_levels = 6
          InitialRefineType = 'magneto12'
       elseif (nI ==10 .and. nJ ==10 .and. nK ==10) then
          x1 = -280.0
          x2 =   40.0
          y1 =  -80.0
          y2 =   80.0
          z1 =  -80.0
          z2 =   80.0
          initial_refine_levels = 6
          InitialRefineType = 'magneto12'
       elseif (nI ==12 .and. nJ ==12 .and. nK ==12) then
          x1 = -288.0
          x2 =   96.0
          y1 =  -96.0
          y2 =   96.0
          z1 =  -96.0
          z2 =   96.0
          initial_refine_levels = 6
          InitialRefineType = 'magneto12'
       else
          write(*,*)'nI, nJ, or nK is not set correctly.  They must be equal '
          write(*,*)'and have a value of 4,6,8,10, or 12 : ',nI,nJ,nK
          call stop_mpi('Please reset nI,nJ,nK in ModSize!!')
       end if

       UseGravity=.false.
       GravityDir=0  
       body1      =.true.
       Rbody      = 2.50
       Rcurrents  = 3.50

       fix_body_level = .true.

       UseNonConservative   = .true.
       nConservCrit         = 1
       if(allocated(TypeConservCrit_I)) deallocate(TypeConservCrit_I)
       allocate( TypeConservCrit_I(nConservCrit) )
       TypeConservCrit_I(1) = 'r'
       rConserv             = 2*rBody

       Bdp_dim              = -31100.0
       THETAtilt       = 0.00
       dt_UpdateB0     = -1.0
       DoUpdateB0      = .false.

       TypeBc_I(east_)        ='outflow'
       TypeBc_I(west_)        ='vary'
       TypeBc_I(south_:top_)  ='fixed'
       TypeBc_I(body1_)       ='ionosphere'

       SW_rho_dim =      5.0    ! n/cc
       SW_T_dim   = 150000.0    ! K
       SW_Ux_dim  =   -400.0    ! km/s
       SW_Uy_dim  =      0.0    ! km/s
       SW_Uz_dim  =      0.0    ! km/s
       SW_Bx_dim  =      0.0    ! nT
       SW_By_dim  =      0.0    ! nT
       SW_Bz_dim  =     -5.0    ! nT

       Body_rho_dim = 10.0    ! n/cc
       Body_T_dim   = 25000.0! K

       nRefineCrit    = 3
       RefineCrit(1)  = 'gradlogP'
       RefineCrit(2)  = 'curlB'
       RefineCrit(3)  = 'Rcurrents'

    end if                                          !^CFG END SIMPLE

  end subroutine set_defaults

  !===========================================================================
  subroutine set_problem_defaults                   !^CFG IF NOT SIMPLE BEGIN


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

    select case(problem_type)                          !^CFG IF CARTESIAN BEGIN
    case(problem_earth,problem_saturn,problem_jupiter,problem_rotation)
       UseNonConservative   = .true.
       nConservCrit         = 1
       if(allocated(TypeConservCrit_I)) deallocate(TypeConservCrit_I)
       allocate( TypeConservCrit_I(nConservCrit) )
       TypeConservCrit_I(1) = 'r'
       rConserv             = 2*rBody
    case default                                       !^CFG END CARTESIAN
       UseNonConservative   = .false.
       if(allocated(TypeConservCrit_I)) deallocate(TypeConservCrit_I)
       nConservCrit         = 0
       rConserv             = -1.
    end select                                         !^CFG IF CARTESIAN

    if(IsStandAlone)then
       select case(problem_type)
       case(problem_earth)
          Bdp_dim = -31100.0
       case(problem_saturn)
          Bdp_dim = 20800.0
       case(problem_jupiter)
          Bdp_dim = 426000.0
       case(problem_comet)
          Bdp_dim = 0.0
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
       TypeBc_I(body1_)='reflect'
    end select

    !   UseRotatingFrame = problem_type == problem_heliosphere
    UseRotatingFrame = .not.UseInertial

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

  end subroutine set_InitialRefineType                      !^CFG END SIMPLE

  !=========================================================================
  subroutine check_options

    ! option and module parameters
    character (len=40) :: Name
    real               :: Version
    logical            :: IsOn


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
    !^CFG IF NOT SIMPLE BEGIN
    if(limiter_type=='mc')then
       call OPTION_FACE(IsOn,Name)
       if(index(Name,'OPTIMIZE')>0)then
          if(iProc==0)then
             write(*,'(a)')NameSub//&
                  ' WARNING: MC limiter is not available any longer !!!'
             if(UseStrict)call stop_mpi( &
                  'Correct PARAM.in or select unoptimized OPTION_FACE!')
             write(*,*)NameSub//' setting beta limiter with parameter 1.2'
          endif
          limiter_type='beta'
          v_limiter_beta_param=1.2
       end if
    end if
    !^CFG END SIMPLE

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
          if(UseStrict)&                             !^CFG IF RUSANOVFLUX
               call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' setting FluxTypeImpl=Rusanov' !^CFG IF RUSANOVFLUX
       end if
       FluxTypeImpl='Rusanov'                        !^CFG IF RUSANOVFLUX
    end select                               !^CFG END IMPLICIT

    if (time_accurate .and. nStage==2) UseBDF2=.true. !^CFG IF IMPLICIT

    ! Make sure periodic boundary conditions are symmetric
    if(any(TypeBc_I(1:2)=='periodic')) TypeBc_I(1:2)='periodic'
    if(any(TypeBc_I(3:4)=='periodic')) TypeBc_I(3:4)='periodic'
    if(any(TypeBc_I(5:6)=='periodic')) TypeBc_I(5:6)='periodic'

    ! Reset initial refine type if it is set to 'default'
    if(InitialRefineType=='default')then
       if(UseSimple)then                             
          InitialRefineType='magneto12'
       else                                          !^CFG IF NOT SIMPLE
          call set_InitialRefineType                 !^CFG IF NOT SIMPLE
       end if
    end if
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

    !^CFG IF DIVBDIFFUSE BEGIN    
    if (UseDivbDiffusion .and. index(optimize_message_pass,'opt') > 0) then
       if(iProc==0 .and. optimize_message_pass /= 'allopt') then
          write(*,'(a)')NameSub//&
               ' WARNING: div B diffusion does not work with'// &
               ' optimize_message_pass='//trim(optimize_message_pass)//' !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' setting optimize_message_pass = all'
       end if
       optimize_message_pass = 'all'
    endif
    !^CFG END DIVBDIFFUSE

    UseCorners= limiter_type=='LSG' .or. limiter_type=='lsg'

    if(prolong_order/=1&
         .and.(index(optimize_message_pass,'old')>0.or.&
         index(optimize_message_pass,'all')>0))&
         call stop_mpi(NameSub// &
         'The prolongation order=2 requires message_pass_dir')

    ! Check test processor
    if(procTEST>nProc)then
       if(iProc==0) write(*,'(a)')NameSub//&
            ' WARNING: procTEST > nProc, setting procTEST=0 !!!'
       procTEST=0
    end if

  end subroutine correct_parameters

  !===========================================================================
  subroutine check_parameters

    ! Make sure the cfl conditions are good for SIMPLE case  !^CFG IF SIMPLE
    if(UseSimple)&                                           !^CFG IF SIMPLE
         cfl=min(cfl,0.80)                                   !^CFG IF SIMPLE

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

    if(UsePointImplicit.and.boris_correction) then     !^CFG IF POINTIMPLICIT BEGIN
       if (iProc == 0)  write(*,'(a)')NameSub// &
            ' WARNING: Point implicit not available for Boris, '// &
            'turning point implicit off !!!'
       if(UseStrict)call stop_mpi('Correct PARAM.in!')
       UsePointImplicit = .false.
    end if                                              !^CFG END POINTIMPLICIT
    if(UseSimple)&                                          !^CFG IF SIMPLE
         boris_cLIGHT_factor=max(boris_cLIGHT_factor,0.02)  !^CFG IF SIMPLE

    ! Finish checks for boris                       !^CFG END BORISCORR

    if(UseSimple)  nStage = max(nStage,nOrder)      !^CFG IF SIMPLE


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
          write(*,'(2a)') NameSub//&
               ' WARNING: PartImplicit scheme does not work with',&
               ' TypeMessagePass=oldopt or allold !!!'
          if(UseStrict)call stop_mpi('Correct PARAM.in!')
          write(*,*)NameSub//' setting optimize_message_pass=allopt'
          optimize_message_pass='allopt'
       end select
    end if
    !Finish checks for implicit                     !^CFG END IMPLICIT

  end subroutine check_parameters

  !===========================================================================
  subroutine set_physics_parameters

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
          dt = DtFixed
          cfl=1.0
       end if
    end if

    ! You have to have the normalizations first
    if (read_new_upstream)  call normalize_upstream_data

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
    !    MaxBoundary=max(MaxBoundary,Top_)           !^CFG IF NOT CELLOUTERBC

    !if(iProc==0)then                                !^CFG IF FACEOUTERBC BEGIN
    !   !These write outs are rather meaningless for the average user
    !   !and the values never change when cell outer bc-s are used 
    !   call write_prefix; write(iUnitOut,*)&
    !        'MinBoundary=',MinBoundary,'MaxBoundary=',MaxBoundary
    !   call write_prefix; write(iUnitOut,*)&
    !        'DoFixOuterBoundary = ',DoFixOuterBoundary
    !   call write_prefix; write(iUnitOut,*)&
    !        'UseExtraBoundary =',UseExtraBoundary,&
    !        'DoFixExtraBoundary =',DoFixExtraBoundary
    !end if                                           !^CFG END FACEOUTERBC 
  end subroutine set_physics_parameters

end subroutine MH_set_parameters
