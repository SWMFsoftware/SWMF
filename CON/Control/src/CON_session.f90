!^CMP COPYRIGHT UM
!
!BOP
!
!MODULE: CON_session - run the framework with a set of parameters
!INTERFACE:
module CON_session

  !USES:
  use CON_comp_param, ONLY: MaxComp, NameComp_I
  use CON_world, ONLY: i_comm, is_proc, is_proc0, i_proc, &
       i_comp, n_comp, use_comp
  use CON_variables, ONLY: UseStrict, DnTiming, lVerbose, iErrorSwmf
  use CON_wrapper, ONLY: set_param_comp, init_session_comp, run_comp
  use CON_couple_all, ONLY: couple_two_comp, couple_all_init
  use CON_coupler, ONLY: &
       check_couple_symm, Couple_CC, nCouple, iCompCoupleOrder_II, &
       DoCoupleOnTime_C
  use CON_io, ONLY : DnShowProgressShort, DnShowProgressLong, &
       SaveRestart, save_restart
  use CON_time, ONLY: iSession, DoTimeAccurate, &
       nStep, nIteration, MaxIteration, DnRun_C, tSimulation, tSimulationMax, &
       CheckStop, DoCheckStopFile, CpuTimeSetup, CpuTimeStart, CpuTimeMax
  use ModFreq, ONLY: is_time_to, FreqType
  use ModMpi, ONLY: MPI_WTIME, MPI_LOGICAL

  implicit none

  save

  private ! except

  public :: init_session ! initialize SWMF with a fixed set of parameters
  public :: do_session   ! advance    SWMF with a fixed set of parameters

  !REVISION HISTORY:
  ! 08/26/03 G.Toth - initial version
  ! 05/20/04 G.Toth - general steady state session model
  ! 08/11/04 G.Toth - removed 'old' and 'parallel' session models
  ! 02/09/05 G.Toth - moved related code from CON_main into init_session.
  ! 09/09/05 G.Toth - removed read_inputs from init_session 
  !                   (to avoid ifort compiler bug)
  !                   moved some code from CON_main into do_session.
  ! 01/20/06 G.Toth - added optional tCoupleExtra_C parameter to do_session
  !EOP

  character (len=*), parameter :: NameMod = 'CON_session'

  !\
  ! Local variable definitions.
  !/
  integer :: lComp, iComp, nComp
  integer :: iCouple, iCompSource, iCompTarget
  integer :: iError
  logical :: IsProc_C(MaxComp)

  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------

contains

  !BOP ======================================================================
  !IROUTINE: init_session - initialize run with a fixed set of input parameters
  !INTERFACE:
  subroutine init_session

    !DESCRIPTION:
    ! Initialize possibly overlapping components for the current session. 
    ! First figure out which components belong to this PE.
    ! Then couple the components in an appropriate order for the first time.
    ! The order is determined by the iCompCoupleOrder\_II array.
    ! Do timings as needed.
    !EOP

    character(len=*), parameter :: NameSub=NameMod//'::init_session'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    DoTestMe = DoTest .and. is_proc0()
    !\
    ! Time execution (timing parameters were read by read_inputs)
    !/
    if(iSession==1)then
       call timing_start('SWMF')
       call timing_start('SETUP')
    end if

    if(is_proc0().and.lVerbose>=0)&
         write(*,*)'----- Starting Session ',iSession,' ------'

    !BOC
    nComp = n_comp()
    !\
    ! Initialize components for this session.
    !/
    do lComp = 1, nComp; iComp = i_comp(lComp)
       call init_session_comp(iComp,iSession,tSimulation)
    end do
    !\
    ! Initialize and broadcast grid descriptors to all the components.
    ! This must involve all PE-s.
    !/
    do lComp = 1, nComp; iComp = i_comp(lComp)
       if(use_comp(iComp)) call set_param_comp(iComp,'GRID')
    end do
    !\
    ! Initialize all couplers. This must involve all PE-s.
    !/
    call couple_all_init
    !\
    ! Figure out which components belong to this PE
    !/
    IsProc_C = .false.
    do lComp = 1, nComp; iComp = i_comp(lComp)
       if(.not.use_comp(iComp)) CYCLE
       IsProc_C(iComp) = is_proc(iComp)
    end do
    !\
    ! Check for unused PE
    !/
    if(.not.any(IsProc_C))then
       write(*,*)NameSub//' WARNING: no component uses iProc=',i_proc()
       if(UseStrict)call CON_stop(NameSub// &
            'SWMF_ERROR: unused PE. Edit LAYOUT.in!')
       RETURN
    end if
    !\
    ! Couple for the first time
    !/
    do iCouple = 1, nCouple

       iCompSource = iCompCoupleOrder_II(1,iCouple)
       iCompTarget = iCompCoupleOrder_II(2,iCouple)

       ! Couple iCompSource --> iCompTarget
       if((IsProc_C(iCompSource) .or. IsProc_C(iCompTarget)) .and. &
            Couple_CC(iCompSource, iCompTarget) % DoThis) &
            call couple_two_comp(iCompSource, iCompTarget, tSimulation)
    end do
    !EOC

    if(iSession==1)then
       call timing_stop('SETUP')
       call timing_stop('SWMF')
       CpuTimeSetup = MPI_WTIME()
       if(DnTiming > -3)call timing_report_total
       if(is_proc0())write(*,*)'Resetting timing counters after setup.'
       call timing_reset('#all',3)
       call timing_start('SWMF')
    end if

  end subroutine init_session

  !BOP ======================================================================
  !IROUTINE: do_session - time loop with a fixed set of input parameters
  !INTERFACE:
  subroutine do_session(IsLastSession, tCoupleExtra_C)

    !INPUT/OUTPUT ARGUMENTS:
    logical, intent(inout) :: IsLastSession ! set it to true if run should stop

    !OPTIONAL INPUT ARGUMENTS:
    real, optional, intent(in) :: tCoupleExtra_C(MaxComp) ! external coupling

    !DESCRIPTION:
    ! This subroutine executes one session.
    ! This is general time looping routine allows overlap in
    ! the layout of the components but allows concurrent execution 
    ! if there is no overlap.
    ! Each component has its own tSimulation\_C(iComp). 
    ! The simulation time for control is defined as the {\bf minimum}
    ! of the component times. Always the component which is lagging behind
    ! should run. Coupling occurs when all the involved components reached or
    ! exceeded the next coupling, save restart or check for stop time.

    !LOCAL VARIABLES:
    real :: tSimulation_C(MaxComp)=-1.0 ! Current time for the component
    real :: tSimulationWait_C(MaxComp)  ! After this time do not progress
    real :: tSimulationLimit_C(MaxComp) ! Time not to be passed by component

    real :: tSimulationCouple, tSimulationLimit

    !EOP

    character(len=*), parameter :: NameSub=NameMod//'::do_session'
    !--------------------------------------------------------------------------

    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(DoTestMe)then
       write(*,*)NameSub,' IsLastSession=',IsLastSession
       if(present(tCoupleExtra_C)) &
            write(*,*)NameSub,'tCoupleExtra_C=',tCoupleExtra_C(1:nComp)
    end if
    
    !BOC
    !\
    ! If no component uses this PE and init_session did not stop
    ! then set tSimulation to the final time and simply return
    !/
    if( .not.any(IsProc_C) ) then
       if(DoTimeAccurate .and. tSimulationMax > 0.0) &
            tSimulation = tSimulationMax
       RETURN
    end if

    !\
    ! Set initial time for all components to be tSimulation
    ! For steady state mode do this only for the first time,
    ! because the SWMF time does not advance but the component time might.
    !/
    if(DoTimeAccurate .and. .not.present(tCoupleExtra_C))then
       tSimulation_C = tSimulation
    else
       where(tSimulation_C < 0.0 .and. IsProc_C) tSimulation_C = tSimulation
    end if

    if(DoTestMe)write(*,*)NameSub,' tSimulation_C=',tSimulation_C(1:nComp)

    TIMELOOP: do

       if(DoTestMe)write(*,*)NameSub,' nIteration, tSimulationMax=',&
            nIteration, tSimulationMax
       !\
       ! Stop this session if stopping conditions are fulfilled
       !/
       if(MaxIteration >= 0 .and. nIteration >= MaxIteration) exit TIMELOOP
       if(DoTimeAccurate .and. tSimulationMax > 0.0 &
            .and. tSimulation >= tSimulationMax) &
            exit TIMELOOP

       !\
       ! Exit from time loop and return if an external coupling should be done
       ! Return is used so that iSession is not modified.
       !/
       if(present(tCoupleExtra_C))then
          if(DoTestMe)write(*,*)NameSub,' checking tCoupleExtra_C'
          if(any(IsProc_C .and. tCoupleExtra_C >= 0.0 &
               .and. tSimulation >= tCoupleExtra_C)) RETURN
       end if

       !\
       ! Check periodically for stop file and cpu time
       !/
       if(is_time_to(CheckStop,nStep,tSimulation,DoTimeAccurate))then
          if(DoTestMe)write(*,*)NameSub,' checking do_stop_now'
          if(do_stop_now())then
             IsLastSession = .true.
             exit TIMELOOP
          end if
       end if

       !\
       ! Advance step number and iteration (in current run)
       ! Inform the TIMING utility about the new time step
       !/
       nStep = nStep+1
       nIteration = nIteration+1
       call timing_step(nStep)

       if(DoTestMe)write(*,*)NameSub,' new nStep, nIteration=',&
            nStep, nIteration

       !\
       ! Calculate next time to synchronize for all local components
       !/
       if(DoTimeAccurate)then
          do lComp=1,nComp
             iComp = i_comp(lComp)
             if(.not.IsProc_C(iComp)) CYCLE

             ! Find the time of the next coupling
             tSimulationCouple = min( &
                  minval(Couple_CC(iComp,:) % tNext, &
                  MASK  =Couple_CC(iComp,:) % DoThis), &
                  minval(Couple_CC(:,iComp) % tNext, &
                  MASK  =Couple_CC(:,iComp) % DoThis))
             
             ! Check for external coupling if present
             if(present(tCoupleExtra_C))then
                if(tCoupleExtra_C(iComp) > 0.0) &
                     tSimulationCouple = &
                     min(tSimulationCouple, tCoupleExtra_C(iComp))
             end if

             if(DoTestMe)write(*,*)NameSub,': iComp, tSimulationCouple=', &
                  iComp, tSimulationCouple

             ! Find the time of next save restart, stop check or end of session
             tSimulationLimit = huge(1.0)
             
             if(SaveRestart % DoThis .and. SaveRestart % Dt > 0) &
                  tSimulationLimit = min(tSimulationLimit, SaveRestart % tNext)

             if(DoTestMe)write(*,*)NameSub,': restart tSimulationLimit=', &
                  tSimulationLimit

             if(CheckStop % DoThis .and. CheckStop % Dt > 0) &
                  tSimulationLimit = min(tSimulationLimit, CheckStop % tNext)

             if(DoTestMe)write(*,*)NameSub,': checkstop tSimulationLimit=', &
                  tSimulationLimit

             if(tSimulationMax > 0)&
                  tSimulationLimit = min(tSimulationLimit, tSimulationMax)

             if(DoTestMe)write(*,*)NameSub,': tmax tSimulationLimit=', &
                  tSimulationLimit

             ! The next wait time is the smaller of coupling and limit times
             tSimulationWait_C(iComp) = min(tSimulationCouple,tSimulationLimit)

             if(DoTestMe)write(*,*)NameSub,': restart tSimulationWait_C=', &
                  tSimulationWait_C(iComp)

             if(DoCoupleOnTime_C(iComp))then
                ! Limit component time to the next waiting time
                tSimulationLimit_C(iComp) = tSimulationWait_C(iComp)
             else
                ! Limit time by save restart, stop check and maximum time only
                tSimulationLimit_C(iComp) = tSimulationLimit
             end if

             if(DoTestMe)write(*,*)NameSub,': final tSimulationLimit_C=', &
                  tSimulationLimit_C(iComp)

          end do
       end if

       if(DoTestMe)write(*,*)NameSub,' advance solution'

       ! Advance solution

       do lComp=1,nComp
          iComp = i_comp(lComp)
          if(.not.IsProc_C(iComp)) CYCLE

          if(DoTimeAccurate)then
             !\
             ! In time accurate mode advance the component(s) which
             ! has/have the smallest physical time == tSimulation
             !/
             if(tSimulation_C(iComp) <= tSimulation) then

                call run_comp(iComp,tSimulation_C(iComp),&
                     tSimulationLimit_C(iComp))

                if(DoTest)write(*,*)NameSub,' run ',NameComp_I(iComp),&
                     ' with tSimulation and Limit=',tSimulation_C(iComp),&
                     tSimulationLimit_C(iComp)
             end if
          else
             !\
             ! In steady state mode advance component every DnRun step
             !/
             if(mod(nStep, DnRun_C(iComp)) == 0) then
                ! tSimulationLimit=Huge since there is no limit on time step
                call run_comp(iComp,tSimulation_C(iComp),Huge(1.0))

                ! There is no progress in time
                !tSimulation_C(iComp) = tSimulation

                if(DoTest)write(*,*)NameSub,' run ',NameComp_I(iComp),&
                     ' at nStep=',nStep
             end if
          end if

       end do

       !\
       ! tSimulation for CON is the minimum of the simulation times of the
       ! components present on this processor
       !/
       if(DoTimeAccurate) tSimulation = min(&
            minval(tSimulation_C,     MASK=IsProc_C), &
            minval(tSimulationWait_C, MASK=IsProc_C))

       !\
       ! Print progress report at given frequency
       !/
       call show_progress

       if(DoTestMe)write(*,*)NameSub,' couple components'
       !\
       ! Couple components as scheduled
       !/
       do iCouple = 1, nCouple

          iCompSource = iCompCoupleOrder_II(1,iCouple)
          iCompTarget = iCompCoupleOrder_II(2,iCouple)

          ! Couple iCompSource --> iCompTarget
          if( (IsProc_C(iCompSource).or.IsProc_C(iCompTarget)) .and. &
               is_time_to(Couple_CC(iCompSource, iCompTarget),&
               nStep, tSimulation, DoTimeAccurate))then
             if(DoTestMe)write(*,*)NameSub,' coupling ',iCompSource,iCompTarget,tSimulation
             call couple_two_comp(iCompSource, iCompTarget, tSimulation)
          end if

       end do

       !\
       ! Save restart files when scheduled
       !/
       if( is_time_to(SaveRestart, nStep, tSimulation, DoTimeAccurate) ) &
            call save_restart

    end do TIMELOOP

    if(.not.IsLastSession)then
       if(is_proc0().and.lVerbose>=0) &
            write(*,*)'----- End of Session   ',iSession,' ------'
       iSession=iSession+1
       if (DnTiming > -2) call timing_report
       call timing_reset_all
    end if

    !EOC

    if(DoTestMe)write(*,*)NameSub,' finished'

  end subroutine do_session

  !BOP =======================================================================
  !IROUTINE: do_stop_now - return true if stop file exists or cpu time is over
  !INTERFACE:
  function do_stop_now() result(DoStopNow)

    !RETURN VALUE:
    logical :: DoStopNow

    !DESCRIPTION:
    ! Check if SWMF.STOP file is present in the run directory.
    ! Also check if cpu time limit is exceeded.
    ! If any of these conditions hold return true, otherwise false
    ! on all PE-s. This subroutine contains an MPI\_bcast, so frequent
    ! checks may affect parallel performance.
    !EOP

    DoStopNow=.false.
    if(is_proc0())then
       if(CpuTimeMax > 0.0 .and. MPI_WTIME()-CpuTimeStart >= CpuTimeMax)&
            then
          write(*,*)'CPU time exceeded:',CpuTimeMax,&
               MPI_WTIME()-CpuTimeStart
          DoStopNow=.true.
       end if
       if(.not.DoStopNow .and. DoCheckStopfile) then
          inquire(file='SWMF.STOP',exist=DoStopNow)
          if (DoStopNow) &
               write(*,*)'SWMF.STOP file exists: recieved stop signal'
       end if
    end if

    call MPI_bcast(DoStopNow,1,MPI_LOGICAL,0,i_comm(),iError)

  end function do_stop_now

  !BOP =======================================================================
  !IROUTINE: show_progress - show how far the simulation progressed
  !INTERFACE:
  subroutine show_progress
    !DESCRIPTION:
    ! Print information about time and time step at the periodicity
    ! given by DnShowProgressShort and DnShowProgressLong which can
    ! be set with the \#PROGRESS command in the input parameter file.
    !EOP

    if(.not.is_proc0()) RETURN

    if(is_proc0() .and. ( &
         (DnShowProgressShort>0 .and. mod(nStep,DnShowProgressShort)==0) .or. &
         (DnShowProgressLong >0 .and.  mod(nStep,DnShowProgressLong)==0))) then
       write(*,'(a,i8,a,g14.6,a,f10.2,a)')          &
            'Progress:',                            &
            nStep,' steps,',                        &
            tSimulation,' s simulation time,',      &
            MPI_WTIME()-CpuTimeSetup,' s CPU time'
    end if

    if( DnTiming > 0 .and. mod(nStep,DnTiming) == 0) then
       call timing_report
    elseif(is_proc0() .and. &
         DnShowProgressLong>0 .and. mod(nStep,DnShowProgressLong)== 0) then
       call timing_tree(2,2)
    end if

  end subroutine show_progress

end module CON_session
