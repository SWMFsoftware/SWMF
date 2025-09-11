!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module CON_session

  ! Execute one session based on the PARAM.in file where the sessions
  ! are separated by #RUN commands.

  use CON_comp_param, ONLY: MaxComp, NameComp_I, PC_, PT_
  use CON_bline, ONLY: BL_
  use CON_world, ONLY: i_comm, is_proc, is_proc0, i_proc, &
       i_comp, n_comp, use_comp, is_thread, world_used, CON_

  use CON_variables, ONLY: UseStrict, DnTiming, lVerbose
  use CON_wrapper, ONLY: set_param_comp, init_session_comp, run_comp
  use CON_couple_all, ONLY: couple_two_comp, couple_all_init
  use CON_coupler, ONLY: &
       Couple_CC, nCouple, iCompCoupleOrder_II, &
       DoCoupleOnTime_C, IsTightCouple_CC
  use CON_io, ONLY: DnShowProgressShort, DnShowProgressLong, &
       SaveRestart, save_restart, IsRestartSaved
  use CON_time, ONLY: iSession, DoTimeAccurate, &
       nStep, nIteration, MaxIteration, DnRun_C, tSimulation, tSimulationMax, &
       CheckStop, DoCheckStopFile, CpuTimeSetup, CpuTimeStart, CpuTimeMax, &
       IsForcedStop, NameCompCheckKill, DoCheckTimeStep, DnCheckTimeStep, &
       nIterationCheck, tSimulationCheck, TimeStepMin, UseEndTime
  use ModFreq, ONLY: is_time_to
  use ModUtilities, ONLY: CON_stop, CON_set_do_test
  use ModMpi, ONLY: MPI_WTIME, MPI_LOGICAL
  use CON_transfer_data, ONLY: transfer_real

  implicit none

  save

  private ! except

  public :: init_session ! initialize SWMF with a fixed set of parameters
  public :: do_session   ! advance    SWMF with a fixed set of parameters

  ! revision history:
  ! 08/26/03 G.Toth - initial version
  ! 05/20/04 G.Toth - general steady state session model
  ! 08/11/04 G.Toth - removed 'old' and 'parallel' session models
  ! 02/09/05 G.Toth - moved related code from CON_main into init_session.
  ! 09/09/05 G.Toth - removed read_inputs from init_session
  !                   (to avoid ifort compiler bug)
  !                   moved some code from CON_main into do_session.
  ! 01/20/06 G.Toth - added optional tCoupleExtra_C parameter to do_session

  character (len=*), parameter :: NameMod = 'CON_session'

  ! Local variables
  integer :: lComp, iComp, nComp
  integer :: iCouple, iCompSource, iCompTarget
  integer :: iError
  logical :: IsProc_C(MaxComp) ! is this PE actively used by a component?

  logical :: DoTest, DoTestMe

contains
  !============================================================================
  subroutine init_session

    ! Initialize possibly overlapping components for the current session.
    ! First figure out which components belong to this PE.
    ! Then couple the components in an appropriate order for the first time.
    ! The order is determined by the iCompCoupleOrder\_II array.
    ! Do timings as needed.

    logical:: UseCore

    character(len=*), parameter:: NameSub = 'init_session'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    DoTestMe = DoTest .and. is_proc0()

    ! Time execution (timing parameters were read by read_inputs)
    if(iSession==1)then
       call timing_start('SWMF')
       call timing_start('SETUP')
    end if

    if(is_proc0().and.lVerbose>=0)&
         write(*,*)'----- Starting Session ',iSession,' ------'

    nComp = n_comp()

    ! Initialize components for this session.
    do lComp = 1, nComp; iComp = i_comp(lComp)
       call init_session_comp(iComp,iSession,tSimulation)
    end do

    ! Initialize and broadcast grid descriptors to all the components.
    ! This must involve all PE-s.
    do lComp = 1, nComp; iComp = i_comp(lComp)
       if(use_comp(iComp)) call set_param_comp(iComp, 'GRID')
    end do

    ! Initialize all couplers. This must involve all PE-s.
    call couple_all_init

    ! Figure out which components belong to this PE
    IsProc_C = .false.
    UseCore  = .false.

    do lComp = 1, nComp; iComp = i_comp(lComp)
       if(.not.use_comp(iComp)) CYCLE
       IsProc_C(iComp) = is_proc(iComp)
       UseCore = UseCore .or. is_thread(iComp)
    end do

    ! Create communicator of active PEs
    call world_used(IsVerbose=.true.)

    if(.not.UseCore)then
       write(*,*)NameSub//' WARNING: no component uses iProc=',i_proc()
       if(UseStrict)call CON_stop(NameSub// &
            'SWMF_ERROR: unused core. Correct layout or number of threads!')
    end if

    ! Processors that are inactive (may be used by OpenMP threads) are done
    if(.not.is_proc(CON_)) RETURN

    ! Couple for the first time
    do iCouple = 1, nCouple

       iCompSource = iCompCoupleOrder_II(1,iCouple)
       iCompTarget = iCompCoupleOrder_II(2,iCouple)

       ! Couple iCompSource --> iCompTarget
       if((IsProc_C(iCompSource) .or. IsProc_C(iCompTarget)) .and. &
            Couple_CC(iCompSource, iCompTarget) % DoThis) &
            call couple_two_comp(iCompSource, iCompTarget, tSimulation)
    end do

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
  !============================================================================
  subroutine do_session(IsLastSession, tCoupleExtra_C)

    logical, intent(inout) :: IsLastSession ! set it to true if run should stop

    real, optional, intent(in) :: tCoupleExtra_C(MaxComp) ! external coupling

    ! This subroutine executes one session.
    ! This is general time looping routine allows overlap in
    ! the layout of the components but allows concurrent execution
    ! if there is no overlap.
    ! Each component has its own tSimulation\_C(iComp).
    ! The simulation time for control is defined as the {\bf minimum}
    ! of the component times. Always the component which is lagging behind
    ! should run. Coupling occurs when all the involved components reached or
    ! exceeded the next coupling, save restart or check for stop time.

    ! local variables

    real :: tSimulation_C(MaxComp)=-1.0 ! Current time for the component
    real :: tSimulationWait_C(MaxComp)  ! After this time do not progress
    real :: tSimulationLimit_C(MaxComp) ! Time not to be passed by component

    real :: tSimulationCouple, tSimulationLimit

    ! Indexes for tight coupling
    integer:: iCompMaster, iCompSlave

    ! Time info transferred from master to slave
    real:: tMaster

    ! Check if a component just ran
    logical:: DoneRun_C(MaxComp)

    ! Check if the run should be killed
    logical:: DoKill

    ! Smallest temporal frequency and a small fraction
    real:: DtTiny = 0.0

    ! If true, save restart before coupling.
    logical:: DoSaveRestartBeforeCoupling

    logical:: DoSaveRestart

    character(len=*), parameter:: NameSub = 'do_session'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(DoTestMe)then
       write(*,*)NameSub,' IsLastSession, tSimulationMax=', &
            IsLastSession, tSimulationMax
       if(present(tCoupleExtra_C)) &
            write(*,*)NameSub,'tCoupleExtra_C=',tCoupleExtra_C(1:nComp)
    end if

    ! If no component uses this PE and init_session did not stop
    ! then set tSimulation to the final time and simply return
    if( .not.is_proc(CON_) ) then
       if(DoTimeAccurate .and. tSimulationMax > 0.0) &
            tSimulation = tSimulationMax

       RETURN
    end if

    ! Set initial time for all components to be tSimulation
    ! For steady state mode do this only for the first time,
    ! because the SWMF time does not advance but the component time might.
    if(DoTimeAccurate .and. .not.present(tCoupleExtra_C))then
       tSimulation_C = tSimulation
    else
       where(tSimulation_C < 0.0 .and. IsProc_C) tSimulation_C = tSimulation
    end if

    if(DoTimeAccurate)then
       ! Calculate the smallest temporal frequency and a tiny fraction of it
       DtTiny = huge(1.0)

       if(DoTestMe)write(*,*)NameSub,' initial DtTiny=', DtTiny

       ! Restart frequency
       if(SaveRestart % DoThis .and. SaveRestart % Dt > 0) &
            DtTiny = SaveRestart % Dt

       if(DoTestMe)write(*,*)NameSub,' restart DtTiny=', DtTiny

       ! Stop check frequency
       if(CheckStop % DoThis .and. CheckStop % Dt > 0) &
            DtTiny = min(DtTiny, CheckStop % Dt)

       if(DoTestMe)write(*,*)NameSub,' check stop DtTiny=', DtTiny

       ! Active coupling frequencies
       DtTiny = min( DtTiny, &
            minval(Couple_CC % Dt, &
            MASK=(Couple_CC % DoThis .and. Couple_CC % Dt > 0.0)))

       if(DoTestMe)write(*,*)NameSub,' couple DtTiny=', DtTiny

       if(DtTiny == huge(1.0)) DtTiny = 0.0

       ! Tiny fraction
       DtTiny = DtTiny * 1e-6

       if(DoTestMe)write(*,*)NameSub,' final DtTiny=', DtTiny

    end if

    if(DoTestMe)then
       write(*,*)NameSub,' nIteration, nStep=', nIteration, nStep
       write(*,*)NameSub,' tSimulation_C=',tSimulation_C(1:nComp)
    end if

    DoSaveRestartBeforeCoupling = .false.
    ! Saving restart before coupling if PC_ or PT_ is involved
    if(use_comp(PC_) .or. (use_comp(PT_).and.PT_/=BL_)) then
       DoSaveRestartBeforeCoupling = .true.
    end if

    TIMELOOP: do

       if(DoTestMe)write(*,*)NameSub,&
            ' nIteration, tSimulation, tSimulationMax=',&
            nIteration, tSimulation, tSimulationMax

       ! Stop this session if stopping conditions are fulfilled
       if(MaxIteration >= 0 .and. nIteration >= MaxIteration) &
            EXIT TIMELOOP
       if(DoTimeAccurate .and. tSimulation + DtTiny >= tSimulationMax) &
            EXIT TIMELOOP

       ! Exit from time loop and return if an external coupling should be done
       ! Return is used so that iSession is not modified.
       if(present(tCoupleExtra_C))then
          if(DoTestMe)write(*,*)NameSub,' checking tCoupleExtra_C'
          if(any(IsProc_C .and. tCoupleExtra_C >= 0.0 &
               .and. tSimulation >= tCoupleExtra_C)) RETURN
       end if

       ! Check periodically for stop file and cpu time
       if(is_time_to(CheckStop, nStep, tSimulation+DtTiny, DoTimeAccurate))then
          if(DoTestMe)write(*,*)NameSub,' checking do_stop_now'
          if(do_stop_now())then
             IsForcedStop  = .true.
             IsLastSession = .true.
             EXIT TIMELOOP
          end if
       end if

       ! Check for kill file unless NameCompCheckKill is set to '!!'
       if(NameCompCheckKill /= '!!')then
          ! The NameCompCheckKill='??' means that all processors check
          DoKill = NameCompCheckKill == '??'
          ! Otherwise only the root of the NameCompCheckKill component checks
          if(.not.DoKill) DoKill = is_proc0(NameCompCheckKill)
          if(DoKill)then
             ! Check if SWMF.KILL file exists
             inquire(file='SWMF.KILL', exist=DoKill)
             if(DoKill) call CON_stop(NameSub//': SWMF.KILL file found')
          end if
       end if

       ! Check if there is sufficient progress
       if(DoTimeAccurate .and. DoCheckTimeStep &
            .and. nIteration == nIterationCheck + DnCheckTimeStep)then
          if(tSimulation - tSimulationCheck < TimeStepMin*DnCheckTimeStep) &
               call CON_stop(NameSub//': advance too slow in time '// &
               'at nIteration, tSimulation, tSimulationCheck, TimeStepMin=', &
               nIteration, tSimulation, tSimulationCheck, TimeStepMin)
          nIterationCheck  = nIteration
          tSimulationCheck = tSimulation
       end if
       ! For a new step, the restart files have not been saved.
       IsRestartSaved = .false.

       ! Advance step number and iteration (in current run)
       ! Inform the TIMING utility about the new time step
       nStep = nStep+1
       nIteration = nIteration+1
       call timing_step(nStep)

       if(DoTestMe)write(*,*)NameSub,' new nStep, nIteration=',&
            nStep, nIteration

       ! Calculate next time to synchronize for all local components
       if(DoTimeAccurate)then
          do lComp = 1, nComp
             iComp = i_comp(lComp)
             if(.not.use_comp(iComp)) CYCLE

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

          ! Loop through potential tight couplings and ensure that
          ! they use the same time step limit
          do iCouple = 1, nCouple

             iCompMaster = iCompCoupleOrder_II(1,iCouple)
             iCompSlave  = iCompCoupleOrder_II(2,iCouple)

             ! Check if there is a master-slave relationship
             if(.not. IsTightCouple_CC(iCompMaster,iCompSlave)) CYCLE

             ! Exclude other processors
             if(.not.(IsProc_C(iCompMaster) .or. IsProc_C(iCompSlave))) CYCLE

             ! Make sure that slaves and master have the same limiting
             ! This loop may need to be iterated if the slave can have
             ! extra couplings
             tMaster = min( &
                  tSimulationLimit_C(iCompMaster), &
                  tSimulationLimit_C(iCompSlave))
             tSimulationLimit_C(iCompMaster) = tMaster
             tSimulationLimit_C(iCompSlave)  = tMaster
          end do

       end if

       if(DoTestMe)write(*,*)NameSub,' advance solution'

       DoneRun_C = .false.

       ! Advance solution
       do lComp = 1, nComp
          iComp = i_comp(lComp)
          if(.not.IsProc_C(iComp)) CYCLE

          if(DoTimeAccurate)then

             ! In time accurate mode advance the component(s) which
             ! has/have the smallest physical time == tSimulation
             if(tSimulation_C(iComp) <= tSimulation) then

                call run_comp(iComp, &
                     tSimulation_C(iComp), tSimulationLimit_C(iComp))

                DoneRun_C(iComp) = .true.

                if(DoTestMe)write(*,*)NameSub,' run ',NameComp_I(iComp),&
                     ' with tSimulation and Limit=',tSimulation_C(iComp),&
                     tSimulationLimit_C(iComp)
             end if
          else
             ! In steady state mode advance component every DnRun step
             if(mod(nStep, DnRun_C(iComp)) == 0) then
                ! tSimulationLimit=Huge since there is no limit on time step
                call run_comp(iComp,tSimulation_C(iComp),Huge(1.0))

                DoneRun_C(iComp) = .true.

                ! There is no progress in time
                ! tSimulation_C(iComp) = tSimulation

                if(DoTestMe)write(*,*)NameSub,' run ',NameComp_I(iComp),&
                     ' at nStep=',nStep
             end if
          end if

       end do

       if(DoTimeAccurate)then
          ! Loop through potential tight couplings
          do iCouple = 1, nCouple

             iCompMaster = iCompCoupleOrder_II(1,iCouple)
             iCompSlave  = iCompCoupleOrder_II(2,iCouple)

             ! Check if there is a master-slave relationship
             if(.not. IsTightCouple_CC(iCompMaster,iCompSlave)) CYCLE

             ! Exclude other processors
             if(.not.(IsProc_C(iCompMaster) .or. IsProc_C(iCompSlave))) CYCLE

             ! Check if the master-slave pair did indeed run
             ! The simulation times for master and slave must be identical
             if(IsProc_C(iCompMaster) .and. .not. DoneRun_C(iCompMaster)) CYCLE
             if(IsProc_C(iCompSlave)  .and. .not. DoneRun_C(iCompSlave)) CYCLE

             ! Force slave simulation time to agree with master
             if(IsProc_C(iCompMaster)) tMaster = tSimulation_C(iCompMaster)
             call transfer_real(iCompMaster, iCompSlave, tMaster, &
                  UseSourceRootOnly = .false.)
             if(IsProc_C(iCompSlave)) tSimulation_C(iCompSlave) = tMaster

             ! Update simulation time on these processors
             tSimulation = min(&
                  minval(tSimulation_C,     MASK=IsProc_C), &
                  minval(tSimulationWait_C, MASK=IsProc_C))

             ! Couple Master to Slave
             call couple_two_comp(iCompMaster, iCompSlave, tSimulation)

             ! Couple Slave to Master if 2-way coupling is requested
             if(Couple_CC(iCompSlave,iCompMaster) % DoThis) &
                  call couple_two_comp(iCompSlave, iCompMaster, tSimulation)

             if(DoTestMe)write(*,*)NameSub, &
                  ' coupling master, slave, time, two-way=', &
                  iCompMaster, iCompSlave, tSimulation, &
                  Couple_CC(iCompSlave,iCompMaster) % DoThis
          end do
       end if

       ! tSimulation for CON is the minimum of the simulation times of the
       ! components present on this processor
       if(DoTimeAccurate) tSimulation = min(&
            minval(tSimulation_C,     MASK=IsProc_C), &
            minval(tSimulationWait_C, MASK=IsProc_C))

       ! Print progress report at given frequency
       call show_progress

       ! Save restart files when scheduled except for the final save
       ! with UseEndTime which overwrites tSimulation.
       DoSaveRestart = .false.
       if(is_time_to(SaveRestart, nStep, tSimulation+DtTiny, DoTimeAccurate) &
            .and. &
            .not. (UseEndTime .and. tSimulation + DtTiny >= tSimulationMax)) &
            then
          DoSaveRestart = .true.
       end if

       if(DoSaveRestartBeforeCoupling .and. DoSaveRestart) then
          call save_restart
          IsRestartSaved = .true.
       end if

       if(DoTestMe)write(*,*)NameSub,' couple components'

       ! Couple components as scheduled
       do iCouple = 1, nCouple

          iCompSource = iCompCoupleOrder_II(1,iCouple)
          iCompTarget = iCompCoupleOrder_II(2,iCouple)

          ! Exclude other processors
          if(.not.(IsProc_C(iCompSource) .or. IsProc_C(iCompTarget))) CYCLE

          if(is_time_to(Couple_CC(iCompSource, iCompTarget),&
               nStep, tSimulation+DtTiny, DoTimeAccurate))then

             if(DoTestMe)write(*,*)NameSub, &
                  ' coupling ',iCompSource,iCompTarget,tSimulation

             call couple_two_comp(iCompSource, iCompTarget, tSimulation)
          end if

       end do

       if(.not.DoSaveRestartBeforeCoupling .and. DoSaveRestart) then
          call save_restart
          IsRestartSaved = .true.
       end if

    end do TIMELOOP

    if(.not.IsLastSession)then
       if(is_proc0().and.lVerbose>=0) &
            write(*,*)'----- End of Session   ',iSession,' ------'
       iSession=iSession+1
       if (DnTiming > -2) call timing_report
       call timing_reset_all
    end if

    if(DoTestMe)write(*,*)NameSub,' finished'

  end subroutine do_session
  !============================================================================
  function do_stop_now() result(DoStopNow)

    !RETURN VALUE:
    logical :: DoStopNow

    ! Check if SWMF.STOP file is present in the run directory.
    ! Also check if cpu time limit is exceeded.
    ! If any of these conditions hold return true, otherwise false
    ! on all PE-s. This subroutine contains an MPI\_bcast, so frequent
    ! checks may affect parallel performance.
    !--------------------------------------------------------------------------
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
  !============================================================================
  subroutine show_progress

    ! Print information about time and time step at the periodicity
    ! given by DnShowProgressShort and DnShowProgressLong which can
    ! be set with the \#PROGRESS command in the input parameter file.

    use CON_time, ONLY: get_time, TimeType

    type(TimeType):: TimeCurrent
    !--------------------------------------------------------------------------
    if(.not.is_proc0()) RETURN

    if(is_proc0() .and. ( &
         (DnShowProgressShort>0 .and. mod(nStep,DnShowProgressShort)==0) .or. &
         (DnShowProgressLong >0 .and.  mod(nStep,DnShowProgressLong)==0))) then

       call get_time(TimeCurrentOut = TimeCurrent)
       write(*,'(a,i8,a,es14.6,a,f10.2,a)')         &
            'Progress:',                            &
            nStep,' steps,',                        &
            tSimulation,' s simulation time,',      &
            MPI_WTIME()-CpuTimeSetup,' s CPU time'//&
            ', Date: ' &
            //TimeCurrent % String(1:8) // '_' // TimeCurrent % String(9:14)
    end if

    if( DnTiming > 0 .and. mod(nStep,DnTiming) == 0) then
       call timing_report
    elseif(is_proc0() .and. &
         DnShowProgressLong>0 .and. mod(nStep,DnShowProgressLong)== 0) then
       call timing_tree(2,2)
    end if

  end subroutine show_progress
  !============================================================================
end module CON_session
!==============================================================================
