!^CMP COPYRIGHT UM
!
!QUOTE:\clearpage
!
!BOP
!
!QUOTE:\section{Execution and Coupling of Components}
!
!MODULE: CON_session - run the framework with a set of parameters
!INTERFACE:
module CON_session

  !USES:
  use CON_world
  use CON_comp_param
  use CON_variables, ONLY: TypeSession, UseStrict
  use CON_wrapper
  use CON_couple_all

  use CON_coupler, ONLY: init_coord_system_all, &
       check_couple_symm, Couple_CC, nCouple, iCompCoupleOrder_II, &
       DoCoupleOnTime_C

  use CON_io, ONLY : DnShowProgressShort, DnShowProgressLong, &
       SaveRestart, save_restart

  use ModUtilities, ONLY: flush_unit

  use CON_time
  use ModFreq
  use ModNumConst
  use ModMpi

  implicit none

  save

  private ! except

  public :: init_session ! initialize SWMF with a fixed set of parameters
  public :: do_session   ! advance    SWMF with a fixed set of parameters

  !REVISION HISTORY:
  ! 08/26/03 G.Toth - initial version
  ! 05/20/04 G.Toth - general steady state session model
  !EOP

  character (len=*), parameter :: NameMod = 'CON_session'

  !\
  ! Local variable definitions.
  !/
  integer :: lComp, iComp, nComp
  integer :: iCouple, iCompSource, iCompTarget
  integer :: iError
  logical :: IsProc_C(MaxComp)      ! for general model

  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------

contains

  !BOP ======================================================================
  !IROUTINE: init_session - initialize run with a fixed set of input parameters
  !INTERFACE:
  subroutine init_session

    !DESCRIPTION:
    ! This subroutine does the session according to the value of TypeSession
    !EOP

    character(len=*), parameter :: NameSub=NameMod//'::init_session'
    !--------------------------------------------------------------------------

    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    !BOC
    nComp = n_comp()
    !\
    ! Initialize components for this session except for IM
    !/
    do lComp = 1,nComp
       iComp = i_comp(lComp)
       call init_session_comp(iComp,iSession,tSimulation)
    end do

    !\
    ! Initialize and broadcast grid descriptors for all the components
    ! This must involve all PE-s.
    !/
    call init_coord_system_all
    do lComp = 1,nComp
       iComp = i_comp(lComp)
       if(use_comp(iComp)) call set_param_comp(iComp,'GRID')
    end do

    !\
    ! Initialize all couplers
    ! This must involve all PE-s.
    !/
    call couple_all_init

    !\
    ! Do the initial coupling according to the value of TypeSession
    !/
    select case(TypeSession)
    case('general')
       call init_session_general
    case('old')                    !^CMP IF GM
       call init_session_old       !^CMP IF GM
    case default
       call CON_stop(NameSub//' SWMF_ERROR invalid TypeSession='//TypeSession)
    end select
    !EOC

  end subroutine init_session

  !BOP ======================================================================
  !IROUTINE: init_session_general - initialize for general layout
  !INTERFACE:
  subroutine init_session_general
    !DESCRIPTION:
    ! Initialize possibly overlapping components for the current session. 
    ! First figure out which components belong to this PE.
    ! Then couple the components in an appropriate order for the first time.
    ! The order is determined by the iCompCoupleOrder\_II array.
    !EOP

    character(len=*), parameter :: NameSub=NameMod//'init_session_general'
    !-----------------------------------------------------------------------
    !\
    ! Figure out which components belong to this PE
    !/
    IsProc_C = .false.
    do lComp = 1,nComp
       iComp = i_comp(lComp)
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
    do iCouple = 1,nCouple

       iCompSource = iCompCoupleOrder_II(1,iCouple)
       iCompTarget = iCompCoupleOrder_II(2,iCouple)

       ! Couple iCompSource --> iCompTarget
       if((IsProc_C(iCompSource) .or. IsProc_C(iCompTarget)) .and. &
            Couple_CC(iCompSource, iCompTarget) % DoThis) &
            call couple_comp(iCompSource, iCompTarget, tSimulation)

    end do

  end subroutine init_session_general

  !^CMP IF GM BEGIN
  !BOP ======================================================================
  !IROUTINE: init_session_old - initialize the backwards compatible way
  !INTERFACE:
  subroutine init_session_old

    !DESCRIPTION:
    ! This subroutine does the initial coupling the old way:
    ! Couple GM and UA to IE\\
    ! solve IE\\
    ! couple IE to GM and UA\\
    ! couple GM and IE to IM\\
    ! couple IM to GM.
    !EOP

    logical :: IsImUninitialized=.true.

    character(len=*), parameter :: NameSub=NameMod//'::init_session_old'
    !--------------------------------------------------------------------------

    call check_couple_symm(IE_,GM_,NameSub)                    !^CMP IF IE
    call check_couple_symm(IM_,GM_,NameSub)                    !^CMP IF IM

    if (use_comp(IH_)) call couple_comp(IH_,GM_,tSimulation)   !^CMP IF IH
  
    !^CMP IF IE BEGIN
    if (use_comp(IE_)) then
       call couple_comp(GM_,IE_,tSimulation)
       if(use_comp(UA_)) call couple_comp(UA_,IE_,tSimulation) !^CMP IF UA
       if (is_proc(IE_) .and. nStep > 0) &
            call run_comp(IE_,tSimulation,tSimulation)
       call couple_comp(IE_,GM_,tSimulation)
       if(use_comp(UA_))call couple_comp(IE_,UA_,tSimulation)  !^CMP IF UA
    end if
    !^CMP END IE

    !^CMP IF IM BEGIN
    if (use_comp(IM_).and.IsImUninitialized) then

       if(DoTestMe) write(*,*)NameSub// &
            ': initializing IM/RCM at start of session.'

       call couple_comp(gm_,im_,tSimulation) ! MH_volume

       call couple_comp(ie_,im_,tSimulation)                 !^CMP IF IE

       call couple_comp(im_,gm_,tSimulation) ! IM_pressure_reset

       IsImUninitialized = .false.
    end if
    !^CMP END IM

  end subroutine init_session_old
  !^CMP END GM

  !BOP ======================================================================
  !IROUTINE: do_session - time loop with a fixed set of input parameters
  !INTERFACE:
  subroutine do_session(IsLastSession)

    !INPUT/OUTPUT ARGUMENTS:
    logical, intent(inout) :: IsLastSession ! set it to true if run should stop

    !DESCRIPTION:
    ! This subroutine does the session according to the value of TypeSession
    !EOP

    character(len=*), parameter :: NameSub=NameMod//'::do_session'
    !--------------------------------------------------------------------------

    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    !BOC
    select case(TypeSession)
    case('general')
       call do_session_general(IsLastSession)
    case('old')                                 !^CMP IF GM
       call do_session_old(IsLastSession)       !^CMP IF GM
    case default
       call CON_stop(NameSub//' SWMF_ERROR invalid TypeSession='//TypeSession)
    end select
    !EOC

  end subroutine do_session

  !BOP ======================================================================
  !IROUTINE: do_session_general - time loop for general layout
  !INTERFACE:
  subroutine do_session_general(IsLastSession)

    !INPUT/OUTPUT ARGUMENTS:
    logical, intent(inout) :: IsLastSession ! set it to true if run should stop

    !DESCRIPTION:
    ! This is the general time looping routine which allows overlap in
    ! the layout of the components but allows concurrent execution. 
    ! Each component has its own tSimulation\_C(iComp). 
    ! The simulation time for control is defined as the {\bf minimum}
    ! of the component times. Always the component which is lagging behind
    ! should run. Coupling occurs when all the involved components reached or
    ! exceeded the next coupling, save restart or check for stop time.
    !EOP

    real :: tSimulation_C(MaxComp)      ! Current time for the component
    real :: tSimulationWait_C(MaxComp)  ! After this time do not progress
    real :: tSimulationLimit_C(MaxComp) ! Time not to be passed by component

    real :: tSimulationCouple, tSimulationLimit

    character(len=*), parameter :: NameSub=NameMod//'::do_session_general'
    !-----------------------------------------------------------------------

    !\
    ! If no component uses this PE and init_session_general did not stop
    ! then simply return
    !/
    if( .not.any(IsProc_C) ) RETURN

    !\
    ! Set initial time for all components to be tSimulation
    !/
    tSimulation_C = tSimulation

    TIMELOOP: do

       !\
       ! Stop this session if stopping conditions are fulfilled
       !/
       if(MaxIteration >= 0 .and. nIteration >= MaxIteration) exit TIMELOOP
       if(DoTimeAccurate .and. tSimulationMax > cZero &
            .and. tSimulation >= tSimulationMax) &
            exit TIMELOOP

       !\
       ! Check periodically for stop file and cpu time
       !/
       if(is_time_to(CheckStop,nStep,tSimulation,DoTimeAccurate))then
          if(do_stop_now())then
             IsLastSession = .true.
             exit TIMELOOP
          end if
       end if

       nStep = nStep+1
       nIteration = nIteration+1
       call timing_step(nStep)

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

             ! Find the time of next save restart, stop check or end of session
             tSimulationLimit = huge(1.0)
             
             if(SaveRestart % DoThis .and. SaveRestart % Dt > 0) &
                  tSimulationLimit = min(tSimulationLimit, SaveRestart % tNext)

             if(CheckStop % DoThis .and. CheckStop % Dt > 0) &
                  tSimulationLimit = min(tSimulationLimit, CheckStop % tNext)

             if(tSimulationMax > 0)&
                  tSimulationLimit = min(tSimulationLimit, tSimulationMax)

             ! The next wait time is the smaller of coupling and limit times
             tSimulationWait_C(iComp) = min(tSimulationCouple,tSimulationLimit)

             if(DoCoupleOnTime_C(iComp))then
                ! Limit component time to the next waiting time
                tSimulationLimit_C(iComp) = tSimulationWait_C(iComp)
             else
                ! Limit time by save restart, stop check and maximum time only
                tSimulationLimit_C(iComp) = tSimulationLimit
             end if

          end do
       end if

       ! Advance solution

       do lComp=1,nComp
          iComp = i_comp(lComp)
          if(.not.IsProc_C(iComp)) CYCLE

          if(DoTimeAccurate)then
             if(tSimulation_C(iComp) < tSimulationWait_C(iComp)) then

                call run_comp(iComp,tSimulation_C(iComp),&
                     tSimulationLimit_C(iComp))

                if(DoTest)write(*,*)NameSub,' run ',NameComp_I(iComp),&
                     ' with tSimulation and Limit=',tSimulation_C(iComp),&
                     tSimulationLimit_C(iComp)
             end if
          else
             if(mod(nStep, DnRun_C(iComp)) == 0) then
                ! tSimulationLimit=Huge since there is no limit on time step
                call run_comp(iComp,tSimulation_C(iComp),Huge(1.0))

                ! There is no progress in time
                tSimulation_C(iComp) = tSimulation

                if(DoTest)write(*,*)NameSub,' run ',NameComp_I(iComp),&
                     ' at nStep=',nStep
             end if
          end if

       end do

       !\
       ! tSimulation for CON is the minimum of the simulation time of the
       ! components present on this processor
       !/
       if(DoTimeAccurate) tSimulation = min(&
            minval(tSimulation_C,     MASK=IsProc_C), &
            minval(tSimulationWait_C, MASK=IsProc_C))

       call show_progress

       !\
       ! Couple components as scheduled
       !/
       do iCouple = 1,nCouple

          iCompSource = iCompCoupleOrder_II(1,iCouple)
          iCompTarget = iCompCoupleOrder_II(2,iCouple)

          ! Couple iCompSource --> iCompTarget
          if( (IsProc_C(iCompSource).or.IsProc_C(iCompTarget)) .and. &
               is_time_to(Couple_CC(iCompSource, iCompTarget),&
               nStep, tSimulation, DoTimeAccurate)) &
               call couple_comp(iCompSource, iCompTarget, tSimulation)

       end do

       !\
       ! Save restart files for the local components when scheduled
       !/
       if( is_time_to(SaveRestart, nStep, tSimulation, DoTimeAccurate) ) then
          do lComp=1,nComp
             iComp = i_comp(lComp)
             if(.not.IsProc_C(iComp)) CYCLE
             call save_restart_comp(iComp, tSimulation)
          end do
       end if

    end do TIMELOOP

  end subroutine do_session_general

  !^CMP IF GM BEGIN
  !BOP ======================================================================
  !IROUTINE: do_session_old - time loop the old way for backwards compatibility
  !INTERFACE:
  subroutine do_session_old(IsLastSession)

    !INPUT/OUTPUT ARGUMENTS:
    logical, intent(inout) :: IsLastSession ! set it to true if run should stop

    !DESCRIPTION:
    ! This is based on BATSRUS main and should give identical results.
    ! GM is the principal component, it controls the flow of time.
    ! IM is coupled in a serial fashion with some time staggering.
    ! IE is also coupled serially, it does the solve between the
    ! GM to IE and IE to GM couplings.
    ! UA is only added experimentally, it was not part of BATSRUS.
    ! IH is only added experimentally, it was not coupled in BATSRUS.
    !EOP

    !^CMP IF IM BEGIN
    real    :: tSimulationIm, tSimulationLimitIm
    real    :: DtCoupleIm
    !^CMP END IM

    real    :: tSimulationUa
    real    :: tSimulationIh
    character(len=*), parameter :: NameSub=NameMod//'::do_session_old'
    !-------------------------------------------------------------------------

    !^CMP IF IM BEGIN
    ! To simplify expressions use this local variable for GM-IM coupling
    ! frequency.
    dtCoupleIM = Couple_CC(IM_,GM_) % Dt
    !^CMP END IM

    ! Simulation time for UA and IH start to be the same as for GM/CON
    tSimulationUa = tSimulation       !^CMP IF UA
    tSimulationIh = tSimulation       !^CMP IF IH

    !\
    ! Begin time step (iterations) loop.
    !/
    TIMELOOP: do

       if(DoTestMe)write(*,*)NameSub,&
            ': Starting TIMELOOP at nStep, tSimulation=',nStep,tSimulation

       !\
       ! Stop this session if stopping conditions are fulfilled
       !/
       if(MaxIteration >= 0 .and. nIteration >= MaxIteration) exit TIMELOOP
       if(DoTimeAccurate .and. tSimulationMax > cZero &
            .and. tSimulation >= tSimulationMax) &
            exit TIMELOOP

       !\
       ! Check periodically for stop file and cpu time
       !/
       if(is_time_to(CheckStop, nStep, tSimulation, DoTimeAccurate))then
          if(do_stop_now())then
             IsLastSession = .true.
             EXIT TIMELOOP
          end if
       end if

       nStep = nStep + 1
       nIteration = nIteration+1

       call timing_step(nStep)

       ! Advance solution

       !^CMP IF UA BEGIN
       ! This is only added experimentally
       if(is_proc(UA_))then
          ! UA can go up to one time step ahead
          do
             if(tSimulationUa > tSimulation) EXIT
             call run_comp(UA_,tSimulationUa, tSimulationMax)
          end do
       end if
       !^CMP END UA

       !^CMP IF IH BEGIN
       ! This is only added experimentally
       if(is_proc(IH_))then
          ! IH can go up to one time step ahead
          !do
          !   if(tSimulationIh > tSimulation) EXIT
          if(tSimulationIh <= tSimulation.or..not.DoTimeAccurate) then
             if(DoTimeAccurate.and.tSimulationMax>cZero)then
                call run_comp(IH_,tSimulationIh,tSimulationMax)
             else
                call run_comp(IH_,tSimulationIh, HUGE(tSimulation))
             end if
          end if
          !end do
       end if
       !^CMP END IH

       ! In the old code GM's time is the time for everyone
       ! Also there was no attempt to make the last time equal to tSimulationMax
       call run_comp(GM_,tSimulation,HUGE(tSimulation))
       call MPI_bcast(tSimulation,1,MPI_REAL,i_proc0(GM_),i_comm(),iError)

       if (is_time_to(SaveRestart, nStep, tSimulation, DoTimeAccurate)) &
            call save_restart

       call show_progress

       ! Couple IH to GM in every step                           !^CMP IF IH
       if (use_comp(IH_)) call couple_comp(ih_,gm_,tSimulation)  !^CMP IF IH

       !^CMP IF IM BEGIN
       !\
       ! Periodically perform IM coupling as required.
       !/
       if (use_comp(IM_)) then
          ! Shift time by dtCoupleIm/2 so that the first coupling
          ! occurs at dtCoupleIm/2.
          if (is_time_to(Couple_CC(IM_,GM_), nStep,&
               tSimulation + 0.5*dtCoupleIM, DoTimeAccurate)) then

             call couple_comp(gm_,im_,tSimulation) ! field line volume

             ! potential from IE to IM               !^CMP IF IE
             call couple_comp(ie_,im_,tSimulation)   !^CMP IF IE

             if(is_proc(IM_))then
                ! Run IM 
                ! from tSimulation-dtCoupleIM/2 to tSimulation+dtCoupleIM/2 + 5
                tSimulationIm      = tSimulation   - 0.5*dtCoupleIM

                tSimulationLimitIm = tSimulationIm + dtCoupleIM

                if(DoTestMe) write(*,*)NameSub,': loop IM_run with',&
                     ' tSimulationIm=',tSimulationIm, &
                     ' tSimulationLimitIm= ',tSimulationLimitIm, &
                     ' dtCoupleIM= ', dtCoupleIM 
                do
                   call run_comp(IM_,tSimulationIm,tSimulationLimitIm)
                   if(tSimulationIm >= tSimulationLimitIm - cTiny) EXIT
                end do
                if(DoTestMe) write(*,*)NameSub,': IM_run finished with',&
                     ' tSimulationIm=',tSimulationIm
                call flush_unit(6)
             end if

             call couple_comp(im_,gm_,tSimulation) ! IM_pressure_reset

          end if
       end if
       !^CMP END IM

       !^CMP IF IE BEGIN
       !\
       ! Periodically perform ionosphere/magnetosphere coupling 
       ! calculations as required.
       !/
       if(is_time_to(Couple_CC(IE_,GM_),nStep,tSimulation,DoTimeAccurate))then
          call couple_comp(gm_,ie_,tSimulation)
          if(use_comp(UA_)) call couple_comp(UA_,IE_,tSimulation) !^CMP IF UA
          if (is_proc(IE_)) call run_comp(IE_,tSimulation,tSimulation)
          call couple_comp(ie_,gm_,tSimulation)
          if(use_comp(UA_))call couple_comp(IE_,UA_,tSimulation)  !^CMP IF UA
       end if
       !^CMP END IE

    end do TIMELOOP

  end subroutine do_session_old
  !^CMP END GM

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
       if(CpuTimeMax > cZero .and. MPI_WTIME()-CpuTimeStart >= CpuTimeMax)&
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

    if( (mod(nStep,DnShowProgressShort)== 0 .or.                 &
         mod(nStep,DnShowProgressLong) == 0 ) .and. is_proc0() ) &
         write(*,'(a,i8,a,g14.6,a,f10.2,a)')                     &
         'Progress:',                                            &
         nStep,' steps,',                                        &
         tSimulation,' s simulation time,',                      &
         MPI_WTIME()-CpuTimeStart,' s CPU time'

    if(  mod(nStep,DnShowProgressLong) == 0 .and. is_proc0() )   &
         call timing_tree(2,2)

  end subroutine show_progress


end module CON_session
