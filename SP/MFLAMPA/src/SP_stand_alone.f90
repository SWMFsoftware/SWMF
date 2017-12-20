!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program MFLAMPA
  use ModKind
  use SP_ModProc,  ONLY: iProc, nProc, iComm
  use ModUtilities, ONLY: remove_file, touch_file
  use SP_ModMain, ONLY: &
       nTiming, IsLastRead, IsStandAlone, &
       iIterGlobal, TimeGlobal, TimeMax, nIterMax, &
       SP_read_param => read_param, &
       SP_check      => check, &
       SP_initialize => initialize, &
       SP_run        => run!, &
       !SP_finalize   => finalize

  use ModReadParam, ONLY: read_file, read_init
  use ModMpi
  
  implicit none

  integer:: iError
  integer:: iSession = 1  
  real(Real8_) :: CpuTimeStart
  logical:: IsFirstSession = .true.
  real:: Time
  !------------------------------------------------
  !\                                                                        
  ! Initialization of MPI/parallel message passing.                        
  !/ 
  call MPI_INIT(iError)
  iComm=MPI_COMM_WORLD
  call MPI_COMM_RANK(iComm, iProc, iError)
  call MPI_COMM_SIZE(iComm, nProc, iError)

  !\
  ! Initialize time which is used to check CPU time
  !/
  CpuTimeStart = MPI_WTIME()

  !\
  ! Delete MFLAMPA.SUCCESS and MFLAMPA.STOP files if found
  !/
  if(iProc==0)then
     call remove_file('MFLAMPA.SUCCESS')
     call remove_file('MFLAMPA.STOP')
  end if
  
  !\
  ! Mark the run as a stand alone
  !/
  IsStandAlone = .true.

  !\
  ! Read input parameter file. Provide the default restart file for #RESTART
  !/
  call read_file('PARAM.in',iComm)
 
  SESSIONLOOP: do
     call read_init('  ', iSessionIn=iSession)

     if(iProc==0)&
         write(*,*)'----- Starting Session ',iSession,' ------'
     !\
     ! Set and check input parameters for this session
     !/
     call SP_read_param('READ')
     call SP_check
     !\
     ! Time execution (timing parameters were set by MH_set_parameters)
     !/
     if(IsFirstSession)then
        call timing_start('MFLAMPA')
        call timing_start('setup')
     end if
     if(IsFirstSession)then
        call SP_initialize
        Time = TimeGlobal
     end if
     if(IsFirstSession)then
        call timing_stop('setup')
        if(nTiming > -3) call timing_report_total
        if(iProc==0) write(*,*)'Resetting timing counters after setup.'
        call timing_reset('#all',3)
     end if

     TIMELOOP: do
        if(stop_condition_true())exit TIMELOOP
        if(is_time_to_stop())exit SESSIONLOOP
        call timing_step(iIterGlobal + 1)
        
        if(TimeMax > 0.0)then
           call SP_run(Time, TimeMax)
        else
           call SP_run(Time, huge(0.0))
        end if

        call show_progress
     end do TIMELOOP

     if(IsLastRead)exit SESSIONLOOP
     if(iProc==0) &
          write(*,*)'----- End of Session   ',iSession,' ------'   
     iSession=iSession+1
     IsFirstSession = .false.
     if (nTiming > -2) call timing_report
     call timing_reset_all
  end do SESSIONLOOP

  if(iProc==0)then
     write(*,*)
     write(*,'(a)')'    Finished Numerical Simulation'
     write(*,'(a)')'    -----------------------------'
  end if

  if (nTiming > -2) call timing_report

  call timing_stop('MFLAMPA')

  if(nTiming > -3)call timing_report_total

  !Finish writing to log file
  !call SP_finalize

  !\
  ! Touch MFLAMPA.SUCCESS
  !/
  if(iProc==0) call touch_file('MFLAMPA.SUCCESS')

  !Finalize MPI
  call MPI_Finalize(iError)
  
contains
  function stop_condition_true() result(IsStopCondition)
    use SP_ModMain, ONLY: nIterMax
    logical :: IsStopCondition
    !--------------------------
    IsStopCondition = .false.

    if(nIterMax >= 0  .and.iIterGlobal >=nIterMax) IsStopCondition = .true.
    if( TimeMax >  0.0.and. TimeGlobal >= TimeMax) IsStopCondition = .true.

  end function stop_condition_true
  !===============================
  function is_time_to_stop() result(IsTimeToStop)
    use SP_ModMain, ONLY: CpuTimeMax, UseStopFile
    logical :: IsTimeToStop
    !---------------------
    IsTimeToStop = .false.

    if(iProc==0)then
       if(CpuTimeMax > 0.0 .and. MPI_WTIME()-CpuTimeStart >= CpuTimeMax)then
          write(*,*)'CPU time exceeded:',CpuTimeMax,MPI_WTIME()-CpuTimeStart
          IsTimeToStop=.true.
       end if
       if(.not.IsTimeToStop .and. UseStopFile) then
          inquire(file='MFLAMPA.STOP',exist=IsTimeToStop)
          if (IsTimeToStop) &
               write(*,*)'MFLAMPA.STOP file exists: received stop signal'
       end if
    end if
    if(nProc==1) RETURN
    call MPI_BCAST(IsTimeToStop,1,MPI_LOGICAL,0,iComm,iError)

  end function is_time_to_stop
  !===================================================================
  subroutine show_progress
    use SP_ModMain, ONLY: UseTiming
    real(Real8_), external :: timing_func_d
    real(Real8_) :: CpuTimeMFLAMPA, CpuTimeAdvance
    integer:: nProgress1 = 0, nProgress2 = 10
    !------------------------------------------------
    !\
    ! Show timing results if required
    !/
    ! Show speed as cells/second/PE/step
    if( UseTiming .and. iProc==0 &
         .and. nProgress1>0 .and. mod(iIterGlobal,nProgress1) == 0 ) then
       CpuTimeMFLAMPA = timing_func_d('sum',1,'MFLAMPA','MFLAMPA')
       CpuTimeAdvance = timing_func_d('sum',1,'advance','MFLAMPA')
       !\
       ! placeholder
       !/
    end if

    ! Show timing tables
    if(nTiming>0.and.mod(iIterGlobal,nTiming)==0) then
       call timing_report
    else if(nProgress2>0.and.mod(iIterGlobal,nProgress2) == 0) then
       call timing_tree(2,2)
    end if

  end subroutine show_progress

end program MFLAMPA

!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  use ModMpi
  use SP_ModProc
  implicit none
  integer:: iError
  !-------------------
  character (len=*), intent(in) :: StringError
  write(*,'(a)')StringError
  
  !Finalize MPI
  call MPI_Finalize(iError)
  
  stop
end subroutine CON_stop
!============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
end subroutine CON_set_do_test
