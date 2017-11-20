!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program BATSRUS

  use BATL_lib, ONLY: &
       test_start, test_stop, lVerbose

  use ModBatsrusMethods, ONLY: &
       BATS_init_session, &
       BATS_setup, &
       BATS_advance, &
       BATS_save_files, &
       BATS_finalize

  use ModKind
  use ModProcMH, ONLY: iComm, iProc, nProc
  use ModUtilities, ONLY: remove_file, touch_file
  use ModMain, ONLY: &
       IsStandAlone, &
       time_accurate, time_loop, time_simulation, t_max, &
       n_step, nIter, iteration_number, &
       IsLastRead, &
       dn_timing, UseTiming, UseTimingAll, iUnitTiming
  use ModSetParameters, ONLY: set_parameters
  use ModWriteProgress, ONLY: write_progress
  use ModRestartFile, ONLY: NameRestartInDir
  use CON_planet, ONLY: init_planet_const, set_planet_defaults

  use ModReadParam

  use ModMpi

  implicit none

  integer :: iSession=1
  integer :: iError
  real(Real8_) :: CpuTimeStart

  !----------------------------------------------------------------------------
  !\
  ! Initialization of MPI/parallel message passing.
  !/
  call MPI_INIT(iError)
  iComm = MPI_COMM_WORLD
  call MPI_COMM_RANK (iComm, iProc, iError)
  call MPI_COMM_SIZE (iComm, nProc, iError)

  !\
  ! Initialize some basic variables for the stand alone code
  !/
  IsStandAlone      = .true.

  !\
  ! Initialize the planetary constant library and set Earth
  ! as the default planet.
  !/
  call init_planet_const
  call set_planet_defaults

  !\
  ! Not yet doing the computation
  !/
  time_loop = .false.

  !\
  ! Announce BATSRUS
  !/
  call write_progress(0)

  !\
  ! Initialize time which is used to check CPU time
  !/
  CpuTimeStart = MPI_WTIME()

  !\
  ! Delete BATSRUS.SUCCESS and BATSRUS.STOP files if found
  !/
  if(iProc==0)then
     call remove_file('BATSRUS.SUCCESS')
     call remove_file('BATSRUS.STOP')
  end if

  !\
  ! Read input parameter file. Provide the default restart file for #RESTART
  !/
  call read_file('PARAM.in',iComm,trim(NameRestartInDir)//'restart.H')

  SESSIONLOOP: do
     call read_init('  ', iSessionIn=iSession)

     if(iProc==0.and.lVerbose>=0)&
          write(*,*)'----- Starting Session ',iSession,' ------'
     !\
     ! Set and check input parameters for this session
     !/
     call set_parameters('READ')
     call set_parameters('CHECK')

     !\
     ! Test string is read, so set the test flags now
     !/
     ! call set_oktest('main',DoTest,DoTest)

     !\
     ! Time execution (timing parameters were set by MH_set_parameters)
     !/
     if(iSession==1)then
        call timing_start('BATSRUS')
        call timing_start('setup')
        call BATS_setup
        call BATS_init_session
        call timing_stop('setup')
        if(dn_timing > -3)call timing_report_total
        if(iProc==0)write(*,*)'Resetting timing counters after setup.'
        call timing_reset('#all',3)
     else
        call BATS_init_session
     end if

     TIMELOOP: do
        !\
        ! Stop this session if stopping conditions are fulfilled
        !/
        if (stop_condition_true()) EXIT TIMELOOP
        if(is_time_to_stop())EXIT SESSIONLOOP

        call timing_step(n_step+1)

        if(time_accurate .and. t_max > 0.0) then
           call BATS_advance(t_max)
        else
           call BATS_advance(huge(0.0))
        end if

        call show_progress

     end do TIMELOOP

     if(IsLastRead)then
        EXIT SESSIONLOOP
     else
        if(iProc==0.and.lVerbose>=0) &
             write(*,*)'----- End of Session   ',iSession,' ------'
        iSession=iSession+1
        if (dn_timing > -2) call timing_report
        call timing_reset_all
     end if

  end do SESSIONLOOP
  time_loop = .false.

  if(iProc==0.and.lVerbose>=0)then
     write(*,*)
     write(*,'(a)')'    Finished Numerical Simulation'
     write(*,'(a)')'    -----------------------------'
     if (time_accurate)   write(*, '(a,e13.5,a,f12.6,a,f12.6,a)') &
          '   Simulated Time T = ',time_simulation, &
          ' (',time_simulation/60.00, &
          ' min, ',time_simulation/3600.00,' hrs)'
  end if
  if (dn_timing > -2) call timing_report

  call BATS_save_files('FINALWITHRESTART')

  if(iProc==0.and.lVerbose>0)then
     write(*,*)
     write(*,'(a)')'    Finished Saving Output Files'
     write(*,'(a)')'    ----------------------------'
  end if

  call timing_stop('BATSRUS')

  if(dn_timing > -3)call timing_report_total

  if(UseTimingAll)close(iUnitTiming)

  call BATS_finalize

  !\
  ! Touch BATSRUS.SUCCESS
  !/
  if(iProc==0) call touch_file('BATSRUS.SUCCESS')

  call MPI_finalize(iError)

contains
  !============================================================================

  function stop_condition_true() result(StopConditionTrue)

    logical :: StopConditionTrue

    !--------------------------------------------------------------------------
    StopConditionTrue = .false.

    if(nIter >= 0 .and. iteration_number >= nIter) StopConditionTrue = .true.
    if(time_accurate .and. t_max > 0.0 &
         .and. time_simulation >= t_max) &
         StopConditionTrue = .true.

  end function stop_condition_true
  !============================================================================

  function is_time_to_stop() result(IsTimeToStop)

    use ModMain, ONLY: cputime_max, check_stopfile

    logical :: IsTimeToStop

    !--------------------------------------------------------------------------
    IsTimeToStop = .false.

    if(iProc==0)then
       if(cputime_max > 0.0 .and. MPI_WTIME()-CpuTimeStart >= cputime_max)then
          write(*,*)'CPU time exceeded:',cputime_max,MPI_WTIME()-CpuTimeStart
          IsTimeToStop=.true.
       end if
       if(.not.IsTimeToStop .and. check_stopfile) then
          inquire(file='BATSRUS.STOP',exist=IsTimeToStop)
          if (IsTimeToStop) &
               write(*,*)'BATSRUS.STOP file exists: recieved stop signal'
       end if
    end if

    call MPI_BCAST(IsTimeToStop,1,MPI_LOGICAL,0,iComm,iError)

  end function is_time_to_stop
  !============================================================================

  subroutine show_progress

    use ModIo,      ONLY: dn_progress1, dn_progress2
    use ModMain,    ONLY: nI, nJ, nK, nBlock, Unused_B, n_step, Dt
    use ModPhysics, ONLY: Si2No_V, UnitT_

    real(Real8_), external :: timing_func_d
    real(Real8_) :: CpuTimeBATSRUS,CpuTimeAdvance

    integer :: nIterExpect, nIterExpectTime

    !\
    ! Show timing results if required
    !/

    ! Show speed as cells/second/PE/step

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'show_progress'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if( UseTiming .and. iProc==0 &
         .and. dn_progress1>0 .and. mod(n_step,dn_progress1) == 0 ) then

       CpuTimeBATSRUS=timing_func_d('sum',3,'BATSRUS','BATSRUS')
       CpuTimeAdvance=timing_func_d('sum',1,'advance','BATSRUS')
       if (.not.time_accurate) then
          write(*,'(a,f9.1,a,f9.1,a,i8)') 'Speed is',&
               nI*nJ*nK*count(.not.Unused_B(1:nBlock)) &
               /max(1.D-10,CpuTimeAdvance),&
               ' c/s/p after',&
               CpuTimeBATSRUS,&
               ' s at N =',n_step
       else
          write(*,'(a,f9.1,a,f9.1,a,i8,a,1p,e10.4,a)') 'Speed is',&
               nI*nJ*nK*count(.not.Unused_B(1:nBlock)) &
               /max(1.D-10,CpuTimeAdvance),&
               ' c/s/p after',&
               CpuTimeBATSRUS,&
               ' s at N =',n_step, ' (', Time_Simulation,' s)'
       endif

    endif

    ! Show timing tables
    if(dn_timing>0.and.mod(iteration_number,dn_timing)==0) then
       call timing_report
    else if(dn_progress2>0.and.mod(n_step,dn_progress2) == 0) then
       call timing_tree(2,2)
    end if

    ! Try to estimate the remaining length of the run
    if(UseTiming .and. iProc==0 &
         .and. dn_progress2>0 .and. mod(n_step,dn_progress2) == 0)then
       nIterExpect = nITER-iteration_number
       if(time_accurate .and. Dt > 0.0 .and. t_max > 0.0)then
          nIterExpectTime = min( real(huge(1)), &
               (t_max - time_simulation)/Dt*Si2No_V(UnitT_) )
          if(nIterExpect < 0)then
             nIterExpect = nIterExpectTime
          else if(nIterExpectTime > 0)then
             nIterExpect = min(nIterExpect, nIterExpectTime)
          endif
       end if
       CpuTimeAdvance=timing_func_d('sum/iter',2,'advance','BATSRUS')
       write (*,'(i6,a,i6,a,f10.2)') iteration_number,' of ', &
            nIterExpect+iteration_number,  &
            ' iterations completed.   Expected time to completion:', &
            CpuTimeAdvance*nIterExpect
       write(*,*)
    end if

    call test_stop(NameSub, DoTest)
  end subroutine show_progress
  !============================================================================

end program BATSRUS
!==============================================================================

! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
subroutine CON_stop(StringError)
  implicit none
  character (len=*), intent(in) :: StringError
  !----------------------------------------------------------------------------
  call stop_mpi(StringError)
end subroutine CON_stop
!==============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe

  !----------------------------------------------------------------------------
end subroutine CON_set_do_test
!==============================================================================

! The subroutines and functions below are defined in srcInterface for SWMF.
subroutine update_grad_density
  !----------------------------------------------------------------------------
  call stop_mpi('ERROR: update_grad_density is for SWMF')
end subroutine update_grad_density
!==============================================================================
subroutine get_from_spher_buffer_grid(Xyz_D,nVar,State_V)
  implicit none
  real,dimension(3),intent(in)::Xyz_D
  integer,intent(in)::nVar
  real,dimension(nVar)::State_V
  !----------------------------------------------------------------------------
  call stop_mpi( &
       'ERROR: get_from_spher_buffer_grid is for SWMF')
end subroutine get_from_spher_buffer_grid
!==============================================================================
subroutine plot_buffer(iFile)
  implicit none
  integer, intent(in)::iFile
  !----------------------------------------------------------------------------
  call stop_mpi( &
       'ERROR: plot_buffer is for SWMF')
end subroutine plot_buffer
!==============================================================================
subroutine read_ih_buffer(y,z,State_V)
  real :: y, z, State_V(8)
  !----------------------------------------------------------------------------
  call stop_mpi('ERROR: read_ih_buffer is for SWMF')
end subroutine read_ih_buffer
!==============================================================================
subroutine read_pw_buffer(FaceCoords_D,nVar,FaceState_V)
  real, intent(in) :: FaceCoords_D(3)
  integer, intent(in) :: nVar
  real, intent(inout) :: FaceState_V(nVar)
  !----------------------------------------------------------------------------
  call stop_mpi('ERROR: read_pw_buffer is for SWMF')
end subroutine read_pw_buffer
!==============================================================================
