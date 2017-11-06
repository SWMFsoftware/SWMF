!******************************************************************************
!
!                               sami3_main.f90
!
!  Contact: Alex Glocer   at alex.glocer-1@nasa.gov,   301-286-9475.
!  Contact: Joe Huba at  huba@nrl.navy.mil
!
!******************************************************************************

program sami3
  use ModSAMI,    ONLY: iProc,nProc,iComm
  use ModMpi
!  use CON_planet, ONLY: init_planet_const, set_planet_defaults
!  use ModPrerunField, ONLY: UsePrerun, read_prerun, read_prerun_IE
 
  implicit none
 
  integer :: iError, iStep
  real    :: DtAdvance
  real    :: DtRestart = 300.0 ! currently set at 300s, should be read in later
  real    :: DtMax = 5.0      ! maximum timestep
  real    :: Time = 0.0
  real    :: TimeMax=10.0
  !---------------------------------------------------------------------------

  !****************************************************************************
  ! Initiallize MPI and get number of processors and rank of given processor
  !****************************************************************************
  
  !---------------------------------------------------------------------------
  call MPI_INIT(iError)
  iComm = MPI_COMM_WORLD
  
  call MPI_COMM_RANK(iComm,iProc,iError)
  call MPI_COMM_SIZE(iComm,nProc,iError)
  
  !****************************************************************************
  ! Read the input file
  !****************************************************************************
  !IsStandAlone=.true.

  ! Initial setup for the SAMI3 model
!  call read_file('PARAM.in',iComm)
!  call read_init('  ',iSessionIn=1,iLineIn=0)
!  call CIMI_set_parameters('READ')

  !if (usePrerun) call read_prerun(t)
  !if (usePrerun .and. iConvect==2) call read_prerun_IE(t)

  !\
  ! Initialize the planetary constant library and set Earth
  ! as the default planet.
  !/
 ! call init_planet_const
 ! call set_planet_defaults
  !****************************************************************************
  ! Initialize the model
  !****************************************************************************
  !call init_mod_cimi
  !call init_mod_field_trace
  

  ! Start Timing
  call timing_active(.true.)
  call timing_step(0)
  call timing_start('SAMI3')

  ! init model
  call timing_start('sami_init')
  call sami_init
  call timing_stop('sami_init')
  
  if (iProc == 0)call timing_report_total
  call timing_reset('#all',3)

  !****************************************************************************
  ! start Timestepping
  !****************************************************************************
  iStep = 0
  TIMELOOP:do
     !report progress on proc 0
     if (iProc==0)write(*,*) 'In Time Loop iStep,Time = ', iStep,Time
     ! If Time exceeds max time then stop advancing
     if (Time >= TimeMax) exit TIMELOOP
     
     ! Set time to advance the model to either the time to reach TimeMax or 
     ! the next restart time
     DtAdvance = min(TimeMax - Time, DtMax)
  
     ! If DtAdvance is too small then just stop advancing
!!!     if (DtAdvance < 1.0e-6) then
!!!        Time = TimeMax
!!!        exit TIMELOOP
!!!     endif
     
     write(*,*) 'calling sami_run at time', Time,'iProc=',iProc,'DtAdvance=',DtAdvance
    
     ! Call sami_run to advance the Timestep
     call timing_step(iStep)
     call timing_start('sami_run')     
     call sami_run(DtAdvance)
     call timing_stop('sami_run')
     
     Time = Time+DtAdvance

!     ! Save restart at DtSaveRestart or TimeMax
!     if (floor((Time+1.0e-5)/DtSaveRestart) /= &
!          floor((Time+1.0e-5-DtAdvance)/DtSaveRestart)) then
!        call sami_write_restart
!     endif
     
     ! Advance the time iteration step
     iStep=iStep+1
  end do TIMELOOP

  ! Save restart at TimeMax
!  call sami_write_restart

  ! finalize SAMI
  call sami_finalize

  ! Finalize timing commands
  call timing_stop('SAMI3')

  if (iProc == 0) then
     write(*,'(a)') 'Finished SAMI run, (reporting timings)'
     write(*,'(a)') '--------------------------------------'
     call timing_report
  endif
  
  ! Finalize MPI
  call MPI_FINALIZE(iError)

end program sami3


!============================================================================
subroutine CON_stop(StringError)
  use ModSAMI,    ONLY: iProc,nProc,iComm!, Time
  use ModMpi
  implicit none
  character (len=*), intent(in) :: StringError
  
  ! Local variables:
  integer :: iError,nError
  !----------------------------------------------------------------------------
  
  write(*,*)'Stopping execution! me=',iProc,' at time=',&!!!,Time,&
       ' with msg:'
  write(*,*)StringError
  call MPI_abort(iComm, nError, iError)
  stop
end subroutine CON_stop
!============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
  DoTest   = .false.
  DoTestMe = .false.
end subroutine CON_set_do_test
