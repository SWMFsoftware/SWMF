!******************************************************************************
!
!                               cimi_main.f90
!  Created on 2 Sept 2014 by Alex Glocer and Mei-Ching Fok, Code 673, NASA GSFC.
!
!  Contact: Alex Glocer   at alex.glocer-1@nasa.gov,   301-286-9475.
!  Contact: Mei-Ching Fok at mei-ching.h.fok@nasa.gov, 301-286-1083.
!
!******************************************************************************

program cimi
  use ModCimiGrid,    ONLY: iProc,nProc,iComm
  use ModCIMI,        ONLY: IsStandalone
  use ModMpi
  use ModCimi,        ONLY: init_mod_cimi, Time
  use ModFieldTrace,  ONLY: init_mod_field_trace
  use ModImTime,      ONLY: TimeMax
  use ModCimiRestart, ONLY: DtSaveRestart,cimi_write_restart
  use ModReadParam
  use ModPrerunField, ONLY: UsePrerun, read_prerun, read_prerun_IE, DtRead
  use ModImSat,       ONLY: UsePrerunSat, read_prerun_sat, DtReadSat, &
       IsFirstWrite
  use ModIeCimi,      ONLY: UseWeimer
  use CON_planet, ONLY: init_planet_const, set_planet_defaults
!  use ModPrerunField, ONLY: UsePrerun, read_prerun, read_prerun_IE
 
  implicit none
 
  integer :: iError, iStep
  real    :: DtAdvance
  real    :: DtRestart = 300.0 ! currently set at 300s, should be read in later
  real    :: DtMax = 60.0      ! maximum timestep
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
  IsStandAlone=.true.

  ! Initial setup for the rbe model
  call read_file('PARAM.in',iComm)
  call read_init('  ',iSessionIn=1,iLineIn=0)
  call CIMI_set_parameters('READ')

  !if (usePrerun) call read_prerun(t)
  !if (usePrerun .and. iConvect==2) call read_prerun_IE(t)

  !\
  ! Initialize the planetary constant library and set Earth
  ! as the default planet.
  !/
  call init_planet_const
  call set_planet_defaults
  !****************************************************************************
  ! Initialize the model
  !****************************************************************************
  call init_mod_cimi
  call init_mod_field_trace
  
  ! Start Timing
  call timing_active(.true.)
  call timing_step(0)
  call timing_start('CIMI')

  !read initial prerun field
  if (UsePrerun) call read_prerun(Time)
  if (UsePrerun .and. .not.UseWeimer) call read_prerun_IE(Time)

  ! init model
  call timing_start('cimi_init')
  call cimi_init
  call timing_stop('cimi_init')
  
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
     if (DtAdvance < 1.0e-6) then
        Time = TimeMax
        exit TIMELOOP
     endif
     
     ! Read next prerun sat file if using prerun sat files and it is
     ! time to read.     
     if (floor((Time+1.0e-5)/DtReadSat) /= &
          floor((Time+1.0e-5+DtAdvance)/DtReadSat)) then
        if (UsePrerunSat) call read_prerun_sat(Time+DtAdvance)
     endif

     ! Call cimi_run to advance the Timestep
     call timing_step(iStep)
     call timing_start('cimi_run')     
     call cimi_run(DtAdvance)
     call timing_stop('cimi_run')
     
     ! Save restart at DtSaveRestart or TimeMax
     if (floor((Time+1.0e-5)/DtSaveRestart) /= &
          floor((Time+1.0e-5-DtAdvance)/DtSaveRestart)) then
        call cimi_write_restart
     endif
     
     ! Read new prerun file if using prerun fields and it is time to read
     if (floor((Time+1.0e-5)/DtRead) /= &
          floor((Time+1.0e-5-DtAdvance)/DtRead)) then 
        if (UsePrerun) call read_prerun(Time)
        if (UsePrerun .and. .not.UseWeimer) call read_prerun_IE(Time)     
     endif

     ! Advance the time iteration step
     iStep=iStep+1
  end do TIMELOOP

  ! Save restart at TimeMax
  call cimi_write_restart

  ! Finalize timing commands
  call timing_stop('CIMI')

  if (iProc == 0) then
     write(*,'(a)') 'Finished CIMI run, (reporting timings)'
     write(*,'(a)') '--------------------------------------'
     call timing_report
  endif
  
  ! Finalize MPI
  call MPI_FINALIZE(iError)

end program cimi
!============================================================================
subroutine CON_stop(StringError)
  use ModCimiGrid,    ONLY: iProc,nProc,iComm
  use ModCimi,        ONLY: Time
  use ModMpi
  implicit none
  character (len=*), intent(in) :: StringError
  
  ! Local variables:
  integer :: iError,nError
  !----------------------------------------------------------------------------
  
  write(*,*)'Stopping execution! me=',iProc,' at time=',Time,&
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
