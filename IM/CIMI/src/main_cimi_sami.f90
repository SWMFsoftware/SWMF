!******************************************************************************
!
!                               main_cimi_sami.f90
!  Created on September 14, 2015 by Alex Glocer, Code 673, NASA GSFC.
!
!  Contact: Alex Glocer   at alex.glocer-1@nasa.gov,   301-286-9475.
!  Contact: Joe Huba      at huba@nrl.navy.mil
!  Contact: Mei-Ching Fok at mei-ching.h.fok@nasa.gov, 301-286-1083.
!
!******************************************************************************

program cimi_sami
  use ModSAMI,        ONLY: iProcSAMI=>iProc,nProcSAMI=>nProc,iCommSAMI=>iComm
  use ModCimiGrid,    ONLY: iProcCIMI=>iProc,nProcCIMI=>nProc,iCommCIMI=>iComm
  use ModCoupleSami,  ONLY: cimi_set_global_mpi,cimi_get_init_for_sami, &
       cimi_send_to_sami,cimi_put_init_from_sami,cimi_get_from_sami
  use ModCoupleCimi,  ONLY: sami_set_global_mpi,sami_put_init_from_cimi, &
       sami_get_from_cimi,sami_send_to_cimi,sami_get_init_for_cimi
  use ModCIMI,        ONLY: IsStandalone,TimeCIMI=>time
  use ModMpi
  use ModCimi,        ONLY: init_mod_cimi
  use ModCimiTrace,  ONLY: init_mod_field_trace
  use ModImTime,      ONLY: TimeMaxCIMI=>TimeMax
  use ModCimiRestart, ONLY: DtSaveRestart,cimi_write_restart
  use ModReadParam
  use ModPrerunField, ONLY: UsePrerun, read_prerun, read_prerun_IE, DtRead
  use ModIeCimi,      ONLY: UseWeimer
  use CON_planet, ONLY: init_planet_const, set_planet_defaults
!  use ModPrerunField, ONLY: UsePrerun, read_prerun, read_prerun_IE
 
  implicit none
 

  ! This is hardcoded for now, CIMI runs on 24 procs and SAMI on 24 so 48 total
  integer,parameter :: nProcCIMItmp = 24, nProcSAMItmp = 24
!!$  integer,parameter :: nProcCIMItmp = 8, nProcSAMItmp = 9
!!$  integer,parameter :: nProcCIMItmp = 1, nProcSAMItmp = 9

  !set some mpi parameters
  integer :: iProcsCIMI_I(nProcCIMItmp),iProcsSAMI_I(nProcSAMItmp)
  integer :: iCommGlobal, nProcGloabl,iProcGlobal
  logical :: IsSamiProc, IsCimiProc
  integer :: CIMI_GROUP,SAMI_GROUP, MPI_GROUP_WORLD
  integer :: nProcGlobal

  !save the index of zero process for CIMI and SAMI
  integer :: iProc0CIMI, iProc0SAMI

  integer :: iError, iStep, iRank
  real    :: DtAdvance
  real    :: DtRestart = 300.0 ! currently set at 300s, should be read in later
  real    :: DtMax = 60.0      ! maximum timestep
  real    :: TimeMax,Time=0.0
  !---------------------------------------------------------------------------

  !****************************************************************************
  ! Initiallize MPI and get number of processors and rank of given processor
  !****************************************************************************
  
  !---------------------------------------------------------------------------
  call MPI_INIT(iError)
  iCommGlobal = MPI_COMM_WORLD
  
  call MPI_COMM_RANK(iCommGlobal,iProcGlobal,iError)
  call MPI_COMM_SIZE(iCommGlobal,nProcGlobal,iError)
  nProcCIMI = nProcCIMItmp
  nProcSAMI = nProcSAMItmp
  
  if (nProcGlobal /= nProcCIMI+nProcSAMI) &
       call con_stop('nProc not equal to nProcCIMI+nProcSAMI')
  
  ! get group id of MPI_COMM_WORLD so we can make sub groups
   call mpi_comm_group(MPI_COMM_WORLD,MPI_GROUP_WORLD,iError)
  
  ! create CIMI and SAMI proc lists and MPI groups
  do iRank = 0,nProcSAMI-1
     iProcsSAMI_I(iRank+1) = iRank
  enddo
  call mpi_group_incl(MPI_GROUP_WORLD,nProcSAMI,iProcsSAMI_I,SAMI_GROUP,iError)
  
  do iRank = nProcSAMI,nProcGlobal-1
     iProcsCIMI_I(iRank-nProcSAMI+1) = iRank
  enddo
  call mpi_group_incl(MPI_GROUP_WORLD,nProcCIMI,iProcsCIMI_I,CIMI_GROUP,iError)
  
  ! create CIMI and SAMI communicators from the proc groups
  call mpi_comm_create(iCommGlobal,SAMI_GROUP,iCommSAMI,iError)
  call mpi_comm_create(iCommGlobal,CIMI_GROUP,iCommCIMI,iError)
  
  ! assign sami and cimi rank
  if (iProcGlobal <nProcSAMI)then
     !set SAMI rank and put CIMI rank to -1
     call MPI_COMM_RANK(iCommSAMI,iProcSAMI,iError)
     iProcCIMI=-1
     IsSamiProc=.true.
     IsCimiProc=.false.
  else
     !set CIMI rank and put SAMI rank to -1
     call MPI_COMM_RANK(iCommCIMI,iProcCIMI,iError)
     iProcSAMI=-1
     IsSamiProc=.false.
     IsCimiProc=.true.
  endif
  
  ! set the zero CIMI and SAMI proc
  iProc0SAMI = 0
  iProc0CIMI = nProcSAMI 

  ! set the global mpi coupling variables in the coupling modules
  if (IsCimiProc) call cimi_set_global_mpi(iProc0Cimi,iProc0Sami, iCommGlobal)
  if (IsSamiProc) call sami_set_global_mpi(iProc0Cimi,iProc0Sami, iCommGlobal)
  !****************************************************************************
  ! Read the input file
  !****************************************************************************
  
  If(IsCimiProc) then
     IsStandAlone=.true.
     
     ! Initial setup for CIMI model and read PARAM.in
     call read_file('PARAM.in',iCommCIMI)
     call read_init('  ',iSessionIn=1,iLineIn=0)
     call CIMI_set_parameters('READ')
     
     !if (usePrerun) call read_prerun(t)
     !if (usePrerun .and. iConvect==2) call read_prerun_IE(t)
     

  endif
  !\
  ! Initialize the planetary constant library and set Earth
  ! as the default planet.
  !/
  call init_planet_const
  call set_planet_defaults
  

  ! TimeMax is set in CIMI PARAM.in file and should be sent to SAMI
  if (iProcGlobal == iProc0CIMI) TimeMax=TimeMaxCIMI
  call MPI_bcast(TimeMax, 1, MPI_REAL, iProc0CIMI, iCommGlobal, iError)

  ! The simulation time is set in CIMI param file and should be sent to SAMI
  if (iProcGlobal == iProc0CIMI) Time=TimeCIMI
  call MPI_bcast(Time, 1, MPI_REAL, iProc0CIMI, iCommGlobal, iError)


  !****************************************************************************
  ! Initialize the model
  !****************************************************************************
  If(IsCimiProc) then
     call init_mod_cimi
     call init_mod_field_trace
  endif
  
  
  ! Set timing procs
  If(IsCimiProc) then
     ! Set the component name and PE number
     call timing_comp_proc('CIMI',iProcGlobal) 
  else
     call timing_comp_proc('SAMI',iProcGlobal) 
  endif

  ! Start Timing
  call timing_active(.true.)
  call timing_step(0)
  call timing_start('CIMI-SAMI')
  
  If(IsCimiProc) then
     !read initial prerun field
     if (UsePrerun) call read_prerun(Time)
     if (UsePrerun .and. .not.UseWeimer) call read_prerun_IE(Time)
  
     ! init CIMI model
     call timing_start('cimi_init')
     call cimi_init
     call timing_stop('cimi_init')
  else
     !init SAMI model
     call timing_start('sami_init')
     call sami_init
     call timing_stop('sami_init')
     
  endif

  call MPI_BARRIER(iCommGlobal,iError)

  if (iProcGlobal==iProc0CIMI) call cimi_get_init_for_sami

  if (IsSamiProc)              call sami_put_init_from_cimi

  
  if (iProcGlobal==iProc0SAMI) call sami_get_init_for_cimi

  if (IsCimiProc)              call cimi_put_init_from_sami

  

  

!  if (iProcGlobal == 0)call timing_report_total
  if (iProcGlobal == iProc0SAMI .or. iProcGlobal == iProc0CIMI) then 
     call timing_report_total
  endif
  call timing_reset('#all',3)
  
  ! couple one time before timestepping
  if(IsCimiProc) then
     call cimi_send_to_sami

     call cimi_get_from_sami(Time)

  endif
  if(IsSamiProc) then
     call sami_get_from_cimi(Time)
     call sami_send_to_cimi
  endif

  !****************************************************************************
  ! start Timestepping
  !****************************************************************************
  iStep = 0
  TIMELOOP:do
     !report progress on proc 0
     if (iProcGlobal==0)write(*,*) 'In Time Loop iStep,Time = ', iStep,Time
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
     
     If(IsCimiProc) then
        ! Call cimi_run to advance the Timestep
        call timing_step(iStep)
        call timing_start('cimi_run')     
        call cimi_run(DtAdvance)
        call timing_stop('cimi_run')
     else
        ! Call cimi_run to advance the Timestep
        call timing_step(iStep)
        call timing_start('sami_run')     
        call sami_run(DtAdvance)
        call timing_stop('sami_run')
        write(*,*) 'Finished sami_run for iProcSAMI = ',iProcSAMI
     endif

     !here we put the coupling
     if(IsCimiProc) then
        call cimi_send_to_sami
        call cimi_get_from_sami(Time)
     endif
     if(IsSamiProc) then
        call sami_get_from_cimi(Time)
        call sami_send_to_cimi
     endif

     ! Save restart at DtSaveRestart or TimeMax
     if (floor((Time+1.0e-5)/DtSaveRestart) /= &
          floor((Time+1.0e-5-DtAdvance)/DtSaveRestart)) then
        If(IsCimiProc) then
           call cimi_write_restart
        else
           !call sami_write_restart
        endif
     endif
     
     If(IsCimiProc) then
        ! Read new prerun file if using prerun fields and it is time to read
        if (floor((Time+1.0e-5)/DtRead) /= &
             floor((Time+1.0e-5-DtAdvance)/DtRead)) then 
           if (UsePrerun) call read_prerun(Time)
           if (UsePrerun .and. .not.UseWeimer) call read_prerun_IE(Time)     
        endif
     endif
     ! Advance the time and iteration step
     iStep=iStep+1
     Time = Time+DtAdvance
  end do TIMELOOP
  write(*,*) 'Finished timestepping for iProcGlobal=', iProcGlobal

  ! Save restart at TimeMax
  If(IsCimiProc) then 
     call cimi_write_restart
  else
     !call sami_write_restart
  endif
  
  ! Finalize timing commands
  call timing_stop('CIMI-SAMI')

  call MPI_BARRIER(iCommGlobal,iError)
!  if (iProcGlobal == 0) then
  if (iProcGlobal == iProc0SAMI) then 
     write(*,'(a)') 'Finished CIMI_SAMI run, (reporting timings SAMI)'
     write(*,'(a)') '--------------------------------------'
     call timing_report
  endif

  if (iProcGlobal == iProc0CIMI) then 
     write(*,'(a)') 'Finished CIMI_SAMI run, (reporting timings CIMI)'
     write(*,'(a)') '--------------------------------------'
     call timing_report
  endif
  
  If(IsSamiProc) call sami_finalize

  ! Finalize MPI
  call MPI_FINALIZE(iError)

end program cimi_sami
!============================================================================
subroutine CON_stop(StringError)
  use ModCimiGrid,    ONLY: iProc,nProc
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
  call MPI_abort(MPI_COMM_WORLD, nError, iError)
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
