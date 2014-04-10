!^CFG COPYRIGHT UM

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!               Space Weather Modeling Framework (SWMF)                !
!    Center for Space Environment Modeling, The University of Michigan !
!-----------------------------------------------------------------------
subroutine CZ_set_param(CompInfo, TypeAction)
  
  !USES:
  use ModPar
  use CON_comp_info,      ONLY: get, put, CompInfoType
  use CON_time,           ONLY: MaxIteration, tSimulationMax, CpuTimeMax
  use ModFSAM,            ONLY: iComm
  use ModInitialization,  ONLY: atmoin, input_param
  use ModControl,         ONLY: itmax, wtlimit, tend
  use ModPhysunits,       ONLY: unit_t
  implicit none

  character (len=*), intent(in)     :: TypeAction ! which action to perform
  type(CompInfoType), intent(inout) :: CompInfo   ! component information
  
  character (len=*), parameter :: NameSub='CZ_set_param'
  integer :: iError, isize
  character (len=100) :: NameCommand, NameItNow

  !------------------------------------------------------------------------
  
  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,                                     &
          Use=.true.,                                       &
          NameVersion='Finite-difference Spherical Anelastic MHD, Y. Fan', &
          Version=1.0)
  case('STDOUT')
     if(myid==0) then
        write(6,'(a)') '------------------------ FSAM starts -------------------------- '
        write(6,'(a,i4,a)') ' Running FSAM on ', nproc, ' CPUs'
        write(6,'(3(a,i4))') ' inmax=', inmax, ' , jnmax=', jnmax, ', knmax=', knmax
        write(6,'(3(a,i4))') ' in = ', in, ', jn = ', jn, ', kn = ', kn
        write(6,'(a)') '--------------------------------------------------------------- '
        write(6,*)
        call flush(6)
     endif
  case('FILEOUT')
     !call get(CompInfo,iUnitOut=iUnitOut)
     !StringPrefix=''
  case('MPI')
     call get(CompInfo, iComm=iComm, iProc=myid, nProc=isize)
     ! check processors
     if(myid==0.and.(isize .ne. nproc)) then
        write(*, *) 'wrong number of processors'
        write(*, *) nproc, ' specified in par.h, ', isize, ' allocated'
        call MPI_ABORT(iComm, 1,iError)
     endif
     call flush(6)
     myid2 = myid/nproc1
     myid1 = mod(myid,nproc1)
     ! input solar atmosphere and define units
     call atmoin
  case('READ')
     ! read PARAM.in
     call input_param
     ! set control parameter for the run
     MaxIteration = itmax
     tSimulationMax = tend*unit_t
     CpuTimeMax = wtlimit
  case('CHECK')
  case('GRID')
  case default
     call CON_stop(NameSub//' CZ_ERROR: invalid TypeAction='//TypeAction)
  end select
  
end subroutine CZ_set_param

!==============================================================================

subroutine CZ_init_session(iSession, TimeSimulation)
  use ModUserSetup,  ONLY: field_init, grid
  use ModBval,       ONLY: bvalv
  use ModInitialization
  use ModGetQtys
  use ModFSAM
  use ModPar
  use ModIoFSAM
  use ModPhysunits
  use ModIteration
  use ModOutfile
  use ModRHS,     ONLY: p_init
  use ModSundry,  ONLY: dt, time, ctime
  use ModControl, ONLY: itnow, tout
  implicit none

  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='CZ_init_session'

  real :: em, ek, ek1, ek2, ek3, eth, anglm, anglmnorm, smeanbot, smeantop
  integer :: iError, ntcond
  !------------------------------------------------------------------------

  ! define grid
  call grid

  ! define blk matrix 
  call blk_init

  ! initialize sundry.h and initial fields 
  call field_init
  if(.not. DoRestart) then
     itnow = 0
     time  = 0.D0
     ctime = time
  else
     call readrst_mpi
     ctime = time
     if(myid==0) write(6,*) 'Restart from restart files at time t = ', time
     call bvalv
  endif
  call p_init
  call nudt

  ! open monitoring file
  ntcond = 0
  call get_global(em,ek,ek1,ek2,ek3,eth,anglm,anglmnorm)
  call get_mean_entropy(smeanbot,smeantop)
  if(myid==0) then
     open(unit=16,file=engfile,status='replace')
     write(16,'(i8,13e23.14,i8)') itnow,time,dt,em,ek,eth, &
          ek+eth+em,smeanbot,smeantop,anglm,anglmnorm, &
          ek1,ek2,ek3,ntcond
     call flush(16)
  endif

  ! output grid and initial fields and set up for future outputs 
  call writegrid
  if(ifile .lt. 0) then
     ifile = time/tout
     tfile = tout*ifile
  else
     tfile=time
  endif
  if(.not.DoRestart) call writedata_mpi
  ifile = ifile + 1
  tfile = tfile + tout

  if(UseImplCond) then
     if(myid==0) write(*,*) 'Initialize Implict Heat Conduction Solver..........'
     ! call init_impl_cond
  else
     if(myid==0) write(*,*) 'Solving Conduction Explicitly......................'
     call nudt_cond
  endif

end subroutine CZ_init_session

!==============================================================================

subroutine CZ_run(TimeSimulation,TimeSimulationLimit)
  use ModInitialization
  use ModGetQtys
  use ModFSAM
  use ModPar
  use ModIoFSAM
  use ModPhysunits
  use ModIteration
  use ModOutfile
  use ModSundry,  ONLY: dt, time, ctime
  use ModControl, ONLY: itnow, tout, itintv, itout
  implicit none

  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded
  real, intent(inout) :: TimeSimulation   ! current time of component
  
  character(len=*), parameter :: NameSub='CZ_run'
  
  real :: em, ek, ek1, ek2, ek3, eth, anglm, anglmnorm, smeanbot, smeantop

  integer :: ntcond
  !------------------------------------------------------------------------

  call nudt
  call pred_corr_step
  if(UseImplCond)then
     !  it = itnow                                                                       
     !  call conduct_impl(ntcond)                         
  else
     call conductstep(ntcond)
  endif

  time = time + dt
  ctime = time
  itnow = itnow + 1
  
  ! monitoring global quantities                                                       
  if(mod(itnow,itintv)==0) then
     call get_global(em,ek,ek1,ek2,ek3,eth,anglm,anglmnorm)
     call get_mean_entropy(smeanbot,smeantop)
     if(myid==0) then
        write(16,'(i8,13e23.14,i8)') itnow,time,dt,em,ek,eth, &
             ek+eth+em,smeanbot,smeantop,anglm,anglmnorm, ek1,ek2,ek3,ntcond
        call flush(16)
     endif
  endif
  
  ! check for output time                                                              
  if(mod(itnow,itout)==itout-1) call writedata_mpi(DoWriteRestart=.true.)
  if(time > tfile) then
     call writedata_mpi
     ifile = ifile + 1
     tfile = tfile + tout
  endif

  TimeSimulation = time*unit_t

end subroutine CZ_run
!===========================================================================

subroutine CZ_finalize(TimeSimulation)

  !USES:
  implicit none

  real,     intent(in) :: TimeSimulation   ! seconds from start time
  character(len=*), parameter :: NameSub='CZ_finalize'

  !-------------------------------------------------------------------------

end subroutine CZ_finalize
!===========================================================================

subroutine CZ_save_restart(TimeSimulation)
  use ModIoFSAM,  ONLY: writedata_mpi
  implicit none

  real,     intent(in) :: TimeSimulation   ! seconds from start time
  character(len=*), parameter :: NameSub='CZ_save_restart'

  !-------------------------------------------------------------------------

  call writedata_mpi(DoWriteRestart=.true.)

end subroutine CZ_save_restart
!===========================================================================

