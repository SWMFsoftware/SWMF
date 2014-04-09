program xfsam
  use ModPar
  use ModOutfile
  use ModControl
  use ModSundry
  use ModGetQtys
  use ModIO
  use ModInitialization
  use ModIteration
  !use ModImplCond
  use ModFSAM
  use ModMpi
  implicit none
  
  integer :: isize, ierr, ntcond, update_rate
  real :: wt0, wtnow, wtelpslc, wtelps, wtprev
  real :: em, ek, eth, anglm, anglmnorm, ek1, ek2, ek3, smeanbot, smeantop
  !---------------------------------------------------------------------------------
  
  IsStandAlone = .true.
  !****************************************************************************
  ! Initiallize MPI and get number of processors and rank of given processor
  !****************************************************************************
  call MPI_Init(ierr)
  iComm = MPI_COMM_WORLD
  
  call MPI_Comm_rank(iComm,myid,ierr)
  call MPI_Comm_size(iComm,isize,ierr)

! check processors
  if(isize .ne. nproc) then
    write(6, *) 'wrong number of processors'
    write(6, *) nproc, ' specified in par.h, ', isize, ' allocated'
    call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
  endif
  call flush(6)

  ! start wall clock
  wt0 = MPI_WTIME()
  
  if(myid==0) then
     write(6,'(a)') '------------------------ FSAM starts -------------------------- '
     write(6,'(a,i4,a)') ' Running FSAM on ', nproc, ' CPUs'
     write(6,'(3(a,i4))') ' inmax = ', inmax, ' , jnmax = ', jnmax, ', knmax = ', knmax
     write(6,'(3(a,i4))') ' in = ', in, ', jn = ', jn, ', kn = ', kn
     write(6,'(a)') '--------------------------------------------------------------- '
     write(6,*) 
     call flush(6)
  endif
  
  myid2 = myid/nproc1
  myid1 = mod(myid,nproc1)
  call initialize
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
  
  !****************************************************************************
  !      main iteration loop
  !****************************************************************************

  ! compute conduction time step (only compute once for fixed ovrth)
  if(UseImplCond) then
     if(myid==0) write(*,*) 'Initialize Implict Heat Conduction Solver..........'
     !     call init_impl_cond
  else
     if(myid==0) write(*,*) 'Solving Conduction Explicitly......................'
     call nudt_cond
  endif

  ! start iteration
  wtnow = MPI_WTIME()
  wtprev = wtnow
  wtelpslc = wtnow - wt0
  call MPI_ALLREDUCE(wtelpslc, wtelps, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
       iComm,ierr)
  
  do while (time .lt. tend .and. itnow .lt. itmax .and. wtelps .lt. wtlimit)
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
     wtnow = MPI_WTIME()
     wtelpslc = wtnow-wt0
     call MPI_ALLREDUCE(wtelpslc,wtelps,1, MPI_DOUBLE_PRECISION,MPI_MAX, &
          iComm,ierr)

     ! monitoring global quantities
     if(mod(itnow,itintv)==0) then
        call get_global(em,ek,ek1,ek2,ek3,eth,anglm,anglmnorm)
        call get_mean_entropy(smeanbot,smeantop)
        if(myid==0) then
           write(16,'(i8,13e23.14,i8)') itnow,time,dt,em,ek,eth, &
                ek+eth+em,smeanbot,smeantop,anglm,anglmnorm, ek1,ek2,ek3,ntcond
           call flush(16)
           update_rate = itintv*inmax*jnmax*knmax/nproc/(wtnow - wtprev)
           write(6,'(a,I7,a,I9)') ' Grid cell updates per CPU sec = ', &
                update_rate, ' at iteration = ', itnow
           wtprev = wtnow
           call flush(6)
        endif
     endif
     
     ! check for output time
     if(mod(itnow,itout)==itout-1) call writedata_mpi(DoWriteRestart=.true.)
     if(time > tfile) then
        call writedata_mpi
        ifile = ifile + 1
        tfile = tfile + tout
     endif
  enddo

  !****************************************************************************
  !      finalize
  !****************************************************************************
  call get_global(em,ek,ek1,ek2,ek3,eth,anglm,anglmnorm)
  call get_mean_entropy(smeanbot,smeantop)
  if(myid==0) then
     write(16,'(i8,13e23.14,i8)') itnow,time,dt,em,ek,eth, &
          ek+eth+em,smeanbot,smeantop,anglm,anglmnorm, ek1,ek2,ek3,ntcond
     close(16) 
     close(6)
  endif
  call writedata_mpi
  call writedata_mpi(DoWriteRestart=.true.)
  call MPI_finalize(ierr)

end program xfsam

!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  use ModPar,      ONLY: myid
  use ModSundry,   ONLY: time
  use ModMpi
  use ModFSAM
  implicit none
  character (len=*), intent(in) :: StringError
  
  ! Local variables:
  integer :: iError, nError
  !----------------------------------------------------------------------------

  write(*,*) 'Stopping execution! me=',myid,' at time=', time, ' with msg:'
  write(*,*) StringError
  call MPI_abort(iComm, nError, iError)
  stop

end subroutine CON_stop
