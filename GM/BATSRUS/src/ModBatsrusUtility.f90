!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!==== Simple subroutines and functions that operate on all used blocks ========
real function integrate_BLK(qnum,qa)

  use BATL_lib, ONLY: &
       test_start, test_stop

  ! Return the volume integral of qa, ie. sum(qa*CellVolume_B)
  ! for all used blocks and true cells
  ! Do for each processor separately if qnum=1, otherwise add them all

  use ModProcMH
  use ModGeometry, ONLY : true_BLK, true_cell
  use BATL_lib, ONLY: IsCartesian, CellVolume_GB, CellVolume_B, Unused_B, &
       nI, nJ, nK, nBlock, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, MaxBlock
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum
  real,    intent(in) :: qa(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

  ! Local variables:
  real    :: qsum, qsum_all
  integer :: iBlock, iError

  logical:: DoTest
  character(len=*), parameter:: NameSub = 'integrate_BLK'
  !----------------------------------------------------------------------------
  call test_start(NameSub, DoTest)

  qsum=0.0

  if(IsCartesian)then
     do iBlock = 1, nBlock
        if(Unused_B(iBlock)) CYCLE
        if(true_BLK(iBlock)) then
           qsum = qsum + CellVolume_B(iBlock)* &
                sum(qa(1:nI,1:nJ,1:nK,iBlock))
        else
           qsum = qsum + CellVolume_B(iBlock)* &
                sum( qa(1:nI,1:nJ,1:nK,iBlock), &
                MASK=true_cell(1:nI,1:nJ,1:nK,iBlock))
        end if
     end do
  else
     do iBlock = 1, nBlock
        if(Unused_B(iBlock)) CYCLE
        if(true_BLK(iBlock)) then
           qsum=qsum + sum(CellVolume_GB(1:nI,1:nJ,1:nK,iBlock)* &
                qa(1:nI,1:nJ,1:nK,iBlock))
        else
           qsum=qsum + sum(CellVolume_GB(1:nI,1:nJ,1:nK,iBlock)* &
                qa(1:nI,1:nJ,1:nK,iBlock), &
                MASK=true_cell(1:nI,1:nJ,1:nK,iBlock))
        end if
     end do
  end if

  if(qnum>1)then
     call MPI_allreduce(qsum, qsum_all, 1,  MPI_REAL, MPI_SUM, &
          iComm, iError)
     integrate_BLK = qsum_all
     if(DoTest)write(*,*)'me,sum,sum_all:',iProc,qsum,qsum_all
  else
     integrate_BLK = qsum
     if(DoTest)write(*,*)'me,qsum:',iProc,qsum
  end if

  call test_stop(NameSub, DoTest)
end function integrate_BLK
!==============================================================================

real function minval_BLK(qnum,qa)

  use BATL_lib, ONLY: &
       test_start, test_stop

  ! Return minval(qa)) corresponding to all used blocks
  ! If qnum<=1, return minval for the processor, otherwise
  ! return the minimum for all processors.

  use ModProcMH
  use ModSize, ONLY: &
       nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nBlock, MaxBlock
  use ModMain, ONLY: Unused_B
  use ModGeometry, ONLY : true_BLK,true_cell
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum

  real, intent(in) :: qa(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

  ! Local variables:
  real    :: qminval, qminval_all
  integer :: iBlock, iError

  logical:: DoTest
  character(len=*), parameter:: NameSub = 'minval_BLK'
  !----------------------------------------------------------------------------
  call test_start(NameSub, DoTest)

  qminval=1.e+30

  do iBlock=1,nBlock
     if(.not.Unused_B(iBlock)) then
        if(true_BLK(iBlock)) then
           qminval=min(qminval, minval(qa(1:nI,1:nJ,1:nK,iBlock)))
        else
           qminval=min(qminval, minval(qa(1:nI,1:nJ,1:nK,iBlock),&
                MASK=true_cell(1:nI,1:nJ,1:nK,iBlock)))
        endif
     end if
  end do

  if(qnum>1)then
     call MPI_allreduce(qminval, qminval_all, 1,  MPI_REAL, MPI_MIN, &
          iComm, iError)
     minval_BLK=qminval_all
     if(DoTest)write(*,*)'me,minval,minval_all:',iProc,qminval,qminval_all
!     if(qminval_all>=1.e+30)call stop_mpi('Error in minval_BLK: NaN-s?')
  else
     minval_BLK=qminval
     if(DoTest)write(*,*)'me,qminval:',iProc,qminval
!     if(qminval>=1.e+30)call stop_mpi('Error in minval_BLK: NaN-s?')
  endif

  call test_stop(NameSub, DoTest)
end function minval_BLK
!==============================================================================

real function maxval_BLK(qnum,qa)

  use BATL_lib, ONLY: &
       test_start, test_stop

  ! Return maxval(qa)) corresponding to all used blocks
  ! If qnum<=1, return maxval for the processor, otherwise
  ! return the maximum for all processors.

  use ModProcMH
  use ModSize, ONLY: &
       nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nBlock, MaxBlock
  use ModMain, ONLY: Unused_B
  use ModGeometry, ONLY: true_BLK,true_cell
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum

  real, intent(in) :: qa(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

  ! Local variables:
  real    :: qmaxval, qmaxval_all
  integer :: iBlock, iError

  logical:: DoTest
  character(len=*), parameter:: NameSub = 'maxval_BLK'
  !----------------------------------------------------------------------------
  call test_start(NameSub, DoTest)

  qmaxval=-1.e+30

  do iBlock=1,nBlock
     if(.not.Unused_B(iBlock)) then
        if(true_BLK(iBlock)) then
           qmaxval=max(qmaxval, maxval(qa(1:nI,1:nJ,1:nK,iBlock)))
        else
           qmaxval=max(qmaxval, maxval(qa(1:nI,1:nJ,1:nK,iBlock),&
                MASK=true_cell(1:nI,1:nJ,1:nK,iBlock)))
        endif
     end if
  end do

  if(qnum>1)then
     call MPI_allreduce(qmaxval, qmaxval_all, 1,  MPI_REAL, MPI_MAX, &
          iComm, iError)
     maxval_BLK=qmaxval_all
     if(DoTest)write(*,*)'me,maxval,maxval_all:',iProc,qmaxval,qmaxval_all
     if(qmaxval_all<=-1.e+30)call stop_mpi('Error in maxval_BLK: NaN-s?')
  else
     maxval_BLK=qmaxval
     if(DoTest)write(*,*)'qmaxval:',qmaxval
     if(qmaxval<=-1.e+30)call stop_mpi('Error in maxval_BLK: NaN-s?')
  endif

  call test_stop(NameSub, DoTest)
end function maxval_BLK
!==============================================================================

real function maxval_loc_BLK(qnum,qa,loc)

  use BATL_lib, ONLY: &
       test_start, test_stop

  ! Return maxval(qa)) corresponding to all used blocks and
  ! also return the location of the maximum value into loc(5)=I,J,K,iBlock,PE
  ! If qnum<=1, return maxval for the processor, otherwise
  ! return the maximum for all processors.

  use ModProcMH
  use ModSize, ONLY: &
       nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nBlock, MaxBlock
  use ModMain, ONLY: Unused_B
  use ModGeometry, ONLY:true_cell
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum

  real, intent(in) :: qa(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

  integer, intent(out):: loc(5)

  ! Local variables:
  real    :: qmaxval, qmaxval_all
  integer :: i,j,k,iBlock, iError

  real, external :: maxval_BLK

  logical:: DoTest
  character(len=*), parameter:: NameSub = 'maxval_loc_BLK'
  !----------------------------------------------------------------------------
  call test_start(NameSub, DoTest)

  qmaxval=maxval_BLK(1,qa)
  if(qnum==1)then
     qmaxval_all=qmaxval
  else
     call MPI_allreduce(qmaxval, qmaxval_all, 1,  MPI_REAL, MPI_MAX, &
          iComm, iError)
  end if

  loc=-1
  if (qmaxval == qmaxval_all) then
     BLKLOOP: do iBlock=1,nBlock
        if(Unused_B(iBlock)) CYCLE
        do k=1,nK; do j=1,nJ; do i=1,nI
           if(.not.true_cell(i,j,k,iBlock)) CYCLE
           if(qa(i,j,k,iBlock)==qmaxval)then
              loc(1)=i; loc(2)=j; loc(3)=k; loc(4)=iBlock; loc(5)=iProc
              EXIT BLKLOOP
           end if
        enddo; enddo; enddo;
     enddo BLKLOOP
  end if

  maxval_loc_BLK=qmaxval_all

  call test_stop(NameSub, DoTest)
end function maxval_loc_BLK
!==============================================================================
real function maxval_loc_abs_BLK(qnum,qa,loc)

  use BATL_lib, ONLY: &
       test_start, test_stop

  ! Return maxval(abs(qa)) corresponding to all used blocks and
  ! also return the location of the maximum value into loc(5)=I,J,K,iBlock,PE
  ! If qnum<=1, return maxval for the processor, otherwise
  ! return the maximum for all processors.

  use ModProcMH
  use ModSize, ONLY: &
       nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nBlock, MaxBlock
  use ModMain, ONLY: Unused_B
  use ModGeometry, ONLY : true_cell
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum

  real, intent(in) :: qa(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

  integer, intent(out):: loc(5)

  ! Local variables:
  real    :: qmaxval, qmaxval_all
  integer :: i,j,k,iBlock, iError

  real, external :: maxval_abs_BLK

  logical:: DoTest
  character(len=*), parameter:: NameSub = 'maxval_loc_abs_BLK'
  !----------------------------------------------------------------------------
  call test_start(NameSub, DoTest)

  qmaxval=maxval_abs_BLK(1,qa)
  if(qnum==1)then
     qmaxval_all=qmaxval
  else
     call MPI_allreduce(qmaxval, qmaxval_all, 1,  MPI_REAL, MPI_MAX, &
          iComm, iError)
  end if

  loc=-1
  if (qmaxval == qmaxval_all) then
     BLKLOOP: do iBlock=1,nBlock
        if(Unused_B(iBlock)) CYCLE
        do k=1,nK; do j=1,nJ; do i=1,nI
           if(.not.true_cell(i,j,k,iBlock)) CYCLE
           if(abs(qa(i,j,k,iBlock))==qmaxval)then
              loc(1)=i; loc(2)=j; loc(3)=k; loc(4)=iBlock; loc(5)=iProc
              EXIT BLKLOOP
           end if
        enddo; enddo; enddo;
     enddo BLKLOOP
  end if

  maxval_loc_abs_BLK=qmaxval_all

  call test_stop(NameSub, DoTest)
end function maxval_loc_abs_BLK
!==============================================================================
real function minval_loc_BLK(qnum,qa,loc)

  use BATL_lib, ONLY: &
       test_start, test_stop

  ! Return minval(qa)) corresponding to all used blocks and
  ! also return the location of the minimum value into loc(5)=I,J,K,iBlock,PE
  ! If qnum<=1, return minval for the processor, otherwise
  ! return the minimum for all processors.

  use ModProcMH
  use ModSize, ONLY: &
       nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nBlock, MaxBlock
  use ModMain, ONLY: Unused_B
  use ModGeometry, ONLY : true_cell
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum

  real, intent(in) :: qa(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

  integer, intent(out):: loc(5)

  ! Local variables:
  real    :: qminval, qminval_all
  integer :: i,j,k,iBlock, iError

  real, external :: minval_BLK

  logical:: DoTest
  character(len=*), parameter:: NameSub = 'minval_loc_BLK'
  !----------------------------------------------------------------------------
  call test_start(NameSub, DoTest)

  qminval=minval_BLK(1,qa)
  if(qnum==1)then
     qminval_all=qminval
  else
     call MPI_allreduce(qminval, qminval_all, 1,  MPI_REAL, MPI_MIN, &
          iComm, iError)
  end if

  loc=-1
  if (qminval == qminval_all) then
     BLKLOOP:do iBlock=1,nBlock
        if(Unused_B(iBlock)) CYCLE
        do k=1,nK; do j=1,nJ; do i=1,nI
           if(.not.true_cell(i,j,k,iBlock)) CYCLE
           if(qa(i,j,k,iBlock)==qminval)then
              loc(1)=i; loc(2)=j; loc(3)=k; loc(4)=iBlock; loc(5)=iProc
              EXIT BLKLOOP
           end if
        enddo; enddo; enddo;
     enddo BLKLOOP
  end if

  minval_loc_BLK=qminval_all

  call test_stop(NameSub, DoTest)
end function minval_loc_BLK
!==============================================================================

real function maxval_abs_BLK(qnum,qa)

  use BATL_lib, ONLY: &
       test_start, test_stop

  ! Return maxval(abs(qa)) corresponding to all used blocks
  ! If qnum<=1, return maxval(abs(qa)) for the processor, otherwise
  ! return the maximum for all processors.

  use ModProcMH
  use ModSize, ONLY: &
       nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nBlock, MaxBlock
  use ModMain, ONLY: Unused_B
  use ModGeometry, ONLY : true_cell,true_BLK
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum

  real, intent(in) :: qa(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

  ! Local variables:
  real    :: qmaxval, qmaxval_all
  integer :: iBlock, iError

  logical:: DoTest
  character(len=*), parameter:: NameSub = 'maxval_abs_BLK'
  !----------------------------------------------------------------------------
  call test_start(NameSub, DoTest)

  qmaxval=-1.0

  do iBlock=1,nBlock
     if(.not.Unused_B(iBlock)) then
        if(true_BLK(iBlock)) then
           qmaxval=max(qmaxval, maxval(abs(qa(1:nI,1:nJ,1:nK,iBlock))))
        else
           qmaxval=max(qmaxval, maxval(abs(qa(1:nI,1:nJ,1:nK,iBlock)),&
                MASK=true_cell(1:nI,1:nJ,1:nK,iBlock)))
        endif
     end if
  end do

  if(qnum>1)then
     call MPI_allreduce(qmaxval, qmaxval_all, 1,  MPI_REAL, MPI_MAX, &
          iComm, iError)
     maxval_abs_BLK=qmaxval_all
     if(DoTest)write(*,*)'me,maxval,maxval_all:',iProc,qmaxval,qmaxval_all
     if(qmaxval_all<0.0) then
        call barrier_mpi
        call stop_mpi('Error in maxval_abs_BLK: negative max!')
     end if
  else
     maxval_abs_BLK=qmaxval
     if(DoTest)write(*,*)'qmaxval:',qmaxval
     if(qmaxval<0.0)call stop_mpi('Error in maxval_abs_BLK: negative max!')
  endif

  call test_stop(NameSub, DoTest)
end function maxval_abs_BLK
!==============================================================================

real function maxval_abs_ALL(qnum,qa)

  use BATL_lib, ONLY: &
       test_start, test_stop

  ! Return maxval(abs(qa)) corresponding to all used blocks
  ! Include ghost cells and .not.true_cell -s too.
  ! If qnum<=1, return maxval(abs(qa)) for the processor, otherwise
  ! return the minimum for all processors.

  use ModProcMH
  use ModSize, ONLY: &
       MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nBlock, MaxBlock
  use ModMain, ONLY: Unused_B
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum

  real, intent(in):: qa(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

  ! Local variables:
  real    :: qmaxval, qmaxval_all
  integer :: iBlock, iError

  logical:: DoTest
  character(len=*), parameter:: NameSub = 'maxval_abs_ALL'
  !----------------------------------------------------------------------------
  call test_start(NameSub, DoTest)

  qmaxval=-1.0

  do iBlock=1,nBlock
     if(.not.Unused_B(iBlock)) then
        qmaxval=max(qmaxval, maxval(abs(qa(:,:,:,iBlock))))
     end if
  end do

  if(qnum>1)then
     call MPI_allreduce(qmaxval, qmaxval_all, 1,  MPI_REAL, MPI_MAX, &
          iComm, iError)
     maxval_abs_ALL=qmaxval_all
     if(DoTest)write(*,*)'me,maxval,maxval_all:',iProc,qmaxval,qmaxval_all
  else
     maxval_abs_ALL=qmaxval
     if(DoTest)write(*,*)'qmaxval:',qmaxval
  endif

  call test_stop(NameSub, DoTest)
end function maxval_abs_ALL

!==============================================================================

subroutine barrier_mpi

  use ModProcMH
  use ModMpi
  implicit none

  ! Local variables:
  integer :: iError
  !----------------------------------------------------------------------------
  call timing_start('barrier')
  call MPI_barrier(iComm, iError)
  call timing_stop('barrier')

end subroutine barrier_mpi
!==============================================================================

subroutine barrier_mpi2(str)

  use ModProcMH
  use ModMpi
  implicit none

  character (len=*), intent(in) :: str

  ! Local variables:
  integer :: iError
  !----------------------------------------------------------------------------
  call timing_start('barrier-'//str)
  call MPI_barrier(iComm, iError)
  call timing_stop('barrier-'//str)

end subroutine barrier_mpi2
!==============================================================================

subroutine stop_mpi(str)

  use ModProcMH
  use ModMain, ONLY : iteration_number, NameThisComp, IsStandAlone, &
       DoWriteCallSequence
  use ModMpi
  implicit none

  character (len=*), intent(in) :: str

  ! Local variables:
  integer :: iError,nError
  !----------------------------------------------------------------------------
  if(IsStandAlone)then
     write(*,*)'Stopping execution! me=',iProc,' at iteration=',&
          iteration_number,' with msg:'
     write(*,*)str

     if(DoWriteCallSequence)then
        write(*,*)'Making floating point exception to write call sequence!'
        write(*,*) sqrt(-1.0-iteration_number)
     end if

     call MPI_abort(iComm, nError, iError)
     stop
  else
     write(*,*)NameThisComp,': stopping execution! me=',iProc,&
          ' at iteration=',iteration_number

     if(DoWriteCallSequence)then
        write(*,*)'Making floating point exception to write call sequence!'
        write(*,*) sqrt(-1.0-iteration_number)
     end if

     call CON_stop(NameThisComp//':'//str)
  end if

end subroutine stop_mpi
!==============================================================================
subroutine error_report(str,value,iErrorIn,show_first)

  use ModProcMH
  use ModMain, ONLY : iteration_number
  use ModIO, ONLY: write_myname
  use ModMpi
  implicit none

  ! Collect global error reports
  ! Reports are identified by an individual string str if iErrorIn<1
  !    and iErrorIn is set to a value > 1 for later use.
  ! If iErrorIn is > 1 to start with, it is used for error identification.
  ! Make statistics of errors based on value
  ! print statistics if str='PRINT'

  ! Parameters:

  ! Maximum number of different types of errors
  integer, parameter :: maxerror=20

  ! Arguments:

  character (LEN=*), intent(in) :: str
  real, intent(in)              :: value
  integer, intent(inout)        :: iErrorIn
  logical, intent(in)           :: show_first

  ! Local variables:

  integer :: iError

  ! Current number of different types of errors
  integer :: nErrors=0

  ! Message, number, occurance, and statistics of errors
  character (LEN=60), dimension(maxerror), save :: error_message
  integer, dimension(maxerror):: &
       error_count=0, error_count_sum, error_count_max,&
       iter_first=100000, iter_last=-1
  real,    dimension(maxerror):: &
       error_min=1e30, error_max=-1e30, &
       error_mean=0., error_last=0., error_last_sum

  character (LEN=60) :: msg

  integer,            dimension(:),   allocatable :: nErrors_all
  character (LEN=60), dimension(:,:), allocatable :: error_message_all
  integer,            dimension(:,:), allocatable :: &
       error_count_all, iter_first_all, iter_last_all
  real,               dimension(:,:), allocatable :: &
       error_min_all, error_max_all, error_mean_all, error_last_all

  integer :: i,i0,ip

  ! Debug
  ! write(*,*)'Error_report me, iErrorIn, value, str=',iProc,iErrorIn,value,str

  !----------------------------------------------------------------------------
  if(str=='PRINT')then
     ! Allocate memory in PROC 0
     allocate(&
          nErrors_all(nProc),&
          error_message_all(maxerror,nProc),&
          error_count_all(maxerror,nProc),&
          iter_first_all(maxerror,nProc),&
          iter_last_all(maxerror,nProc),&
          error_min_all(maxerror,nProc),&
          error_max_all(maxerror,nProc),&
          error_mean_all(maxerror,nProc),&
          error_last_all(maxerror,nProc))

     ! Collect the error reports
     call MPI_gather(nErrors, 1, MPI_INTEGER, &
          nErrors_all, 1, MPI_INTEGER, 0, iComm, iError)
     call MPI_gather(error_message, 60*maxerror, MPI_CHARACTER, &
          error_message_all, 60*maxerror, MPI_CHARACTER, 0, iComm,iError)
     call MPI_gather(error_count, maxerror, MPI_INTEGER, &
          error_count_all, maxerror, MPI_INTEGER, 0, iComm, iError)
     call MPI_gather(iter_first, maxerror, MPI_INTEGER, &
          iter_first_all, maxerror, MPI_INTEGER, 0, iComm, iError)
     call MPI_gather(iter_last, maxerror, MPI_INTEGER, &
          iter_last_all, maxerror, MPI_INTEGER, 0, iComm, iError)
     call MPI_gather(error_min, maxerror, MPI_REAL, &
          error_min_all, maxerror, MPI_REAL, 0, iComm, iError)
     call MPI_gather(error_max, maxerror, MPI_REAL, &
          error_max_all, maxerror, MPI_REAL, 0, iComm, iError)
     call MPI_gather(error_mean, maxerror, MPI_REAL, &
          error_mean_all, maxerror, MPI_REAL, 0, iComm, iError)
     call MPI_gather(error_last, maxerror, MPI_REAL, &
          error_last_all, maxerror, MPI_REAL, 0, iComm, iError)

     ! Analyze errors in PROC 0
     if(iProc==0)then
        nErrors=0
        do ip=1,nProc
           do i=1,nErrors_all(ip)
              msg=error_message_all(i,ip)
              i0=1
              do
                 if(i0>nErrors)then
                    nErrors=i0
                    error_message(i0)=msg
                    error_count_max(i0)=error_count_all(i,ip)
                    error_count_sum(i0)=error_count_all(i,ip)
                    iter_first(i0)=iter_first_all(i,ip)
                    iter_last(i0)=iter_last_all(i,ip)
                    error_min(i0)=error_min_all(i,ip)
                    error_max(i0)=error_max_all(i,ip)
                    error_mean(i0)=error_mean_all(i,ip)
                    error_last(i0)=error_last_all(i,ip)
                    error_last_sum(i0)=error_last_all(i,ip)
                    EXIT
                 end if
                 if(error_message(i0)==msg)then

                    error_mean(i0)=&
                         (error_mean_all(i,ip)*error_count_all(i,ip)+&
                         error_mean(i0)*error_count_sum(i0)) &
                         /(error_count_all(i,ip)+error_count_sum(i0))

                    if(iter_last(i0)<iter_last_all(i,ip))then
                       error_last(i0)=error_last_all(i,ip)
                       error_last_sum(i0)=error_last_all(i,ip)
                       iter_last(i0)=iter_last_all(i,ip)
                    elseif(iter_last(i0)==iter_last_all(i,ip))then
                       error_last_sum(i0)=error_last_sum(i0)+&
                            error_last_all(i,ip)
                    end if

                    error_count_sum(i0)=&
                         error_count_all(i,ip)+error_count_sum(i0)
                    error_count_max(i0)=&
                         max(error_count_all(i,ip),error_count_max(i0))
                    iter_first(i0)=&
                         min(iter_first_all(i,ip),iter_first(i0))
                    error_min(i0)=&
                         min(error_min_all(i,ip),error_min(i0))
                    error_max(i0)=&
                         max(error_max_all(i,ip),error_min(i0))
                    EXIT
                 end if
                 i0=i0+1
              end do ! i0
           end do ! error types
        end do ! processors

        ! Report errors
        if(nErrors==0)then
           call write_myname; write(*,*)'error report: no errors...'
        else
           do i=1,nErrors
              call write_myname
              write(*,'(a,a)')'Error_report for ',trim(error_message(i))
              call write_myname
              write(*,*)'OCCURED first=',iter_first(i),&
                   ' last=',iter_last(i),&
                   ' count_max=',error_count_max(i),&
                   ' count_sum=',error_count_sum(i)
              call write_myname
              write(*,*)'VALUES min=',error_min(i),' max=',error_max(i),&
                   ' mean=',error_mean(i),' last=',error_last(i),&
                   ' last_sum=',error_last_sum(i)
              call write_myname; write(*,*)
           end do
        end if

     end if ! iProc==0

     deallocate(nErrors_all,error_message_all,error_count_all,&
          iter_first_all,iter_last_all,error_min_all,error_max_all,&
          error_mean_all,error_last_all)

     RETURN

  end if ! PRINT

  if(iErrorIn<1 .or. iErrorIn>nErrors) then
     ! Determine iErrorIn based on str
     iErrorIn=1
     do
        if(iErrorIn>nErrors)then
           ! First occurance of this error type
           nErrors=iErrorIn
           EXIT
        end if
        if(error_message(iErrorIn)==str)EXIT
        iErrorIn=iErrorIn+1
     end do
  end if

  i=iErrorIn

  error_count(i)=error_count(i)+1

  iter_last(i)=iteration_number
  error_last(i)=value

  if(error_count(i)==1)then
     if(show_first)then
        call write_myname;
        write(*,*)'First error for ',str,' (PE=',iProc,&
          ') at iter=',iteration_number,' with value=',value
     end if
     error_message(i)=str
     iter_first(i)=iteration_number
     error_min(i)=value
     error_max(i)=value
     error_mean(i)=value
  else
     error_min(i)=min(error_min(i),value)
     error_max(i)=max(error_max(i),value)
     error_mean(i)=(error_mean(i)*(error_count(i)-1)+value)/error_count(i)
  end if

end subroutine error_report
!==============================================================================

subroutine test_error_report

  use ModProcMH
  use ModMain, ONLY : iteration_number
  implicit none

  integer:: ierr1=-1, ierr2=-1, ierr3=-1

  ! Test error_report
  !----------------------------------------------------------------------------
  select case(iProc)
  case(0)
     iteration_number=1
     call error_report('negative pressure',-1.,ierr1,.true.)
     call error_report('negative pressure',-2.,ierr1,.true.)
     call error_report('energy correction',0.1,ierr2,.false.)
     iteration_number=2
     call error_report('energy correction',0.2,ierr2,.false.)
     iteration_number=3
     call error_report('negative pressure',-6.,ierr1,.true.)
     call error_report('only PE 0',100.,ierr3,.true.)
     call error_report('energy correction',0.6,ierr2,.false.)
  case(1)
     iteration_number=1
     call error_report('only PE 1',200.,ierr3,.true.)
     call error_report('energy correction',0.01,ierr2,.false.)
     iteration_number=2
     call error_report('energy correction',0.02,ierr2,.false.)
     iteration_number=3
     call error_report('energy correction',0.06,ierr2,.false.)
  end select

  call error_report('PRINT',0.,iErr1,.true.)

end subroutine test_error_report
!==============================================================================

subroutine get_date_time_start(iTime_I)

  use ModMain,        ONLY : StartTime
  use ModTimeConvert, ONLY : time_real_to_int

  implicit none
  integer, intent(out) :: iTime_I(7)

  !----------------------------------------------------------------------------
  call time_real_to_int(StartTime,iTime_I)

end subroutine get_date_time_start
!==============================================================================

subroutine get_date_time(iTime_I)

  use ModMain,        ONLY : StartTime, Time_Simulation
  use ModTimeConvert, ONLY : time_real_to_int

  implicit none
  integer, intent(out) :: iTime_I(7)

  !----------------------------------------------------------------------------
  call time_real_to_int(StartTime+Time_Simulation,iTime_I)

end subroutine get_date_time
!==============================================================================

subroutine get_time_string

  use ModIO,   ONLY: StringDateOrTime, NameMaxTimeUnit
  use ModMain, ONLY: StartTime, Time_Simulation
  use ModTimeConvert, ONLY: TimeType, time_real_to_int
  implicit none

  integer:: i
  type(TimeType):: Time

  ! This is the value if the time is too large
  !----------------------------------------------------------------------------
  StringDateOrTime = '99999999'
  select case(NameMaxTimeUnit)
  case('hour')
     if(Time_Simulation < 10000.0*3600) &
          write(StringDateOrTime,'(i4.4,i2.2,i2.2)') &
          int(                            Time_Simulation/3600.), &
          int((Time_Simulation-(3600.*int(Time_Simulation/3600.)))/60.), &
          int( Time_Simulation-(  60.*int(Time_Simulation/  60.)))
  case('hr')
     if(Time_Simulation < 100.0*3600) &
          write(StringDateOrTime,'(i2.2,i2.2,f4.1)') &
          int(                            Time_Simulation/3600.), &
          int((Time_Simulation-(3600.*int(Time_Simulation/3600.)))/60.), &
          Time_Simulation-(  60.*int(Time_Simulation/  60.))
  case('minute')
     if(Time_Simulation < 100.0*60) &
          write(StringDateOrTime,'(i2.2,f6.3)') &
          int(Time_Simulation/60.), &
          Time_Simulation-(60.*int(Time_Simulation/60.))
  case('second')
     if(Time_Simulation < 100.0) &
          write(StringDateOrTime,'(f8.5)') Time_Simulation
  case('millisecond')
     if(Time_Simulation < 1.0) &
          write(StringDateOrTime,'(f8.4)') Time_Simulation*1e3
  case('microsecond')
     if(Time_Simulation < 1e-3) &
          write(StringDateOrTime,'(f8.4)') Time_Simulation*1e6
  case('nanosecond')
     if(Time_Simulation < 1e-6) &
          write(StringDateOrTime,'(f8.4)') Time_Simulation*1e9
  case default
     ! Could not find unit
     StringDateOrTime = ''
  end select

  if(StringDateOrTime /= '')then
     ! The time tag is 8-character long for the above cases
     ! Replaces spaces with 0-s up to 8 characters
     do i=1,8
        if(StringDateOrTime(i:i)==' ') StringDateOrTime(i:i)='0'
     end do
     RETURN
  end if

  ! Convert current date and time into string Time % String
  Time % Time = StartTime + Time_Simulation
  call time_real_to_int(Time)

  ! Select part of the string
  select case(NameMaxTimeUnit)
  case('date')
     StringDateOrTime = Time % String(1:14)
  case('year')
     StringDateOrTime = Time % String(1:8)
  case('yr')
     StringDateOrTime = Time % String(3:10)
  case('month')
     StringDateOrTime = Time % String(5:12)
  case('day')
     StringDateOrTime = Time % String(7:14)
  case default
     ! the unit is wrong, but what can we do? Let's write something.
     StringDateOrTime = '00000000'
  end select

end subroutine get_time_string
!==============================================================================

subroutine get_iVar(NameVar, iVar)

  use ModMain,       ONLY: NameVarLower_V
  use ModVarIndexes, ONLY: NameFluid_I, nVar
  use ModMultiFluid, ONLY: extract_fluid_name, iFluid
  use ModUtilities,  ONLY: lower_case

  character(len=10), intent(inout) :: NameVar
  integer,intent(out)              :: iVar

  integer :: iVarLoop, iError

  ! Initialize iVar
  character(len=*), parameter:: NameSub = 'get_iVar'
  !----------------------------------------------------------------------------
  iVar = -1

  ! In case NameVar is in upper case...
  call lower_case(NameVar)

  ! Remove the fluid name first
  call extract_fluid_name(NameVar)

  ! In case the user specifies ux/uy/uz instead of mx/my/mz
  select case(NameVar)
  case('ux')
     NameVar = 'mx'
  case('uy')
     NameVar = 'my'
  case('uz')
     NameVar = 'mz'
  end select

  ! Put back the fluid name
  if (iFluid /= 1) NameVar=trim(NameFluid_I(iFluid))//trim(NameVar)

  ! The first character in NameFluid_I is in upper case...
  call lower_case(NameVar)

  ! Find NameVar in NameVarLower_V if it is there
  do iVarLoop =1,nVar
     if (NameVar /= NameVarLower_V(iVarLoop)) CYCLE
     iVar = iVarLoop
     EXIT
  end do

  if(iVar < 0) then
     ! Try reading iVar as an index.
     read(NameVar,*,iostat=iError) iVar

     ! If it is not an integer and not one of the element in NameVarLower_V
     if (iError /= 0) call stop_mpi(NameSub//': unknown NameVar =' &
          //trim(NameVar))
  end if

  ! iVar must be within 1 and nVar
  if (iVar < 0 .or. iVar > nVar) call stop_mpi(NameSub//': check NameVar, ' &
       //'iVar is not within 1 and nVar???')

end subroutine get_iVar
!==============================================================================
