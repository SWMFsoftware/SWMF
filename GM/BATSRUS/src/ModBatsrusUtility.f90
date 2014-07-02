!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!==== Simple subroutines and functions that operate on all used blocks ========
!=============================================================================
real function integrate_BLK(qnum,qa)              

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
  integer :: iBLK, iError

  logical :: oktest, oktest_me
  !---------------------------------------------------------------------------
  call set_oktest('integrate_BLK',oktest, oktest_me)

  qsum=0.0

  if(IsCartesian)then
     do iBLK = 1, nBlock
        if(Unused_B(iBLK)) CYCLE
        if(true_BLK(iBLK)) then
           qsum = qsum + CellVolume_B(iBLK)* &
                sum(qa(1:nI,1:nJ,1:nK,iBLK))
        else
           qsum = qsum + CellVolume_B(iBLK)* &
                sum( qa(1:nI,1:nJ,1:nK,iBLK), &
                MASK=true_cell(1:nI,1:nJ,1:nK,iBLK))
        end if
     end do
  else
     do iBLK = 1, nBlock
        if(Unused_B(iBLK)) CYCLE
        if(true_BLK(iBLK)) then
           qsum=qsum + sum(CellVolume_GB(1:nI,1:nJ,1:nK,iBLK)* &
                qa(1:nI,1:nJ,1:nK,iBLK))
        else
           qsum=qsum + sum(CellVolume_GB(1:nI,1:nJ,1:nK,iBLK)* &
                qa(1:nI,1:nJ,1:nK,iBLK), &
                MASK=true_cell(1:nI,1:nJ,1:nK,iBLK))
        end if
     end do
  end if

  if(qnum>1)then
     call MPI_allreduce(qsum, qsum_all, 1,  MPI_REAL, MPI_SUM, &
          iComm, iError)
     integrate_BLK = qsum_all
     if(oktest)write(*,*)'me,sum,sum_all:',iProc,qsum,qsum_all
  else
     integrate_BLK = qsum
     if(oktest)write(*,*)'me,qsum:',iProc,qsum
  end if

end function integrate_BLK    

!=============================================================================

real function minval_BLK(qnum,qa)

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
  integer :: iBLK, iError

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('minval_BLK',oktest, oktest_me)

  qminval=1.e+30

  do iBLK=1,nBlock
     if(.not.Unused_B(iBLK)) then
        if(true_BLK(iBLK)) then
           qminval=min(qminval, minval(qa(1:nI,1:nJ,1:nK,iBLK)))
        else
           qminval=min(qminval, minval(qa(1:nI,1:nJ,1:nK,iBLK),&
                MASK=true_cell(1:nI,1:nJ,1:nK,iBLK)))
        endif
     end if
  end do

  if(qnum>1)then
     call MPI_allreduce(qminval, qminval_all, 1,  MPI_REAL, MPI_MIN, &
          iComm, iError)
     minval_BLK=qminval_all
     if(oktest)write(*,*)'me,minval,minval_all:',iProc,qminval,qminval_all
!     if(qminval_all>=1.e+30)call stop_mpi('Error in minval_BLK: NaN-s?')
  else
     minval_BLK=qminval
     if(oktest)write(*,*)'me,qminval:',iProc,qminval
!     if(qminval>=1.e+30)call stop_mpi('Error in minval_BLK: NaN-s?')
  endif

end function minval_BLK

!=============================================================================

real function maxval_BLK(qnum,qa)

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
  integer :: iBLK, iError

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('maxval_BLK',oktest, oktest_me)

  qmaxval=-1.e+30

  do iBLK=1,nBlock
     if(.not.Unused_B(iBLK)) then
        if(true_BLK(iBLK)) then
           qmaxval=max(qmaxval, maxval(qa(1:nI,1:nJ,1:nK,iBLK)))
        else
           qmaxval=max(qmaxval, maxval(qa(1:nI,1:nJ,1:nK,iBLK),&
                MASK=true_cell(1:nI,1:nJ,1:nK,iBLK)))
        endif
     end if
  end do

  if(qnum>1)then
     call MPI_allreduce(qmaxval, qmaxval_all, 1,  MPI_REAL, MPI_MAX, &
          iComm, iError)
     maxval_BLK=qmaxval_all
     if(oktest)write(*,*)'me,maxval,maxval_all:',iProc,qmaxval,qmaxval_all
     if(qmaxval_all<=-1.e+30)call stop_mpi('Error in maxval_BLK: NaN-s?')
  else
     maxval_BLK=qmaxval
     if(oktest)write(*,*)'qmaxval:',qmaxval
     if(qmaxval<=-1.e+30)call stop_mpi('Error in maxval_BLK: NaN-s?')
  endif

end function maxval_BLK

!=============================================================================
real function maxval_loc_BLK(qnum,qa,loc)

  ! Return maxval(qa)) corresponding to all used blocks and 
  ! also return the location of the maximum value into loc(5)=I,J,K,IBLK,PE
  ! If qnum<=1, return maxval for the processor, otherwise
  ! return the maximum for all processors.

  use ModProcMH
  use ModSize, ONLY: &
       nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nBlock, MaxBlock
  use ModMain, ONLY: Unused_B
  use ModGeometry, ONLY :true_cell
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum

  real, intent(in) :: qa(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

  integer, intent(out):: loc(5)

  ! Local variables:
  real    :: qmaxval, qmaxval_all
  integer :: i,j,k,iBLK, iError

  real, external :: maxval_BLK

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('maxval_loc_BLK',oktest, oktest_me)

  qmaxval=maxval_BLK(1,qa)
  if(qnum==1)then
     qmaxval_all=qmaxval
  else
     call MPI_allreduce(qmaxval, qmaxval_all, 1,  MPI_REAL, MPI_MAX, &
          iComm, iError)
  end if

  loc=-1
  if (qmaxval == qmaxval_all) then
     BLKLOOP: do iBLK=1,nBlock
        if(Unused_B(iBLK)) CYCLE
        do k=1,nK; do j=1,nJ; do i=1,nI
           if(.not.true_cell(i,j,k,iBLK)) CYCLE
           if(qa(i,j,k,iBLK)==qmaxval)then
              loc(1)=i; loc(2)=j; loc(3)=k; loc(4)=iBLK; loc(5)=iProc
              EXIT BLKLOOP
           end if
        enddo; enddo; enddo; 
     enddo BLKLOOP
  end if

  maxval_loc_BLK=qmaxval_all

end function maxval_loc_BLK
!=============================================================================
real function maxval_loc_abs_BLK(qnum,qa,loc)

  ! Return maxval(abs(qa)) corresponding to all used blocks and 
  ! also return the location of the maximum value into loc(5)=I,J,K,IBLK,PE
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
  integer :: i,j,k,iBLK, iError

  real, external :: maxval_abs_BLK

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('maxval_loc_abs_BLK',oktest, oktest_me)

  qmaxval=maxval_abs_BLK(1,qa)
  if(qnum==1)then
     qmaxval_all=qmaxval
  else
     call MPI_allreduce(qmaxval, qmaxval_all, 1,  MPI_REAL, MPI_MAX, &
          iComm, iError)
  end if

  loc=-1
  if (qmaxval == qmaxval_all) then
     BLKLOOP: do iBLK=1,nBlock
        if(Unused_B(iBLK)) CYCLE
        do k=1,nK; do j=1,nJ; do i=1,nI
           if(.not.true_cell(i,j,k,iBLK)) CYCLE
           if(abs(qa(i,j,k,iBLK))==qmaxval)then
              loc(1)=i; loc(2)=j; loc(3)=k; loc(4)=iBLK; loc(5)=iProc
              EXIT BLKLOOP
           end if
        enddo; enddo; enddo; 
     enddo BLKLOOP
  end if

  maxval_loc_abs_BLK=qmaxval_all

end function maxval_loc_abs_BLK
!=============================================================================
real function minval_loc_BLK(qnum,qa,loc)

  ! Return minval(qa)) corresponding to all used blocks and 
  ! also return the location of the minimum value into loc(5)=I,J,K,IBLK,PE
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
  integer :: i,j,k,iBLK, iError

  real, external :: minval_BLK

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('minval_loc_BLK',oktest, oktest_me)

  qminval=minval_BLK(1,qa)
  if(qnum==1)then
     qminval_all=qminval
  else
     call MPI_allreduce(qminval, qminval_all, 1,  MPI_REAL, MPI_MIN, &
          iComm, iError)
  end if

  loc=-1
  if (qminval == qminval_all) then
     BLKLOOP:do iBLK=1,nBlock
        if(Unused_B(iBLK)) CYCLE
        do k=1,nK; do j=1,nJ; do i=1,nI
           if(.not.true_cell(i,j,k,iBLK)) CYCLE
           if(qa(i,j,k,iBLK)==qminval)then
              loc(1)=i; loc(2)=j; loc(3)=k; loc(4)=iBLK; loc(5)=iProc
              EXIT BLKLOOP
           end if
        enddo; enddo; enddo; 
     enddo BLKLOOP
  end if

  minval_loc_BLK=qminval_all

end function minval_loc_BLK

!=============================================================================

real function maxval_abs_BLK(qnum,qa)

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
  integer :: iBLK, iError

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('maxval_abs_BLK',oktest, oktest_me)

  qmaxval=-1.0

  do iBLK=1,nBlock
     if(.not.Unused_B(iBLK)) then
        if(true_BLK(iBLK)) then
           qmaxval=max(qmaxval, maxval(abs(qa(1:nI,1:nJ,1:nK,iBLK))))
        else
           qmaxval=max(qmaxval, maxval(abs(qa(1:nI,1:nJ,1:nK,iBLK)),&
                MASK=true_cell(1:nI,1:nJ,1:nK,iBLK)))
        endif
     end if
  end do

  if(qnum>1)then
     call MPI_allreduce(qmaxval, qmaxval_all, 1,  MPI_REAL, MPI_MAX, &
          iComm, iError)
     maxval_abs_BLK=qmaxval_all
     if(oktest)write(*,*)'me,maxval,maxval_all:',iProc,qmaxval,qmaxval_all
     if(qmaxval_all<0.0) then
        call barrier_mpi
        call stop_mpi('Error in maxval_abs_BLK: negative max!')
     end if
  else
     maxval_abs_BLK=qmaxval
     if(oktest)write(*,*)'qmaxval:',qmaxval
     if(qmaxval<0.0)call stop_mpi('Error in maxval_abs_BLK: negative max!')
  endif

end function maxval_abs_BLK

!=============================================================================

real function maxval_abs_ALL(qnum,qa)

  ! Return maxval(abs(qa)) corresponding to all used blocks
  ! Include ghost cells and .not.true_cell -s too.
  ! If qnum<=1, return maxval(abs(qa)) for the processor, otherwise
  ! return the minimum for all processors.

  use ModProcMH
  use ModSize, ONLY: &
       nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nBlock, MaxBlock
  use ModMain, ONLY: Unused_B
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum

  real, intent(in):: qa(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

  ! Local variables:
  real    :: qmaxval, qmaxval_all
  integer :: iBLK, iError

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('maxval_abs_ALL',oktest, oktest_me)

  qmaxval=-1.0

  do iBLK=1,nBlock
     if(.not.Unused_B(iBLK)) then
        qmaxval=max(qmaxval, maxval(abs(qa(:,:,:,iBLK))))
     end if
  end do

  if(qnum>1)then
     call MPI_allreduce(qmaxval, qmaxval_all, 1,  MPI_REAL, MPI_MAX, &
          iComm, iError)
     maxval_abs_ALL=qmaxval_all
     if(oktest)write(*,*)'me,maxval,maxval_all:',iProc,qmaxval,qmaxval_all
  else
     maxval_abs_ALL=qmaxval
     if(oktest)write(*,*)'qmaxval:',qmaxval
  endif

end function maxval_abs_ALL

!=============================================================================
real function test_cell_value(qa,Imin,Imax,Jmin,Jmax,Kmin,Kmax)

  ! Find the value at the test cell location Itest, Jtest, Ktest,BLKtest
  ! PROCtest and then broadcast it to all PROC.

  use ModProcMH
  use ModMain, ONLY : nBLK,PROCtest,Itest,Jtest,Ktest,BLKtest
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: Imin,Imax,Jmin,Jmax,Kmin,Kmax
  real, dimension(Imin:Imax,Jmin:Jmax,Kmin:Kmax,nBLK), &
       intent(in) :: qa

  ! Local variables:
  integer :: iError
  real    :: qval

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('test_cell_value',oktest, oktest_me)

  qval=0.0

  if (PROCtest == iProc) qval = qa(Itest,Jtest,Ktest,BLKtest)

  call MPI_Bcast(qval,1,MPI_REAL,PROCtest,iComm,iError)

  test_cell_value=qval

  if(oktest)write(*,*)'i,j,k,BLK,PROC,qval:',Itest,Jtest,Ktest,BLKtest, &
       PROCtest,qval

end function test_cell_value

!=============================================================================

subroutine xyz_to_spherical(x, y, z, r, Phi, Colatitude)

  use ModNumConst, ONLY: cTwoPi
  implicit none

  real, intent(in)  :: x, y, z
  real, intent(out) :: r, Phi, Colatitude
  !---------------------------------------------------------------------------
  r = sqrt(x**2 + y**2 + z**2)

  if(r == 0.0)then
     Phi      = 0.0
     Colatitude = 0.0
     RETURN
  end if

  ! get the phi(longitude relative to +x) and colatitude
  if (x == 0.0 .and. y == 0.0) then
     Phi = 0.0
  else
     Phi = modulo(atan2(y, x), cTwoPi)
  end if
  Colatitude = acos(z/R)

end subroutine xyz_to_spherical

!=============================================================================

subroutine set_oktest(str,oktest,oktest_me)

  use ModProcMH
  use ModMain, ONLY : iteration_number,Ttest,iterTEST,PROCtest,lVerbose, &
       time_accurate,time_simulation,test_string
  implicit none

  character (len=*) :: str
  integer, external :: index_mine
  logical :: oktest, oktest_me
  !----------------------------------------------------------------------------

  if(iteration_number>iterTEST .or. &
       (time_accurate .and. time_simulation>Ttest))then
     oktest=index_mine(' '//test_string,' '//str//' ')>0
     oktest_me = oktest .and. iProc==PROCtest
     if(oktest_me)then
        write(*,*)str,' at iter=',iteration_number
     else if(lVerbose>=100)then
        write(*,*)str,' CALLED by me=',iProc,' at iter=',iteration_number
     else if(iProc==PROCtest.and.lVerbose>=10)then
        write(*,*)str,' CALLED at iter=',iteration_number
     endif
  else
     oktest    = .false.
     oktest_me = .false.
  end if

end subroutine set_oktest

!=============================================================================

integer function index_mine(str1,str2)

  implicit none

  character (len=*), intent(in) :: str1, str2

  index_mine=index(str1,str2)

end function index_mine

!=============================================================================

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

!=============================================================================

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

!=============================================================================

subroutine stop_mpi(str)

  use ModProcMH
  use ModMain, ONLY : iteration_number,NameThisComp,IsStandAlone
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
     call MPI_abort(iComm, nError, iError)
     stop
  else
     write(*,*)NameThisComp,': stopping execution! me=',iProc,&
          ' at iteration=',iteration_number
     call CON_stop(NameThisComp//':'//str)
  end if

end subroutine stop_mpi

!============================================================================

subroutine alloc_check(iError,String)

  implicit none

  integer, intent(in) :: iError
  character (len=*), intent(in) :: String
  !----------------------------------------------------------------------------

  if (iError>0) call stop_mpi("Allocation error for "//String)

end subroutine alloc_check

!==========================================================================

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

  !--------------------------------------------------------------------------

  !Debug
  !write(*,*)'Error_report me, iErrorIn, value, str=',iProc,iErrorIn,value,str

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
                    exit
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
                    exit
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

     return

  end if ! PRINT

  if(iErrorIn<1 .or. iErrorIn>nErrors) then
     ! Determine iErrorIn based on str
     iErrorIn=1
     do
        if(iErrorIn>nErrors)then
           ! First occurance of this error type
           nErrors=iErrorIn
           exit
        end if
        if(error_message(iErrorIn)==str)exit
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
subroutine join_str(n, String_I, String)

  ! Append the n strings in the string array String_I to the end of String
  ! separated by spaces.
  implicit none

  integer, intent(in)            :: n
  character(len=*), intent(in)   :: String_I(n)
  character(len=*), intent(inout):: String

  integer :: i
  !---------------------------------------------------------------------------
  do i = 1, n
     String = trim(String) // ' ' // String_I(i)
  end do

end subroutine join_str
!==============================================================================
subroutine find_test_cell

  ! Find cell indices corresponding to Xtest, Ytest, Ztest coordinates
  ! or print out cell coordinates corresponding to Itest, Jtest, Ktest, ...

  use ModProcMH
  use ModMain
  use ModParallel, ONLY: NOBLK, neiLEV, neiPE, neiBLK
  use ModGeometry, ONLY: r_BLK
  use BATL_lib, ONLY: Xyz_DGB, CellSize_DB, CellFace_DFB, CellVolume_GB, &
       IsCartesian, MaxDim, find_grid_block
  use ModMpi

  implicit none

  logical :: DoBroadcast
  integer :: IjkTest_D(MaxDim), iDir, iError

  character(len=*), parameter:: NameSub = 'find_test_cell'
  !----------------------------------------------------------------------------

  DoBroadcast = .false.

  if(.not.coord_test)then
     if(iProc == ProcTest)then
        if(1 <= BLKtest .and. BLKtest <= nBlock)then
           if(Unused_B(BLKtest))then
              if(lVerbose>0) write(*,*)'Test cell is in an unused block'
           else
              XyzTestCell_D = Xyz_DGB(:,Itest,Jtest,Ktest,BLKtest)
              DoBroadcast = .true.
           end if
        else
           if(lVerbose>0) write(*,*)'BLKtest=',BLKtest,&
                ' is out of 1..nBlock=',nBlock
        end if
     end if
     call MPI_Bcast(DoBroadcast,1,MPI_LOGICAL,PROCtest,iComm,iError)
     if (.not. DoBroadcast) RETURN

  else   ! if a coord_test
     DoBroadcast = .true.

     call find_grid_block( (/ xTest, yTest, zTest /), &
          ProcTest, BlkTest, IjkTest_D)

     if(ProcTest < 0)then

        if(iProc==0)then
           write(*,*) NameSub,' xTest, yTest, zTest=', xTest, yTest, zTest
           write(*,*) NameSub,' WARNING test point was not found! Setting defaults.'
        end if

        iTest = 1
        jTest = 1
        kTest = 1
        BlkTest  = 1
        ProcTest = 0
        xTest = 0.0
        yTest = 0.0
        zTest = 0.0
        XyzTestCell_D = 0.0
     else        
        iTest = IjkTest_D(1)
        jTest = IjkTest_D(2)
        kTest = Ijktest_D(3)
        if(iProc==ProcTest) XyzTestCell_D = Xyz_DGB(:,Itest,Jtest,Ktest,BLKtest)
     end if

  end if
  if (DoBroadcast) then

     call MPI_Bcast(XyzTestCell_D, 3, MPI_REAL, ProcTest, iComm, iError)

     if(iProc == ProcTest .and. UseTestCell .and. lVerbose>0)then
        write(*,*)
        write(*,*)'Selected test cell:'
        write(*,'(a,i4,a,i4,a,i4,a,i8,a,i5)')&
             'I=',Itest,' J=',Jtest,' K=',Ktest,&
             ' BLK=',BLKtest,' PE=',PROCtest
        write(*,'(a,3f12.5,a,f12.5)') &
             'x,y,z=',Xyz_DGB(:,Itest,Jtest,Ktest,BLKtest), &
             ' r=',r_BLK(iTest,jTest,kTest,BLKtest)
        write(*,'(a,3f12.5,a,f12.5)') &
             ' CellSize_D=',CellSize_DB(:,BLKtest),&
             ' CellVolume=',CellVolume_GB(iTest,jTest,kTest,BLKtest)
        if(.not.IsCartesian) write(*,'(a,3f12.5)') &
             ' CellFace_D=',CellFace_DFB(:,iTest,jTest,kTest,BLKtest)
  	do iDir = 1, 6
  	   select case(neiLEV(iDir,BLKtest))
           case(0,1)
              write(*,'(a,i2,a,i2,a,i5,a,i8)')&
                   'idir=',idir,' neiLEV=',neiLEV(idir,BLKtest),&
                   ' neiPE=',neiPE(1,idir,BLKtest),&
                   ' neiBLK=',neiBLK(1,idir,BLKtest)
           case(-1)
              write(*,'(a,i2,a,i2,a,4i5,a,4i8)')&
                   'idir=',idir,' neiLEV=',neiLEV(idir,BLKtest),&
                   ' neiPE=',neiPE(:,idir,BLKtest),&
                   ' neiBLK=',neiBLK(:,idir,BLKtest)
  	   case(NOBLK)
  	      write(*,'(a,i2,a,i5)')&
                   'idir=',iDir,' neiLEV=',neiLEV(iDir,BLKtest)
  	   end select
  	end do
  	write(*,*)

     end if

  end if

end subroutine find_test_cell

!=============================================================================

subroutine xyz_to_peblk(x, y, z, iPe, iBlock, DoFindCell, i, j, k)

  ! Find the processor (iPe) and block (iBlock) for position x, y, z
  ! If DoFindCell is true then the cell indexes i, j, k are returned too.

  use BATL_lib, ONLY: MaxDim, find_grid_block

  implicit none

  real,    intent(in) :: x, y, z
  integer, intent(out):: iPE, iBlock
  logical, intent(in) :: DoFindCell
  integer, intent(out):: i, j, k

  ! Variables for BATL
  integer:: iCell_D(MaxDim)

  character(len=*), parameter:: NameSub = 'xyz_to_peblk'
  !----------------------------------------------------------------------

  call find_grid_block( (/x,y,z/), iPe, iBlock, iCell_D)
  i = iCell_D(1)
  j = iCell_D(2)
  k = iCell_D(3)

end subroutine xyz_to_peblk
!=============================================================================
subroutine fill_edge_corner(Array_G)

  ! Fill edges and corners for Array_G (e.g. for bilinear interpolation)

  use ModSize, ONLY: nI, nJ, nK

  implicit none
  real, intent(inout) :: Array_G(0:nI+1,0:nJ+1,0:nK+1)
  !---------------------------------------------------------------------------
  ! Edges in the` K direction
  Array_G(0,0,:) = &
       Array_G(1,0,:)     + Array_G(0,1,:)     - Array_G(1,1,:)
  Array_G(nI+1,0,:) = &
       Array_G(nI,0,:)    + Array_G(nI+1,1,:)  - Array_G(nI,1,:)
  Array_G(0,nJ+1,:) = &
       Array_G(0,nJ,:)    + Array_G(1,nJ+1,:)  - Array_G(1,nJ,:)
  Array_G(nI+1,nJ+1,:) = &
       Array_G(nI+1,nJ,:) + Array_G(nI,nJ+1,:) - Array_G(nI,nJ,:)
  
  ! Edges in the J direction
  Array_G(0,:,0) = &
       Array_G(1,:,0)     + Array_G(0,:,1)     - Array_G(1,:,1)
  Array_G(nI+1,:,0) = &
       Array_G(nI,:,0)    + Array_G(nI+1,:,1)  - Array_G(nI,:,1)
  Array_G(0,:,nK+1) = &
       Array_G(1,:,nK+1)  + Array_G(0,:,nK)    - Array_G(1,:,nK)
  Array_G(nI+1,:,nK+1) = &
       Array_G(nI,:,nK+1) + Array_G(nI+1,:,nK) - Array_G(nI,:,nK)
  
  ! Edges in the I direction including the corners
  Array_G(:,0,0) = &
       Array_G(:,1,0)     + Array_G(:,0,1)     - Array_G(:,1,1)
  Array_G(:,nJ+1,0) = &
       Array_G(:,nJ,0)    + Array_G(:,nJ+1,1)  - Array_G(:,nJ,1)
  Array_G(:,0,nK+1)=&
       Array_G(:,1,nK+1)  + Array_G(:,0,nK)    - Array_G(:,1,nK)
  Array_G(:,nJ+1,nK+1) = &
       Array_G(:,nJ,nK+1) + Array_G(:,nJ+1,nK) - Array_G(:,nJ,nK)

end subroutine fill_edge_corner

!=========================================================================
subroutine get_date_time_start(iTime_I)
  
  use ModMain,        ONLY : StartTime
  use ModTimeConvert, ONLY : time_real_to_int

  implicit none
  integer, intent(out) :: iTime_I(7)

  call time_real_to_int(StartTime,iTime_I)

end subroutine get_date_time_start

!=========================================================================
subroutine get_date_time(iTime_I)
  
  use ModMain,        ONLY : StartTime, Time_Simulation
  use ModTimeConvert, ONLY : time_real_to_int

  implicit none
  integer, intent(out) :: iTime_I(7)

  call time_real_to_int(StartTime+Time_Simulation,iTime_I)

end subroutine get_date_time

!=========================================================================
subroutine get_time_string

  use ModIO,   ONLY: StringDateOrTime, NameMaxTimeUnit
  use ModMain, ONLY: StartTime, Time_Simulation
  use ModTimeConvert, ONLY: TimeType, time_real_to_int
  implicit none

  integer:: i
  type(TimeType):: Time
  !---------------------------------------------------------------------------

  ! This is the value if the time is too large
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

