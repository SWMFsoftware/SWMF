!^CFG COPYRIGHT UM
!==== Simple subroutines and functions that operate on all used blocks ========
!^CFG IF PROJECTION BEGIN
subroutine set_BLK(qa,qb)

  ! Set qa=qb for all used blocks, where qb is a scalar

  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModGeometry, ONLY : true_BLK,true_cell
  implicit none

  ! Arguments

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(out) :: qa
  real, intent(in) :: qb

  ! Local variables:
  integer:: iBLK

  !---------------------------------------------------------------------------

  do iBLK=1,nBlockMax
     if(.not.unusedBLK(iBLK)) then
        if(true_BLK(iBLK))then
           qa(1:nI,1:nJ,1:nK,iBLK)=qb
        else
           where(true_cell(1:nI,1:nJ,1:nK,iBLK)) &
                qa(1:nI,1:nJ,1:nK,iBLK)=qb
        end if
     end if
  end do

end subroutine set_BLK

!=============================================================================
subroutine eq_BLK(qa,qb)

  ! Do qa=qb for all used blocks

  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModGeometry, ONLY : true_BLK,true_cell
  implicit none

  ! Arguments

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK):: qa,qb
  intent(out) :: qa
  intent(in)  :: qb

  ! Local variables:
  integer:: iBLK

  !---------------------------------------------------------------------------

  do iBLK=1,nBlockMax
     if(.not.unusedBLK(iBLK)) then
        if(true_BLK(iBLK))then
           qa(1:nI,1:nJ,1:nK,iBLK)= qb(1:nI,1:nJ,1:nK,iBLK)
        else
           where(true_cell(1:nI,1:nJ,1:nK,iBLK)) &
                qa(1:nI,1:nJ,1:nK,iBLK)= qb(1:nI,1:nJ,1:nK,iBLK)
        end if

     end if
  end do

end subroutine eq_BLK

!=============================================================================
subroutine add_BLK(qa,qb)

  ! Do qa=qa+qb for all used blocks

  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModGeometry, ONLY : true_BLK,true_cell
  implicit none

  ! Arguments

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK) :: qa, qb
  intent(inout) :: qa
  intent(in)    :: qb

  ! Local variables:
  integer:: iBLK

  !---------------------------------------------------------------------------

  do iBLK=1,nBlockMax
     if(.not.unusedBLK(iBLK)) then
        if(true_BLK(iBLK))then
           qa(1:nI,1:nJ,1:nK,iBLK)= &
                qa(1:nI,1:nJ,1:nK,iBLK)+qb(1:nI,1:nJ,1:nK,iBLK)
        else
           where(true_cell(1:nI,1:nJ,1:nK,iBLK)) &
                qa(1:nI,1:nJ,1:nK,iBLK)= &
                qa(1:nI,1:nJ,1:nK,iBLK)+qb(1:nI,1:nJ,1:nK,iBLK)
        end if
     end if
  end do

end subroutine add_BLK

!=============================================================================
subroutine sub_BLK(qa,qb)

  ! Do qa=qa-qb for all used blocks

  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModGeometry, ONLY : true_BLK,true_cell
  implicit none

  ! Arguments

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK) :: qa, qb
  intent(inout) :: qa
  intent(in)    :: qb

  ! Local variables:
  integer:: iBLK

  !---------------------------------------------------------------------------

  do iBLK=1,nBlockMax
     if(.not.unusedBLK(iBLK)) then
        if(true_BLK(iBLK))then
           qa(1:nI,1:nJ,1:nK,iBLK)= &
                qa(1:nI,1:nJ,1:nK,iBLK)-qb(1:nI,1:nJ,1:nK,iBLK)
        else
           where(true_cell(1:nI,1:nJ,1:nK,iBLK)) &
                qa(1:nI,1:nJ,1:nK,iBLK)= &
                qa(1:nI,1:nJ,1:nK,iBLK)-qb(1:nI,1:nJ,1:nK,iBLK)
        end if
     end if
  end do

end subroutine sub_BLK

!=============================================================================
subroutine eq_plus_BLK(qa,qb,qc)

  ! Do qa=qb+qc for all used blocks

  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModGeometry, ONLY : true_BLK,true_cell
  implicit none

  ! Arguments

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK) :: qa,qb,qc
  intent(out) :: qa
  intent(in)  :: qb,qc

  ! Local variables:
  integer:: iBLK

  !---------------------------------------------------------------------------

  do iBLK=1,nBlockMax
     if(.not.unusedBLK(iBLK)) then
        if(true_BLK(iBLK))then
           qa(1:nI,1:nJ,1:nK,iBLK)= &
                qb(1:nI,1:nJ,1:nK,iBLK)+qc(1:nI,1:nJ,1:nK,iBLK)
        else
           where(true_cell(1:nI,1:nJ,1:nK,iBLK)) &
                qa(1:nI,1:nJ,1:nK,iBLK)= &
                qb(1:nI,1:nJ,1:nK,iBLK)+qc(1:nI,1:nJ,1:nK,iBLK)
        end if
     end if
  end do


end subroutine eq_plus_BLK

!=============================================================================
subroutine add_times_BLK(qa,qb,qc)

  ! Do qa=qa+qb*qc for all used blocks, where qb is a scalar

  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModGeometry, ONLY : true_BLK,true_cell
  implicit none

  ! Arguments

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK) :: qa,qc
  intent(inout) :: qa
  intent(in)    :: qc

  real, intent(in) :: qb

  ! Local variables:
  integer:: iBLK

  !---------------------------------------------------------------------------

  do iBLK=1,nBlockMax
     if(.not.unusedBLK(iBLK)) then
        if(true_BLK(iBLK))then
           qa(1:nI,1:nJ,1:nK,iBLK)= &
                qa(1:nI,1:nJ,1:nK,iBLK)+qb*qc(1:nI,1:nJ,1:nK,iBLK)
        else
           where(true_cell(1:nI,1:nJ,1:nK,iBLK)) &
                qa(1:nI,1:nJ,1:nK,iBLK)= &
                qa(1:nI,1:nJ,1:nK,iBLK)+qb*qc(1:nI,1:nJ,1:nK,iBLK)
        end if
     end if
  end do


end subroutine add_times_BLK

!=============================================================================
subroutine eq_plus_times_BLK(qa,qb,qc,qd)

  ! Do qa=qb+qc*qd for all used blocks, where qc is a scalar

  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModGeometry, ONLY : true_BLK,true_cell
  implicit none

  ! Arguments

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK) :: qa,qb,qd
  intent(inout) :: qa
  intent(in)    :: qb
  intent(inout) :: qd

  real, intent(in) :: qc

  ! Local variables:
  integer:: iBLK

  !---------------------------------------------------------------------------

  do iBLK=1,nBlockMax
     if(.not.unusedBLK(iBLK)) then
        if(true_BLK(iBLK))then
           qa(1:nI,1:nJ,1:nK,iBLK)= &
                qb(1:nI,1:nJ,1:nK,iBLK)+qc*qd(1:nI,1:nJ,1:nK,iBLK)
        else
           where(true_cell(1:nI,1:nJ,1:nK,iBLK)) &
                qa(1:nI,1:nJ,1:nK,iBLK)= &
                qb(1:nI,1:nJ,1:nK,iBLK)+qc*qd(1:nI,1:nJ,1:nK,iBLK)
        end if
     end if
  end do

end subroutine eq_plus_times_BLK

!=============================================================================
real function dot_product_BLK(qa,qb)

  ! Return qa.qb=sum(qa*qb) for all used blocks

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModGeometry, ONLY : true_BLK,true_cell
  use ModMpi
  implicit none

  ! Arguments

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(in) :: qa,qb

  ! Local variables:
  real    :: qproduct, qproduct_all
  integer :: iBLK, iError

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('dot_product_BLK',oktest, oktest_me)

  qproduct=0.0

  do iBLK=1,nBlockMax
     if(.not.unusedBLK(iBLK)) then
        if(true_BLK(iBLK)) then
           qproduct=qproduct + &
                sum(qa(1:nI,1:nJ,1:nK,iBLK)*qb(1:nI,1:nJ,1:nK,iBLK))
        else
           qproduct=qproduct + &
                sum(qa(1:nI,1:nJ,1:nK,iBLK)*qb(1:nI,1:nJ,1:nK,iBLK),&
                MASK=true_cell(1:nI,1:nJ,1:nK,iBLK))
        end if
     end if
  end do

  if(nProc>1)then
     call MPI_allreduce(qproduct, qproduct_all, 1,  MPI_REAL, MPI_SUM, &
          iComm, iError)
     dot_product_BLK=qproduct_all
     if(oktest)write(*,*)'me,product,product_all:',&
          iProc,qproduct,qproduct_all
  else
     dot_product_BLK=qproduct
     if(oktest)write(*,*)'me,qproduct:',iProc,qproduct
  end if

end function dot_product_BLK

!=============================================================================
real function sum_BLK(qnum,qa)

  ! Return sum(qa) for all used blocks and true cells
  ! Do for each processor separately if qnum=1, otherwise add them all

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModGeometry, ONLY : true_BLK,true_cell
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum
  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(in) :: qa

  ! Local variables:
  real    :: qsum, qsum_all
  integer :: iBLK, iError

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('sum_BLK',oktest, oktest_me)

  qsum=0.0

  do iBLK=1,nBlockMax
     if(.not.unusedBLK(iBLK)) then
        if(true_BLK(iBLK)) then
           qsum=qsum + sum(qa(1:nI,1:nJ,1:nK,iBLK))
        else
           qsum=qsum + sum(qa(1:nI,1:nJ,1:nK,iBLK), &
                MASK=true_cell(1:nI,1:nJ,1:nK,iBLK))
        end if
     end if
  end do

  if(qnum>1)then
     call MPI_allreduce(qsum, qsum_all, 1,  MPI_REAL, MPI_SUM, &
          iComm, iError)
     sum_BLK=qsum_all
     if(oktest)write(*,*)'me,sum,sum_all:',iProc,qsum,qsum_all
  else
     sum_BLK=qsum
     if(oktest)write(*,*)'me,qsum:',iProc,qsum
  end if

end function sum_BLK
!^CFG END PROJECTION
!=============================================================================
real function integrate_BLK(qnum,qa)              

  ! Return the volume integral of qa, ie. sum(qa*cV_BLK) 
  ! for all used blocks and true cells
  ! Do for each processor separately if qnum=1, otherwise add them all

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModGeometry, ONLY : cV_BLK,true_BLK,true_cell 
  use ModGeometry, ONLY : UseCovariant      
  use ModMpi
  implicit none 

  ! Arguments

  integer, intent(in) :: qnum
  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(in) :: qa

  ! Local variables:
  real    :: qsum, qsum_all
  integer :: iBLK, iError

  logical :: oktest, oktest_me
  real,external:: integrate_BLK_covar  
  !---------------------------------------------------------------------------
  if(UseCovariant)then
     integrate_BLK=integrate_BLK_covar(qnum,qa)
     return
  end if                               
  call set_oktest('integrate_BLK',oktest, oktest_me)

  qsum=0.0
                                                     
  do iBLK=1,nBlockMax
     if(.not.unusedBLK(iBLK)) then
        if(true_BLK(iBLK)) then
           qsum=qsum + sum(qa(1:nI,1:nJ,1:nK,iBLK)*&
                cV_BLK(iBLK))
        else
           qsum=qsum + sum(qa(1:nI,1:nJ,1:nK,iBLK)*&
                cV_BLK(iBLK), &
                MASK=true_cell(1:nI,1:nJ,1:nK,iBLK))
        end if
     end if
  end do
                                                    
  if(qnum>1)then
     call MPI_allreduce(qsum, qsum_all, 1,  MPI_REAL, MPI_SUM, &
          iComm, iError)
     integrate_BLK=qsum_all
     if(oktest)write(*,*)'me,sum,sum_all:',iProc,qsum,qsum_all
  else
     integrate_BLK=qsum
     if(oktest)write(*,*)'me,qsum:',iProc,qsum
  end if
end function integrate_BLK    



!=============================================================================

real function minval_BLK(qnum,qa)

  ! Return minval(qa)) corresponding to all used blocks
  ! If qnum<=1, return minval for the processor, otherwise
  ! return the minimum for all processors.

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModGeometry, ONLY : true_BLK,true_cell
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(in) :: qa

  ! Local variables:
  real    :: qminval, qminval_all
  integer :: iBLK, iError

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('minval_BLK',oktest, oktest_me)

  qminval=1.e+30

  do iBLK=1,nBlockMax
     if(.not.unusedBLK(iBLK)) then
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
!     if(qminval_all>=1.e+30)call stop_mpi('Error in minval_BLK: huge min!')
  else
     minval_BLK=qminval
     if(oktest)write(*,*)'me,qminval:',iProc,qminval
!     if(qminval>=1.e+30)call stop_mpi('Error in minval_BLK: huge min!')
  endif

end function minval_BLK

!=============================================================================

real function maxval_BLK(qnum,qa)

  ! Return maxval(qa)) corresponding to all used blocks
  ! If qnum<=1, return maxval for the processor, otherwise
  ! return the maximum for all processors.

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModGeometry, ONLY : true_BLK,true_cell
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(in) :: qa

  ! Local variables:
  real    :: qmaxval, qmaxval_all
  integer :: iBLK, iError

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('maxval_BLK',oktest, oktest_me)

  qmaxval=-1.e+30

  do iBLK=1,nBlockMax
     if(.not.unusedBLK(iBLK)) then
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
     if(qmaxval_all<=-1.e+30)call stop_mpi('Error in maxval_BLK: tiny max!')
  else
     maxval_BLK=qmaxval
     if(oktest)write(*,*)'qmaxval:',qmaxval
     if(qmaxval<=-1.e+30)call stop_mpi('Error in maxval_BLK: tiny max!')
  endif

end function maxval_BLK

!=============================================================================
real function maxval_loc_BLK(qnum,qa,loc)

  ! Return maxval(qa)) corresponding to all used blocks and 
  ! also return the location of the maximum value into loc(5)=I,J,K,IBLK,PE
  ! If qnum<=1, return maxval for the processor, otherwise
  ! return the maximum for all processors.

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModGeometry, ONLY : true_BLK,true_cell
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(in) :: qa

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
     BLKLOOP: do iBLK=1,nBlockMax
        if(unusedBLK(iBLK)) CYCLE
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
  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModGeometry, ONLY : true_cell
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(in) :: qa

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
     BLKLOOP: do iBLK=1,nBlockMax
        if(unusedBLK(iBLK)) CYCLE
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
  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModGeometry, ONLY : true_cell
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(in) :: qa

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
     BLKLOOP:do iBLK=1,nBlockMax
        if(unusedBLK(iBLK)) CYCLE
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
  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModGeometry, ONLY : true_cell,true_BLK
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(in) :: qa

  ! Local variables:
  real    :: qmaxval, qmaxval_all
  integer :: iBLK, iError

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('maxval_abs_BLK',oktest, oktest_me)

  qmaxval=-1.0

  do iBLK=1,nBlockMax
     if(.not.unusedBLK(iBLK)) then
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
  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum

  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), &
       intent(in) :: qa

  ! Local variables:
  real    :: qmaxval, qmaxval_all
  integer :: iBLK, iError

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('maxval_abs_ALL',oktest, oktest_me)

  qmaxval=-1.0

  do iBLK=1,nBlockMax
     if(.not.unusedBLK(iBLK)) then
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

subroutine xyz_to_spherical(x, y, z, r, Phi, Latitude)

  use ModNumConst, ONLY: cTwoPi
  implicit none

  real, intent(in)  :: x, y, z
  real, intent(out) :: r, Phi, Latitude
  !---------------------------------------------------------------------------
  r = sqrt(x**2 + y**2 + z**2)

  if(r == 0.0)then
     Phi      = 0.0
     Latitude = 0.0
     RETURN
  end if

  ! get the phi(longitude relative to +x) and latitude (relative to equator)
  if (x == 0.0 .and. y == 0.0) then
     Phi = 0.0
  else
     Phi = modulo(atan2(y, x), cTwoPi)
  end if
  Latitude = acos(z/R)

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
  integer, parameter :: maxerror=100

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
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,r_BLK,dx_BLK
  use ModParallel, ONLY : NOBLK, neiLEV,neiPE,neiBLK
  use ModAdvance,  ONLY : tmp1_BLK
  use ModMpi
  implicit none

  real :: qdist, qdist_min
  logical :: pass_message
  integer :: loc(5), idir, iError, iProcTestMe
  real, external :: minval_loc_BLK
  !----------------------------------------------------------------------------

  pass_message = .false.

  if(.not.coord_test)then
     if(iProc==PROCtest)then
        if(1<=BLKtest.and.BLKtest<=nBlockMax)then
           if(unusedBLK(BLKtest))then
              write(*,*)'Test cell is in an unused block'
           else
              Xtest_mod = x_BLK(Itest,Jtest,Ktest,BLKtest)
              Ytest_mod = y_BLK(Itest,Jtest,Ktest,BLKtest)
              Ztest_mod = z_BLK(Itest,Jtest,Ktest,BLKtest)
              pass_message = .true.
           end if
        else
           write(*,*)'BLKtest=',BLKtest,' is out of 1..nBlockMax=',&
                nBlockMax
        end if
     end if
     call MPI_Bcast(pass_message,1,MPI_LOGICAL,PROCtest,iComm,iError)
     if (.not. pass_message) return

  else   ! if a coord_test

     pass_message = .true.

     tmp1_BLK(1:nI,1:nJ,1:nK,1:nBlockMax)=&
          abs(x_BLK(1:nI,1:nJ,1:nK,1:nBlockMax)-Xtest)+&
          abs(y_BLK(1:nI,1:nJ,1:nK,1:nBlockMax)-Ytest)+&
          abs(z_BLK(1:nI,1:nJ,1:nK,1:nBlockMax)-Ztest)

     qdist=minval_loc_BLK(nProc,tmp1_BLK,loc)

     !!! write(*,*)'minval=',qdist,' loc=',loc

     Itest=loc(1)
     Jtest=loc(2)
     Ktest=loc(3)
     BLKtest=loc(4)
     iProcTestMe=loc(5)

     ! Tell everyone which processor contains the test cell
     ! The others have -1 so MPI_MAX behaves like a broadcast.
     call MPI_allreduce(iProcTestMe,PROCtest,1,MPI_INTEGER,MPI_MAX,&
          iComm,iError)

     if(iProc==ProcTest)then
        Xtest_mod = x_BLK(Itest,Jtest,Ktest,BLKtest)
        Ytest_mod = y_BLK(Itest,Jtest,Ktest,BLKtest)
        Ztest_mod = z_BLK(Itest,Jtest,Ktest,BLKtest)
     end if
  end if

  if (pass_message) then

     call MPI_Bcast(Xtest_mod,1,MPI_REAL,PROCtest,iComm,iError)
     call MPI_Bcast(Ytest_mod,1,MPI_REAL,PROCtest,iComm,iError)
     call MPI_Bcast(Ztest_mod,1,MPI_REAL,PROCtest,iComm,iError)

     if(iProc==PROCtest .and. UseTestCell .and. lVerbose>0)then
        write(*,*)
        write(*,*)'Selected test cell:'
        write(*,'(a,i4,a,i4,a,i4,a,i8,a,i5)')&
             'I=',Itest,' J=',Jtest,' K=',Ktest,&
             ' BLK=',BLKtest,' PE=',PROCtest
        write(*,'(a,f12.5,a,f12.5,a,f12.5,a,f12.5)') &
             'x=',x_BLK(Itest,Jtest,Ktest,BLKtest),&
             ' y=',y_BLK(Itest,Jtest,Ktest,BLKtest),&
             ' z=',z_BLK(Itest,Jtest,Ktest,BLKtest),&
             ' r=',r_BLK(iTest,jTest,kTest,BLKtest),&
             ' dx=',dx_BLK(BLKtest)

  	do idir=1,6
  	   select case(neiLEV(idir,BLKtest))
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
                   'idir=',idir,' neiLEV=',neiLEV(idir,BLKtest)
  	   end select
  	end do
  	write(*,*)

     end if

  end if

end subroutine find_test_cell

!=============================================================================

subroutine xyz_to_peblk(x,y,z,iPe,iBlock,DoFindCell,iCell,jCell,kCell)

  ! The programm returns the value of iPE and iBlock for
  ! the given Xyz values. If DoFindIjk=.true., 
  ! the i,j,k values are returned too

  use ModParallel,ONLY : proc_dims
  use ModOctree, ONLY: adaptive_block_ptr, octree_roots
  use ModSize, ONLY: nIJK_D
  use ModGeometry, ONLY : UseCovariant         
  use ModGeometry, ONLY : XyzMin_D, XyzMax_D
  use ModNumConst
  implicit none

  real, intent(in) :: x,y,z
  integer, intent(out) :: iPE,iBlock
  logical, intent(in) :: DoFindCell
  integer, intent(out):: iCell,jCell,kCell

  type(adaptive_block_ptr):: Octree
  real,dimension(3) :: Xyz_D,DXyz_D,XyzCorner_D,XyzCenter_D
  integer,dimension(3)::IjkRoot_D
  logical,dimension(3):: IsLowerThanCenter_D
  character(len=*), parameter:: NameSub = 'xyz_to_peblk'
  !----------------------------------------------------------------------

  nullify(Octree % ptr)

  ! Perform the coordinate transformation, if needed
                     
  if(UseCovariant)then                       
     call xyz_to_gen((/x,y,z/),Xyz_D)
  else                                       
     Xyz_D = (/x,y,z/)
  end if                                     

  !Check if we are within the domain:
  if(  any(Xyz_D(1:3) < XyzMin_D(1:3)-cTiny) .or. &
       any(Xyz_D(1:3) > XyzMax_D(1:3)+cTiny) )then
     write(*,*)NameSub,' Xyz_D   =',Xyz_D
     write(*,*)NameSub,' XyzMin_D=',XyzMin_D
     write(*,*)NameSub,' XyzMax_D=',XyzMax_D
     write(*,*)NameSub,' x,y,z   =',x,y,z
     if(UseCovariant)write(*,*)NameSub,' UseCovariant is TRUE'
     call stop_mpi(NameSub//': the point is out of the domain')
  end if

  !Find the octree root
  Dxyz_D      = (XyzMax_D-XyzMin_D)/proc_dims
  IjkRoot_D   = int((Xyz_D-XyzMin_D)/DXyz_D)

  ! Make sure that we remain within the index range
  ! The int function takes care of values below XyzMin_D
  ! We fix indices that are too large if any Xyz_D >= XyzMax_D
  IjkRoot_D   = min(IjkRoot_D, proc_dims-1)

  XyzCorner_D = XyzMin_D + DXyz_D*IjkRoot_D

  Octree % ptr => &
       octree_roots(IjkRoot_D(1)+1,IjkRoot_D(2)+1,IjkRoot_D(3)+1) % ptr

  ! Descend the octree to find the block containing the point
  do
     if(Octree % ptr % used) then
        iPE    = octree % ptr % PE
        iBlock = octree % ptr % BLK
        if(DoFindCell)then
           DXyz_D = DXyz_D/nIJK_D
           iCell  = int((Xyz_D(1)-XyzCorner_D(1))/DXyz_D(1))+1
           jCell  = int((Xyz_D(2)-XyzCorner_D(2))/DXyz_D(2))+1
           kCell  = int((Xyz_D(3)-XyzCorner_D(3))/DXyz_D(3))+1
        end if
        EXIT
     else
        DXyz_D = 0.5*DXyz_D
        XyzCenter_D = XyzCorner_D + DXyz_D
        IsLowerThanCenter_D = Xyz_D < XyzCenter_D
        if(IsLowerThanCenter_D(2))then
           if(.not.IsLowerThanCenter_D(3))then
              XyzCorner_D(3) = XyzCenter_D(3)
              if(IsLowerThanCenter_D(1))then
                 Octree % ptr => Octree % ptr % child(1)%ptr
              else
                 XyzCorner_D(1) = XyzCenter_D(1)
                 Octree % ptr => Octree % ptr % child(2)%ptr
              end if
           else
              if(.not.IsLowerThanCenter_D(1))then
                 XyzCorner_D(1) = XyzCenter_D(1)
                 Octree % ptr => Octree % ptr % child(3)%ptr
              else
                 Octree % ptr => Octree % ptr % child(4)%ptr
              end if
           end if
        else
           XyzCorner_D(2) = XyzCenter_D(2)
           if(IsLowerThanCenter_D(3))then
              if(IsLowerThanCenter_D(1))then
                 Octree % ptr => Octree % ptr % child(5)%ptr
              else
                 XyzCorner_D(1)=XyzCenter_D(1)
                 Octree % ptr => Octree % ptr % child(6)%ptr
              end if
           else
              XyzCorner_D(3) = XyzCenter_D(3)
              if(.not.IsLowerThanCenter_D(1))then
                 XyzCorner_D(1)=XyzCenter_D(1)
                 Octree % ptr => Octree % ptr % child(7)%ptr
              else
                 Octree % ptr => Octree % ptr % child(8)%ptr
              end if
           end if
        end if
     end if
  end do
end subroutine xyz_to_peblk
!=============================================================================
subroutine fill_edge_corner(Array_G)

  ! Fill edges and corners for Array_G (e.g. for bilinear interpolation)

  use ModSize, ONLY: nI, nJ, nK

  implicit none
  real, intent(inout) :: Array_G(0:nI+1,0:nJ+1,0:nK+1)
  !---------------------------------------------------------------------------
  ! Edges in the K direction
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

