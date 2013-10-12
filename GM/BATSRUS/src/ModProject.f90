!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!This code is a copyright protected software (c) 2002- University of Michigan


module ModProject

  use ModSize, ONLY: nI, nJ, nK, nBlock, &
       MinI, MaxI, MinJ, MaxJ, MinK, MaxK, MaxBlock

  use ModProcMH, ONLY: iProc, nProc, iComm
  use ModMain, ONLY: Unused_B
  use ModGeometry, ONLY: true_BLK, true_cell
  use ModMpi

  implicit none
  save
  !\
  ! Parameters for projection scheme:
  !
  ! proj_method determines the iterative scheme to be used
  ! if proj_method='cg'       then Conjugate Gradient method is used
  !                           (this is only good for symmetric matrices!)
  !    proj_method='bicgstab' then BiConjugate Gradient method is used
  !
  ! proj_typestop determines the stopping condition for the iterations.
  ! if proj_typestop='rel' then sum((div B)**2) is reduced by a factor 
  !    proj_divbcoeff
  ! if proj_typestop='max' then |div B| is kept below a limit determeined by 
  !    proj_divbcoeff and
  !    proj_divbconst 
  !    if proj_divbcoeff>0 then 
  !       |div B|< proj_divbcoeff * divbmax_0  
  !          where divbmax_0 is the max(|div B|) after the first time step
  !    otherwise
  !       |div B|< proj_divbconst
  ! 
  ! proj_matvecmax is an upper limit on matrix-vector products for iterations
  !
  ! proj_boundtype determines the boundary conditions for the Poisson equation
  ! if proj_boundtype='zero' : phi=0 in ghost cells outside of comput. domain
  !
  !/

  character (len=10):: proj_method   ='cg        '
  character (len=3)::  proj_typestop ='rel'
  real ::              proj_divbcoeff=0.1
  real ::              proj_divbconst=0.0
  integer ::           proj_matvecmax=50
  character (len=10):: proj_boundtype='zero'

  ! Counter for matrix vector multiplications, and for errors of solver

  integer :: nmatvectotal, poissonerror

  ! Minimum value for divB (a small number)

  real, parameter :: divbmin=1E-10

  ! Maximum value for divB (absolute or relative)
  real :: DivbMax

contains
  subroutine set_BLK(qa,qb)

    ! Set qa=qb for all used blocks, where qb is a scalar

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), intent(out) :: qa
    real, intent(in) :: qb

    ! Local variables:
    integer:: iBLK

    !---------------------------------------------------------------------------

    do iBLK=1,nBlock
       if(.not.Unused_B(iBLK)) then
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

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock):: qa,qb
    intent(out) :: qa
    intent(in)  :: qb

    ! Local variables:
    integer:: iBLK

    !---------------------------------------------------------------------------

    do iBLK=1,nBlock
       if(.not.Unused_B(iBLK)) then
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

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock) :: qa, qb
    intent(inout) :: qa
    intent(in)    :: qb

    ! Local variables:
    integer:: iBLK

    !---------------------------------------------------------------------------

    do iBLK=1,nBlock
       if(.not.Unused_B(iBLK)) then
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

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock) :: qa, qb
    intent(inout) :: qa
    intent(in)    :: qb

    ! Local variables:
    integer:: iBLK

    !---------------------------------------------------------------------------

    do iBLK=1,nBlock
       if(.not.Unused_B(iBLK)) then
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

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock) :: qa,qb,qc
    intent(out) :: qa
    intent(in)  :: qb,qc

    ! Local variables:
    integer:: iBLK

    !---------------------------------------------------------------------------

    do iBLK=1,nBlock
       if(.not.Unused_B(iBLK)) then
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

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock) :: qa,qc
    intent(inout) :: qa
    intent(in)    :: qc

    real, intent(in) :: qb

    ! Local variables:
    integer:: iBLK

    !---------------------------------------------------------------------------

    do iBLK=1,nBlock
       if(.not.Unused_B(iBLK)) then
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

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock) :: qa,qb,qd
    intent(inout) :: qa
    intent(in)    :: qb
    intent(inout) :: qd

    real, intent(in) :: qc

    ! Local variables:
    integer:: iBLK

    !---------------------------------------------------------------------------

    do iBLK=1,nBlock
       if(.not.Unused_B(iBLK)) then
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

    ! Arguments

    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), intent(in) :: qa,qb

    ! Local variables:
    real    :: qproduct, qproduct_all
    integer :: iBLK, iError

    logical :: oktest, oktest_me

    !---------------------------------------------------------------------------

    call set_oktest('dot_product_BLK',oktest, oktest_me)

    qproduct=0.0

    do iBLK=1,nBlock
       if(.not.Unused_B(iBLK)) then
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

    ! Arguments

    integer, intent(in) :: qnum
    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), intent(in) :: qa

    ! Local variables:
    real    :: qsum, qsum_all
    integer :: iBLK, iError

    logical :: oktest, oktest_me

    !---------------------------------------------------------------------------

    call set_oktest('sum_BLK',oktest, oktest_me)

    qsum=0.0

    do iBLK=1,nBlock
       if(.not.Unused_B(iBLK)) then
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

end module ModProject
