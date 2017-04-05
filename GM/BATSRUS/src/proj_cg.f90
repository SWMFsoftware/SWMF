!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!This code is a copyright protected software (c) 2002- University of Michigan


! The CG-algorithm is implemented as shown on page 12 of the thesis
! "Preconditioning for sparse matrices with applications."
! Auke van der Ploeg, University of Groningen, 1994.
! Rewritten to F90 by G. Toth based on the F77 subroutine in src/conjgrad.f 
! Rewritten by G. Toth for BATS-R-US code to be used with MPI, 2000.
!
! This subroutine determines the solution of A.QX=RHS, where
! the matrix-vector multiplication with A is performed by 
! the subroutine 'proj_matvec'.
!
! !! If the matrix is not symmetric positive definite, CG is likely to fail.
!
subroutine proj_cg(rhs,qx,iter,tol,typestop,info)
  use ModProcMH
  use ModMain, ONLY :nBLK
  use ModAdvance, ONLY : tmp1_BLK,tmp2_BLK
  use ModProject
  implicit none

  ! Arguments

  real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nBLK), &
       intent(inout) :: rhs
  !        on input:  right-hand side vector.
  !        on output: residual vector.

  real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nBLK), &
       intent(out):: qx
  !        on output: solution vector.

  integer, intent(inout) :: iter
  !       on input:  maximum number of iterations to be performed.
  !       on output: actual  number of iterations done.

  real, intent(inout) :: tol
  !       on input:  required (relative) 2-norm or maximum norm of residual
  !       on output: achieved (relative) 2-norm or maximum norm of residual

  character (len=3), intent(in) :: typestop
  !       Determine stopping criterion (||.|| denotes the 2-norm):
  !       typestop='rel'    -- relative stopping crit.: ||res|| <= tol*||res0||
  !       typestop='abs'    -- absolute stopping crit.: ||res|| <= tol
  !       typestop='max'    -- maximum  stopping crit.: max(abs(res)) <= tol

  integer, intent(out)   :: info
  !       Gives reason for returning:
  !     abs(info)=  0 - solution found satisfying given tolerance.
  !                 1 - iteration aborted due to division by very small value.
  !                 2 - no convergence within maximum number of iterations.
  !                 3 - initial guess satisfies the stopping criterion.
  !    sign(info)=  + - residual decreased
  !                 - - residual did not reduce

  ! Local variables

  integer ::itr,matv
  real :: rho,rhonew,res,res0,bet,alf,assumedzero

  logical :: oktest, oktest_me

  ! external Functions:

  real :: maxval_abs_BLK

  !---------------------------------------------------------------------------

  call set_oktest('proj_cg',oktest, oktest_me)

  if(oktest_me)write(*,*)'Proj_CG: tol, iter=',tol,iter

  !Debug
  !call show_BLK('rhs',rhs)

  assumedzero=1.E-16; itr=0; matv=0
  call set_BLK(qx,0.0)

  if (typestop/='rel'.and.typestop/='abs'.and.typestop/='max') &
       call stop_mpi('Error in CG: typestop='//typestop// &
       ' should be one of rel/abs/max.')

  if(oktest_me) write(*,*)'n gives the number of CG-iterations.'

  ! rho=||rhs||
  rho=sqrt(dot_product_BLK(rhs,rhs))

  res0=rho
  if (typestop=='max') res0=maxval_abs_BLK(nProc,rhs)
  !DEBUG
  !write(*,*)'res0:',res0,' rhs:',rhs

  res=res0

  assumedzero = assumedzero*res0

  if (oktest_me) then
     if (typestop=='max') then
        write(*,*)'n:',itr,' Maximum norm initial residual:',res0
     else
        write(*,*)'n:',itr,' 2-norm initial residual:',res0
     end IF
  end if

  if (res0<divbmin .or. (typestop/='rel'.and.res0<=tol)) then
     info = 3
  else
     ! Initialize rho and tmp1_BLK=Z
     rho=rho*rho
     call eq_BLK(tmp1_BLK,rhs)

     ! Do iteration
     do
        ! AZ=A.Z
        call proj_matvec(tmp1_BLK,tmp2_BLK)
        matv=matv+1

        !Debug
        !write(*,*)'Z =tmp1_BLK=',tmp1_BLK(:,:,:,1:2)
        !write(*,*)'AZ=tmp2_BLK=',tmp2_BLK(:,:,:,1:2)
        !call stop_mpi('Debug')

        ! alf=A.AZ
        alf=dot_product_BLK(tmp1_BLK,tmp2_BLK)

        if(oktest_me)write(*,*)'alf=',alf
        !Debug
        !call show_BLK('Z',tmp1_BLK)
        !call show_BLK('AZ',tmp2_BLK)
        !call stop_mpi('Debug')
        if (abs(alf)<=assumedzero**2) then
           info = 1
           exit
        end if
        alf=rho/alf
        if(oktest_me)write(*,*)'alf=',alf

        call add_times_BLK(qx,alf,tmp1_BLK)
        call add_times_BLK(rhs,-alf,tmp2_BLK)

        ! rhonew=||rhs||
        rhonew=sqrt(dot_product_BLK(rhs,rhs))
        if(oktest_me)write(*,*)'rhonew=',rhonew

        select case(typestop)
        case('max')
           res=maxval_abs_BLK(nProc,rhs)

        case('rel')
           res=rhonew/res0

        case('abs')
           res=rhonew
        end select
        rhonew=rhonew*rhonew

        itr=itr+1
        if (oktest_me) &
             write(*,*)'n:',itr,' ',typestop,'. norm of residual:',res

        if (res<=tol) then
           info = 0
           exit
        end if
        if (itr>=iter) then
           info = 2
           exit
        end if
        if (rho<=assumedzero**2) then
           info = 1
           exit
        end if

        bet=rhonew/rho
        if(oktest_me)write(*,*)'bet=',bet

        call eq_plus_times_BLK(tmp1_BLK,rhs,bet,tmp1_BLK)

        if (oktest_me) write(*,*)'alf,bet,rho,rhonew:',alf,bet,rho,rhonew
        rho=rhonew
     end do
  end if

  ! return number of iterations and achieved residual
  iter=itr
  tol =res
  if((typestop=='rel'.and.res>1.0).or.(typestop/='rel'.and.res>res0))info=-info

  ! report results
  if(oktest_me)then
     write(*,*)'Total Number of CG-iterations:',itr
     write(*,*)'Number of matrix-vector mult.:',matv
     select case(abs(info))
     case(0)
        write(*,*)'Successful iteration, norm of res. is:',tol
     case(1)
        write(*,*)'Iteration aborted due to division by a'
        write(*,*)'very small value.' 
     case(2)
        write(*,*)'Stopping crit. not fulfilled within given'
        write(*,*)'maximum number of iterations.'
     case(3)
        write(*,*)'Initial guess for the solution satisfies'
        write(*,*)'given stopping criterion.' 
     case default
        write(*,*)'Impossible value for info:',info
     end select
     if(info<0)write(*,*)'The residual did not reduce'
  endif

end subroutine proj_cg
