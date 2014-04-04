!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!=============================================================================
subroutine impl_matvec(x_I, y_I, n)

  ! Calculate y=P_L.A.P_R.x for the iterative solver, where 
  ! P_L and P_R are the left and right preconditioner matrices,
  ! A = I - beta*dt*dR/dw, and R is the residual from dw/dt = R(w).
  !
  ! The multiplication by A is done in a matrix free fashion.

  use ModImplicit
  use ModLinearSolver, ONLY: precond_left_multiblock, precond_right_multiblock
  implicit none

  integer, intent(in):: n
  real, intent(in)   :: x_I(n)
  real, intent(out)  :: y_I(n)

  logical :: DoTest, DoTestMe
  character(len=*), parameter:: NameSub = 'impl_matvec'
  !----------------------------------------------------------------------------
  call set_oktest(NameSub, DoTest, DoTestMe)

  if(DoTestMe)write(*,*)NameSub,': initial n,sum(x**2)=', n, sum(x_I**2)

  if(ImplParam%DoPrecond)then

     ! y = P_R.x, where P_R = I, U^{-1}, or U^{-1}L^{-1}
     ! for left, symmetric and right preconditioning, respectively
     y_I = x_I
     call precond_right_multiblock(ImplParam, &
          nVar, nDim, nI, nJ, nK, nBlock, MAT, y_I)

     ! y = A.y
     call impl_matvec_free(y_I, y_I)

     ! y = P_L.y, where P_L==U^{-1}.L^{-1}, L^{-1}, or I
     ! for left, symmetric, and right preconditioning, respectively
     call precond_left_multiblock(ImplParam, &
          nVar, nDim, nI, nJ, nK, nBlock, MAT, y_I)
  else
     ! y = A.y
     call impl_matvec_free(x_I, y_I)
  end if

  if(DoTestMe)write(*,*)'impl_matvec_prec final n, sum(y**2)=', n, sum(y_I**2)

end subroutine impl_matvec

!=============================================================================
subroutine impl_matvec_free(x_I, y_I)

  ! Calculate y=L.x for the iterative solver, matrix-free 
  ! where L= I - beta*dt*dR/dw   (dt=dt_implicit)
  !
  ! One sided derivative:
  !----------------------
  ! ImplEps_VGB = Impl_VGB + eps*x            ! perturbation
  !
  ! ImplEps_VGB'=ImplEps_VGB + R(ImplEps_VGB,dtexpl)  ! advance ImplEps_VGB
  !
  ! dR/dw.x = (R(w+eps*x)-R(w))/eps 
  !          = [(ImplEps_VGB'-ImplEps_VGB) - (Impl_VGB'-Impl_VGB)]/eps/dtexpl
  !          = (ImplEps_VGB'-Impl_VGB')/eps/dtexpl - x/dtexpl
  !
  ! L.x = dx - beta*dt*dR/dw.x 
  !      = (1 + beta*dtcoeff)*x - beta*dtcoeff*(ImplEps_VGB' - w')/eps
  !
  ! where w=Impl_VGB, w'=w+R_low, beta=ImplCoeff, eps=sqrt(JacobianEps)/||x||
  ! instead of eps=(JacobianEps)^(1/2)*(Impl_VGB.x)/(x.x) suggested by Keyes

  use ModProcMH
  use ModMain, ONLY : Itest, Jtest, Ktest, VARtest
  use ModImplicit
  use ModMpi
  implicit none

  real, intent(in)   :: x_I(nImpl)
  ! Sometimes this subroutine called with the same array in both arguments
  ! that's why the intent of y cannot be set to out.
  real, intent(inout):: y_I(nImpl)

  real, allocatable, save:: ImplEps_VCB(:,:,:,:,:)
  integer:: n, i, j, k, iVar, iBlock, iError
  real:: Eps, xNorm, xNormTotal, Coef1, Coef2

  logical :: DoTest, DoTestMe
  character(len=*), parameter:: NameSub = 'impl_matvec_free'
  !----------------------------------------------------------------------------
  call set_oktest(NameSub, DoTest, DoTestMe)

  call timing_start(NameSub)

  if(.not.allocated(ImplEps_VCB)) &
       allocate(ImplEps_VCB(nVar,nI,nJ,nK,MaxImplBLK))

  xNorm = sum(x_I**2)
  call MPI_allreduce(xNorm, xNormTotal, 1, MPI_REAL, MPI_SUM,iComm,iError)

  if(DoTestMe)write(*,*) NameSub,': initial n,sum(x**2),xNormTotal=', &
       nImpl, xNorm, xNormTotal

  xNorm = sqrt(xNormTotal/nimpl_total)

  if(xNorm < SmallDouble) xNorm = 1.0

  Eps = sqrt(JacobianEps)/xNorm

  if(DoTestMe)write(*,*)'Eps, xNorm =',Eps,xNorm

  n=0
  do iBlock=1,nImplBLK; do k=1,nK; do j=1,nJ; do i=1,nI; do iVar=1,nVar
     n = n + 1
     ImplEps_VCB(iVar,i,j,k,iBlock) = Impl_VGB(iVar,i,j,k,iBlock) &
          + Eps*x_I(n)*wnrm(iVar)
  enddo; enddo; enddo; enddo; enddo

  ! Advance ImplEps_VCB:low order,  no dt, don't subtract
  call get_residual(.true., .false., .false., ImplEps_VCB, ImplEps_VCB) 

  ! Calculate y = L.x = (1 + beta*dtcoeff)*x 
  !                       - beta*dtcoeff*(ImplEps_VCB' - Impl_VGB')/eps
  ! where ImplEps_VCB  = Impl_VGB + eps*x, 
  !       ImplEps_VCB' = ImplEps_VCB + dt*R(ImplEps_VCB) and 
  !       Impl_VGB'    = Impl_VGB + dt*R(w)
  ! y = x + beta*dtcoeff*x 
  !       - beta*dtcoeff*(Impl_VGB + eps*x + R(ImplEps_VCB) - w - R(w))/eps
  !   = x - beta*dtcoeff*(R(ImplEps_VCB)-R(Impl_VGB))/eps 
  !   = x - beta*dt*dR/dU*x

  Coef1 = 1 + ImplCoeff*dtcoeff
  Coef2 = ImplCoeff*dtcoeff/Eps

  if(DoTestMe)write(*,*)'dtcoeff,ImplCoeff,Coef1,Coef2=', &
       dtcoeff,ImplCoeff,Coef1,Coef2

  n=0
  do iBlock=1,nImplBLK; do k=1,nK; do j=1,nJ; do i=1,nI; do iVar=1,nVar
     n=n+1
     y_I(n) = Coef1*x_I(n) - Coef2*(ImplEps_VCB(iVar,i,j,k,iBlock) &
          - Impl_VGB(iVar,i,j,k,iBlock) &
          - ResImpl_VCB(iVar,i,j,k,iBlock))/wnrm(iVar)
  enddo; enddo; enddo; enddo; enddo

  call timing_stop(NameSub)

  if(DoTestMe)write(*,*) NameSub,': final n,sum(y**2)=', nImpl, sum(y_I**2)

end subroutine impl_matvec_free

