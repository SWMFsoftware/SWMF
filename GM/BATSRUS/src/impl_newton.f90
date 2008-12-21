!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
subroutine impl_newton_init

  ! initialization for NR

  use ModProcMH
  use ModMain, ONLY : Itest,Jtest,Ktest,VARtest,n_step,dt,nOrder, &
       UseGrayDiffusion
  use ModAdvance, ONLY : FluxType
  use ModImplicit
  use ModMpi
  use ModGrayDiffusion, ONLY: IsNewTimestepGrayDiffusion
  implicit none

  integer :: i,j,k,n,iw,implBLK,iBLK, iError
  real :: coef1, coef2, q1, q2, q3

  logical :: oktest, oktest_me
  !---------------------------------------------------------------------------
  call set_oktest('impl_newton',oktest,oktest_me)

  ! Calculate high and low order residuals
  ! RES_expl= dtexpl * R

  if(UseSemiImplicit)then
     call get_semi_impl_rhs(w_k, RES_expl)
     !!! We don't really need RES_impl, but it is easier if it is set
     RES_impl = RES_expl
  else
     if(UseGrayDiffusion) IsNewTimestepGrayDiffusion = .true.

     !                not low,  dt,  subtract
     call get_residual(.false.,.true.,.true.,w_k(1:nI,1:nJ,1:nK,:,:),RES_expl)

     if(UseGrayDiffusion) IsNewTimestepGrayDiffusion = .false.

     if (nOrder==nOrder_Impl .and. FluxType==FluxTypeImpl) then
        ! If R_low=R then RES_impl = RES_expl
        RES_impl(:,:,:,:,1:nImplBLK)=RES_expl(:,:,:,:,1:nImplBLK)
     else
        ! RES_impl = dtexpl * R_low
        !                  low,  no dt, subtract
        call get_residual(.true.,.false.,.true.,w_k(1:nI,1:nJ,1:nK,:,:),RES_impl) 
     endif
  end if
  if(oktest_me.and.nImplBLK>0)write(*,*)'RES_expl,RES_impl(test)=',&
       RES_expl(Itest,Jtest,Ktest,VARtest,implBLKtest),&
       RES_impl(Itest,Jtest,Ktest,VARtest,implBLKtest)

  ! Calculate rhs used for NewtonIter=1
  n=0
  if(UseBDF2.and.n_step==n_prev+1)then
     ! Collect RHS terms from Eq 8 in Paper implvac
     ! Newton-Raphson iteration. The BDF2 scheme implies
     ! beta+alpha=1 and beta=(dt_n+dt_n-1)/(2*dt_n+dt_n-1)
     coef1=ImplCoeff*dtcoeff
     coef2=(1-ImplCoeff)*dt/dt_prev
     do implBLK=1,nImplBLK; do k=1,nK; do j=1,nJ; do i=1,nI; do iw=1,nw
        iBLK=impl2iBLK(implBLK)
        n=n+1
        ! For 1st Newton iteration
        ! RHS = dt*(beta*R + alpha*(w_n-w_n-1)/dt_n-1)/wnrm 
        rhs(n)=(coef1*RES_expl(i,j,k,iw,implBLK) &
             +coef2*(w_k(i,j,k,iw,implBLK)-w_prev(i,j,k,iw,iBLK)))/wnrm(iw)
     end do; end do; enddo; enddo; enddo
  else
     do implBLK=1,nImplBLK; do k=1,nK; do j=1,nJ; do i=1,nI; do iw=1,nw
        n=n+1
        ! RHS = dt*R/wnrm for the first iteration
        rhs(n)=RES_expl(i,j,k,iw,implBLK)*dtcoeff/wnrm(iw)

        !!!DEBUG
        !if(i==Itest.and.j==Jtest.and.k==Ktest)write(*,*)'iw,RES_expl,rhs=',&
        !     iw,RES_expl(i,j,k,iw),rhs(n)

     end do; end do; enddo; enddo; enddo
  endif

  if(UseNewton .or. UseConservativeImplicit)then
     ! Calculate RHS0 used for RHS when NewtonIter>1
     n=0
     do implBLK=1,nImplBLK; do k=1,nK; do j=1,nJ; do i=1,nI; do iw=1,nw
        n=n+1
        !RHS0 = [dt*(R - beta*R_low) + w_n]/wnrm 
        !     = RHS + [-beta*dt*R_low + w_n]/wnrm
        rhs0(n) = rhs(n) + (- ImplCoeff*dtcoeff*RES_impl(i,j,k,iw,implBLK) &
             + w_k(i,j,k,iw,implBLK))/wnrm(iw)
     end do; end do; enddo; enddo; enddo
  endif

  if(oktest)then
     call MPI_allreduce(sum(RES_impl(:,:,:,:,1:nImplBLK)**2),q1,&
          1,MPI_REAL,MPI_SUM,iComm,iError)
     call MPI_allreduce(sum(RES_expl(:,:,:,:,1:nImplBLK)**2),q2,&
          1,MPI_REAL,MPI_SUM,iComm,iError)
     call MPI_allreduce(sum(rhs(1:nimpl)**2),q3,&
          1,MPI_REAL,MPI_SUM,iComm,iError)

     if(oktest_me)write(*,*)'Sum RES_expl**2,RES_impl**2,rhs**2:',q1,q2,q3
  end if

  ! Initial guess for dw = w_n+1 - w_n
  non0dw=.true.
  select case(KrylovInitType)
  case('explicit')
     ! w_n+1-w_n = dt * R_n
     dw(1:nimpl)=rhs(1:nimpl)
  case('scaled')
     ! Like explicit, but amplitude reduced
     ! w_n+1-w_n = dtexpl * R_n
     dw(1:nimpl)=rhs(1:nimpl)/dtcoeff
  case('nul')
     ! w_n+1-w_n = 0
     dw(1:nimpl)=0.0
     non0dw=.false.
  case('old')
  case default
     call stop_mpi('Unknown type for KrylovInitType='//KrylovInitType)
  end select

end subroutine impl_newton_init
!=============================================================================

subroutine impl_newton_loop

  use ModProcMH
  use ModMain, ONLY : Itest,Jtest,Ktest,VARtest
  use ModImplicit
  use ModMpi
  implicit none

  integer :: i,j,k,iw,implBLK,n, iError
  real    :: q1
  logical :: oktest, oktest_me
  !---------------------------------------------------------------------------
  call set_oktest('impl_newton',oktest,oktest_me)

  ! Caculate RHS for 2nd or later Newton iteration
  n=0
  do implBLK=1,nImplBLK; do k=1,nK; do j=1,nJ; do i=1,nI; do iw=1,nw
     n=n+1
     ! RHS = (dt*R_n - beta*dt*R_n_low + w_n + beta*dt*R_k_low - w_k)/wnrm
     ! use: RHS0 and RES_impl = dtexpl * R_k_low
     rhs(n)= rhs0(n)+(ImplCoeff*dtcoeff*RES_impl(i,j,k,iw,implBLK) &
          - w_k(i,j,k,iw,implBLK))/wnrm(iw)
  enddo; enddo; enddo; enddo; enddo

  if(oktest)then
     call MPI_allreduce(sum(rhs(1:nimpl)**2),q1,1,MPI_REAL,MPI_SUM,&
          iComm,iError)
     if(oktest_me)then
        write(*,*)'norm of rhs:',sqrt(q1/nimpl_total)
        if(nImplBLK>0)write(*,*)'rhs,rhs0,RES_impl,w_k(test)=',&
             rhs(implVARtest),rhs0(implVARtest),               &
             RES_impl(Itest,Jtest,Ktest,VARtest,implBLKtest),  &
             w_k(Itest,Jtest,Ktest,VARtest,implBLKtest)
     end if
  end if

  ! Initial guess for dw is always zero in later NR iterations
  dw(1:nimpl)=0.0
  non0dw=.false.

end subroutine impl_newton_loop

!=============================================================================
subroutine impl_newton_update(dwnrm, converged)

  ! Update w: w_k+1 = w_k + coeff*dw  with coeff from backtracking
  ! such that F(w_k+1) <= F(w_k) if possible

  use ModProcMH
  use ModMain, ONLY : nOrder,time_accurate
  use ModAdvance, ONLY : FluxType
  use ModGeometry, ONLY : true_cell
  use ModImplicit
  use ModMpi
  implicit none

  real,    intent(out):: dwnrm
  logical, intent(out):: converged

  integer:: i, j, k, iw, implBLK, n, itry, iError
  real:: coeff, resold, dwnrm_local, wnrm2, resexpl2

  logical :: oktest,oktest_me
  !---------------------------------------------------------------------------

  call set_oktest('impl_newton',oktest,oktest_me)

  if(UseNewton)then
     ! Calculate progress in NR scheme to set linear solver accuracy
     ! dwnrm=||w_k+1 - w_k||/||w_n||
     dwnrm_local=sum(dw(1:nimpl)**2)
     call MPI_allreduce(dwnrm_local,dwnrm,1,MPI_REAL,MPI_SUM,iComm,&
          iError)
     dwnrm=sqrt(dwnrm/nimpl_total)
     converged = dwnrm<KrylovErrorMax
     if(oktest_me)write(*,*)'dwnrm:',dwnrm
  else
     converged = .true.
  endif

  ! Initial guess for coeff is limited to avoid non-physical w for pseudo-time
  !if((.not.time_accurate).and.impldwlimit<bigdouble)then
  !  coeff=min(1.0,impldwlimit/(maxval(abs(dw(1:nimpl)))+smalldouble))
  !else
  coeff=1.0
  !endif

  ! For steady state calculations try changing coeff to reduce the residual
  itry=0
  resold=residual
  do
     itry=itry+1
     if(oktest_me)write(*,*)'itry, coeff:',itry,coeff

     ! w=w+dw for all true cells
     n=0
     do implBLK=1,nImplBLK; do k=1,nK; do j=1,nJ; do i=1,nI; do iw=1,nw
        n=n+1
        if(true_cell(i,j,k,impl2iBLK(implBLK)))&
             w_k(i,j,k,iw,implBLK)=w_k(i,j,k,iw,implBLK)+coeff*dw(n)*wnrm(iw)
     enddo; enddo; enddo; enddo; enddo

     if(UseConservativeImplicit .or. .not.Converged) then
        if(UseSemiImplicit)then
           call get_semi_impl_rhs(w_k, RES_impl)
        else
           !calculate low order residual RES_impl = dtexpl*RES_low_k+1
           !                  low,   no dt, subtract
           call get_residual(.true.,.false.,.true.,w_k(1:nI,1:nJ,1:nK,:,:),&
                RES_impl)
        end if
     end if

     ! Do not backtrack in a time accurate calculation or
     ! if Newton-Raphson converged or no Newton-Raphson is done
     if (time_accurate .or. converged) EXIT

     ! calculate high order residual RES_expl = dt*R(w_k+1)
     if ( (nOrder==nOrder_Impl .and. FluxType==FluxTypeImpl) &
          .or. UseSemiImplicit) then
        RES_expl(:,:,:,:,1:nImplBLK)=RES_impl(:,:,:,:,1:nImplBLK)
     else
        !                 not low, no dt, subtract
        call get_residual(.false.,.false.,.true.,w_k(1:nI,1:nJ,1:nK,:,:),&
             RES_expl)
     endif

     ! Calculate norm of high order residual
     residual=0.0
     do iw=1,nw
        call MPI_allreduce(sum(w_k(1:nI,1:nJ,1:nK,iw,1:nImplBLK)**2),&
             wnrm2,   1,MPI_REAL,MPI_SUM,iComm,iError)
        call MPI_allreduce(sum(RES_expl(1:nI,1:nJ,1:nK,iw,1:nImplBLK)**2),&
             resexpl2,1,MPI_REAL,MPI_SUM,iComm,iError)

        if(wnrm2<smalldouble)wnrm2=1.0
        residual = residual + resexpl2/wnrm2
     enddo
     residual=sqrt(residual/nw)
     if(oktest_me)write(*,*)'resold,residual:',resold,residual

     ! Exit if backtracked towards steady-state or giving up
     if(residual<=resold .or. itry>3)exit
     coeff=coeff*0.5
  end do

end subroutine impl_newton_update

!==============================================================================
subroutine impl_newton_conserve

  ! Replace the final Newton iterate w_k with a flux based conservative update

  use ModImplicit
  use ModGeometry, ONLY : true_cell
  implicit none
  integer :: i,j,k,iw,n,implBLK
  !---------------------------------------------------------------------------

  ! w = Rhs0*wNrm + ImplCoeff*DtCoeff*ResImpl
  !
  ! Rhs0 is the normalized iteration independent part of the right hand side,
  ! which is calculated in impl_newton_init.
  !
  ! wNrm converts each variable from the normalized (second norm=1) 
  ! units to the units used in the explicit part BATSRUS
  !
  ! ResImpl is the (low order) residual obtained from the final newton iterate
  ! with a DtExpl time step. It is calculated in impl_newton_update.
  !
  ! DtCoeff = Dt/DtExpl is used to convert to the implicit time step

  n=0
  do ImplBlk=1,nImplBlk; do k=1,nK; do j=1,nJ; do i=1,nI; do iW=1,nW
     n=n+1
     if(true_cell(i,j,k,impl2iBLK(ImplBlk))) &
          w_k(i,j,k,iW,ImplBlk) = &
          Rhs0(n)*wNrm(iW) + ImplCoeff*DtCoeff*Res_Impl(i,j,k,iW,ImplBlk)
  enddo; enddo; enddo; enddo; enddo

end subroutine impl_newton_conserve
