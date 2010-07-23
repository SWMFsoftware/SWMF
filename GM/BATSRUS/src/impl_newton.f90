!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
subroutine impl_newton_init

  ! initialization for NR

  use ModProcMH
  use ModMain, ONLY : Itest,Jtest,Ktest,VARtest,n_step,dt,nOrder, &
       UseRadDiffusion, UseLaserPackage
  use ModAdvance, ONLY : FluxType
  use ModImplicit
  use ModMpi
  use ModRadDiffusion, ONLY: IsNewTimestepRadDiffusion
  implicit none

  integer :: i,j,k,n,iw,implBLK,iBLK, iError
  real :: coef1, coef2, q1, q2, q3

  logical :: oktest, oktest_me
  !---------------------------------------------------------------------------
  call set_oktest('impl_newton',oktest,oktest_me)

  ! Calculate high and low order residuals
  ! ResExpl_VCB= dtexpl * R

  if(UseSemiImplicit)then
     call get_semi_impl_rhs(Impl_VGB, ResExpl_VCB)

     !\
     ! If the laser package is used, this is the place to apply it
     !/
     if(UseLaserPackage) call add_laser_energy_deposition

     !!! We don't really need ResImpl_VCB, but it is easier if it is set
     ResImpl_VCB(:,:,:,:,1:nImplBLK) = ResExpl_VCB(:,:,:,:,1:nImplBLK)
  else
     if(UseRadDiffusion) IsNewTimestepRadDiffusion = .true.

     !                not low,  dt,  subtract
     call get_residual(.false.,.true.,.true., &
          Impl_VGB(:,1:nI,1:nJ,1:nK,:),ResExpl_VCB)

     if(UseRadDiffusion) IsNewTimestepRadDiffusion = .false.

     if (nOrder==nOrder_Impl .and. FluxType==FluxTypeImpl) then
        ! If R_low=R then ResImpl_VCB = ResExpl_VCB
        ResImpl_VCB(:,:,:,:,1:nImplBLK) = ResExpl_VCB(:,:,:,:,1:nImplBLK)
     else
        ! ResImpl_VCB = dtexpl * R_low
        !                  low,  no dt, subtract
        call get_residual(.true.,.false.,.true., &
             Impl_VGB(:,1:nI,1:nJ,1:nK,:), ResImpl_VCB) 
     endif
  end if
  if(oktest_me.and.nImplBLK>0)write(*,*)'ResExpl_VCB,ResImpl_VCB(test)=',&
       ResExpl_VCB(VARtest,Itest,Jtest,Ktest,implBLKtest),&
       ResImpl_VCB(VARtest,Itest,Jtest,Ktest,implBLKtest)

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
        rhs(n)=(coef1*ResExpl_VCB(iw,i,j,k,implBLK) &
             + coef2*(Impl_VGB(iw,i,j,k,implBLK) &
             -        ImplOld_VCB(iw,i,j,k,iBLK)))/wnrm(iw)
     end do; end do; enddo; enddo; enddo

  else
     do implBLK = 1, nImplBLK; do k=1,nK; do j=1,nJ; do i=1,nI
        do iw=1,nVarSemi
           n=n+1
           ! RHS = dt*R/wnrm for the first iteration
           rhs(n)=ResExpl_VCB(iw,i,j,k,implBLK)*dtcoeff/wnrm(iw)

     end do; end do; enddo; enddo; enddo

  endif

  if(UseNewton .or. UseConservativeImplicit)then
     ! Calculate RHS0 used for RHS when NewtonIter>1
     n=0
     do implBLK=1,nImplBLK; do k=1,nK; do j=1,nJ; do i=1,nI; do iw=1,nw
        n=n+1
        !RHS0 = [dt*(R - beta*R_low) + w_n]/wnrm 
        !     = RHS + [-beta*dt*R_low + w_n]/wnrm
        rhs0(n) = rhs(n) + (- ImplCoeff*dtcoeff*ResImpl_VCB(iw,i,j,k,implBLK) &
             + Impl_VGB(iw,i,j,k,implBLK))/wnrm(iw)
     end do; end do; enddo; enddo; enddo
  endif

  if(oktest)then
     call MPI_allreduce(sum(ResImpl_VCB(:,:,:,:,1:nImplBLK)**2),q1,&
          1,MPI_REAL,MPI_SUM,iComm,iError)
     call MPI_allreduce(sum(ResExpl_VCB(:,:,:,:,1:nImplBLK)**2),q2,&
          1,MPI_REAL,MPI_SUM,iComm,iError)
     call MPI_allreduce(sum(rhs(1:nimpl)**2),q3,&
          1,MPI_REAL,MPI_SUM,iComm,iError)

     if(oktest_me)write(*,*)'Sum ResExpl_VCB**2,ResImpl_VCB**2,rhs**2:', &
          q1, q2, q3
  end if

  ! Initial guess for dw = w_n+1 - w_n
  non0dw=.true.
  select case(KrylovInitType)
  case('explicit')
     ! w_n+1-w_n = dt * R_n
     dw(1:nimpl) = rhs(1:nimpl)
  case('scaled')
     ! Like explicit, but amplitude reduced
     ! w_n+1-w_n = dtexpl * R_n
     dw(1:nimpl) = rhs(1:nimpl)/dtcoeff
  case('nul')
     ! w_n+1-w_n = 0
     dw(1:nimpl) = 0.0
     non0dw = .false.
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
     ! RHS = (dt*R_n - beta*dt*R_n_low + w_n + beta*dt*R_k_low - Impl_VGB)/wnrm
     ! use: RHS0 and ResImpl_VCB = dtexpl * R_k_low
     rhs(n)= rhs0(n)+(ImplCoeff*dtcoeff*ResImpl_VCB(iw,i,j,k,implBLK) &
          - Impl_VGB(iw,i,j,k,implBLK))/wnrm(iw)
  enddo; enddo; enddo; enddo; enddo

  if(oktest)then
     call MPI_allreduce(sum(rhs(1:nimpl)**2),q1,1,MPI_REAL,MPI_SUM,&
          iComm,iError)
     if(oktest_me)then
        write(*,*)'norm of rhs:',sqrt(q1/nimpl_total)
        if(nImplBLK>0)write(*,*)'rhs,rhs0,ResImpl_VCB,Impl_VGB(test)=',&
             rhs(implVARtest),rhs0(implVARtest),               &
             ResImpl_VCB(Ktest,VARtest,Itest,Jtest,implBLKtest),  &
             Impl_VGB(VARtest,Itest,Jtest,Ktest,implBLKtest)
     end if
  end if

  ! Initial guess for dw is always zero in later NR iterations
  dw(1:nimpl)=0.0
  non0dw=.false.

end subroutine impl_newton_loop

!=============================================================================
subroutine impl_newton_update(dwnrm, converged)

  ! Update Impl_VGB(k+1) = Impl_VGB(k) + coeff*dw  with coeff from backtracking
  ! such that F(Impl_VGB+1) <= F(Impl_VGB) if possible

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
     ! dwnrm = ||Impl_VGB(k+1) - Impl_VGB(k)||/||w_n||
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
     do implBLK=1,nImplBLK; do k=1,nK; do j=1,nJ; do i=1,nI;
        do iw = iVarSemiMin, iVarSemiMax
        n=n+1
        if(true_cell(i,j,k,impl2iBLK(implBLK)))&
             Impl_VGB(iw,i,j,k,implBLK) = Impl_VGB(iw,i,j,k,implBLK) &
             + coeff*dw(n)*wnrm(iw)
     enddo; enddo; enddo; enddo; enddo

     if(UseConservativeImplicit .or. .not.Converged) then
        if(UseSemiImplicit)then
           call get_semi_impl_rhs(Impl_VGB, ResImpl_VCB)
        else
           !calculate low order residual ResImpl_VCB = dtexpl*RES_low(k+1)
           !                  low,   no dt, subtract
           call get_residual(.true., .false., .true., &
                Impl_VGB(:,1:nI,1:nJ,1:nK,:), ResImpl_VCB)
        end if
     end if

     ! Do not backtrack in a time accurate calculation or
     ! if Newton-Raphson converged or no Newton-Raphson is done
     if (time_accurate .or. converged) EXIT

     ! calculate high order residual ResExpl_VCB = dt*R(Impl_VGB(k+1))
     if ( (nOrder==nOrder_Impl .and. FluxType==FluxTypeImpl) &
          .or. UseSemiImplicit) then
        ResExpl_VCB(:,:,:,:,1:nImplBLK)=ResImpl_VCB(:,:,:,:,1:nImplBLK)
     else
        !                 not low, no dt, subtract
        call get_residual(.false.,.false.,.true.,Impl_VGB(:,1:nI,1:nJ,1:nK,:),&
             ResExpl_VCB)
     endif

     ! Calculate norm of high order residual
     residual = 0.0
     do iw = 1, nw
        call MPI_allreduce(sum(Impl_VGB(iw,1:nI,1:nJ,1:nK,1:nImplBLK)**2),&
             wnrm2,   1,MPI_REAL,MPI_SUM,iComm,iError)
        call MPI_allreduce(sum(ResExpl_VCB(iw,1:nI,1:nJ,1:nK,1:nImplBLK)**2),&
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

  ! Replace the final Newton iterate Impl_VGB with a flux based conservative update

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
          Impl_VGB(iW,i,j,k,ImplBlk) = &
          Rhs0(n)*wNrm(iW) + ImplCoeff*DtCoeff*ResImpl_VCB(iW,i,j,k,ImplBlk)
  enddo; enddo; enddo; enddo; enddo

end subroutine impl_newton_conserve
