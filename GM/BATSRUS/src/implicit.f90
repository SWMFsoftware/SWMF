!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
subroutine OPTION_IMPLICIT(on,name)

  logical, intent(out) :: on
  character (len=40), intent(out) :: name

  on  =.true.
  name='IMPLICIT SCHEME Toth 1.3'

end subroutine OPTION_IMPLICIT

!=============================================================================
subroutine advance_impl

  ! The implicit schemes used in this module are described in detail in
  !
  ! Keppens, Toth, Botchev, van der Ploeg, 
  ! J. Int. Num. Methods in Fluids, 30, 335-352, 1999
  !
  ! Equation numbers below refer to this paper unless stated otherwise.
  !
  ! We solve the MHD equation written as 
  !
  !    dw/dt = R(t)                                                 (1)
  !
  ! by one of the following implicit schemes:
  !
  ! If UseBDF2 is false (and in any case for the 1st time step):
  ! 
  !    w^n+1 = w^n + dt^n*[R^n + ImplCoeff*(R_low^n+1 - R_low^n)]   (4)
  !
  ! where ImplCoeff is a fixed parameter in the [0,1] range. 
  ! Here R is a high order while R_imp is a possibly low order discretization.
  ! The low order scheme is typically the first order Rusanov scheme. 
  !
  ! If UseBDF2 is true (except for the 1st time step):
  !
  !    w^n+1 = w^n + dt^n*[ ImplCoeff*(R^n + R_low^n+1 - R_low^n)
  !
  !                        + (1-ImplCoeff)*(w^n - w^n-1)/dt^n-1]    (8)
  !
  ! where
  !
  !    ImplCoeff = (dt^n + dt^n-1)/(2*dt^n + dt^n-1) 
  !
  ! provides second order time accuracy.
  !
  ! A Newton iteration is used to solve the discrete equations (4) or (8):
  !
  !    w^(k=0) = w^n
  !
  ! Solve
  !
  !    (I - dt*ImplCoeff*dR_low/dw).dw = ImplCoeff*dt*R(w^n) 
  !
  !            + (1-ImplCoeff)*(w^n - w^n-1)*dt/dt^n-1 
  !
  !            + ImplCoeff*dt*[R_low(w^k) - R_low^n] + (w^n - w^k)
  !
  ! for the increment dw. Terms in the second line are only included for BDF2, 
  ! while terms in the third line are zero for the first k=0 iteration.
  !
  ! In each iteration update the iterate as
  !
  !    w^k+1 = w^k + dw
  !
  ! At most NewtonIterMax iterations are done. When the Newton iteration
  ! is finished, the solution is updated as
  !
  !    w^n+1 = w^k
  !
  ! In each iteration the linear problem is solved by a Krylov type iterative 
  ! method.
  !
  ! We use get_residual(.false.,...) to calculate R_expl
  ! and    get_residual(.true.,....) to calculate R_impl

  use ModProcMH
  use ModMain
  use ModVarIndexes
  use ModAdvance, ONLY : State_VGB, E_BLK, StateOld_VCB, E_o_BLK,UseUpdateCheck
  use ModPhysics, ONLY : gm1
  use ModImplicit
  use ModAMR, ONLY : unusedBlock_BP
  use ModNumConst
  use ModMpi
  implicit none

  integer :: n,iVar, i, j, k, iw, implBLK, iBLK, KrylovMatVec, info, NewtonIter
  integer :: iError
  real    :: KrylovError, dwnrm, local_wnrm(nw), coef1, coef2

  logical:: converged
  character (LEN=3) :: typestop

  real*8  :: time_before,time_before_all
  logical :: oktest, oktest_me, DoTestKrylov, DoTestKrylovMe

  logical :: UseUpdateCheckOrig

  external impl_matvec

  !----------------------------------------------------------------------------

  call set_oktest('implicit',oktest,oktest_me)

  UseUpdateCheckOrig = UseUpdateCheck
  UseUpdateCheck = .false.

  ! Store unusedBLK into skippedBLK
  skippedBLK(1:nBlockMax)=unusedBLK(1:nBlockMax)

  ! Advance explicitly treated blocks if any
  if(UsePartImplicit .and. (.not. &
       all(skippedBLK(1:nBlockMax).or.implicitBLK(1:nBlockMax))))then

     if(UseBDF2)then

        ! Save the current state into the previous state for BDF2 scheme
        ! This is needed for explicit blocks, because they may become
        ! implicit in the next time step...
        do iBLK=1,nBlock
           if(skippedBLK(iBLK).or.implicitBLK(iBLK))CYCLE
           do iVar=1,nVar
              if(iVar==E_)then
                 w_prev(:,:,:,E_    ,iBLK)=    E_BLK(1:nI,1:nJ,1:nK,iBLK)
              else
                 w_prev(:,:,:,iVar,iBLK)=  &
                      State_VGB(iVar,1:nI,1:nJ,1:nK,iBLK)
              end if
           end do
        end do
     end if

     ! Select unusedBLK = not explicit blocks = skipped or implicit blocks
     iNewDecomposition=mod(iNewDecomposition+1, 10000)
     unusedBLK(1:nBlockMax) = skippedBLK(1:nBlockMax) .or. &
          implicitBLK(1:nBlockMax)
     call MPI_ALLGATHER(unusedBLK,      nBLK, MPI_LOGICAL, &
          unusedBlock_BP, nBLK, MPI_LOGICAL, iComm, iError)

     ! advance explicit blocks, calc timestep but do not exchange messages
     ! after the last stage (it may exchange messages in the first stage)
     if(.not.UseDtFixed)cfl=explCFL
     call advance_expl(.true.,.false.) 

     ! update ghost cells for the implicit blocks to time level n+1
     iNewDecomposition=mod(iNewDecomposition+1, 10000)
     unusedBLK(1:nBlockMax) = skippedBLK(1:nBlockMax)
     call MPI_ALLGATHER(unusedBLK,      nBLK, MPI_LOGICAL, &
          unusedBlock_BP, nBLK, MPI_LOGICAL, iComm, iError)

     call exchange_messages

     ! The implicit scheme is only applied on implicit blocks
     iNewDecomposition=mod(iNewDecomposition+1, 10000)
     unusedBLK(1:nBlockMax) = .not.implicitBLK(1:nBlockMax)
     call MPI_ALLGATHER(unusedBLK,      nBLK, MPI_LOGICAL, &
          unusedBlock_BP, nBLK, MPI_LOGICAL, iComm, iError)
  end if

  !\
  ! Advance implicitly treated blocks
  !/

  ! Use implicit time step
  if(.not.UseDtFixed)cfl=implCFL

  ! Initialize some variables in ModImplicit
  call implicit_init

  ! Get initial iterate from current state
  call explicit2implicit(0,nI+1,0,nJ+1,0,nK+1,w_k)
  if(oktest_me)write(*,*)'Implicit initial nImplBLK=',nImplBLK
  if(oktest_me.and.nImplBLK>0)write(*,*)'w_k=',&
       w_k(Itest,Jtest,Ktest,VARtest,implBLKtest)

  wnrm(1:nw)=-1.0
  ! Global norm of current w_(k=0) = w_n
  do iw=1,nw
     local_wnrm(iw)=sum(w_k(1:nI,1:nJ,1:nK,iw,1:nImplBLK)**2)
  end do
  call MPI_allreduce(local_wnrm, wnrm, nw, MPI_REAL, MPI_SUM, iComm,iError)

  call MPI_allreduce(nimpl,nimpl_total, 1,MPI_INTEGER,MPI_SUM,iComm,iError)
  wnrm=sqrt(wnrm/(nimpl_total/nw))
  where(wnrm < smalldouble) wnrm =1.0

  if(oktest_me)write(*,*)'wnrm:',nimpl_total,wnrm(bat2vac)

  if(UseDtFixed)then
     if(oktest_me)write(*,*)'call getdt_courant'
     call getdt_courant(dtexpl)
     dtexpl=cHalf*dtexpl
     dtcoeff=dt/dtexpl
  else
     if(oktest_me)write(*,*)'no call of getdt_courant'
     dtcoeff=implCFL/cHalf
  endif

  if (UseBDF2.and.n_step==n_prev+1) then
     ! For 3 level BDF2 scheme set beta=ImplCoeff if previous state is known
     ImplCoeff = (dt+dt_prev)/(2*dt+dt_prev)
  else
     ImplCoeff = ImplCoeff0
  end if

  if(oktest_me.and.time_accurate)&
       write(*,*)'dtcoeff,dtexpl,dt:',dtcoeff,dtexpl,dt
  if(oktest_me.and.UseBDF2)write(*,*)'n_prev,dt_prev,ImplCoeff:',&
       n_prev,dt_prev,ImplCoeff

  ! Initialize right hand side and dw
  call impl_newton_init

  ! Save previous timestep for 3 level scheme
  if(UseBDF2)then
     n_prev  = n_step
     dt_prev = dt
  endif
  ! Save the current state into w_prev for BDF2 scheme and/or update check
  if(UseBDF2 .or. UseUpdateCheckOrig)then
     ! The implicit blocks haven't been updated, so save current state
     do implBLK=1,nImplBlk
        iBLK=impl2iBLK(implBLK)
        w_prev(:,:,:,1:nw,iBLK) = w_k(1:nI,1:nJ,1:nK,1:nw,implBLK)
     end do
  endif

  ! Newton-Raphson iteration and iterative linear solver
  dwnrm=bigdouble
  NewtonIter=0
  do
     NewtonIter=NewtonIter+1;
     if(oktest_me)write(*,*)'NewtonIter:',NewtonIter
     if(NewtonIter>NewtonIterMax)then
        write(*,*)'Newton-Raphson failed to converge NewtonIter=',NewtonIter
        if(time_accurate)call stop_mpi('Newton-Raphson failed to converge')
        exit
     endif
     nnewton=nnewton+1

     ! Calculate Jacobian matrix if required
     if(JacobianType/='free'.and.(NewtonIter==1.or.NewMatrix))then

        if(NewtonIter>1)then
           ! Update ghost cells for w_k, because it is needed by impl_jacobian
           call implicit2explicit(w_k(1:nI,1:nJ,1:nK,:,:))
           call exchange_messages
           call explicit2implicit(0,nI+1,0,nJ+1,0,nK+1,w_k)
        end if

        call timing_start('impl_jacobian')
        do implBLK=1,nImplBLK
           call impl_jacobian(implBLK,MAT(1,1,1,1,1,1,implBLK))
        end do
        call timing_stop('impl_jacobian')

        if(oktest)then
           call MPI_reduce(sum(MAT(:,:,:,:,:,:,1:nImplBLK)**2),coef1,1,MPI_REAL,&
                MPI_SUM,PROCtest,iComm,iError)
           if(oktest_me)write(*,*)'sum(MAT**2)=',coef1
        end if

     endif

     ! Update rhs and initial dw if required
     if (NewtonIter>1) call impl_newton_loop

     if(oktest_me.and.nImplBLK>0)write(*,*)'initial dw(test), rhs(test)=',&
          dw(implVARtest),rhs(implVARtest)

     ! Precondition matrix if required
     if(JacobianType=='prec'.and.(NewtonIter==1.or.NewMatrix))then
        do implBLK=1,nImplBLK
           ! Preconditioning: MAT --> LU
           call prehepta(nIJK,nw,nI,nI*nJ,-GustafssonPar,&
                MAT(1,1,1,1,1,1,implBLK),&
                MAT(1,1,1,1,1,2,implBLK),&
                MAT(1,1,1,1,1,3,implBLK),&
                MAT(1,1,1,1,1,4,implBLK),&
                MAT(1,1,1,1,1,5,implBLK),&
                MAT(1,1,1,1,1,6,implBLK),&
                MAT(1,1,1,1,1,7,implBLK))

           ! rhs --> P_L.rhs, where P_L=U^{-1}.L^{-1}, L^{-1}, or I
           ! for left, symmetric, and right preconditioning, respectively
           if(PrecondSide/='right')then
              call Lhepta(nIJK,nw,nI,nI*nJ,&
                   rhs(nwIJK*(implBLK-1)+1),&
                   MAT(1,1,1,1,1,1,implBLK),&
                   MAT(1,1,1,1,1,2,implBLK),&
                   MAT(1,1,1,1,1,4,implBLK),&
                   MAT(1,1,1,1,1,6,implBLK))
              if(PrecondSide=='left') &
                   call Uhepta(.true.,nIJK,nw,nI,nI*nJ,&
                   rhs(nwIJK*(implBLK-1)+1),  &
                   MAT(1,1,1,1,1,3,implBLK),  &   ! +i diagonal
                   MAT(1,1,1,1,1,5,implBLK),  &   ! +j
                   MAT(1,1,1,1,1,7,implBLK))      ! +k
           end if

           ! Initial guess x --> P_R^{-1}.x where P_R^{-1} = I, U, LU for
           ! left, symmetric and right preconditioning, respectively
           ! Multiplication with LU is NOT implemented
           if(non0dw .and. PrecondSide=='symmetric') &
                call Uhepta(.false.,nIJK,nw,nI,nI*nJ,&
                dw(nwIJK*(implBLK-1)+1),   &
                MAT(1,1,1,1,1,3,implBLK),  &   ! +i diagonal
                MAT(1,1,1,1,1,5,implBLK),  &   ! +j
                MAT(1,1,1,1,1,7,implBLK))      ! +k
        end do

        if(oktest)then
           call MPI_reduce(sum(MAT(:,:,:,:,:,:,1:nImplBLK)**2),coef1,1,MPI_REAL,&
                MPI_SUM,procTEST,iComm,iError)
           call MPI_reduce(sum(rhs(1:nimpl)**2),coef2,1,MPI_REAL,&
                MPI_SUM,procTEST,iComm,iError)
           if(oktest_me)then
              write(*,*)'preconditioned sum(MAT**2), sum(rhs**2)=',coef1,coef2
              if(nImplBLK>0)&
                   write(*,*)'preconditioned dw(test)   , rhs(test)  =',&
                   dw(implVARtest),  rhs(implVARtest)
           endif
        end if
     endif

     ! Set tolerance and stopping conditions for iterative solver

     if(UseNewton)then
        ! link the inner iterative solver with the outer NR
        ! require inner loop more accurate than outer loop
        KrylovError=0.1
        typestop='rel'
     else if(JacobianType=='prec')then
        ! No normalization is needed for preconditioned solvers (???!!!)
        KrylovError=KrylovErrorMax
        typestop='rel'
     else
        ! Normalize KrylovError by dt(expl) for non-preconditioned solvers
        ! Distinguish between steady state and time accurate cases
        if(time_accurate)then
           KrylovError=KrylovErrorMax
        else
           KrylovError=KrylovErrorMax*dtcoeff
        endif
        typestop='abs'
     endif

     KrylovMatVec=KrylovMatvecMax

     if(oktest_me)write(*,*)'Before ',KrylovType,' KrylovMatVec,KrylovError:',&
          KrylovMatVec,KrylovError

     ! Solve linear problem

     call set_oktest('krylov',DoTestKrylov,DoTestKrylovMe)
     call timing_start('krylov solver')
     select case(KrylovType)                       !^CFG IF NOT SIMPLE BEGIN
     case('bicgstab','BiCGSTAB')
        call bicgstab(impl_matvec,rhs,dw,non0dw,nimpl,&
             KrylovError,typestop,KrylovMatVec,info,DoTestKrylovMe)
     case('GMRES','gmres')                         !^CFG END SIMPLE
        call gmres(impl_matvec,rhs,dw,non0dw,nimpl,nKrylovVector, &
             KrylovError,typestop,KrylovMatVec,info,DoTestKrylovMe)
     case default                                  !^CFG IF NOT SIMPLE BEGIN
        call stop_mpi('ERROR: Unknown TypeKrylov='//KrylovType)
     end select                                    !^CFG END SIMPLE
     call timing_stop('krylov solver')

     if(oktest_me.and.nImplBLK>0)write(*,*)'solution dw(test)=',dw(implVARtest)

     ! Postprocessing: dw = P_R.dw' where P_R = I, U^{-1}, U^{-1}L^{-1} for 
     ! left, symmetric and right preconditioning, respectively
     if(JacobianType=='prec' .and. PrecondSide/='left')then
        do implBLK=1,nImplBLK
           if(PrecondSide=='right') &
                call Lhepta(nIJK,nw,nI,nI*nJ,&
                dw(nwIJK*(implBLK-1)+1) ,&
                MAT(1,1,1,1,1,1,implBLK),&   ! Main diagonal
                MAT(1,1,1,1,1,2,implBLK),&   ! -i
                MAT(1,1,1,1,1,4,implBLK),&   ! -j
                MAT(1,1,1,1,1,6,implBLK))    ! -k
           call Uhepta(.true.,nIJK,nw,nI,nI*nJ,&
                dw(nwIJK*(implBLK-1)+1),   &
                MAT(1,1,1,1,1,3,implBLK),  &   ! +i diagonal
                MAT(1,1,1,1,1,5,implBLK),  &   ! +j
                MAT(1,1,1,1,1,7,implBLK))      ! +k
        end do
        if(oktest_me.and.nImplBLK>0)&
             write(*,*)'final     dw(test)=',dw(implVARtest)
     end if
     if(oktest_me)write(*,*)'After KrylovMatVec,info,KrylovError:',&
          KrylovMatVec,info,KrylovError

     if(oktest_me.and.info/=0)write(*,*) &
          'Advance_Impl warning: no convergence, info:',info

     ! Update w: w_k+1 = w_k + coeff*dw  with coeff=1 or coeff<1 from
     ! backtracking (for steady state only) based on reducing the residual 
     ! ||RES_expl(w_k+1)||<=||RES_expl(w_k)||. 
     ! Also calculates RES_impl=dtexpl*R_low_k+1 and logical converged.
     call impl_newton_update(dwnrm, converged)

     if(oktest_me)write(*,*)'dwnrm, converged=',dwnrm, converged
     if(converged)exit
  enddo

  ! Make the update conservative
  if(UseConservativeImplicit)call impl_newton_conserve

  ! Put back implicit result into the explicit code
  call implicit2explicit(w_k(1:nI,1:nJ,1:nK,:,:))

  ! Make explicit part available again for partially explicit scheme
  if(UsePartImplicit)then
     ! Restore unusedBLK
     iNewDecomposition=mod(iNewDecomposition-3, 10000)
     unusedBLK(1:nBlockMax)=skippedBLK(1:nBlockMax)
     call MPI_ALLGATHER(unusedBLK,      nBLK, MPI_LOGICAL, &
          unusedBlock_BP, nBLK, MPI_LOGICAL, iComm, iError)

  endif

  ! Exchange messages, so ghost cells of all blocks are updated
  call exchange_messages

  if(oktest_me)write(*,*)'Advance_impl: nmatvec=',nmatvec
  if(oktest_me.and.nImplBLK>0)write(*,*)'Advance_impl: new w=',&
       w_k(Itest,Jtest,Ktest,VARtest,implBLKtest)
  if(UseNewton.and.oktest_me)write(*,*)'Final NewtonIter, dwnrm=',&
       NewtonIter, dwnrm

  if(UseUpdateCheckOrig)then
     UseUpdateCheck = .true.

     ! Restore _o_ variables in the implicit blocks
     do implBLK=1,nImplBlk
        iBLK=impl2iBLK(implBLK)
        do iVar=1,nVar
           if(iVar==E_)then
              E_o_BLK(:,:,:,iBLK)=w_prev(:,:,:,E_    ,iBLK)
           else
              StateOld_VCB(iVar,:,:,:,iBLK)=w_prev(:,:,:,iVar,iBLK)
           end if
        end do
        StateOld_VCB(P_,:,:,:,iBLK)     = gm1*(E_o_BLK(:,:,:,iBLK)-cHalf*(&
             (StateOld_VCB(rhoUx_,:,:,:,iBLK)**2                              &
             +StateOld_VCB(rhoUy_,:,:,:,iBLK)**2                              &
             +StateOld_VCB(rhoUz_,:,:,:,iBLK)**2)/                            &
             StateOld_VCB(rho_,:,:,:,iBLK)                                    &
             +StateOld_VCB(Bx_,:,:,:,iBLK)**2                                 &
             +StateOld_VCB(By_,:,:,:,iBLK)**2                                 &
             +StateOld_VCB(Bz_,:,:,:,iBLK)**2)                           )
     end do

     call update_check(nStage)
  endif

end subroutine advance_impl
