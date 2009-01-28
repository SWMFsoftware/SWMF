!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT

!=============================================================================
subroutine advance_impl

  ! The implicit schemes used in this module are described in detail in
  !
  ! G. Toth, D. L. De Zeeuw, T. I. Gombosi, K. G. Powell, 2006,
  !  Journal of Computational Physics, 217, 722-758, 
  !  doi:10.1016/j.jcp.2006.01.029
  !
  ! and
  !
  ! Keppens, Toth, Botchev, van der Ploeg, 
  ! J. Int. Num. Methods in Fluids, 30, 335-352, 1999
  !
  ! Equation numbers below refer to the latter paper unless stated otherwise.
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
  use ModMultifluid, ONLY: select_fluid, iFluid, nFluid, iP
  use ModAdvance, ONLY : State_VGB, Energy_GBI, StateOld_VCB, EnergyOld_CBI, &
       time_BlK, &
       tmp1_BLK, UseUpdateCheck, iTypeAdvance_B, iTypeAdvance_BP, &
       SkippedBlock_, ExplBlock_, ImplBlock_
  use ModPhysics, ONLY : No2Si_V, UnitT_
  use ModImplicit
  use ModPointImplicit, ONLY: UsePointImplicit
  use ModAMR, ONLY : UnusedBlock_BP
  use ModNumConst
  use ModLinearSolver, ONLY: gmres, bicgstab, cg, prehepta, Uhepta, Lhepta
  use ModMpi
  use ModEnergy, ONLY: calc_old_pressure

  implicit none

  real, external :: minval_BLK

  integer :: iVar, iw, implBLK, iBLK, KrylovMatVec, info, NewtonIter
  integer :: iError
  real    :: KrylovError, dwnrm, local_wnrm(nw), coef1, coef2

  logical:: converged
  character (LEN=3) :: typestop

  real :: TimeSimulationOrig

  logical :: DoTest, DoTestMe, DoTestKrylov, DoTestKrylovMe

  logical :: UseUpdateCheckOrig, UsePointImplicitOrig

  real    :: pRhoRelativeMin

  character(len=15) :: NameSub = 'MH_advance_impl'

  external impl_matvec

  !----------------------------------------------------------------------------

  NameSub(1:2) = NameThisComp
  call set_oktest('implicit',DoTest,DoTestMe) 
  if(DoTestMe) write(*,*)NameSub,' starting at step=',n_step

  ! Initialize some variables in ModImplicit
  call implicit_init

  ! Get initial iterate from current state
  call explicit2implicit(0,nI+1,0,nJ+1,0,nK+1,Impl_VGB)

  if(DoTestMe)write(*,*)NameSub,': nImplBLK=',nImplBLK
  if(DoTestMe.and.nImplBLK>0)write(*,*)NameSub,': Impl_VGB=',&
       Impl_VGB(:,iTest,jTest,kTest,implBLKtest)

  call MPI_allreduce(nimpl,nimpl_total, 1,MPI_INTEGER,MPI_SUM,iComm,iError)
  if(UseSemiImplicit)then
     wnrm = 1.0
  else
     wnrm(1:nw)=-1.0
     ! Global norm of current w_(k=0) = w_n
     do iw=1,nw
        local_wnrm(iw)=sum(Impl_VGB(iw,1:nI,1:nJ,1:nK,1:nImplBLK)**2)
     end do
     call MPI_allreduce(local_wnrm, wnrm, nw, MPI_REAL, MPI_SUM, iComm,iError)
     wnrm=sqrt(wnrm/(nimpl_total/nw))
     where(wnrm < smalldouble) wnrm =1.0
  end if

  if(DoTestMe)write(*,*)NameSub,': nimpltot, wnrm=',nimpl_total,wnrm

  TimeSimulationOrig = Time_Simulation
  UseUpdateCheckOrig = UseUpdateCheck
  UseUpdateCheck = .false.

  ! Advance explicitly treated blocks if any
  if(UsePartImplicit .and. nBlockExplALL > 0)then

     if(DoTestMe)write(*,*)NameSub,': advance explicit blocks'

     if(UseBDF2)then

        ! Save the current state into the previous state for BDF2 scheme
        ! This is needed for explicit blocks, because they may become
        ! implicit in the next time step...
        do iBLK=1,nBlock
           if(iTypeAdvance_B(iBLK) /= ExplBlock_)CYCLE
           ImplOld_VCB(:,:,:,:,iBLK) = State_VGB(:,1:nI,1:nJ,1:nK,iBLK)
           ! Overwrite pressure with energy
           do iFluid = 1, nFluid
              call select_fluid
              ImplOld_VCB(iP,:,:,:,iBLK) = &
                   Energy_GBI(1:nI,1:nJ,1:nK,iBLK,iFluid)
           end do
        end do
     end if

     if(.not.UsePartImplicit2)then
        ! Select unusedBLK = not explicit blocks
        iNewDecomposition=mod(iNewDecomposition+1, 10000)
        UnusedBlock_BP(1:nBlockMax,:) = &
             iTypeAdvance_BP(1:nBlockMax,:) /= ExplBlock_
        UnusedBLK(1:nBlockMax) = UnusedBlock_BP(1:nBlockMax,iProc)
     end if

     ! advance explicit blocks, calc timestep 
     if(.not.UseDtFixed)cfl=ExplCfl
     call advance_expl(.true.) 

     if(.not.UsePartImplicit2)then
        ! update ghost cells for the implicit blocks to time level n+1
        iNewDecomposition=mod(iNewDecomposition+1, 10000)
        UnusedBlock_BP(1:nBlockMax,:) = &
             iTypeAdvance_BP(1:nBlockMax,:) == SkippedBlock_
        UnusedBLK(1:nBlockMax) = UnusedBlock_BP(1:nBlockMax,iProc)
     end if

     call exchange_messages

     ! The implicit scheme is only applied on implicit blocks
     iNewDecomposition=mod(iNewDecomposition+1, 10000)
     UnusedBlock_BP(1:nBlockMax,:) = &
          iTypeAdvance_BP(1:nBlockMax,:) /= ImplBlock_
     UnusedBLK(1:nBlockMax) = UnusedBlock_BP(1:nBlockMax,iProc)
  end if

  !\
  ! Advance implicitly treated blocks
  !/

  ! Switch off point implicit scheme while advancing the implicit blocks
  UsePointImplicitOrig = UsePointImplicit
  UsePointImplicit = .false.

  if(UseSemiImplicit)then
     ! time step is set by the explicit scheme
     dtexpl = dt
     dtcoeff = 1.0
     ImplCfl = Cfl
  else
     ! Use implicit time step
     if(.not.UseDtFixed)Cfl = ImplCfl

     if(UseDtFixed)then
        if(DoTestMe)write(*,*)NameSub,': call getdt_courant'
        call getdt_courant(dtexpl)
        dtexpl=cHalf*dtexpl
        dtcoeff=dt/dtexpl
     else
        if(DoTestMe)write(*,*)NameSub,': no call of getdt_courant'
        dtcoeff=implCFL/cHalf
     endif
  end if

  if (UseBDF2.and.n_step==n_prev+1) then
     ! For 3 level BDF2 scheme set beta=ImplCoeff if previous state is known
     ImplCoeff = (dt+dt_prev)/(2*dt+dt_prev)
  else
     ImplCoeff = ImplCoeff0
  end if

  ! Advance time to level n+1 in case there is explicit time dependence:
  !   R(U^n+1,t^n+1) = R(U^n,t^n+1) + dR/dU(U^n,t^n+1).(U^n+1 - U^n)
  ! so the Jacobian should be evaliated at t^n+1
  ! Semi-implicit scheme has already advanced time in advance_expl
  if(.not.UseSemiImplicit) &
       Time_Simulation = TimeSimulationOrig + Dt*No2Si_V(UnitT_)

  if(DoTestMe.and.time_accurate)&
       write(*,*)NameSub,': dtcoeff,dtexpl,dt=',dtcoeff,dtexpl,dt
  if(DoTestMe.and.UseBDF2)write(*,*)NameSub,': n_prev,dt_prev,ImplCoeff=',&
       n_prev,dt_prev,ImplCoeff

  ! Initialize right hand side and dw
  call impl_newton_init

  ! Save previous timestep for 3 level scheme
  if(UseBDF2)then
     n_prev  = n_step
     dt_prev = dt
  endif

  ! Save the current state into ImplOld_VCB so that StateOld_VCB can be restored
  ! The implicit blocks haven't been updated, so save current state
  do implBLK=1,nImplBlk
     iBLK=impl2iBLK(implBLK)
     ImplOld_VCB(:,:,:,:,iBLK) = Impl_VGB(:,1:nI,1:nJ,1:nK,implBLK)
  end do

  ! Newton-Raphson iteration and iterative linear solver
  dwnrm=bigdouble
  NewtonIter=0
  do
     NewtonIter=NewtonIter+1;
     if(DoTestMe)write(*,*)NameSub,': NewtonIter=',NewtonIter
     if(NewtonIter>NewtonIterMax)then
        write(*,*)'Newton-Raphson failed to converge NewtonIter=',NewtonIter
        if(time_accurate)call stop_mpi('Newton-Raphson failed to converge')
        exit
     endif
     nnewton=nnewton+1

     ! Calculate Jacobian matrix if required
     if(JacobianType/='free'.and.(NewtonIter==1.or.NewMatrix))then

        !!! need to be changed for semi-implicit
        if(NewtonIter>1)then
           ! Update ghost cells for Impl_VGB, 
           ! because it is needed by impl_jacobian
           call implicit2explicit(Impl_VGB(:,1:nI,1:nJ,1:nK,:))
           call exchange_messages
           call explicit2implicit(0,nI+1,0,nJ+1,0,nK+1,Impl_VGB)
        end if

        call timing_start('impl_jacobian')
        if(UseSemiImplicit)then
           call get_semi_impl_jacobian
        else
           do implBLK=1,nImplBLK
              call impl_jacobian(implBLK,MAT(1,1,1,1,1,1,implBLK))
           end do
        end if
        call timing_stop('impl_jacobian')

        if(DoTest)then
           call MPI_reduce(sum(MAT(:,:,:,:,:,:,1:nImplBLK)**2),coef1,1,&
                MPI_REAL,MPI_SUM,PROCtest,iComm,iError)
           if(DoTestMe)write(*,*)NameSub,': sum(MAT**2)=',coef1
        end if

     endif

     ! Update rhs and initial dw if required
     if (NewtonIter>1) call impl_newton_loop

     if(DoTestMe.and.nImplBLK>0)write(*,*)NameSub,&
          ': initial dw(test), rhs(test)=',dw(implVARtest),rhs(implVARtest)

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

        if(DoTest)then
           call MPI_reduce(sum(MAT(:,:,:,:,:,:,1:nImplBLK)**2),coef1,1,&
                MPI_REAL,MPI_SUM,procTEST,iComm,iError)
           call MPI_reduce(sum(rhs(1:nimpl)**2),coef2,1,MPI_REAL,&
                MPI_SUM,procTEST,iComm,iError)
           if(DoTestMe)then
              write(*,*)NameSub,': preconditioned sum(MAT**2), sum(rhs**2)=',&
                   coef1,coef2
              if(nImplBLK>0)&
                   write(*,*)NameSub,&
                   ': preconditioned dw(test)   , rhs(test)  =',&
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

     if(DoTestMe)write(*,*)NameSub,': Before ',KrylovType,&
          ' KrylovMatVec,KrylovError:',KrylovMatVec,KrylovError

     ! Solve linear problem

     call set_oktest('krylov',DoTestKrylov,DoTestKrylovMe)
     call timing_start('krylov solver')
     select case(KrylovType)
     case('bicgstab','BiCGSTAB')
        call bicgstab(impl_matvec,rhs,dw,non0dw,nimpl,&
             KrylovError,typestop,KrylovMatVec,info,DoTestKrylovMe,iComm)
     case('GMRES','gmres')
        call gmres(impl_matvec,rhs,dw,non0dw,nimpl,nKrylovVector, &
             KrylovError,typestop,KrylovMatVec,info,DoTestKrylovMe,iComm)
     case('CG','cg')
        call cg(impl_matvec,rhs,dw,non0dw,nimpl,&
             KrylovError,typestop,KrylovMatVec,info,DoTestKrylovMe,iComm)        
     case default
        call stop_mpi('ERROR: Unknown TypeKrylov='//KrylovType)
     end select
     call timing_stop('krylov solver')

     if(DoTestMe.and.nImplBLK>0)write(*,*)NameSub,&
          ': solution dw(test)=',dw(implVARtest)

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
        if(DoTestMe.and.nImplBLK>0)&
             write(*,*)NameSub,': final     dw(test)=',dw(implVARtest)
     end if
     if(DoTestMe.or.DoTestKrylovMe)write(*,*)NameSub,&
          ': After KrylovMatVec,info,KrylovError=',&
          KrylovMatVec,info,KrylovError

     if(DoTestMe.and.info/=0)write(*,*) NameSub, &
          ' warning: no convergence, info:',info

     ! Update w: Impl_VGB(k+1) = Impl_VGB(k) + coeff*dw  
     ! with coeff=1 or coeff<1 from backtracking (for steady state only) 
     ! based on reducing the residual 
     ! ||ResExpl_VCB(Impl_VGB+1)|| <= ||ResExpl_VCB(Impl_VGB)||. 
     ! Also calculates ResImpl_VCB=dtexpl*R_loImpl_VGB+1 and logical converged.
     call impl_newton_update(dwnrm, converged)

     if(DoTestMe.and.UseNewton) &
          write(*,*)NameSub,': dwnrm, converged=',dwnrm, converged
     if(converged)exit
  enddo

  ! Make the update conservative
  if(UseConservativeImplicit)call impl_newton_conserve

  ! Put back implicit result into the explicit code
  call implicit2explicit(Impl_VGB(:,1:nI,1:nJ,1:nK,:))

  ! Make explicit part available again for partially explicit scheme
  if(UsePartImplicit)then
     ! Restore unusedBLK
     if(.not.UsePartImplicit2)then
        iNewDecomposition=mod(iNewDecomposition-3, 10000)
     else
        iNewDecomposition=mod(iNewDecomposition-1, 10000)
     end if
     UnusedBlock_BP(1:nBlockMax,:) = &
          iTypeAdvance_BP(1:nBlockMax,:) == SkippedBlock_
     UnusedBLK(1:nBlockMax) = UnusedBlock_BP(1:nBlockMax,iProc)
  endif

  ! Exchange messages, so ghost cells of all blocks are updated
  call exchange_messages

  if(DoTestMe)write(*,*) NameSub,': nmatvec=',nmatvec
  if(DoTestMe.and.nImplBLK>0)write(*,*)NameSub,': new w=',&
       Impl_VGB(VARtest,Itest,Jtest,Ktest,implBLKtest)
  if(UseNewton.and.DoTestMe)write(*,*)NameSub,': final NewtonIter, dwnrm=',&
       NewtonIter, dwnrm

  if(.not.UseSemiImplicit)then
     ! Restore StateOld and E_o_BLK in the implicit blocks
     do implBLK=1,nImplBlk
        iBLK=impl2iBLK(implBLK)
        StateOld_VCB(:,:,:,:,iBLK)=ImplOld_VCB(:,:,:,:,iBLK)
        do iFluid = 1, nFluid
           call select_fluid
           EnergyOld_CBI(:,:,:,iBLK,iFluid) = ImplOld_VCB(iP,:,:,:,iBLK)
        end do
        call calc_old_pressure(iBlk) ! restore StateOld_VCB(P_...)
     end do
  end if

  if(UseUpdateCheckOrig .and. time_accurate .and. UseDtFixed)then
     
     ! Calculate the largest relative drop in density or pressure
     do iBLK=1,nBlock
        if(UnusedBlk(iBLK)) CYCLE
        tmp1_BLK(1:nI,1:nJ,1:nK,iBLK)=&
             min( &
             State_VGB(P_,1:nI,1:nJ,1:nK,iBLK) / &
             StateOld_VCB(P_,1:nI,1:nJ,1:nK,iBLK), &
             State_VGB(Rho_,1:nI,1:nJ,1:nK,iBLK) / &
             StateOld_VCB(Rho_,1:nI,1:nJ,1:nK,iBLK) &
             )
     end do
     pRhoRelativeMin = minval_BLK(nProc,tmp1_BLK)

     if(pRhoRelativeMin < RejectStepLevel .or. info/=0)then
        ! Redo step if pressure decreased below RejectStepLevel
        ! or the Krylov iteration failed.
        Dt = 0.0
        ! Do not use previous step in BDF2 scheme
        n_prev = -1
        ! Reset the state variable, the energy and set time_BLK variable to 0
        do iBLK = 1,nBlock
           if(UnusedBlk(iBLK)) CYCLE
           State_VGB(:,1:nI,1:nJ,1:nK,iBLK)  = StateOld_VCB(:,:,:,:,iBLK)
           Energy_GBI(1:nI,1:nJ,1:nK,iBLK,:) = EnergyOld_CBI(:,:,:,iBLK,:)
           time_BLK(1:nI,1:nJ,1:nK,iBLK)     = 0.0
        end do
        ! Reduce next time step
        DtFixed = RejectStepFactor*DtFixed
     elseif(pRhoRelativeMin < ReduceStepLevel)then
        ! Reduce next time step if pressure is reduced below ReduceStepLevel
        DtFixed = ReduceStepFactor*DtFixed
     elseif(pRhoRelativeMin > IncreaseStepLevel .and. Dt == DtFixed)then
        ! Increase next time step if pressure remained above IncreaseStepLevel
        ! and the last step was taken with DtFixed. Do not exceed DtFixedOrig
        DtFixed = min(DtFixedOrig, DtFixed*IncreaseStepFactor)
     end if

     if(DoTestMe) write(*,*) NameSub,': pRelMin,Dt,DtFixed=',&
          pRhoRelativeMin,Dt*No2Si_V(UnitT_), DtFixed*No2Si_V(UnitT_)
  endif

  ! Advance time by Dt
  if(.not.UseSemiImplicit) &
       Time_Simulation = TimeSimulationOrig + Dt*No2Si_V(UnitT_)

  UseUpdateCheck   = UseUpdateCheckOrig
  UsePointImplicit = UsePointImplicitOrig

end subroutine advance_impl
