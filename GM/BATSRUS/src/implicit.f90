!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

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

  use ModImplicit
  use ModProcMH, ONLY: iComm, nProc
  use ModMain, ONLY: nBlockMax, nBlockExplAll, time_accurate, &
       n_step, time_simulation, dt, UseDtFixed, DtFixed, DtFixedOrig, Cfl, &
       iNewDecomposition, NameThisComp, &
       test_string, iTest, jTest, kTest, BlkTest, ProcTest, VarTest
  use ModVarIndexes, ONLY: Rho_
  use ModMultifluid, ONLY: select_fluid, iFluid, nFluid, iP
  use ModAdvance, ONLY : State_VGB, Energy_GBI, StateOld_VCB, EnergyOld_CBI, &
       time_BlK, tmp1_BLK, iTypeAdvance_B, iTypeAdvance_BP, &
       SkippedBlock_, ExplBlock_, ImplBlock_, UseUpdateCheck, DoFixAxis
  use ModPhysics, ONLY : No2Si_V, UnitT_
  use ModPointImplicit, ONLY: UsePointImplicit
  use ModLinearSolver, ONLY: solve_linear_multiblock
  use ModEnergy, ONLY: calc_old_pressure, calc_old_energy
  use ModImplHypre, ONLY: hypre_initialize
  use ModMessagePass, ONLY: exchange_messages
  use ModResistivity, ONLY: init_impl_resistivity
  use BATL_lib, ONLY: Unused_B, Unused_BP, Xyz_DGB
  use BATL_size, ONLY: j0_, nJp1_, k0_, nKp1_
  use ModMpi

  implicit none

  real, external :: minval_BLK, minval_loc_BLK

  integer :: iw, implBLK, iBLK, KrylovMatVec, info
  integer :: NewtonIter
  integer :: iError, iError1
  real    :: KrylovError, dwnrm, local_wnrm(nw), coef1

  logical:: converged

  real :: TimeSimulationOrig

  logical :: DoTest, DoTestMe, DoTestKrylov, DoTestKrylovMe

  logical :: UseUpdateCheckOrig, UsePointImplicitOrig, DoFixAxisOrig

  real    :: pRhoRelativeMin

  external impl_matvec

  integer :: i, j, k, iBlock, iLoc_I(5)

  character(len=15) :: NameSub = 'MH_advance_impl'
  !----------------------------------------------------------------------------
  NameSub(1:2) = NameThisComp
  call set_oktest('implicit',DoTest,DoTestMe) 
  if(DoTestMe) write(*,*)NameSub,' starting at step=',n_step

  ! Initialize some variables in ModImplicit
  call implicit_init

  ! Get initial iterate from current state
  call explicit2implicit(0,nI+1,j0_,nJp1_,k0_,nKp1_,Impl_VGB)

  if(DoTestMe)write(*,*)NameSub,': nImplBLK=',nImplBLK
  if(DoTestMe.and.nImplBLK>0)write(*,*)NameSub,': Impl_VGB=',&
       Impl_VGB(:,iTest,jTest,kTest,implBLKtest)

  call MPI_allreduce(nimpl,nimpl_total, 1,MPI_INTEGER,MPI_SUM,iComm,iError)
  wnrm(1:nw)=-1.0
  ! Global norm of current w_(k=0) = w_n
  do iw=1,nw
     local_wnrm(iw)=sum(Impl_VGB(iw,1:nI,1:nJ,1:nK,1:nImplBLK)**2)
  end do
  call MPI_allreduce(local_wnrm, wnrm, nw, MPI_REAL, MPI_SUM, iComm,iError)
  wnrm=sqrt(wnrm/(nimpl_total/nw))
  where(wnrm < smalldouble) wnrm =1.0

  if(DoTestMe)write(*,*)NameSub,': nimpltot, wnrm=',nimpl_total,wnrm

  TimeSimulationOrig = Time_Simulation
  UseUpdateCheckOrig = UseUpdateCheck
  UseUpdateCheck     = .false.

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

           if(.not. UseImplicitEnergy) CYCLE
           ! Overwrite pressure with energy
           do iFluid = 1, nFluid
              call select_fluid
              ImplOld_VCB(iP,:,:,:,iBLK) = &
                   Energy_GBI(1:nI,1:nJ,1:nK,iBLK,iFluid)
           end do
        end do
     end if

     if(.not.UsePartImplicit2)then
        ! Select Unused_B = not explicit blocks
        iNewDecomposition=mod(iNewDecomposition+1, 10000)
        Unused_BP(1:nBlockMax,:) = &
             iTypeAdvance_BP(1:nBlockMax,:) /= ExplBlock_
        Unused_B(1:nBlockMax) = Unused_BP(1:nBlockMax,iProc)
     end if

     ! advance explicit blocks, calc timestep 
     if(.not.UseDtFixed)cfl=ExplCfl
     call advance_expl(.true., -1) 

     if(.not.UsePartImplicit2)then
        ! update ghost cells for the implicit blocks to time level n+1
        iNewDecomposition=mod(iNewDecomposition+1, 10000)
        Unused_BP(1:nBlockMax,:) = &
             iTypeAdvance_BP(1:nBlockMax,:) == SkippedBlock_
        Unused_B(1:nBlockMax) = Unused_BP(1:nBlockMax,iProc)
     end if

     call exchange_messages

     ! The implicit scheme is only applied on implicit blocks
     iNewDecomposition=mod(iNewDecomposition+1, 10000)
     Unused_BP(1:nBlockMax,:) = &
          iTypeAdvance_BP(1:nBlockMax,:) /= ImplBlock_
     Unused_B(1:nBlockMax) = Unused_BP(1:nBlockMax,iProc)
  end if

  !\
  ! Advance implicitly treated blocks
  !/

  ! Switch off point implicit scheme while advancing the implicit blocks
  UsePointImplicitOrig = UsePointImplicit
  UsePointImplicit = .false.

  ! Switch off merging the cells around the poles during the implicit solve
  DoFixAxisOrig      = DoFixAxis
  DoFixAxis          = .false.

  ! Use implicit time step
  if(.not.UseDtFixed)Cfl = ImplCfl

  if(UseDtFixed)then
     if(DoTestMe)write(*,*)NameSub,': call getdt_courant'
     call getdt_courant(dtexpl)
     dtexpl = 0.5*dtexpl
     dtcoeff = dt/dtexpl
  else
     if(DoTestMe)write(*,*)NameSub,': no call of getdt_courant'
     dtcoeff = implCFL/0.5
  endif

  if (UseBDF2.and.n_step==n_prev+1) then
     ! For 3 level BDF2 scheme set beta=ImplCoeff if previous state is known
     ImplCoeff = (dt+dt_prev)/(2*dt+dt_prev)
  else
     ImplCoeff = ImplCoeff0
  end if

  ! Advance time to level n+1 in case there is explicit time dependence:
  !   R(U^n+1,t^n+1) = R(U^n,t^n+1) + dR/dU(U^n,t^n+1).(U^n+1 - U^n)
  ! so the Jacobian should be evaliated at t^n+1

  Time_Simulation = TimeSimulationOrig + Dt*No2Si_V(UnitT_)

  if(DoTestMe.and.time_accurate)&
       write(*,*)NameSub,': dtcoeff,dtexpl,dt=',dtcoeff,dtexpl,dt
  if(DoTestMe.and.UseBDF2)write(*,*)NameSub,': n_prev,dt_prev,ImplCoeff=',&
       n_prev,dt_prev,ImplCoeff

  if(.not.UseBDF2)then
     ! Save the current state into ImplOld_VCB so that StateOld_VCB 
     ! can be restored. 
     ! The implicit blocks haven't been updated, so save current state
     do implBLK=1,nImplBlk
        iBLK=impl2iBLK(implBLK)
        ImplOld_VCB(:,:,:,:,iBLK) = Impl_VGB(:,1:nI,1:nJ,1:nK,implBLK)
     end do
  end if

  ! Initialize right hand side and dw. Uses ImplOld_VCB for BDF2 scheme.
  call impl_newton_init

  ! Save previous timestep for 3 level scheme
  if(UseBDF2)then
     n_prev  = n_step
     dt_prev = dt

     ! Save the current state into ImplOld_VCB so that StateOld_VCB 
     ! can be restored. 
     ! The implicit blocks haven't been updated, so save current state
     do implBLK=1,nImplBlk
        iBLK=impl2iBLK(implBLK)
        ImplOld_VCB(:,:,:,:,iBLK) = Impl_VGB(:,1:nI,1:nJ,1:nK,implBLK)
     end do
  endif

  ! Newton-Raphson iteration and iterative linear solver
  dwnrm = bigdouble
  NewtonIter = 0
  do
     NewtonIter = NewtonIter+1;
     if(DoTestMe)write(*,*)NameSub,': NewtonIter=',NewtonIter
     if(NewtonIter > NewtonIterMax)then
        write(*,*)'Newton-Raphson failed to converge NewtonIter=',NewtonIter
        if(time_accurate)call stop_mpi('Newton-Raphson failed to converge')
        exit
     endif
     nnewton=nnewton+1

     ! Calculate Jacobian matrix if required
     if(ImplParam%DoPrecond .and. (NewtonIter==1 .or. NewMatrix))then

        if(NewtonIter>1)then
           ! Update ghost cells for Impl_VGB, 
           ! because it is needed by impl_jacobian
           call implicit2explicit(Impl_VGB(:,1:nI,1:nJ,1:nK,:))
           call exchange_messages
           call explicit2implicit(0,nI+1,j0_,nJp1_,k0_,nKp1_,Impl_VGB)
        end if

        call timing_start('impl_jacobian')

        ! Initialize variables for preconditioner calculation
        call init_impl_resistivity

        ! Calculate approximate dR/dU matrix
        do implBLK = 1, nImplBLK
           call impl_jacobian(implBLK,MAT(1,1,1,1,1,1,implBLK))
        end do
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

     ! solve implicit system

     ! For Newton solver the outer loop has to converge, 
     ! the inner loop only needs to reduce the error somewhat.
     if(UseNewton) ImplParam%KrylovErrorMax = 0.1

     call set_oktest('krylov', DoTestKrylov, DoTestKrylovMe)

     call solve_linear_multiblock(ImplParam, &
          nVar, nDim, nI, nJ, nK, nImplBlk, iComm, &
          impl_matvec, Rhs, Dw, info, DoTestKrylovMe, MAT)

     if(DoTestMe .and. nImplBLK>0)&
          write(*,*)NameSub,': final     dw(test)=',dw(implVARtest)

     if(info /= 0 .and. iProc == 0 .and. time_accurate) then
        call error_report('Krylov solver failure, Krylov error', &
             KrylovError, iError1, .true.)
     end if

     ! Update w: Impl_VGB(k+1) = Impl_VGB(k) + coeff*dw  
     ! with coeff=1 or coeff<1 from backtracking (for steady state only) 
     ! based on reducing the residual 
     ! ||ResExpl_VCB(Impl_VGB+1)|| <= ||ResExpl_VCB(Impl_VGB)||. 
     ! Also calculates ResImpl_VCB=dtexpl*R_loImpl_VGB+1 
     ! and logical converged.
     call impl_newton_update(dwnrm, converged)

     if(DoTestMe.and.UseNewton) &
          write(*,*)NameSub,': dwnrm, converged=',dwnrm, converged

     if(converged) EXIT
  enddo ! Newton iteration

  ! Make the update conservative
  if(UseConservativeImplicit)call impl_newton_conserve

  ! Put back implicit result into the explicit code
  call implicit2explicit(Impl_VGB(:,1:nI,1:nJ,1:nK,:))

  if(DoFixAxisOrig)call fix_axis_cells

  ! Make explicit part available again for partially explicit scheme
  if(UsePartImplicit)then
     ! Restore Unused_B
     if(.not.UsePartImplicit2)then
        iNewDecomposition=mod(iNewDecomposition-3, 10000)
     else
        iNewDecomposition=mod(iNewDecomposition-1, 10000)
     end if
     Unused_BP(1:nBlockMax,:) = &
          iTypeAdvance_BP(1:nBlockMax,:) == SkippedBlock_
     Unused_B(1:nBlockMax) = Unused_BP(1:nBlockMax,iProc)
  endif

  ! Exchange messages, so ghost cells of all blocks are updated
  call exchange_messages

  if(DoTestMe)write(*,*) NameSub,': nmatvec=',nmatvec
  if(DoTestMe.and.nImplBLK>0)write(*,*)NameSub,': new w=',&
       Impl_VGB(VARtest,Itest,Jtest,Ktest,implBLKtest)
  if(UseNewton.and.DoTestMe)write(*,*)NameSub,': final NewtonIter, dwnrm=',&
       NewtonIter, dwnrm

  ! Restore StateOld and EnergyOld in the implicit blocks
  do implBLK=1,nImplBlk
     iBLK=impl2iBLK(implBLK)
     StateOld_VCB(:,:,:,:,iBLK) = ImplOld_VCB(:,:,:,:,iBLK)

     if(UseImplicitEnergy) then
        do iFluid = 1, nFluid
           call select_fluid
           EnergyOld_CBI(:,:,:,iBLK,iFluid) = ImplOld_VCB(iP,:,:,:,iBLK)
        end do
        call calc_old_pressure(iBlk) ! restore StateOld_VCB(P_...)
     else
        call calc_old_energy(iBlk) ! restore EnergyOld_CBI
     end if
  end do

  if(UseUpdateCheckOrig .and. time_accurate .and. UseDtFixed)then

     ! Calculate the largest relative drop in density or pressure
     do iBLK = 1, nBlock
        if(Unused_B(iBLK)) CYCLE
        ! Check p and rho
        tmp1_BLK(1:nI,1:nJ,1:nK,iBLK)=&
             min(State_VGB(P_,1:nI,1:nJ,1:nK,iBLK) / &
             StateOld_VCB(P_,1:nI,1:nJ,1:nK,iBLK), &
             State_VGB(Rho_,1:nI,1:nJ,1:nK,iBLK) / &
             StateOld_VCB(Rho_,1:nI,1:nJ,1:nK,iBLK) )
     end do

     if(index(Test_String, 'updatecheck') > 0)then
        pRhoRelativeMin = minval_loc_BLK(nProc, tmp1_BLK, iLoc_I)
        if(iLoc_I(5) == iProc)then
           i = iLoc_I(1); j = iLoc_I(2); k = iLoc_I(3); iBlock = iLoc_I(4)
           write(*,*) 'pRhoRelativeMin is at i,j,k,iBlock,iProc = ',iLoc_I
           write(*,*) 'x,y,z =', Xyz_DGB(:,i,j,k,iBlock)
           write(*,*) 'RhoOld,pOld=', StateOld_VCB((/Rho_,P_/),i,j,k,iBlock)
           write(*,*) 'RhoNew,pNew=', State_VGB((/Rho_,P_/),i,j,k,iBlock)
           write(*,*) 'pRhoRelativeMin=', pRhoRelativeMin
        end if
     else
        pRhoRelativeMin = minval_BLK(nProc,tmp1_BLK)
     end if
     if(pRhoRelativeMin < RejectStepLevel .or. info/=0)then
        ! Redo step if pressure decreased below RejectStepLevel
        ! or the Krylov iteration failed.
        Dt = 0.0
        ! Do not use previous step in BDF2 scheme
        n_prev = -1
        ! Reset the state variable, the energy and set time_BLK variable to 0
        do iBLK = 1,nBlock
           if(Unused_B(iBLK)) CYCLE
           State_VGB(:,1:nI,1:nJ,1:nK,iBLK)  = StateOld_VCB(:,:,:,:,iBLK)
           Energy_GBI(1:nI,1:nJ,1:nK,iBLK,:) = EnergyOld_CBI(:,:,:,iBLK,:)
           time_BLK(1:nI,1:nJ,1:nK,iBLK)     = 0.0
        end do
        ! Reduce next time step
        DtFixed = RejectStepFactor*DtFixed
        if(index(Test_String, 'updatecheck') > 0) write(*,*) NameSub, &
             ': RejectStepLevel, info, DtFixed=', RejectStepLevel,info, DtFixed
     elseif(pRhoRelativeMin < ReduceStepLevel)then
        ! Reduce next time step if pressure is reduced below ReduceStepLevel
        DtFixed = ReduceStepFactor*DtFixed
        if(index(Test_String, 'updatecheck') > 0) write(*,*) NameSub, &
             ': ReduceStepLevel, DtFixed=', ReduceStepLevel, DtFixed
     elseif(pRhoRelativeMin > IncreaseStepLevel .and. Dt == DtFixed)then
        ! Increase next time step if pressure remained above IncreaseStepLevel
        ! and the last step was taken with DtFixed. Do not exceed DtFixedOrig
        DtFixed = min(DtFixedOrig, DtFixed*IncreaseStepFactor)
        if(index(Test_String, 'updatecheck') > 0) write(*,*) NameSub, &
             ': IncreaseStepLevel, DtFixed=', IncreaseStepLevel, DtFixed
     end if

     if(DoTestMe) write(*,*) NameSub,': pRelMin,Dt,DtFixed=',&
          pRhoRelativeMin,Dt*No2Si_V(UnitT_), DtFixed*No2Si_V(UnitT_)
  endif

  ! Advance time by Dt
  Time_Simulation = TimeSimulationOrig + Dt*No2Si_V(UnitT_)

  ! Restore logicals
  UseUpdateCheck   = UseUpdateCheckOrig
  UsePointImplicit = UsePointImplicitOrig
  DoFixAxis        = DoFixAxisOrig

end subroutine advance_impl

  
