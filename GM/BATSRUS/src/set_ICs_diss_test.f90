!^CFG COPYRIGHT UM
!^CFG FILE DISSFLUX
subroutine set_ICs_diss_test
  use ModPhysics,  ONLY:TypeProblemDiss
  implicit none
  
  select case(TypeProblemDiss)
  case('heat_test1')
     call heat_test_1
  case('heat_test2')
     call heat_test_2
  case('heat_test3')
     call heat_test_3
  case('heat_test4')
     call heat_test_4
  case('resist_test1')
     call resist_test_1
  case('mag_rec1')
     call mag_rec_1
  case('mag_rec2')
     call mag_rec_2
  case('mag_rec3')
     call mag_rec_3
  case('loop_rec1')
     call loop_rec_1
  case('spicules1')
     call spicules_1
  case('spicules2')
     call spicules_2
  case('spicules3')
     call spicules_3
  case default
     call stop_mpi('Unknown TypeProblemDiss='//TypeProblemDiss)
  end select
end subroutine set_ICs_diss_test
!***********************************************************************
subroutine heat_test_1
  use ModProcMH,   ONLY:iProc 
  use ModMain,     ONLY:nI,nJ,nK,gcn,globalBLK,nBLK,&
       BLKtest,unusedBLK,UseGravity
  use ModVarIndexes,ONLY: rho_,rhoUx_,rhoUy_,&
       rhoUz_,Bx_,By_,Bz_,P_
  use ModGeometry, ONLY:x_BLK,y_BLK,x1,x2,y1,y2
  use ModAdvance,  ONLY:State_VGB
  use ModPhysics,  ONLY:TypeProblemDiss,EpsilonDiss,DeltaDiss, &
       ThetaDiss,RhoDifDiss,yShiftDiss,BZ0Diss,scaleHeightDiss, &
       scaleFactorDiss,Gbody,Grav0Diss,Beta0Diss,Length0Diss, &
       Time0Diss,Rho0Diss,Tem0Diss
  use ModNumConst, ONLY:cOne,cTwo,cFour,cHalf, &
       cZero,cTiny,cHundred,cPi
  implicit none

  integer:: i,j,k, iBLK
  real:: yyy
  !\
  ! This experiment is aimed to test heat conduction: Test1
  !/
  do iBLK=1,nBLK
     if (unusedBLK(iBLK)) CYCLE
     do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
        yyy = y_BLK(i,j,k,iBLK)/y2
        State_VGB(rho_,i,j,k,iBLK)   = cOne+ &
             RhoDifDiss*(cOne-tanh(cTwo*(cOne+cFour)*yyy)) 
        State_VGB(rhoUx_,i,j,k,iBLK) = cZero
        State_VGB(rhoUy_,i,j,k,iBLK) = cZero
        State_VGB(rhoUz_,i,j,k,iBLK) = cZero
        State_VGB(Bx_,i,j,k,iBLK)    = cZero
        State_VGB(By_,i,j,k,iBLK)    = cOne
        State_VGB(Bz_,i,j,k,iBLK)    = cZero
        State_VGB(P_,i,j,k,iBLK)     = Beta0Diss/cTwo
     end do; end do; end do
  end do
  if(globalBLK==BLKtest.and.iProc==0)then
     write(*,*) ''
     write(*,*) '>>>>>>>>>>ILIAs heat_test_1<<<<<<<<<<'
     write(*,*) ''
     write(*,*) 'Grav0Diss   = ',Grav0Diss
     write(*,*) 'Beta0Diss   = ',Beta0Diss
     write(*,*) 'Length0Diss = ',Length0Diss
     write(*,*) 'Time0Diss   = ',Time0Diss
     write(*,*) 'Rho0Diss    = ',Rho0Diss
     write(*,*) 'Tem0Diss    = ',Tem0Diss
     write(*,*) ''
     write(*,*) 'ThetaDiss   = ',ThetaDiss
     write(*,*) 'DeltaDiss   = ',DeltaDiss
     write(*,*) 'EpsilonDiss = ',EpsilonDiss
     write(*,*) 'RhoDifDiss  = ',RhoDifDiss
     write(*,*) 'yShiftDiss  = ',yShiftDiss
     write(*,*) 'BZ0Diss     = ',BZ0Diss
     write(*,*) 'scaleHeightDiss = ',scaleHeightDiss
     write(*,*) 'scaleFactorDiss = ',scaleFactorDiss
     write(*,*) 'BZ0Diss  =',BZ0Diss
     write(*,*) 'min(Rho),max(Rho) = ',minval(State_VGB(rho_,:,:,:,:)),&
          maxval(State_VGB(rho_,:,:,:,:))
     write(*,*) 'min(p)  ,max(p)   = ',minval(State_VGB(P_,:,:,:,:)),&
          maxval(State_VGB(P_,:,:,:,:))
     write(*,*) 'min(Px) ,max(Px)  = ',minval(State_VGB(rhoUx_,:,:,:,:)),&
          maxval(State_VGB(rhoUx_,:,:,:,:))
     write(*,*) 'min(Py) ,max(Py)  = ',minval(State_VGB(rhoUy_,:,:,:,:)),&
          maxval(State_VGB(rhoUy_,:,:,:,:))
     write(*,*) 'min(Pz) ,max(Pz)  = ',minval(State_VGB(rhoUz_,:,:,:,:)),&
          maxval(State_VGB(rhoUz_,:,:,:,:))
     write(*,*) 'min(Bx) ,max(Bx)  = ',minval(State_VGB(Bx_,:,:,:,:)),&
          maxval(State_VGB(Bx_,:,:,:,:))
     write(*,*) 'min(By) ,max(By)  = ',minval(State_VGB(By_,:,:,:,:)),&
          maxval(State_VGB(By_,:,:,:,:))
     write(*,*) 'min(Bz) ,max(Bz)  = ',minval(State_VGB(Bz_,:,:,:,:)),&
          maxval(State_VGB(Bz_,:,:,:,:))
     write(*,*) ''
     write(*,*) '>>>>>>>>>>                 <<<<<<<<<<'
     write(*,*) ''
  endif
end subroutine heat_test_1
!***********************************************************************
subroutine heat_test_2
  use ModProcMH,   ONLY:iProc 
  use ModMain,     ONLY:nI,nJ,nK,gcn,globalBLK,nBLK,&
       BLKtest,unusedBLK,UseGravity
  use ModVarIndexes,ONLY:rho_,P_,rhoUx_,rhoUy_,&
       rhoUz_,Bx_,By_,Bz_
  use ModGeometry, ONLY:x_BLK,y_BLK,x1,x2,y1,y2
  use ModAdvance,  ONLY:State_VGB
  use ModPhysics,  ONLY:TypeProblemDiss,EpsilonDiss,DeltaDiss, &
       ThetaDiss,RhoDifDiss,yShiftDiss,BZ0Diss,scaleHeightDiss, &
       scaleFactorDiss,Gbody,Grav0Diss,Beta0Diss,Length0Diss, &
       Time0Diss,Rho0Diss,Tem0Diss
  use ModNumConst, ONLY:cOne,cTwo,cFour,cHalf, &
       cZero,cTiny,cHundred,cPi
  implicit none
  
  integer:: i,j,k, iBLK
  real:: xxx,yyy
  !\
  ! This experiment is aimed to test heat conduction: Test2
  !/
  do iBLK=1,nBLK
     if (unusedBLK(iBLK)) CYCLE
     do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
        xxx = x_BLK(i,j,k,iBLK)
        yyy = y_BLK(i,j,k,iBLK)+yShiftDiss*y2
        !
        State_VGB(rhoUx_,i,j,k,iBLK) = cZero
        State_VGB(rhoUy_,i,j,k,iBLK) = cZero
        State_VGB(rhoUz_,i,j,k,iBLK) = cZero
        State_VGB(Bx_,i,j,k,iBLK)    = cZero
        State_VGB(By_,i,j,k,iBLK)    = tanh(DeltaDiss*xxx)
        State_VGB(Bz_,i,j,k,iBLK)    = cZero
        State_VGB(rho_,i,j,k,iBLK)   = (cOne+RhoDifDiss*(cOne-tanh(cTwo*     &
                                (cOne+cFour)*yyy)))*((cOne+Beta0Diss- &
                                State_VGB(By_,i,j,k,iBLK)**2)/Beta0Diss)**   &
                                ThetaDiss
        State_VGB(P_,i,j,k,iBLK)     = cHalf*(cOne+Beta0Diss-State_VGB(By_,i,j,k,iBLK)**2)
     end do; end do; end do
  end do
  if(globalBLK==BLKtest.and.iProc==0)then
     write(*,*) ''
     write(*,*) '>>>>>>>>>>ILIAs heat_test_2<<<<<<<<<<'
     write(*,*) ''
     write(*,*) 'Grav0Diss   = ',Grav0Diss
     write(*,*) 'Beta0Diss   = ',Beta0Diss
     write(*,*) 'Length0Diss = ',Length0Diss
     write(*,*) 'Time0Diss   = ',Time0Diss
     write(*,*) 'Rho0Diss    = ',Rho0Diss
     write(*,*) 'Tem0Diss    = ',Tem0Diss
     write(*,*) ''
     write(*,*) 'ThetaDiss   = ',ThetaDiss
     write(*,*) 'DeltaDiss   = ',DeltaDiss
     write(*,*) 'EpsilonDiss = ',EpsilonDiss
     write(*,*) 'RhoDifDiss  = ',RhoDifDiss
     write(*,*) 'yShiftDiss  = ',yShiftDiss
     write(*,*) 'BZ0Diss     = ',BZ0Diss
     write(*,*) 'scaleHeightDiss = ',scaleHeightDiss
     write(*,*) 'scaleFactorDiss = ',scaleFactorDiss
     write(*,*) 'min(Rho),max(Rho) = ',minval(State_VGB(rho_,:,:,:,:))&
          ,maxval(State_VGB(rho_,:,:,:,:))
     write(*,*) 'min(p)  ,max(p)   = ',minval(State_VGB(P_,:,:,:,:)),maxval(State_VGB(P_,:,:,:,:))
     write(*,*) 'min(Px) ,max(Px)  = ',minval(State_VGB(rhoUx_,:,:,:,:)),maxval(State_VGB(rhoUx_,:,:,:,:))
     write(*,*) 'min(Py) ,max(Py)  = ',minval(State_VGB(rhoUy_,:,:,:,:)),maxval(State_VGB(rhoUy_,:,:,:,:))
     write(*,*) 'min(Pz) ,max(Pz)  = ',minval(State_VGB(rhoUz_,:,:,:,:)),maxval(State_VGB(rhoUz_,:,:,:,:))
     write(*,*) 'min(Bx) ,max(Bx)  = ',minval(State_VGB(Bx_,:,:,:,:)),maxval(State_VGB(Bx_,:,:,:,:))
     write(*,*) 'min(By) ,max(By)  = ',minval(State_VGB(By_,:,:,:,:)),maxval(State_VGB(By_,:,:,:,:))
     write(*,*) 'min(Bz) ,max(Bz)  = ',minval(State_VGB(Bz_,:,:,:,:)),maxval(State_VGB(Bz_,:,:,:,:))
     write(*,*) ''
     write(*,*) '>>>>>>>>>>                 <<<<<<<<<<'
     write(*,*) ''
  endif
end subroutine heat_test_2
!***********************************************************************
subroutine heat_test_3
  use ModProcMH,   ONLY:iProc 
  use ModMain,     ONLY:nI,nJ,nK,gcn,globalBLK,nBLK,&
       BLKtest,unusedBLK,UseGravity
  use ModVarIndexes,ONLY:rho_,&
       rhoUx_,rhoUy_,rhoUz_,Bx_,By_,Bz_,P_
  use ModGeometry, ONLY:x_BLK,y_BLK,x1,x2,y1,y2
  use ModAdvance,  ONLY:State_VGB
  use ModPhysics,  ONLY:TypeProblemDiss,EpsilonDiss,DeltaDiss, &
       ThetaDiss,RhoDifDiss,yShiftDiss,BZ0Diss,scaleHeightDiss, &
       scaleFactorDiss,Gbody,Grav0Diss,Beta0Diss,Length0Diss, &
       Time0Diss,Rho0Diss,Tem0Diss
  use ModNumConst, ONLY:cOne,cTwo,cFour,cHalf, &
       cZero,cTiny,cHundred,cPi
  implicit none
  
  integer:: i,j,k, iBLK
  real:: yyy
  !\
  ! This experiment is aimed to test heat conduction: Test3
  ! The exact solution of the problem is known:::
  ! T = To+Tdif*{1-cos[0.5*cPi*(y+y2)/y2]*exp[-Kappao*cPi^2*t/(4*y2^2)]}
  !/
  do iBLK=1,nBLK
     if (unusedBLK(iBLK)) CYCLE
     do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
        yyy = y_BLK(i,j,k,iBLK)/y2
        !
        State_VGB(rho_,i,j,k,iBLK)   = cOne
        State_VGB(rhoUx_,i,j,k,iBLK) = cZero
        State_VGB(rhoUy_,i,j,k,iBLK) = cZero
        State_VGB(rhoUz_,i,j,k,iBLK) = cZero
        State_VGB(Bx_,i,j,k,iBLK)    = cZero
        State_VGB(By_,i,j,k,iBLK)    = cOne
        State_VGB(Bz_,i,j,k,iBLK)    = cZero
        State_VGB(P_,i,j,k,iBLK)     = (cOne+RhoDifDiss*(cOne- &
             cos(cHalf*cPi*(yyy+cOne))))/(cTwo*(cOne+cFour))
     end do; end do; end do
  end do
  if(globalBLK==BLKtest.and.iProc==0)then
     write(*,*) ''
     write(*,*) '>>>>>>>>>>ILIAs heat_test_3<<<<<<<<<<'
     write(*,*) ''
     write(*,*) 'Grav0Diss   = ',Grav0Diss
     write(*,*) 'Beta0Diss   = ',Beta0Diss
     write(*,*) 'Length0Diss = ',Length0Diss
     write(*,*) 'Time0Diss   = ',Time0Diss
     write(*,*) 'Rho0Diss    = ',Rho0Diss
     write(*,*) 'Tem0Diss    = ',Tem0Diss
     write(*,*) ''
     write(*,*) 'ThetaDiss   = ',ThetaDiss
     write(*,*) 'DeltaDiss   = ',DeltaDiss
     write(*,*) 'EpsilonDiss = ',EpsilonDiss
     write(*,*) 'RhoDifDiss  = ',RhoDifDiss
     write(*,*) 'yShiftDiss  = ',yShiftDiss
     write(*,*) 'BZ0Diss     = ',BZ0Diss
     write(*,*) 'scaleHeightDiss = ',scaleHeightDiss
     write(*,*) 'scaleFactorDiss = ',scaleFactorDiss
     write(*,*) 'min(Rho),max(Rho) = ',minval(State_VGB(rho_,:,:,:,:)),maxval(State_VGB(rho_,:,:,:,:))
     write(*,*) 'min(p)  ,max(p)   = ',minval(State_VGB(P_,:,:,:,:)),maxval(State_VGB(P_,:,:,:,:))
     write(*,*) 'min(Px) ,max(Px)  = ',minval(State_VGB(rhoUx_,:,:,:,:)),maxval(State_VGB(rhoUx_,:,:,:,:))
     write(*,*) 'min(Py) ,max(Py)  = ',minval(State_VGB(rhoUy_,:,:,:,:)),maxval(State_VGB(rhoUy_,:,:,:,:))
     write(*,*) 'min(Pz) ,max(Pz)  = ',minval(State_VGB(rhoUz_,:,:,:,:)),maxval(State_VGB(rhoUz_,:,:,:,:))
     write(*,*) 'min(Bx) ,max(Bx)  = ',minval(State_VGB(Bx_,:,:,:,:)),maxval(State_VGB(Bx_,:,:,:,:))
     write(*,*) 'min(By) ,max(By)  = ',minval(State_VGB(By_,:,:,:,:)),maxval(State_VGB(By_,:,:,:,:))
     write(*,*) 'min(Bz) ,max(Bz)  = ',minval(State_VGB(Bz_,:,:,:,:)),maxval(State_VGB(Bz_,:,:,:,:))
     write(*,*) ''
     write(*,*) '>>>>>>>>>>                 <<<<<<<<<<'
     write(*,*) ''
  endif
end subroutine heat_test_3
!***********************************************************************
subroutine heat_test_4
  use ModProcMH,   ONLY:iProc
  use ModMain,     ONLY:nI,nJ,nK,gcn,globalBLK,&
       nBLK,BLKtest,unusedBLK,UseGravity
  use ModVarIndexes,ONLY:rho_,&
       rhoUx_,rhoUy_,rhoUz_,Bx_, By_, Bz_, P_
  use ModGeometry, ONLY:x_BLK,y_BLK,x1,x2,y1,y2
  use ModAdvance,  ONLY:State_VGB
  use ModPhysics,  ONLY:TypeProblemDiss,EpsilonDiss,DeltaDiss, &
       ThetaDiss,RhoDifDiss,yShiftDiss,BZ0Diss,scaleHeightDiss, &
       scaleFactorDiss,Gbody,Grav0Diss,Beta0Diss,Length0Diss, &
       Time0Diss,Rho0Diss,Tem0Diss,Kappa0Heat,ExponentHeat
  use ModNumConst, ONLY:cOne,cTwo,cFour,cHalf, &
       cZero,cTiny,cHundred,cPi
  use ModConst,    ONLY:cBoltzmann,cProtonMass
  use ModMpi
  implicit none
  
  integer:: i,j,k
  real:: xxx,yyy
  real:: TExact0Diss,TExactDiss
  real:: t0Diss,q0Diss,Zeta0Diss
  !\
  ! This experiment is aimed to test heat conduction: Test4
  ! The exact solutions is known:: Zeldovich & Rayser 1966, pp.510-519
  !/
  call set_parameters_diss_test
  !
  do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
     xxx = x_BLK(i,j,k,globalBLK)
     yyy = y_BLK(i,j,k,globalBLK)
     !
     State_VGB(rho_,i,j,k,globalBLK)   = cOne
     State_VGB(rhoUx_,i,j,k,globalBLK) = cZero
     State_VGB(rhoUy_,i,j,k,globalBLK) = cZero
     State_VGB(rhoUz_,i,j,k,globalBLK) = cZero
     State_VGB(Bx_,i,j,k,globalBLK)    = cZero
     State_VGB(By_,i,j,k,globalBLK)    = cOne ! cZero for test4_2
     State_VGB(Bz_,i,j,k,globalBLK)    = cZero
     !
     call compute_T_diss_test
     !
     State_VGB(P_,i,j,k,globalBLK)     = TExactDiss
  end do; end do; end do
  !
  if(globalBLK==BLKtest.and.iProc==0)then
     write(*,*) ''
     write(*,*) '>>>>>>>>>>                 <<<<<<<<<<'
     write(*,*) ''
     write(*,*) 'Grav0Diss   = ',Grav0Diss
     write(*,*) 'Beta0Diss   = ',Beta0Diss
     write(*,*) 'Length0Diss = ',Length0Diss
     write(*,*) 'Time0Diss   = ',Time0Diss
     write(*,*) 'Rho0Diss    = ',Rho0Diss
     write(*,*) 'Tem0Diss    = ',Tem0Diss
     write(*,*) ''
     write(*,*) 'ThetaDiss   = ',ThetaDiss
     write(*,*) 'DeltaDiss   = ',DeltaDiss
     write(*,*) 'EpsilonDiss = ',EpsilonDiss
     write(*,*) 'RhoDifDiss  = ',RhoDifDiss
     write(*,*) 'yShiftDiss  = ',yShiftDiss
     write(*,*) 'BZ0Diss     = ',BZ0Diss
     write(*,*) 'scaleHeightDiss = ',scaleHeightDiss
     write(*,*) 'scaleFactorDiss = ',scaleFactorDiss
     write(*,*) 'min(Rho),max(Rho) = ',minval(State_VGB(rho_,:,:,:,:)),maxval(State_VGB(rho_,:,:,:,:))
     write(*,*) 'min(p)  ,max(p)   = ',minval(State_VGB(P_,:,:,:,:)),maxval(State_VGB(P_,:,:,:,:))
     write(*,*) 'min(Px) ,max(Px)  = ',minval(State_VGB(rhoUx_,:,:,:,:)),maxval(State_VGB(rhoUx_,:,:,:,:))
     write(*,*) 'min(Py) ,max(Py)  = ',minval(State_VGB(rhoUy_,:,:,:,:)),maxval(State_VGB(rhoUy_,:,:,:,:))
     write(*,*) 'min(Pz) ,max(Pz)  = ',minval(State_VGB(rhoUz_,:,:,:,:)),maxval(State_VGB(rhoUz_,:,:,:,:))
     write(*,*) 'min(Bx) ,max(Bx)  = ',minval(State_VGB(Bx_,:,:,:,:)),maxval(State_VGB(Bx_,:,:,:,:))
     write(*,*) 'min(By) ,max(By)  = ',minval(State_VGB(By_,:,:,:,:)),maxval(State_VGB(By_,:,:,:,:))
     write(*,*) 'min(Bz) ,max(Bz)  = ',minval(State_VGB(Bz_,:,:,:,:)),maxval(State_VGB(Bz_,:,:,:,:))
     write(*,*) ''
     write(*,*) '>>>>>>>>>>                 <<<<<<<<<<'
     write(*,*) ''
  endif
  
Contains
  
  subroutine compute_T_diss_test
    use ModMpi
    implicit none
    
    real:: tDiss,nDiss,qDiss,aDiss
    real:: ZetaDiss,FuncZetaDiss
    !\
    ! For the test run:: N[x,y]=500, y2=4.0, tend=1.0, tstep=0.05
    !/
    tDiss = t0Diss 
    nDiss = ExponentHeat 
    qDiss = q0Diss 
    aDiss = Kappa0Heat
    !\
    ! Compute the Temperature profile at t = Time0Diss
    !/  
    ZetaDiss = sqrt(xxx**2+yyy**2)/ & 
         (aDiss*tDiss*qDiss**nDiss)**(cOne/(nDiss+cTwo))
    !\
    ! Compute FuncZetaDiss(ZetaDiss(R)) and TExactDiss(ZetaDiss(R))
    !/
    if (ZetaDiss.lt.Zeta0Diss) then
       FuncZetaDiss = (nDiss*Zeta0Diss**cTwo/(cTwo* &
            (nDiss+cTwo)))**(cOne/nDiss)*(cOne-     &
            (ZetaDiss/Zeta0Diss)**cTwo)**(cOne/nDiss)
    else
       FuncZetaDiss = cZero
    end if
    TExactDiss = TExact0Diss+(qDiss**cTwo/(aDiss*tDiss))** &
         (cOne/(nDiss+cTwo))*FuncZetaDiss
  end subroutine compute_T_diss_test
  
  subroutine set_parameters_diss_test
    use ModMpi
    implicit none
    
    !\
    ! Define some parameters
    !/
    t0Diss      = cTwo*cTiny
    TExact0Diss = cTiny
    q0Diss      = cHalf
    if (abs(ExponentHeat-cHalf).lt.cTiny)             &
         Zeta0Diss = 2.4796721935
    if (abs(ExponentHeat-cOne).lt.cTiny)              &
         Zeta0Diss = 1.6509636641
    if (abs(ExponentHeat-(cOne+cHalf)).lt.cTiny)      &
         Zeta0Diss = 1.3133586645
    if (abs(ExponentHeat-cTwo).lt.cTiny)              &
         Zeta0Diss = 1.1283791065
    if (abs(ExponentHeat-(cTwo+cHalf)).lt.cTiny)      &
         Zeta0Diss = 1.0115311146
    if (abs(ExponentHeat-(cTwo+cOne)).lt.cTiny)       &
         Zeta0Diss = 0.9310743809
    if (abs(ExponentHeat-(cTwo+cOne+cHalf)).lt.cTiny) &
         Zeta0Diss = 0.8723382950
    if (abs(ExponentHeat-cFour).lt.cTiny)             &
         Zeta0Diss = 0.8275982141
    if (abs(ExponentHeat-(cFour+cHalf)).lt.cTiny)     &
         Zeta0Diss = 0.7923995256
    if (abs(ExponentHeat-(cFour+cOne)).lt.cTiny)      &
         Zeta0Diss = 0.7639933228
    if(globalBLK==BLKtest.and.iProc==0)then
       write(*,*) ''
       write(*,*) '>>>>>>>>>>ILIAs heat_test_4<<<<<<<<<<'
       write(*,*) ''
       write(*,*) 't0    = ',t0Diss
       write(*,*) 'T0    = ',TExact0Diss
       write(*,*) 'Q0    = ',q0Diss
       write(*,*) 'Zeta0 = ',Zeta0Diss
       write(*,*) ''
       write(*,*) '>>>>>>>>>>                 <<<<<<<<<<'
       write(*,*) ''
    endif
  end subroutine set_parameters_diss_test
  
end subroutine heat_test_4
!***********************************************************************
subroutine resist_test_1
  use ModProcMH,   ONLY:iProc 
  use ModMain,     ONLY:nI,nJ,nK,gcn,globalBLK,nBLK,&
       BLKtest,unusedBLK,UseGravity
  use ModVarIndexes,ONLY:rho_,&
       rhoUx_,rhoUy_,rhoUz_,Bx_, By_, Bz_, P_
  use ModGeometry, ONLY:x_BLK,y_BLK,x1,x2,y1,y2
  use ModAdvance,  ONLY:State_VGB
  use ModPhysics,  ONLY:TypeProblemDiss,EpsilonDiss,DeltaDiss, &
       ThetaDiss,RhoDifDiss,yShiftDiss,BZ0Diss,scaleHeightDiss, &
       scaleFactorDiss,Gbody,Grav0Diss,Beta0Diss,Length0Diss, &
       Time0Diss,Rho0Diss,Tem0Diss,Eta0Resist
  use ModNumConst, ONLY:cOne,cTwo,cFour,cHalf, &
       cZero,cTiny,cHundred,cPi
  use ModMpi
  implicit none
  
  logical:: UseExactFormResist = .true.
  integer:: i,j,k
  real:: xxx,t0Diss,Zeta0Diss
  real:: BExactDiss,BExact0Diss
  !\
  ! This experiment is aimed to test constant resistivity: Test1
  ! The exact solutions is known:: Priest 1982, pp.97-98
  !/
  call set_parameters_diss_test
  !
  do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
     xxx = x_BLK(i,j,k,globalBLK)
     !
     State_VGB(rho_,i,j,k,globalBLK)   = cOne
     State_VGB(rhoUx_,i,j,k,globalBLK) = cZero
     State_VGB(rhoUy_,i,j,k,globalBLK) = cZero
     State_VGB(rhoUz_,i,j,k,globalBLK) = cZero
     State_VGB(P_,i,j,k,globalBLK)     = cOne
     State_VGB(Bx_,i,j,k,globalBLK)    = cZero
     !
     call compute_B_diss_test
     !
     State_VGB(By_,i,j,k,globalBLK)    = BExactDiss
     State_VGB(Bz_,i,j,k,globalBLK)    = cZero
     
  end do; end do; end do
  !
  if(globalBLK==BLKtest.and.iProc==0)then
     write(*,*) ''
     write(*,*) '>>>>>>>>>>                 <<<<<<<<<<'
     write(*,*) ''
     write(*,*) 'Grav0Diss   = ',Grav0Diss
     write(*,*) 'Beta0Diss   = ',Beta0Diss
     write(*,*) 'Length0Diss = ',Length0Diss
     write(*,*) 'Time0Diss   = ',Time0Diss
     write(*,*) 'Rho0Diss    = ',Rho0Diss
     write(*,*) 'Tem0Diss    = ',Tem0Diss
     write(*,*) ''
     write(*,*) 'ThetaDiss   = ',ThetaDiss
     write(*,*) 'DeltaDiss   = ',DeltaDiss
     write(*,*) 'EpsilonDiss = ',EpsilonDiss
     write(*,*) 'RhoDifDiss  = ',RhoDifDiss
     write(*,*) 'yShiftDiss  = ',yShiftDiss
     write(*,*) 'BZ0Diss     = ',BZ0Diss
     write(*,*) 'scaleHeightDiss = ',scaleHeightDiss
     write(*,*) 'scaleFactorDiss = ',scaleFactorDiss
     write(*,*) 'min(Rho),max(Rho) = ',minval(State_VGB(rho_,:,:,:,:)),maxval(State_VGB(rho_,:,:,:,:))
     write(*,*) 'min(p)  ,max(p)   = ',minval(State_VGB(P_,:,:,:,:)),maxval(State_VGB(P_,:,:,:,:))
     write(*,*) 'min(Px) ,max(Px)  = ',minval(State_VGB(rhoUx_,:,:,:,:)),maxval(State_VGB(rhoUx_,:,:,:,:))
     write(*,*) 'min(Py) ,max(Py)  = ',minval(State_VGB(rhoUy_,:,:,:,:)),maxval(State_VGB(rhoUy_,:,:,:,:))
     write(*,*) 'min(Pz) ,max(Pz)  = ',minval(State_VGB(rhoUz_,:,:,:,:)),maxval(State_VGB(rhoUz_,:,:,:,:))
     write(*,*) 'min(Bx) ,max(Bx)  = ',minval(State_VGB(Bx_,:,:,:,:)),maxval(State_VGB(Bx_,:,:,:,:))
     write(*,*) 'min(By) ,max(By)  = ',minval(State_VGB(By_,:,:,:,:)),maxval(State_VGB(By_,:,:,:,:))
     write(*,*) 'min(Bz) ,max(Bz)  = ',minval(State_VGB(Bz_,:,:,:,:)),maxval(State_VGB(Bz_,:,:,:,:))
     write(*,*) ''
     write(*,*) '>>>>>>>>>>                 <<<<<<<<<<'
     write(*,*) ''
  endif
  
Contains
  
  subroutine compute_B_diss_test
    use ModMpi
    implicit none
    
    !\
    ! For the test run:: N[x,y]=128,y2=0.25,tend=0.0015,tstep=0.00005
    !/
    !\
    ! Compute Zeta0Diss = xxx/sqrt(cFour*Eta0Resist*t0Diss)::
    !/  
    Zeta0Diss = xxx/sqrt(cFour*Eta0Resist*t0Diss) 
    !\
    ! Compute BsssResistTest at time t = 0 
    ! (if UseExactFormResist == .true.) 
    !          OR at t = t0Diss 
    ! (if UseExactFormResist == .false.)
    !/
    if (UseExactFormResist) then
       if (xxx.lt.cZero) then
          BExactDiss = -BExact0Diss
       else
          BExactDiss = BExact0Diss
       end if
    else
       if (abs(Zeta0Diss).le.cOne) then
          BExactDiss = BExact0Diss*Zeta0Diss
       else if (Zeta0Diss.gt.cOne) then
          BExactDiss = BExact0Diss
       else
          BExactDiss = -BExact0Diss
       end if
    end if
  end subroutine compute_B_diss_test
  
  subroutine set_parameters_diss_test
    use ModMpi
    implicit none
    
    !\
    ! Define some parameters
    !/
    t0Diss      = cFour*cTiny
    BExact0Diss = cOne
    if(globalBLK==BLKtest.and.iProc==0)then
       write(*,*) ''
       write(*,*) '>>>>>>>>>>ILIAs resist_test_1<<<<<<<<<<'
       write(*,*) ''
       write(*,*) 't0   = ',t0Diss
       write(*,*) 'B0   = ',BExact0Diss
       write(*,*) 'Eta0 = ',Eta0Resist

       write(*,*) ''
       write(*,*) '>>>>>>>>>>                 <<<<<<<<<<'
       write(*,*) ''
    endif
  end subroutine set_parameters_diss_test
  
end subroutine resist_test_1
!***********************************************************************
subroutine mag_rec_1
  use ModProcMH,   ONLY:iProc 
  use ModMain,     ONLY:nI,nJ,nK,gcn,globalBLK,nBLK,&
       BLKtest,unusedBLK,UseGravity
  use ModVarIndexes,ONLY:rho_,&
       rhoUx_,rhoUy_,rhoUz_,Bx_, By_, Bz_, P_
  use ModGeometry, ONLY:x_BLK,y_BLK,x1,x2,y1,y2
  use ModAdvance,  ONLY:State_VGB
  use ModPhysics,  ONLY:TypeProblemDiss,EpsilonDiss,DeltaDiss, &
       ThetaDiss,RhoDifDiss,yShiftDiss,BZ0Diss,scaleHeightDiss, &
       scaleFactorDiss,Gbody,Grav0Diss,Beta0Diss,Length0Diss, &
       Time0Diss,Rho0Diss,Tem0Diss
  use ModNumConst, ONLY:cOne,cTwo,cFour,cHalf, &
       cZero,cTiny,cHundred,cPi
  implicit none

  integer i,j,k, iBLK
  real xxx,yyy
  real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK) :: &
       pGrav1Diss,pGrav2Diss
  !\
  ! This experiment is aimed to test magnetic reconnection: Test1
  !/
  ! Set uniform gravity acceleration
  if (UseGravity) then
     Gbody = -Grav0Diss
  else
     Gbody = cZero
  end if
  !
  do iBLK=1,nBLK
     if (unusedBLK(iBLK)) CYCLE
     do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
        xxx = x_BLK(i,j,k,iBLK)
        yyy = y_BLK(i,j,k,iBLK)+yShiftDiss*y2
        !
        State_VGB(rhoUx_,i,j,k,iBLK)  = cZero
        State_VGB(rhoUy_,i,j,k,iBLK)  = cZero
        State_VGB(rhoUz_,i,j,k,iBLK)  = cZero
        State_VGB(Bx_,i,j,k,iBLK)     = cZero
        State_VGB(By_,i,j,k,iBLK)     = tanh(DeltaDiss*xxx) 
        State_VGB(Bz_,i,j,k,iBLK)     = cZero
        State_VGB(rho_,i,j,k,iBLK)    = (cOne+RhoDifDiss*(cOne-tanh(EpsilonDiss*yyy)))*    &
                                 ((cOne+Beta0Diss-State_VGB(By_,i,j,k,iBLK)**2)/           &
                                 Beta0Diss)**ThetaDiss
        pGrav1Diss(i,j,k,iBLK) = -cTwo*Gbody*((RhoDifDiss+cOne)*yyy-(RhoDifDiss/    &
                                 EpsilonDiss)*log(cosh(EpsilonDiss*yyy)))*((cOne+   &
                                 Beta0Diss-State_VGB(By_,i,j,k,iBLK)**2)/Beta0Diss)**      &
                                 ThetaDiss
        pGrav2Diss(i,j,k,iBLK) = -cTwo*Gbody*((RhoDifDiss+cOne)*(y2+yShiftDiss*y2)- &
                                 (RhoDifDiss/EpsilonDiss)*log(cosh(EpsilonDiss*(y2+ &
                                 yShiftDiss*y2))))*((cOne+Beta0Diss-                &
                                 State_VGB(By_,i,j,k,iBLK)**2)/Beta0Diss)**ThetaDiss
        State_VGB(P_,i,j,k,iBLK)      = (cOne+Beta0Diss-State_VGB(By_,i,j,k,iBLK)**2+             &
                                 pGrav2Diss(i,j,k,iBLK)-pGrav1Diss(i,j,k,iBLK))/cTwo
     end do; end do; end do
  end do
  if(globalBLK==BLKtest.and.iProc==0)then
     write(*,*) ''
     write(*,*) '>>>>>>>>>>ILIAs mag_rec_1<<<<<<<<<<'
     write(*,*) ''
     write(*,*) 'Grav0Diss   = ',Gbody
     write(*,*) 'Beta0Diss   = ',Beta0Diss
     write(*,*) 'Length0Diss = ',Length0Diss
     write(*,*) 'Time0Diss   = ',Time0Diss
     write(*,*) 'Rho0Diss    = ',Rho0Diss
     write(*,*) 'Tem0Diss    = ',Tem0Diss
     write(*,*) ''
     write(*,*) 'ThetaDiss   = ',ThetaDiss
     write(*,*) 'DeltaDiss   = ',DeltaDiss
     write(*,*) 'EpsilonDiss = ',EpsilonDiss
     write(*,*) 'RhoDifDiss  = ',RhoDifDiss
     write(*,*) 'yShiftDiss  = ',yShiftDiss
     write(*,*) 'BZ0Diss     = ',BZ0Diss
     write(*,*) 'scaleHeightDiss = ',scaleHeightDiss
     write(*,*) 'scaleFactorDiss = ',scaleFactorDiss
     write(*,*) 'min(Rho),max(Rho) = ',minval(State_VGB(rho_,:,:,:,:)),maxval(State_VGB(rho_,:,:,:,:))
     write(*,*) 'min(p)  ,max(p)   = ',minval(State_VGB(P_,:,:,:,:)),maxval(State_VGB(P_,:,:,:,:))
     write(*,*) 'min(Px) ,max(Px)  = ',minval(State_VGB(rhoUx_,:,:,:,:)),maxval(State_VGB(rhoUx_,:,:,:,:))
     write(*,*) 'min(Py) ,max(Py)  = ',minval(State_VGB(rhoUy_,:,:,:,:)),maxval(State_VGB(rhoUy_,:,:,:,:))
     write(*,*) 'min(Pz) ,max(Pz)  = ',minval(State_VGB(rhoUz_,:,:,:,:)),maxval(State_VGB(rhoUz_,:,:,:,:))
     write(*,*) 'min(Bx) ,max(Bx)  = ',minval(State_VGB(Bx_,:,:,:,:)),maxval(State_VGB(Bx_,:,:,:,:))
     write(*,*) 'min(By) ,max(By)  = ',minval(State_VGB(By_,:,:,:,:)),maxval(State_VGB(By_,:,:,:,:))
     write(*,*) 'min(Bz) ,max(Bz)  = ',minval(State_VGB(Bz_,:,:,:,:)),maxval(State_VGB(Bz_,:,:,:,:))
     write(*,*) ''
     write(*,*) '>>>>>>>>>>                 <<<<<<<<<<'
     write(*,*) ''
  endif
end subroutine mag_rec_1
!***********************************************************************
subroutine mag_rec_2
  use ModProcMH,   ONLY:iProc 
  use ModMain,     ONLY:nI,nJ,nK,gcn,globalBLK,nBLK,&
       BLKtest,unusedBLK,UseGravity
  use ModVarIndexes,ONLY:rho_,&
       rhoUx_,rhoUy_,rhoUz_,Bx_, By_, Bz_, P_
  use ModGeometry, ONLY:x_BLK,y_BLK,x1,x2,y1,y2
  use ModAdvance,  ONLY:State_VGB
  use ModPhysics,  ONLY:TypeProblemDiss,EpsilonDiss,DeltaDiss, &
       ThetaDiss,RhoDifDiss,yShiftDiss,BZ0Diss,scaleHeightDiss, &
       scaleFactorDiss,Gbody,Grav0Diss,Beta0Diss,Length0Diss, &
       Time0Diss,Rho0Diss,Tem0Diss
  use ModNumConst, ONLY:cOne,cTwo,cFour,cHalf, &
       cZero,cTiny,cHundred,cPi
  implicit none

  integer i,j,k,l, iBLK
  real xxx,yyy
  real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK) :: &
       ProfileDiss
  real, dimension(7) :: PolyIndexDiss
  !\
  ! This experiment is aimed to test magnetic reconnection: Test2
  !/
  ! Set uniform gravity acceleration
  if (UseGravity) then
     Gbody = -Grav0Diss
  else
     Gbody = cZero
  end if
  !
  PolyIndexDiss(1) =  cOne 
  PolyIndexDiss(2) =  2.9540557E-05 
  PolyIndexDiss(3) = -8.2382945E-11
  PolyIndexDiss(4) =  1.5141202E-16
  PolyIndexDiss(5) = -1.5301333E-22
  PolyIndexDiss(6) =  7.8460268E-29
  PolyIndexDiss(7) = -1.5958496E-35

  do iBLK=1,nBLK
     if (unusedBLK(iBLK)) CYCLE
     do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
        xxx = x_BLK(i,j,k,iBLK)
        yyy = y_BLK(i,j,k,iBLK)
        do l=1,7
           ProfileDiss(i,j,k,iBLK) = PolyIndexDiss(l)* &
                (1.4402537E+6*yyy)**(l-1)
        end do
        State_VGB(rhoUx_,i,j,k,iBLK) = cZero
        State_VGB(rhoUy_,i,j,k,iBLK) = cZero
        State_VGB(rhoUz_,i,j,k,iBLK) = cZero
        State_VGB(Bx_,i,j,k,iBLK)    = cZero
        State_VGB(By_,i,j,k,iBLK)    = tanh(DeltaDiss*xxx) 
        State_VGB(Bz_,i,j,k,iBLK)    = cZero
        State_VGB(rho_,i,j,k,iBLK)   = (cOne/ProfileDiss(i,j,k,iBLK))*          &
                                ((cOne+Beta0Diss-State_VGB(By_,i,j,k,iBLK)**2)/ &
                                Beta0Diss)**ThetaDiss
        State_VGB(P_,i,j,k,iBLK)     = (cOne+Beta0Diss-State_VGB(By_,i,j,k,iBLK)**2)/cTwo
     end do; end do; end do
  end do
  if(globalBLK==BLKtest.and.iProc==0)then
     write(*,*) ''
     write(*,*) '>>>>>>>>>>ILIAs mag_rec_2<<<<<<<<<<'
     write(*,*) ''
     write(*,*) 'Grav0Diss   = ',Gbody
     write(*,*) 'Beta0Diss   = ',Beta0Diss
     write(*,*) 'Length0Diss = ',Length0Diss
     write(*,*) 'Time0Diss   = ',Time0Diss
     write(*,*) 'Rho0Diss    = ',Rho0Diss
     write(*,*) 'Tem0Diss    = ',Tem0Diss
     write(*,*) ''
     write(*,*) 'ThetaDiss   = ',ThetaDiss
     write(*,*) 'DeltaDiss   = ',DeltaDiss
     write(*,*) 'EpsilonDiss = ',EpsilonDiss
     write(*,*) 'RhoDifDiss  = ',RhoDifDiss
     write(*,*) 'yShiftDiss  = ',yShiftDiss
     write(*,*) 'BZ0Diss     = ',BZ0Diss
     write(*,*) 'scaleHeightDiss = ',scaleHeightDiss
     write(*,*) 'scaleFactorDiss = ',scaleFactorDiss
     write(*,*) 'min(Rho),max(Rho) = ',minval(State_VGB(rho_,:,:,:,:)),&
          maxval(State_VGB(rho_,:,:,:,:))
     write(*,*) 'min(p)  ,max(p)   = ',minval(State_VGB(P_,:,:,:,:)),&
          maxval(State_VGB(P_,:,:,:,:))
     write(*,*) 'min(Px) ,max(Px)  = ',minval(State_VGB(rhoUx_,:,:,:,:)),&
          maxval(State_VGB(rhoUx_,:,:,:,:))
     write(*,*) 'min(Py) ,max(Py)  = ',minval(State_VGB(rhoUy_,:,:,:,:)),&
          maxval(State_VGB(rhoUy_,:,:,:,:))
     write(*,*) 'min(Pz) ,max(Pz)  = ',minval(State_VGB(rhoUz_,:,:,:,:)),&
          maxval(State_VGB(rhoUz_,:,:,:,:))
     write(*,*) 'min(Bx) ,max(Bx)  = ',minval(State_VGB(Bx_,:,:,:,:)),&
          maxval(State_VGB(Bx_,:,:,:,:))
     write(*,*) 'min(By) ,max(By)  = ',minval(State_VGB(By_,:,:,:,:)),&
          maxval(State_VGB(By_,:,:,:,:))
     write(*,*) 'min(Bz) ,max(Bz)  = ',minval(State_VGB(Bz_,:,:,:,:)),&
          maxval(State_VGB(Bz_,:,:,:,:))
     write(*,*) ''
     write(*,*) '>>>>>>>>>>                 <<<<<<<<<<'
     write(*,*) ''
  endif
end subroutine mag_rec_2
!***********************************************************************
subroutine mag_rec_3
  use ModProcMH,   ONLY:iProc 
  use ModMain,     ONLY:nI,nJ,nK,gcn,globalBLK,nBLK,&
       BLKtest,unusedBLK,UseGravity
  use ModVarIndexes,ONLY:rho_,&
       rhoUx_,rhoUy_,rhoUz_,Bx_, By_, Bz_, P_
  use ModGeometry, ONLY:x_BLK,y_BLK,x1,x2,y1,y2
  use ModAdvance,  ONLY:State_VGB
  use ModPhysics,  ONLY:TypeProblemDiss,EpsilonDiss,DeltaDiss, &
       ThetaDiss,RhoDifDiss,yShiftDiss,BZ0Diss,scaleHeightDiss, &
       scaleFactorDiss,Gbody,Grav0Diss,Beta0Diss,Length0Diss, &
       Time0Diss,Rho0Diss,Tem0Diss
   use ModNumConst, ONLY:cOne,cTwo,cFour,cHalf, &
       cZero,cTiny,cHundred,cPi
  use ModMpi
  implicit none
  
  integer:: i,j,k,iBLK
  real:: xxx,yyy
  real:: bDiss,cDiss,mDiss
  real:: BX0RefDiss,BY0RefDiss
  real:: Rho0RefDiss,p0RefDiss,p1RefDiss
  real:: scrARR1,scrARR2,scrARR3,scrARR4
  !\
  ! This experiment is aimed to test magnetic reconnection: Test3
  !/
  !\
  ! Calculate the reference values of some physical cuantities
  !/
  call calculate_ref_values_diss_test
  !
  do iBLK=1,nBLK
     if (unusedBLK(iBLK)) CYCLE
     do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
        xxx  = x_BLK(i,j,k,iBLK)
        yyy  = y_BLK(i,j,k,iBLK)
        scrARR1 = log(cosh(bDiss*xxx))
        scrARR2 = exp(-yyy/cDiss)
        scrARR3 = (cOne/bDiss/cDiss)*scrARR1*scrARR2
        scrARR4 = tanh(bDiss*xxx)*scrARR2
        !
        State_VGB(rhoUx_,i,j,k,iBLK) = cZero
        State_VGB(rhoUy_,i,j,k,iBLK) = cZero
        State_VGB(rhoUz_,i,j,k,iBLK) = cZero
        State_VGB(Bx_,i,j,k,iBLK)    = mDiss*scrARR3
        State_VGB(By_,i,j,k,iBLK)    = mDiss*scrARR4
        State_VGB(Bz_,i,j,k,iBLK)    = BZ0Diss
        State_VGB(rho_,i,j,k,iBLK)   = (Rho0RefDiss+(mDiss**2/cDiss/Grav0Diss)*        &
                                (cOne+scrARR1)/cosh(bDiss*xxx)**2)*             &
                                scrARR2**2
        State_VGB(P_,i,j,k,iBLK)     = -(p0RefDiss+cHalf*(mDiss/bDiss/cDiss)**2*       &
                                (scrARR1**2-(bDiss*cDiss/cosh(bDiss*xxx))**2))* &
                                scrARR2**2+p1RefDiss+Beta0Diss/cTwo
     end do; end do; end do
  end do
  if(globalBLK==BLKtest.and.iProc==0)then
     write(*,*) ''
     write(*,*) '>>>>>>>>>>ILIAs mag_rec_3<<<<<<<<<<'
     write(*,*) ''
     write(*,*) 'Grav0Diss   = ',Gbody
     write(*,*) 'Beta0Diss   = ',Beta0Diss
     write(*,*) 'Length0Diss = ',Length0Diss
     write(*,*) 'Time0Diss   = ',Time0Diss
     write(*,*) 'Rho0Diss    = ',Rho0Diss
     write(*,*) 'Tem0Diss    = ',Tem0Diss
     write(*,*) ''
     write(*,*) 'ThetaDiss   = ',ThetaDiss
     write(*,*) 'DeltaDiss   = ',bDiss
     write(*,*) 'EpsilonDiss = ',EpsilonDiss
     write(*,*) 'RhoDifDiss  = ',RhoDifDiss
     write(*,*) 'yShiftDiss  = ',yShiftDiss
     write(*,*) 'BZ0Diss     = ',BZ0Diss
     write(*,*) '1/|B|       = ',mDiss
     write(*,*) 'scaleFactorDiss   = ',scaleFactorDiss
     write(*,*) '1/scaleHeightDiss = ',cDiss
     write(*,*) 'min(Rho),max(Rho) = ',minval(State_VGB(rho_,:,:,:,:)),&
          maxval(State_VGB(rho_,:,:,:,:))
     write(*,*) 'min(p)  ,max(p)   = ',minval(State_VGB(P_,:,:,:,:)),&
          maxval(State_VGB(P_,:,:,:,:))
     write(*,*) 'min(Px) ,max(Px)  = ',minval(State_VGB(rhoUx_,:,:,:,:)),&
          maxval(State_VGB(rhoUx_,:,:,:,:))
     write(*,*) 'min(Py) ,max(Py)  = ',minval(State_VGB(rhoUy_,:,:,:,:)),&
          maxval(State_VGB(rhoUy_,:,:,:,:))
     write(*,*) 'min(Pz) ,max(Pz)  = ',minval(State_VGB(rhoUz_,:,:,:,:)),&
          maxval(State_VGB(rhoUz_,:,:,:,:))
     write(*,*) 'min(Bx) ,max(Bx)  = ',minval(State_VGB(Bx_,:,:,:,:)),&
          maxval(State_VGB(Bx_,:,:,:,:))
     write(*,*) 'min(By) ,max(By)  = ',minval(State_VGB(By_,:,:,:,:)),&
          maxval(State_VGB(By_,:,:,:,:))
     write(*,*) 'min(Bz) ,max(Bz)  = ',minval(State_VGB(Bz_,:,:,:,:)),&
          maxval(State_VGB(Bz_,:,:,:,:))
     write(*,*) ''
     write(*,*) '>>>>>>>>>>                 <<<<<<<<<<'
     write(*,*) ''
  endif

Contains
  
  subroutine calculate_ref_values_diss_test
    use ModMpi
    implicit none
    
    !\
    ! Set uniform gravity acceleration:: 
    !/
    if (UseGravity) then
       Gbody = -Grav0Diss
    else
       Gbody = cZero
    end if
    !\
    ! Redefine some parameters::
    ! Follow the notations in Roussev etal. 2002b
    !/
    cDiss = scaleFactorDiss*y2/scaleHeightDiss
    bDiss = DeltaDiss
    !\
    ! Compute the magnitude of the `B' field at the 
    ! reference location -- (x1,y1)::
    !/
    BX0RefDiss = (cOne/bDiss/cDiss)*log(cosh(bDiss*x1))* &
                 exp(-y1/cDiss)
    BY0RefDiss = tanh(bDiss*x1)*exp(-y1/cDiss)
    mDiss      = cOne/sqrt(BX0RefDiss**2+BY0RefDiss**2+  &
                 BZ0Diss**2)
    !\
    ! Compute the reference `rho' and `p' at (x1,y1)::
    !/
    Rho0RefDiss = exp(cTwo*y1/cDiss)
    p0RefDiss   = cHalf*Rho0RefDiss*Grav0Diss*cDiss
    p1RefDiss   = (p0RefDiss+cHalf*(mDiss/bDiss/cDiss)**2* &
                  (log(cosh(bDiss*x1))**2-(bDiss*cDiss/    &
                  cosh(bDiss*x1))**2))*exp(-cTwo*y1/cDiss)
  end subroutine calculate_ref_values_diss_test
  
end subroutine mag_rec_3
!***********************************************************************
subroutine loop_rec_1
  use ModProcMH,   ONLY:iProc 
  use ModMain,     ONLY:nI,nJ,nK,gcn,globalBLK,nBLK,&
       BLKtest,unusedBLK,UseGravity
  use ModVarIndexes,ONLY:rho_,&
       rhoUx_,rhoUy_,rhoUz_,Bx_, By_, Bz_, P_
  use ModGeometry, ONLY:x_BLK,y_BLK,x1,x2,y1,y2
  use ModAdvance,  ONLY:State_VGB
  use ModPhysics,  ONLY:TypeProblemDiss,EpsilonDiss,DeltaDiss, &
       ThetaDiss,RhoDifDiss,yShiftDiss,BZ0Diss,scaleHeightDiss, &
       scaleFactorDiss,Gbody,Grav0Diss,Beta0Diss,Length0Diss, &
       Time0Diss,Rho0Diss,Tem0Diss
  use ModNumConst, ONLY:cOne,cTwo,cFour,cHalf, &
       cZero,cTiny,cHundred,cPi
  implicit none

  integer i,j,k, iBLK
  real xxx,yyy
  real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK) :: &
       pGravDiss
  !\
  ! This experiment is aimed to test magnetic reconnection: Test1
  !/
  ! Set uniform gravity acceleration
  if (UseGravity) then
     Gbody = -Grav0Diss
  else
     Gbody = cZero
  end if
  !
  do iBLK=1,nBLK
     if (unusedBLK(iBLK)) CYCLE
     do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
        xxx = x_BLK(i,j,k,iBLK)
        yyy = y_BLK(i,j,k,iBLK)
        !
        State_VGB(rhoUx_,i,j,k,iBLK) = cZero
        State_VGB(rhoUy_,i,j,k,iBLK) = cZero
        State_VGB(rhoUz_,i,j,k,iBLK) = cZero
        State_VGB(Bx_,i,j,k,iBLK)    = -(cTwo*scaleHeightDiss/cPi)*cos(cPi*xxx/cTwo)*      &
                                exp(-scaleHeightDiss*yyy)/1.53718
        State_VGB(By_,i,j,k,iBLK)    = sin(cPi*xxx/cTwo)*exp(-scaleHeightDiss*yyy)/1.53718
        State_VGB(Bz_,i,j,k,iBLK)    = cZero
        State_VGB(rho_,i,j,k,iBLK)   = cOne+RhoDifDiss*(cOne-tanh(EpsilonDiss*(yyy+        &
                                yShiftDiss*y2)))
        pGravDiss(i,j,k,iBLK) = -cTwo*Gbody*((RhoDifDiss+cOne)*(yyy+yShiftDiss*y2)- &
                                (RhoDifDiss/EpsilonDiss)*log(cosh(EpsilonDiss*(yyy+ &
                                yShiftDiss*y2))))
        State_VGB(P_,i,j,k,iBLK)     = (Beta0Diss*exp(cTwo*scaleHeightDiss*(yyy+           &
                                yShiftDiss*y2))-pGravDiss(i,j,k,iBLK))/cTwo
     end do; end do; end do
  end do
end subroutine loop_rec_1
!***********************************************************************
subroutine spicules_1
  use ModProcMH,   ONLY:iProc 
  use ModMain,     ONLY:nI,nJ,nK,gcn,globalBLK,nBLK,&
       BLKtest,unusedBLK,UseGravity
  use ModVarIndexes,ONLY:rho_,&
       rhoUx_,rhoUy_,rhoUz_,Bx_, By_, Bz_, P_
  use ModGeometry, ONLY:x_BLK,y_BLK,x1,x2,y1,y2
  use ModAdvance,  ONLY:State_VGB
  use ModPhysics,  ONLY:TypeProblemDiss,EpsilonDiss,DeltaDiss, &
       ThetaDiss,RhoDifDiss,yShiftDiss,BZ0Diss,scaleHeightDiss, &
       scaleFactorDiss,Gbody,Grav0Diss,Beta0Diss,Length0Diss, &
       Time0Diss,Rho0Diss,Tem0Diss
  use ModNumConst, ONLY:cOne,cTwo,cFour,cHalf, &
       cZero,cTiny,cHundred,cPi
  implicit none
  
  integer i,j,k, iBLK
  real xxx,yyy
  real B2InvDiss,B2MinDiss,pGrav0Diss
  real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK) :: &
       pGravDiss,scrARR1,scrARR2
  !\
  ! This experiment is aimed to test magnetic reconnection: Test1
  !/
  ! Set uniform gravity acceleration
  if (UseGravity) then
     Gbody = -Grav0Diss
  else
     Gbody = cZero
  end if
  !
  do iBLK=1,nBLK
     if (unusedBLK(iBLK)) CYCLE
     do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
        xxx = x_BLK(i,j,k,iBLK)
        yyy = y_BLK(i,j,k,iBLK)
        !
        scrARR1(i,j,k,iBLK)   = sin(cPi-cPi*scaleHeightDiss*xxx/x2)*           &
                                exp(-cPi*scaleHeightDiss*yyy/y2)
        scrARR2(i,j,k,iBLK)   = -cos(cPi-cPi*scaleHeightDiss*xxx/x2)*          &
                                exp(-cPi*scaleHeightDiss*yyy/y2)
        pGravDiss(i,j,k,iBLK) = -Gbody*((RhoDifDiss+cOne)*(yyy+yShiftDiss*y2)- &
                                (RhoDifDiss/EpsilonDiss)*log(cosh(EpsilonDiss* &
                                (yyy+yShiftDiss*y2))))/(cOne+cTwo*RhoDifDiss)
     end do; end do; end do
  end do
  B2InvDiss  = cOne/sqrt(maxval(scrARR1**2+scrARR2**2))
  B2MinDiss  = minval(scrARR1**2+scrARR2**2)
  pGrav0Diss = -Gbody*((RhoDifDiss+cOne)*(y2+yShiftDiss*y2)- &
       (RhoDifDiss/EpsilonDiss)*log(cosh(EpsilonDiss*(y2+    &
       yShiftDiss*y2))))/(cOne+cTwo*RhoDifDiss)
  do iBLK=1,nBLK
     if (unusedBLK(iBLK)) CYCLE
     do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
        xxx = x_BLK(i,j,k,iBLK)
        yyy = y_BLK(i,j,k,iBLK)+yShiftDiss*y2
        !
        State_VGB(rhoUx_,i,j,k,iBLK) = cZero
        State_VGB(rhoUy_,i,j,k,iBLK) = cZero
        State_VGB(rhoUz_,i,j,k,iBLK) = cZero
        State_VGB(Bx_,i,j,k,iBLK)    = B2InvDiss*scrARR1(i,j,k,iBLK)
        State_VGB(By_,i,j,k,iBLK)    = B2InvDiss*scrARR2(i,j,k,iBLK)
        State_VGB(Bz_,i,j,k,iBLK)    = cZero
        State_VGB(rho_,i,j,k,iBLK)   = (cOne+RhoDifDiss*(cOne-tanh(EpsilonDiss*yyy)))/ &
                                (cOne+cTwo*RhoDifDiss)
        State_VGB(P_,i,j,k,iBLK)     = (Beta0Diss*B2MinDiss/cTwo+pGrav0Diss-           &
                                pGravDiss(i,j,k,iBLK))
     end do; end do; end do
  end do
end subroutine spicules_1
!***********************************************************************
subroutine spicules_2
  use ModProcMH,   ONLY:iProc 
  use ModMain,     ONLY:nI,nJ,nK,gcn,globalBLK,nBLK,&
       BLKtest,unusedBLK,UseGravity
  use ModVarIndexes,ONLY:rho_,&
       rhoUx_,rhoUy_,rhoUz_,Bx_, By_, Bz_, P_
  use ModGeometry, ONLY:x_BLK,y_BLK,x1,x2,y1,y2
  use ModAdvance,  ONLY:State_VGB
  use ModPhysics,  ONLY:TypeProblemDiss,EpsilonDiss,DeltaDiss, &
       ThetaDiss,RhoDifDiss,yShiftDiss,BZ0Diss,scaleHeightDiss, &
       scaleFactorDiss,Gbody,Grav0Diss,Beta0Diss,Length0Diss, &
       Time0Diss,Rho0Diss,Tem0Diss
  use ModNumConst, ONLY:cOne,cTwo,cFour,cHalf, &
       cZero,cTiny,cHundred,cPi
  implicit none

  integer i,j,k, iBLK
  real xxx,yyy
  real B2InvDiss,pGrav0Diss
  real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK) :: &
       pGravDiss,scrARR1,scrARR2
  !\
  ! This experiment is aimed to test magnetic reconnection: Test1
  !/
  ! Set uniform gravity acceleration
  if (UseGravity) then
     Gbody = -Grav0Diss
  else
     Gbody = cZero
  end if
  !
  do iBLK=1,nBLK
     if (unusedBLK(iBLK)) CYCLE
     do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
        xxx = x_BLK(i,j,k,iBLK)
        yyy = y_BLK(i,j,k,iBLK)
        !
        scrARR1(i,j,k,iBLK)   = sinh(xxx)*cos(yyy)
        scrARR2(i,j,k,iBLK)   = -cosh(xxx)*sin(yyy)+cOne+cHalf
        pGravDiss(i,j,k,iBLK) = -cTwo*Gbody*((RhoDifDiss+cOne)* &
                             (yyy+yShiftDiss*y2)-(RhoDifDiss/   &
                             EpsilonDiss)*log(cosh(EpsilonDiss* &
                             (yyy+yShiftDiss*y2))))
     end do; end do; end do
  end do
  B2InvDiss  = cOne/sqrt(maxval(scrARR1**2+scrARR2**2))
  pGrav0Diss = -Gbody*((RhoDifDiss+cOne)*(y2+yShiftDiss*y2)- &
       (RhoDifDiss/EpsilonDiss)*log(cosh(EpsilonDiss*(y2+    &
       yShiftDiss*y2))))/(cOne+cTwo*RhoDifDiss)
  do iBLK=1,nBLK
     if (unusedBLK(iBLK)) CYCLE
     do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
        xxx = x_BLK(i,j,k,iBLK)
        yyy = y_BLK(i,j,k,iBLK)+yShiftDiss*y2
        !        
        State_VGB(rhoUx_,i,j,k,iBLK) = cZero
        State_VGB(rhoUy_,i,j,k,iBLK) = cZero
        State_VGB(rhoUz_,i,j,k,iBLK) = cZero
        State_VGB(Bx_,i,j,k,iBLK)    = B2InvDiss*scrARR1(i,j,k,iBLK)
        State_VGB(By_,i,j,k,iBLK)    = B2InvDiss*scrARR2(i,j,k,iBLK)
        State_VGB(Bz_,i,j,k,iBLK)    = cZero
        State_VGB(rho_,i,j,k,iBLK)   = cOne+RhoDifDiss*(cOne-tanh(EpsilonDiss*yyy))
        State_VGB(P_,i,j,k,iBLK)     = (Beta0Diss+pGrav0Diss-pGravDiss(i,j,k,iBLK))/cTwo
     end do; end do; end do
  end do
end subroutine spicules_2
!***********************************************************************
subroutine spicules_3
  use ModProcMH,   ONLY:iProc 
  use ModMain,     ONLY:nI,nJ,nK,gcn,globalBLK,nBLK,&
       BLKtest,unusedBLK,UseGravity
  use ModVarIndexes,ONLY:rho_,&
       rhoUx_,rhoUy_,rhoUz_,Bx_, By_, Bz_, P_
  use ModGeometry, ONLY:x_BLK,y_BLK,x1,x2,y1,y2
  use ModAdvance,  ONLY:State_VGB
  use ModPhysics,  ONLY:TypeProblemDiss,EpsilonDiss,DeltaDiss, &
       ThetaDiss,RhoDifDiss,yShiftDiss,BZ0Diss,scaleHeightDiss, &
       scaleFactorDiss,Gbody,Grav0Diss,Beta0Diss,Length0Diss, &
       Time0Diss,Rho0Diss,Tem0Diss
  use ModNumConst, ONLY:cOne,cTwo,cFour,cHalf, &
       cZero,cTiny,cHundred,cPi
  implicit none
  
  integer i,j,k, iBLK
  real xxx,yyy,pGrav0Diss
  real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK) :: &
       pGravDiss,scrARR1,scrARR2
  !\
  ! This experiment is aimed to test magnetic reconnection: Test1
  !/
  ! Set uniform gravity acceleration
  if (UseGravity) then
     Gbody = -Grav0Diss
  else
     Gbody = cZero
  end if
  !
  do iBLK=1,nBLK
     if (unusedBLK(iBLK)) CYCLE
     do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
        xxx = x_BLK(i,j,k,iBLK)
        yyy = y_BLK(i,j,k,iBLK)+yShiftDiss*y2
        !
        pGravDiss(i,j,k,iBLK) = -cTwo*Gbody*((RhoDifDiss+cOne)* &
             yyy-(RhoDifDiss/EpsilonDiss)*log(cosh(EpsilonDiss*yyy)))
     end do; end do; end do
  end do
  pGrav0Diss = -Gbody*((RhoDifDiss+cOne)*(y2+yShiftDiss*y2)- &
       (RhoDifDiss/EpsilonDiss)*log(cosh(EpsilonDiss*(y2+    &
       yShiftDiss*y2))))/(cOne+cTwo*RhoDifDiss)
  do iBLK=1,nBLK
     if (unusedBLK(iBLK)) CYCLE
     do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
        xxx = x_BLK(i,j,k,iBLK)
        yyy = y_BLK(i,j,k,iBLK)+yShiftDiss*y2
        !
        State_VGB(rhoUx_,i,j,k,iBLK) = cZero
        State_VGB(rhoUy_,i,j,k,iBLK) = cZero
        State_VGB(rhoUz_,i,j,k,iBLK) = cZero
        State_VGB(Bx_,i,j,k,iBLK)    = cZero
        State_VGB(By_,i,j,k,iBLK)    = cOne
        State_VGB(Bz_,i,j,k,iBLK)    = cZero
        State_VGB(rho_,i,j,k,iBLK)   = cOne+RhoDifDiss*(cOne-tanh(EpsilonDiss*yyy))
        State_VGB(P_,i,j,k,iBLK)     = (Beta0Diss+pGrav0Diss-pGravDiss(i,j,k,iBLK))/cTwo
     end do; end do; end do
  end do
end subroutine spicules_3
!***********************************************************************
