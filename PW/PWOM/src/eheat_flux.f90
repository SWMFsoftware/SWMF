!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
SUBROUTINE PW_eheat_flux
  use ModCommonVariables
  use ModPWOM,ONLY: UseEleHeat
  REAL C1(MaxGrid),C2(MaxGrid),C3(MaxGrid),C4(MaxGrid),&
       D1(MaxGrid),D2(MaxGrid),&
       D3(MaxGrid),D4(MaxGrid)
  REAL YL(MaxGrid),YK(MaxGrid)

  !----------------------------------------------------------------------------
  if(.not.UseEleHeat) return
  C2(1)=H0*(HeatCon_GI(2,nIon)-HeatCon_GI(0,nIon))
  C2(NDIM)=H0*(HeatCon_GI(nDim+1,nIon)-HeatCon_GI(NDIM2,nIon))
  
  !C2 = d/dr Kappa 
  DO K=2,NDIM2
     C2(K)=H0*(HeatCon_GI(K+1,nIon)-HeatCon_GI(K-1,nIon))
  enddo
  ! C1 = Kappa/Rho
  ! C2 = 1/Rho * d/dr Kappa 
  DO K=1,NDIM
     C1(K)=HeatCon_GI(K,nIon)/State_GV(K,RhoE_)
     C2(K)=C2(K)/State_GV(K,RhoE_)
  enddo
  C3(1)=H0*(State_GV(2,uE_)-State_GV(0,uE_))
  C3(NDIM)=H0*(State_GV(nDim+1,uE_)-State_GV(nDim2,uE_))
  ! C3 = du/dr 
  DO  K=2,NDIM2
     C3(K)=H0*(State_GV(K+1,uE_)-State_GV(K-1,uE_))
  enddo
  ! C2 = 1/Rho * d/dr Kappa - u + A'/A * Kappa/Rho
  ! C3 = -(gamma-1)*( A'/A*u + du/dr )
  DO K=1,NDIM
     C2(K)=C2(K)-State_GV(K,uE_)
     C2(K)=C2(K)+DAREA(K)*C1(K)
     C3(K)=-GMIN1*(DAREA(K)*State_GV(K,uE_)+C3(K))
  enddo
  
  XHLP=EXP(-(TIME-300.)**2/2./150./150.)
  ! C3 = -(gamma-1)*( A'/A*u + du/dr ) - S/Rho
  ! C4 = (gamma-1)/R * E/Rho
  DO K=1,NDIM
     C3(K)=C3(K)-Source_CV(K,RhoE_)/State_GV(K,RhoE_)
     C4(K)=HLPE0*(Source_CV(K,pE_)+XHLP*QELECT(K))/State_GV(K,RhoE_)
  enddo
  ! XX1 = 0.5/dr**2 * Kappa/Rho
  ! XX2 = 0.25/dr * (1/Rho * d/dr Kappa - u + A'/A * Kappa/Rho)
  ! XX3 = 1/dr**2 * (Kappa/Rho )
  ! XX4 = 0.5 *( -(gamma-1)*( A'/A*u + du/dr ) - S/Rho ) - 1/dr**2 * (Kappa/Rho ) 
  ! D1 = -0.5/dr**2 * Kappa/Rho  - 0.25/dr * (1/Rho * d/dr Kappa - u + A'/A * Kappa/Rho) 
  ! D2 = 1/dt - 0.5 *( -(gamma-1)*( A'/A*u + du/dr ) - S/Rho ) + 1/dr**2 * (Kappa/Rho )
  ! D3 = -0.5/dr**2 * Kappa/Rho  + 0.25/dr * (1/Rho * d/dr Kappa - u + A'/A * Kappa/Rho)
  ! D4 = 1/dt + 0.5 *( -(gamma-1)*( A'/A*u + du/dr ) - S/Rho ) - 1/dr**2 * (Kappa/Rho ) 
  DO K=1,NDIM
     XX1=H3*C1(K)
     XX2=H4*C2(K)
     XX3=H2*C1(K)
     XX4=0.5*C3(K)-XX3
     D1(K)=-XX1-XX2
     D2(K)=H1E2-XX4
     D3(K)=-XX1+XX2
     D4(K)=H1E2+XX4
  enddo
  
  D4(1)=&
       C4(1)-D3(1)*State_GV(0,Te_)+D4(1)*State_GV(1,Te_)&
       -D1(1)*State_GV(2,Te_)
  D4(NDIM)=C4(NDIM)-D3(NDIM)*State_GV(nDim2,Te_)+ &
       D4(NDIM)*State_GV(nDim,Te_)-D1(NDIM)*State_GV(nDim+1,Te_)
  !D4 = (gamma-1)/R * E/Rho 
  !     -[-0.5/dr**2 * Kappa/Rho  + 0.25/dr * (1/Rho * d/dr Kappa - u + A'/A * Kappa/Rho)]T(k-1)
  !     +[1/dt + 0.5 *( -(gamma-1)*( A'/A*u + du/dr ) - S/Rho ) - 1/dr**2 * (Kappa/Rho ) ]T(k)
  !     -[-0.5/dr**2 * Kappa/Rho  - 0.25/dr * (1/Rho * d/dr Kappa - u + A'/A * Kappa/Rho)]T(k+1)
  !   = (gamma-1)/R * E/Rho + .5*Kappa/rho * d^2T/dr^2 - .5(gamma-1)(A'/A*u + du/dr)T 
  !     - .5*S/rho * T + .5*(1/rho * dk/dr - u + A'/A*Kappa/rho) * dT/dr + T(k)/dt
  !   = (gamma-1)/R*E/Rho - .5*S/Rho *T + T(k)/dt + 1/(2*Rho*A) * d/dr(A*Kappa*dT/dr) 
  !     - (gamma-1)/2A * T * d(Au)/dr - u/2 * dT/dr
  DO K=2,NDIM2
     D4(K)=C4(K)-D3(K)*State_GV(K-1,Te_)+D4(K)*State_GV(K,Te_)-&
          D1(K)*State_GV(K+1,Te_)   
  enddo
  ! We have a tridiagonal system: 
  ! D3(k)T(k-1) + D2(k)T(k) + D1(k)T(k+1) = D4(k)
  ! Use a tridiagonal matrix algoritm to get T(k)
  ! Use a forward sweep to get D3' = 0 and D2' = 1 (' is updated coef)
  ! D3'= YL and D4'=YK. 
  YL(1)=D1(1)/D2(1)
  YK(1)=(D4(1)-D3(1)*State_GV(0,Te_))/D2(1)
  
  DO  K=2,NDIM
     XX1=D2(K)-D3(K)*YL(K-1)
     YK(K)=(D4(K)-D3(K)*YK(K-1))/XX1
     
     YL(K)=D1(K)/XX1
  enddo

  ! Use a back subsitution to get T
  State_GV(nDim,Te_)=YK(NDIM)-YL(NDIM)*State_GV(nDim+1,Te_)
  State_GV(nDim,pE_)=RGAS_I(nIon)*State_GV(nDim,RhoE_)*State_GV(nDim,Te_)
  
  DO K=NDIM2,1,-1
     
     State_GV(K,Te_)=YK(K)-YL(K)*State_GV(K+1,Te_)
     State_GV(K,pE_)=RGAS_I(nIon)*State_GV(K,RhoE_)*State_GV(K,Te_)
     
  enddo
  
  RETURN
end SUBROUTINE PW_eheat_flux
