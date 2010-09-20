
SUBROUTINE PW_calc_efield(nCell,State_GV)
  use ModCommonVariables,only: nDim,nIon,iRho_I,iU_I,iP_I,MassElecIon_I,&
       Source_CV,uE_,pE_,DrBnd,Mass_I,Efield,GRAVTY,&
       dArea,CurrMx,MaxGrid,Te_,nVar, NamePlanet, Ion1_,AltD
  use ModPWOM,    ONLY:UseCentrifugal
  use ModPwWaves, ONLY: UseWaveAcceleration, WaveAcceleration_C, &
                        calc_wave_acceleration
  implicit none
  
  integer,intent(in) :: nCell
  real,   intent(in) :: State_GV(-1:nCell+2,nVar)
  
  integer :: iIon, K
  real    :: PR_GI(0:MaxGrid,nIon),&
             PR(0:MaxGrid),EZ(0:MaxGrid),eNumDens,Dz,Ez1,Dnom
  real, allocatable :: aCentrifugal_C(:)
  
  !----------------------------------------------------------------------------
  
  if (.not.allocated(aCentrifugal_C)) allocate (aCentrifugal_C(nDim))

  PR(:)=0.0
  EZ(:)=0.0
  
  do iIon=1,nIon
     do K=0,nDim+1
        !PR_species = rho * u^2
        PR_GI(K,iIon)=&
             State_GV(K,iU_I(iIon))**2*State_GV(K,iRho_I(iIon))
        
        if (iIon /= nIon) &
             PR(K)=PR(K)-MassElecIon_I(iIon)*(State_GV(K,iP_I(iIon))+PR_GI(K,iIon))
     enddo
  enddo
  
  !PR = sum( M_e/M_species *(P_thermal(k) + rho(k) * u(k)^2))
  PR(0:nDim+1)=State_GV(0:nDim+1,pE_)+PR_GI(0:nDim+1,nIon)+PR(0:nDim+1)
  
  EZ(1:nDim)=Source_CV(1:nDim,uE_)
  
  do iIon=1,nIon-1
     DO K=1,NDIM
        
        !     dM_e/dt-sum(dM_i/dt * m_e/m_i) :: dM_i/dt = momentum collision term
        !     sum(m_e/m_i * (u_e - u_i) * S_i :: S_i is the mass density source
        EZ(K) = EZ(K)&
             -MassElecIon_I(iIon)*Source_CV(K,iU_I(iIon))&
             +MassElecIon_I(iIon)&
             *(State_GV(K,uE_)-State_GV(K,iU_I(iIon)))*Source_CV(K,iRho_I(iIon))

     enddo
  enddo

  DZ=0.5/DRBND
  DO k=1,NDIM
     
     !dnom = q*n_e
     eNumDens=0.0
     do iIon=1,nIon-1
        eNumDens=eNumDens+State_GV(k,iRho_I(iIon))/Mass_I(iIon)
     enddo
     
     DNOM = 4.803242E-10*eNumDens
     
     !      A'/A * [u_e^2*rho_e - sum(m_e/m_i * u_i^2*rho_i)]
     EZ1=PR_GI(k,nIon)
     do iIon=1,nIon-1
        EZ1=EZ1-MassElecIon_I(iIon)*PR_GI(k,iIon)
     enddo
     EZ1=-DAREA(k)*EZ1
     
     !     E_z = -1/(e*n_e) * d/dz [ sum(m_e/m_s(p_s + u_s^2*rho_s)) 
     !           -sum_ions ( m_e/m_i*( (u_e-u_i)*S_i - dM_i/dt ) ) - dM_e/dt ] 
     !           -1/(e*n_e) * A'/A * ( sum( m_e/m_s * u_s^2 * rho_s ) )
     ! note: derivatives are central difference


     
     EFIELD(K) = ( EZ1 + DZ*(EZ(k+1)-EZ(k-1)) - DZ*(PR(k+1)-PR(k-1)) )&
          /DNOM
     
     !         EFIELD(K) = -DZ*(PR(k+1)-PR(k-1))/DNOM
     
     
  enddo
  EFIELD(NDIM)=2.13E-7*CURRMX/State_GV(NDIM,Te_)**1.5
  
  IF (EFIELD(NDIM).LT.0.) EFIELD(NDIM)=0.
  
  !     Update the momentum and energy source terms now that the 
  !     Efield is calculated.
  
  do iIon=1,nIon-1
     
     aCentrifugal_C(:) = 0.0
     if (UseCentrifugal) then
        call calc_centrifugal(nDim, State_GV(1:nDim,iU_I(iIon)),aCentrifugal_C)
     endif

     if (UseWaveAcceleration) then
        call calc_wave_acceleration
     endif

     DO K=1,NDIM
        Source_CV(K,iU_I(iIon))=Source_CV(K,iU_I(iIon))&
             +Source_CV(K,iRho_I(iIon))*State_GV(K,iU_I(iIon))&
             +State_GV(K,iRho_I(iIon))&
             *(aCentrifugal_C(K)+WaveAcceleration_C(K) &
              +GRAVTY(K)+4.803242E-10*EFIELD(K)/Mass_I(iIon))
        
        Source_CV(K,iP_I(iIon))=Source_CV(K,iP_I(iIon)) &
             +State_GV(K,iU_I(iIon))&
             *(Source_CV(K,iU_I(iIon))-0.5*State_GV(K,iU_I(iIon))*Source_CV(K,iRho_I(iIon)))
     enddo
  enddo

  if (allocated(aCentrifugal_C)) deallocate (aCentrifugal_C)  
  
  RETURN
  
END SUBROUTINE PW_CALC_EFIELD
