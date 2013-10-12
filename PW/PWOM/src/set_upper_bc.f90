!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine PW_set_upper_bc
  use ModCommonVariables
  use ModGmPressure
  use ModPWOM, ONLY: iLine
  use ModConst,ONLY: cBoltzmann

  real :: ScaleHeight_I(nIon-1)
  !----------------------------------------------------------------------------

  if (TypeSolver == 'Godunov') then
     
     do iIon = 1,nIon-1
        State_GV(nDim+1:nDim+2,iU_I(iIon)) = State_GV(nDim,iU_I(iIon))
        State_GV(nDim+1:nDim+2,iT_I(iIon)) = State_GV(nDim,iT_I(iIon))
        State_GV(nDim+1:nDim+2,iP_I(iIon)) = &
          (State_GV(nDim-1,iP_I(iIon))+State_GV(nDim,iP_I(iIon))&
           -State_GV(nDim-2,iP_I(iIon)))*0.8
        State_GV(nDim+1:nDim+2,iRho_I(iIon)) = &
          State_GV(nDim+1,iP_I(iIon)) / RGAS_I(iIon) / State_GV(nDim+1,iT_I(iIon))

!        State_GV(-1:0,iU_I(iIon))=&
!             State_GV(1,iU_I(iIon))*State_GV(1,iRho_I(iIon))/State_GV(-1:0,iRho_I(iIon))
     enddo
     
  ENDIF
  
  if (TypeSolver == 'Rusanov') then



     do iIon = 1,nIon-1
        State_GV(nDim+1:nDim+2,iU_I(iIon)) = State_GV(nDim,iU_I(iIon))
        State_GV(nDim+1:nDim+2,iT_I(iIon)) = State_GV(nDim,iT_I(iIon))

        ! cBoltzmann [cgs] = 1.0e7 * cBoltzmann [SI] 
        ScaleHeight_I(iIon) =&
             1.0e7*cBoltzmann*(State_GV(nDim,iT_I(iIon))+State_GV(nDim,Te_))&
             /(abs(Gravty(nDim))*Mass_I(iIon))
        
        State_GV(nDim+1:nDim+2,iP_I(iIon)) = &
             State_GV(nDim,iP_I(iIon))*exp(-DrBnd/ScaleHeight_I(iIon))
        
        State_GV(nDim+1,iRho_I(iIon)) = &
          State_GV(nDim+1,iP_I(iIon)) / RGAS_I(iIon) / State_GV(nDim+1,iT_I(iIon))

        State_GV(nDim+2,iRho_I(iIon)) = &
             State_GV(nDim+2,iP_I(iIon)) / RGAS_I(iIon) / State_GV(nDim+2,iT_I(iIon))
     enddo
  endif
  
  
      
  
  XHTM=1.+ELFXIN*EXP(-(TIME-300.)**2/2./150./150.)
  
  ! Top boundary electron temperature is calculated from
  ! Te(nDim) + Tgradient*Dx where Tgradient = Etop/kappa
  ! and Etop is the electron heat flux or proportional to that
  State_GV(nDim+1,Te_)=State_GV(nDim,Te_)+XHTM*ETOP/State_GV(nDim,Te_)**2.5
  State_GV(nDim+2,Te_)=State_GV(nDim+1,Te_)+XHTM*ETOP/State_GV(nDim,Te_)**2.5
  
  ! If using GM --> PW coupling, overwrite ghostcell pressure
  ! with fraction of MHD pressure using 
  ! P_total = P_oxyg+P_e + smaller terms.
  ! n_e ~= n_oxyg hence P_oxyg ~= T_oxyg/(T_oxyg+T_e) * P_total
  if (UseGmToPw) then
     State_GV(nDim+1,iP_I(Ion1_)) = &
         (  State_GV(nDim+1,iT_I(Ion1_))&
            / (State_GV(nDim+1,iT_I(Ion1_))+State_GV(nDim+1,Te_))&
          )*p_I(iLine)
     State_GV(nDim+1,iRho_I(Ion1_)) = &
          State_GV(nDim+1,iP_I(Ion1_)) / RGAS_I(Ion1_) / State_GV(nDim+1,iT_I(Ion1_))
  endif
  
  ! n_e   = n_H+  + n_O+ + n_He+
  ! rho_e = rho_H*(m_e/m_H) + rho_O*(m_e/m_O) + rho_He*(m_e/m_He)
  ! Set electron density and velocity and pressure

  State_GV(nDim+1:nDim+2,RhoE_)=0.0
  State_GV(nDim+1:nDim+2,uE_)  =0.0
  
  do iIon=1,nIon-1
     State_GV(nDim+1:nDim+2,RhoE_) = &
          State_GV(nDim+1:nDim+2,RhoE_)+MassElecIon_I(iIon)*State_GV(nDim+1:nDim+2,iRho_I(iIon))
     State_GV(nDim+1:nDim+2,uE_)= &
          State_GV(nDim+1:nDim+2,uE_)+ &
          ( MassElecIon_I(iIon)*State_GV(nDim+1:nDim+2,iRho_I(iIon))&
          *State_GV(nDim+1:nDim+2,iU_I(iIon)) )
  enddo
  State_GV(nDim+1:nDim+2,uE_)=&
       (State_GV(nDim+1:nDim+2,uE_) -1.8965E-18*CURRMX)/State_GV(nDim+1:nDim+2,RhoE_)


  State_GV(nDim+1:nDim+2,pE_)=&
       RGAS_I(nIon)*State_GV(nDim+1:nDim+2,Te_)*State_GV(nDim+1:nDim+2,RhoE_)  
  

  !these are sound speeds
  do iIon=1,nIon
     SoundSpeed_GI(nDim+1,iIon) = &
          SQRT(&
           GAMMA*State_GV(nDim+1,iP_I(iIon)) / State_GV(nDim+1,iRho_I(iIon))&
          )
  enddo

  call  PW_upper_heat_conduction

  RETURN
END SUBROUTINE PW_set_upper_bc
