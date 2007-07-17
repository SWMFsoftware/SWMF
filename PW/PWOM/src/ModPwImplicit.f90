module ModPwImplicit

  implicit none

  contains
  !============================================================================

  subroutine PW_calc_residual(nOrder, DtIn, nCell, nVar, StateIn_GV, Resid_CV)

    use ModCommonVariables, ONLY: iRho_I,iU_I,iP_I,iT_I,Source_CV,RGAS_I,nDim, &
         MassElecIon_I,RhoE_,uE_,nIon,CURR
    implicit none
    integer, intent(in) :: nOrder, nCell, nVar
    real,    intent(in) :: DtIn
    real,    intent(in) :: StateIn_GV(-1:nCell+2, nVar)
    real,    intent(out):: Resid_CV(nCell, nVar)

    real, allocatable    :: NewState_GV(:,:)
    integer              :: iIon,k
    !---------------------------------------------------------------------------

    if (.not. allocated(NewState_GV)) allocate(NewState_GV(-1:nCell+2, nVar))

    NewState_GV = StateIn_GV
!    if (nOrder == 1) then
!       do iIon=1,nIon-1
!          call solver(iIon,DtIn, StateIn_GV(:,iRho_I(iIon):iT_I(iIon)),&
!               Source_CV(:,iRho_I(iIon)),Source_CV(:,iP_I(iIon)),&
!               Source_CV(:,iU_I(iIon)),&
!               RGAS_I(iIon),NewState_GV(:,iRho_I(iIon):iT_I(iIon)))
!       enddo
!
!    else if (nOrder == 2) then
       do iIon=1,nIon-1
          call rusanov_solver(iIon,RGAS_I(iIon),DtIn,   &
               StateIn_GV(:,iRho_I(iIon):iT_I(iIon)),&
               Source_CV(:,iRho_I(iIon)), Source_CV(:,iU_I(iIon)),&
               Source_CV(:,iP_I(iIon)),  &
               NewState_GV(:,iRho_I(iIon):iT_I(iIon)))
       enddo
!    endif

    Resid_CV(1:nDim,:) = NewState_GV(1:nDim,:)

    deallocate(NewState_GV)

    ! Set electron density and velocity
    Resid_CV(1:nDim,RhoE_)=0.0
    Resid_CV(1:nDim,uE_)  =0.0
    do k=1,nDim
       do iIon=1,nIon-1
          Resid_CV(K,RhoE_) = &
               Resid_CV(k,RhoE_)+MassElecIon_I(iIon)*Resid_CV(K,iRho_I(iIon))
          Resid_CV(K,uE_)= &
               Resid_CV(k,uE_)+ &
               (MassElecIon_I(iIon)*Resid_CV(K,iRho_I(iIon))&
               *Resid_CV(K,iU_I(iIon)))
       enddo

       Resid_CV(K,uE_)=(Resid_CV(K,uE_) -1.8965E-18*CURR(K))/Resid_CV(K,RhoE_)
    enddo

  end subroutine PW_calc_residual

  !=============================================================================

  subroutine PW_update_boundary(nCell, nVar, State_GV)

    use ModCommonVariables, ONLY: iRho_I,iU_I,iP_I,iT_I,Source_CV,RGAS_I,nDim, &
         MassElecIon_I,RhoE_,uE_,Te_,Pe_,nIon,DrBnd,&
         Gravty,ELFXIN,Mass_I,TIME,CURRMX,GAMMA,ETOP,&
         Ion1_
    use ModGmPressure
    use ModPWOM, ONLY: iLine
    use ModConst,ONLY: cBoltzmann
    implicit none
    integer, intent(in)    :: nCell, nVar
    real,    intent(inout) :: State_GV(-1:nCell+2, nVar)

    real    :: ScaleHeight_I(nIon-1),XHTM
    integer :: iIon
    !---------------------------------------------------------------------------
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

    XHTM=1.+ELFXIN*EXP(-(TIME-300.)**2/2./150./150.)

    ! Top boundary electron temperature is calculated from
    ! Te(nDim) + Tgradient*Dx where Tgradient = Etop/kappa
    ! and Etop is the electron heat flux or proportional to that
    State_GV(nDim+1,Te_)=State_GV(nDim,Te_)+XHTM*ETOP/State_GV(nDim,Te_)**2.5

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

    State_GV(nDim+1,RhoE_)=0.0
    State_GV(nDim+1,uE_)  =0.0

    do iIon=1,nIon-1
       State_GV(nDim+1,RhoE_) = &
            State_GV(nDim+1,RhoE_)+MassElecIon_I(iIon)*State_GV(nDim+1,iRho_I(iIon))
       State_GV(nDim+1,uE_)= &
            State_GV(nDim+1,uE_)+ &
            (MassElecIon_I(iIon)*State_GV(nDim+1,iRho_I(iIon))*State_GV(nDim+1,iU_I(iIon)))
    enddo
    State_GV(nDim+1,uE_)=(State_GV(nDim+1,uE_) -1.8965E-18*CURRMX)/State_GV(nDim+1,RhoE_)

    State_GV(nDim+1,pE_)=RGAS_I(nIon)*State_GV(nDim+1,Te_)*State_GV(nDim+1,RhoE_)  


    !  !these are sound speeds
    !  do iIon=1,nIon
    !     SoundSpeed_GI(nDim+1,iIon) = &
    !          SQRT(&
    !           GAMMA*State_GV(nDim+1,iP_I(iIon)) / State_GV(nDim+1,iRho_I(iIon))&
    !          )
    !  enddo

    call  PW_upper_heat_conduction

  end subroutine PW_update_boundary
  
  subroutine PW_implicit_update
    
    use PW_ModImplicit,     ONLY: PW_calc_residual, PW_update_boundary
    use ModLinearSolver,    ONLY: implicit_solver
    use ModCommonVariables, ONLY: Dt,nDim,nVar,State_GV
    
    implicit none
    real  :: ImplPar = 0.005, DtExpl, DtImpl
    
    !--------------------------------------------------------------------------
    DtImpl = Dt
    DtExpl = Dt * 0.1
    
    call implicit_solver(ImplPar, DtExpl, DtImpl, nDim, nVar, &
         State_GV(-1:nDim+2,1:nVar), &
         PW_calc_residual, PW_update_boundary)
    
  end subroutine PW_implicit_update
  
  
end module ModPwImplicit
!==============================================================================
