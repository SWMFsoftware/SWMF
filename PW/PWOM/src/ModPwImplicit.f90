module ModPwImplicit

  implicit none

  real :: Te
  real, allocatable :: StateOrig_GV(:,:)
  contains
  !============================================================================

  subroutine PW_calc_residual(nOrder, DtIn, nCell, nVar, StateIn_GV, Resid_CV)

    use ModCommonVariables, ONLY: iRho_I,iU_I,iP_I,iT_I,Source_CV,RGAS_I, &
         MassElecIon_I,nIon,RhoE_,uE_,Te_,pE_,Curr
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
    
    !Update source terms. First fill in complete state array

    ! Put the updated ReducedState back into the State array,
    ! and get temperature
    do iIon = 1,nIon-1
       StateOrig_GV(:,iRho_I(iIon))= StateIn_GV(:,iRho(iIon))
       StateOrig_GV(:,iU_I  (iIon))= StateIn_GV(:,iU  (iIon))
       StateOrig_GV(:,iP_I  (iIon))= StateIn_GV(:,iP  (iIon))
       StateOrig_GV(:,iT_I  (iIon))= &
            StateOrig_GV(:,iP_I(iIon))/Rgas_I(iIon)/StateOrig_GV(:,iRho_I(iIon))
    enddo
        
    ! Set electron density and velocity
    StateOrig_GV(1:nCell,RhoE_)=0.0
    StateOrig_GV(1:nCell,uE_)  =0.0
    do k=1,nCell
       do iIon=1,nIon-1
          StateOrig_GV(k,RhoE_) = &
               StateOrig_GV(k,RhoE_)+MassElecIon_I(iIon)*StateOrig_GV(k,iRho_I(iIon))
          StateOrig_GV(k,uE_)= &
               StateOrig_GV(k,uE_)+ &
               (MassElecIon_I(iIon)*StateOrig_GV(k,iRho_I(iIon))&
               *StateOrig_GV(k,iU_I(iIon)))
       enddo

       StateOrig_GV(k,uE_)=(StateOrig_GV(k,uE_) -1.8965E-18*CURR(k))&
                               /StateOrig_GV(k,RhoE_)
    enddo
    
    call COLLIS(nCell,StateOrig_GV)

    call PW_calc_efield(nCell,StateOrig_GV)         


!    Source_CV(:,:) = 0.0
    
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
          call rusanov_solver(iIon,nCell,RGAS_I(iIon),DtIn,   &
               StateIn_GV(-1:nCell+2,iRho(iIon):iP(iIon)),&
               Source_CV(1:nCell,iRho_I(iIon)), Source_CV(1:nCell,iU_I(iIon)),&
               Source_CV(1:nCell,iP_I(iIon)),  &
               NewState_GV(-1:nCell+2,iRho(iIon):iP(iIon)))
       enddo
!    endif


    Resid_CV(1:nCell,:) = NewState_GV(1:nCell,:)

    deallocate(NewState_GV)

    Resid_CV(1:nCell,:)=Resid_CV(1:nCell,:)-StateIn_GV(1:nCell,:)

  end subroutine PW_calc_residual

  !=============================================================================

  subroutine PW_update_boundary(nCell, nVar, State_GV)

    use ModCommonVariables, ONLY: iRho_I,iU_I,iP_I,iT_I,Source_CV,RGAS_I, &
         MassElecIon_I,RhoE_,uE_,Te_,Pe_,nIon,DrBnd,&
         Gravty,ELFXIN,Mass_I,TIME,CURRMX,GAMMA,ETOP,&
         Ion1_
    use ModGmPressure
    use ModPWOM, ONLY: iLine
    use ModConst,ONLY: cBoltzmann
    implicit none
    integer, intent(in)    :: nCell, nVar
    real,    intent(inout) :: State_GV(-1:nCell+2, nVar)

    real    :: ScaleHeight_I(nIon-1),XHTM,T_C(nIon)
    integer :: iIon
    !---------------------------------------------------------------------------
    do iIon = 1,nIon-1
       T_C(iIon) = &
            State_GV(nCell,iP(iIon))/Rgas_I(iIon)/State_GV(nCell,iRho(iIon))

       State_GV(nCell+1:nCell+2,iU(iIon)) = State_GV(nCell,iU(iIon))
              
       ! cBoltzmann [cgs] = 1.0e7 * cBoltzmann [SI] 
       ScaleHeight_I(iIon) =&
            1.0e7*cBoltzmann*( T_C(iIon)+Te )&
            /(abs(Gravty(nCell))*Mass_I(iIon))

       State_GV(nCell+1:nCell+2,iP(iIon)) = &
            State_GV(nCell,iP(iIon))*exp(-DrBnd/ScaleHeight_I(iIon))

       State_GV(nCell+1,iRho(iIon)) = &
            State_GV(nCell+1,iP(iIon)) / Rgas_I(iIon) / T_C(iIon)

       State_GV(nCell+2,iRho(iIon)) = &
            State_GV(nCell+2,iP(iIon)) / Rgas_I(iIon) / T_C(iIon)
    enddo


    ! If using GM --> PW coupling, overwrite ghostcell pressure
    ! with fraction of MHD pressure using 
    ! P_total = P_oxyg+P_e + smaller terms.
    ! n_e ~= n_oxyg hence P_oxyg ~= T_oxyg/(T_oxyg+T_e) * P_total
    if (UseGmToPw) then
       State_GV(nCell+1,iP(Ion1_)) = &
            (  T_C(Ion1_)&
            / ( T_C(Ion1_) + Te )&
            )*p_I(iLine)
       State_GV(nCell+1,iRho(Ion1_)) = &
            State_GV(nCell+1,iP(Ion1_)) / Rgas_I(Ion1_) / T_C(Ion1_)
    endif

!    call  PW_upper_heat_conduction

  end subroutine PW_update_boundary
  
  !============================================================================
  subroutine PW_implicit_update
    
    use ModLinearSolver,    ONLY: implicit_solver
    use ModCommonVariables, ONLY: Dt,nDim,nVar,State_GV,nIon,&
                                  RhoE_,uE_,iRho_I,iU_I,iP_I,iT_I,&
                                  MassElecIon_I,&
                                  Curr,nCell => nDim,Rgas_I
    use ModPWOM, ONLY: Beta
    
    real                :: ImplPar = 0.005, DtExpl, DtImpl
    real,allocatable    :: ReducedState_GV(:,:)
    integer :: k,iIon, nVarReduced
    !--------------------------------------------------------------------------
    
    Beta = 1.5
    
    DtImpl = Dt
    DtExpl = 0.01
    nVarReduced = 3*(nIon-1)

    if (.not.allocated(ReducedState_GV)) then
       allocate(ReducedState_GV(-1:nCell+2,nVarReduced),&
            StateOrig_GV(-1:nCell+2,nVar))
    endif

    StateOrig_GV(-1:nCell+2,:) =  State_GV(-1:nCell+2,:)
    ! Create ReducedState array of only independent variables
    do iIon = 1,nIon-1
       ReducedState_GV(-1:nCell+2,iRho(iIon)) = State_GV(-1:nCell+2,iRho_I(iIon))
       ReducedState_GV(-1:nCell+2,iU  (iIon)) = State_GV(-1:nCell+2,iU_I  (iIon))
       ReducedState_GV(-1:nCell+2,iP  (iIon)) = State_GV(-1:nCell+2,iP_I  (iIon))
    enddo
    
    Te = State_GV(nCell,iT_I(nIon))
    
    ! Get the implicit update for the independent variables 
    call implicit_solver(ImplPar, DtExpl, DtImpl, nDim, nVarReduced, &
         ReducedState_GV(-1:nDim+2,:), &
         PW_calc_residual, PW_update_boundary)
    
    ! Put the updated ReducedState back into the State array,
    ! and get temperature
    do iIon = 1,nIon-1
       State_GV(-1:nCell+2,iRho_I(iIon)) = &
            ReducedState_GV(-1:nCell+2,iRho(iIon))
       State_GV(-1:nCell+2,iU_I  (iIon)) = &
            ReducedState_GV(-1:nCell+2,iU  (iIon))
       State_GV(-1:nCell+2,iP_I  (iIon)) = &
            ReducedState_GV(-1:nCell+2,iP  (iIon))
       State_GV(-1:nCell+2,iT_I  (iIon)) = &
            State_GV(-1:nCell+2,iP_I  (iIon))&
            /Rgas_I(iIon)/State_GV(-1:nCell+2,iRho_I(iIon))
    enddo
        
    ! Set electron density and velocity
    State_GV(1:nCell,RhoE_)=0.0
    State_GV(1:nCell,uE_)  =0.0
    do k=1,nCell
       do iIon=1,nIon-1
          State_GV(k,RhoE_) = &
               State_GV(k,RhoE_)+MassElecIon_I(iIon)*State_GV(k,iRho_I(iIon))
          State_GV(k,uE_)= &
               State_GV(k,uE_)+ &
               (MassElecIon_I(iIon)*State_GV(k,iRho_I(iIon))&
               *State_GV(k,iU_I(iIon)))
       enddo

       State_GV(k,uE_)=(State_GV(k,uE_) -1.8965E-18*CURR(k))/State_GV(k,RhoE_)
    enddo
    
    deallocate(ReducedState_GV,StateOrig_GV)
  end subroutine PW_implicit_update
  !============================================================================

  integer function iRho(iIon)
    integer, intent(in) :: iIon
    
    !--------------------------------------------------------------------------
    
    iRho = 3*(iIon-1) + 1
    return
  end function iRho
  
  !============================================================================

  integer function iU(iIon)
    integer, intent(in) :: iIon
    
    !--------------------------------------------------------------------------
    
    iU = 3*(iIon-1) + 2
    return
  end function iU
  
  !============================================================================

  integer function iP(iIon)
    integer, intent(in) :: iIon
    
    !--------------------------------------------------------------------------
    
    iP = 3*(iIon-1) + 3
    return
  end function iP


end module ModPwImplicit
!==============================================================================
