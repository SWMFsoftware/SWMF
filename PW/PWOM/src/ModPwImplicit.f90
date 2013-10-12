!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModPwImplicit
  use ModCommonPlanet, ONLY: iRho,iU,iP,iTe_
  implicit none


  real, allocatable :: StateOrig_GV(:,:),StateFull_GV(:,:)
  contains
  !============================================================================

  subroutine PW_calc_residual(nOrder, DtIn, nCell, nVar, StateIn_GV, Resid_CV)

    use ModCommonVariables, ONLY: iRho_I,iU_I,iP_I,iT_I,Source_CV,RGAS_I, &
         MassElecIon_I,nIon,RhoE_,uE_,Te_,pE_,Curr,CURRMX,HeatCon_GI,&
         nVarFull=>nVar,TypeSolver
    use ModPWOM, ONLY: Beta, BetaIn

    implicit none
    integer, intent(in) :: nOrder, nCell, nVar
    real,    intent(in) :: DtIn
    real,    intent(in) :: StateIn_GV(-1:nCell+2, nVar)
    real,    intent(out):: Resid_CV(nCell, nVar)

    real, allocatable    :: NewState_GV(:,:),NewStateFull_GV(:,:)
    integer              :: iIon,k
    !--------------------------------------------------------------------------

    if (.not. allocated(NewState_GV)) &
         allocate(NewState_GV(-1:nCell+2, nVar))
    if (.not. allocated(NewStateFull_GV)) &
         allocate(NewStateFull_GV(-1:nCell+2, nVarFull))

    NewState_GV = StateIn_GV
    
    !Update source terms. First fill in complete state array

    ! Put the updated ReducedState back into the State array,
    ! and get temperature
    do iIon = 1,nIon-1
       StateFull_GV(:,iRho_I(iIon))= StateIn_GV(:,iRho(iIon))
       StateFull_GV(:,iU_I  (iIon))= StateIn_GV(:,iU  (iIon))
       StateFull_GV(:,iP_I  (iIon))= StateIn_GV(:,iP  (iIon))
       StateFull_GV(:,iT_I  (iIon))= &
            StateFull_GV(:,iP_I(iIon))/Rgas_I(iIon)/StateFull_GV(:,iRho_I(iIon))
    enddo
    StateFull_GV(:,iT_I(nIon))= StateIn_GV(:,iTe_)
    
    ! Set electron density and velocity and pressure
    StateFull_GV(1:nCell+2,RhoE_)=0.0
    StateFull_GV(1:nCell+2,uE_)  =0.0
    StateFull_GV(-1:0,RhoE_)=StateOrig_GV(-1:0,RhoE_)
    StateFull_GV(-1:0,uE_)  =StateOrig_GV(-1:0,uE_)
    do k=1,nCell+2
       do iIon=1,nIon-1
          StateFull_GV(k,RhoE_) = &
               StateFull_GV(k,RhoE_)+MassElecIon_I(iIon)*StateFull_GV(k,iRho_I(iIon))
          StateFull_GV(k,uE_)= &
               StateFull_GV(k,uE_)+ &
               (MassElecIon_I(iIon)*StateFull_GV(k,iRho_I(iIon))&
               *StateFull_GV(k,iU_I(iIon)))
       enddo
       if (k <= nCell) then
          StateFull_GV(k,uE_)=(StateFull_GV(k,uE_) -1.8965E-18*CURR(k))&
               /StateFull_GV(k,RhoE_)
       elseif (k > nCell) then
          StateFull_GV(k,uE_)=(StateFull_GV(k,uE_) -1.8965E-18*CURRMX)&
               /StateFull_GV(k,RhoE_)
       endif
    enddo

    StateFull_GV(:,iP_I(nIon))= &
         Rgas_I(nIon)*StateFull_GV(:,iRho_I(nIon))*StateIn_GV(:,iTe_)
    
    call COLLIS(nCell,StateFull_GV)
    
    call PW_calc_efield(nCell,StateFull_GV)         

    if(nOrder == 1)then
       Beta = 0.0
    else
       Beta = BetaIn
    end if
    
    if (TypeSolver == 'Rusanov') then
       do iIon=1,nIon-1
          call rusanov_solver(iIon,nCell,RGAS_I(iIon),DtIn,   &
               StateIn_GV(-1:nCell+2,iRho(iIon):iP(iIon)),&
               Source_CV(1:nCell,iRho_I(iIon)), Source_CV(1:nCell,iU_I(iIon)),&
               Source_CV(1:nCell,iP_I(iIon)),&
               HeatCon_GI(0:nCell+1,iIon), &
               NewState_GV(-1:nCell+2,iRho(iIon):iP(iIon)))
       enddo
    elseif (TypeSolver == 'Godunov') then
       do iIon=1,nIon-1
          CALL Solver(iIon,nCell,DtIn,&
               StateFull_GV(-1:nCell+2,iRho_I(iIon):iT_I(iIon)),&
               Source_CV(1:nCell,iRho_I(iIon)),Source_CV(1:nCell,iP_I(iIon)),&
               Source_CV(1:nCell,iU_I(iIon)),&
               RGAS_I(iIon),HeatCon_GI(0:nCell+1,iIon),&
               NewStateFull_GV(-1:nCell+2,iRho_I(iIon):iT_I(iIon)))
          NewState_GV(-1:nCell+2,iRho(iIon):iP(iIon))=&
               NewStateFull_GV(-1:nCell+2,iRho_I(iIon):iP_I(iIon))
       enddo
    else
       write(*,*) 'ERROR_PW: TypeSolver not recognized'
    endif
    call PW_eheat_flux_explicit(nCell,RGAS_I(nIon),DtIn,      &
         StateFull_GV(-1:nCell+2,iRho_I(nIon):iP_I(nIon)),       &
         Source_CV(1:nCell,iRho_I(nIon)), Source_CV(1:nCell,iU_I(nIon)),&
         Source_CV(1:nCell,iP_I(nIon)),                           &
         HeatCon_GI(0:nCell+1,nIon),                        &
         NewState_GV(-1:nCell+2,iTe_))
    
    Resid_CV(1:nCell,:) = NewState_GV(1:nCell,:)

    deallocate(NewState_GV)

    Resid_CV(1:nCell,:)=Resid_CV(1:nCell,:)-StateIn_GV(1:nCell,:)

  end subroutine PW_calc_residual

  !============================================================================

  subroutine PW_update_boundary(nCell, nVar, State_GV)

    use ModCommonVariables, ONLY: iRho_I,iU_I,iP_I,iT_I,Source_CV,RGAS_I, &
         MassElecIon_I,RhoE_,uE_,Te_,Pe_,nIon,DrBnd,&
         Gravty,ELFXIN,Mass_I,CURRMX,GAMMA,ETOP,&
         Ion1_,AR12top,AR23top
    use ModGmPressure
    use ModPWOM, ONLY: iLine
    use ModConst,ONLY: cBoltzmann
    implicit none
    integer, intent(in)    :: nCell, nVar
    real,    intent(inout) :: State_GV(-1:nCell+2, nVar)

    real    :: ScaleHeight_I(nIon-1),XHTM,T_C(nIon)
    integer :: iIon
    !--------------------------------------------------------------------------
    do iIon = 1,nIon-1
       
       ! Lower boundary is fixed, copy it from StateOrig_GV
       State_GV(-1:0,iRho(iIon):iP(iIon)) = &
            StateOrig_GV(-1:0,iRho_I(iIon):iP_I(iIon))
      
      State_GV(-1:0,iTe_) = &
           StateOrig_GV(-1:0,iT_I(nIon))

       ! Upper boundary requires temperature and scaleheight
       T_C(iIon) = &
            State_GV(nCell,iP(iIon))/Rgas_I(iIon)/State_GV(nCell,iRho(iIon))

       !State_GV(nCell+1:nCell+2,iU(iIon)) = State_GV(nCell,iU(iIon))
              
       ! cBoltzmann [cgs] = 1.0e7 * cBoltzmann [SI] 
       ScaleHeight_I(iIon) =&
            1.0e7*cBoltzmann*( T_C(iIon)+State_GV(nCell,iTe_) )&
            /(abs(Gravty(nCell))*Mass_I(iIon))

       !State_GV(nCell+1:nCell+2,iP(iIon)) = &
       !     State_GV(nCell,iP(iIon))*0.99999!*exp(-DrBnd/ScaleHeight_I(iIon))

!       State_GV(nCell+1,iP(iIon)) = &
!            State_GV(nCell  ,iP(iIon))*exp(-DrBnd/ScaleHeight_I(iIon))*.9
!       State_GV(nCell+2,iP(iIon)) = State_GV(nCell+1,iP(iIon))
!       State_GV(nCell+2,iP(iIon)) = &
!            (2.0*State_GV(nCell+1  ,iP(iIon))-State_GV(nCell,iP(iIon)))*0.9

!       State_GV(nCell+1,iP(iIon)) = &
!            (2.0*State_GV(nCell  ,iP(iIon))-State_GV(nCell-1,iP(iIon)))*0.9
!       State_GV(nCell+2,iP(iIon)) = &
!            (2.0*State_GV(nCell+1  ,iP(iIon))-State_GV(nCell,iP(iIon)))*0.9
!       State_GV(nCell+1,iP(iIon)) = &
!            (2.0*State_GV(nCell  ,iP(iIon))-State_GV(nCell-1,iP(iIon)))*0.99
!       State_GV(nCell+2,iP(iIon)) = State_GV(nCell+1,iP(iIon))
       State_GV(nCell+1,iP(iIon)) = &
            (2.0*State_GV(nCell  ,iP(iIon))-State_GV(nCell-1,iP(iIon)))*0.8
       State_GV(nCell+2,iP(iIon)) = State_GV(nCell+1,iP(iIon))
       

!       State_GV(nCell+2,iP(iIon)) = &
!            State_GV(nCell+1,iP(iIon))*.99999999999!*exp(-DrBnd/ScaleHeight_I(iIon))*0.1
       
       !set density bv to satisfy EOS
       State_GV(nCell+1,iRho(iIon)) = &
            State_GV(nCell+1,iP(iIon)) / Rgas_I(iIon) / T_C(iIon)

       State_GV(nCell+2,iRho(iIon)) = &
            State_GV(nCell+2,iP(iIon)) / Rgas_I(iIon) / T_C(iIon)

       !set velocity bc to conserve flux Rho(n+1) * U(n+1) = Rho(n) * U(n)
       !State_GV(nCell+1,iU(iIon)) = 2.0*State_GV(nCell,iU(iIon))-State_GV(nCell-1,iU(iIon))
       !State_GV(nCell+2,iU(iIon)) = 2.0*State_GV(nCell+1,iU(iIon))-State_GV(nCell,iU(iIon))
       State_GV(nCell+1:nCell+2,iU(iIon)) = State_GV(nCell,iU(iIon))
!       State_GV(nCell+1,iU(iIon)) = &
!            AR12top(1)*State_GV(nCell,iRho(iIon))* State_GV(nCell,iU(iIon))&
!            /(AR23top(1)*State_GV(nCell+1,iRho(iIon)))
!       State_GV(nCell+2,iU(iIon)) = &
!            AR12top(2)*State_GV(nCell+1,iRho(iIon))* State_GV(nCell+1,iU(iIon))&
!            /(AR23top(2)*State_GV(nCell+2,iRho(iIon)))
       
    enddo
    
    ! Set Te boundary
      
      State_GV(nCell+1,iTe_)=&
           State_GV(nCell,iTe_)+ETOP/State_GV(nCell,iTe_)**2.5
      State_GV(nCell+2,iTe_)=&
           State_GV(nCell+1,iTe_)+ETOP/State_GV(nCell,iTe_)**2.5

    ! If using GM --> PW coupling, overwrite ghostcell pressure
    ! with fraction of MHD pressure using 
    ! P_total = P_oxyg+P_e + smaller terms.
    ! n_e ~= n_oxyg hence P_oxyg ~= T_oxyg/(T_oxyg+T_e) * P_total
    if (UseGmToPw) then
       State_GV(nCell+1,iP(Ion1_)) = &
            (  T_C(Ion1_)&
            / ( T_C(Ion1_) + State_GV(nCell,iTe_) )&
            )*p_I(iLine)
       State_GV(nCell+1,iRho(Ion1_)) = &
            State_GV(nCell+1,iP(Ion1_)) / Rgas_I(Ion1_) / T_C(Ion1_)
    endif

!    call  PW_upper_heat_conduction

  end subroutine PW_update_boundary
  
  !============================================================================
  subroutine PW_implicit_update
    
    use ModLinearSolver,    ONLY: implicit_solver
    use ModCommonVariables, ONLY: Dt,DrBnd,nDim,nVar,State_GV,nIon,&
                                  RhoE_,uE_,iRho_I,iU_I,iP_I,iT_I,&
                                  MassElecIon_I,&
                                  Curr,nCell => nDim,Rgas_I,HeatCon_GI,Dt
    
    real                :: ImplPar = 1.0, DtExpl, DtImpl
    real,allocatable    :: ReducedState_GV(:,:)
    integer :: k,iIon, nVarReduced
    !--------------------------------------------------------------------------
    DtImpl = Dt
    DtExpl = 1.0e-3*(DrBnd/2.0e6)
!godunov    DtExpl = 1.0e-2
    !DtExpl = dx^2/2/max(kappa) 
    
   
!    DtExpl = &
!         0.5*DrBnd**2.0 &
!         *minval(State_GV(1:nCell,iRho_I(nIon))/HeatCon_GI(1:nCell,nIon))

    !write(*,*)'DtImpl, DtExpl=',DtImpl, DtExpl

    nVarReduced = 3*(nIon-1)+1

    if (.not.allocated(ReducedState_GV)) then
       allocate(ReducedState_GV(-1:nCell+2,nVarReduced),&
            StateOrig_GV(-1:nCell+2,nVar),StateFull_GV(-1:nCell+2,nVar))
    endif

    StateOrig_GV(-1:nCell+2,:) =  State_GV(-1:nCell+2,:)
    ! Create ReducedState array of only independent variables
    do iIon = 1,nIon-1
       ReducedState_GV(-1:nCell+2,iRho(iIon)) = State_GV(-1:nCell+2,iRho_I(iIon))
       ReducedState_GV(-1:nCell+2,iU  (iIon)) = State_GV(-1:nCell+2,iU_I  (iIon))
       ReducedState_GV(-1:nCell+2,iP  (iIon)) = State_GV(-1:nCell+2,iP_I  (iIon))
    enddo
    ReducedState_GV(-1:nCell+2,iTe_) = State_GV(-1:nCell+2,iT_I(nIon))

    ! Get the implicit update for the independent variables 
    call implicit_solver(ImplPar, DtImpl, DtExpl, nDim, nVarReduced, &
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
    State_GV(-1:nCell+2,iT_I(nIon)) = ReducedState_GV(-1:nCell+2,iTe_)
    
    ! Set electron density and velocity and pressure
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
    State_GV(-1:nCell+2,iP_I  (nIon)) = &
            State_GV(-1:nCell+2,iT_I(nIon))&
            *Rgas_I(iIon)*State_GV(-1:nCell+2,iRho_I(nIon))

    deallocate(ReducedState_GV,StateOrig_GV,StateFull_GV)
  end subroutine PW_implicit_update



end module ModPwImplicit
!==============================================================================
