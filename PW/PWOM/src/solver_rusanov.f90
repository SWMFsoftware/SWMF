!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!\
! -----------------------------------------------------------------------------
! Here we update using the rusanov method
! ----------------------------------------------------------------------------
!/

subroutine rusanov_solver(iIon, nCell,&
     Rgas, DtIn,&
     OldState_GV,&
     OldeState_GV,&
     RhoSource_C, RhoUSource_C, eSource_C,&
     HeatCon_G,UpdateState_GV)
  
  use ModCommonVariables
  use ModPWOM, ONLY: IsFullyImplicit,UseIonHeat, TypeFlux 
  implicit none

  integer, intent(in)      :: iIon,nCell
  real, intent(in)         :: Rgas,DtIn
  real, intent(in)         :: OldState_GV(-1:nCell+2,3)
  real, intent(in)         :: OldeState_GV(-1:nCell+2,3)
  real, dimension(nCell), intent(in)  :: RhoSource_C, RhoUSource_C, eSource_C
  real, intent(in)  :: HeatCon_G(0:nCell+1)
  real,          intent(out) :: UpdateState_GV(-1:nCell+2,3)
  
  
  integer,parameter    :: Rho_=1,U_=2,P_=3,T_=4
  real, allocatable    :: LeftRho_F(:), LeftU_F(:), LeftP_F(:)
  real, allocatable    :: LefteRho_F(:), LefteP_F(:)
  real, allocatable    :: RightRho_F(:), RightU_F(:), RightP_F(:)
  real, allocatable    :: RighteRho_F(:), RighteP_F(:)
  real, allocatable    :: RhoFlux_F(:), RhoUFlux_F(:), eFlux_F(:)
  real, allocatable    :: T_G(:),Heat1_G(:),HeatSource_C(:),eTotalSource_C(:),&
                          Heat1_F(:),GradT_F(:),Diffusion_C(:),&
                          Conduction_F(:)
  integer :: i,iCell
  !---------------------------------------------------------------------------


  if (.not.allocated(LeftRho_F)) then
     allocate(LeftRho_F(nCell+1), LeftU_F(nCell+1), LeftP_F(nCell+1),&
          RightRho_F(nCell+1), RightU_F(nCell+1), RightP_F(nCell+1),&
          LefteRho_F(nCell+1), LefteP_F(nCell+1),&
          RighteRho_F(nCell+1), RighteP_F(nCell+1),&
          RhoFlux_F(nCell+1), RhoUFlux_F(nCell+1), eFlux_F(nCell+1),&
          T_G(-1:nCell+2),Heat1_G(-1:nCell+1),Heat1_F(0:nCell),     &
          HeatSource_C(nCell),eTotalSource_C(nCell),GradT_F(0:nCell),&
          Diffusion_C(nCell),Conduction_F(0:nCell))
  endif
  eTotalSource_C = eSource_C

  ! get temperature from pressure and density, and calculate heat flow
  ! source term for the fully implicit case.
  T_G(-1:nCell+2)=OldState_GV(-1:nCell+2,3)/Rgas/OldState_GV(-1:nCell+2,1)
  
  if (IsFullyImplicit .and. UseIonHeat)then
     ! Get Temperature Gradient
     do iCell = 0,nCell
        GradT_F(iCell)=(T_G(iCell+1)-T_G(iCell)) / DrBnd
     enddo
     ! Get facevalues of Kappa*Grad(T)
     do iCell = 0,nCell
        Conduction_F(iCell) = (HeatCon_G(iCell+1)+HeatCon_G(iCell)) / 2.0 &
             * GradT_F(iCell)
     enddo
     ! Get heat diffusion term
     do iCell = 1,nCell
        Diffusion_C(iCell) = Rgas/Gmin1 *&
           (Ar23(iCell)*Conduction_F(iCell)-Ar12(iCell)*Conduction_F(iCell-1))&
             / (CellVolume_C(iCell))
     enddo
     ! Add heat flux to energy source

     eTotalSource_C = eSource_C+Diffusion_C
  endif
  
  ! Save the old state
  UpdateState_GV=OldState_GV

  ! get the face values
  call calc_facevalues(nCell,OldState_GV(-1:nCell+2,Rho_),LeftRho_F,RightRho_F)
  call calc_facevalues(nCell,OldeState_GV(-1:nCell+2,Rho_),LefteRho_F,RighteRho_F)

  call calc_facevalues(nCell, OldState_GV(-1:nCell+2,U_), LeftU_F, RightU_F)

  call calc_facevalues(nCell, OldState_GV(-1:nCell+2,P_), LeftP_F, RightP_F)
  call calc_facevalues(nCell, OldeState_GV(-1:nCell+2,P_), LefteP_F, RighteP_F)

  !\
  ! Select the intercell flux for the Approximate Riemann Solver
  !/
  select case (TypeFlux)
     case('LaxFriedrichs')
        !get the LaxFriedrichs flux
        call laxfriedrichs_flux( DtIn,nCell,&
             LeftRho_F, LeftU_F, LeftP_F, &
             RightRho_F, RightU_F, RightP_F, &
             LefteRho_F, LefteP_F, &
             RighteRho_F, RighteP_F, &
             RhoFlux_F, RhoUFlux_F, eFlux_F)

     case('Rusanov')
        !get the rusanov flux
        call rusanov_flux( DtIn,nCell,&
             LeftRho_F, LeftU_F, LeftP_F, &
             RightRho_F, RightU_F, RightP_F, &
             LefteRho_F, LefteP_F, &
             RighteRho_F, RighteP_F, &
             RhoFlux_F, RhoUFlux_F, eFlux_F)

     case('HLL')
        !get the rusanov flux
        call HLL_flux( DtIn,nCell,&
             LeftRho_F, LeftU_F, LeftP_F, &
             RightRho_F, RightU_F, RightP_F, &
             LefteRho_F, LefteP_F, &
             RighteRho_F, RighteP_F, &
             RhoFlux_F, RhoUFlux_F, eFlux_F)
     end select


  ! update the cells one timestep 
  call update_state( iIon,nCell,&
       Rgas, DtIn,&
       OldState_GV(1:nCell,Rho_), OldState_GV(1:nCell,U_), &
       OldState_GV(1:nCell,P_), &
       LeftRho_F, RightRho_F,LeftU_F, RightU_F,LeftP_F, RightP_F, &
       RhoFlux_F, RhoUFlux_F, eFlux_F, &
       RhoSource_C, RhoUSource_C, eTotalSource_C,&
       UpdateState_GV(1:nCell,Rho_), UpdateState_GV(1:nCell,U_), &
       UpdateState_GV(1:nCell,P_), T_G(1:nCell))

      deallocate(LeftRho_F, LeftU_F, LeftP_F,&
          RightRho_F, RightU_F, RightP_F,&
          LefteRho_F, LefteP_F,&
          RighteRho_F, RighteP_F,&
          RhoFlux_F, RhoUFlux_F, eFlux_F,T_G,Heat1_G,Heat1_F,HeatSource_C,&
          eTotalSource_C,GradT_F,Diffusion_C,Conduction_F)

end subroutine rusanov_solver

!==============================================================================

subroutine laxfriedrichs_flux( DtIn, nCell,&
     LeftRho_F, LeftU_F, LeftP_F, &
     RightRho_F, RightU_F, RightP_F, &
     LefteRho_F, LefteP_F, &
     RighteRho_F, RighteP_F, &
     RhoFlux_F, RhoUFlux_F, eFlux_F)

  use ModCommonVariables
  implicit none

  real, intent(in)   :: DtIn
  integer,intent(in) :: nCell
  real, dimension(nCell+1), intent(in) :: LeftRho_F, LeftU_F, LeftP_F, &
                                          RightRho_F, RightU_F, RightP_F

  real, dimension(nCell+1), intent(in) :: LefteRho_F, LefteP_F, &
                                          RighteRho_F, RighteP_F
  real, dimension(nCell+1), intent(out) :: RhoFlux_F, RhoUFlux_F, eFlux_F
  real    :: GammaOverGammaMinus1, InvGammaMinus1, RightE, LeftE, Coeff
  integer :: i
  real :: aIonLEft, aIonRight,aElecLeft,aElecRight,aplusLeft,aplusRight
  !----------------------------------------------------------------------------

  InvGammaMinus1       = 1.0/(Gamma - 1.0)
  GammaOverGammaMinus1 = Gamma/(Gamma - 1.0)

!
!  Coeff = 0.5/dtr1
  
  do i=1,nCell+1
     ! Take half of the maximum of left and right hydro propagation speeds
     Coeff = 0.47*DRBND/DtIn
     
     RhoFlux_F(i)  = 0.5*(LeftRho_F(i) *LeftU_F(i) &
          +               RightRho_F(i)*RightU_F(i))  &
          - Coeff * (RightRho_F(i) - LeftRho_F(i))

     RhoUFlux_F(i) = 0.5*(LeftRho_F(i) *LeftU_F(i)**2  + LeftP_F(i) &
          +               RightRho_F(i)*RightU_F(i)**2 + RightP_F(i))  &
          - Coeff * (RightRho_F(i)*RightU_F(i) - LeftRho_F(i)*LeftU_F(i))

     
     RightE = 0.5*RightRho_F(i)*RightU_F(i)**2 + InvGammaMinus1*RightP_F(i)
     LeftE  = 0.5*LeftRho_F(i)*LeftU_F(i)**2   + InvGammaMinus1*LeftP_F(i)

     eFlux_F(i)    = 0.5*(0.5*LeftRho_F(i) *LeftU_F(i)**3  &
          + GammaOverGammaMinus1*LeftU_F(i)* LeftP_F(i) &
          +               0.5*RightRho_F(i)*RightU_F(i)**3 &
          + GammaOverGammaMinus1*RightU_F(i)*RightP_F(i))  &
          - Coeff * (RightE - LeftE)
  end do
  
end subroutine laxfriedrichs_flux

!==============================================================================
!==============================================================================

subroutine rusanov_flux( DtIn, nCell,&
     LeftRho_F, LeftU_F, LeftP_F, &
     RightRho_F, RightU_F, RightP_F, &
     LefteRho_F, LefteP_F, &
     RighteRho_F, RighteP_F, &
     RhoFlux_F, RhoUFlux_F, eFlux_F)

  use ModCommonVariables
  implicit none

  real, intent(in)   :: DtIn
  integer,intent(in) :: nCell
  real, dimension(nCell+1), intent(in) :: LeftRho_F, LeftU_F, LeftP_F, &
                                          RightRho_F, RightU_F, RightP_F

  real, dimension(nCell+1), intent(in) :: LefteRho_F, LefteP_F, &
                                          RighteRho_F, RighteP_F
  real, dimension(nCell+1), intent(out) :: RhoFlux_F, RhoUFlux_F, eFlux_F
  real    :: GammaOverGammaMinus1, InvGammaMinus1, RightE, LeftE, Coeff
  integer :: i
  real :: aIonLEft, aIonRight,aElecLeft,aElecRight,aplusLeft,aplusRight
  !----------------------------------------------------------------------------

  InvGammaMinus1       = 1.0/(Gamma - 1.0)
  GammaOverGammaMinus1 = Gamma/(Gamma - 1.0)

  do i=1,nCell+1
     ! Find the max characteristic speeds from the ion and electron 
     ! accustic speeds. From Gombosi and Rasmussen 1991 (eq 64)
     aIonLeft=sqrt(LeftP_F(i)/LeftRho_F(i))
     aIonRight=sqrt(RightP_F(i)/RightRho_F(i))
     aElecLeft=sqrt(LefteP_F(i)/LefteRho_F(i))
     aElecRight=sqrt(RighteP_F(i)/RighteRho_F(i))
     
     ! Find the max characteristic speeds 
     aplusLeft = sqrt((aElecLeft**2.0+6.0*aIonLeft**2.0 &
          +sqrt(aElecLeft**4.0+24.0*aIonLeft**4))*0.5)
     aplusRight = sqrt((aElecRight**2.0+6.0*aIonRight**2.0 &
          +sqrt(aElecRight**4.0+24.0*aIonRight**4))*0.5)
     
     Coeff = &
          0.5*max( abs(LeftU_F(i)) + aplusLeft, &
          abs(RightU_F(i)) + aplusRight)

     ! Get the intercell fluxes
     RhoFlux_F(i)  = 0.5*(LeftRho_F(i) *LeftU_F(i) &
          +               RightRho_F(i)*RightU_F(i))  &
          - Coeff * (RightRho_F(i) - LeftRho_F(i))

     RhoUFlux_F(i) = 0.5*(LeftRho_F(i) *LeftU_F(i)**2  + LeftP_F(i) &
          +               RightRho_F(i)*RightU_F(i)**2 + RightP_F(i))  &
          - Coeff * (RightRho_F(i)*RightU_F(i) - LeftRho_F(i)*LeftU_F(i))

     
     RightE = 0.5*RightRho_F(i)*RightU_F(i)**2 + InvGammaMinus1*RightP_F(i)
     LeftE  = 0.5*LeftRho_F(i)*LeftU_F(i)**2   + InvGammaMinus1*LeftP_F(i)

     eFlux_F(i)    = 0.5*(0.5*LeftRho_F(i) *LeftU_F(i)**3  &
          + GammaOverGammaMinus1*LeftU_F(i)* LeftP_F(i) &
          +               0.5*RightRho_F(i)*RightU_F(i)**3 &
          + GammaOverGammaMinus1*RightU_F(i)*RightP_F(i))  &
          - Coeff * (RightE - LeftE)
  end do
  
end subroutine rusanov_flux

!==============================================================================
!==============================================================================

subroutine HLL_flux( DtIn, nCell,&
     LeftRho_F, LeftU_F, LeftP_F, &
     RightRho_F, RightU_F, RightP_F, &
     LefteRho_F, LefteP_F, &
     RighteRho_F, RighteP_F, &
     RhoFlux_F, RhoUFlux_F, eFlux_F)

  use ModCommonVariables
  implicit none

  real, intent(in)   :: DtIn
  integer,intent(in) :: nCell
  real, dimension(nCell+1), intent(in) :: LeftRho_F, LeftU_F, LeftP_F, &
                                          RightRho_F, RightU_F, RightP_F

  real, dimension(nCell+1), intent(in) :: LefteRho_F, LefteP_F, &
                                          RighteRho_F, RighteP_F
  real, dimension(nCell+1), intent(out) :: RhoFlux_F, RhoUFlux_F, eFlux_F
  real    :: GammaOverGammaMinus1, InvGammaMinus1, RightE, LeftE, Coeff
  integer :: i
  real :: aIonLEft, aIonRight,aElecLeft,aElecRight,aplusLeft,aplusRight
  real :: Sleft, Sright, FL,FR,UL,UR
  !----------------------------------------------------------------------------

  InvGammaMinus1       = 1.0/(Gamma - 1.0)
  GammaOverGammaMinus1 = Gamma/(Gamma - 1.0)

  do i=1,nCell+1
     
     
     ! From Gombosi and Rasmussen 1991 (eq 64)
     aIonLeft=sqrt(LeftP_F(i)/LeftRho_F(i))
     aIonRight=sqrt(RightP_F(i)/RightRho_F(i))
     aElecLeft=sqrt(LefteP_F(i)/LefteRho_F(i))
     aElecRight=sqrt(RighteP_F(i)/RighteRho_F(i))
     
     aplusLeft = sqrt((aElecLeft**2.0+6.0*aIonLeft**2.0 &
          +sqrt(aElecLeft**4.0+24.0*aIonLeft**4))*0.5)
     aplusRight = sqrt((aElecRight**2.0+6.0*aIonRight**2.0 &
          +sqrt(aElecRight**4.0+24.0*aIonRight**4))*0.5)

     ! set the Left and right signal speeds from Toro p. 328
     Sleft = min(LeftU_F(i)- aplusLeft,RightU_F(i)- aplusRight)
     Sright= max(LeftU_F(i)+ aplusLeft,RightU_F(i)+ aplusRight)

     ! get HLL fluxes (Harten, Lax and van Leer, 1983)    
     ! Toro p. 321
     UR=RightRho_F(i)
     UL=LeftRho_F(i)
     FL=LeftRho_F(i) *LeftU_F(i)
     FR=RightRho_F(i)*RightU_F(i)
     
     if(Sleft>=0.0) then
        RhoFlux_F(i)  = FL
     elseif(Sleft<=0.0 .and. Sright>=0.0)then
        RhoFlux_F(i)  = (Sright*FL-Sleft*FR+Sleft*Sright*(UR-UL))/(Sright-Sleft)
     elseif(Sright<=0.0) then
        RhoFlux_F(i)  = FR
     endif

     UR=RightRho_F(i)*RightU_F(i)
     UL=LeftRho_F(i)*LeftU_F(i)
     FL=LeftRho_F(i) *LeftU_F(i)**2  + LeftP_F(i)
     FR=RightRho_F(i)*RightU_F(i)**2 + RightP_F(i)

     if(Sleft>=0.0) then
        RhoUFlux_F(i) = FL
     elseif(Sleft<=0.0 .and. Sright>=0.0)then
        RhoUFlux_F(i) = (Sright*FL-Sleft*FR+Sleft*Sright*(UR-UL))/(Sright-Sleft)
     elseif(Sright<=0.0) then
        RhoUFlux_F(i) = FR
     endif

     RightE = 0.5*RightRho_F(i)*RightU_F(i)**2 + InvGammaMinus1*RightP_F(i)
     LeftE  = 0.5*LeftRho_F(i)*LeftU_F(i)**2   + InvGammaMinus1*LeftP_F(i)

     UR=RightE
     UL=LeftE
     FL=0.5*LeftRho_F(i) *LeftU_F(i)**3  &
          + GammaOverGammaMinus1*LeftU_F(i)* LeftP_F(i)
     FR= 0.5*RightRho_F(i)*RightU_F(i)**3 &
          + GammaOverGammaMinus1*RightU_F(i)*RightP_F(i)
     

     if(Sleft>=0.0) then
         eFlux_F(i) = FL
     elseif(Sleft<=0.0 .and. Sright>=0.0)then
         eFlux_F(i) = (Sright*FL-Sleft*FR+Sleft*Sright*(UR-UL))/(Sright-Sleft)
     elseif(Sright<=0.0) then
         eFlux_F(i) = FR
     endif

  end do
  
end subroutine HLL_flux

!==============================================================================

subroutine update_state( iIon,nCell,&
     Rgas,DtIn,&
     OldRho_C, OldU_C, OldP_C, &
     LeftRho_F, RightRho_F,LeftU_F, RightU_F,LeftP_F, RightP_F, &
     RhoFlux_F, RhoUFlux_F, eFlux_F, &
     RhoSource_C, RhoUSource_C, eSource_C,&
     NewRho_C, NewU_C, NewP_C, NewT_C)

  use ModCommonVariables
  use ModPWOM, ONLY: IsPointImplicit, IsPointImplicitAll
  implicit none

  integer, intent(in)                   :: iIon,nCell  
  real, intent(in)                      :: Rgas,DtIn 
  real, dimension(nCell), intent(in)  :: OldRho_C, OldU_C, OldP_C
  real, dimension(nCell+1), intent(in):: LeftRho_F, LeftU_F, LeftP_F
  real, dimension(nCell+1), intent(in):: RightRho_F, RightU_F, RightP_F
  real, dimension(nCell+1), intent(in):: RhoFlux_F, RhoUFlux_F, eFlux_F
  real, dimension(nCell), intent(in)  :: RhoSource_C, RhoUSource_C, eSource_C
  real, dimension(nCell), intent(out) :: NewRho_C, NewU_C, NewP_C,NewT_C
  real, allocatable :: OldT_C(:), NewRhoU(:), NewE(:)
  real :: InvGammaMinus1,SumCollisionFreq,SumHTxCol,SumMFxCol
  real, allocatable :: NewRhoUSource_C(:), ExplicitRhoU_C(:)
  real :: CollisionNeutral1, CollisionNeutral2, &
          CollisionNeutral3,Collisionneutral4
  real :: CollisionCoefT1, CollisionCoefM1, CollisionCoefT2, CollisionCoefM2,&
          CollisionCoefT3, CollisionCoefM3, CollisionCoefT4, CollisionCoefM4
  real :: ModifyESource
  integer :: i,iSpecies,MinSpecies
  !----------------------------------------------------------------------------
  if (.not. allocated(NewRhoU)) then
     allocate(NewRhoU(nCell), NewE(nCell),NewRhoUSource_C(nCell), &
          ExplicitRhoU_C(nCell),OldT_C(nCell))
  endif
  OldT_C(1:nCell)=OldP_C(1:nCell)/Rgas/OldRho_C(1:nCell)
  
  InvGammaMinus1 = 1.0/(Gamma - 1.0)

  do i=1,nCell
 !    ImpliciteSource_C(i)=0.0
     NewRho_C(i) = OldRho_C(i) &
          - DtIn * (Ar23(i)*RhoFlux_F(i+1)-Ar12(i)*RhoFlux_F(i)) &
          / CellVolume_C(i) &
          + DtIn*RhoSource_C(i)
     NewRhoU(i) = OldRho_C(i)*OldU_C(i) &
          - DtIn * (Ar23(i)*RhoUFlux_F(i+1)-Ar12(i)*RhoUFlux_F(i)) &
          / CellVolume_C(i) &
          + DtIn*Darea(i)*OldP_C(i) &
          + DtIn*RhoUSource_C(i)

     NewE(i) = 0.5*OldRho_C(i)*OldU_C(i)**2 + InvGammaMinus1*OldP_C(i) &
          - DtIn * (Ar23(i)*eFlux_F(i+1)-Ar12(i)*eFlux_F(i)) &
          / CellVolume_C(i) &
          + DtIn*eSource_C(i) 
 
     NewU_C(i) = NewRhoU(i)/NewRho_C(i)

     NewP_C(i) = gMin1*(NewE(i) - 0.5 * NewU_C(i)**2 * NewRho_C(i))

     If (IsPointImplicit) then
        !If is ImplictAll = T then treat ions and neutral contributions
        !with point implicit. Else Just treat neutrals with point implicit
        If (IsPointImplicitAll) then
           MinSpecies = 1
        else
           MinSpecies = 5
        endif
        SumCollisionFreq = 0.0
        SumHTxCol        = 0.0
        SumMFxCol        = 0.0
        do iSpecies=MinSpecies,nSpecies
           if (iSpecies /= iIon) then
              SumCollisionFreq = &
                   SumCollisionFreq+CollisionFreq_IIC(iIon,iSpecies,i)
              SumHTxCol        = SumHTxCol &
                   + HeatFlowCoef_II(iIon,iSpecies)*CollisionFreq_IIC(iIon,iSpecies,i)
              SumMFxCol        = SumMFxCol &
                   + MassFracCoef_II(iIon,iSpecies)*CollisionFreq_IIC(iIon,iSpecies,i)
           endif
        enddo
        
        !eliminate contribution of implicit sources in RhoU
        ExplicitRhoU_C(i) = &
             NewRhoU(i)+DtIn*OldU_C(i)*OldRho_C(i)*SumCollisionFreq
        

        NewRhoU(i)= ExplicitRhoU_C(i) / &
             (1+DtIn*SumCollisionFreq)
        
        NewU_C(i) = NewRhoU(i)/NewRho_C(i)

        !When making p implicit we have to take the cf*cm*Ti term out of 
        !esource. Do this by adding ModifyEsource*dt to esource.
        ModifyESource = OldRho_C(i)* &
             SumHTxCol*OldT_C(i)

        NewE(i)   = NewE(i) + DtIn*ModifyESource

        NewP_C(i) = gMin1*(NewE(i) - 0.5 * NewU_C(i)**2 * NewRho_C(i))
        
        NewP_C(i) = NewP_C(i) / (1+DtIn*2.*SumMFxCol)

     endif
     
     
     NewT_C(i)=NewP_C(i)/Rgas/NewRho_C(i)
  end do
  deallocate(NewRhoU, NewE,NewRhoUSource_C, ExplicitRhoU_C,OldT_C)
end subroutine update_state
