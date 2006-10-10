!\
! -----------------------------------------------------------------------------
! Here we update using the rusanov method
! ----------------------------------------------------------------------------
!/

subroutine rusanov_solver(NameVariable, &
     Rgas, &
     OldRho_C, OldU_C, OldP_C, OldT_C, &
     RhoBottom, UBottom, PBottom, &
     RhoTop, UTop, Ptop, &
     RhoTop2, UTop2, Ptop2, &
     RhoSource_C, RhoUSource_C, eSource_C, &
     NewRho_C, NewU_C, NewP_C, NewT_C)
  
  use ModCommonVariables
  implicit none

  character(len=*), intent(in)          :: NameVariable
  real, intent(in)                      :: Rgas
  real, dimension(MaxGrid), intent(in)  :: OldRho_C, OldU_C, OldP_C, OldT_C
  real, intent(in)                      :: RhoBottom, UBottom, PBottom
  real, intent(in)                      :: RhoTop, UTop, PTop
  real, intent(in)                      :: RhoTop2, UTop2, PTop2
  real, dimension(MaxGrid), intent(in)  :: RhoSource_C, RhoUSource_C, eSource_C
  real, dimension(MaxGrid), intent(out) :: NewRho_C, NewU_C, NewP_C, NewT_C
  
  real, dimension(-1:MaxGrid+2) :: OldRho_G, OldU_G, OldP_G
  real, dimension(MaxGrid+1)    :: LeftRho_F, LeftU_F, LeftP_F
  real, dimension(MaxGrid+1)    :: RightRho_F, RightU_F, RightP_F
  real, dimension(MaxGrid+1)    :: RhoFlux_F, RhoUFlux_F, eFlux_F

  logical:: DoTest,DoTest2,DoTest3

  !---------------------------------------------------------------------------


  DoTest  = index(NameVariable,'H') > 0


  if (NamePlanet == 'Earth ') then
     DoTest  = index(NameVariable,'H') > 0
     DoTest2 = index(NameVariable,'O') > 0
     DoTest3 = index(NameVariable,'e') > 0
  endif




!  write(*,*)'rusanov solver starting for variable ',NameVariable, DoTest,DoTest2,DoTest3
 
  ! add real ghost cells for use in getting face values
  OldRho_G(1:nDim)      = OldRho_C(1:nDim)
  OldRho_G(-1:0)        = RhoBottom
  OldRho_G(nDim+1)      = RhoTop
  OldRho_G(nDim+2)      = RhoTop2

  OldU_G(1:nDim)        = OldU_C(1:nDim)
  OldU_G(-1:0)          = UBottom
  OldU_G(nDim+1)        = UTop
  OldU_G(nDim+2)        = UTop2

  OldP_G(1:nDim)        = OldP_C(1:nDim)
  OldP_G(-1:0)          = PBottom
  OldP_G(nDim+1)        = PTop
  OldP_G(nDim+2)        = PTop2


  ! get the face values
  call calc_facevalues(nDim, OldRho_G, LeftRho_F, RightRho_F)
  call calc_facevalues(nDim, OldU_G, LeftU_F, RightU_F)
  call calc_facevalues(nDim, OldP_G, LeftP_F, RightP_F)

  
  
  !get the rusanov flux
  call rusanov_flux( &
       LeftRho_F, LeftU_F, LeftP_F, &
       RightRho_F, RightU_F, RightP_F, &
       RhoFlux_F, RhoUFlux_F, eFlux_F)




  ! update the cells one timestep 
  call update_state( DoTest, DoTest2, DoTest3,&
       Rgas, &
       OldRho_C, OldU_C, OldP_C, OldT_C, &
       LeftRho_F, RightRho_F,LeftU_F, RightU_F,LeftP_F, RightP_F, &
       RhoFlux_F, RhoUFlux_F, eFlux_F, &
       RhoSource_C, RhoUSource_C, eSource_C, &
       NewRho_C, NewU_C, NewP_C, NewT_C)

!  if(DoTest)then
!!     write(*,*)'OldRho=',OldRho_C(nDim-1:nDim)
!!     write(*,*)'DHYD  =',DHYD(nDim-1:nDim)
!!     write(*,*)'OldP  =',OldP_C(nDim-1:nDim)
!!     write(*,*)'PHYD  =',PHYD(nDim-1:nDim)
!!     write(*,*)'THYD  =',THYD(nDim-1:nDim)
!!     write(*,*)'LeftP =',LeftP_F(nDim-1:nDim+1)
!!     write(*,*)'RightP=',RightP_F(nDim-1:nDim+1)
!!     write(*,*)'RhoFlux =',RhoFlux_F(nDim-1:nDim+1)
!!     write(*,*)'RhoUFlux=',RhoUFlux_F(nDim-1:nDim+1)
!!     write(*,*)'eFlux   =',eFlux_F(nDim-1:nDim+1)
!!     write(*,*)'LeftRho_F  =',LeftRho_F(nDim-1:nDim+1)
!!     write(*,*)'RightRho_F =',RightRho_F(nDim-1:nDim+1)
!!     write(*,*)'LeftU_F    =',LeftU_F(nDim-1:nDim+1)
!!     write(*,*)'RightU_F   =',RightU_F(nDim-1:nDim+1)
!!     write(*,*)'OldU_G     =',OldU_G(nDim-1:nDim+2)
!!     write(*,*)'LeftP_F    =',LeftP_F(nDim-1:nDim+1)
!!     write(*,*)'RightP_F   =',RightP_F(nDim-1:nDim+1)
!!     write(*,*)'NewRho=',NewRho_C(nDim-1:nDim)
!!     write(*,*)'NewP  =',NewP_C(nDim-1:nDim)
!!     write(*,*)'NewT  =',NewT_C(nDim-1:nDim)
!
!     write(*,*)'OldRho=',OldRho_C(1:3)
!     write(*,*)'DHYD  =',DHYD(1:3)
!     write(*,*)'OldP  =',OldP_C(1:3)
!     write(*,*)'PHYD  =',PHYD(1:3)
!     write(*,*)'THYD  =',THYD(1:3)
!     write(*,*)'LeftP =',LeftP_F(1:3)
!     write(*,*)'RightP=',RightP_F(1:3)
!     write(*,*)'RhoFlux =',RhoFlux_F(1:3)
!     write(*,*)'RhoUFlux=',RhoUFlux_F(1:3)
!     write(*,*)'eFlux   =',eFlux_F(1:3)
!     write(*,*)'LeftRho_F  =',LeftRho_F(1:3)
!     write(*,*)'RightRho_F =',RightRho_F(1:3)
!     write(*,*)'LeftU_F    =',LeftU_F(1:3)
!     write(*,*)'RightU_F   =',RightU_F(1:3)
!     write(*,*)'OldU_G     =',OldU_G(1:3)
!     write(*,*)'LeftP_F    =',LeftP_F(1:3)
!     write(*,*)'RightP_F   =',RightP_F(1:3)
!     write(*,*)'NewRho=',NewRho_C(1:3)
!     write(*,*)'NewP  =',NewP_C(1:3)
!     write(*,*)'NewT  =',NewT_C(1:3)
!
!  end if

end subroutine rusanov_solver

!==============================================================================

subroutine rusanov_flux( &
     LeftRho_F, LeftU_F, LeftP_F, &
     RightRho_F, RightU_F, RightP_F, &
     RhoFlux_F, RhoUFlux_F, eFlux_F)

  use ModCommonVariables
  implicit none


  real, dimension(MaxGrid+1), intent(in) :: &
       LeftRho_F, LeftU_F, LeftP_F, &
       RightRho_F, RightU_F, RightP_F
  real, dimension(MaxGrid+1), intent(out) :: RhoFlux_F, RhoUFlux_F, eFlux_F
  real    :: GammaOverGammaMinus1, InvGammaMinus1, RightE, LeftE, Coeff
  integer :: i
  !----------------------------------------------------------------------------
  InvGammaMinus1       = 1.0/(Gamma - 1.0)
  GammaOverGammaMinus1 = Gamma/(Gamma - 1.0)

  Coeff = 0.47/dtr1
!  Coeff = 0.5/dtr1
  
  do i=1,nDim+1
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

subroutine update_state( DoTest, DoTest2, DoTest3,&
     Rgas,&
     OldRho_C, OldU_C, OldP_C, OldT_C, &
     LeftRho_F, RightRho_F,LeftU_F, RightU_F,LeftP_F, RightP_F, &
     RhoFlux_F, RhoUFlux_F, eFlux_F, &
     RhoSource_C, RhoUSource_C, eSource_C, &
     NewRho_C, NewU_C, NewP_C, NewT_C)

  use ModCommonVariables
  implicit none
  logical, intent(in)                   :: DoTest,DoTest2,DoTest3  
  real, intent(in)                      :: Rgas 
  real, dimension(MaxGrid), intent(in)  :: OldRho_C, OldU_C, OldP_C, OldT_C
  real, dimension(MaxGrid+1), intent(in):: LeftRho_F, LeftU_F, LeftP_F
  real, dimension(MaxGrid+1), intent(in):: RightRho_F, RightU_F, RightP_F
  real, dimension(MaxGrid+1), intent(in):: RhoFlux_F, RhoUFlux_F, eFlux_F
  real, dimension(MaxGrid), intent(in)  :: RhoSource_C, RhoUSource_C, eSource_C
  real, dimension(MaxGrid), intent(out) :: NewRho_C, NewU_C, NewP_C,NewT_C
  real, dimension(MaxGrid) :: NewRhoU, NewE
  real :: InvGammaMinus1
  real, dimension(MaxGrid) :: NewRhoUSource_C, ExplicitRhoU_C
  real :: CollisionNeutral1, CollisionNeutral2, &
          CollisionNeutral3,Collisionneutral4
  real :: CollisionCoefT1, CollisionCoefM1, CollisionCoefT2, CollisionCoefM2,&
          CollisionCoefT3, CollisionCoefM3, CollisionCoefT4, CollisionCoefM4
  real :: ModifyESource
  integer :: i
  !----------------------------------------------------------------------------
  InvGammaMinus1 = 1.0/(Gamma - 1.0)

  do i=1,nDim
 !    ImpliciteSource_C(i)=0.0
     NewRho_C(i) = OldRho_C(i) &
          - Dt * (Ar23(i)*RhoFlux_F(i+1)-Ar12(i)*RhoFlux_F(i)) &
          / CellVolume_C(i) &
          + Dt*RhoSource_C(i)
     NewRhoU(i) = OldRho_C(i)*OldU_C(i) &
          - Dt * (Ar23(i)*RhoUFlux_F(i+1)-Ar12(i)*RhoUFlux_F(i)) &
          / CellVolume_C(i) &
          + Dt*Darea(i)*OldP_C(i) &
          + Dt*RhoUSource_C(i)

     NewE(i) = 0.5*OldRho_C(i)*OldU_C(i)**2 + InvGammaMinus1*OldP_C(i) &
          - Dt * (Ar23(i)*eFlux_F(i+1)-Ar12(i)*eFlux_F(i)) &
          / CellVolume_C(i) &
          + Dt*eSource_C(i) 
 
     NewU_C(i) = NewRhoU(i)/NewRho_C(i)

     NewP_C(i) = gMin1*(NewE(i) - 0.5 * NewU_C(i)**2 * NewRho_C(i))

     If (IsImplicit .and. (NamePlanet .eq. 'Saturn')) then
!        !get rid of the contribution to the terms that are implicit
        If (DoTest) then
           CollisionNeutral1=CFHpH(i)
           CollisionNeutral2=CFHPH2(i)
           CollisionCoefT1  =CTHpH
           CollisionCoefM1  =XMSH/(XMSH+XMSH)
           CollisionCoefT2  =CTHpH2
           CollisionCoefM2  =XMSH/(XMSH+2.*XMSH)

           CollisionNeutral3=CFHpH3p(i)
           CollisionCoefT3  =CTHpH3p
           CollisionCoefM3  =XMSH/(XMSH+XMSO)

           CollisionNeutral4=CFHpEL(i)
           CollisionCoefT4  =CTHpEL
           CollisionCoefM4  =XMSH/(XMSH+XMSE)

        else
           CollisionNeutral1=CFH3pH(i)
           CollisionNeutral2=CFH3pH2(i)
           CollisionCoefT1  =CTH3pH
           CollisionCoefM1  =XMSO/(XMSO+XMSH)
           CollisionCoefT2  =CTH3pH2
           CollisionCoefM2  =XMSO/(XMSO+2.*XMSH)

           CollisionNeutral3=CFH3pHp(i)
           CollisionCoefT3  =CTH3pHp
           CollisionCoefM3  =XMSO/(XMSO+XMSH)

           CollisionNeutral4=CFH3pEL(i)
           CollisionCoefT4  =CTH3pEL
           CollisionCoefM4  =XMSO/(XMSO+XMSE)

        endif
        !eliminate contribution of implicit sources in RhoU
        ExplicitRhoU_C(i) = NewRhoU(i)+Dt*OldU_C(i)*OldRho_C(i)* &
             ( CollisionNeutral1 + CollisionNeutral2 &
             + CollisionNeutral3 + CollisionNeutral4)
        
        NewRhoU(i)= ExplicitRhoU_C(i) / &
             (1+Dt*(CollisionNeutral1 + CollisionNeutral2 &
             +      CollisionNeutral3 + CollisionNeutral4))
        
        NewU_C(i) = NewRhoU(i)/NewRho_C(i)
        !When making p implicit we have to take the cf*cm*Ti term out of 
        !esource. Do this by adding ModifyEsource*dt to esource.
        ModifyESource = OldRho_C(i)* &
             (CollisionNeutral1*CollisionCoefT1  &
             +CollisionNeutral2*CollisionCoefT2  &
             +CollisionNeutral3*CollisionCoefT3  &
             +CollisionNeutral4*CollisionCoefT4)*OldT_C(i)
!
!        if (i .lt. 3 .and. dotest)write(*,*) ModifyESource, &
!             -tHyd(i)*dHyd(i)*(CTHpH*CFHpH(I)+CTHpH2*CFHpH2(I))
!        
!        if (i .lt. 4)write(*,*) i,newe(i)
        NewE(i)=NewE(i) + Dt*ModifyESource
!        if (i .lt. 4)write(*,*) i,newe(i)

        NewP_C(i) = gMin1*(NewE(i) - 0.5 * NewU_C(i)**2 * NewRho_C(i))
!        if (i .lt. 3) write(*,*) i,(1+Dt*3.* &
!             (CollisionNeutral1*CollisionCoefM1 + &
!             CollisionNeutral1*CollisionCoefM1)),NewP_C(i)
        
        NewP_C(i) = NewP_C(i) / (1+Dt*2.* &
             (CollisionNeutral1*CollisionCoefM1  &
             +CollisionNeutral2*CollisionCoefM2  & 
             +CollisionNeutral3*CollisionCoefM3  &
             +CollisionNeutral4*CollisionCoefM4))

!        if (i .lt. 4)write(*,*) NewP_C(i)/gMin1 +.5*NewU_C(i)**2*NewRho_C(i)
     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     If (IsImplicit .and. (NamePlanet .eq. 'Earth ')) then
!        !get rid of the contribution to the terms that are implicit
        If (DoTest .and. .not. DoTest3) then
           CollisionNeutral1=CFHN2(i)
           CollisionNeutral2=CFHO(i)
           CollisionCoefT1  =CTHN2
           CollisionCoefM1  = XMSH/(XMSH+28.0*XAMU)
           CollisionCoefT2  =CTHO
           CollisionCoefM2  =XMSH/(XMSH+XMSO)

           CollisionNeutral3=CFHO2(i)
           CollisionCoefT3  =CTHO2
           CollisionCoefM3  =XMSH/(XMSH+2.0*XMSO)
                             
           CollisionNeutral4=CFHHL(i)
           CollisionCoefT4  =CTHHL
           CollisionCoefM4  =XMSH/(XMSH+XMSHE)

        else if (DoTest2) then
           CollisionNeutral1=CFOXN2(i)
           CollisionNeutral2=CFOXO(i)
           CollisionCoefT1  =CTOXN2
           CollisionCoefM1  =XMSO/(XMSO+28.0*XAMU)
           CollisionCoefT2  =CTOXO
           CollisionCoefM2  =XMSO/(XMSO+XMSO)

           CollisionNeutral3=CFOXO2(i)
           CollisionCoefT3  =CTOXO2
           CollisionCoefM3  =xmso/(xmso+2.0*xmso)
                             
           CollisionNeutral4=CFOXHL(i)
           CollisionCoefT4  =CTOXHL
           CollisionCoefM4  =xmso/(xmso+xmshe)

        else
           CollisionNeutral1=CFHeN2(i)
           CollisionNeutral2=CFHeO(i)
           CollisionCoefT1  =CTHeN2
           CollisionCoefM1  =XMSHe/(XMSHe+28.0*XAMU)
           CollisionCoefT2  =CTHeO
           CollisionCoefM2  =XMSHe/(XMSHe+XMSO)

           CollisionNeutral3=CFHeO2(i)
           CollisionCoefT3  =CTHeO2
           CollisionCoefM3  =xmsHe/(xmsHe+2.0*xmso)
                             
           CollisionNeutral4=CFHeHE(i)
           CollisionCoefT4  =CTHeHE
           CollisionCoefM4  =xmsHe/(xmsHe+xmshe)

        endif
        !eliminate contribution of implicit sources in RhoU
        ExplicitRhoU_C(i) = NewRhoU(i)+Dt*OldU_C(i)*OldRho_C(i)* &
             ( CollisionNeutral1 + CollisionNeutral2 &
             + CollisionNeutral3 + CollisionNeutral4)
        
        NewRhoU(i)= ExplicitRhoU_C(i) / &
             (1+Dt*(CollisionNeutral1 + CollisionNeutral2 &
             +      CollisionNeutral3 + CollisionNeutral4))
        
        NewU_C(i) = NewRhoU(i)/NewRho_C(i)
        !When making p implicit we have to take the cf*cm*Ti term out of 
        !esource. Do this by adding ModifyEsource*dt to esource.
        ModifyESource = OldRho_C(i)* &
             (CollisionNeutral1*CollisionCoefT1  &
             +CollisionNeutral2*CollisionCoefT2  &
             +CollisionNeutral3*CollisionCoefT3  &
             +CollisionNeutral4*CollisionCoefT4)*OldT_C(i)
!
!        if (i .lt. 3 .and. dotest)write(*,*) ModifyESource, &
!             -tHyd(i)*dHyd(i)*(CTHpH*CFHpH(I)+CTHpH2*CFHpH2(I))
!        
!        if (i .lt. 4)write(*,*) i,newe(i)
        NewE(i)=NewE(i) + Dt*ModifyESource
!        if (i .lt. 4)write(*,*) i,newe(i)

        NewP_C(i) = gMin1*(NewE(i) - 0.5 * NewU_C(i)**2 * NewRho_C(i))
!        if (i .lt. 3) write(*,*) i,(1+Dt*3.* &
!             (CollisionNeutral1*CollisionCoefM1 + &
!             CollisionNeutral1*CollisionCoefM1)),NewP_C(i)
        
        NewP_C(i) = NewP_C(i) / (1+Dt*2.* &
             (CollisionNeutral1*CollisionCoefM1  &
             +CollisionNeutral2*CollisionCoefM2  & 
             +CollisionNeutral3*CollisionCoefM3  &
             +CollisionNeutral4*CollisionCoefM4))

!        if (i .lt. 4)write(*,*) NewP_C(i)/gMin1 +.5*NewU_C(i)**2*NewRho_C(i)
     endif

 
     
     NewT_C(i)=NewP_C(i)/Rgas/NewRho_C(i)
  end do

end subroutine update_state
