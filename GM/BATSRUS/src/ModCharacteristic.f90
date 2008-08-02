!^CFG COPYRIGHT UM
module ModCharacteristicMhd

  use ModCoordTransform, ONLY: cross_product
  use ModVarIndexes
  use ModPhysics,ONLY:g,inv_gm1
  use ModNumConst

  implicit none

  private
  public:: get_dissipation_flux_mhd

  ! Named MHD wave indexes
  integer, parameter :: EntropyW_=Rho_, AlfvenRW_=Ux_, AlfvenLW_=Uy_, &
       SlowRW_=Uz_, FastRW_=Bx_, SlowLW_=By_, FastLW_=Bz_, DivBW_=nVar

  real, parameter :: cTolerance2=cTolerance**2

contains
  !============================================================================
  subroutine generate_tangent12(Normal_D,Tangent1_D,Tangent2_D)

    real,dimension(3),intent(in) ::Normal_D
    real,dimension(3),intent(out)::Tangent1_D,Tangent2_D
    !------------------------------------------------------------------------
    integer,dimension(1)::iMinAbs
    iMinAbs=minloc(abs(Normal_D))

    !Construct the vector along one of the coordinate axis which is 
    !farthest from the direction of Normal_D

    Tangent2_D=cZero;Tangent2_D(iMinAbs(1))=cOne 
    Tangent1_D=-cross_product(Normal_D,Tangent2_D)
    Tangent1_D=Tangent1_D/sqrt(sum(Tangent1_D**2))
    Tangent2_D=cross_product(Normal_D,Tangent1_D)

  end subroutine generate_tangent12

  !============================================================================

  function primitive_from_pseudochar(PseudoChar_V,RhoInv,XH)

    !Variable transformtaions: 
    !make sense only for the increments, not for variables

    real,dimension(nVar)::primitive_from_pseudochar
    real,dimension(nVar),intent(in)::PseudoChar_V
    real,intent(in)::RhoInv,XH
    !-------------------------------------------------------------------------
    primitive_from_pseudochar=PseudoChar_V
    primitive_from_pseudochar(Ux_:Uz_)=&
         primitive_from_pseudochar(Ux_:Uz_)*RhoInv
    primitive_from_pseudochar(P_)=&
         primitive_from_pseudochar(P_)-&
         primitive_from_pseudochar(Rho_)*XH

  end function primitive_from_pseudochar
  !============================================================================
  function flux_from_pseudochar(PseudoChar_V,U_D,B_D,XH)

    !Again: "flux" means the flux difference recovered from the jumps in
    !pseudocharacterictic variables

    real,dimension(nVar+1)::flux_from_pseudochar
    real,dimension(nVar)::PseudoChar_V(nVar)
    real,dimension(3),intent(in)::U_D,B_D
    real,intent(in)::XH
    flux_from_pseudochar(1:nVar)=PseudoChar_V
    flux_from_pseudochar(RhoUx_:RhoUz_)=&
         PseudoChar_V(RhoUx_:RhoUz_)+&
         PseudoChar_V(Rho_)*U_D
    flux_from_pseudochar(nVar+1)=PseudoChar_V(P_)*inv_gm1+&
         sum(PseudoChar_V(RhoUx_:RhoUz_)*U_D)+&
         sum(PseudoChar_V(Bx_:Bz_)*B_D)+&
         PseudoChar_V(Rho_)*&
         (cHalf*sum(U_D**2)+(cOne-inv_gm1)*XH)
    flux_from_pseudochar(P_)=PseudoChar_V(P_) - XH*PseudoChar_V(Rho_)

  end function flux_from_pseudochar

  !===========================================================================

  subroutine decompose_state(&
       Normal_D,            &
       StateL_V,StateR_V,   &
       B0_D,                &
       Eigenvector_VV,      &
       DeltaWave_V,         &   !Wave amplitudes, dimensionless
       RhoH,XH,UH_D,B1H_D,  &   !Transformation coefficients
       Eigenvalue_V,        & 
       EigenvalueL_V,       &
       EigenvalueR_V)

    real, intent(in):: Normal_D(3)
    real,intent(in) :: StateL_V(nVar), StateR_V(nVar) !Primitive vars:Rho,u,B,P
    real,intent(in) :: B0_D(3)
    real,intent(out):: Eigenvector_VV(nVar,nVar) !Rho, RhoU, B, P+XH*Rho
    real,intent(out):: DeltaWave_V(nVar)         !Dimensionless
    real,intent(out):: RhoH, UH_D(3), B1H_D(3), XH
    real,optional,intent(out),dimension(nVar-1):: &
         Eigenvalue_V, EigenvalueL_V, EigenvalueR_V

    ! Jump in the state
    real, dimension(nVar) :: dState_V
    ! Jumps in variables
    real::dRho,dBn,dBt1,dBt2,dUnRho,dUt1Rho,dUt2Rho,dP

    !Vectors of the coordinate system associated with the averaged nagnetic
    !field: Tangent1_D is the vector of a unity length directed towards the 
    !tangential part of the averaged magnetic field, 
    !Tangent2_D=Normal_D x Tangent1_D
    real,dimension(3):: Tangent1_D, Tangent2_D

    ! B0Face
    real:: B0n, B0Tangent_D(3)

    ! Left face
    real :: RhoL,UL_D(3),BL_D(3)
    real :: BnL,aL,CsL,CaL,CfL

    ! Right face
    real :: RhoR,UR_D(3),BR_D(3)
    real :: BnR,aR,CsR,CaR,CfR

    ! Average (hat)
    real :: BH_D(3),XnH
    real :: BnH,aH,CsH,CaH,CfH

    !Used in averaging
    real :: WeightInv , ScalarAvr
    integer::iScalar

    !Normalization coefficient for left eigenvectors
    real:: NormCoef

    !More face variables

    real :: AlphaS, AlphaF, SignBnH

    real :: RhoInvL , RhoInvR , RhoInvH
    real :: RhoSqrtH, RhoSqrtL, RhoSqrtR

    real ::BTang2 !Reusable tangential magnetic field squared

    real ::Tmp    
    !--------------------------------------------------------------------------

    dState_V=StateL_V-StateR_V
    ! Scalar variables
    RhoL    =  StateL_V(rho_)
    RhoInvL =  cOne/RhoL
    aL      =  g * StateL_V(P_ ) * RhoInvL

    RhoR    =  StateR_V(rho_)
    RhoInvR =  cOne/RhoR
    aR      =  g * StateR_V(P_ ) * RhoInvR

    !\
    ! Set some values that are reused over and over
    !/

    RhoSqrtL   =sqrt(RhoL)
    RhoSqrtR   =sqrt(RhoR)


    !Srart averaging
    RhoH       = RhoSqrtL * RhoSqrtR
    !\
    ! Set some values that are reused over and over
    !/
    RhoInvH     = cOne/RhoH
    RhoSqrtH    = sqrt(RhoH)    

    WeightInv=cOne/(RhoSqrtL + RhoSqrtR)
    !Average velocity:

    UL_D=StateL_V(Ux_:Uz_)
    UR_D=StateR_V(Ux_:Uz_)

    UH_D=WeightInv * (RhoSqrtL * UL_D + RhoSqrtR * UR_D)


    BL_D=StateL_V(Bx_:Bz_)
    BR_D=StateR_V(Bx_:Bz_)

    BnL         =sum(Normal_D*BL_D)
    BnR         =sum(Normal_D*BR_D)

    !Leave only tangential components in the magnetic field vectors:

    BL_D=BL_D-BnL*Normal_D; BR_D=BR_D-BnR*Normal_D

    !Average them with the density dependent weights and add 
    !the arithmetic average of the normal components:
    B1H_D=WeightInv * (RhoSqrtL * BR_D + RhoSqrtR * BL_D)&
         + cHalf * (BnL+BnR) * Normal_D

    XnH = 0.25*RhoInvH*(BnR-BnL)**2 
    XH  =cHalf*WeightInv**2*sum((BL_D-BR_D)**2)

    !Average the speed of sound
    aH   =WeightInv * ( RhoSqrtL* aL +  RhoSqrtR* aR) +&
         g * XH + (g-cOne) * (XnH + cHalf * sum(dState_V(Ux_:Uz_)**2)*&
         RhoH  * WeightInv**2)

    !Below B1H is used only in the transformation matrix for the 
    !conservative variables. Add B0 field 

    BH_D=B1H_D+B0_D
    BnH=sum(Normal_D*BH_D)

    !Leave only tangential components in the magnetic field 
    BH_D=BH_D-BnH*Normal_D


    call get_characteristic_speeds(aH,RhoInvH,RhoSqrtH,BnH,BH_D,CsH,CaH,CfH)

    ! The components of eigenvectors for fast- and slow- sounds depend 
    ! on sgn(BnH) and Alphas 
    SignBnH     = sign(cOne,BnH)
    Tmp = CfH**2 - CsH**2
    if (Tmp > cTolerance2) then
       AlphaF = max(cZero,(aH**2  - CsH**2)/Tmp)
       AlphaS = max(cZero,(CfH**2 - aH**2 )/Tmp)

       AlphaF = sqrt(AlphaF)
       AlphaS = sqrt(AlphaS)
    else if (BnH**2 * RhoInvH <= aH**2 ) then
       AlphaF = cOne
       AlphaS = cZero
    else
       AlphaF = cZero
       AlphaS = cOne
    endif

    !Set the vectors of direction
    if(BTang2<cTolerance2)then
       !The direction of Tangent1_D vector is set more or less arbitralily
       call generate_tangent12(Normal_D,Tangent1_D,Tangent2_D)
    else
       !In non-degenerated case Tangent1_D is along the 
       !averaged magnetic field
       Tangent1_D=BH_D/sqrt(BTang2)
       Tangent2_D=cross_product(Normal_D,Tangent1_D)
    end if

    !Calculate jumps
    dRho=dState_V(rho_)
    dUnRho =sum(Normal_D  *dState_V(Ux_:Uz_))*RhoH
    dUt1Rho=sum(Tangent1_D*dState_V(Ux_:Uz_))*RhoH
    dUt2Rho=sum(Tangent2_D*dState_V(Ux_:Uz_))*RhoH
    dBn =sum(Normal_D  *dState_V(Bx_:Bz_))
    dBt1=sum(Tangent1_D*dState_V(Bx_:Bz_))
    dBt2=sum(Tangent2_D*dState_V(Bx_:Bz_))
    dP  =dState_V(p_)+XH*dRho

    !Calculate wave amplitudes and eigenvectors
    Eigenvector_VV=cZero; DeltaWave_V=cZero ; NormCoef=cHalf/(RhoH*aH*aH)
    !---------------------------------------------------------------------! 
    DeltaWave_V(EntropyW_)     = dRho*RhoInvH-dP*2*NormCoef

    Eigenvector_VV(rho_,EntropyW_)= RhoH
    !---------------------------------------------------------------------!

    DeltaWave_V(AlfvenRW_)  = NormCoef*(-aH*dUt2Rho + SignBnH*RhoSqrtH*aH*dBt2)


    Eigenvector_VV(RhoUx_:RhoUz_,AlfvenRW_)=            -RhoH*aH*Tangent2_D
    Eigenvector_VV(Bx_:Bz_,      AlfvenRW_)= SignBnH*RhoSqrtH*aH*Tangent2_D
    !---------------------------------------------------------------------!
    DeltaWave_V(AlfvenLW_)  = NormCoef*( aH*dUt2Rho + SignBnH*RhoSqrtH*aH*dBt2)

    Eigenvector_VV(RhoUx_:RhoUz_,AlfvenLW_)=             RhoH*aH*Tangent2_D
    Eigenvector_VV(Bx_:Bz_,      AlfvenLW_)= SignBnH*RhoSqrtH*aH*Tangent2_D
    !---------------------------------------------------------------------!
    DeltaWave_V(SlowRW_)    = NormCoef*(AlphaS*( dP        +CsH*dUnRho) + &
         AlphaF*( CfH*dUt1Rho*SignBnH &
         -RhoSqrtH*aH*dBt1))

    Eigenvector_VV(rho_         ,SlowRW_)=  RhoH*AlphaS
    Eigenvector_VV(RhoUx_:RhoUz_,SlowRW_)=  RhoH*(AlphaS*CsH*Normal_D  &
         + AlphaF*CfH*SignBnH*Tangent1_D)
    Eigenvector_VV(Bx_:Bz_      ,SlowRW_)= -RhoSqrtH*AlphaF*aH*Tangent1_D
    Eigenvector_VV(P_           ,SlowRW_)=  RhoH*AlphaS*aH**2
    !------------------------------------------------------------------------! 
    DeltaWave_V(FastRW_)    = NormCoef*(AlphaF*( dP        +CfH*dUnRho) + &
         AlphaS*(-CsH*dUt1Rho*SignBnH &
         +RhoSqrtH*aH*dBt1))

    Eigenvector_VV(rho_         ,FastRW_)=  RhoH*AlphaF
    Eigenvector_VV(RhoUx_:RhoUz_,FastRW_)=  RhoH*(AlphaF*CfH*Normal_D  &
         - AlphaS*CsH*SignBnH*Tangent1_D)
    Eigenvector_VV(Bx_:Bz_      ,FastRW_)=  RhoSqrtH*AlphaS*aH*Tangent1_D
    Eigenvector_VV(P_           ,FastRW_)=  RhoH*AlphaF*aH**2
    !------------------------------------------------------------------------!
    DeltaWave_V(SlowLW_)    = NormCoef*(AlphaS*( dP        -CsH*dUnRho) + &
         AlphaF*(-CfH*dUt1Rho*SignBnH &
         -RhoSqrtH*aH*dBt1))

    Eigenvector_VV(rho_         ,SlowLW_)=  RhoH*AlphaS
    Eigenvector_VV(RhoUx_:RhoUz_,SlowLW_)=  RhoH*(-AlphaS*CsH*Normal_D  &
         - AlphaF*CfH*SignBnH*Tangent1_D)
    Eigenvector_VV(Bx_:Bz_      ,SlowLW_)= -RhoSqrtH*AlphaF*aH*Tangent1_D
    Eigenvector_VV(P_           ,SlowLW_)=  RhoH*AlphaS*aH**2 

    !------------------------------------------------------------------------! 
    DeltaWave_V(FastLW_)    = NormCoef*(AlphaF*( dP        -CfH*dUnRho) + &
         AlphaS*( CsH*dUt1Rho*SignBnH &
         +RhoSqrtH*aH*dBt1))

    Eigenvector_VV(rho_         ,FastLW_)=  RhoH*AlphaF
    Eigenvector_VV(RhoUx_:RhoUz_,FastLW_)=  RhoH*(-AlphaF*CfH*Normal_D  &
         + AlphaS*CsH*SignBnH*Tangent1_D)
    Eigenvector_VV(Bx_:Bz_      ,FastLW_)=  RhoSqrtH*AlphaS*aH*Tangent1_D
    Eigenvector_VV(P_           ,FastLW_)=  RhoH*AlphaF*aH**2 
    !------------------------------------------------------------------------!

    DeltaWave_V(DivBW_)     = dBn

    Eigenvector_VV(Bx_:Bz_,DivBW_)      = Normal_D
    !------------------------------------------------------------------------!
    do iScalar = ScalarFirst_, ScalarLast_
       ScalarAvr=WeightInv*(StateL_V(iScalar)*RhoSqrtR+&
            StateR_V(iScalar)*RhoSqrtL)
       Tmp=ScalarAvr*RhoInvH

       DeltaWave_V(iScalar)=dState_V(iScalar)-Tmp*dRho

       Eigenvector_VV(iScalar,iScalar)  = cOne 

       !All these waves advect the density, hence, the passive scalar

       Eigenvector_VV(iScalar,EntropyW_) = Eigenvector_VV(rho_,EntropyW_)*Tmp
       Eigenvector_VV(iScalar,SlowRW_ )  = Eigenvector_VV(rho_,SlowRW_  )*Tmp
       Eigenvector_VV(iScalar,FastRW_ )  = Eigenvector_VV(rho_,FastRW_  )*Tmp
       Eigenvector_VV(iScalar,SlowLW_ )  = Eigenvector_VV(rho_,SlowLW_  )*Tmp
       Eigenvector_VV(iScalar,FastLW_ )  = Eigenvector_VV(rho_,FastLW_  )*Tmp

    end do


    if(.not.present(Eigenvalue_V))return
    !Calculate eigenvalues:add the normal velocity first
    Eigenvalue_V =sum(UH_D*Normal_D)
    call set_eigenvalues(Eigenvalue_V,CsH,CaH,CfH)

    if(.not.present(EigenvalueL_V))return
    EigenvalueL_V=sum(UL_D*Normal_D)
    EigenvalueR_V=sum(UR_D*Normal_D)


    !For left and right stetes we need first to add the B0 field
    !Split it for normal and tangential components:

    B0n=sum(Normal_D*B0_D)
    B0Tangent_D = B0_D - B0n*Normal_D

    BnL = BnL + B0n; BL_D = BL_D + B0Tangent_D
    call get_characteristic_speeds(aL,RhoInvL,RhoSqrtL,&
         BnH,BL_D,CsL,CaL,CfL)

    !CaL=CaL*sign(cOne,BnL)*SignBnH
    call set_eigenvalues(EigenvalueL_V,CsL,CaL,CfL)

    BnR=BnR+B0n; BR_D = BR_D+B0Tangent_D
    call get_characteristic_speeds(aR,RhoInvR,RhoSqrtR,&
         BnH,BR_D,CsR,CaR,CfR)
    !CaR=CaR*sign(cOne,BnR)*SignBnH
    call set_eigenvalues(EigenvalueR_V,CsR,CaR,CfR)

  contains
    !-------------------------------------------------------------------------!
    subroutine set_eigenvalues(Value_V,Cs,Ca,Cf)
      real,dimension(nVar-1),intent(inout)::Value_V
      real,intent(in):: Cs, Ca, Cf
      !\
      ! Eigenvalues
      !/
      Value_V(AlfvenRW_) = Value_V(AlfvenRW_) + Ca
      Value_V(AlfvenLW_) = Value_V(AlfvenLW_) - Ca
      Value_V(SlowRW_)   = Value_V(SlowRW_  ) + Cs
      Value_V(FastRW_)   = Value_V(FastRW_  ) + Cf
      Value_V(SlowLW_)   = Value_V(SlowLW_  ) - Cs
      Value_V(FastLW_)   = Value_V(FastLW_  ) - Cf
    end subroutine set_eigenvalues
    !-------------------------------------------------------------------------!
    subroutine get_characteristic_speeds(A,         &
         RhoInv,    &
         RhoSqrt,   &
         Bn,        &
         BTang_D,   &
         Cs,        &
         Ca,        &
         Cf)

      real,intent(inout)::A !In: speed of sound squared, out speed of sound
      real,intent(in)   ::RhoInv,RhoSqrt,Bn,BTang_D(3)
      real,intent(out)  ::Cs,Ca,Cf
      !-----------------------------------------------------------------------
      A = sqrt(A) !Speed of sound

      !BTang2 is reused while constructing the coordinate system
      BTang2 = sum(BTang_D**2) 
      Tmp=BTang2*RhoInv

      Ca=abs(Bn)/RhoSqrt                                 !Alfven speed
      Cf=cHalf*(sqrt((A-Ca)**2+Tmp)+sqrt((A+Ca)**2+Tmp)) !Fast magnetossonic
      Cs=Ca * A/Cf                                       !Slow magnetosonic

    end subroutine get_characteristic_speeds

  end subroutine decompose_state

  !===========================================================================

  subroutine get_fixed_abs_eigenvalue(&
       Eigenvalue_V,&
       EigenvalueL_V,&
       EigenvalueR_V,&
       LambdaB0, &
       EigenvalueFixed_V,CMax,IsBoundary)

    use ModMain, ONLY: Climit

    real,intent(in) ,dimension(nVar-1)::Eigenvalue_V
    real,intent(in) ,dimension(nVar-1)  ::EigenvalueL_V,EigenvalueR_V
    real,intent(in) ::LambdaB0
    real,intent(out),dimension(nVar-1)  ::EigenvalueFixed_V
    real,intent(out)::cMax
    logical,intent(in)::IsBoundary
    integer::iWave
    real::Eps_V(nVar-1),Lambda
    !-------------------------------------------------------------------------
    Eps_V=max(LambdaB0,abs(Eigenvalue_V(FastRW_)-Eigenvalue_V(FastLW_))*0.05)
    do iWave=1,nVar-1
       Lambda=Eigenvalue_V(iWave)
       Eps_V(iWave)=max(EigenvalueR_V(iWave)-Lambda,&
            Lambda -EigenvalueL_V(iWave),Eps_V(iWave))
       EigenvalueFixed_V(iWave)=max(abs(Lambda),Eps_V(iWave))!+&
       !cHalf*min(Lambda**2/Eps_V(iWave)-Eps_V(iWave),cZero)
    end do

    cMax = max(EigenvalueFixed_V(FastRW_),EigenvalueFixed_V(FastLW_))

    if(Climit > 0.0) EigenvalueFixed_V = min(Climit, EigenvalueFixed_V)

  end subroutine get_fixed_abs_eigenvalue

  !===========================================================================!

  subroutine get_dissipation_flux_mhd(Dir_D, StateL_V, StateR_V, &
       B0_D, DeltaB0_D, uL_D, uR_D, DeltaBnL,DeltaBnR, IsBoundary, DoTest, &
       DissipationFlux_V, cMax, Un)

    real,    intent(in):: Dir_D(3)
    real,    intent(in):: StateL_V(nVar), StateR_V(nVar)
    real,    intent(in):: B0_D(3), DeltaB0_D(3), uL_D(3), uR_D(3)
    real,    intent(in):: DeltaBnL,DeltaBnR
    logical, intent(in):: IsBoundary, DoTest
    real,   intent(out):: DissipationFlux_V(nVar+1)
    real,   intent(out):: cMax, Un

    real,dimension(nVar,nVar)::Eigenvector_VV
    real,dimension(nVar)     ::DeltaWave_V
    real,dimension(nVar-1)   :: Eigenvalue_V
    real,dimension(nVar-1)     ::   EigenvalueL_V,  &
         EigenvalueR_V
    real :: RhoH,UH_D(3),B1H_D(3),XH,UnL,UnR,LambdaB0
    real,dimension(nVar)  ::EigenvalueFixed_V,FluxPseudoChar_V
    integer::iWave!,iDir_D(3)
    !--------------------------------------------------------------------------

    call decompose_state(Dir_D, StateL_V,StateR_V, B0_D, &
         Eigenvector_VV,        &
         DeltaWave_V,           &    !Wave amplitudes, dimensionless
         RhoH, XH, UH_D, B1H_D, &   !Transformation coefficients
         Eigenvalue_V, EigenvalueL_V, EigenvalueR_V)

    LambdaB0=sqrt((sum(DeltaB0_D**2))/RhoH)
    call get_fixed_abs_eigenvalue(&
         Eigenvalue_V,&
         EigenvalueL_V,&
         EigenvalueR_V,&
         LambdaB0,&
         EigenvalueFixed_V(1:nVar-1),&
         CMax,IsBoundary)

    UnL= sum(uL_D*Dir_D); UnR=sum(uR_D*Dir_D)
    EigenvalueFixed_V(DivBW_)=max(abs(UnL),abs(UnR))
    cMax=max(cMax,EigenvalueFixed_V(DivBW_))
    if(IsBoundary)EigenvalueFixed_V=cMax
    FluxPseudoChar_V=cZero

    do iWave=1,nVar-1
       FluxPseudoChar_V = FluxPseudoChar_V + &
            Eigenvector_VV(:,iWave)*EigenvalueFixed_V(iWave)*&
            DeltaWave_V(iWave)
    end do

    FluxPseudoChar_V = FluxPseudoChar_V+&
         Eigenvector_VV(:,DivBW_)*&
         (UnL*DeltaBnR+UnR*DeltaBnL+&
         EigenvalueFixed_V(DivBW_)*(&
         DeltaWave_V(DivBW_)+DeltaBnL-DeltaBnR))

    DissipationFlux_V = 0.5*&
         flux_from_pseudochar(FluxPseudoChar_V,UH_D,B1H_D,XH)

    Un = sum(UH_D*Dir_D)

  end subroutine get_dissipation_flux_mhd

end module ModCharacteristicMhd
