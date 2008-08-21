Module ModTurbulence
  use ModConst
  implicit none
  SAVE
  logical::DoInitSpectrum=.true.
  logical::UseTurbulentSpectrum=.true.
  logical::UseAdvectionWithAlfvenSpeed=.false.
  logical::DoOutputGamma=.true.
  integer::iXOutputGamma=600
  real::DispersionOutput
  real,allocatable::Gamma_I(:)
  real,allocatable::IPlus_IX(:,:),IMinus_IX(:,:),IC(:)
  real,allocatable::VA_I(:)           !Alfven speed

  real,allocatable::DispersionA_I(:)  !Rate of advection in k space, neglecting 
                                      !a spacial advection with the Alfven speed

  real,allocatable::DispersionPlus_I(:)  !Rate of advection in k space, for I_+ wave 
  real,allocatable::DispersionMinus_I(:)  !Rate of advection in k space, for I_- wave

  real,allocatable:: RhoCompression_I(:) ! This is the ratio of densities powered 3/2 


    !-----------------------------------------------------------------------------------!
    !          Grid in the momentum space                                               !
    !iP     0     1                         nP   nP+1                                   !
    !       |     |    ....                 |     |                                     !
    !P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))              !
    !             |    Grid in k-space      |     |                                     !
    !K/B         KMax                      KMin                                         !
    !ik    nP+1   nP                        1     0                                     !
    !-----------------------------------------------------------------------------------!


  real::Kmin     !cElectronCharge/(PInjection*exp(real(nP)*DeltaLnP)) 
  real::KMax     !cElectronCharge/(PInjection*exp(DeltaLnP))
  real::DeltaLnK !DeltaLnP
  real::cInjectionHere
  real,external,private:: momentum_to_energy,momentum_to_kinetic_energy
  real,external,private:: energy_to_momentum,kinetic_energy_to_momentum
  real,external,private:: energy_in

  real,allocatable,private::AK_II(:,:)
  real,allocatable,private::BK_II(:,:)

  real,allocatable,private::CFL_I(:)

  integer,allocatable::CorrectionMode(:)

  !the intensity of the back travelling wave in the initial condition
  real,parameter::Alpha=cOne/10.0
  real,parameter::Lambda0=4.0/10.0  ![AU]
 
contains
  subroutine assign_kolmogorov_spectrum(iXStart,iXLast,iKStart,iKLast,X_DI,B_I)
    integer,intent(in)::iXStart,iXLast,iKStart,iKLast
    real,intent(in)::X_DI(1:3,iXStart:iXLast),B_I(iXStart:iXLast)
    integer::iX,iK
    real::K,KR0
    !--------------------
    do iX=iXStart,iXLast
       !We use the formulae:
       !I_{+}(k)=(1-\Alpha)*IC/K^{5/3} and I_{-}=\Alpha*I_{+}/(1-\Alpha),
       !where IC=54*B^2/(7*\pi* \Lambda_0* K0^{1/3}*(r/1 [AU])) and \Lambda_0=0.4 [AU], and K0=|e|B/1[GeV],
       !in CGS units
       
       !SI units
       KR0=cElectronCharge*B_I(iX)*cLightSpeed/energy_in('GeV')
       IC(iX)=(216./7.)*(B_I(iX)**2)/(Lambda0*cMu)*&  ! /1 AU
            (KR0**(-(1.0/3)))/&                        ! *1 AU
            sqrt(sum(X_DI(:,iX)**2))
       do iK=iKStart,iKLast
          K=B_I(iX)*Kmin*exp(real(iK-1)*DeltaLnK)
          IPlus_IX( iK,iX)=(cOne-Alpha)*IC(iX)/(K**(5.0/3))
          Iminus_IX(iK,iX)=Alpha/(cOne-Alpha)*IPlus_IX(iK,iX) 
       end do
    end do
  end subroutine assign_kolmogorov_spectrum
  !=======================Initial spectrum of turbulence=========================================!
  !We recover the initial spectrum of turbulence from the spatial distribution of
  !the diffusion coefficient and its dependence on the particle energy.
  subroutine initSpectrum(nX,nP,X_DI,B_I,pInjection,DeltaLnP,iShock,&
       cInjection,AlfvenMach)
    integer,intent(in)::nX,nP
    real,dimension(1:3,1:nX),intent(in)::X_DI     !Coordinates of Lagrangian Meshes [m]
    real,dimension(1:nX),intent(in)::B_I          !Magnetic field intensity         [T]
    real,intent(in)::pInjection,DeltaLnP,cInjection,AlfvenMach
    integer, intent(in):: iShock

    integer::iX,iK
    real::K
    real::ICOld
    real::R,RSh

    !-----------------------------------------------------------------------------------!
    allocate(IPlus_IX(0:nP+1,1:nX),IMinus_IX(0:nP+1,1:nX),IC(1:nX),CorrectionMode(1:nX)) 
    IPlus_IX=cZero;IMinus_IX=cZero
    allocate(VA_I(nX)); VA_I=cZero
    allocate(DispersionA_I(1:nX)); DispersionA_I=cZero
    if(UseAdvectionWithAlfvenSpeed)then
       allocate(DispersionPlus_I(1:nX)); DispersionPlus_I=cZero
       allocate(DispersionMinus_I(1:nX));DispersionMinus_I=cZero
       allocate(CFL_I(1:nX));CFL_I=cOne
    end if
    allocate(RhoCompression_I(1:nX));RhoCompression_I=cZero
    if(DoOutputGamma)allocate(Gamma_I(nP))
    !-----------------------------------------------------------------------------------!

    !-----------------------------------------------------------------------------------!
    !          Grid in the momentum space                                               !
    !iP     0     1                         nP   nP+1                                   !
    !       |     |    ....                 |     |                                     !
    !P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))              !
    !             |    Grid in k-space      |     |                                     !
    !K/B         KMax                      KMin                                         !
    !ik    nP+1   nP                        1     0                                     !
    !-----------------------------------------------------------------------------------!


    KMin=cElectronCharge/(PInjection*exp(real(nP)*DeltaLnP))   
    KMax=cElectronCharge/(PInjection*exp(DeltaLnP))
    DeltaLnK=DeltaLnP
    cInjectionHere=cInjection
    call assign_kolmogorov_spectrum(0,nP+1,1,nX,X_DI,B_I)
    CorrectionMode=1.0
   
    RSh=sqrt(sum(X_DI(:,iShock)**2))
    do iX=1,nX
       R=sqrt(sum(X_DI(:,iX)**2))
       if (R<(cOne+cOne/10.0)*RSh.and.AlfvenMach>1.0) then
          ! In this part of the spectrum another equation governs the diffusion
          ICOld=IC(iX)
          IC(iX)=B_I(iX)**2*10.0*CInjection*AlfvenMach*&
               min(cOne,R/Rsh/(cOne-cOne/10.0))/(cMu*3.0)
          if(IC(iX)>ICOld)then
             do iK=0,nP+1
                K=B_I(iX)*Kmin*exp(real(iK-1)*DeltaLnK) 
                IPlus_IX( iK,iX)=(cOne-Alpha)*IC(iX)/K
                IMinus_IX(iK,iX)=Alpha/(cOne-Alpha)*IPlus_IX(iK,iX)
             end do
             CorrectionMode(iX)=2
          else
             !Do not change if the newly calculated value is less than the old one
             IC(iX)=ICOld
          end if
       end if
    end do
    
    DoInitSpectrum = .false.
   
  end subroutine initSpectrum

  !======================================================================!

  !=======================Set advection rate in k-space==================!
  subroutine set_wave_advection_rates(nX,B_I,BOld_I,Rho_I,RhoOld_I,X_DI,&
       DeltaLnP,Dt,DtReduction)
    integer,intent(in)::nX
    real,intent(in)::B_I(1:nX),BOld_I(1:nX),Rho_I(1:nX),RhoOld_I(1:nX)
    real,intent(in)::X_DI(1:3,1:nX)
    real,intent(in)::DeltaLnP, Dt
    real,intent(out)::DtReduction
    !------------------------------
    integer::iLoop
   
    !Calculate alfven speed
    VA_I(1:nX)=B_I(1:nX)/sqrt(cMu*cProtonMass*Rho_I(1:nX))
    
    !Calculate the wave increment and update the wave spectra
    RhoCompression_I(1:nX)=(3.0*cHalf)*log(Rho_I(1:nX)/RhoOld_I(1:nX))

    !Contribution to the advection in k_space from the Lagrangian derivatives:
    ! Dispersion = Dt * [ (D ln {\rho})/(Dt) - 2 (D ln { B })/(D t) ] + ...
 
    DispersionA_I(1:nX)=log(Rho_I(   1:nX)*BOld_I(1:nX)**2/&
                           (RhoOld_I(1:nX)*B_I(   1:nX)**2)) 
   
    if(UseAdvectionWithAlfvenSpeed)then

       CFL_I(1 )=Dt*VA_I(1 )/sqrt(sum( (X_DI(:, 2)-X_DI(:,   1))**2 ))
       CFL_I(nX)=Dt*VA_I(nX)/sqrt(sum( (X_DI(:,nX)-X_DI(:,nX-1))**2 ))

       do iLoop=2,nX-1
          CFL_I(iLoop)=Dt*VA_I(iLoop)/&
               (cHalf*sqrt(sum((X_DI(:,iLoop+1)-X_DI(:,iLoop-1))**2)))   
       end do

       !Now add the contribution from the spacial derivative:
       ! ... + Dt * [ \partial ((1/2) ln {rho})/(\partial s)-2\partial( ln{B})/(\partial s) ]  
       !for the "plus" wave 
       
       DispersionPlus_I(1)=DispersionA_I(1  )+CFL_I(1 )*&
            (cHalf*log(Rho_I(2)/Rho_I(1    ))-2.0*log(B_I(2 )/B_I(1   )))
       DispersionPlus_I(nX)=DispersionA_I(nX)+CFL_I(nX)*&
            (cHalf*log(Rho_I(nX)/Rho_I(nX-1))-2.0*log(B_I(nX)/B_I(nX-1)))
       
       do iLoop=2,nX-1
          DispersionPlus_I(iLoop)=DispersionA_I(iLoop)+CFL_I(iLoop)*&
               (0.250*log(Rho_I(iLoop+1)/Rho_I(iLoop-1))-&
                         log(B_I(  iLoop+1)/B_I(  iLoop-1)))
       end do
      
       !Now add the contribution from the spacial derivative:
       ! ... - Dt * [ \partial ((1/2) ln {rho})/(\partial s)-2\partial( ln{B})/(\partial s) ]  
       !for the "minus" wave 

       DispersionMinus_I(1:nX)=2.0*DispersionA_I(   1:nX)-&
                                    DispersionPlus_I(1:nX)
      
       DispersionPlus_I( 1:nX)=DispersionPlus_I( 1:nX)/DeltaLnP
       DispersionMinus_I(1:nX)=DispersionMinus_I(1:nX)/DeltaLnP
       
       DtReduction=max(maxval( CFL_I(1:nX)),&
            maxval(abs(DispersionPlus_I( 1:nX))),&
            maxval(abs(DispersionMinus_I(1:nX))))
    else
       DispersionA_I(1:nX)=DispersionA_I(1:nX)/DeltaLnP
   
       DtReduction = maxval(abs(DispersionA_I(1:nX)))
    end if
  end subroutine set_wave_advection_rates
  !=======================================================================!
  subroutine reduce_advection_rates(nStep)
    integer,intent(in)::nStep
    if(UseAdvectionWithAlfvenSpeed)then
       DispersionPlus_I=DispersionPlus_I/real(nStep)
       DispersionMinus_I=DispersionMinus_I/real(nStep)
       CFL_I=CFL_I/real(nStep)
    else
       DispersionA_I=DispersionA_I/real(nStep)
    end if
    RhoCompression_I=RhoCompression_I/real(nStep)
  end subroutine reduce_advection_rates
  
  !=====The coefficient of a spatial diffusion along the magnetic field===!
  !Given by a formula:

  !D_{xx}=4*v*B^2/(cMu)*(AK_I(k_{res})-k_{res}^2*BK_I(k_{res})),

  !where AK_I=\int_{k_{res}}^\infty{d(\log k)/(k^2*(I_{+}(k) + I_{-}(k) )) }
  !and   BK_I=\int_{k_{res}}^\infty{d(\log k)/(k^4*(I_{+}(k) + I_{-}(k) )) }

  real function Dxx(P,B,iX,NameParticle)
    !--------------------------------------------------------------------
    real,intent(in)   :: P,B        !Momentum, magnetic field intensity, both in SI
    integer,intent(in):: iX         !Node number
    character(LEN=*),intent(in)::NameParticle
    !--------------------------------------------------------------------
    real              :: KR         !The resonant wave number 
    integer           :: iK
 



    !-----------------------------------------------------------------------------------!
    !          Grid in the momentum space                                               !
    !iP     0     1                         nP   nP+1                                   !
    !       |     |    ....                 |     |                                     !
    !P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))              !
    !             |    Grid in k-space      |     |                                     !
    !K/B         KMax                      KMin                                         !
    !ik    nP+1   nP                        1     0                                     !
    !-----------------------------------------------------------------------------------!


    !Calculate the resonant wave number
    KR=cElectronCharge*B/P

    !The integer representative for the wave number to take the partial sums 
    iK=1+nint(log(KR/(B*Kmin)) / DeltaLnK)  
    !Remember that  Kmin=cElectronCharge/(PInjection*exp(real(nP)*DeltaLnP)) 

    !Calculate D_{xx}: KRes-dependent part
    Dxx = (4.0*B**2*P*cLightSpeed**2)/&
         (momentum_to_energy(P,NameParticle)*cMu)*&
         (AK_II(iK,iX)-BK_II(iK,iX)*(KR**2))

  end function Dxx
  !=======================================================================!
  !=Calculate the partial sums as described above: 
  !    \sum_{k_{res}}^\infty d(log K)/(k^2*(I_+ + I_-))  (This is AK_I) 
  !and \sum_{k_{res}}^\infty d(log K)/(k^4*(I_+ + I_-))  (This is BK_I)
  subroutine set_dxx(nX,nP,B_I)
    integer,intent(in)::nX,nP
    real,intent(in),dimension(nX)::B_I  !The magnetic field intensity
    !-----------------------------------------------------------------!
    integer:: iX,iK
    real:: F01,F02,F11,F12
    real:: K0,K1,expdLogK
    real::SpectralIndexAtKMax

    if (.not.allocated(AK_II)) allocate(AK_II(nP,nX),BK_II(nP,nX)) 

    do iX=1,nX
       select case(CorrectionMode(iX))
       case(1)
          SpectralIndexAtKMax=cOne+2.0*(1.0/3)
       case(2)
          SpectralIndexAtKMax=cOne
       end select
       



    !-----------------------------------------------------------------------------------!
    !          Grid in the momentum space                                               !
    !iP     0     1                         nP   nP+1                                   !
    !       |     |    ....                 |     |                                     !
    !P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))              !
    !             |    Grid in k-space      |     |                                     !
    !K/B         KMax                      KMin                                         !
    !ik    nP+1   nP                        1     0                                     !
    !-----------------------------------------------------------------------------------!



       K0=B_I(iX)*Kmax

       !The integrand for AK_I
       F01=cOne/(K0**2)/(IPlus_IX(nP,iX)+IMinus_IX(nP,iX))

       !The integrand for BK_I 
       F02=cOne/(K0**4)/(IPlus_IX(nP,iX)+IMinus_IX(nP,iX))

       ExpDLogK=exp(DeltaLnK)

       !As the starting values for AK_I and BK_I at the maximal momentum,
       !solve the integrals from K_{max} up to \infty, assuming the
       !power law spectrum of turbulence at K>K_{max}
       AK_II(nP,iX)=F01/(2.0-SpectralIndexAtKMax)
       BK_II(nP,iX)=F02/(4.0-SpectralIndexAtKMax)

       do iK=nP-1,1,-1
          !We calculate the partial sums for a set of the wave number values.
          !The integral is taken from KRes up to infinity, so we start from 
          !the maximal wave number and add the contributions from each of 
          !the wave number intervals.

          !The lower value of the wave number
          K1=K0/expdLogK 
          
          !The integrands at the lower value of the wave number
          F11=cOne/(K1**2)/(IPlus_IX(iK,iX)+IMinus_IX(iK,iX)) 
          F12=cOne/(K1**4)/(IPlus_IX(iK,iX)+IMinus_IX(iK,iX))

          !Calculate the new partial sums

          AK_II(iK,iX)=AK_II(iK+1,iX)+cHalf*(F01+F11)*DeltaLnK   
          BK_II(iK,iX)=BK_II(iK+1,iX)+cHalf*(F02+F12)*DeltaLnK 

          ! current values saved as the initial values for the next step in the loop
          K0=K1; F01=F11; F02=F12
       end do
    end do
  end subroutine set_dxx
  
  !=====================Update the Turbulence spectrum===============!
  !This is the subroutine which advances the wave spectrum through
  !a time step, by solving the equations
  !dI_+/dt=\gamma_+ I_+
  !dI_-/dt=\gamma_- I_-
  !==================================================================!
  
  subroutine UpdateSpectrum(nX,nP,pInjection,DeltaLnP,&
       X_DI,F_II,B_I,Rho_I,SP_Dt)
    use ModLogAdvection
    use ModLinAdvection
    !-----------------------------------------------------------------
    
    integer,intent(in)::nX,nP                    !The number of points and the number of the energy intervals
    real,intent(in)::pInjection,DeltaLnP,SP_Dt  
    real,dimension(3,nX),intent(in)::X_DI         !Coordinates of the Lagrangian points
    real,dimension(0:nP+1,1:nX),intent(in)::F_II  !The distribution function
    real,dimension(nX),intent(in)::B_I,Rho_I      !The magnetic field intensity and the particle 
                                                  !number density, [T], [m^{-3}]

    integer::iX,iK,iP
    real:: PRes       !Resonant value of the particle momentum, for a given K

    real:: Gamma      !The increment for the I_+ wave, multiplied by (I_+  +  I_-). Below 
                      !the dynamical equations for I_+ and I_- are reformulated in terms of \gamma
                      !and integrated exactly assuming the constant value of \gamma (NOT the increment)
                      !through the time step 

    real:: K          !Wave number
    real:: DsPlus,DsMinus      !Forward and backward values for Ds 
    real:: DfDs0,P0,DfDs1,P1   !Spatial derivatives of the distribution function at given values of the 
                               !momentum
    real:: A(nP+1),B(nP+1)     !Partial sums in the integral for Gamma

    !Variables used in the spectral model for turbulence
    
    real:: ExpRhoCompression

    real::C1,C2,C3             !Miscellaneous 
    
    do iX=1,nX

       !Advection in k space:
       if(UseAdvectionWithAlfvenSpeed)then
          call advance_log_advection(DispersionPlus_I( iX),nP,1,1,&
                                     IPlus_IX( :,iX),IsConservative=.true.,&
                                     DeltaLnP=DeltaLnP)
          call advance_log_advection(DispersionMinus_I(iX),nP,1,1,&
                                     IMinus_IX(:,iX),IsConservative=.true.,&
                                     DeltaLnP=DeltaLnP)
       else
 
          call advance_log_advection(DispersionA_I(iX),nP,0,0,&
                                     IPlus_IX( 1:nP,iX),IsConservative=.true.,&
                                     DeltaLnP=DeltaLnP)
          call advance_log_advection(DispersionA_I(iX),nP,0,0,&
                                     IMinus_IX(1:nP,iX),IsConservative=.true.,&
                                     DeltaLnP=DeltaLnP)
       end if



    !-----------------------------------------------------------------------------------!
    !          Grid in the momentum space                                               !
    !iP     0     1                         nP   nP+1                                   !
    !       |     |    ....                 |     |                                     !
    !P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))              !
    !             |    Grid in k-space      |     |                                     !
    !K/B         KMax                      KMin                                         !
    !ik    nP+1   nP                        1     0                                     !
    !-----------------------------------------------------------------------------------!



       !Calculate the partial sums in the integral for \Gamma
       !Nullify the sums

       A(nP+1)=cZero; B(nP+1)=cZero

       !Calculate the spatial derivatives of the distribution function at the maximal energy

       if (iX==1) then

          DsPlus=max(sqrt(sum( (X_DI(:,iX+1)-X_DI(:,iX))**2 )),cTiny)
          !Use the forward spatial derivative

          DfDs0=(F_II(nP,iX+1)-F_II(nP,iX))/DsPlus

       else if (iX==nX) then

          DsMinus=max(sqrt(sum( (X_DI(:,iX)-X_DI(:,iX-1))**2 )),cTiny)
          !Use the backward spatial derivative

          DfDs0=(F_II(nP,iX)-F_II(nP,iX-1))/DsMinus

       else

          DsPlus =max(sqrt(sum( (X_DI(:,iX+1)-X_DI(:,iX  ))**2 )),cTiny)
          DsMinus=max(sqrt(sum( (X_DI(:,iX  )-X_DI(:,iX-1))**2 )),cTiny)
          !Use the average between the forward and backward spatial derivatives

          DfDs0=cHalf*((F_II(nP,iX+1)-F_II(nP,iX  ))/DsPlus+&
                       (F_II(nP,iX  )-F_II(nP,iX-1))/DsMinus)

       end if

       !The momentum at the maximal energy

       P0=PInjection*exp(real(nP+1)*DeltaLnP)  

       !We calculate the partial sums for a set of the momentum values.
       !The integral is taken from pRes up to infinity, so we start from the
       !maximal energy and add the contributions from each of the energy intervals.

       do iP=nP,1,-1

          !Calculate the momentum at the lower energy
          P1=PInjection*exp(real(iP)*DeltaLnP)  

          !Calculate the distribution function gradient at the lower energy

          if (iX==1) then
             DfDs1=(F_II(iP,iX+1)-F_II(iP,iX  ))/DsPlus
          else if (ix==nX) then
             DfDs1=(F_II(iP,iX  )-F_II(iP,iX-1))/DsMinus
          else
             DfDs1=cHalf*((F_II(iP,iX+1)-F_II(iP,iX  ))/DsPlus+&
                          (F_II(iP,iX  )-F_II(iP,iX-1))/DsMinus)
          end if

          ! here are the parts for \gamma being integrated 
          ! For I_+ the increment equals \frac{\gamma}{(I_{+}+I_{-})}
          ! where \gamma =-4.0*(cPi**2)*Va/k *Integral/cProtonMass 
          ! Integral=\int_{p_{res}(k)}^{\infty}{\frac{dp}{p}*p^4
          ! (p_{res}(k)-\frac{p_{res}^3}{p^2})\frac{\partial{f}}{\partial{s}}
          !
          ! A=p_{res}^{k}=\int{d(ln{p})*(p^{4}*{\frac{\partial{f}}{\partial{s}}}}
          ! B=p_{res}^{3}\int{d(ln{p})*p^2*{\frac{\partial{f}}{\partial{s}}

          A(iP)=A(iP+1)+cHalf*(DfDs0*(P0**4)+DfDs1*(P1**4))*DeltaLnP
          B(iP)=B(iP+1)+cHalf*(DfDs0*(P0**2)+DfDs1*(P1**2))*DeltaLnP

          !Save the values for the lower energy to re-use them as the values for
          !the higher energy end of the next interval (in CYCLING the momentum
          !DECREASES)

          P0=P1; DfDs0=DfDs1

       end do

       !Calculate the wave increment and update the wave spectra
       ExpRhoCompression=exp(RhoCompression_I(iX))
   
       call assign_kolmogorov_spectrum(iX,iX,0,0,X_DI(:,iX:iX),B_I(iX:iX))
       !IPlus_IX(    0,iX) = IPlus_IX(  0,iX)*ExpRhoCompression
       !IMinus_IX(   0,iX) = IMinus_IX( 0,iX)*ExpRhoCompression
       do iK=1,nP 




    !-----------------------------------------------------------------------------------!
    !          Grid in the momentum space                                               !
    !iP     0     1                         nP   nP+1                                   !
    !       |     |    ....                 |     |                                     !
    !P      P_inj P_inj*exp(\Delta (Ln P))  P_Max P_Max*exp(\Delta (Ln P))              !
    !             |    Grid in k-space      |     |                                     !
    !K/B         KMax                      KMin                                         !
    !ik    nP+1   nP                        1     0                                     !
    !-----------------------------------------------------------------------------------!


          K=B_I(iX)*KMin*exp(real(iK-1)*DeltaLnP)   !The wave number

          PRes=cElectronCharge*B_I(iX)/K      !The resonant momentum, for a given K
 

          iP=nP+1-iK

          
          ! here is the dynamic \gamma evaluated

          Gamma=-4.0*2.0*(cPi**2)*VA_I(iX)/K*&
               (PRes*A(iP)-(PRes**3)*B(iP))/    &
               cProtonMass 
          if(DoOutputGamma.and.iX==iXOutputGamma)Gamma_I(iK)=Gamma

          !We need to integrate the two coupled equations:
          !
          !DI_+/Dt =(3/2)* (D ln rho/Dt) * I_+  +  \gamma/(I_+ + I_-) * I_+
          !and
          !DI_-/Dt =(3/2)* (D ln rho/Dt)) *I_-  -  \gamma/(I_+ + I_-) * I_-
          !see Eq.(\ref{eq:Lagrangian})
          !
          !Make a sibstitution 
          !I_+ = I_+(new) * rho^{3/2},    I_- = I_-(new) * rho^{3/2}
          !
          !The modified equations do not involve the density derivative on the right 
          !hand side and are given below. The solution will be expressed in terms of
          !I_+ * I_- and I_+  - I_-, in which we multiply the "old" values of intinsity 
          !by the compression ratio to account for the contribution from the density 
          !derivative.
          !

          C1=(IPlus_IX(iK,iX)*IMinus_IX(iK,iX))*ExpRhoCompression**2 ! Compression  factor comes twice
          C2=(IPlus_IX(iK,iX)-IMinus_IX(iK,iX))*ExpRhoCompression
          
          !The solution of the equations
          !dI_+/dt=(\gamma/(I_+ + I_-)*I_+
          !and
          ! !dI_-/dt=-(\gamma/(I_+  +  I_-)*I_-
          !reads: (I_+   -   I_-)^{n+1}=(I_+   -   I_-)^n + \gamma*dt,
          !so calculate the right hand side:

          C3=gamma*SP_Dt+C2

          !Having in mind the conservation of the c1= I_-*I_+ product, solve the
          !quadratic equation for I_+ :   I_+^2-c3*I_+ - c1=0.

          !To avoid the loss in accuracy, occuring at C3>>C1, calculate that root of the
          !quadratic equation which may be obataibed as a total of two POSITIVE numbers
          
          if(C3>0)then
             
             Iplus_IX(iK,iX)=(C3+sqrt(C3**2+4.0*C1))/2.0 

             !Now we take again into account the conservation of the product I_+*I_-
             !which is the consequence of the relationship \gamma_-=-\gamma_+

             Iminus_IX(iK,iX)=C1/IPlus_IX(iK,iX)
          else
             
             IMinus_IX(iK,iX)=(-C3+sqrt(C3**2+4.0*C1))/2.0 

             !Now we take again into account the conservation of the product I_+*I_-
             !which is the consequence of the relationship \gamma_-=-\gamma_+

             IPlus_IX(iK,iX)=C1/IMinus_IX(iK,iX)
          end if

          if(IPlus_IX(iK,iX)<=cZero)then
             write(*,*)'IPlus(iK,iX) is non-positive, IPlus,iK,iX=',&
                  IPlus_IX(iK,iX),iK,iX,c1,c2,c3,gamma,SP_Dt
             stop
          end if

          if(IMinus_IX(iK,iX)<=cZero)then
             write(*,*)'IMinus(iK,iX) is non-positive, IMinus,iK,iX=',&
                  IMinus_IX(iK,iX),iK,iX,c1,c2,c3,gamma,SP_Dt
             stop
          end if

       end do             !cycling iK
       IPlus_IX( nP+1,iX) = IPlus_IX( nP,iX)
       IMinus_IX(nP+1,iX) = IMinus_IX(nP,iX)

    end do            !cycling iX

    if(UseAdvectionWithAlfvenSpeed)then
       do iK=0,nP+1
          IPlus_IX(iK,1:nX)=VA_I(1:nX)*IPlus_IX(iK,1:nX)/B_I( 1:nX)
          call advance_lin_advection_plus(CFL_I( 1:nX),&
               nX,0,0,IPlus_IX(iK,1:nX))
          IPlus_IX(iK,1:nX)=B_I( 1:nX)*IPlus_IX(iK,1:nX)/VA_I(1:nX)

          IMinus_IX(iK,1:nX)=VA_I(1:nX)*IMinus_IX(iK,1:nX)/B_I( 1:nX)
          call advance_lin_advection_minus(CFL_I(1:nX),&
               nX,0,0,IMinus_IX(iK,1:nX))
          IMinus_IX(iK,1:nX)=B_I( 1:nX)*IMinus_IX(iK,1:nX)/VA_I(1:nX)
       end do
    end if

  end subroutine UpdateSpectrum

  end Module ModTurbulence
