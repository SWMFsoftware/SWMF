!^CFG COPYRIGHT UM
module ModStatSum
  use ModStatSumMix, ZAv => ZAvMix, EAv => EAvMix, iIter=>iIterMix, &
       Te => TeMix, iIterTe=>iIterTeMix
  use ModIonizPotential
  use ModAtomicMass,ONLY : nZMax
  use ModConst
  implicit none
  SAVE
  
 

Contains
  !==========================================================================
  !Set the element and its Ionization Potentials
  subroutine set_element( nZIn)
    integer,intent(in) :: nZIn
    !--------------------------!
    call set_mixture(nMixIn=1, nZIn_I=(/nZIn/), ConcentrationIn_I=(/1.0/))
  end subroutine set_element

  !=========================================================================
  ! Find the final values of ZAv and the ion populations, from temperature 
  ! and heavy particle density
  subroutine set_ionization_equilibrium(TeIn, NaIn, IsDegenerated )
    ! Concentration of heavy particles (atoms+ions) in the plasma 
    ! (# of particles per m^3):
    real, intent(in)             ::  NaIn,& ![1/m^3]
         TeIn !electron temperature [eV] 
    logical,optional,intent(out) :: IsDegenerated
   
    call set_ionization_equilibrium_in_mix(TeIn, NaIn, IsDegenerated )
  end subroutine set_ionization_equilibrium

  !=======================================  
  subroutine set_temperature(Uin, NaIn,IsDegenerated)
    real,intent(in) :: Uin,& !Average internal energy per atomic unit [eV]
         NaIn !Density of heavy particles [# of particles/m^3]
    logical,intent(out),optional::IsDegenerated
    !------------------------------------------!
    call set_temperature_in_mix(Uin, NaIn,IsDegenerated)
  end subroutine set_temperature
   

  !============================================
  subroutine pressure_to_temperature(PToNaRatio,NaIn,IsDegenerated)
    real,intent(in) :: PToNaRatio,& !Presure divided by Na [eV]
         NaIn !Density of heavy particles [# of particles/m^3]
    logical,intent(out),optional::IsDegenerated

    !accuracy of internal energy needed [(% deviation)/100]
    real,parameter :: ToleranceP = 0.001 

    ! The difference between the given pressure and the calculated one
    real :: PDeviation,& 
         TolerancePeV ! The required accuracy of P in eV*Na

    !Thermodinamical derivative (\partial P/\partial T)
    !at constant Na, divided by Na
    real :: NaInvDPOvDTAtCNa  
    !------------------------------------

    Na = NaIn

    !It is difficult to make an initial guess about temperature
    !The suggested estimate optimizes the number of idle iterations:
    if(.not.UsePreviousTe) Te = max(1.50* PToNaRatio, IonizPotential_II(1,1) ) * 0.1

    if(PToNaRatio<=0.03*IonizPotential_II(1,1))then
       Te = PToNaRatio; ZAv=0.0; EAv=0.0; DeltaZ2Av=0.0 
       if(present(IsDegenerated))IsDegenerated=.false.
       return
    end if

    TolerancePeV = ToleranceP * PToNaRatio
    iIterTe = 0
    PDeviation = 2.*TolerancePeV

    !Use Newton-Rapson iterations to get a better approximation of Te:

    iterations: do while(abs(PDeviation) >TolerancePeV .and. iIterTe<=20)

       !Find the populations for the trial Te
       call set_ionization_equilibrium(Te, Na, IsDegenerated) 
       PDeviation = pressure()/(Na*cEV) - PToNaRatio

       call get_thermodyn_derivatives(NaInvDPOvDTAtCNaOut=NaInvDPOvDTAtCNa)

       !Calculate the improved value of Te, limiting the iterations so 
       !they can't jump too far out, if the initial guess for Te is bad.
       Te = min(2.0*Te, max(0.5*Te, Te - PDeviation/NaInvDPOvDTAtCNa))  

       iIterTe = iIterTe+1  !Global variable, which is accessible from outside

    end do iterations
  end subroutine pressure_to_temperature
  !============================================
  ! Calculate the pressure in the plasma [Pa]
  ! Can only be called after set_ionization_equilibrium has executed
  real function pressure()
    pressure = (1+Zav)*Na*Te*cEV
  end function pressure

  !============================================
  ! calculate the average internal energy per atomic unit [eV]
  ! Can only be called after set_ionization_equilibrium has executed 
  real function internal_energy()
    internal_energy = 1.50* Te *(1+ZAv) + EAv
  end function internal_energy

  !==================================
  !Calculate the specific heat capacity at constant volume 
  !(derivative of internal energy wrt Te) from temperature:
  !Can only be called after set_ionization_equilibrium has executed
  
  real function heat_capacity()
    real :: TeInv,   & !The inverse of the electron temperature [1/eV]
         ETeInvAv,& ! < Ei/Te> (Ei - energy levels, Te - electron temperature [eV])
         DeltaETeInv2Av,&	! <(delta Ei/Te)^2>
         ETeInvDeltaZAv ! <delta i * Ei/Te>

    ! Array of ionization energy levels of ions divided by the temperature in eV
    real,dimension(0:nZMax) :: ETeInv_I 
    !------------------
  
    if( ZAv<=cTiny )then
       heat_capacity=1.50;return
    end if

    !calculate the values of the variables defined above:
    TeInv = cOne/Te
    ETeInv_I(iZMin_I(1):iZMax_I(1)) = IonizEnergyNeutral_II( iZMin_I(1):iZMax_I(1),1 )*TeInv
    ETeInvAv              = EAv*TeInv

    !Calculate the missing deviators
    DeltaETeInv2Av        = &
         sum( Population_II(iZMin_I(1):iZMax_I(1),1) * (ETeInv_I(iZMin_I(1):iZMax_I(1))-ETeInvAv)**2 )
    ETeInvDeltaZAv        = &
         sum( Population_II(iZMin_I(1):iZMax_I(1),1) * ETeInv_I(iZMin_I(1):iZMax_I(1)) * &
         (N_I(iZMin_I(1):iZMax_I(1))-ZAv) )

    ! calculate the heat capacity:
    heat_capacity = 1.50*(cOne +ZAv) + DeltaETeInv2Av &
         +( 3.0*ZAv*(0.750*DeltaZ2Av + ETeInvDeltaZAv) &
         - ETeInvDeltaZAv**2 ) / (ZAv + DeltaZ2Av)

  end function heat_capacity
  !=======================================!

  subroutine get_thermodyn_derivatives(GammaOut,GammaSOut,GammaMaxOut,CvOut,NaInvDPOvDTAtCNaOut)
    real,optional,intent(out)::GammaOut    !1+P/UInternal
    real,optional,intent(out)::GammaSOut   !The speed of sound squared*Rho/P
    real,optional,intent(out)::GammaMaxOut !max(Gamma,GammaS)
    real,optional,intent(out)::CvOut       !The same as heat_capacity()
    real,optional,intent(out)::NaInvDPOvDTAtCNaOut
    real::Gamma,GammaS,Cv
    !Thermodynamic derivatives: use abbreviations:
    !Ov means over 
    !AtCrho means at constant  rho

    !Derivatives of Z:
    real::TDZOvDTAtCNa  !Te*(\partial Z/\partial Te)_{Na=const}
    real::NaDZOvDNaAtCT !Na*(\partial Z/\partial Na)_{Te=const}

    !Derivatives of pressure:
    real::RhoOvPDPOvDRhoAtCT !(\rho/P)*(\partial P/\Partial \rho)_{T=const}
    real::NaInvDPOvDTAtCNa   !(1/Na)*(\partial P/\partial T)_{\rho=const}

    !Misc:
    real :: TeInv,   & !The inverse of the electron temperature [1/eV]
         ETeInvAv,& ! < Ei/Te> (Ei - energy levels, Te - electron temperature [eV])
         DeltaETeInv2Av,&	! <(delta Ei/Te)^2>
         ETeInvDeltaZAv ! <delta i * Ei/Te>

    ! Array of ionization energy levels of ions divided by the temperature in eV
    real,dimension(0:nZMax) :: ETeInv_I
    real,parameter::Gamma0=5.0/3.0
    !--------------------------------------!

    if(ZAv==0.0)then
       if(present(GammaOut))   GammaOut=Gamma0
       if(present(GammaSOut))  GammaSOut=Gamma0
       if(present(GammaMaxOut))GammaMaxOut=Gamma0
       if(present(CvOut))CvOut=1.50
       if(present(NaInvDPOvDTAtCNaOut))NaInvDPOvDTAtCNaOut = 1.0 
       return
    end if
    if(present(GammaOut))GammaOut = 1.0 + Te * (1.0 + ZAv)/&
         (1.50 * (1.0 + ZAv) * Te + EAv)

    !Derivatives at constant T:
    NaDZOvDNaAtCT = - DeltaZ2Av/(1.0 + DeltaZ2Av/ZAv)   
    RhoOvPDPOvDRhoAtCT  = 1.0 + NaDZOvDNaAtCT/(1.0 + Zav) 

    !calculate the values of the variables defined above:
    TeInv = cOne/Te
    ETeInv_I(iZMin_I(1):iZMax_I(1)) = IonizEnergyNeutral_II( iZMin_I(1):iZMax_I(1),1 )*TeInv
    ETeInvAv              = EAv*TeInv

    !Calculate the missing deviator
    ETeInvDeltaZAv        = &
         sum( Population_II(iZMin_I(1):iZMax_I(1),1) * ETeInv_I(iZMin_I(1):iZMax_I(1)) * &
         (N_I(iZMin_I(1):iZMax_I(1))-ZAv) )

    !Calculate derivatives at const Na:
    TDZOvDTAtCNa =(ETeInvDeltaZAv + 1.50 * DeltaZ2Av)/(1.0 + DeltaZ2Av / ZAv)   
    NaInvDPOvDTAtCNa = 1.0 + ZAv + TDZOvDTAtCNa
    if(present(NaInvDPOvDTAtCNaOut)) NaInvDPOvDTAtCNaOut = NaInvDPOvDTAtCNa


    !Calculate another missing deviator
    DeltaETeInv2Av        = &
         sum( Population_II(iZMin_I(1):iZMax_I(1),1) * (ETeInv_I(iZMin_I(1):iZMax_I(1))-ETeInvAv)**2 )

    !Calculate Cv: it may be among the output variables and also is used to find Cs
    Cv=1.50 * NaInvDPOvDTAtCNa + DeltaETeInv2Av +  ETeInvDeltaZAv *&
         (1.50 -  TDZOvDTAtCNa/ZAv)
    
    if(present(CvOut)) CvOut = Cv
    
    !Define GammaS in terms the speed of sound: GammaS*P/\rho=(\partial P/\partial \rho)_{s=const}
    !Use a formula 3.72 from R.P.Drake,{\it High-Energy Density Physics}, Springer, Berlin-Heidelberg-NY, 2006
    ! and apply the thermodinamic entity: 
    !(\partial \epsilon/\partial \rho)_{T=const}-P/\rho^2)=T(\partial P/\partial T)_{\rho=const)/\rho^2
    
    GammaS =  RhoOvPDPOvDRhoAtCT +  NaInvDPOvDTAtCNa **2 /(Cv * (1.0 + ZAv))
    
    if(present(GammaSOut))GammaSOut = GammaS

    if(present(GammaMaxOut))GammaMaxOut = max( GammaS, 1.0 + Te * (1.0 + ZAv)/&
         (1.50 * (1.0 + ZAv) * Te + EAv))

  end subroutine get_thermodyn_derivatives

  !=======================================!
  ! Calculating the Z average values from populations
  real function z_averaged()
    z_averaged = ZAv
  end function z_averaged

  !=======================================!
  ! Calculating the Z^2 average values from populations
  real function z2_averaged()
    !The average of square is the averaged squared plus dispersion squared
    z2_averaged = DeltaZ2Av + ZAv*ZAv 
  end function z2_averaged

  !==================================
  !Calculate the average ionization energy from neutral atoms of the ions
  real function E_averaged()
    E_averaged = EAv

  end function E_averaged

end module ModStatSum
