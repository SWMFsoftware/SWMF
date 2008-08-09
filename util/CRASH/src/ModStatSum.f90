!^CFG COPYRIGHT UM
module ModStatSum
  use ModIonizPotential
  use ModAtomicMass,ONLY : nZMax
  use ModConst
  implicit none
  SAVE
  logical::DoInit = .true.
  integer::iIter     !To provide the output, if needed

  integer :: nZ=-1                      !Atomic number of element in question
  
  integer :: iZMin  ! Numbers of the ionization states, such that the population
  integer :: iZMax  ! of ion states with iZ<iZMin or iZ>iZMax is negligible.
  
  ! Array of ionization potentials - energy needed to create 
  ! i-level ion from (i-1)-level ion
  real,dimension(1:nZMax) :: IonizPotential_I

  ! Array of energies needed to create i-level ion from a neutral atom
  real,dimension(0:nZMax) :: IonizEnergyNeutral_I 
  
  ! Array of the populations of ions
  real,dimension(0:nZMax) :: Population_I

  ! array of consecutive integers (with type real)
  real,dimension(0:nZMax) :: N_I 

  ! array of natural logarithms of consecutive integers
  real,dimension(nZMax) :: LogN_I 
  
  !The combination of fundamental constants, used to calculate the electron
  !statistical weight: electron statistical weight at the temperature of 1eV
  !and the conscentration of 1 particle per m^3
  real,private :: eWight1eV1m3          ! 2/(Lambda^3)

  !Input parameters (note though that the temperature may be not
  !directly assigned and should be found from the internal energy 
  !in this case):
  real :: Te=1.0 ! the electron temperature [eV] = ([cKToeV] * [Te in K])
  real,private::Na  ! The density of heavy particles in the plasma [#/m^3]

  !Averages:
  real :: ZAv,& ! the average charge per ion - <Z> (elementary charge units)
          EAv   ! The average ionization energy level of ions

  !Deviator <(\delta i)^2>:
  real,private :: DeltaZ2Av ! The value of <(delta i)^2>=<i^2>-<i>^2

  integer :: iIterTe !Temperature iterations counter

  real,parameter::c0=0.0

  private::mod_init
Contains
  !=========================================================================
  !Severl initialazations of constants, such as
  !calculating the natural logarithms of the first nZMax integers
  subroutine mod_init
    integer:: iZ  !Used for loops
    real   :: DeBroglieInv
	!-----------------
    LogN_I = (/(log(real(iZ)), iZ = 1,nZMax)/)
    N_I    = (/(real(iZ), iZ = 0,nZMax)/)

    DeBroglieInv = sqrt(cTwoPi*(cElectronMass/cPlanckH)*(cEV/cPlanckH)) 
    !*sqrt(cBoltzmann/cEV * T) - temperature in eV

    eWight1eV1m3 = 2*DeBroglieInv**3 ! 2/(Lambda^3)
    DoInit=.false.
  end subroutine mod_init
 
  !==========================================================================
  !Set the element and its Ionization Potentials
  subroutine set_element( nZIn)
    integer,intent(in) :: nZIn
    integer            :: iZ   ! for loop
    !--------------------------!
    if(DoInit)call mod_init
    if(nZIn==nZ)return
    nZ = nZIn
    call get_ioniz_potential(nZ,IonizPotential_I(1:nZ))

	IonizEnergyNeutral_I(0) = 0.0
	do iZ = 1,nZ
	   IonizEnergyNeutral_I(iZ) = IonizEnergyNeutral_I(iZ-1) +&
                IonizPotential_I(iZ)
	end do
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
    real :: lnC1,&  ! natural log C1 
	    TeInv   ! The inverse of the electron temperature

    real,parameter :: ToleranceZ = 0.001 !Accuracy of Z needed
    !---------------------------------------------------------
    
    Te = TeIn
    Na = NaIn

    !At too low temperature the ionization degree is set to 0
    if(Te<=0.02*IonizPotential_I(1))then
       ZAv=c0; EAv=c0; DeltaZ2Av=c0; return
    end if	

    TeInv = cOne / TeIn        ! 1/kT; units: [1/eV]
    lnC1  = log(eWight1eV1m3 * sqrt(TeIn)*TeIn / Na)
    call set_Z
    EAv = sum(Population_I(iZMin:iZmax)*IonizEnergyNeutral_I(iZMin:iZMax))
    if( present(IsDegenerated) ) IsDegenerated = lnC1 -log(ZAv) < 2.0
  contains
    !=====================================
    ! Calculating Z averaged iteratively
    subroutine set_Z
      real    :: ZTrial  ! The trial values of Z for iterations
      real    :: Z1, Z2  !<z> and <z^2> 
      integer,dimension(1) :: InitZ ! The initial approximation of Z
      !-----------------------------
      
      ! First approximate the value of Z by finding for what i=Z 
      ! the derivative of the populations sequence~0 (population is maximum):
      InitZ = minloc( abs(lnC1 - LogN_I(1:nZ) - IonizPotential_I(1:nZ)*TeInv) )
                           
      if(InitZ(1)==1)then
         !Find ZAv in the case when Z~0
         ZAv  = min( real(InitZ(1) ), &
              exp(cHalf*( lnC1 -IonizPotential_I(1)*TeInv) ) )
      else
         ZAv  = real(InitZ(1)) -cHalf
      end if

      ! Use Newton's method to iteratively get a better approximation of Z:
      iIter  =  0
      ZTrial = -ToleranceZ
      iterations: do while (abs(ZAv-ZTrial) >= ToleranceZ .and. iIter<10)
         ZTrial = ZAv
         call set_population(lnC1 - log(ZTrial))
         Z1  = sum(Population_I(iZMin:iZMax)*N_I(iZMin:iZMax))
         DeltaZ2Av  = sum( Population_I(iZMin:iZmax) * (N_I(iZMin:iZMax)-Z1)**2 ) 
         ZAv = ZTrial - (ZTrial - Z1)/(cOne + DeltaZ2Av/ZTrial)
         iIter = iIter+1
      end do iterations

      ZAv=Z1
    end subroutine set_Z

    !==============================================
    ! Finding the populations of the ion states
    subroutine set_population(GeLog)
      real, intent(in) :: GeLog   ! Natural logarithm of the electron stat weight:
                                  !  log(1/(Ne*lambda^3)) !<*>yv:calc.it.ind
      real :: StatSumTermMax,StatSumTermMin
      real,dimension(0:nZMax) :: StatSumTermLog_I
      integer,dimension(1) :: iZDominant    !Most populated ion state
      ! ln(1.0e-2), to truncate terms of the statistical sum, which a factor of 
      ! 1e-2 less than the maximal one:
      real, parameter :: StatSumToleranceLog = 7.0 
      
      integer :: iZ
      real    :: PITotal
      !--------------------------------------!

      ! First, make the sequence of ln(StatSumTerm) values; let ln(P0)=0 )
      StatSumTermLog_I(0) = c0	
      do iZ = 1, nZ               
         !Fill up the sequence using the following equation:
         StatSumTermLog_I(iZ)  =  StatSumTermLog_I(iZ-1)   &
                                - IonizPotential_I(iZ  )*TeInv + GeLog
      end do

      ! Find the location of that maximum value
      iZDominant = maxloc(StatSumTermLog_I(0:nZ))-1 
      
      StatSumTermMax = StatSumTermLog_I(iZDominant(1))
      
      StatSumTermMin = StatSumTermMax -StatSumToleranceLog 
      
      
      ! Find the lower boundary of the array 
      ! below which the values of Pi can be neglected
      
      iZMin = count( StatSumTermLog_I(0:iZDominant(1)) < StatSumTermMin ) 
      
      
      !Find the similar upper boundary
      iZMax = max(nZ - count(StatSumTermLog_I(iZDominant(1):nZ) < StatSumTermMin),1)
      
      !Initialize the population array to zeros
      Population_I(0:nZ) = c0
      
      
      !Convert the array into the Pi values from ln(Pi)
      Population_I(iZMin:iZMax) = exp(StatSumTermLog_I(iZMin:iZMax)-StatSumTermMax)
      
      !Add up all the values of Pi found so far
      PITotal = sum(Population_I(iZMin:iZMax))	

      !Normalize the Pi-s so that their sum =1
      Population_I(iZMin:iZMax) = Population_I(iZMin:iZMax)/PITotal 
      
    end subroutine set_population

  end subroutine set_ionization_equilibrium

  !=======================================  
  subroutine set_temperature(Uin, NaIn,IsDegenerated)
    real,intent(in) :: Uin,& !Average internal energy per atomic unit [eV]
	               NaIn !Density of heavy particles [# of particles/m^3]
    logical,intent(out),optional::IsDegenerated

    !accuracy of internal energy needed [(% deviation)/100]
    real,parameter :: ToleranceU = 0.001 

    ! The difference between the given internal energy and the calculated one
    real :: UDeviation,& 
            ToleranceUeV ! The required accuracy of U in eV
    !-------------------------
    Na = NaIn
    if(UIn<=0.03*IonizPotential_I(1))then
       Te=UIn/1.50; ZAv=c0; EAv=c0; DeltaZ2Av=c0; return
    end if
    ToleranceUeV = ToleranceU * Uin
    iIterTe = 0
    !Use Newton-Rapson iterations to get a better approximation of Te:
    UDeviation = 2.*ToleranceUeV
    iterations: do while(abs(UDeviation) >ToleranceUeV .and. iIterTe<=20)
       !Find the populations for the trial Te
       call set_ionization_equilibrium(Te, Na, IsDegenerated) 
       UDeviation = internal_energy()-Uin

       !Calculate the improved value of Te, limiting the iterations so 
       !they can't jump too far out, if the initial guess for Te is bad.
       Te = min(2.0*Te, max(0.5*Te, Te - UDeviation/heat_capacity()))  

       iIterTe = iIterTe+1  !Global variable, which is accessible from outside
    end do iterations
  end subroutine set_temperature
  
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
    if(ZAv==c0)then
       heat_capacity=1.50;return
    end if
    
    !calculate the values of the variables defined above:
    TeInv = cOne/Te
    ETeInv_I(iZMin:iZMax) = IonizEnergyNeutral_I( iZMin:iZMax )*TeInv
    ETeInvAv              = EAv*TeInv
    
    !Calculate the missing deviators
    DeltaETeInv2Av        = &
         sum( Population_I(iZMin:iZmax) * (ETeInv_I(iZMin:iZmax)-ETeInvAv)**2 )
    ETeInvDeltaZAv        = &
         sum( Population_I(iZMin:iZMax) * ETeInv_I(iZMin:iZmax) * &
                           (N_I(iZMin:iZMax)-ZAv) )

    ! calculate the heat capacity:
    heat_capacity = 1.50*(cOne +ZAv) + DeltaETeInv2Av &
	           +( 3.0*ZAv*(0.750*DeltaZ2Av + ETeInvDeltaZAv) &
                   - ETeInvDeltaZAv**2 ) / (ZAv + DeltaZ2Av)

  end function heat_capacity 
  !=======================================!
  subroutine get_termodyn_derivatives(GammaOut,GammaSOut,GammaMaxOut,CvOut)
    real,optional,intent(out)::GammaOut    !1+P/UInternal
    real,optional,intent(out)::GammaSOut   !The speed of sound squared*Rho/P
    real,optional,intent(out)::GammaMaxOut !max(Gamma,GammaS)
    real,optional,intent(out)::CvOut       !The same as heat_capacity()
    
    real::Gamma,GammaS,Cv
    !Thermodynamic derivatives: use abbreviations:
    !Ov means over 
    !AtCrho means at constant  rho

    !Derivatives of Z:
    real::TDZOvDTAtCNa  !Te*(\partial Z/\partial Te)_{Na=const}
    real::NaDZOvDNaAtCT !Na*(\partial Z/\partial Na)_{Te=const}

    !Derivatives of pressure:
    real::RhoOvPDPOvDRhoAtCT !(\rho/P)*(\partial P/\Partial \rho)_{T=const}
    real::NaInvDPOvDTAtCNa   !(1/Na)*(\partial P/\partial T)

    !Misc:
    real :: TeInv,   & !The inverse of the electron temperature [1/eV]
            ETeInvAv,& ! < Ei/Te> (Ei - energy levels, Te - electron temperature [eV])
            DeltaETeInv2Av,&	! <(delta Ei/Te)^2>
   	    ETeInvDeltaZAv ! <delta i * Ei/Te>

    ! Array of ionization energy levels of ions divided by the temperature in eV
    real,dimension(0:nZMax) :: ETeInv_I
    real,parameter::Gamma0=5.0/3.0
    !--------------------------------------!
    if(ZAv==c0)then
       if(present(GammaOut))   GammaOut=Gamma0
       if(present(GammaSOut))  GammaSOut=Gamma0
       if(present(GammaMaxOut))GammaMaxOut=Gamma0
       if(present(CvOut))CvOut=1.50
       return
    end if
    if(present(GammaOut))GammaOut = 1.0 + Te * (1.0 + ZAv)/&
         (1.50 * (1.0 + ZAv) * Te + EAv)
    
    !Derivatives at constant T:
    NaDZOvDNaAtCT = - DeltaZ2Av/(1.0 + DeltaZ2Av/ZAv)   
    RhoOvPDPOvDRhoAtCT  = 1.0 + NaDZOvDNaAtCT/(1.0 + Zav) 

    !calculate the values of the variables defined above:
    TeInv = cOne/Te
    ETeInv_I(iZMin:iZMax) = IonizEnergyNeutral_I( iZMin:iZMax )*TeInv
    ETeInvAv              = EAv*TeInv
    
    !Calculate the missing deviator
    ETeInvDeltaZAv        = &
         sum( Population_I(iZMin:iZMax) * ETeInv_I(iZMin:iZmax) * &
                           (N_I(iZMin:iZMax)-ZAv) )

    !Calculate derivatives at const Na:
    TDZOvDTAtCNa =(ETeInvDeltaZAv + 1.50 * DeltaZ2Av)/(1.0 + DeltaZ2Av / ZAv)   
    NaInvDPOvDTAtCNa = 1.0 + ZAv + TDZOvDTAtCNa

    !Calculate another missing deviator
    DeltaETeInv2Av        = &
         sum( Population_I(iZMin:iZmax) * (ETeInv_I(iZMin:iZmax)-ETeInvAv)**2 )
    
    !Calculate Cv: it may be among the output variables and also is used to find Cs
    Cv=1.50 * NaInvDPOvDTAtCNa + DeltaETeInv2Av +  ETeInvDeltaZAv *&
                                                 (1.50 -  TDZOvDTAtCNa/ZAv)
    if(present(CvOut))CvOut=Cv
    GammaS =  RhoOvPDPOvDRhoAtCT +  NaInvDPOvDTAtCNa **2 /(Cv * (1.0 + ZAv))
    if(present(GammaSOut))GammaSOut=GammaS
    if(present(GammaMaxOut))GammaMaxOut=max(GammaS,1.0 + Te * (1.0 + ZAv)/&
         (1.50 * (1.0 + ZAv) * Te + EAv))
  end subroutine get_termodyn_derivatives

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
