!^CFG COPYRIGHT UM
module ModStatSumMix
  use ModIonizPotential
  use ModAtomicMass,ONLY : nZMax
  use ModConst
  implicit none
  SAVE
  logical::DoInit = .true.
  integer::iIterMix     !To provide the output, if needed
  integer,parameter :: nMixMax = 6
  integer :: nMix
  integer,dimension(nMixMax) :: nZ_I=0                      !Atomic number of element in question
  
  integer,dimension(nMixMax) :: iZMin_I  ! Numbers of the ionization states, such that the population
  integer,dimension(nMixMax) :: iZMax_I  ! of ion states with iZ<iZMin or iZ>iZMax is negligible.
  
  real,dimension(nMixMax) :: Concentration_I=0.0 ! relative concentrations of the elements in the mixture (part of whole comprised by the element)
  
  ! Array of ionization potentials - energy needed to create i-level ion from (i-1)-level ion
  real,dimension(1:nZMax,nMixMax) :: IonizPotential_II

  ! Array of energies needed to create i-level ion from a neutral atom
  real,dimension(0:nZMax,nMixMax) :: IonizEnergyNeutral_II 
  
  real,dimension(0:nZMax,nMixMax) :: Population_II ! Array of the populations of ions
  real,dimension(0:nZMax) :: N_I ! array of consecutive integers (with type real)

  real,dimension(nZMax) :: LogN_I ! array of natural logarithms of consecutive integers
  
  real :: C0          ! 2/(Lambda^3)
  real :: ZAvMix,&  ! the average charge per ion - <Z> (elementary charge units)
		  EAvMix,&    ! The average ionization energy level of ions
		  DeltaZ2Av,& ! The value of <(delta i)^2>=<i^2>-<i>^2
		  TeMix=1.0,&   ! the electron temperature [eV] (cBoltzmann in eV * Te in Kelvin)
          Na          ! The density of heavy particles in the plasma
  integer :: iIterTeMix !Temperature iterations counter
  private::mod_init,Na,DeltaZ2Av,C0,DoInit
  !private :: z_averaged, z2_averaged, E_averaged
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

    C0 = 2*DeBroglieInv**3 ! 2/(Lambda^3)
    DoInit=.false.
  end subroutine mod_init
  !Set the element and its Ionization Potentials
  !==========================================================================
  subroutine set_elements( nZIn_I)
    integer,dimension(nMixMax),intent(in) :: nZIn_I
    integer            :: iZ, iMix  ! for loops
    !--------------------------!
    if(DoInit)call mod_init
    if(nZIn_I==nZ_I)return
	nMix = count(nZIn_I(1:nMixMax)>0)
	nZ_I(1:nMix) = nZIn_I
	do iMix=1,nMix
       call get_ioniz_potential(nZ_I(iMix),IonizPotential_I(1:nZ,iMix))
	end do

	IonizEnergyNeutral_II(0,1:nMix) = 0.0
	do iMix=1,nMix
		do iZ = 1,nZ_I(iMix)
			IonizEnergyNeutral_I(iZ,iMix) = IonizEnergyNeutral_I(iZ-1,iMix) + IonizPotential_I(iZ,iMix)
		end do
	end do
  end subroutine set_element
  
  !=========================================================================
		  
  
   ! Find the final values of ZAv and the ion populations from Temperature and heavy particle density
  
  subroutine set_ionization_equilibrium(TeIn, NaIn, ConcentrationIn_I, IsDegenerated )
    ! Concentration of heavy particles (atoms+ions) in the plasma 
    ! (# of particles per m^3):

    real, intent(in)             ::  NaIn,& ![1/m^3]
	                             TeIn !electron temperature [eV]
	real,dimension(nMixMax),intent(in) :: ConcentrationIn_I ! relative concentrations of elements in plasma 
	                                                        !(same order as the elements entered for set_elements)
    logical,optional,intent(out) :: IsDegenerated
    real :: lnC1,&  ! natural log C1 
	    TeInv   ! The inverse of the electron temperature	[1/eV]
												 
    real,parameter :: ToleranceZ = 0.001 !Accuracy of Z needed
    !---------------------------------------------------------
    
    TeMix = TeIn
    Na = NaIn
	Concentration_I(1:nMix) = ConcentrationIn_I	
    TeInv = cOne / TeIn        ! 1/kT; units: [1/eV]
    lnC1  = log(C0 * sqrt(TeIn)*TeIn / Na)
    call set_Z()
    EAvMix = EMix_averaged()	
    if( present(IsDegenerated) ) IsDegenerated = lnC1 -log(ZAvMix) < 2.0

  contains

    ! Calculating Z averaged iteratively
    subroutine set_Z()
      real    :: ZTrial, Z1 ! The trial values of Z for iterations
      real,dimension(iMixMax,1) :: InitZ ! The initial approximation of Z
	  integer :: iMix
      !=====================================
      ! First approximate the value of Z by finding for what i=Z 
      ! the derivative of the populations sequence~0 (population is maximum):
	  initial: do iMix = 1,nMix
		InitZ(iMix) = real(minloc( abs(lnC1 - LogN_I(1:nZ) - IonizPotential_I(1:nZ,iMix)*TeInv) ))
                           !Find ZAv in the case when Z~0
		if(InitZ(iMix,1)==1.0)then
			InitZ(iMix,1)  = min( 1.0, exp(cHalf*( lnC1 -IonizPotential_I(1,iMix)*TeInv) ) )
		else
			InitZ(iMix,1)  = InitZ(iMix, 1) -cHalf
		end if
      end do initial
	  !Average the initial approximation across elements:
	  ZAvMix = sum(Concentration_I(1:nMix) * InitZ(1:nMix,1))
	  
      ! Use Newton's method to iteratively get a better approximation of Z:
      iIterMix  =  0
      ZTrial = -ToleranceZ
      iterations: do while (abs(ZAvMix-ZTrial) >= ToleranceZ .and. iIterMix<10)
         ZTrial = ZAvMix
         call set_populations(lnC1 - log(ZTrial))
         Z1  = zMix_averaged()
		 DeltaZ2Av = 0.0
		 do iMix=1,nMix
             DeltaZ2Av = DeltaZ2Av + Concentration_I(iMix) * sum( Population_I(iZMin_I(iMix):iZMax_I(iMix), iMix) * (N_I(iZMin_I(iMix):iZMax_I(iMix))-Z1)**2 ) 
		 end do
         ZAvMix = ZTrial - (ZTrial - Z1)/(cOne + DeltaZ2Av/ZTrial)
         iIterMix = iIterMix+1
      end do iterations
      ZAvMix=Z1
    end subroutine set_Z

    !==============================================
    ! Finding the populations of the ion states
    subroutine set_populations(GeLog)
      real, intent(in) :: GeLog   ! Natural logarithm of the electron stat weight:
                                  !  log(1/(Ne*lambda^3)) !<*>yv:calc.it.ind
      real :: StatSumTermMax,StatSumTermMin
      real,dimension(0:nZMax) :: StatSumTermLog_I
      integer,dimension(1) :: iZDominant    !Most populated ion state
      ! ln(1.0e-2), to truncate terms of the statistical sum, which a factor of 
      ! 1e-2 less than the maximal one:
      real, parameter :: StatSumToleranceLog = 7.0 
      
      integer :: iZ, iMix
      real    :: PITotal
      !--------------------------------------!
     mixture: do iMix=1,nMix
      ! First, make the sequence of ln(StatSumTerm) values; let ln(P0)=0 )
      StatSumTermLog_I(0) = 0.	
      do iZ = 1, nZ               !Fill up the sequence using the following equation:
         StatSumTermLog_I(iZ)  =  StatSumTermLog_I(iZ-1)                          &
                                - IonizPotential_I(iZ,iMix)*TeInv + GeLog
      end do

      ! Find the location of that maximum value
      iZDominant = maxloc(StatSumTermLog_I(0:nZ))-1 
      
      StatSumTermMax = StatSumTermLog_I(iZDominant(1))
      
      StatSumTermMin = StatSumTermMax -StatSumToleranceLog 
      
      
      ! Find the lower boundary of the array 
      ! below which the values of Pi can be neglected
      
      iZMin_I(iMix) = count( StatSumTermLog_I(0:iZDominant(1)) < StatSumTermMin ) 
      
      
      !Find the similar upper boundary
      iZMax_I(iMix) = max(nZ - count(StatSumTermLog_I(iZDominant(1):nZ) < StatSumTermMin),1)
      
      !Initialize the population array to zeros
      Population_I(0:nZ, iMix) = cZero 
      
      
      !Convert the array into the Pi values from ln(Pi)
      Population_I(iZMin_I(iMix):iZMax_I(iMix), iMix) = exp(StatSumTermLog_I(iZMin_I(iMix):iZMax_I(iMix))-StatSumTermMax)
      
      
      PITotal = sum(Population_I(iZMin_I(iMix):iZMax_I(iMix),iMix))	!Add up all the values of Pi found so far
      !Normalize the Pi-s so that their sum =1
      Population_I(iZMin_I(iMix):iZMax_I(iMix),iMix) = Population_I(iZMin_I(iMix):iZMax_I(iMix),iMix)/PITotal 
     end do mixture  
    end subroutine set_populations

  end subroutine set_ionization_equilibrium

!=======================================  
subroutine set_temperature(Uin, NaIn,ConcentrationIn_I,IsDegenerated)
    real,intent(in) :: Uin,& !Average internal energy per atomic unit [eV]
	               NaIn !Density of heavy particles [# of particles/m^3]
    real,dimension(nMixMax),intent(in) :: ConcentrationIn_I
    logical,intent(out),optional::IsDegenerated
    real,parameter :: ToleranceU = 0.001 !accuracy of internal energy needed [(% deviation)/100]
    real :: UDeviation,& ! The difference between the given internal energy and the calculated one
            ToleranceUeV ! The required accuracy of U in eV
    !-------------------------
    Na = NaIn
	Concentration_I(1:nMix) = ConcentrationIn_I
    ToleranceUeV = ToleranceU * Uin
    iIterTeMix = 0
    !Use Newton-Rapson iterations to get a better approximation of Te:
    UDeviation = 2.*ToleranceUeV
    iterations: do while(abs(UDeviation) >ToleranceUeV .and. iIterTeMix<=20)
       !Find the populations for the trial Te
       call set_ionization_equilibrium(TeMix, Na, Concentration_I, IsDegenerated) 
       UDeviation = internal_energy_mix()-Uin

       !Calculate the improved value of Te, limiting the iterations so 
       !they can't jump too far out
       TeMix = min(2.0*TeMix, max(0.5*TeMix, TeMix - UDeviation/heat_capacity_mix()))  
       iIterTeMix = iIterTeMix+1
    end do iterations
  end subroutine set_temperature
  
  !============================================
  ! Calculate the pressure in the plasma [Pa]
  ! Can only be called after set_ionization_equilibrium has executed
  real function pressure_mix()
     pressure = (1+Zav)*Na*Te*cEV
  end function pressure_mix

  !============================================
  ! calculate the average internal energy per atomic unit [eV]
  ! Can only be called after set_ionization_equilibrium has executed 
  real function internal_energy_mix()
	 internal_energy = 1.50*Te*(1+ZAv) + EAvMix
  end function internal_energy_mix

  !==================================
  !Calculate the specific heat capacity at constant volume 
  !(derivative of internal energy wrt Te) from temperature:
  !Can only be called after set_ionization_equilibrium has executed
  real function heat_capacity_mix()
    real :: TeInv,   & !The inverse of the electron temperature [1/eV]
            ETeInvAv,& ! < Ei/Te> (Ei - energy levels, Te - electron temperature [eV])
            DeltaETeInv2Av,&	! <(delta Ei/Te)^2>
   	    ETeInvDeltaZAv ! <delta i * Ei/Te>

    ! Array of ionization energy levels of ions divided by the temperature in eV
    real,dimension(0:nZMax) :: ETeInv_I 
    !------------------
    !calculate the values of the variables defined above:
    TeInv = cOne/Te
    ETeInv_I(iZMin_I(iMix):iZMax_I(iMix)) = IonizEnergyNeutral_I( iZMin_I(iMix):iZMax_I(iMix) )*TeInv
    ETeInvAv              = EAv*TeInv
    DeltaETeInv2Av        = sum( Population_I(iZMin_I(iMix):iZMax_I(iMix)) * (ETeInv_I(iZMin_I(iMix):iZMax_I(iMix))-ETeInvAv)**2 )
    ETeInvDeltaZAv   = sum( Population_I(iZMin_I(iMix):iZMax_I(iMix)) * ETeInv_I(iZMin_I(iMix):iZMax_I(iMix)) * &
                           (N_I(iZMin_I(iMix):iZMax_I(iMix))-ZAv) )
    ! calculate the heat capacity:
    heat_capacity = 1.50*(cOne +ZAv) + DeltaETeInv2Av &
	           +( 3.0*ZAv*(0.750*DeltaZ2Av + ETeInvDeltaZAv) &
                   - ETeInvDeltaZAv**2 ) / (ZAv + DeltaZ2Av)

  end function heat_capacity_mix ! ^^^^^^ /\ iZMin_I(iMix) >=1 /\ 


  !=======================================!
  ! Calculating the Z average values from populations
  real function zMix_averaged()
    integer :: iMix
	zMix_averaged = 0.0
	do iMix = 1,nMix
       zMix_averaged = zMix_averaged + Concentration_I(iMix) * sum(Population_I(iZMin_I(iMix):iZMax_I(iMix),iMix)*N_I(iZMin_I(iMix):iZMax_I(iMix)))
	end do
  end function z_averaged

   !==================================
  !Calculate the average ionization energy from neutral atoms of the ions
  real function EMix_averaged()
	integer :: iMix
	EMix_averaged = 0.0
	do iMix = 1,nMix
       EMix_averaged = EMix_averaged + Concentration_I(iMix) * sum(Population_I(iZMin_I(iMix):iZMax_I(iMix), iMix) * &
	            IonizEnergyNeutral_I(iZMin_I(iMix):iZMax_I(iMix),iMix))
	end do
  end function EMix_averaged				  
end module ModStatSum
