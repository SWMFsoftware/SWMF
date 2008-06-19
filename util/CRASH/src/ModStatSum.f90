!^CFG COPYRIGHT UM
module ModStatSum
  use ModIonizPotential
  use ModAtomicMass,ONLY : nZMax
  use ModConst

  implicit none  
  integer :: nZ                         !Atomic number of element in question
  integer,dimension(1) :: iZDominant    !Most populated ion state
  integer :: iZMin  !Numbers of the ionization states, such that the population
  integer :: iZMax  !of ion states with iZ<iZMin or iZ>iZMax is negligible.
  
  !Array of ionization potentials - energy needed to create i-level ion from (i-1)-level ion
  real,dimension(1:nZMax) :: IonizPotential_I

  !Array of energies needed to create i-level ion from a neutral atom
  real,dimension(1:nZMax) :: IonizEnergyNeutral_I 
  
  real,dimension(0:nZMax) :: Population_I, StatSumTermLog_I

  real,dimension(nZMax) :: LogN_I,& !array of natural logarithms of consecutive integers
			   N_I !array of consecutive integers (with type real)
  
  real :: C0          ! 2/(Lambda^3)
  real :: ZAv,&       ! the average charge per ion - <Z> (elementary charge units)
          Z2Av,&      ! the average of the squares of the ion charges 
          UInternal,& ! the average internal energy in the plasma
          EAv,&       ! The average energy level of ions
          Cv          ! The heat capacity at constant volume for the plasma
  real,save :: Te = 1.! the electron temperature [eV] (cBoltzmann in eV * Te in Kelvin)
  
Contains
  !=========================================================================
  !Severl initialazations of constants, such as
  !calculating the natural logarithms of the first nZMax integers
  subroutine mod_init
    integer:: iZ  !Used for loops
    real   :: DeBroglieInv
	!-----------------
    LogN_I = (/(log(real(iZ)), iZ = 1,nZMax)/)
    N_I    = (/(real(iZ), iZ = 1,nZMax)/)

    DeBroglieInv=sqrt(cTwoPi*(cElectronMass/cPlanckH)*(cEV/cPlanckH)) 
    !*sqrt(cBoltzmann/cEV * T) - temperature in eV

    C0 = cTwo*DeBroglieInv**3 ! 2/(Lambda^3)
  end subroutine mod_init
  !Set the element and its Ionization Potentials
  !==========================================================================
  subroutine set_element( nZIn)
    integer,intent(in) :: nZIn
	integer :: iZ !for loop
    !--------------------------!
    nZ = nZIn
    call get_ioniz_potential(nZ,IonizPotential_I(1:nZ))
	IonizEnergyNeutral_I(1) = IonizPotential_I(1)
	do iZ = 2,nZ
	   IonizEnergyNeutral_I(iZ) = IonizEnergyNeutral_I(iZ-1) + IonizPotential_I(iZ)
	end do
  end subroutine set_element
  
  !=========================================================================
		  
  
  ! Find the final values of ZAv and the ion populations from Temperature and heavy particle density
  
  subroutine set_ionization_equilibrium(TeIn, Na)
    ! Concentration of heavy particles (atoms+ions) in the plasma 
    ! (# of particles per m^3):
    real, intent(in)::   Na,& ![1/m^3]
	                 TeIn !electron temperature [eV] 
    real :: lnC1,&     ! natural log C1 
	    TeInv ! The inverse of the electron temperature	[1/eV]
												 
    real,parameter :: ToleranceZ = 0.001 !Accuracy of Z needed
    !---------------------------------------------------------
	
    TeInv = cOne / TeIn        ! 1/kT; units: [1/eV]
    lnC1 = log(C0 * sqrt(TeIn)*TeIn / Na)
    call set_Z()	

  contains

    ! Calculating Z averaged iteratively
    subroutine set_Z()
      real    :: ZTrial, Z1 !The trial values of Z for iterations
      integer,dimension(1) :: InitZ !The initial approximation of Z
      integer :: iIter

      !=====================================
      ! First approximate the value of Z by finding for what i=Z 
      ! the derivative of the populations sequence~0 (population is maximum):
      InitZ = minloc(abs(lnC1 - LogN_I(1:nZ) - IonizPotential_I(1:nZ)*TeInv))
                           !Find ZAv in the case when Z~0
      if(InitZ(1)==1)then
         ZAv  = min(real(InitZ(1)),exp(cHalf*(lnC1-IonizPotential_I(1)*TeInv)))
      else
         ZAv  = real(InitZ(1))-cHalf
      end if

      ! Use Newton's method to iteratively get a better approximation of Z:
      iIter  =  0
      ZTrial = -ToleranceZ
      iterations: do while (abs(ZAv-ZTrial) >= ToleranceZ .and. iIter<10)
         ZTrial = ZAv
         call set_population(lnC1 - log(ZTrial))
         Z1  = z_averaged()
         Z2Av  = z2_averaged()
         ZAv = ZTrial - (ZTrial - Z1)/(cOne + (Z2Av - Z1*Z1)/ZTrial)
         iIter = iIter+1
      end do iterations
    end subroutine set_Z

    !==============================================
    ! Finding the populations of the ion states
    subroutine set_population(GeLog)
      real, intent(in) :: GeLog   ! Natural logarithm of the electron stat weight:
                                  !  log(1/(Ne*lambda^3)) !<*>yv:calc.it.ind
      real :: StatSumTermMax,StatSumTermMin


      ! ln(1.0e-2), to truncate terms of the statistical sum, which a factor of 
      ! 1e-2 less than the maximal one:
      real, parameter :: StatSumToleranceLog = 4.6 
      
      integer :: iZ
      real    :: PITotal
      !--------------------------------------!

      ! First, make the sequence of ln(StatSumTerm) values; let ln(P0)=0 )
      StatSumTermLog_I(0) = 0.	
      do iZ = 1, nZ               !Fill up the sequence using the following equation:
         StatSumTermLog_I(iZ)  =  StatSumTermLog_I(iZ-1)                          &
                                - IonizPotential_I(iZ  )*TeInv + GeLog
      end do

      ! Find the location of that maximum value
      iZDominant = maxloc(StatSumTermLog_I(0:nZ))-1 
      
      StatSumTermMax = StatSumTermLog_I(iZDominant(1))
      
      StatSumTermMin = StatSumTermMax -StatSumToleranceLog 
      
      
      !Find the lower boundary of the array 
      !below which the values of Pi can be neglected
      
      iZMin = count( StatSumTermLog_I(0:iZDominant(1)) < StatSumTermMin) 
      
      
      !Find the similar upper boundary
      iZMax = max(nZ - count(StatSumTermLog_I(iZDominant(1):nZ) < StatSumTermMin),1)
      
      !Initialize the population array to zeros
      Population_I(0:nZ) = cZero 
      
      
      !Convert the array into the Pi values from ln(Pi)
      Population_I(iZMin:iZMax) = exp(StatSumTermLog_I(iZMin:iZMax)-StatSumTermMax)
      
      
      PITotal = sum(Population_I(iZMin:iZMax))	!Add up all the values of Pi found so far
      !Normalize the Pi-s so that their sum =1
      Population_I(iZMin:iZMax) = Population_I(iZMin:iZMax)/PITotal 
      
    end subroutine set_population

  end subroutine set_ionization_equilibrium

  !=======================================

  subroutine set_temp(Uin, Na)
    real,intent(in) :: Uin,& !Average internal energy per atomic unit [eV]
	               Na !Density of heavy particles [# of particles/m^3]
    integer :: iIter !iteration counter
    real,parameter :: ToleranceU = 0.001 !accuracy of internal energy needed [(% deviation)/100]
    real :: UDeviation,& !The difference between the given internal energy and the calculated one
            ToleranceUeV
    !-------------------------
    ToleranceUeV = ToleranceU * Uin
    iIter = 0
    UInternal = Uin
    !For initial approximation, use the value for Te found last time
    !Use Newton-Rapson iterations to get a better approximation of Te:
    !UDeviation = ToleranceU	  
    iterations: do 
       call set_ionization_equilibrium(Te, Na) !Find the populations for the trial Te
       EAv = E_averaged() !Find the average of the ionization energry levels of the ions
       UDeviation = internal_energy()-Uin 

       !The exit condition for the loop:
       !(has to be here because it is based on UDeviation)
       if (abs(UDeviation) < ToleranceUeV .or. iIter>10) exit iterations
       
       Cv = heat_capacity() !Find the heat capacity at const. V
       Te = Te - UDeviation/Cv !Calculate the improved value of Te
       iIter = iIter+1
    end do iterations
  end subroutine set_temp

  !============================================
  !calculate the average internal energy per atomic unit from 
  real function internal_energy()
	 internal_energy = 1.50*Te*(1+ZAv) + EAv
  end function internal_energy

  !==================================
  !Calculate the specific heat capacity at constant volume 
  !(derivative of internal energy wrt Te) from temperature:
  real function heat_capacity()
    real :: TeInv,& !The inverse of the electron temperature [1/eV]
            ETeInvAv,&          ! <Ei/Te> (Ei - energy levels, Te - electron temperature [eV])
            DeltaETeInv2Av,&	 ! <(Ei/Te)^2> - <Ei/Te>^2
	    DeltaZ2Av,&          ! <i^2>-<i>^2
   	    DeltaZDeltaETeInvAv ! <i*Ei/Te> - <i><Ei/Te>

    ! Array of energy levels of ions divided by the temperature in eV
    real,dimension(1:nZMax) :: ETeInv_I 
    !------------------
    !calculate the values of the variables defined above:
    TeInv = cOne/Te
    ETeInv_I(iZMin:iZMax) = IonizEnergyNeutral_I(iZMin:iZMax)*TeInv
    ETeInvAv = EAv*TeInv
    DeltaETeInv2Av = sum(Population_I(iZmin:iZmax)*ETeInv_I(iZmin:iZmax)**2)&
                         & - ETeInvAv**2
    DeltaZ2Av = Z2Av - ZAv*ZAv
    DeltaZDeltaETeInvAv = sum(Population_I(iZMin:iZMax)*ETeInv_I(iZmin:iZmax)*N_I(iZMin:iZMax))&
	                 & - ZAv * ETeInvAv
    !calculate the heat capacity:
    heat_capacity = 1.50*(1+ZAv) + DeltaETeInv2Av + &
	          & (3*ZAv*(0.75*DeltaZ2Av + DeltaZDeltaETeInvAv) - DeltaZDeltaETeInvAv**2)/&
			  & (ZAv + DeltaZ2Av)

  end function heat_capacity 

  !=======================================!
  ! Calculating the Z average values from populations
  real function z_averaged()
    z_averaged = sum(Population_I(iZMin:iZMax)*N_I(iZMin:iZMax))
  end function z_averaged

  !=======================================!
  ! Calculating the Z^2 average values from populations
  real function z2_averaged()
    z2_averaged = sum(Population_I(iZMin:iZMax)*N_I(iZMin:iZMax)**2)
  end function z2_averaged

  !==================================
  !Calculate the average ionization energy from neutral atoms of the ions
  real function E_averaged()
     E_averaged = sum(Population_I(iZmin:iZmax)*IonizEnergyNeutral_I(iZmin:iZmax))
  end function E_averaged				  
end module ModStatSum
