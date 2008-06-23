   !^CFG COPYRIGHT UM
module ModStatSum
  use ModIonizPotential
  use ModAtomicMass,ONLY : nZMax
  use ModConst
  implicit none
  SAVE
  integer :: nZ=-1                      !Atomic number of element in question
  
  integer :: iZMin  !Numbers of the ionization states, such that the population
  integer :: iZMax  !of ion states with iZ<iZMin or iZ>iZMax is negligible.
  
  !Array of ionization potentials - energy needed to create i-level ion from (i-1)-level ion
  real,dimension(1:nZMax) :: IonizPotential_I

  !Array of energies needed to create i-level ion from a neutral atom
  real,dimension(1:nZMax) :: IonizEnergyNeutral_I 
  
  real,dimension(0:nZMax) :: Population_I,& !Array of the populations of ions
                             N_I !array of consecutive integers (with type real)

  real,dimension(nZMax) :: LogN_I !array of natural logarithms of consecutive integers
  
  real :: C0          ! 2/(Lambda^3)
  real :: ZAv,&       ! the average charge per ion - <Z> (elementary charge units)
          EAv,&       ! The average ionization energy level of ions
          Te = 1.,&   ! the electron temperature [eV] (cBoltzmann in eV * Te in Kelvin)
          Na          ! The density of heavy particles in the plasma
  
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

    DeBroglieInv=sqrt(cTwoPi*(cElectronMass/cPlanckH)*(cEV/cPlanckH)) 
    !*sqrt(cBoltzmann/cEV * T) - temperature in eV

    C0 = cTwo*DeBroglieInv**3 ! 2/(Lambda^3)
  end subroutine mod_init
  !Set the element and its Ionization Potentials
  !==========================================================================
  subroutine set_element( nZIn)
    integer,intent(in) :: nZIn
    integer            :: iZ   ! for loop
    !--------------------------!
    if(nZIn==nZ)return
    nZ = nZIn
    call get_ioniz_potential(nZ,IonizPotential_I(1:nZ))

	IonizEnergyNeutral_I(1) = IonizPotential_I(1)
	do iZ = 2,nZ
	   IonizEnergyNeutral_I(iZ) = IonizEnergyNeutral_I(iZ-1) + IonizPotential_I(iZ)
	end do
  end subroutine set_element
  
  !=========================================================================
		  
  
  ! Find the final values of ZAv and the ion populations from Temperature and heavy particle density
  
  subroutine set_ionization_equilibrium(TeIn, NaIn, IsDegenerated )
    ! Concentration of heavy particles (atoms+ions) in the plasma 
    ! (# of particles per m^3):

    real, intent(in)             ::  NaIn,& ![1/m^3]
	                             TeIn !electron temperature [eV] 
    logical,optional,intent(out) :: IsDegenerated
    real :: lnC1,&  ! natural log C1 
	    TeInv   ! The inverse of the electron temperature	[1/eV]
												 
    real,parameter :: ToleranceZ = 0.001 !Accuracy of Z needed
    !---------------------------------------------------------
    
    Te = TeIn
    Na = NaIn	
    TeInv = cOne / TeIn        ! 1/kT; units: [1/eV]
    lnC1  = log(C0 * sqrt(TeIn)*TeIn / Na)
    call set_Z()
    EAv = E_averaged()	
    if( present(IsDegenerated) ) IsDegenerated = lnC1 -log(ZAv) < 2.0

  contains

    ! Calculating Z averaged iteratively
    subroutine set_Z()
      real    :: ZTrial, Z1, Z2 ! The trial values of Z for iterations
      integer,dimension(1) :: InitZ ! The initial approximation of Z
      integer :: iIter
      !=====================================
      ! First approximate the value of Z by finding for what i=Z 
      ! the derivative of the populations sequence~0 (population is maximum):
      InitZ = minloc( abs(lnC1 - LogN_I(1:nZ) - IonizPotential_I(1:nZ)*TeInv) )
                           !Find ZAv in the case when Z~0
      if(InitZ(1)==1)then
         ZAv  = min( real(InitZ(1) ), exp(cHalf*( lnC1 -IonizPotential_I(1)*TeInv) ) )
      else
         ZAv  = real(InitZ(1)) -cHalf
      end if

      ! Use Newton's method to iteratively get a better approximation of Z:
      iIter  =  0
      ZTrial = -ToleranceZ
      iterations: do while (abs(ZAv-ZTrial) >= ToleranceZ .and. iIter<10)
         ZTrial = ZAv
         call set_population(lnC1 - log(ZTrial))
         Z1  = z_averaged()
         Z2  = z2_averaged()
         ZAv = ZTrial - (ZTrial - Z1)/(cOne + (Z2 - Z1*Z1)/ZTrial)
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
      
      iZMin = count( StatSumTermLog_I(0:iZDominant(1)) < StatSumTermMin ) 
      
      
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
subroutine set_temperature(Uin, NaIn)
    real,intent(in) :: Uin,& !Average internal energy per atomic unit [eV]
	               NaIn !Density of heavy particles [# of particles/m^3]
    integer :: iIter !iteration counter
    real,parameter :: ToleranceU = 0.001 !accuracy of internal energy needed [(% deviation)/100]
    real :: UDeviation,& !The difference between the given internal energy and the calculated one
            ToleranceUeV ! The required accuracy of U in eV
    !-------------------------
    Na = NaIn
    ToleranceUeV = ToleranceU * Uin
    iIter = 0
    !For initial approximation, use the value for Te found last time
    !Use Newton-Rapson iterations to get a better approximation of Te:
    !UDeviation = ToleranceU	  
    iterations: do 
       call set_ionization_equilibrium(Te, Na) !Find the populations for the trial Te
       UDeviation = internal_energy()-Uin 

       !The exit condition for the loop:
       !(has to be here because it is based on UDeviation)
       if (abs(UDeviation) < ToleranceUeV .or. iIter>10) exit iterations
       
       Te = Te - UDeviation/heat_capacity() !Calculate the improved value of Te
       iIter = iIter+1
    end do iterations
  end subroutine set_temperature
  
  !============================================
  !Calculate the pressure in the plasma [Pa]
  !Can only be called after set_ionization_equilibrium has executed
  real function pressure()
     pressure = (1+Zav)*Na*Te*cEV
  end function pressure

  !============================================
  !calculate the average internal energy per atomic unit [eV]
  !Can only be called after set_ionization_equilibrium has executed 
  real function internal_energy()
	 internal_energy = 1.50*Te*(1+ZAv) + EAv
  end function internal_energy

  !==================================
  !Calculate the specific heat capacity at constant volume 
  !(derivative of internal energy wrt Te) from temperature:
  !Can only be called after set_ionization_equilibrium has executed
  real function heat_capacity()
    real :: TeInv,& !The inverse of the electron temperature [1/eV]
            ETeInvAv,&          ! < Ei/Te> (Ei - energy levels, Te - electron temperature [eV])
            DeltaETeInv2Av,&	! <(Ei/Te)^2> - <Ei/Te>^2
	      DeltaZ2Av,&         ! <i^2>-<i>^2
   	      DeltaZDeltaETeInvAv ! <i*Ei/Te> - <i><Ei/Te>

    ! Array of ionization energy levels of ions divided by the temperature in eV
    real,dimension(1:nZMax) :: ETeInv_I 
    !------------------
    !calculate the values of the variables defined above:
    TeInv = cOne/Te
    ETeInv_I(iZMin:iZMax) = IonizEnergyNeutral_I( iZMin:iZMax )*TeInv
    ETeInvAv              = EAv*TeInv
    DeltaETeInv2Av        = sum( Population_I(iZMin:iZmax) * (ETeInv_I(iZMin:iZmax)-ETeInvAv)**2 )
    DeltaZ2Av             = sum( Population_I(iZMin:iZmax) * (N_I(iZMin:iZMax)-Zav)**2 )
    DeltaZDeltaETeInvAv   = sum( Population_I(iZMin:iZMax) * (ETeInv_I(iZMin:iZmax)-ETeInvAv) * &
                           (N_I(iZMin:iZMax)-ZAv) )
    !calculate the heat capacity:
    heat_capacity = 1.50*(cOne +ZAv) + DeltaETeInv2Av &
	           +( 3.0*ZAv*(0.750*DeltaZ2Av + DeltaZDeltaETeInvAv) &
                   - DeltaZDeltaETeInvAv**2 ) / (ZAv + DeltaZ2Av)

  end function heat_capacity ! ^^^^^^ /\ iZmin >=1 /\ 


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
     E_averaged = sum(Population_I(iZMin:iZmax)*IonizEnergyNeutral_I(iZMin:iZmax))
  end function E_averaged				  
end module ModStatSum
