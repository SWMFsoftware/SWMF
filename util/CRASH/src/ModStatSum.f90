
module ModStatSum
  use ModIonizPotential
  use ModAtomicMass,ONLY : nZMax
  use ModConst

  implicit none  
  integer :: nZ                         !Atomic number of element in question
  integer,dimension(1) :: iZDominant    !Most populated ion state
  integer :: iZMin  !Numbers of the ionization states, such that the population
  integer :: iZMax  !of ion states with iZ<iZMin or iZ>iZMax is negligible.
  
 
  real,dimension(1:nZMax) :: IonizPotential_I,&	!array of ionization potentials - energy needed to create i-level ion from (i-1)-level ion
							 EnergyLevel_I !array of energy levels of ions - energy needed to create i-level ion from a neutral atom
  real,dimension(0:nZMax) :: Population_I, StatSumTermLog_I

  real,parameter :: cBoltzmannEVInv = cEV/cBoltzmann !Inverse of the Boltzmann constant in [K/eV]
  real,dimension(nZMax) :: LogN_I
  real :: C0 ! 2/(Lambda^3)
  real :: TeInv,& ! the inverse of the electron temperature [1/eV] (1/(cBoltzmann in eV * Te in Kelvin)
		  Zav  ! the average charge per ion - <Z> (elementary charge units)
Contains
  !=========================================================================
  !Calculates the natural logarithms of the first nZMax integers
  subroutine mod_init
    integer:: iZ  !Used for loops
    real   :: DeBroglie
    LogN_I = (/(log(real(iZ)), iZ = 1,nZMax)/)
    DeBroglie=sqrt(cTwoPi*(cElectronMass/cPlanckH)*(cBoltzmann/cPlanckH))
    C0 = cTwo*DeBroglie**3 ! 2/(Lambda^3)
  end subroutine mod_init
  !Set the element and its Ionization Potentials
  !==========================================================================
  subroutine set_element( nZIn)
    integer,intent(in) :: nZIn
	integer :: iZ !for loop
    !--------------------------!
    nZ = nZIn
    call get_ioniz_potential(nZ,IonizPotential_I(1:nZ))
	EnergyLevel_I(1) = IonizPotential_I(1)
	do iZ = 2,nZ
	   EnergyLevel_I(iZ) = EnergyLevel_I(iZ-1) + IonizPotential_I(iZ)
	end do
	!EnergyLevel_I(1:nZ) = (/(sum(IonizPotential_I(1:iZ)), iZ = 1,nZ)/)
  end subroutine set_element
  
  !=========================================================================
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
      write (*,*)"_set_pop: terMin:",  StatSumTermMin  ," terMax:",  StatSumTermMax    


    !Find the lower boundary of the array 
    !below which the values of Pi can be neglected

    iZMin = count( StatSumTermLog_I(0:iZDominant(1)) < StatSumTermMin) 


    !Find the similar upper boundary
    iZMax = max(nZ - count(StatSumTermLog_I(iZDominant(1):nZ) < StatSumTermMin),1)

    !Get rid of all the negligible values
    Population_I(0:nZ) = cZero 


    !Convert the array into the Pi values from ln(Pi)
    Population_I(iZMin:iZMax) = exp(StatSumTermLog_I(iZMin:iZMax)-StatSumTermMax)


    PITotal = sum(Population_I(iZMin:iZMax))	!Add up all the values of Pi found so far
    !Normalize the Pi-s so that their sum =1
    Population_I(iZMin:iZMax) = Population_I(iZMin:iZMax)/PITotal 

  end subroutine set_population


  !=======================================!
  ! Calculating the Z average values
  real function z_averaged()
    integer::iLoop
    !-------------!
    z_averaged=0
    do iLoop=iZMin,iZMax
       z_averaged = z_averaged + Population_I(iLoop)*real(iLoop)
    end do
  end function z_averaged

  !=======================================!
  ! Calculating the Z^2 average values
  real function z2_averaged()
    integer::iLoop
    z2_averaged=0
    do iLoop=iZMin,iZMax
       z2_averaged = z2_averaged + Population_I(iLoop) * real(iLoop)**2
    end do
  end function z2_averaged
!=======================================!
		  
  
  ! Find the final values of Zav and the ion populations
  
  subroutine set_ionization_equilibrium(Na, Te)
    ! Concentration of heavy particles (atoms+ions) in the plasma 
    ! (# of particles per m^3)
    real, intent(in)::   Na   ! {1/m^3}  
    real, intent(in)::   Te   ! Electron temperature (Kelvin)
    real :: lnC1     ! natural log C1 
												 
!    real,parameter :: ToleranceZ = 0.01 !Accuracy of Z needed
    real,parameter :: ToleranceZ = 0.001 !Accuracy of Z needed


    !==========================================
    write(*,*)'Start set_ionization_equilibrium', Na,Te

    TeInv = cBoltzmannEVInv / Te        ! 1/kT; units: [1/eV]
    lnC1 = log(C0 * sqrt(Te)*Te / Na)
    call set_Z()	
  contains
    ! Calculating Z averaged iteratively
    subroutine set_Z()
      real    :: ZTrial, Z1, Z2 !The trial values of Z for iterations
      integer,dimension(1) :: InitZ !The initial approximation of Z
      integer :: iIter

      !=====================================
      ! First approximate the value of Z by finding for what i=Z 
      ! the derivative of the populations sequence~0 (population is maximum):
      InitZ = minloc(abs(lnC1 - LogN_I(1:nZ) - IonizPotential_I(1:nZ)*TeInv))
                           !Find Zav in the case when Z~0
      if(InitZ(1)==1)then
         Zav  = min(real(InitZ(1)),exp(cHalf*(lnC1-IonizPotential_I(1)*TeInv)))
      else
         Zav  = real(InitZ(1))-cHalf
      end if

      write(*,*) "Initial approximation of Z:", Zav

      ! Use Newton's method to iteratively get a better approximation of Z:
      iIter  =  0
      ZTrial = -ToleranceZ
      iterations: do while (abs(Zav-ZTrial) >= ToleranceZ.and.iIter<10)
         ZTrial = Zav
         call set_population(lnC1 - log(ZTrial))
         Z1  = z_averaged()
         Z2  = z2_averaged()
         Zav = ZTrial - (ZTrial - Z1)/(cOne + (Z2 - Z1*Z1)/ZTrial)
         iIter = iIter+1
      end do iterations
      write (*,*) "Iterations done:", iIter
      write (*,*) "Final Zav = ", Zav
    end subroutine set_Z
  end subroutine set_ionization_equilibrium

  !=======================================
  real function internal_energy()
      real :: Uparticle
	  integer :: iZ	 
	  !-------------
	  Uparticle = 1.50/TeInv
	  internal_energy = Uparticle * (1+Zav) + sum(Population_I(iZmin:iZmax)*EnergyLevel_I(iZmin:iZmax))
  end function internal_energy

end module ModStatSum
