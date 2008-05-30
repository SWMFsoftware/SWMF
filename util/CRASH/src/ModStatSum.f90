module ModStatSum
  use ModIonizPotential
  use ModAtomicMass,ONLY : nZMax
  use ModConst
  implicit none  
  integer :: nZ                        !Atomic number of element in question
  integer,dimension(1)::iZDominant    !Most populated ion state
  integer :: iZMin  !Numbers of the ionization states, such that &
  integer :: iZMax  !the population of ion states with iZ<iZMin or &
  integer :: iZ     !iZ>iZMax is negligible.
  real    :: PiSum  !The sum of all non-normalized Pi values  
  real,dimension(nZMax) :: IonizPotential_I
  real,dimension(0:nZMax) :: Population_I
  !real,parameter :: cBoltzmann = 1.3806503e-23 !(m^2*kg)/(s^2*K) - Bolzmann constant (only for testing purposes)


Contains
  !Set the element and its Ionization Potentials
  subroutine set_element(nZIn)
    integer,intent(in) :: nZIn
    nZ=nZIn
    call get_ioniz_potential(nZ,IonizPotential_I(1:nZ))
  end subroutine set_element

  !Finding the populations of the ions
  subroutine set_population(TeInv,GeLog)
         real,intent(in)::TeInv !the inverse of the electron temperature
         real,intent(in)::GeLog !Natural logarithm of the electron stat weight: log(1/(Ne*lambda^3))
         Population_I(0) = 0.	!First, make the sequence of ln(Pi) values; let P0=1 ( ln(P0)=0 )
		 do iZ = 1, nZ  !Fill up the sequence using the following equation:
		   Population_I(iZ) = Population_I(iZ-1) - IonizPotential_I(iZ)*TeInv/cBoltzmann + GeLog
		 end do
		 Population_I(1:nZ) = Population_I(1:nZ) - maxval(Population_I)	!Normalize the array to the maximum value of Pi
		 iZDominant = maxloc(Population_I) !Find the location of that maximum value
		 iZMin = count(Population_I(1:iZDominant(1)) <= -5.) + 1 !Find the lower boundary of the array 
		                                                       !below which the values of Pi can be neglected
		 iZMax = nZ - count(Population_I(iZDominant(1):nZ) <= -5.) !Find the similar upper boundary
		 Population_I(iZMin:iZMax) = exp(Population_I(iZMin:iZMax)) !Convert the array back into the Pi values from ln(Pi)
		 Population_I(1:iZMin-1) = 0. !Get rid of all the negligible values
		 Population_I(iZMax+1:nZMax) = 0.
		 PiSum = sum(Population_I(iZMin:iZMax))	!Add up all the values of Pi found so far
		 Population_I(iZMin:iZMax) = Population_I(iZMin:iZMax)/PiSum !Normalize the Pi-s so that their sum =1
		                        
         !The program sets iZMin and iZMax, such that the ionization state
         !population is negligibly small for iZ<iZMin or iZ>iZMax. It also
         !sets the ionization state populations, Population_I(iZMin:iZMax),
         !such that:
         !log(Population_I(I)/Population_I(I-1))=
         !=-IonizPotential_I(I)*TeInv/kB + GeLog
         !and sum(Pupulation_I(iZMin:iZMax)=1
  end subroutine set_population

  !Calculating the Z average values
  real function z_averaged()
    integer::iLoop
    z_averaged=0
    do iLoop=iZMin,iZMax
       z_averaged = z_averaged + Population_I(iLoop)*iLoop
    end do
  end function z_averaged

  !Calculating the Z^2 average values
  real function z2_averaged()
    integer::iLoop
    z2_averaged=0
    do iLoop=iZMin,iZMax
       z2_averaged = z2_averaged + Population_I(iLoop) * real(iLoop)**2
    end do
  end function z2_averaged

end module ModStatSum
