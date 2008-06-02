module ModStatSum
  use ModIonizPotential
  use ModAtomicMass,ONLY : nZMax
  use ModConst
  implicit none  
  integer :: nZ                        !Atomic number of element in question
  integer,dimension(1)::iZDominant    !Most populated ion state
  integer :: iZMin  !Numbers of the ionization states, such that the population
  integer :: iZMax  !of ion states with iZ<iZMin or iZ>iZMax is negligible.
  
 
  real,dimension(1:nZMax) :: IonizPotential_I
  real,dimension(0:nZMax) :: Population_I, StatSumTermLog_I


Contains
  !Set the element and its Ionization Potentials
  !==========================================================================
  subroutine set_element( nZIn)
    integer,intent(in) :: nZIn
    !--------------------------!
    nZ=nZIn
    call get_ioniz_potential(nZ,IonizPotential_I(1:nZ))
  end subroutine set_element
  
  !=========================================================================
  ! Finding the populations of the ion states
  subroutine set_population(TeInv, GeLog)
    real, intent(in) :: TeInv   ! the inverse of the electron temperature, [eV]
    real, intent(in) :: GeLog   ! Natural logarithm of the electron stat weight:
                                !  log(1/(Ne*lambda^3)) !<*>yv:calc.it.ind
    real :: StatSumTermMax,StatSumTermMin


    ! ln(1.0e-2), to truncate terms of the statistical sum, which a factor of 
    ! 1e-2 less than the maximal one: 
    real, parameter :: StatSumToleranceLog = 4.6 

    integer :: iZ
    real :: PITotal
    !--------------------------------------!

    ! First, make the sequence of ln(StatSumTerm) values; let ln(P0)=0 )
    StatSumTermLog_I(0) = 0.	
    do iZ = 1, nZ               !Fill up the sequence using the following equation:
       StatSumTermLog_I(iZ)  =  StatSumTermLog_I(iZ-1)                          &
                              - IonizPotential_I(iZ  )*TeInv + GeLog
     write (*,*)"_set_pop: ", iZ, StatSumTermLog_I(iZ) 
    end do

	
    iZDominant = maxloc(StatSumTermLog_I(0:nZ)) !Find the location of that maximum value
    !debuT       
      write (*,*)"_set_pop:", iZdominant
    
    StatSumTermMax = StatSumTermLog_I(iZDominant(1))
      write (*,*)"_set_pop:", iZdominant

    StatSumTermMin = StatSumTermMax -StatSumToleranceLog
      write (*,*)"_set_pop: iZdom", iZdominant    


                     !Find the lower boundary of the array 
                     !below which the values of Pi can be neglected
    iZMin = count( StatSumTermLog_I(0:iZDominant(1)) <StatSumTermMin) 
      write (*,*)"_set_pop: iZmin", iZmin


                     !Find the similar upper boundary
    iZMax = nZ - count(StatSumTermLog_I(iZDominant(1):nZ) <StatSumTermMin)
      write (*,*)"_set_pop: iZmax", iZmax


                     !Get rid of all the negligible values
    Population_I(0:nZ) = cZero 

    !Convert the array into the Pi values from ln(Pi)
    Population_I(iZMin:iZMax) = exp(StatSumTermLog_I(iZMin:iZMax)-StatSumTermMax)
 
    PITotal = sum(Population_I(iZMin:iZMax))	!Add up all the values of Pi found so far
    Population_I(iZMin:iZMax) = Population_I(iZMin:iZMax)/PITotal !Normalize the Pi-s so that their sum =1
		                        
  end subroutine set_population
  !=======================================!
  !Calculating the Z average values
  real function z_averaged()
    integer::iLoop
    !-------------!
    z_averaged=0
    do iLoop=iZMin,iZMax
       z_averaged = z_averaged + Population_I(iLoop)*real(iLoop)
    end do
  end function z_averaged

  !=======================================!
  !Calculating the Z^2 average values
  real function z2_averaged()
    integer::iLoop
    z2_averaged=0
    do iLoop=iZMin,iZMax
       z2_averaged = z2_averaged + Population_I(iLoop) * real(iLoop)**2
    end do
  end function z2_averaged
  !=======================================!
end module ModStatSum
