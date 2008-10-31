!^CFG COPYRIGHT UM
module ModStatSumMix
  use ModIonizPotential
  use ModAtomicMass,ONLY : nZMax
  use ModConst
  implicit none
  SAVE

  logical :: DoInit = .true.             !Set to false after the initialization

  integer :: nMix = -1
  integer, parameter :: nMixMax = 6

  integer,private,dimension(nMixMax) :: nZ_I = -1  !Atomic numbers of elements in the mizture 
  
  !\
  !For each component in the mixture:
  !/
  
  integer,dimension(nMixMax) :: iZMin_I  ! Numbers of the ionization states, such that the population
  integer,dimension(nMixMax) :: iZMax_I  ! of ion states with iZ<iZMin or iZ>iZMax is negligible.
  
  ! relative concentrations of the elements in the mixture (part of whole comprised by the element)
  real,dimension(nMixMax) :: Concentration_I=0.0 
  
  ! Array of ionization potentials - energy needed to create i-level ion from (i-1)-level ion
  real,dimension(1:nZMax,nMixMax) :: IonizPotential_II

  ! Array of energies needed to create i-level ion from a neutral atom
  real,dimension(0:nZMax,nMixMax) :: IonizEnergyNeutral_II 
  
  !\
  !The distribution over the ion states, for each of the elements:
  !/
  real,dimension(0:nZMax,nMixMax) :: Population_II ! Array of the populations of ions

  
  !The combination of fundamental constants, used to calculate the electron
  !statistical weight: electron statistical weight at the temperature of 1eV
  !and the conscentration of 1 particle per m^3
  real,private :: eWight1eV1m3           ! 2/(Lambda^3)


  !Input parameters (note though that the temperature may be not
  !directly assigned and should be found from the internal energy 
  !or pressure in this case):
  real :: TeMix=1.0 ! the electron temperature [eV] = ([cKToeV] * [Te in K])
  real ::Na  ! The density of heavy particles in the plasma [#/m^3]


  !Averages:
  integer,parameter :: DeltaI_ = 1, DeltaE_=2
  real :: ZAvMix,&  ! 1 the average charge per ion - <Z> (elementary charge units)
          EAvMix          ! 2 The average ionization energy level of ions

  !Deviators <(\delta i)^2>, <delta i delta E> and <(delta e)^2:
  real :: DeltaZ2Av ! The value of <(delta i)^2>=<i^2>-<i>^2


  integer :: iIterMix     !To provide the output for the convergence efficiency, if needed
  integer :: iIterTeMix !Temperature iterations counter


  !Auxiliary arrys of numerical constants:
  real,dimension(0:nZMax) :: N_I ! array of consecutive integers (with type real)
  real,dimension(nZMax)   :: LogN_I ! array of natural logarithms of consecutive integers
 
  logical :: UsePreviousTe = .true.

  private :: mod_init

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
  !Set the element and its Ionization Potentials
  !==========================================================================
  subroutine set_mixture(nMixIn, nZIn_I, ConcentrationIn_I)
    integer,intent(in)::nMixIn
    integer,dimension(nMixIn),intent(in) :: nZIn_I
    real,dimension(nMixIn),intent(in) :: ConcentrationIn_I

    integer            :: iZ, iMix  ! for loops
    !--------------------------!
    
    if(DoInit)call mod_init
    
    if(nMixIn==nMix)then
       if( all( nZIn_I( 1:nMix ) == nZ_I( 1:nMix ) ).and. &
            all( ConcentrationIn_I( 1:nMix ) == Concentration_I( 1:nMix ) ) )return
    end if
    
    nMix = nMixIn
    nZ_I( 1:nMix ) = nZIn_I( 1:nMix )
    do iMix=1,nMix
       call get_ioniz_potential(nZ_I(iMix),IonizPotential_II(1:nZ_I(iMix),iMix))
    end do

    IonizEnergyNeutral_II(0,1:nMix) = 0.0
    do iMix=1,nMix
       do iZ = 1,nZ_I(iMix)
          IonizEnergyNeutral_II(iZ,iMix) = IonizEnergyNeutral_II(iZ-1,iMix) + IonizPotential_II(iZ,iMix)
       end do
    end do
    
    Concentration_I( 1:nMix ) =  ConcentrationIn_I( 1:nMix )
    if(abs (sum(Concentration_I( 1:nMix ) ) - 1.0 ) > cTiny)then
       write(*,*)'Wrong input - the total of relative concentrations differs from 1.0: ', Concentration_I( 1:nMix ) 
       call CON_stop('Stop')
    end if
  end subroutine set_mixture
  
  !=========================================================================
		  
  
   ! Find the final values of ZAv and the ion populations from Temperature and heavy particle density
  
  subroutine set_ionization_equilibrium_in_mix(TeIn, NaIn, IsDegenerated )
    ! Concentration of heavy particles (atoms+ions) in the plasma 
    ! (# of particles per m^3):
    real, intent(in)             ::  NaIn,& ![1/m^3]
	                             TeIn !electron temperature [eV]
    logical,optional,intent(out) :: IsDegenerated


    !Internal variables
    real :: lnC1,&  ! natural log C1 
	    TeInv   ! The inverse of the electron temperature	[1/eV]
												 
    real,parameter :: ToleranceZ = 0.001 !Accuracy of Z needed
    !---------------------------------------------------------
    
    TeMix = TeIn
    Na = NaIn

    if( TeMix <= 0.02 * minval( IonizPotential_II( 1,1:nMix) ) )then
       
       ZAvMix=0.0; EAvMix=0.0; DeltaZ2Av=0.0
       
       if(present(IsDegenerated))IsDegenerated=.false.
       return
    end if

	
    TeInv = cOne / TeIn        ! 1/kT; units: [1/eV]
    lnC1  = log( eWight1eV1m3  * sqrt(TeIn)*TeIn / Na)

    call set_Z

    EAvMix = EMix_averaged()	

    if( present(IsDegenerated) ) IsDegenerated = lnC1 -log(ZAvMix) < 2.0

  contains

    ! Calculating Z averaged iteratively
    subroutine set_Z
      !Internal variables
      real    :: ZTrial, Z1              ! The trial values of Z for iterations
      real,dimension(nMixMax) :: InitZ_I ! The initial approximation of Z
      real :: ZJ                         ! Average Z for each component
      integer :: iMix,iZ(1)              ! Misc
      !----------------------------------------------------------
      ! First approximate the value of Z by finding for what i=Z 
      ! the derivative of the populations sequence~0 (population is maximum):
      initial: do iMix = 1,nMix

         iZ = minloc( abs(lnC1 - LogN_I(1:nZ_I(iMix)) - &
              IonizPotential_II(1:nZ_I(iMix),iMix)*TeInv) )
        
        
         if(iZ(1)==1)then
            !Find ZAv in the case when Z<=1
            InitZ_I(iMix)  = min( 1.0, exp(cHalf*( lnC1 -IonizPotential_II(1,iMix)*TeInv) ) )
         else
            !Apply the above estimate
            InitZ_I(iMix)  = real(iZ( 1)) -cHalf
        
         end if
      
      end do initial
      
      !Average the initial approximation across the mixture compenents:
      ZAvMix = sum(Concentration_I(1:nMix) * InitZ_I(1:nMix))
      
      ! Use Newton's method to iteratively get a better approximation of Z:
      ! Organize the iteration loop:
      iIterMix  =  0
      ZTrial = -ToleranceZ
      iterations: do while (abs(ZAvMix-ZTrial) >= ToleranceZ .and. iIterMix<10)
         
         !Set a trial value for Z
         ZTrial = ZAvMix

         !Set the ion population with the electron stqtistical weight being expressed
         !in terms of ZTrial:
         call set_populations(lnC1 - log(ZTrial))
         
         !Take an average over the ion populations:
         Z1  = zMix_averaged()
		
         DeltaZ2Av = 0.0
         do iMix=1,nMix
            !Calculate averaged Z for a given component of the mixture
            ZJ = sum( Population_II(iZMin_I(iMix):iZMax_I(iMix), iMix) * N_I(iZMin_I(iMix):iZMax_I(iMix)))
            
            !Calculate a << (\delta i)^2 >>
            DeltaZ2Av = DeltaZ2Av + Concentration_I(iMix) * &
                 sum( Population_II(iZMin_I(iMix):iZMax_I(iMix), iMix) * (N_I(iZMin_I(iMix):iZMax_I(iMix))-ZJ)**2 ) 
         end do
         
         ZAvMix = ZTrial - (ZTrial - Z1)/(cOne + DeltaZ2Av/ZTrial)
         
         iIterMix = iIterMix+1
      end do iterations
      ZAvMix=Z1
    end subroutine set_Z

    !==============================================
    ! Finding the populations of the ion states
    subroutine set_populations(GeLogIn)
      real, intent(in) :: GeLogIn ! Natural logarithm of the electron stat weight:
                                  !  log(1/(Ne*lambda^3)) 
  
      real :: StatSumTermMax,StatSumTermMin
      
      real,dimension(0:nZMax) :: StatSumTermLog_I
      
      integer,dimension(1) :: iZDominant    !Most populated ion state

      ! ln(1.0e-3), to truncate terms of the statistical sum, which a factor of 
      ! 1e-3 less than the maximal one:
      real, parameter :: StatSumToleranceLog = 7.0 
      
      integer :: iZ, iMix  !Loop variable
      real    :: PITotal   !Normalization constant, to set the populations total to equal 1.
      real    :: GeLog
      !--------------------------------------!

      ! The present version does not stop simulation
      ! when the Fermi degeneracy occurs (the electron 
      ! statictical weight is not large), but it sets the
      ! low bound for the weight to be e^2\approx 8
      
      GeLog = max( 2.0, GeLogIn )

      mixture: do iMix=1,nMix

         ! First, make the sequence of ln(StatSumTerm) values; let ln(P0)=0 )
         StatSumTermLog_I(0) = 0.0

         do iZ = 1, nZ_I(iMix)       !Fill up the sequence using the following equation:

            StatSumTermLog_I(iZ)  =  StatSumTermLog_I(iZ-1)                          &
                                   - IonizPotential_II(iZ,iMix)*TeInv + GeLog
         end do

         ! Find the location of that maximum value
         iZDominant = maxloc(StatSumTermLog_I(0:nZ_I(iMix)))-1 
      
         StatSumTermMax = StatSumTermLog_I(iZDominant(1))
      
         StatSumTermMin = StatSumTermMax -StatSumToleranceLog 
      
      
         ! Find the lower boundary of the array 
         ! below which the values of Pi can be neglected
         
         iZMin_I(iMix) = count( StatSumTermLog_I(0:iZDominant(1)) < StatSumTermMin ) 
      
      
         !Find the similar upper boundary
         iZMax_I(iMix) = max(nZ_I(iMix) - count(StatSumTermLog_I( iZDominant(1):nZ_I(iMix) ) < StatSumTermMin),1)
      
         !Initialize the population array to zeros
         Population_II(0:nZ_I(iMix), iMix) = 0.0
      
      
         !Convert the array into the Pi values from ln(Pi)
         Population_II(iZMin_I(iMix):iZMax_I(iMix), iMix) = exp(StatSumTermLog_I(iZMin_I(iMix):iZMax_I(iMix))-StatSumTermMax)
      
         !Normalize the Pi-s so that their sum =1

         PITotal = sum(Population_II(iZMin_I(iMix):iZMax_I(iMix),iMix))	!Add up all the values of Pi found so far         
         Population_II(iZMin_I(iMix):iZMax_I(iMix),iMix) = Population_II(iZMin_I(iMix):iZMax_I(iMix),iMix)/PITotal
 
      end do mixture
    end subroutine set_populations

  end subroutine set_ionization_equilibrium_in_mix

  !=======================================  
  !Find a temperature from the internal energy per atom:
  
  subroutine set_temperature_in_mix(Uin, NaIn,IsDegenerated)
    
    real,intent(in) :: Uin,& !Average internal energy per atomic unit [eV]
                       NaIn !Density of heavy particles [# of particles/m^3]
    logical,intent(out),optional::IsDegenerated
    
    real,parameter :: ToleranceU = 0.001 !accuracy of internal energy needed [(% deviation)/100]
    real :: UDeviation,& ! The difference between the given internal energy and the calculated one
            ToleranceUeV ! The required accuracy of U in eV
    !-------------------------
    Na = NaIn

    if( .not.UsePreviousTe ) then
       if(nMix>1) then 
          !To keep the backward compatibility, because the choice for mixture was different
          TeMix = UIn / 1.5
       else
          TeMix = max(UIn,IonizPotential_II(1,1)) * 0.1
       end if
    end if
    if( UIn <= 0.03 * minval(IonizPotential_II(1,1:nMix))) then

       TeMix=UIn/1.50 ;  ZAvMix=0.0 ;  EAvMix=0.0; DeltaZ2Av=0.0
       
       if(present(IsDegenerated))IsDegenerated=.false.
       return
    end if

    ToleranceUeV = ToleranceU * Uin
    iIterTeMix = 0
    
    !Use Newton-Rapson iterations to get a better approximation of Te:
    UDeviation = 2.*ToleranceUeV
    iterations: do while(abs(UDeviation) >ToleranceUeV .and. iIterTeMix<=20)
       
       !Find the populations for the trial Te
       
       call set_ionization_equilibrium_in_mix(TeMix, Na, IsDegenerated) 
       UDeviation = internal_energy_mix()-Uin

       !Calculate the improved value of Te, limiting the iterations so 
       !they can't jump too far out
       
       TeMix = min(2.0*TeMix, max(0.5*TeMix, TeMix - UDeviation/heat_capacity_mix()))  
    
       iIterTeMix = iIterTeMix+1
    end do iterations

  end subroutine set_temperature_in_mix
  
  !============================================
  ! Calculate the pressure in the plasma [Pa]
  ! Can only be called after set_ionization_equilibrium has executed
  
  real function pressure_mix()
    
    pressure_mix = (1+ZavMix)*Na*TeMix*cEV

  end function pressure_mix

  !============================================
  ! calculate the average internal energy per atomic unit [eV]
  ! Can only be called after set_ionization_equilibrium has executed 
  
  real function internal_energy_mix()
    
    internal_energy_mix = 1.50*TeMix*(1+ZAvMix) + EAvMix
 
  end function internal_energy_mix

  !==================================
  !Calculate the specific heat capacity at constant volume 
  !(derivative of internal energy wrt Te) from temperature:
  !Can only be called after set_ionization_equilibrium has executed
  real function heat_capacity_mix()
    real :: TeInv,    &         ! The inverse of the electron temperature [1/eV]
            ETeInvAvJ,&         ! < Ei/Te>_J (Ei - energy levels, Te - electron temperature [eV])
            ZJ,       &         ! Z_J, averaged Z for J component
            DeltaETeInv2Av,&	! <(delta Ei/Te)^2>
   	    ETeInvDeltaZAv      ! <delta i * Ei/Te>

    ! Array of ionization energy levels of ions divided by the temperature in eV
    real,dimension(0:nZMax) :: ETeInv_I
    integer :: iMix
    !------------------

    if( ZAvMix <= cTiny )then
       heat_capacity_mix = 1.50; return
    end if

    !calculate the values of the variables defined above:
 
    TeInv = cOne/TeMix
    DeltaETeInv2Av = 0.0
    ETeInvDeltaZAv = 0.0
    
    do iMix = 1, nMix

       ETeInv_I(iZMin_I(iMix):iZMax_I(iMix)) = IonizEnergyNeutral_II( iZMin_I(iMix):iZMax_I(iMix),iMix )*TeInv

       !Calculate average vaues of E/T and Z, for this component:

       ETeInvAvJ   = sum(Population_II(iZMin_I(iMix):iZMax_I(iMix),iMix )* ETeInv_I(iZMin_I(iMix):iZMax_I(iMix)))
       ZJ          = sum(Population_II(iZMin_I(iMix):iZMax_I(iMix),iMix )* N_I(iZMin_I(iMix):iZMax_I(iMix)))

       !Calculate contributions to bi-linear deviators
       DeltaETeInv2Av   = DeltaETeInv2Av + Concentration_I(iMix)*&
            sum( Population_II(iZMin_I(iMix):iZMax_I(iMix), iMix) * (ETeInv_I(iZMin_I(iMix):iZMax_I(iMix))-ETeInvAvJ)**2 )
       ETeInvDeltaZAv   = ETeInvDeltaZAv + Concentration_I(iMix)*&
            sum( Population_II(iZMin_I(iMix):iZMax_I(iMix), iMix) * ETeInv_I(iZMin_I(iMix):iZMax_I(iMix)) * &
                           (N_I(iZMin_I(iMix):iZMax_I(iMix))-ZJ) )
    end do
    
    ! calculate the heat capacity:
    heat_capacity_mix = 1.50*(cOne +ZAvMix) + DeltaETeInv2Av &
	           +( 3.0*ZAvMix*(0.750*DeltaZ2Av + ETeInvDeltaZAv) &
                   - ETeInvDeltaZAv**2 ) / (ZAvMix + DeltaZ2Av)

  end function heat_capacity_mix !


  !=======================================!
  ! Calculating the Z average values from populations

   real function zMix_averaged()
    integer :: iMix
    zMix_averaged = 0.0
    do iMix = 1,nMix
       zMix_averaged = zMix_averaged + &
            Concentration_I(iMix) * sum( Population_II(iZMin_I(iMix):iZMax_I(iMix) , iMix ) * N_I(iZMin_I(iMix):iZMax_I(iMix) ))
    end do
  end function zMix_averaged

  !==================================
  !Calculate the average ionization energy from neutral atoms of the ions
  
  real function EMix_averaged()
    integer :: iMix
    EMix_averaged = 0.0
    do iMix = 1,nMix
       EMix_averaged = EMix_averaged + Concentration_I(iMix) * sum(Population_II(iZMin_I(iMix):iZMax_I(iMix), iMix) * &
            IonizEnergyNeutral_II(iZMin_I(iMix):iZMax_I(iMix),iMix))
    end do
  end function EMix_averaged				  

end module ModStatSumMix
