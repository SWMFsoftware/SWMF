module ModHeidiInput

  use ModReadParam, ONLY: read_file, read_init, i_session_read, &
       read_line, read_command, read_var
  use ModUtilities, ONLY: split_string

  implicit none

  character(len=*), parameter:: NameMod = 'ModHeidiInput'
  ! GIPHT BEGIN DECLARATIONS

  ! >>> Testing <<<
  ! "#STARTTIME", "#SETREALTIME"
  integer :: iYear = 2002
  integer :: iMonth = 4
  integer :: iDay = 17
  integer :: iHour = 0
  integer :: iMinute = 0
  integer :: iSecond = 0
  ! "#STOP"
  real :: tSimulationMax = 0.0
  ! "#TIMESTEP"
  real :: TimeStep = 20.0
  ! "#GRID"
  integer :: nRadialGrid = 20
  integer :: nPhiGrid = 24
  integer :: nEnergyGrid = 42
  integer :: nPitchGrid = 71
  integer :: nPointFieldLine = 24
  ! "#STORM"
  character(len=100) :: TypeStorm
  ! "#INNERBOUNDARY"
  real :: Height = 1e6
  ! "#ENERGYSETUP"
  real :: EnergyLowerBoundary
  real :: LowestEnergyCellWidth
  real :: GrowthMultiplier
  ! "#SPECIES"
  logical :: UseElectron = .true.
  logical :: UseHPlus = .true.
  logical :: UseHePlus = .true.
  logical :: UseOPlus = .true.
  logical :: UseIon = .true.
  logical :: UsePhotoelElectron = .false.
  logical :: UsePlasmaSheetElectron = .false.
  ! "#INDICES"
  integer :: WhichKp
  real :: KpIndex
  real :: ApIndex
  real :: SunspotAverage
  ! "#BOUNDARY"
  character(len=100) :: TypeBoundary
  ! "#INITIAL"
  character(len=100) :: TypeInitial
  real :: MaxwellianScallingFactor
  real :: CharacteristicEnergy
  ! "#OUTPUT"
  logical :: DoSaveDistributionFunctionEverywhere = .true.
  logical :: DoSaveEquatorialDistributionFunction = .false.
  logical :: DoSaveEnergyDeposition = .true.
  logical :: DoSaveTotalPrecipitationFlux = .false.
  logical :: DoSaveDifferentialPrecipitationFlux = .false.
  logical :: DoSaveParticleEnergyLosses = .false.
  logical :: DoSaveThermalPlasmaDensity = .true.
  logical :: DoSaveCflForAdvection = .false.
  logical :: DoSaveDriftVelocities = .true.
  logical :: DoSaveEvsLDistributions = .false.
  logical :: DoSaveParticleLifetimes = .false.
  logical :: DoSavePressureDensityDst = .true.
  logical :: DoSaveUnformatted = .true.
  logical :: DoSaveContinuousSourcesLosses = .true.
  logical :: DoSaveNightsideBCDistribution = .true.
  logical :: DoSaveDifferentialNumberFlux = .false.
  ! "#OUTPUTINFO"
  character(len=20) :: NameRunIn = 'test1'
  real :: nFrequency
  ! "#INJECTIONFREQUENCY"
  real :: iFrequency
  ! "#CONVECTION"
  character(len=100) :: TypeConvection
  ! "#INITIALTHERMALPLASMA"
  logical :: DoReadDGCPM = .true.
  ! "#SOLARWIND"
  logical :: DoReadSolarWind = .true.
  ! "#PITCHANGLE"
  logical :: UseConstantStepPA = .false.
  ! "#INCLUDEWAVES"
  logical :: UseWaves = .false.
  ! "#BFIELD"
  character(len=20) :: TypeBField = 'analytic'
  ! "#SAVERESTART"
  logical           :: DoSaveRestart = .true.
  real              :: DtSaveRestart = 40.0
  character(len=20) :: TypeFile = 'ascii'
 


  ! GIPHT END DECLARATIONS

contains

  subroutine set_parameters

    use ModProcIM, ONLY: iComm
    use ModHeidiIO, ONLY: IsFramework

    character (len=100) :: NameCommand, StringPart_I(100)
    integer :: iSession = -1, nStringPart
    logical :: UseStrict
    ! GIPHT BEGIN INDEXES

    ! GIPHT END INDEXES
    character(len=*), parameter:: NameSub = NameMod//'::set_parameters'
    !-------------------------------------------------------------------------
    if(.not.IsFramework .and. iSession == -1)then
       call read_file('PARAM.in', iComm)
       call read_init('  ',iSessionIn=1, iLineIn=0)
    end if

    iSession = i_session_read()

    write(*,*) NameSub,' starting for iSession=',iSession

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       ! GIPHT BEGIN COMMANDS
       select case(NameCommand)

          ! >>> Testing <<<

       case("#STARTTIME", "#SETREALTIME")
          if(.not.is_first_session() .or. IsFramework)CYCLE
          call read_var('iYear', iYear)
          call read_var('iMonth', iMonth)
          call read_var('iDay', iDay)
          call read_var('iHour', iHour)
          call read_var('iMinute', iMinute)
          call read_var('iSecond', iSecond)
       case("#STOP")
          if(.not.is_first_session() .or. IsFramework)CYCLE
          call read_var('tSimulationMax', tSimulationMax)
       case("#TIMESTEP")
          call read_var('dtmax', TimeStep)
       case("#GRID")
          call read_var('nRadialGrid', nRadialGrid)
          call read_var('nPhiGrid', nPhiGrid)
          call read_var('nEnergyGrid', nEnergyGrid)
          call read_var('nPitchGrid', nPitchGrid)
          call read_var('nPointFieldLine', nPointFieldLine)
       case("#STORM")
          call read_var('TypeStorm', TypeStorm)
       case("#INNERBOUNDARY")
          call read_var('Height', Height)
       case("#ENERGYSETUP")
          call read_var('EnergyLowerBoundary', EnergyLowerBoundary)
          call read_var('LowestEnergyCellWidth', LowestEnergyCellWidth)
          call read_var('GrowthMultiplier', GrowthMultiplier)
       case("#SPECIES")
          call read_var('UseElectron', UseElectron)
          call read_var('UseHPlus', UseHPlus)
          call read_var('UseHePlus', UseHePlus)
          call read_var('UseOPlus', UseOPlus)
          call read_var('UseIon', UseIon)
          call read_var('UsePhotoelElectron', UsePhotoelElectron)
          call read_var('UsePlasmaSheetElectron', UsePlasmaSheetElectron)
       case("#INDICES")
          call read_var('WhichKp', WhichKp)
          call read_var('KpIndex', KpIndex)
          call read_var('ApIndex', ApIndex)
          call read_var('SunspotAverage', SunspotAverage)
       case("#BOUNDARY")
          call read_var('TypeBoundary', TypeBoundary, IsLowerCase=.true.)
       case("#INITIAL")
          call read_var('TypeInitial', TypeInitial, IsLowerCase=.true.)
          call read_var('MaxwellianScallingFactor', MaxwellianScallingFactor)
          call read_var('CharacteristicEnergy', CharacteristicEnergy)
       case("#OUTPUT")
          call read_var('DoSaveDistributionFunctionEverywhere', DoSaveDistributionFunctionEverywhere)
          call read_var('DoSaveEquatorialDistributionFunction', DoSaveEquatorialDistributionFunction)
          call read_var('DoSaveEnergyDeposition', DoSaveEnergyDeposition)
          call read_var('DoSaveTotalPrecipitationFlux', DoSaveTotalPrecipitationFlux)
          call read_var('DoSaveDifferentialPrecipitationFlux', DoSaveDifferentialPrecipitationFlux)
          call read_var('DoSaveParticleEnergyLosses', DoSaveParticleEnergyLosses)
          call read_var('DoSaveThermalPlasmaDensity', DoSaveThermalPlasmaDensity)
          call read_var('DoSaveCflForAdvection', DoSaveCflForAdvection)
          call read_var('DoSaveDriftVelocities', DoSaveDriftVelocities)
          call read_var('DoSaveEvsLDistributions', DoSaveEvsLDistributions)
          call read_var('DoSaveParticleLifetimes', DoSaveParticleLifetimes)
          call read_var('DoSavePressureDensityDst', DoSavePressureDensityDst)
          call read_var('DoSaveUnformatted', DoSaveUnformatted)
          call read_var('DoSaveContinuousSourcesLosses', DoSaveContinuousSourcesLosses)
          call read_var('DoSaveNightsideBCDistribution', DoSaveNightsideBCDistribution)
          call read_var('DoSaveDifferentialNumberFlux', DoSaveDifferentialNumberFlux)
       case("#OUTPUTINFO")
          call read_var('NameRun', NameRunIn)
          call read_var('nFrequency', nFrequency)
       case("#INJECTIONFREQUENCY")
          call read_var('iFrequency', iFrequency)
       case("#CONVECTION")
          call read_var('TypeConvection', TypeConvection, IsLowerCase=.true.)
       case("#BFIELD")
          call read_var('TypeBField', TypeBField, IsLowerCase=.true.)   
       case("#INITIALTHERMALPLASMA")
          call read_var('DoReadDGCPM', DoReadDGCPM)
       case("#SOLARWIND")
          call read_var('DoReadSolarWind', DoReadSolarWind)
       case("#PITCHANGLE")
          call read_var('UseConstantStepPA', UseConstantStepPA)
       case("#INCLUDEWAVES")
          call read_var('UseWaves', UseWaves)
       case("#SAVERESTART")
          call read_var('DoSaveRestart',DoSaveRestart)
          call read_var('DtSaveRestart',DtSaveRestart)
          call read_var('TypeFile',TypeFile)


          ! GIPHT END COMMANDS
       case default
          !if(iProc==0) then
          write(*,*) NameSub // ' WARNING: unknown command ' // &
               trim(NameCommand),' !',DtSaveRestart
          if(UseStrict)call CON_stop('Correct PARAM.in!')
          !end if
       end select
    end do

    call set_heidi_variables

  contains
    !==========================================================================
    subroutine set_heidi_variables
      use ModHeidiIO
      use ModHeidiSize
      use ModHeidiMain, ONLY: itherminit

      implicit none
      integer :: iSpecies
      !-----------------------------------------
      year                  = iYear 
      month                 = iMonth
      day                   = iDay
      ut                    = iHour 
      iMinute               = 0.0
      iSecond               = 0.0

      tmax                  = tSimulationMax
      dtmax                 = TimeStep 
      io                    = nRadialGrid 
      jo                    = nPhiGrid
      ko                    = nEnergyGrid
      lo                    = nPitchGrid 
      iso                   = nPointFieldLine

      hmin                  = Height  
      elb                   = EnergyLowerBoundary
      swe                   = LowestEnergyCellWidth
      rw                    = GrowthMultiplier

      ikp                   = WhichKp 
      kp                    = KpIndex  
      Ap                    = ApIndex   
      r                     = SunspotAverage  

      tinj                  = iFrequency
      tint                  = nFrequency
      NameRun               = NameRunIn 

      Ab                    = MaxwellianScallingFactor
      Eob                   = CharacteristicEnergy

      !Convert from strings to integers
      if (TypeStorm == ('major'))         istorm = 1
      if (TypeStorm == ('moderate'))      istorm = 2
      if (TypeStorm == ('test'))          istorm = 3
     

      !Convert from logicals to integers 
      if (DoSaveDistributionFunctionEverywhere)         ires(1) = 1
      if (DoSaveEquatorialDistributionFunction)         ires(2) = 1 
      if (DoSaveEnergyDeposition)                       ires(3) = 1
      if (DoSaveTotalPrecipitationFlux)                 ires(4) = 1
      if (DoSaveDifferentialPrecipitationFlux)          ires(5) = 1
      if (DoSaveParticleEnergyLosses)                   ires(6) = 1
      if (DoSaveThermalPlasmaDensity)                   ires(7) = 1
      if (DoSaveCflForAdvection)                        ires(8) = 1
      if (DoSaveDriftVelocities)                        ires(9) = 1
      if (DoSaveEvsLDistributions)                      ires(10) = 1
      if (DoSaveParticleLifetimes)                      ires(11) = 1
      if (DoSavePressureDensityDst)                     ires(12) = 1
      if (DoSaveUnformatted)                            ires(13) = 1
      if (DoSaveContinuousSourcesLosses)                ires(14) = 1
      if (DoSaveNightsideBCDistribution)                ires(15) = 1
      if (.not. DoSaveDifferentialNumberFlux)           ifac = 0
      if (DOSaveDifferentialNumberFlux)                 ifac = 1
      

      !Set Convection 
      select case(TypeConvection)
      case('0',  'none')
         iA = 0
      case('1',  'kpvsmaynardchen')
         IA = 1  
      case ('2', 'mbivs')
         IA = 2
      case ('3', 'mcilwain')
         IA = 3
      case ('4', 'mcilwainpluscpcp')
         IA = 4
      case ('5','kpvsplusburkewygant')
         IA = 5
      case ('6','mcilwainplusburkewygant')
         IA = 6
      case ('7','mcilwaincpcpburkewygant')
         IA = 7
      case ('8','kpvsselfconsintent')
         IA = 8
      case ('9','mcilwainselfconsistent')
         IA = 9
      case ('10','mcilwaincpcpselfconsistent')
         IA = 10
      case ('11','unshieldedvswithpen')
         IA = 11
      case ('12','unshieldedvsnopen')
         IA = 12
      case ('13','w96sc')
         IA = 13
      case ('14','w96')
         IA = 14
      case ('15','amiesc')
         IA = 15 
      case ('16','amie')
         IA = 16
      case ('20','readfromfile')
         IA = 20
      case ('21','amiepot')
         IA = 21
      case ('22','w96pot')
         IA = 22
      case ('23','foster')
         IA = 23
      end select
      
      if (DoReadSolarWind)         isw  = 1
      if (UseWaves)                iwpi = 1
      if (UseConstantStepPA)       ipa  = 0
      if (.not. UseConstantStepPA) ipa  = 1
          
      if (DoReadDGCPM)            Itherminit = 1
      if (UseIon)                 ist = 0
      if (UsePhotoelElectron)     ist = 1
      if (UsePlasmaSheetElectron) ist = 2
             
      if (UseElectron)          SCALC(1) = 1
      if (UseHPlus)             SCALC(2) = 1 
      if (UseHePlus)            SCALC(3) = 1
      if (UseOPlus)             SCALC(4) = 1

      ! Set initial conditions
      select case(TypeInitial)
      case('0',  'none')
         do iSpecies = 1, 4
            ini(iSpecies) = 0
         end do
      case('1',  'maxwellian')
         do iSpecies = 1, 4
            ini(iSpecies) = 1
         end do
      case('2',  'gaussian')
         do iSpecies = 1, 4
            ini(iSpecies) = 2
         end do
      case('3',  'frominputfile')
         do iSpecies = 1, 4
            ini(iSpecies) = 3
         end do
      case('4',  'quietrc')
         do iSpecies = 1, 4
            ini(iSpecies) = 4
         end do
      case('5',  'fromfile')
         do iSpecies = 1, 4
            ini(iSpecies) = 5
         end do
      case('6',  'psinject')
         do iSpecies = 1, 4
            ini(iSpecies) = 6
         end do
      case('7',  'fromrestart')
         do iSpecies = 1, 4
            ini(iSpecies) = 7
         end do
      end select
  
        ! Set boundary conditions
      select case(TypeBoundary)
      case('0',  'none')
         do iSpecies = 1, 4
            ibc(iSpecies) = 0
         end do
      case('1',  'maxwellian')
         do iSpecies = 1, 4
            ibc(iSpecies) = 1
         end do
      case('2',  'gaussian')
         do iSpecies = 1, 4
            ibc(iSpecies) = 2
         end do
      case('3',  'frominputfile')
         do iSpecies = 1, 4
            ibc(iSpecies) = 3
         end do
      case('4',  'quietrc')
         do iSpecies = 1, 4
            ibc(iSpecies) = 4
         end do
      case('5',  'fromfile')
         do iSpecies = 1, 4
            ibc(iSpecies) = 5
         end do
      case('6',  'psinject')
         do iSpecies = 1, 4
            ibc(iSpecies) = 6
         end do
      case('7',  'fromrestart')
         do iSpecies = 1, 4
            ibc(iSpecies) = 7
         end do
      end select
      

    end subroutine set_heidi_variables

    !==========================================================================
    logical function is_first_session()
      is_first_session = iSession == 1
      if(iSession > 1)then
         ! if(iProc==0) then
         write(*,*) NameSub // ' WARNING: command ',trim(NameCommand), &
              ' can be used in first session only!'
         if(UseStrict)call CON_stop('Correct PARAM.in!')
         ! end if
      end if
    end function is_first_session

  end subroutine set_parameters

end module ModHeidiInput
