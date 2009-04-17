module ModHeidiInput

  use ModReadParam, ONLY: i_session_read, read_line, read_command, read_var
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
  ! "#TIMESIMULATION"
  real :: tSimulation = 0.0
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
  integer :: BCElectron = 0
  integer :: BCHPlus = 0
  integer :: BCHePlus = 0
  integer :: BCOPlus = 0
  ! "#INITIAL"
  integer :: InitialElectron = 0
  integer :: InitialHPlus = 0
  integer :: InitialHePlus = 0
  integer :: InitialOPlus = 0
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
  character(len=10) :: NameRun
  real :: nFrequency
  ! "#INJECTIONFREQUENCY"
  real :: iFrequency
  ! "#CONVECTION"
  logical :: IsNoConvection = .false.
  logical :: IsKpVSMaynardChen = .false.
  logical :: IsMBIVS = .false.
  logical :: IsMcIlwain = .false.
  logical :: IsMcIlwainPlusCPCP = .false.
  logical :: IsKpVSPlusBurkeWygant = .false.
  logical :: IsMcIlwainPlusBurkeWygant = .false.
  logical :: IsMcIlwainCPCPBurkeWygant = .false.
  logical :: IsKpVSSelfConsintent = .false.
  logical :: IsMcIlwainSelfConsistent = .false.
  logical :: IsMcIlwainCPCPSelfConsistent = .false.
  logical :: IsUnshieldedVSwithPen = .false.
  logical :: IsUnshieldedVSnoPen = .false.
  logical :: IsW96SC = .false.
  logical :: IsW96 = .true.
  logical :: IsAMIESC = .false.
  logical :: IsAMIE = .false.
  logical :: IsReadFromFile = .false.
  logical :: IsAMIEPot = .false.
  logical :: IsW96Pot = .false.
  logical :: IsFoster = .false.
  ! "#INITIALTHERMALPLASMA"
  logical :: DoReadDGCPM = .true.
  ! "#SOLARWIND"
  logical :: DoReadSolarWind = .true.
  ! "#PITCHANGLE"
  logical :: UseConstantStepPA = .false.
  ! "#INCLUDEWAVES"
  logical :: UseWaves = .false.
  ! GIPHT END DECLARATIONS

contains

  subroutine set_parameters

    character (len=100) :: NameCommand, StringPart_I(100)
    integer :: iSession, nStringPart
    logical :: UseStrict
    ! GIPHT BEGIN INDEXES

    ! GIPHT END INDEXES
    character(len=*), parameter:: NameSub = NameMod//'::set_parameters'
    !-------------------------------------------------------------------------
    iSession = i_session_read()

    write(*,*) NameSub,' starting for iSession=',iSession

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       ! GIPHT BEGIN COMMANDS
       select case(NameCommand)

          ! >>> Testing <<<

       case("#STARTTIME", "#SETREALTIME")
          if(.not.is_first_session())CYCLE
          call read_var('iYear', iYear)
          call read_var('iMonth', iMonth)
          call read_var('iDay', iDay)
          call read_var('iHour', iHour)
          call read_var('iMinute', iMinute)
          call read_var('iSecond', iSecond)
       case("#TIMESIMULATION")
          if(.not.is_first_session())CYCLE
          call read_var('tSimulation', tSimulation)
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
          call read_var('BCElectron', BCElectron)
          call read_var('BCHPlus', BCHPlus)
          call read_var('BCHePlus', BCHePlus)
          call read_var('BCOPlus', BCOPlus)
       case("#INITIAL")
          call read_var('InitialElectron', InitialElectron)
          call read_var('InitialHPlus', InitialHPlus)
          call read_var('InitialHePlus', InitialHePlus)
          call read_var('InitialOPlus', InitialOPlus)
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
          call read_var('NameRun', NameRun)
          call read_var('nFrequency', nFrequency)
       case("#INJECTIONFREQUENCY")
          call read_var('iFrequency', iFrequency)
       case("#CONVECTION")
          call read_var('IsNoConvection', IsNoConvection)
          call read_var('IsKpVSMaynardChen', IsKpVSMaynardChen)
          call read_var('IsMBIVS', IsMBIVS)
          call read_var('IsMcIlwain', IsMcIlwain)
          call read_var('IsMcIlwainPlusCPCP', IsMcIlwainPlusCPCP)
          call read_var('IsKpVSPlusBurkeWygant', IsKpVSPlusBurkeWygant)
          call read_var('IsMcIlwainPlusBurkeWygant', IsMcIlwainPlusBurkeWygant)
          call read_var('IsMcIlwainCPCPBurkeWygant', IsMcIlwainCPCPBurkeWygant)
          call read_var('IsKpVSSelfConsintent', IsKpVSSelfConsintent)
          call read_var('IsMcIlwainSelfConsistent', IsMcIlwainSelfConsistent)
          call read_var('IsMcIlwainCPCPSelfConsistent', IsMcIlwainCPCPSelfConsistent)
          call read_var('IsUnshieldedVSwithPen', IsUnshieldedVSwithPen)
          call read_var('IsUnshieldedVSnoPen', IsUnshieldedVSnoPen)
          call read_var('IsW96SC', IsW96SC)
          call read_var('IsW96', IsW96)
          call read_var('IsAMIESC', IsAMIESC)
          call read_var('IsAMIE', IsAMIE)
          call read_var('IsReadFromFile', IsReadFromFile)
          call read_var('IsAMIEPot', IsAMIEPot)
          call read_var('IsW96Pot', IsW96Pot)
          call read_var('IsFoster', IsFoster)
       case("#INITIALTHERMALPLASMA")
          call read_var('DoReadDGCPM', DoReadDGCPM)
       case("#SOLARWIND")
          call read_var('DoReadSolarWind', DoReadSolarWind)
       case("#PITCHANGLE")
          call read_var('UseConstantStepPA', UseConstantStepPA)
       case("#INCLUDEWAVES")
          call read_var('UseWaves', UseWaves)


          ! GIPHT END COMMANDS
       case default
          !if(iProc==0) then
          write(*,*) NameSub // ' WARNING: unknown command ' // &
               trim(NameCommand),' !'
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

      tmax                  = tSimulation
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
      name                  = NameRun 

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
      if (IsNoConvection)               IA = 0
      if (IsKpVSMaynardChen)            IA = 1 
      if (IsMBIVS)                      IA = 2
      if (IsMcIlwain)                   IA = 3
      if (IsMcIlwainPlusCPCP )          IA = 4
      if (IsKpVSPlusBurkeWygant)        IA = 5
      if (IsMcIlwainPlusBurkeWygant)    IA = 6
      if (IsMcIlwainCPCPBurkeWygant)    IA = 7
      if (IsKpVSSelfConsintent)         IA = 8
      if (IsMcIlwainSelfConsistent)     IA = 9
      if (IsMcIlwainCPCPSelfConsistent) IA = 10
      if (IsUnshieldedVSwithPen)        IA = 11
      if (IsUnshieldedVSnoPen )         IA = 12
      if (IsW96SC)                      IA = 13
      if (IsW96)                        IA = 14
      if (IsAMIESC)                     IA = 15
      if (IsAMIE)                       IA = 16
      if (IsReadFromFile)               IA = 20
      if (IsAMIEPot)                    IA = 21
      if (IsW96Pot)                     IA = 22
      if (IsFoster)                     IA = 23
      
      
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

      ini(1) = InitialElectron 
      ini(2) = InitialHPlus 
      ini(3) = InitialHePlus
      ini(4) = InitialOPlus

      ibc(1) = BCElectron 
      ibc(2) = BCHPlus
      ibc(3) = BCHePlus
      ibc(4) = BCOPlus


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
