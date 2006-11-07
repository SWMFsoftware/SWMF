!^CFG COPYRIGHT UM
Module ModConst
  use ModNumConst
  use ModKind

  implicit none

  save

  !\  
  ! Physical and solar astronomical constants.
  !
  ! All constants for planets, satellites, comets and
  ! other astronomical bodyies are found in ModPlanetConst
  !/
  !\
  ! Physical constants
  !/

  real,parameter::  cBoltzmann  = 13.807/cE24
  ! Boltzmann constant  1.3807E-23 [J/K]

  real, parameter :: cProtonMass = 1.6726/cE27
  ! Proton mass  1.6726E-27[kg]

  real, parameter :: cElectronMass = 9.1094/cE30/cE1
  ! Electron mass  9.1094E-31[kg]

  real, parameter :: cLightSpeed =  0.29979*cE9
  ! Speed of light (true) 2.9979E+08 [m/s]

  real, parameter :: cMu = cPi*0.4/cE6
  ! Vacuum permeability Pi*4.0E-7[H/m]

  real, parameter :: cEps = cOne/cMu/(cLightSpeed**2)
  ! Vacuum permittivity 8.8542E-12[F/m]
  
  real,parameter :: cVacuumResistivity=cLightSpeed*cMu
  ! Vacuum resistance 120*cPi [Ohm]

  real, parameter :: cElectronCharge  = 0.16022/cE18
  ! Fundamental charge 1.6022E-19 [C (coulomb)]

  real, parameter :: cAvogadro =0.6022045*cE24
  ! Number of particles per mole 6.022045E+23[mole^-1]

  real, parameter :: cGravitation = 66.726/cE12
  ! Gravitation constant 6.6726E-11  (NRL 1994)

  !\
  ! Units for energy.
  !/
  real, parameter :: cEV  = cElectronCharge
  real, parameter :: cKEV = cE3 * cEV
  real, parameter :: cMEV = cE3 * cKEV
  real, parameter :: cGEV = cE3 * cMEV
  real, parameter :: cTEV = cE3 * cGEV

  real, parameter :: cEVToK  =  cEV / cBoltzmann
  real, parameter :: cKEVToK = cKEV / cBoltzmann
  real, parameter :: cMEVToK = cMEV / cBoltzmann

  !\
  ! Here RME stands for Rest Mass Energy.
  !/
  real, parameter :: cRMEProton   = cProtonMass   * cLightSpeed**2
  real, parameter :: cRMEElectron = cElectronMass * cLightSpeed**2

  !\
  ! Non-relativistic formulae for gyrofrequencies:
  ! Gyrofrequency = cGyroParticle * |B|
  !/
  real, parameter :: cGyroProton   = cElectronCharge / cProtonMass
  real, parameter :: cGyroElectron = cElectronCharge / cElectronMass
  !\
  ! Relativistic formulae for gyrofrequencies:
  ! Gyrofrequency = cGyroRel * |B| / Energy
  !/
  real, parameter :: cGyroRel = cElectronCharge * cLightSpeed**2
  !\
  ! Formula for gyroradius:
  ! Gyroradius = cGyroRadius * momentum / |B|
  !/
  real, parameter :: cGyroRadius = cOne / cElectronCharge



  real(Real8_), parameter :: cSecondPerYear   = 31536000.0
  real(Real8_), parameter :: cSecondPerDay    =    86400.0
  real(Real8_), parameter :: cSecondPerHour   =     3600.0
  real(Real8_), parameter :: cSecondPerMinute =       60.0

  !\
  ! Solar Astronomical constants
  !/

  real,parameter :: cAU = 1.4959787E11

  real,parameter :: rSun              = 0.696E9                ! [ m]
  real,parameter :: mSun              = 1.99E30                ! [kg]
  real,parameter :: RotationPeriodSun = 25.38 * 24.0 * 3600.0  ! [ s]



end module ModConst
!====================================================================
!====================================================================
real function momentum_to_energy(Momentum,NameParticle)
  use ModConst
  implicit none
  real,intent(in):: Momentum
  character(LEN=*),intent(in):: NameParticle
  select case(NameParticle)
  case('e','Electron','electron','ELECTRON')
     momentum_to_energy=sqrt((Momentum*cLightSpeed)**2+&
          cRMEElectron**2)
  case('p','Proton','proton','PROTON')
     momentum_to_energy=sqrt((Momentum*cLightSpeed)**2+&
          cRMEProton**2)
  case default
     call CON_stop(&
          'Do not know the rest mass energy for '//NameParticle)
  end select
end function momentum_to_energy
!====================================================================
real function momentum_to_kinetic_energy(Momentum,NameParticle)
  use ModConst
  implicit none
  real,intent(in):: Momentum
  character(LEN=*),intent(in):: NameParticle
  select case(NameParticle)
  case('e','Electron','electron','ELECTRON')
     momentum_to_kinetic_energy=(Momentum*cLightSpeed)**2/(&
          sqrt((Momentum*cLightSpeed)**2+cRMEElectron**2) +&
          cRMEElectron)
  case('p','Proton','proton','PROTON')
     momentum_to_kinetic_energy=(Momentum*cLightSpeed)**2/(&
          sqrt((Momentum*cLightSpeed)**2+cRMEProton**2)   +&
          cRMEProton)
  case default
     call CON_stop(&
          'Do not know the rest mass energy for '//NameParticle)
  end select
end function momentum_to_kinetic_energy
!====================================================================
real function energy_to_momentum(Energy,NameParticle)
  use ModConst
  implicit none
  real,intent(in):: Energy
  character(LEN=*),intent(in):: NameParticle
  select case(NameParticle)
  case('e','Electron','electron','ELECTRON')
     energy_to_momentum=sqrt(Energy**2-cRMEElectron**2)/&
          cLightSpeed
  case('p','Proton','proton','PROTON')
     energy_to_momentum=sqrt(Energy**2 - cRMEProton**2)/&
          cLightSpeed
  case default
     call CON_stop(&
          'Do not know the rest mass energy for '//NameParticle)
  end select
end function energy_to_momentum
!====================================================================
real function kinetic_energy_to_momentum(Energy,NameParticle)
  use ModConst
  implicit none
  real,intent(in):: Energy
  character(LEN=*),intent(in):: NameParticle
  select case(NameParticle)
  case('e','Electron','electron','ELECTRON')
     kinetic_energy_to_momentum=sqrt(&
          Energy*(Energy+cTwo*cRMEElectron))/cLightSpeed
  case('p','Proton','proton','PROTON')
     kinetic_energy_to_momentum=sqrt(&
          Energy*(Energy+cTwo*cRMEProton  ))/cLightSpeed
  case default
     call CON_stop(&
          'Do not know the rest mass energy for '//NameParticle)
  end select
end function kinetic_energy_to_momentum
!====================================================================
real function energy_in(NameEnergyUnit)
  use ModConst
  implicit none
  character(LEN=*),intent(in):: NameEnergyUnit
  select case (NameEnergyUnit)
  case('J','j','Joule','joule','JOULE')
     energy_in=cOne   ! Do Nothing
  case('K','k','Kelvin','kelvin','KELVIN')
     energy_in=cBoltzmann
  case('ev','eV','EV','Ev')
     energy_in=cEV
  case('kEV','KeV','KEV','kev')
     energy_in=cKEV
  case('mEV','MeV','MEV','mev')
     energy_in=cMEV
  case('gEV','GeV','GEV','gev')
     energy_in=cGEV
  case('tEV','TeV','TEV','tev')
     energy_in=cTEV
  case default
     call CON_stop(&
          'We do not support energy units like: '//&
          'ton of trotil equivalent, '//&
          NameEnergyUnit//' and many other.')
  end select
end function energy_in
!====================================================================

