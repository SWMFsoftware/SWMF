!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================!
module SP_ModUnit
  use SP_ModGrid, ONLY: nVar, LagrID_, FluxMax_
  use ModConst,   ONLY: &
       gen_kin_energy_to_momentum=>kinetic_energy_to_momentum, &
       gen_momentum_to_kin_energy=>momentum_to_kinetic_energy, &
       gen_momentum_to_energy    =>momentum_to_energy             
  implicit none
  SAVE
  private !Except
  public :: kinetic_energy_to_momentum, momentum_to_energy, &
       momentum_to_kinetic_energy, NameEUnit
  !\
  ! unit of SEP energy is also applicable for ion temperature
  character(len=*), parameter :: NameEUnit = 'kev'
  real, public                :: UnitEnergy
  ! simulated particles, used for converting momentum to energy
  ! is back. For different sorts of ions the sqrt(A) factor needs
  ! to be used in these relations.
  character(len=*), parameter :: NameParticle = 'proton'
  !/
  character(len=6), public, parameter:: NameVarUnit_V(LagrID_:FluxMax_) = (/&
       'none  ', &
       'RSun  ', &
       'RSun  ', &
       'RSun  ', &
       'amu/m3', &
       NameEUnit//'   ', &
       'm/s   ', &
       'm/s   ', &
       'm/s   ', &
       'T     ', &
       'T     ', &
       'T     ', &
       'J/m3  ', &
       'J/m3  ', &
       'RSun  ', &
       'RSun  ', &
       'RSun  ', &
       'm/s   ', &
       'T     ', &
       'none  ', &
       'amu/m3', &
       'T     ', &
       'p.f.u.', &
       'p.f.u.', &
       'p.f.u.', &
       'p.f.u.', &
       'p.f.u.', &
       'p.f.u.', &
       'p.f.u.', &
       '??????'  /)
contains
  real function kinetic_energy_to_momentum(Energy)
    real, intent(in) :: Energy
    !---------------
    kinetic_energy_to_momentum = &
         gen_kin_energy_to_momentum(Energy, NameParticle)
  end function kinetic_energy_to_momentum
  !======================================
  real function momentum_to_energy(Momentum)
    real, intent(in) :: Momentum
    !---------------
    momentum_to_energy = &
         gen_momentum_to_energy(Momentum, NameParticle)
  end function momentum_to_energy
  !======================================
  real function momentum_to_kinetic_energy(Momentum)
    real, intent(in) :: Momentum
    !---------------
    momentum_to_kinetic_energy = &
         gen_momentum_to_kin_energy(Momentum, NameParticle)
  end function momentum_to_kinetic_energy
end module SP_ModUnit
