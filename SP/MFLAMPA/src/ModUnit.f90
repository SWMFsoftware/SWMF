!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================!
module SP_ModUnit
  ! Unit for particle energy,  energy to momentum conversion
  ! for proton, names for all units
  ! Dec.24 2017 Sokolov & Borovikov.
  use SP_ModGrid, ONLY: nVar, LagrID_, FluxMax_
  use ModConst,   ONLY: energy_in                            , &
       gen_kin_energy_to_momentum=>kinetic_energy_to_momentum, &
       gen_momentum_to_kin_energy=>momentum_to_kinetic_energy, &
       gen_momentum_to_energy    =>momentum_to_energy             
  implicit none
  SAVE
  private !Except
  ! public members
  public :: init               ! Initialize
  public :: read_param         ! Read particle energy unit
  public :: UnitParticleEnergy ! Energy unit is SI
  public :: NameVarUnit_V      ! Units for state vector components
  ! Convert particle momentum to energy or kinetic energy and
  ! kinetic energy to momnetum, for proton
  public :: kinetic_energy_to_momentum, momentum_to_energy, &
       momentum_to_kinetic_energy
  !\
  ! unit of SEP energy is also applicable for ion temperature
  character(len=3)            :: NameEnergyUnit = 'kev'
  real                        :: UnitParticleEnergy ! In SI
  ! simulated particles, used for converting momentum to energy
  ! is back. For different sorts of ions the sqrt(A) factor needs
  ! to be used in these relations.
  character(len=*), parameter :: NameParticle = 'proton'
  !/
  character(len=6)            :: NameVarUnit_V(LagrID_:FluxMax_) = (/&
       'none  ', &
       'RSun  ', &
       'RSun  ', &
       'RSun  ', &
       'amu/m3', &
       'kev   ', &
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
       '??????'  /) ! TBD
  logical :: DoInit = .true.
contains
   subroutine read_param(NameCommand)
    use ModReadParam, ONLY: read_var
    use ModUtilities, ONLY: lower_case
    use SP_ModGrid  , ONLY: T_
    character(len=*), intent(in):: NameCommand ! From PARAM.in  
    character(len=*), parameter :: NameSub='SP:read_param_unit'
    !----------------------------------------------------------
    select case(NameCommand)
    case('#PARTICLEENERGYUNIT')
       !Read unit to be used for particle energy: eV, keV, GeV
       call read_var('ParticleEnergyUnit',NameEnergyUnit)
       call lower_case(NameEnergyUnit)
       NameVarUnit_V(T_) = NameEnergyUnit//'   '
    case default
       call CON_stop(NameSub//'Unknown command '//NameCommand)
    end select
  end subroutine read_param
  !==============
  subroutine init
    character(len=*), parameter :: NameSub='SP:init_unit'
    !------------------
    if(.not.DoInit)RETURN
    DoInit = .false.
     ! account for units of energy
    UnitParticleEnergy = energy_in(NameEnergyUnit)
  end subroutine init
  !==============
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
  !==============================
  real function momentum_to_kinetic_energy(Momentum)
    real, intent(in) :: Momentum
    !---------------
    momentum_to_kinetic_energy = &
         gen_momentum_to_kin_energy(Momentum, NameParticle)
  end function momentum_to_kinetic_energy
end module SP_ModUnit
