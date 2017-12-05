!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModExtraVariables

  ! Define indexes that are not used in most equation modules

  ! Set impossible value for indexes
  integer, parameter :: Ppar_ = 1, iPparIon_I(1) = 1
  integer, parameter :: Pe_ = 1, Pepar_ = 1
  integer, parameter :: Hyp_ = 1
  integer, parameter :: Erad_ = 1, ExtraEint_ = 1
  integer, parameter :: Ew_  = 1
  integer, parameter :: SignB_ = 1
  integer, parameter :: Te0_   = 1
  integer, parameter :: Ehot_ = 1

  ! Electric field and hyperbolic scalar field
  integer, parameter :: Ex_ = 1, Ey_ = 2, Ez_ = 3
  integer, parameter :: HypE_ = 1

  ! Number of electron fluids
  integer, parameter :: nElectronFluid = 0

  ! The named index range for frequency bins in multi-group
  integer, parameter :: nWave = 1
  integer, parameter :: WaveFirst_ = 1, WaveLast_ = 1

  ! The named index range for material levels
  integer, parameter :: nMaterial = 1
  integer, parameter :: MaterialFirst_ = 1, MaterialLast_ = 1

  ! The named index range for species and their masses
  integer, parameter :: SpeciesFirst_ = 1, SpeciesLast_ = 1
  real               :: MassSpecies_V(SpeciesFirst_:SpeciesLast_)

end module ModExtraVariables
!==============================================================================
