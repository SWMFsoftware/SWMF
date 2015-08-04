!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModExtraVariables
  
  ! Define indexes that are not used in most equation modules

  ! Set impossible value for indexes
  integer, parameter :: Pe_ = 1, Ppar_ = 1, iPparIon_I(1) = 1
  integer, parameter :: Hyp_ = 1, Erad_ = 1, ExtraEint_ = 1
  integer, parameter :: Ew_  = 1
  integer, parameter :: SignB_ = 1
  integer, parameter :: Te0_   = 1
  integer, parameter :: Ehot_ = 1

  ! The named index range for frequency bins in multi-group
  integer, parameter :: nWave = 1
  integer, parameter :: WaveFirst_ = 1, WaveLast_ = 1

  ! The named index range for material levels
  integer, parameter :: nMaterial = 1
  integer, parameter :: MaterialFirst_ = 1, MaterialLast_ = 1

end module ModExtraVariables
