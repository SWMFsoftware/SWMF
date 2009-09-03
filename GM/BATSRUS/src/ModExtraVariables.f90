module ModExtraVariables
  
  ! Define indexes that are not used in most equation modules

  ! Set impossible value for indexes
  integer, parameter :: Pe_ = 1, Ppar_ = 1, Hyp_ = 1, Erad_ = 1
  integer, parameter :: ExtraEint_ = 1, Ee_ = 1

  

  ! The named index range for frequency bins in multi-group
  integer, parameter :: nWave = 1
  integer, parameter :: WaveFirst_ = 1, WaveLast_ = 1
  
end module ModExtraVariables
