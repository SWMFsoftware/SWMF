module ModExtraVariables
  
  ! Define indexes that are not used in most equation modules

  ! Set impossible value for indexes
  integer, parameter :: Pe_ = 1, Ppar_ = 1, Hyp_ = 1, Erad_ = 1
  integer, parameter :: ExtraEint_ = 1, Ee_ = 1

  

  ! The named index range for frequency bins in multi-group
  integer, parameter :: nWave = 0
  integer, parameter :: WaveFirst_ = 1, WaveLast_ = 0
  

  integer,private :: i

  character(len=*),dimension(0:99), parameter :: NameNumber_I=(/ &
       ( ACHAR(IACHAR('0')+i/10)//ACHAR(IACHAR('0')+mod(i,10)),i=0,99)&
       /)

end module ModExtraVariables
