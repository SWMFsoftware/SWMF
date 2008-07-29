
module ModHeidiSatellites

  integer, parameter :: nMaxSatellites = 50
  integer            :: nImSats = 0
  logical            :: DoWriteSats = .false.
  real               :: SatLoc_3I(3,2,nMaxSatellites)
  character(len=100) :: NameSat_I(nMaxSatellites)

end module ModHeidiSatellites
