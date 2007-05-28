
module ModSatellites

  use ModInputs, only : iCharLen_
  use ModIoUnit, only : UnitTmp_

  implicit none

  integer, parameter :: nMaxSats = 20
  integer, parameter :: nMaxSatInputLines = 10000
  integer, parameter :: nMaxSatPos = 50

  integer :: nSats = 0

  integer :: iSatUnit = UnitTmp_

  character (len=iCharLen_) :: cSatFileName(nMaxSats)

  Integer, parameter :: dblprec = selected_real_kind(14,200)

  real (kind=dblprec) :: SatTime(nMaxSats, nMaxSatInputLines)
  real                :: SatPos(nMaxSats, 3, nMaxSatPos, nMaxSatInputLines)
  real                :: SatCurrentPos(nMaxSats, 3, nMaxSatPos)
  real                :: SatDtPlot(nMaxSats)
  integer             :: nSatPos(nMaxSats, nMaxSatInputLines)
  integer             :: nSatLines(nMaxSats)
  integer             :: iSatCurrentIndex(nMaxSats)

  real                :: CurrentSatellitePosition(3)
  character (len=8)   :: CurrentSatelliteName

end module ModSatellites
