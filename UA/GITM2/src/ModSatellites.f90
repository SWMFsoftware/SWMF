
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

  real (kind=dblprec), allocatable :: SatTime(:,:)
  real, allocatable                :: SatPos(:,:,:,:)
  integer, allocatable             :: nSatPos(:,:)
  real                :: SatCurrentPos(nMaxSats, 3, nMaxSatPos)
  real                :: SatDtPlot(nMaxSats)
  integer             :: nSatLines(nMaxSats)
  integer             :: iSatCurrentIndex(nMaxSats)

  real                :: CurrentSatellitePosition(3)
  character (len=8)   :: CurrentSatelliteName

contains

  subroutine init_mod_satellites

    if(allocated(SatTime)) return
    allocate( &
         SatTime(nMaxSats, nMaxSatInputLines), &
         SatPos(nMaxSats, 3, nMaxSatPos, nMaxSatInputLines), &
         nSatPos(nMaxSats, nMaxSatInputLines))

  end subroutine init_mod_satellites

end module ModSatellites
