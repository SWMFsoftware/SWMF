!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModSatellites

  use ModInputs, only : iCharLen_
  use ModIoUnit, only : UnitTmp_

  implicit none

  integer, parameter :: nMaxSats = 20
  integer, parameter :: nMaxSatInputLines = 100000 !Asad: increased for RCAC
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

  integer :: CurrSat, nRCMRSat
  integer, dimension(nMaxSats) :: RCMRSat

  real, allocatable :: SatDat(:,:)
  real, allocatable :: SatCurrentDat(:)
  real, allocatable :: SatAltDat(:)

contains

  subroutine init_mod_satellites
    ! Asad: Added allocation for new variables SatDat, SatCurrentDat,
    !       and SatAltDat

    if(allocated(SatTime)) return
    allocate( &
         SatTime(nMaxSats, nMaxSatInputLines), &
         SatPos(nMaxSats, 3, nMaxSatPos, nMaxSatInputLines), &
         nSatPos(nMaxSats, nMaxSatInputLines), &
         SatDat(nMaxSats, nMaxSatInputLines), &
         SatCurrentDat(nMaxSats), SatAltDat(nMaxSats))

  end subroutine init_mod_satellites

end module ModSatellites
