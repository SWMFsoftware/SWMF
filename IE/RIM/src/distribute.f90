
subroutine distribute

  use ModRIM
  use ModMpi
  use ModProcIE
  use ModParamRIM, only: iDebugLevel

  implicit none

  integer :: iError, iSize

  write(*,*) "mm(OuterMagJrAll) : ", minval(OuterMagJrAll), maxval(OuterMagJrAll)

  if (maxval(OuterMagJrAll) > -1.0e31) then
     call rearrange(OuterMagJrAll, OuterMagJr)
  else
     OuterMagJr = 0.0
  endif

  if (maxval(InnerMagJrAll) > -1.0e31) then
     call rearrange(InnerMagJrAll, InnerMagJr)
  else
     InnerMagJr = 0.0
  endif

  if (maxval(IonoJrAll) > -1.0e31) then
     call rearrange(IonoJrAll, IonoJr)
  else
     IonoJr = 0.0
  endif

  Jr = IonoJr + InnerMagJr + OuterMagJr

  if (maxval(OuterMagInvBAll) > -1.0e31) then
     call rearrange(OuterMagInvBAll, OuterMagInvB)
     call rearrange(OuterMagRhoAll, OuterMagRho)
     call rearrange(OuterMagPAll, OuterMagP)
     OuterMagT = -1.0
     where (OuterMagP > 0.0) OuterMagT = OuterMagP/OuterMagRho
  endif

contains

  subroutine rearrange(ValueAll, ValueLocal)

    real, intent(out) :: ValueLocal(0:nLons+1, nLats)
    real, intent(in)  :: ValueAll(nLats+2, nLons*nProc+1)

    integer :: iLat, iLon, iLatFrom, iLonFrom, iError, iLatTo, iLonTo

    !\
    ! We need to reverse the latitude, shift the longitudes, and swap the
    ! index order....
    !/

    ! Skip the North and South pole, since RIM doesn't actually use this
    ! information

    do iLat = 1, nLats
       ! Need to reverse latitudes
       iLatFrom = iLat+1
       iLatTo   = nLats-(iLat-1)
       do iLon = 0, nLons+1
          ! Need to shift longitudes
          iLonTo   = iLon
          iLonFrom = mod(iLon + iProc*nLons + nLons*nProc/2, nLons*nProc)
          if (iLonFrom == 0) iLonFrom = nLons*nProc
          ValueLocal(iLonTo, iLatTo) = ValueAll(iLatFrom, iLonFrom)
       enddo
    enddo

  end subroutine rearrange

end subroutine distribute
