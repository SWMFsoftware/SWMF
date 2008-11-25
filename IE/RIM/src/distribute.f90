
subroutine distribute

  use ModRIM
  use ModMpi
  use ModProcIE
  use ModParamRIM, only: iDebugLevel
  use ModNumConst

  implicit none

  integer :: iError, iSize

  if (iDebugLevel > 1) then 
     write(*,*) "RIM==> mm(OuterMagJrAll):", &
          minval(OuterMagJrAll), maxval(OuterMagJrAll)
     write(*,*) "RIM==> mm(InnerMagJrAll):", &
          minval(InnerMagJrAll), maxval(InnerMagJrAll)
     write(*,*) "RIM==> mm(IonoJrAll):", &
          minval(IonoJrAll), maxval(IonoJrAll)
  endif

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
     ! This creates full domain arrays of these, but rearranged

     OuterMagInvBAllR = -1.0
     OuterMagRhoAllR = -1.0
     OuterMagPAllR = -1.0
     OuterMagTAllR = -1.0
     LatitudeAllR = -1.0
     
     call rearrangeR(OuterMagInvBAll, OuterMagInvBAllR)
     call rearrangeR(OuterMagRhoAll, OuterMagRhoAllR)
     call rearrangeR(OuterMagPAll, OuterMagPAllR)
     where (OuterMagRhoAllR> 0.0) OuterMagTAllR = OuterMagPAllR/OuterMagRhoAllR
     call rearrangeR(LatitudeAll, LatitudeAllR)
     LatitudeAllR = cPi/2 - LatitudeAllR
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

  subroutine rearrangeR(ValueAll, ValueLocal)

    real, intent(out) :: ValueLocal(0:nLons*nProc+1, nLats)
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
       do iLon = 0, nLons*nProc+1
          ! Need to shift longitudes
          iLonTo   = iLon
          iLonFrom = mod(iLon + nLons*nProc/2, nLons*nProc)
          if (iLonFrom == 0) iLonFrom = nLons*nProc
          ValueLocal(iLonTo, iLatTo) = ValueAll(iLatFrom, iLonFrom)
       enddo
    enddo

  end subroutine rearrangeR
  
end subroutine distribute
