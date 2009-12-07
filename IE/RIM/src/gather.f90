
subroutine gather

  use ModRIM
  use ModMpi
  use ModProcIE
  use ModParamRIM, only: iDebugLevel

  implicit none

  integer :: iError, iSize

  ! Let's try message passing some of this stuff to Processor 0

  iSize = (nLats+2) * (nProc*nLons+1)
  call rearrange(Potential, PotentialAll)

  if (iProc == 0) then
     cpcpn = maxval(PotentialAll(1:nLats/2,:)) - &
             minval(PotentialAll(1:nLats/2,:))
     cpcps = maxval(PotentialAll(nLats/2+1:nLats+2,:)) - &
             minval(PotentialAll(nLats/2+1:nLats+2,:))
  endif

  if (iDebugLevel > 1) &
     write(*,*) "RIM==> CPCP of total solution : ", &
          (maxval(PotentialAll) - minval(PotentialAll))/1000.0, " kV"
  if (iDebugLevel > 0) &
     write(*,"(a,1p,2e12.3,a)") " RIM=> CPCP South/North : ", &
     cpcps/1000.0, cpcpn/1000.0, " kV"

  call rearrange(JouleHeating,JouleHeatingAll)
  call rearrange(AveE, AveEAll)
  call rearrange(Eflux,EfluxAll)

contains

  subroutine rearrange(ValueLocal, ValueAll)

    real, intent(in)  :: ValueLocal(0:nLons+1, nLats)
    real, intent(out) :: ValueAll(nLats+2, nLons*nProc+1)

    integer :: iLat, iLon, iLatFrom, iLonFrom, iError, iLatTo, iLonTo

    iError = 0

    !\
    ! We need to reverse the latitude, shift the longitudes, and swap the
    ! index order....
    !/

    ! Skip the North and South pole, since this needs to be done on a case by
    ! case basis

    ValueAll = -1.0e32

    do iLat = 1, nLats
       ! Need to reverse latitudes
       iLatTo   = iLat+1
       iLatFrom = nLats-(iLat-1)
       do iLon = 1, nLons
          ! Need to shift longitudes
          iLonFrom = iLon
          iLonTo   = mod(iLon + nLons*iProc + nLonsAll/2, nLonsAll)
          if (iLonTo == 0) iLonTo = nLonsAll
          ValueAll(iLatTo, iLonTo) = ValueLocal(iLonFrom, iLatFrom)
       enddo
    enddo

    localVar = ValueAll

    call MPI_REDUCE(localVar, ValueAll, iSize, MPI_REAL, MPI_MAX, &
         0, iComm, iError)

    if (iProc == 0) then
       ! Pole is the average of the surrounding points.
       ValueAll(1,:)       = &
            sum(ValueAll(      2,1:nLons*nProc))/(nLons*nProc)
       ValueAll(nLats+2,:) = &
            sum(ValueAll(nLats+1,1:nLons*nProc))/(nLons*nProc)
       ! Fill in overlap point in longitude
       ValueAll(:,nLons*nProc+1) = ValueAll(:,1)
    endif

  end subroutine rearrange

end subroutine gather
