
subroutine distribute

  use ModRIM
  use ModMpi
  use ModProcIE
  use ModParamRIM, only: iDebugLevel
  use ModNumConst

  implicit none

  real :: lw
  integer :: iError, iSize, iLM, iLat, iLatGood, iLon, nSmooth
  integer :: iLonFrom, iSubLon, iSubLat
  real, dimension(:), allocatable    :: SmoothVar
  integer, dimension(:), allocatable :: iSmooth

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

  if (.not.allocated(SmoothVar)) then
     allocate(SmoothVar(nLonsAll+1))
     SmoothVar = 0.0
     allocate(iSmooth(nLonsAll+1))
     iSmooth = 0
  endif

  InnerMagEFlux = -1.0
  InnerMagAveE  = -1.0
  if (maxval(InnerMagJrAll) > -1.0e31) then
     ! If we are getting solution from the RCM, only the North values are
     ! supplied.  Therefore, before we rearrage everything, we have to 
     ! map the Northern hemisphere solution to the Southern hemisphere
     if (maxval(InnerMagJrAll(nLats/2+1:nLats+1,:)) <= 0.0) then
        do iLat = 1, nLats/2+1
           iLM = nLats+2 - iLat + 1
           InnerMagJrAll(iLM,:)    = InnerMagJrAll(iLat,:)
           InnerMagEFluxAll(iLM,:) = InnerMagEFluxAll(iLat,:)
           InnerMagAveEAll(iLM,:)  = InnerMagAveEAll(iLat,:)
        enddo
     endif

     nSmooth = 2
     do iLat = nSmooth+1, nLats-nSmooth
        do iLon = 1, nLonsAll+1
           SmoothVar(iLon) = 0.0
           iSmooth(iLon) = 0
           do iSubLon = iLon-nSmooth, iLon+nSmooth
              iLonFrom = iSubLon
              if (iLonFrom <= 0) iLonFrom = iLonFrom + nLonsAll
              if (iLonFrom > nLonsAll) iLonFrom = iLonFrom - nLonsAll
              if (OuterMagRhoAll(iLat,iLonFrom) > 0.0) then
                 SmoothVar(iLon) = SmoothVar(iLon) + &
                      InnerMagJrAll(iLat,iLonFrom)
                 iSmooth(iLon) = iSmooth(iLon) + 1
              endif
           enddo
           do iSubLat = iLat-nSmooth, iLat+nSmooth
              if (OuterMagRhoAll(iSubLat,iLon) > 0.0) then
                 SmoothVar(iLon) = SmoothVar(iLon) + &
                      InnerMagJrAll(iSubLat,iLon)
                 iSmooth(iLon) = iSmooth(iLon) + 1
              endif
           enddo
        enddo
        do iLon = 1, nLonsAll+1
           if (iSmooth(iLon) > 0) then
              InnerMagJrAll(iLat,iLon) = SmoothVar(iLon)/(2*iSmooth(iLon)+1)
           endif
        enddo
     enddo

     nSmooth = 2
     do iLat = nSmooth+1, nLats-nSmooth
        do iLon = 1, nLonsAll+1
           SmoothVar(iLon) = 0.0
           iSmooth(iLon) = 0
           do iSubLon = iLon-nSmooth, iLon+nSmooth
              iLonFrom = iSubLon
              if (iLonFrom <= 0) iLonFrom = iLonFrom + nLonsAll
              if (iLonFrom > nLonsAll) iLonFrom = iLonFrom - nLonsAll
              if (OuterMagRhoAll(iLat,iLonFrom) > 0.0) then
                 SmoothVar(iLon) = SmoothVar(iLon) + &
                      InnerMagEFluxAll(iLat,iLonFrom)
                 iSmooth(iLon) = iSmooth(iLon) + 1
              endif
           enddo
           do iSubLat = iLat-nSmooth, iLat+nSmooth
              if (OuterMagRhoAll(iSubLat,iLon) > 0.0) then
                 SmoothVar(iLon) = SmoothVar(iLon) + &
                      InnerMagEFluxAll(iSubLat,iLon)
                 iSmooth(iLon) = iSmooth(iLon) + 1
              endif
           enddo
        enddo
        do iLon = 1, nLonsAll+1
           if (iSmooth(iLon) > 0) then
              InnerMagEFluxAll(iLat,iLon) = SmoothVar(iLon)/(2*iSmooth(iLon)+1)
           endif
        enddo
     enddo

     nSmooth = 2
     do iLat = nSmooth+1, nLats-nSmooth
        do iLon = 1, nLonsAll+1
           SmoothVar(iLon) = 0.0
           iSmooth(iLon) = 0
           do iSubLon = iLon-nSmooth, iLon+nSmooth
              iLonFrom = iSubLon
              if (iLonFrom <= 0) iLonFrom = iLonFrom + nLonsAll
              if (iLonFrom > nLonsAll) iLonFrom = iLonFrom - nLonsAll
              if (OuterMagRhoAll(iLat,iLonFrom) > 0.0) then
                 SmoothVar(iLon) = SmoothVar(iLon) + &
                      InnerMagAveEAll(iLat,iLonFrom)
                 iSmooth(iLon) = iSmooth(iLon) + 1
              endif
           enddo
           do iSubLat = iLat-nSmooth, iLat+nSmooth
              if (OuterMagRhoAll(iSubLat,iLon) > 0.0) then
                 SmoothVar(iLon) = SmoothVar(iLon) + &
                      InnerMagAveEAll(iSubLat,iLon)
                 iSmooth(iLon) = iSmooth(iLon) + 1
              endif
           enddo
        enddo
        do iLon = 1, nLonsAll+1
           if (iSmooth(iLon) > 0) then
              InnerMagAveEAll(iLat,iLon) = SmoothVar(iLon)/(2*iSmooth(iLon)+1)
           endif
        enddo
     enddo

     call rearrange(InnerMagJrAll, InnerMagJr)
     if (maxval(InnerMagEFluxAll) > 0.0) then
        call rearrange(InnerMagEFluxAll, InnerMagEFlux)
        call rearrange(InnerMagAveEAll,  InnerMagAveE)
        where(InnerMagAveE < 0.1) InnerMagAveE = 0.1
        where(InnerMagEFlux < 0.1) InnerMagEFlux = 0.1
     endif
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

     lw = 10.0*3.1415/180.0
     do iLon = 1, nLonsAll+1
        iLatGood = 1
        do while (OuterMagRhoAll(iLatGood,iLon) == 0.0 .and. iLatGood<nLats/2)
           iLatGood = iLatGood + 1
        enddo
        if (iLatGood > 1) then
           OuterMagRhoAll(1:iLatGood-1,iLon) = OuterMagRhoAll(iLatGood,iLon)* &
                exp(-abs(LatitudeAll(1:iLatGood-1,iLon)-LatitudeAll(iLatGood,iLon))/lw)
           OuterMagPAll(1:iLatGood-1,iLon) = OuterMagPAll(iLatGood,iLon) * &
                exp(-abs(LatitudeAll(1:iLatGood-1,iLon)-LatitudeAll(iLatGood,iLon))/lw)
        endif
     enddo
 
     do iLon = 1, nLonsAll+1
        iLatGood = nLats+2
        do while (OuterMagRhoAll(iLatGood,iLon) == 0.0 .and. iLatGood > nLats/2)
           iLatGood = iLatGood - 1
        enddo
        OuterMagRhoAll(iLatGood+1:nLats+2,iLon) = OuterMagRhoAll(iLatGood,iLon) * &
                exp(-abs(LatitudeAll(iLatGood+1:nLats+2,iLon)-LatitudeAll(iLatGood,iLon))/lw)
        OuterMagPAll(iLatGood+1:nLats+2,iLon) = OuterMagPAll(iLatGood,iLon) * &
                exp(-abs(LatitudeAll(iLatGood+1:nLats+2,iLon)-LatitudeAll(iLatGood,iLon))/lw)
     enddo

     nSmooth = 4
     do iLat = 1, nLats
        do iLon = 1, nLonsAll+1
           SmoothVar(iLon) = 0.0
           do iSubLon = iLon-nSmooth, iLon+nSmooth
              iLonFrom = iSubLon
              if (iLonFrom <= 0) iLonFrom = iLonFrom + nLonsAll
              if (iLonFrom > nLonsAll) iLonFrom = iLonFrom - nLonsAll
              SmoothVar(iLon) = SmoothVar(iLon) + &
                   OuterMagRhoAll(iLat,iLonFrom)
           enddo
        enddo
        OuterMagRhoAll(iLat,:) = SmoothVar/(2*nSmooth+1)
     enddo

     do iLat = 1, nLats
        do iLon = 1, nLonsAll+1
           SmoothVar(iLon) = 0.0
           do iSubLon = iLon-nSmooth, iLon+nSmooth
              iLonFrom = iSubLon
              if (iLonFrom <= 0) iLonFrom = iLonFrom + nLonsAll
              if (iLonFrom > nLonsAll) iLonFrom = iLonFrom - nLonsAll
              SmoothVar(iLon) = SmoothVar(iLon) + &
                   OuterMagPAll(iLat,iLonFrom)
           enddo
        enddo
        OuterMagPAll(iLat,:) = SmoothVar/(2*nSmooth+1)
     enddo

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
