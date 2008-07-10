
program PostProcess

  implicit none

  character (len=80) :: FileName, cLine
  character (len=6)  :: cBlock
  integer :: iStart, iError

  integer, parameter :: iInputUnit_  = 13
  integer, parameter :: iOutputUnit_ = 14
  integer, parameter :: nMaxVars = 200

  logical :: IsThere, DoWrite, UseGhostCells

  integer :: iYear,iMonth,iDay,iHour,iMinute,iSecond,iMilli
  integer :: nBlocksLon, nBlocksLat, nBlocksAlt
  integer :: nVars, nLons, nLats, nAlts, nLonsTotal, nLatsTotal, nAltsTotal
  integer :: nGhostLats, nGhostLons, nGhostAlts
  integer :: iLon, iLat, iAlt, iBlockLon, iBlockLat, iBlockAlt
  integer :: iiLon, iiLat, iiAlt
  real    :: Version
  real    :: Data(nMaxVars)

  real, allocatable :: AllData(:,:,:,:)

  character (len=40) :: Variables(nMaxVars)

  integer :: i, iVar, iBlock

  nGhostLats = 2
  nGhostLons = 2
  nGhostAlts = 2

  write(*,*) "Enter file group to process (can include .header) : "
  read(5,*) FileName

  iStart = index(FileName,".header")-1
  if (iStart > 0) then 
  else
     iStart = index(FileName," ")-1
  endif

  inquire(file=FileName(1:iStart)//".header",EXIST=IsThere)
  if (.not.IsThere) then
     write(*,*) "Could not find header file : ",FileName(1:iStart)//".header"
     stop
  endif

  open(iInputUnit_, file=FileName(1:iStart)//".header",status="old")

  iError = 0

  nBlocksAlt = 1
  nBlocksLat = 1
  nBlocksLon = 1

  UseGhostCells = .true.

  do while (iError == 0)

     read(iInputUnit_,'(a)',iostat=iError) cLine

     if (index(cLine,"BLOCKS") > 0) then
        read(iInputUnit_,*) nBlocksAlt
        read(iInputUnit_,*) nBlocksLat
        read(iInputUnit_,*) nBlocksLon
     endif

     if (index(cLine,"NUMERICAL") > 0) then
        read(iInputUnit_,*) nVars
        read(iInputUnit_,*) nAlts
        read(iInputUnit_,*) nLats
        read(iInputUnit_,*) nLons
     endif

     if (index(cLine,"TIME") > 0) then
        read(iInputUnit_,*) iYear
        read(iInputUnit_,*) iMonth
        read(iInputUnit_,*) iDay
        read(iInputUnit_,*) iHour
        read(iInputUnit_,*) iMinute
        read(iInputUnit_,*) iSecond
        read(iInputUnit_,*) iMilli
     endif

     if (index(cLine,"VERSION") > 0) then
        read(iInputUnit_,*) Version
     endif

     if (index(cLine,"VARIABLE") > 0) then
        do iVar = 1, nVars
           read(iInputUnit_,"(i7,a40)") i, Variables(iVar)
        enddo
     endif

     if (index(cLine,"NGHOSTCELLS") > 0) then
        read(iInputUnit_,*) nGhostAlts
        read(iInputUnit_,*) nGhostLats
        read(iInputUnit_,*) nGhostLons
     endif

     if (index(cLine,"NO GHOSTCELLS") > 0) then
        UseGhostCells = .false.
        nGhostLats = 0
        nGhostLons = 0
        nGhostAlts = 0
     endif

  enddo

  close(iInputUnit_)

  ! If it is 1D, then we know that we have to have 1,1,1,1.....
  if (nLons == 1 .and. nLats == 1) then
     nBlocksLon = 1
     nBlocksLat = 1
  endif

  write(*,*) "Inputs :"
  write(*,*) "  nBlocksLon, nBlocksLat, nBlocksAlt : ", &
       nBlocksLon, nBlocksLat, nBlocksAlt
  write(*,*) "  nLons,      nLats,      nAlts      : ", nLons, nLats, nAlts
  if (UseGhostCells) then
     nLonsTotal = nBlocksLon*(nLons-nGhostLons*2)+nGhostLons*2
!     write(*,*) nBlocksLon, nLons, nGhostLons
     nLatsTotal = nBlocksLat*(nLats-nGhostLats*2)+nGhostLats*2
!     write(*,*) nBlocksLat, nLats, nGhostLats, nLatsTotal
     nAltsTotal = nBlocksAlt*(nAlts-nGhostAlts*2)+nGhostAlts*2
  else
     nLonsTotal = nBlocksLon*nLons
     nLatsTotal = nBlocksLat*nLats
     nAltsTotal = nBlocksAlt*nAlts
  endif
  write(*,*) "  nLonsTotal, nLatsTotal, nAltsTotal : ", &
       nLonsTotal, nLatsTotal, nAltsTotal, " (Predicted Values)"


  allocate(AllData(nLonsTotal,nLatsTotal,nAltsTotal,nVars))

  open(iOutputUnit_,file=FileName(1:iStart)//".bin",&
       status="unknown",form="unformatted")

  write(iOutputUnit_) Version
  write(iOutputUnit_) nLonsTotal, nLatsTotal, nAltsTotal
  write(iOutputUnit_) nVars
  do iVar=1,nVars
     write(iOutputUnit_) Variables(iVar)
  enddo

  write(iOutputUnit_) iYear, iMonth, iDay, iHour, iMinute, iSecond, iMilli

  nAltsTotal = 0
  nLatsTotal = 0
  nLonsTotal = 0

  iBlock = 1
  do iBlockAlt = 1, nBlocksAlt
     do iBlockLat = 1, nBlocksLat
        do iBlockLon = 1, nBlocksLon

           write(cBlock,'(a2,i4.4)') ".b",iBlock

           inquire(file=FileName(1:iStart)//cBlock,EXIST=IsThere)
           if (.not.IsThere) then
              write(*,*) "Must be a satellite file...."
              cBlock = ".sat"
           endif

           write(*,*) "Reading File : ",FileName(1:iStart)//cBlock
           open(iInputUnit_, file=FileName(1:iStart)//cBlock, &
                status="old", form="unformatted")

           do iAlt = 1, nAlts
              do iLat = 1, nLats
                 do iLon = 1, nLons

                    read(iInputUnit_) Data(1:nVars)

                    DoWrite = .true.

                    !! Only write ghostcells on the edges

                    if (UseGhostCells) then

                       if (iBlockLon /= 1 .and. iLon <= nGhostLons) &
                            DoWrite = .false.
                       if (iBlockLon /= nBlocksLon .and. &
                            iLon > nLons-nGhostLons) DoWrite = .false.

                       if (iBlockLat /= 1 .and. iLat <= nGhostLats) &
                            DoWrite = .false.
                       if (iBlockLat /= nBlocksLat .and. &
                            iLat > nLats-nGhostLats) DoWrite = .false.

                       if (iBlockAlt /= 1 .and. iAlt <= nGhostAlts) &
                            DoWrite = .false.
                       if (iBlockAlt /= nBlocksAlt .and. &
                            iAlt > nAlts-nGhostAlts) DoWrite = .false.

                    endif

                    if (DoWrite) then

                       if (UseGhostCells) then

                          iiLon = (iBlockLon-1)*(nLons-nGhostLons*2) + iLon
                          iiLat = (iBlockLat-1)*(nLats-nGhostLats*2) + iLat
                          iiAlt = (iBlockAlt-1)*(nAlts-nGhostAlts*2) + iAlt

                       else

                          iiLon = (iBlockLon-1)*(nLons) + iLon
                          iiLat = (iBlockLat-1)*(nLats) + iLat
                          iiAlt = (iBlockAlt-1)*(nAlts) + iAlt

                       endif

                       AllData(iiLon,iiLat,iiAlt,1:nVars) = Data(1:nVars)

                       if (iiLon > nLonsTotal) nLonsTotal = iiLon
                       if (iiLat > nLatsTotal) nLatsTotal = iiLat
                       if (iiAlt > nAltsTotal) nAltsTotal = iiAlt


                    endif

                 enddo
              enddo
           enddo

           close(iInputUnit_)

           iBlock = iBlock + 1

        enddo
     enddo
  enddo

  do iVar = 1, nVars
     write(iOutputUnit_) AllData(:,:,:,iVar)
  enddo

  write(*,*) "Output :"
  write(*,*) "  nLons, nLats, nAlts : ", nLonsTotal, nLatsTotal, nAltsTotal

  close(iOutputUnit_)

end program PostProcess
