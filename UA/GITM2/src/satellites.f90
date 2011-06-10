
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

subroutine read_satellites(iError)

  use ModSatellites
  use ModGITM, only: iUp_, iEast_, iNorth_
  use ModInputs, only: iDebugLevel
  use ModConstants

  implicit none

  integer, intent(out)      :: iError
  integer                   :: iSat, iMax, i
  logical                   :: IsStartFound, IsFine
  character (len=iCharLen_) :: cLine
  real (kind = dblprec)     :: OldTime, NewTime
  integer                   :: itime(7)
  real                      :: pos(3)
  !---------------------------------------------------------------------
  call init_mod_satellites

  iError = 0
  nSatLines = 0
  nSatPos = 0
  iSatCurrentIndex = 0

  do iSat = 1, nSats

     if (iDebugLevel > 2) write(*,*) "Reading Satellite File : ",cSatFileName(iSat), iSat

     open(unit=iSatUnit,file=cSatFileName(iSat),status="old",iostat = iError)

     if (iError /= 0) then
        write(*,*) "Error opening satellite file : ",cSatFileName(iSat)
        return
     endif

     IsStartFound = .false.

     do while (.not. IsStartFound)
        read(iSatUnit, *, iostat = iError) cLine
        if (iError /= 0) IsStartFound = .true.
        if (index(cline,"#START") > 0) IsStartFound = .true.
     enddo

     if (iError /= 0) then
        write(*,*) "Error finding #START in satellite file : ",&
             cSatFileName(iSat)
        close(iSatUnit)
        return
     endif

     OldTime = 0.0
     IsFine = .false.

     do while (iError == 0)

        read(iSatUnit,*,iostat=iError) iTime, Pos

        if (iError == 0) then

           IsFine = .true.

           if (Pos(iEast_) < 0.0) Pos(iEast_) = Pos(iEast_) + 360.0
           Pos(iEast_)  = Pos(iEast_)*pi/180.0
           Pos(iNorth_) = Pos(iNorth_)*pi/180.0
           Pos(iUp_)    = Pos(iUp_)*1000.0

           call time_int_to_real(iTime, NewTime)

           if (NewTime /= OldTime) then
              nSatLines(iSat) = nSatLines(iSat) + 1
              if (nSatLines(iSat) > nMaxSatInputLines) then
                 write(*,*) "Too Many Lines in satfile : ",cSatFileName(iSat)
                 iError = 1
              else
                 nSatPos(iSat,nSatLines(iSat)) = &
                      nSatPos(iSat,nSatLines(iSat)) + 1
                 SatPos(iSat,1,nSatPos(iSat,nSatLines(iSat)),nSatLines(iSat))=Pos(1)
                 SatPos(iSat,2,nSatPos(iSat,nSatLines(iSat)),nSatLines(iSat))=Pos(2)
                 SatPos(iSat,3,nSatPos(iSat,nSatLines(iSat)),nSatLines(iSat))=Pos(3)
                 SatTime(iSat,nSatLines(iSat)) = NewTime
              endif
           else
              nSatPos(iSat,nSatLines(iSat)) = nSatPos(iSat,nSatLines(iSat)) + 1
              if (nSatPos(iSat,nSatLines(iSat)) > nMaxSatPos) then
                 write(*,*) "Too Many Positions in satfile : ", &
                      cSatFileName(iSat)
                 write(*,*) "Line : ",nSatLines(iSat)
                 iError = 1
              else
                 SatPos(iSat,:,nSatPos(iSat,nSatLines(iSat)),nSatLines(iSat))=&
                      Pos
              endif
           endif

           OldTime = NewTime

        endif

     enddo

     close(iSatUnit)
     if (IsFine) iError = 0

  enddo

  if (iDebugLevel > 1) then

     write(*,*) "Number of Satellite Files Read : ",nSats
     do iSat = 1, nSats
        write(*,*) iSat,". Number of Different times : ",nSatLines(iSat)
        iMax = 1
        do i=1,nSatLines(iSat)
           if (nSatPos(iSat,i) > iMax) iMax = nSatPos(iSat,i)
        enddo
        write(*,*) iSat,". Max number of Positions : ",iMax
     enddo

  endif

end subroutine read_satellites

!----------------------------------------------------------------------
! Find locations for current time
!----------------------------------------------------------------------

subroutine move_satellites

  use ModSatellites
  use ModGITM, only: nBlocks, dt
  use ModTime, only: CurrentTime, tSimulation

  implicit none

  integer :: iSat, iPos, iLine = 1, i, iBlock
  real    :: r
  character (len=8)         :: cName
  character (len=3)         :: cPos
  character (len=iCharLen_) :: cTmp

  do iSat = 1, nSats

     if (floor((tSimulation-dt)/SatDtPlot(iSat)) /= &
          floor((tsimulation)/SatDtPlot(iSat))) then

        if (iSatCurrentIndex(iSat) == 0) then
           if ( CurrentTime >= SatTime(iSat,1) .and. &
                CurrentTime <= SatTime(iSat,nSatLines(iSat))) then
              iLine = 1
              if (CurrentTime == SatTime(iSat,nSatLines(iSat))) then
                 iLine = nSatLines(iSat)
              else
                 do while (SatTime(iSat,iLine) < CurrentTime)
                    iLine = iLine + 1
                 enddo
                 iLine = iLine - 1
              endif
              iSatCurrentIndex(iSat) = iLine
           endif
        else
           iLine = iSatCurrentIndex(iSat)
           if (SatTime(iSat,nSatLines(iSat)) >= CurrentTime) then
              if (CurrentTime == SatTime(iSat,nSatLines(iSat))) then
                 iLine = nSatLines(iSat)
              else
                 do while (SatTime(iSat,iLine) < CurrentTime)
                    iLine = iLine + 1
                 enddo
                 iLine = iLine - 1
              endif
           else
              iLine = 0
           endif
           iSatCurrentIndex(iSat) = iLine
        endif

        if (iSatCurrentIndex(iSat) > 0) then

           if (iSatCurrentIndex(iSat) < nSatLines(iSat)) then
              iLine = iSatCurrentIndex(iSat)
              r =  1-(CurrentTime         - SatTime(iSat,iLine)) / &
                   (SatTime(iSat,iLine+1) - SatTime(iSat,iLine))
           else
              iLine = iSatCurrentIndex(iSat) - 1
              r = 1-(CurrentTime          - SatTime(iSat,iLine)) / &
                   (SatTime(iSat,iLine+1) - SatTime(iSat,iLine))
           endif

           do iPos = 1, nSatPos(iSat,iLine)
              do i=1,3
                 if (i == 1 .and. &
                      SatPos(iSat, i, iPos, iLine) > 300 .and. &
                      SatPos(iSat, i, iPos, iLine+1) < 60) then
                    SatCurrentPos(iSat, i, iPos) = &
                         (  r)*SatPos(iSat, i, iPos, iLine) + &
                         (1-r)*(SatPos(iSat, i, iPos, iLine+1)+360.0)
                 else
                    SatCurrentPos(iSat, i, iPos) = &
                         (  r)*SatPos(iSat, i, iPos, iLine) + &
                         (1-r)*SatPos(iSat, i, iPos, iLine+1)
                 endif

!                write(*,*) "Satellite : ",iSat, iSatCurrentIndex(iSat), &
!                     CurrentTime - SatTime(iSat,iSatCurrentIndex(iSat)), &
!                     iPos, i, SatCurrentPos(iSat, i, iPos)
              enddo
              write(cPos,'(i3.3)') iPos
              cTmp  = cSatFileName(iSat)
              cName = cTmp(1:4)//"_"//cPos
              CurrentSatellitePosition = SatCurrentPos(iSat,:,iPos)
              CurrentSatelliteName     = cName
              do iBlock = 1, nBlocks
!                 call output_1d("UA/data/",cName,iBlock, &
!                      SatCurrentPos(iSat,:,iPos))
                 call output("UA/data/",iBlock, -1)
              enddo
           enddo

        endif

     endif

  enddo

end subroutine move_satellites

