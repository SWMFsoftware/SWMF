
!----------------------------------------------------------------
!
!----------------------------------------------------------------

integer function bad_outputtype()

  use ModInputs, only : OutputType, nOutputTypes

  implicit none

  integer :: iOutputType
  logical :: IsFound

  do iOutputType = 1, nOutputTypes

     IsFound = .false.

     if (OutputType(iOutputType) == '3DAUR')     IsFound = .true.
     if (OutputType(iOutputType) == '3DSRC')     IsFound = .true.
     if (OutputType(iOutputType) == '3DALL')     IsFound = .true.
     if (OutputType(iOutputType) == '1DALL')     IsFound = .true.
     if (OutputType(iOutputType) == '3DNEUTRAL') IsFound = .true.
     if (OutputType(iOutputType) == '3DION')     IsFound = .true.
     if (OutputType(iOutputType) == '2DGEO')     IsFound = .true.
     if (OutputType(iOutputType) == '2DMAG')     IsFound = .true.

     if (.not. IsFound) then
        bad_outputtype = iOutputType
        return
     endif

  enddo

  bad_outputtype = 0
  return

end function bad_outputtype

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_1d(dir, cName, iBlock, Position)

  use ModGITM
  use ModEUV
  use ModTime
  use ModInputs
  use ModSources
  use ModConstants

  implicit none

  character (len=*), intent(in) :: dir
  character (len=*), intent(in) :: cName
  integer, intent(in)           :: iBlock
  real, intent(in)              :: Position(3)

  character (len=14) :: cTime
  integer :: iLon,iLat,iAlt, nvars_to_write, nlines, iBLK, i
  integer :: iiLon, iiLat, iiAlt, nGCs
  logical :: done

  real :: LatFind, LonFind
  real :: rLon, rLat

  character (len=2) :: cYear, cMonth, cDay, cHour, cMinute, cSecond

  LatFind = Position(iNorth_)
  LonFind = Position(iEast_)

  if ((Longitude(0,iBlock)+Longitude(1,iBlock))/2 <=LonFind .and. &
       (Longitude(nLons,iBlock)+Longitude(nLons+1,iBlock))/2 >LonFind) then
     if ((Latitude(0,iBlock)+Latitude(1,iBlock))/2 <=LatFind .and. &
          (Latitude(nLats,iBlock)+Latitude(nLats+1,iBlock))/2 >LatFind) then

        iiLat = -1
        iiLon = -1
        do iLon = 0,nLons
           if (Longitude(iLon,iBlock) <= LonFind .and. &
                Longitude(iLon+1,iBlock) > LonFind) then
              iiLon = iLon
              rLon = 1.0 - (LonFind - Longitude(iLon,iBlock)) / &
                   (Longitude(iLon+1,iBlock) - Longitude(iLon,iBlock))
           endif
        enddo

        do iLat = 0,nLats
           if (Latitude(iLat,iBlock) <= LatFind .and. &
                Latitude(iLat+1,iBlock) > LatFind) then
              iiLat = iLat
              rLat = 1.0 - (LatFind - Latitude(iLat,iBlock)) / &
                   (Latitude(iLat+1,iBlock) - Latitude(iLat,iBlock))
           endif
        enddo

     else
        return
     endif
  else 
     return
  endif

  if (iProc == 0 .and. iBlock == 1) &
       write(*,'(a,i7,i5,5i3)') &
       "Writing Output files at iStep : ",iStep, iTimeArray(1:6)

  call calc_physics(iBlock)
  call chapman_integrals(iBlock)
  call calc_rates(iBlock)
  if (.not. Is1D) call calc_efield(iBlock)

  !! construct naming strings

  call i2s(mod(iTimeArray(1),100), cYear, 2)
  call i2s(iTimeArray(2), cMonth, 2)
  call i2s(iTimeArray(3), cDay, 2)
  call i2s(iTimeArray(4), cHour, 2)
  call i2s(iTimeArray(5), cMinute, 2)
  call i2s(iTimeArray(6), cSecond, 2)

  cTime = "t"//cYear//cMonth//cDay//"_"//cHour//cMinute//cSecond

  !! open file
  open(unit=iOutputUnit_, &
       file=dir//"/"//cName//"_"//cTime//"."//"dat",&
       status="unknown")

  write(iOutputUnit_,*) ""
  write(iOutputUnit_,*) "TIME"
  do i=1,7
     write(iOutputUnit_,*) iTimeArray(i)
  enddo
  write(iOutputUnit_,*) ""

  call output_1dall(iiLon, iiLat, iBlock, rLon, rLat, iOutputUnit_)

  close(unit=iOutputUnit_)

end subroutine output_1d

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output(dir, iBlock, iOutputType)

  use ModSphereInterface, only:iStartBlk
  use ModGITM
  use ModEUV
  use ModTime
  use ModInputs
  use ModSources

  implicit none

  character (len=*), intent(in) :: dir
  integer, intent(in) :: iBlock
  integer, intent(in) :: iOutputType

  character (len=5) :: proc_str,blk_str
  character (len=14) :: cTime
  character (len=100) :: output_format
  integer :: iLon,iLat,iAlt, nvars_to_write, nlines, iBLK
  integer :: iiLon, iiLat, iiAlt, nGCs
  logical :: done

  real :: LatFind = 2.0*pi/180.0, LonFind = 10.0*pi/180.0
  real :: rLon, rLat

  character (len=2) :: cYear, cMonth, cDay, cHour, cMinute, cSecond

  !! construct naming strings

  if (index(OutputType(iOutputType), "1D") > 0) then
     if ((Longitude(0,iBlock)+Longitude(1,iBlock))/2 <=LonFind .and. &
          (Longitude(nLons,iBlock)+Longitude(nLons+1,iBlock))/2 >LonFind) then
        if ((Latitude(0,iBlock)+Latitude(1,iBlock))/2 <=LatFind .and. &
             (Latitude(nLats,iBlock)+Latitude(nLats+1,iBlock))/2 >LatFind) then

           iiLat = -1
           iiLon = -1
           do iLon = 0,nLons
              if (Longitude(iLon,iBlock) <= LonFind .and. &
                   Longitude(iLon+1,iBlock) > LonFind) then
                 iiLon = iLon
                 rLon = 1.0 - (LonFind - Longitude(iLon,iBlock)) / &
                      (Longitude(iLon+1,iBlock) - Longitude(iLon,iBlock))
              endif
           enddo

           do iLat = 0,nLats
              if (Latitude(iLat,iBlock) <= LatFind .and. &
                   Latitude(iLat+1,iBlock) > LatFind) then
                 iiLat = iLat
                 rLat = 1.0 - (LatFind - Latitude(iLat,iBlock)) / &
                      (Latitude(iLat+1,iBlock) - Latitude(iLat,iBlock))
              endif
           enddo

        else
           return
        endif
     else 
        return
     endif
  endif

  if (iProc == 0 .and. iBlock == 1) &
       write(*,'(a,i7,i5,5i3)') &
       "Writing Output files at iStep : ",iStep, iTimeArray(1:6)

  call calc_physics(iBlock)
  call chapman_integrals(iBlock)
  call calc_rates(iBlock)
  if (.not. Is1D) call calc_efield(iBlock)

  iBLK = iStartBLK + iBlock

  write(blk_str,'(a1,i4.4)') "b",iBLK

  call i2s(mod(iTimeArray(1),100), cYear, 2)
  call i2s(iTimeArray(2), cMonth, 2)
  call i2s(iTimeArray(3), cDay, 2)
  call i2s(iTimeArray(4), cHour, 2)
  call i2s(iTimeArray(5), cMinute, 2)
  call i2s(iTimeArray(6), cSecond, 2)

csecond="00"
  cTime = "t"//cYear//cMonth//cDay//"_"//cHour//cMinute//cSecond

  select case (OutputType(iOutputType))

  case ('3DALL')

     !! open file
     open(unit=iOutputUnit_, &
          file=dir//"/"//blk_str//"_"//cTime//"."//OutputType(iOutputType),&
          status="unknown")

     call write_head_blocks
     call write_head_time
     call write_head_version
     call output_3dall

     !! close file
     close(unit=iOutputUnit_)

  case ('3DSRC')

     !! open file
     open(unit=iOutputUnit_, &
          file=dir//"/"//blk_str//"_"//cTime//"."//OutputType(iOutputType),&
          status="unknown")

     call write_head_blocks
     call write_head_time
     call write_head_version
     call output_3dsrc

     !! close file
     close(unit=iOutputUnit_)

  case ('3DAUR')

     !! open file
     open(unit=iOutputUnit_, &
          file=dir//"/"//blk_str//"_"//cTime//"."//OutputType(iOutputType),&
          status="unknown")

     call write_head_blocks
     call write_head_time
     call write_head_version
     call output_3daur

     !! close file
     close(unit=iOutputUnit_)

  case ('2DGEO')

     !! open file
     open(unit=iOutputUnit_, &
          file=dir//"/"//blk_str//"_"//cTime//"."//OutputType(iOutputType),&
          status="unknown")

     call write_head_blocks
     call write_head_time
     call write_head_version
     call output_2d

     !! close file
     close(unit=iOutputUnit_)

  case ('2DMAG')

     if (iBLK == 1) then

        !! open file
        open(unit=iOutputUnit_, &
             file=dir//"/"//blk_str//"_"//cTime//"."//OutputType(iOutputType),&
             status="unknown")

        call write_head_time
        call write_head_version
        call output_2delectro

        !! close file
        close(unit=iOutputUnit_)

     endif

  case ('1DALL')

     !! open file
     open(unit=iOutputUnit_, &
          file=dir//"/"//blk_str//"_"//cTime//"."//OutputType(iOutputType),&
          status="unknown")

     call write_head_time
     call write_head_version
     call output_1dall(iiLon, iiLat, iBlock, rLon, rLat, iOutputUnit_)

     !! close file
     close(unit=iOutputUnit_)

  case ('3DION')

     !! open file
     open(unit=iOutputUnit_, &
          file=dir//"/"//blk_str//"_"//cTime//"."//OutputType(iOutputType),&
          status="unknown")

     call write_head_blocks
     call write_head_time
     call write_head_version
     call output_3dion

     !! close file
     close(unit=iOutputUnit_)

  end select

contains

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine write_head_blocks

    write(iOutputUnit_,*) "BLOCKS"
    write(iOutputUnit_,"(I7,A)") 1, " nBlocksAlt"
    write(iOutputUnit_,"(I7,A)") nBlocksLat, " nBlocksLat"
    write(iOutputUnit_,"(I7,A)") nBlocksLon, " nBlocksLon"
    write(iOutputUnit_,*) ""

  end subroutine write_head_blocks


  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine write_head_time

    write(iOutputUnit_,*) "TIME"
    write(iOutputUnit_,"(I7,A)") iTimeArray(1), " Year"
    write(iOutputUnit_,"(I7,A)") iTimeArray(2), " Month"
    write(iOutputUnit_,"(I7,A)") iTimeArray(3), " Day"
    write(iOutputUnit_,"(I7,A)") iTimeArray(4), " Hour"
    write(iOutputUnit_,"(I7,A)") iTimeArray(5), " Minute"
    write(iOutputUnit_,"(I7,A)") iTimeArray(6), " Second"
    write(iOutputUnit_,*) ""

  end subroutine write_head_time

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine write_head_version

    write(iOutputUnit_,*) "VERSION"
    write(iOutputUnit_,*) 2.4
    write(iOutputUnit_,*) ""

  end subroutine write_head_version

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine output_2d

    use ModElectrodynamics

    nvars_to_write = 12
    write(output_format,"('(1p,',I2,'E11.3)')") nvars_to_write

    write(iOutputUnit_,*) "NUMERICAL VALUES"
    write(iOutputUnit_,"(I7,A)") nvars_to_write, " nvars"
    write(iOutputUnit_,"(I7,A)") 1, " nAltitudes"
    write(iOutputUnit_,"(I7,A)") nLats, " nLatitude"
    write(iOutputUnit_,"(I7,A)") nLons+1, " nLongitudes"

    write(iOutputUnit_,*) ""

    write(iOutputUnit_,*) "VARIABLE LIST"
    write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
    write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
    write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
    write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Potential"
    write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Pedersen Conductance"
    write(iOutputUnit_,"(I7,A1,a)")  6, " ", "Hall Conductance"
    write(iOutputUnit_,"(I7,A1,a)")  7, " ", "Electron_Average_Energy"
    write(iOutputUnit_,"(I7,A1,a)")  8, " ", "Electron_Energy_Flux"
    write(iOutputUnit_,"(I7,A1,a)")  9, " ", "DivJuAlt"
    write(iOutputUnit_,"(I7,A1,a)") 10, " ", "Pedersen FL Conductance"
    write(iOutputUnit_,"(I7,A1,a)") 11, " ", "Hall FL Conductance"
    write(iOutputUnit_,"(I7,A1,a)") 12, " ", "FL Length"
    write(iOutputUnit_,*) ""
    write(iOutputUnit_,*) "BEGIN"

    iAlt = 1
    do iLat=1,nLats
       do iLon=1,nLons+1
          write(iOutputUnit_,output_format)       &
               Longitude(iLon,iBlock), &
               Latitude(iLat,iBlock),&
               Altitude(iAlt),&
               Potential(iLon,iLat,iAlt,iBlock), &
               PedersenConductance(iLon,iLat,iBlock), &
               HallConductance(iLon,iLat,iBlock), &
               ElectronAverageEnergy(iLon,iLat), &
               ElectronEnergyFlux(iLon,iLat), &
               DivJuAlt(iLon,iLat), &
               PedersenFieldLine(iLon, iLat), &
               HallFieldLine(iLon, iLat), &
               LengthFieldLine(iLon, iLat)
       enddo
    enddo

  end subroutine output_2d

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine output_2delectro

    use ModElectrodynamics
    use ModConstants, only:Pi

    nvars_to_write = 7
    write(output_format,"('(1p,',I2,'E11.3)')") nvars_to_write

    write(iOutputUnit_,*) "BLOCKS"
    write(iOutputUnit_,"(I7,A)") 1, " nBlocksAlt"
    write(iOutputUnit_,"(I7,A)") 1, " nBlocksLat"
    write(iOutputUnit_,"(I7,A)") 1, " nBlocksLon"
    write(iOutputUnit_,*) ""

    write(iOutputUnit_,*) "NUMERICAL VALUES"
    write(iOutputUnit_,"(I7,A)") nvars_to_write, " nvars"
    write(iOutputUnit_,"(I7,A)") 1, " nAltitudes"
    write(iOutputUnit_,"(I7,A)") nMagLats, " nLatitude"
    write(iOutputUnit_,"(I7,A)") nMagLons+1, " nLongitudes"

    write(iOutputUnit_,*) ""

    write(iOutputUnit_,*) "VARIABLE LIST"
    write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
    write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
    write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
    write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Pedersen Conductance"
    write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Hall Conductance"
    write(iOutputUnit_,"(I7,A1,a)")  6, " ", "DivJuAlt"
    write(iOutputUnit_,"(I7,A1,a)")  6, " ", "Field Line Length"
    write(iOutputUnit_,*) ""
    write(iOutputUnit_,*) "BEGIN"

    iAlt = 1
    do iLat=1,nMagLats
       do iLon=1,nMagLons+1
          write(iOutputUnit_,output_format)       &
               MagLocTimeMC(iLon,iLat)*Pi/12.0, &
               MagLatMC(iLon,iLat)*Pi/180.0,Altitude(iAlt),&
               SigmaPedersenMC(iLon,iLat), &
               SigmaHallMC(iLon,iLat), &
               DivJuAltMC(iLon,iLat), &
               LengthMC(iLon,iLat)

       enddo
    enddo

  end subroutine output_2delectro

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine output_3dall

    nvars_to_write = 33+nSpecies
    write(output_format,"('(1p,',I2,'E11.3)')") nvars_to_write

    if (Is1D) then
       nGCs = 0
    else
       nGCs = 2
    endif

    write(iOutputUnit_,*) "NUMERICAL VALUES"
    write(iOutputUnit_,"(I7,6A)") nvars_to_write, " nvars"
    write(iOutputUnit_,"(I7,7A)") nAlts+4, " nAltitudes"
    write(iOutputUnit_,"(I7,7A)") nLats+nGCs*2, " nLatitude"
    write(iOutputUnit_,"(I7,7A)") nLons+nGCs*2, " nLongitudes"
    write(iOutputUnit_,*) ""

    write(iOutputUnit_,*) "VARIABLE LIST"
    write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
    write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
    write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
    write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Rho"
    write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Temperature"
    write(iOutputUnit_,"(I7,A1,a)")  6, " ", "[O]"
    write(iOutputUnit_,"(I7,A1,a)")  7, " ", "[O2]"
    write(iOutputUnit_,"(I7,A1,a)")  8, " ", "[N2]"
    write(iOutputUnit_,"(I7,A1,a)")  9, " ", "[NO]"
    write(iOutputUnit_,"(I7,A1,a)") 10, " ", "[N_4S]"
    write(iOutputUnit_,"(I7,A1,a)") 11, " ", "[N_2D]"
    write(iOutputUnit_,"(I7,A1,a)") 12, " ", "Vn (east)"
    write(iOutputUnit_,"(I7,A1,a)") 13, " ", "Vn (north)"
    write(iOutputUnit_,"(I7,A1,a)") 14, " ", "Vn (up)"
    write(iOutputUnit_,"(I7,A1,a)") 15, " ", "[O+]"
    write(iOutputUnit_,"(I7,A1,a)") 16, " ", "[N+]"
    write(iOutputUnit_,"(I7,A1,a)") 17, " ", "[O2+]"
    write(iOutputUnit_,"(I7,A1,a)") 18, " ", "[N2+]"
    write(iOutputUnit_,"(I7,A1,a)") 19, " ", "[NO+]"
    write(iOutputUnit_,"(I7,A1,a)") 20, " ", "[e-]"
    write(iOutputUnit_,"(I7,A1,a)") 21, " ", "eTemperature"
    write(iOutputUnit_,"(I7,A1,a)") 22, " ", "iTemperature"
    write(iOutputUnit_,"(I7,A1,a)") 23, " ", "Vi (east)"
    write(iOutputUnit_,"(I7,A1,a)") 24, " ", "Vi (north)"
    write(iOutputUnit_,"(I7,A1,a)") 25, " ", "Vi (up)"
    write(iOutputUnit_,"(I7,A1,a)") 26, " ", "Potential (kV)"
    write(iOutputUnit_,"(I7,A1,a)") 27, " ", "ExB(east)"
    write(iOutputUnit_,"(I7,A1,a)") 28, " ", "ExB (north)"   
    write(iOutputUnit_,"(I7,A1,a)") 29, " ", "ExB (up)" 
    write(iOutputUnit_,"(I7,A1,a)") 30, " ", "EUV_IonRate (O+)" 
    write(iOutputUnit_,"(I7,A1,a)") 31, " ", "Aurora_IonRate (O+)" 
    write(iOutputUnit_,"(I7,A1,a)") 32, " ", "Joule heating (W/m3)" 
    write(iOutputUnit_,"(I7,A1,a)") 33, " ", "Heating rate (K/s)" 

    write(iOutputUnit_,"(I7,A1,a)") 34, " ", "Vn (up,O)"
    write(iOutputUnit_,"(I7,A1,a)") 35, " ", "Vn (up,O2)"
    if (nSpecies >= iN2_) &
         write(iOutputUnit_,"(I7,A1,a)") 36, " ", "Vn (up,N2)"
    if (nSpecies >= iN_4S_) &
         write(iOutputUnit_,"(I7,A1,a)") 37, " ", "Vn (up,N)"

    write(iOutputUnit_,*) ""

    write(iOutputUnit_,*) "BEGIN"

    do iAlt=-1,nAlts+2
       iiAlt = max(min(iAlt,nAlts),1)
       do iLat=1-nGCs,nLats+nGCs
          iiLat = min(max(iLat,1),nLats)
          do iLon=1-nGCs,nLons+nGCs
             iiLon = min(max(iLon,1),nLons)
             write(iOutputUnit_,output_format)       &
                  Longitude(iLon,iBlock), &
                  Latitude(iLat,iBlock), &
                  altitude(iAlt),&
                  Rho(iLon,iLat,iAlt,iBlock),&
                  Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt),&
                  NDensityS(iLon,iLat,iAlt,iO_,iBlock), &
                  NDensityS(iLon,iLat,iAlt,iO2_,iBlock), &
                  NDensityS(iLon,iLat,iAlt,iN2_,iBlock), &
                  NDensityS(iLon,iLat,iAlt,iNO_,iBlock), &
                  NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock), &
                  NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock), &
                  velocity(iLon,iLat,iAlt,iEast_,iBlock) , &
                  velocity(iLon,iLat,iAlt,iNorth_,iBlock), &
                  velocity(iLon,iLat,iAlt,iUp_,iBlock)   , &
                  IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock),&
                  IDensityS(iLon,iLat,iAlt,iNP_,iBlock) , &
                  IDensityS(iLon,iLat,iAlt,iO2P_,iBlock), &
                  IDensityS(iLon,iLat,iAlt,iN2P_,iBlock), &
                  IDensityS(iLon,iLat,iAlt,iNOP_,iBlock), &
                  IDensityS(iLon,iLat,iAlt,ie_,iBlock), &
                  eTemperature(iLon,iLat,iAlt,iBlock)  ,&
                  ITemperature(iLon,iLat,iAlt,iBlock)  ,&
                  Ivelocity(iLon,iLat,iAlt,iEast_:iUp_,iBlock), &
                  Potential(iiLon,iiLat,iiAlt,iBlock), &
                  ExB(iLon,iLat,iAlt,iEast_:iUp_),& 
                  EuvIonRateS(iLon,iLat,iAlt,iO_4SP_,iBlock)+&
                  EuvIonRateS(iLon,iLat,iAlt,iO_2PP_,iBlock)+&
                  EuvIonRateS(iLon,iLat,iAlt,iO_2DP_,iBlock),&
                  AuroralIonRateS(iLon,iLat,iAlt,iO_,iBlock),& 
                  JouleHeating(iiLon,iilat,iiAlt) * 1.38e-23 *1.5 * &
                  ( NDensityS(iLon,iLat,iAlt,iO_,iBlock) + &
                    NDensityS(iLon,iLat,iAlt,iO2_,iBlock) + &
                    NDensityS(iLon,iLat,iAlt,iN2_,iBlock) ), &
                  JouleHeating(iiLon,iilat,iiAlt) , &
                    VerticalVelocity(iLon,iLat,iAlt,1:nSpecies,iBlock) 
          enddo
       enddo
    enddo

  end subroutine output_3dall

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine output_3dsrc

    nvars_to_write = 15
    write(output_format,"('(1p,',I2,'E11.3)')") nvars_to_write

    if (Is1D) then
       nGCs = 0
    else
       nGCs = 2
    endif

    write(iOutputUnit_,*) "NUMERICAL VALUES"
    write(iOutputUnit_,"(I7,6A)") nvars_to_write, " nvars"
    write(iOutputUnit_,"(I7,7A)") nAlts+4, " nAltitudes"
    write(iOutputUnit_,"(I7,7A)") nLats+nGCs*2, " nLatitude"
    write(iOutputUnit_,"(I7,7A)") nLons+nGCs*2, " nLongitudes"
    write(iOutputUnit_,*) ""

    write(iOutputUnit_,*) "VARIABLE LIST"
    write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
    write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
    write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
    write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Temperature"
    write(iOutputUnit_,"(I7,A1,a)")  5, " ", "EUV Heating"
    write(iOutputUnit_,"(I7,A1,a)")  6, " ", "NO Cooling"
    write(iOutputUnit_,"(I7,A1,a)")  7, " ", "O Cooling"
    write(iOutputUnit_,"(I7,A1,a)")  8, " ", "Auroral Heating"
    write(iOutputUnit_,"(I7,A1,a)")  9, " ", "Ion Precipitation Heating"
    write(iOutputUnit_,"(I7,A1,a)") 10, " ", "Joule Heating"
    write(iOutputUnit_,"(I7,A1,a)") 11, " ", "Conduction"
    write(iOutputUnit_,"(I7,A1,a)") 12, " ", "Chemical Heating"
    write(iOutputUnit_,"(I7,A1,a)") 13, " ", "Vertical Advection"
    write(iOutputUnit_,"(I7,A1,a)") 14, " ", "Horizontal Advection"
    write(iOutputUnit_,"(I7,A1,a)") 15, " ", "TempUnit"

    write(iOutputUnit_,*) ""

    write(iOutputUnit_,*) "BEGIN"

    do iAlt=-1,nAlts+2
       iiAlt = max(min(iAlt,nAlts),1)
       do iLat=1-nGCs,nLats+nGCs
          iiLat = min(max(iLat,1),nLats)
          do iLon=1-nGCs,nLons+nGCs
             iiLon = min(max(iLon,1),nLons)
             write(iOutputUnit_,output_format)       &
                  Longitude(iLon,iBlock), &
                  Latitude(iLat,iBlock), &
                  altitude(iAlt),&
                  Temperature(iLon,iLat,iAlt,iBlock),&
                  EuvHeating(iiLon, iiLat, iiAlt, iBlock)*dt, &
                  -NOCooling(iiLon, iiLat, iiAlt)*dt, &
                  -OCooling(iiLon, iiLat, iiAlt)*dt, &
                  AuroralHeating(iiLon, iiLat, iiAlt)*dt, &
                  IonPrecipHeating(iiLon, iiLat, iiAlt)*dt, &
                  JouleHeating(iiLon, iiLat, iiAlt)*dt, &
                  Conduction(iiLon, iiLat, iiAlt), &
                  ChemicalHeatingRate(iiLon, iiLat, iiAlt), & 
                  VerticalTempSource(iiLon, iiLat, iiAlt), &
                  HorizontalTempSource(iiLon, iiLat, iiAlt), &
                  TempUnit(iLon, iLat, iAlt)

          enddo
       enddo
    enddo

  end subroutine output_3dsrc

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine output_3daur

    nvars_to_write = 9
    write(output_format,"('(1p,',I2,'E11.3)')") nvars_to_write

    if (Is1D) then
       nGCs = 0
    else
       nGCs = 2
    endif

    write(iOutputUnit_,*) "NUMERICAL VALUES"
    write(iOutputUnit_,"(I7,6A)") nvars_to_write, " nvars"
    write(iOutputUnit_,"(I7,7A)") nAlts, " nAltitudes"
    write(iOutputUnit_,"(I7,7A)") nLats+nGCs*2, " nLatitude"
    write(iOutputUnit_,"(I7,7A)") nLons+nGCs*2, " nLongitudes"
    write(iOutputUnit_,*) ""

    write(iOutputUnit_,*) "VARIABLE LIST"
    write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
    write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
    write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
    write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Auroral Ion Rate O"
    write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Auroral Ion Rate O2"
    write(iOutputUnit_,"(I7,A1,a)")  6, " ", "Auroral Ion Rate N2"
    write(iOutputUnit_,"(I7,A1,a)")  7, " ", "Ion Precip Ion Rate O"
    write(iOutputUnit_,"(I7,A1,a)")  8, " ", "Ion Precip Ion Rate O2"
    write(iOutputUnit_,"(I7,A1,a)")  9, " ", "Ion Precip Ion Rate N2"

    write(iOutputUnit_,*) ""

    write(iOutputUnit_,*) "BEGIN"

    do iAlt=1,nAlts
       do iLat=1-nGCs,nLats+nGCs
          iiLat = min(max(iLat,1),nLats)
          do iLon=1-nGCs,nLons+nGCs
             iiLon = min(max(iLon,1),nLons)
             write(iOutputUnit_,output_format)       &
                  Longitude(iLon,iBlock), &
                  Latitude(iLat,iBlock), &
                  altitude(iAlt),&
                  AuroralIonRateS(iLon,iLat,iAlt,iO_,iBlock),&
                  AuroralIonRateS(iLon,iLat,iAlt,iO2_,iBlock),&
                  AuroralIonRateS(iLon,iLat,iAlt,iN2_,iBlock),&
                  IonPrecipIonRateS(iLon,iLat,iAlt,iO_,iBlock),&
                  IonPrecipIonRateS(iLon,iLat,iAlt,iO2_,iBlock),&
                  IonPrecipIonRateS(iLon,iLat,iAlt,iN2_,iBlock)

          enddo
       enddo
    enddo

  end subroutine output_3daur

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine output_3dion

    nvars_to_write = 15
    write(output_format,"('(1p,',I2,'E11.3)')") nvars_to_write

    if (Is1D) then
       nGCs = 0
    else
       nGCs = 2
    endif

    write(iOutputUnit_,*) "NUMERICAL VALUES"
    write(iOutputUnit_,"(I7,6A)") nvars_to_write, " nvars"
    write(iOutputUnit_,"(I7,7A)") nAlts+4, " nAltitudes"
    write(iOutputUnit_,"(I7,7A)") nLats+nGCs*2, " nLatitude"
    write(iOutputUnit_,"(I7,7A)") nLons+nGCs*2, " nLongitudes"
    write(iOutputUnit_,*) ""

    write(iOutputUnit_,*) "VARIABLE LIST"
    write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
    write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
    write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
    write(iOutputUnit_,"(I7,A1,a)")  4, " ", "[O4S+]"
    write(iOutputUnit_,"(I7,A1,a)")  5, " ", "[O2+]"
    write(iOutputUnit_,"(I7,A1,a)")  6, " ", "[N2+]"
    write(iOutputUnit_,"(I7,A1,a)")  7, " ", "[N+]"
    write(iOutputUnit_,"(I7,A1,a)")  8, " ", "[NO+]"
    write(iOutputUnit_,"(I7,A1,a)")  9, " ", "[e-]"
    write(iOutputUnit_,"(I7,A1,a)") 10, " ", "eTemperature"
    write(iOutputUnit_,"(I7,A1,a)") 11, " ", "iTemperature"
    write(iOutputUnit_,"(I7,A1,a)") 12, " ", "Vi (east)"
    write(iOutputUnit_,"(I7,A1,a)") 13, " ", "Vi (north)"
    write(iOutputUnit_,"(I7,A1,a)") 14, " ", "Vi (up)"
    write(iOutputUnit_,"(I7,A1,a)") 15, " ", "Potential (kV)"
    write(iOutputUnit_,*) ""

    write(iOutputUnit_,*) "BEGIN"

    do iAlt=-1,nAlts+2
       iiAlt = max(min(iAlt,nAlts),1)
       do iLat=1-nGCs,nLats+nGCs
          iiLat = min(max(iLat,1),nLats)
          do iLon=1-nGCs,nLons+nGCs
             iiLon = min(max(iLon,1),nLons)
             write(iOutputUnit_,output_format)       &
                  Longitude(iLon,iBlock), &
                  Latitude(iLat,iBlock), &
                  altitude(iAlt),&
                  IDensityS(iLon,iLat,iAlt,iO_4SP_,iBlock), &
                  IDensityS(iLon,iLat,iAlt,iO2P_,iBlock), &
                  IDensityS(iLon,iLat,iAlt,iN2P_,iBlock), &
                  IDensityS(iLon,iLat,iAlt,iNP_,iBlock), &
                  IDensityS(iLon,iLat,iAlt,iNOP_,iBlock), &
                  IDensityS(iLon,iLat,iAlt,ie_,iBlock), &
                  eTemperature(iLon,iLat,iAlt,iBlock),&
                  ITemperature(iLon,iLat,iAlt,iBlock),&
                  Ivelocity(iLon,iLat,iAlt,iEast_,iBlock), &
                  Ivelocity(iLon,iLat,iAlt,iNorth_,iBlock), &
                  Ivelocity(iLon,iLat,iAlt,iUp_,iBlock), &
                  Potential(iiLon,iiLat,iiAlt,iBlock)
          enddo
       enddo
    enddo

  end subroutine output_3dion

end subroutine output

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

subroutine output_1dall(iiLon, iiLat, iBlock, rLon, rLat, iUnit)

  use ModGITM
  use ModSources, only: JouleHeating
  implicit none

  integer, intent(in) :: iiLat, iiLon, iBlock, iUnit
  real, intent(in)    :: rLon, rLat

  integer, parameter :: nVars = 28
  real :: Vars(nVars)
  character (len=100) :: output_format
  integer :: iAlt, iiAlt

  write(output_format,"('(1p,',I2,'E11.3)')") nVars

  write(iUnit,*) "BLOCKS"
  write(iUnit,"(I7,A)") 1, " nBlocksAlt"
  write(iUnit,"(I7,A)") 1, " nBlocksLat"
  write(iUnit,"(I7,A)") 1, " nBlocksLon"
  write(iUnit,*) ""

  write(iUnit,*) "NUMERICAL VALUES"
  write(iUnit,"(I7,6A)") nVars, " nvars"
  write(iUnit,"(I7,7A)") nAlts+4, " nAltitudes"
  write(iUnit,"(I7,7A)") 1, " nLatitude"
  write(iUnit,"(I7,7A)") 1, " nLongitudes"
  write(iUnit,*) ""
  write(iUnit,*) "VARIABLE LIST"
  write(iUnit,"(I7,A1,a)")  1, " ", "Longitude"
  write(iUnit,"(I7,A1,a)")  2, " ", "Latitude"
  write(iUnit,"(I7,A1,a)")  3, " ", "Altitude"
  write(iUnit,"(I7,A1,a)")  4, " ", "Rho"
  write(iUnit,"(I7,A1,a)")  5, " ", "Temperature"
  write(iUnit,"(I7,A1,a)")  6, " ", "[O]"
  write(iUnit,"(I7,A1,a)")  7, " ", "[O2]"
  write(iUnit,"(I7,A1,a)")  8, " ", "[N2]"
  write(iUnit,"(I7,A1,a)")  9, " ", "[NO]"
  write(iUnit,"(I7,A1,a)") 10, " ", "[N_4S]"
  write(iUnit,"(I7,A1,a)") 11, " ", "[N_2D]"
  write(iUnit,"(I7,A1,a)") 12, " ", "Vn (east)"
  write(iUnit,"(I7,A1,a)") 13, " ", "Vn (north)"
  write(iUnit,"(I7,A1,a)") 14, " ", "Vn (up)"
  write(iUnit,"(I7,A1,a)") 15, " ", "[O+]"
  write(iUnit,"(I7,A1,a)") 16, " ", "[N+]"
  write(iUnit,"(I7,A1,a)") 17, " ", "[O2+]"
  write(iUnit,"(I7,A1,a)") 18, " ", "[N2+]"
  write(iUnit,"(I7,A1,a)") 19, " ", "[NO+]"
  write(iUnit,"(I7,A1,a)") 20, " ", "[e-]"
  write(iUnit,"(I7,A1,a)") 21, " ", "eTemperature"
  write(iUnit,"(I7,A1,a)") 22, " ", "iTemperature"
  write(iUnit,"(I7,A1,a)") 23, " ", "Vi (east)"
  write(iUnit,"(I7,A1,a)") 24, " ", "Vi (north)"
  write(iUnit,"(I7,A1,a)") 25, " ", "Vi (up)"
  write(iUnit,"(I7,A1,a)") 26, " ", "Potential (kV)"
  write(iUnit,"(I7,A1,a)") 27, " ", "Heating rate (K/s)" 
  write(iUnit,"(I7,A1,a)") 28, " ", "Heating rate (K/s)"

  write(iUnit,*) ""
  write(iUnit,*) "BEGIN"

  do iAlt=-1,nAlts+2

     iiAlt = max(min(iAlt,nAlts),1)

     Vars(1) = &
          rLon*Longitude(iiLon,iBlock)+(1-rLon)*Longitude(iiLon+1,iBlock)
     Vars(2) = &
          rLat*Latitude(iiLat,iBlock)+(1-rLat)*Latitude(iiLat+1,iBlock)
     Vars(3) = altitude(iAlt)
     Vars(4) = &
          (  rLon)*(  rLat)*Rho(iiLon  ,iiLat  ,iAlt,iBlock) + &
          (1-rLon)*(  rLat)*Rho(iiLon+1,iiLat  ,iAlt,iBlock) + &
          (  rLon)*(1-rLat)*Rho(iiLon  ,iiLat+1,iAlt,iBlock) + &
          (1-rLon)*(1-rLat)*Rho(iiLon+1,iiLat+1,iAlt,iBlock)
     Vars(5) = TempUnit(iiLon,iiLat,iAlt) * ( &
          (  rLon)*(  rLat)*Temperature(iiLon,iiLat,iAlt,iBlock) + &
          (1-rLon)*(  rLat)*Temperature(iiLon+1,iiLat,iAlt,iBlock) + &
          (  rLon)*(1-rLat)*Temperature(iiLon,iiLat+1,iAlt,iBlock) + &
          (1-rLon)*(1-rLat)*Temperature(iiLon+1,iiLat+1,iAlt,iBlock))
     Vars(6) = &
          (  rLon)*(  rLat)*NDensityS(iiLon  ,iiLat  ,iAlt,iO_,iBlock) + &
          (1-rLon)*(  rLat)*NDensityS(iiLon+1,iiLat  ,iAlt,iO_,iBlock) + &
          (  rLon)*(1-rLat)*NDensityS(iiLon  ,iiLat+1,iAlt,iO_,iBlock) + &
          (1-rLon)*(1-rLat)*NDensityS(iiLon+1,iiLat+1,iAlt,iO_,iBlock)
     Vars(7) = &
          (  rLon)*(  rLat)*NDensityS(iiLon  ,iiLat  ,iAlt,iO2_,iBlock) + &
          (1-rLon)*(  rLat)*NDensityS(iiLon+1,iiLat  ,iAlt,iO2_,iBlock) + &
          (  rLon)*(1-rLat)*NDensityS(iiLon  ,iiLat+1,iAlt,iO2_,iBlock) + &
          (1-rLon)*(1-rLat)*NDensityS(iiLon+1,iiLat+1,iAlt,iO2_,iBlock)
     Vars(8) = &
          (  rLon)*(  rLat)*NDensityS(iiLon  ,iiLat  ,iAlt,iN2_,iBlock) + &
          (1-rLon)*(  rLat)*NDensityS(iiLon+1,iiLat  ,iAlt,iN2_,iBlock) + &
          (  rLon)*(1-rLat)*NDensityS(iiLon  ,iiLat+1,iAlt,iN2_,iBlock) + &
          (1-rLon)*(1-rLat)*NDensityS(iiLon+1,iiLat+1,iAlt,iN2_,iBlock)
     Vars(9) = &
          (  rLon)*(  rLat)*NDensityS(iiLon  ,iiLat  ,iAlt,iNO_,iBlock) + &
          (1-rLon)*(  rLat)*NDensityS(iiLon+1,iiLat  ,iAlt,iNO_,iBlock) + &
          (  rLon)*(1-rLat)*NDensityS(iiLon  ,iiLat+1,iAlt,iNO_,iBlock) + &
          (1-rLon)*(1-rLat)*NDensityS(iiLon+1,iiLat+1,iAlt,iNO_,iBlock)
     Vars(10) = &
          (  rLon)*(  rLat)*NDensityS(iiLon  ,iiLat  ,iAlt,iN_4S_,iBlock) + &
          (1-rLon)*(  rLat)*NDensityS(iiLon+1,iiLat  ,iAlt,iN_4S_,iBlock) + &
          (  rLon)*(1-rLat)*NDensityS(iiLon  ,iiLat+1,iAlt,iN_4S_,iBlock) + &
          (1-rLon)*(1-rLat)*NDensityS(iiLon+1,iiLat+1,iAlt,iN_4S_,iBlock)
     Vars(11) = &
          (  rLon)*(  rLat)*NDensityS(iiLon  ,iiLat  ,iAlt,iN_2D_,iBlock) + &
          (1-rLon)*(  rLat)*NDensityS(iiLon+1,iiLat  ,iAlt,iN_2D_,iBlock) + &
          (  rLon)*(1-rLat)*NDensityS(iiLon  ,iiLat+1,iAlt,iN_2D_,iBlock) + &
          (1-rLon)*(1-rLat)*NDensityS(iiLon+1,iiLat+1,iAlt,iN_2D_,iBlock)
     Vars(12) = &
          (  rLon)*(  rLat)*Velocity(iiLon  ,iiLat  ,iAlt,iEast_,iBlock) + &
          (1-rLon)*(  rLat)*Velocity(iiLon+1,iiLat  ,iAlt,iEast_,iBlock) + &
          (  rLon)*(1-rLat)*Velocity(iiLon  ,iiLat+1,iAlt,iEast_,iBlock) + &
          (1-rLon)*(1-rLat)*Velocity(iiLon+1,iiLat+1,iAlt,iEast_,iBlock)
     Vars(13) = &
          (  rLon)*(  rLat)*Velocity(iiLon  ,iiLat  ,iAlt,iNorth_,iBlock) + &
          (1-rLon)*(  rLat)*Velocity(iiLon+1,iiLat  ,iAlt,iNorth_,iBlock) + &
          (  rLon)*(1-rLat)*Velocity(iiLon  ,iiLat+1,iAlt,iNorth_,iBlock) + &
          (1-rLon)*(1-rLat)*Velocity(iiLon+1,iiLat+1,iAlt,iNorth_,iBlock)
     Vars(14) = &
          (  rLon)*(  rLat)*Velocity(iiLon  ,iiLat  ,iAlt,iUp_,iBlock) + &
          (1-rLon)*(  rLat)*Velocity(iiLon+1,iiLat  ,iAlt,iUp_,iBlock) + &
          (  rLon)*(1-rLat)*Velocity(iiLon  ,iiLat+1,iAlt,iUp_,iBlock) + &
          (1-rLon)*(1-rLat)*Velocity(iiLon+1,iiLat+1,iAlt,iUp_,iBlock)
     Vars(15) = &
          (  rLon)*(  rLat)*IDensityS(iiLon  ,iiLat  ,iAlt,iO_4SP_,iBlock) + &
          (1-rLon)*(  rLat)*IDensityS(iiLon+1,iiLat  ,iAlt,iO_4SP_,iBlock) + &
          (  rLon)*(1-rLat)*IDensityS(iiLon  ,iiLat+1,iAlt,iO_4SP_,iBlock) + &
          (1-rLon)*(1-rLat)*IDensityS(iiLon+1,iiLat+1,iAlt,iO_4SP_,iBlock)
     Vars(16) = &
          (  rLon)*(  rLat)*IDensityS(iiLon  ,iiLat  ,iAlt,iNP_,iBlock) + &
          (1-rLon)*(  rLat)*IDensityS(iiLon+1,iiLat  ,iAlt,iNP_,iBlock) + &
          (  rLon)*(1-rLat)*IDensityS(iiLon  ,iiLat+1,iAlt,iNP_,iBlock) + &
          (1-rLon)*(1-rLat)*IDensityS(iiLon+1,iiLat+1,iAlt,iNP_,iBlock)
     Vars(17) = &
          (  rLon)*(  rLat)*IDensityS(iiLon  ,iiLat  ,iAlt,iO2P_,iBlock) + &
          (1-rLon)*(  rLat)*IDensityS(iiLon+1,iiLat  ,iAlt,iO2P_,iBlock) + &
          (  rLon)*(1-rLat)*IDensityS(iiLon  ,iiLat+1,iAlt,iO2P_,iBlock) + &
          (1-rLon)*(1-rLat)*IDensityS(iiLon+1,iiLat+1,iAlt,iO2P_,iBlock)
     Vars(18) = &
          (  rLon)*(  rLat)*IDensityS(iiLon  ,iiLat  ,iAlt,iN2P_,iBlock) + &
          (1-rLon)*(  rLat)*IDensityS(iiLon+1,iiLat  ,iAlt,iN2P_,iBlock) + &
          (  rLon)*(1-rLat)*IDensityS(iiLon  ,iiLat+1,iAlt,iN2P_,iBlock) + &
          (1-rLon)*(1-rLat)*IDensityS(iiLon+1,iiLat+1,iAlt,iN2P_,iBlock)
     Vars(19) = &
          (  rLon)*(  rLat)*IDensityS(iiLon  ,iiLat  ,iAlt,iNO_,iBlock) + &
          (1-rLon)*(  rLat)*IDensityS(iiLon+1,iiLat  ,iAlt,iNO_,iBlock) + &
          (  rLon)*(1-rLat)*IDensityS(iiLon  ,iiLat+1,iAlt,iNO_,iBlock) + &
          (1-rLon)*(1-rLat)*IDensityS(iiLon+1,iiLat+1,iAlt,iNO_,iBlock)
     Vars(20) = &
          (  rLon)*(  rLat)*IDensityS(iiLon  ,iiLat  ,iAlt,ie_,iBlock) + &
          (1-rLon)*(  rLat)*IDensityS(iiLon+1,iiLat  ,iAlt,ie_,iBlock) + &
          (  rLon)*(1-rLat)*IDensityS(iiLon  ,iiLat+1,iAlt,ie_,iBlock) + &
          (1-rLon)*(1-rLat)*IDensityS(iiLon+1,iiLat+1,iAlt,ie_,iBlock)
     Vars(21) = &
          (  rLon)*(  rLat)*eTemperature(iiLon,iiLat,iAlt,iBlock) + &
          (1-rLon)*(  rLat)*eTemperature(iiLon+1,iiLat,iAlt,iBlock) + &
          (  rLon)*(1-rLat)*eTemperature(iiLon,iiLat+1,iAlt,iBlock) + &
          (1-rLon)*(1-rLat)*eTemperature(iiLon+1,iiLat+1,iAlt,iBlock)
     Vars(22) = &
          (  rLon)*(  rLat)*ITemperature(iiLon,iiLat,iAlt,iBlock) + &
          (1-rLon)*(  rLat)*ITemperature(iiLon+1,iiLat,iAlt,iBlock) + &
          (  rLon)*(1-rLat)*ITemperature(iiLon,iiLat+1,iAlt,iBlock) + &
          (1-rLon)*(1-rLat)*ITemperature(iiLon+1,iiLat+1,iAlt,iBlock)
     Vars(23) = &
          (  rLon)*(  rLat)*IVelocity(iiLon  ,iiLat  ,iAlt,iEast_,iBlock) + &
          (1-rLon)*(  rLat)*IVelocity(iiLon+1,iiLat  ,iAlt,iEast_,iBlock) + &
          (  rLon)*(1-rLat)*IVelocity(iiLon  ,iiLat+1,iAlt,iEast_,iBlock) + &
          (1-rLon)*(1-rLat)*IVelocity(iiLon+1,iiLat+1,iAlt,iEast_,iBlock)
     Vars(24) = &
          (  rLon)*(  rLat)*IVelocity(iiLon  ,iiLat  ,iAlt,iNorth_,iBlock) + &
          (1-rLon)*(  rLat)*IVelocity(iiLon+1,iiLat  ,iAlt,iNorth_,iBlock) + &
          (  rLon)*(1-rLat)*IVelocity(iiLon  ,iiLat+1,iAlt,iNorth_,iBlock) + &
          (1-rLon)*(1-rLat)*IVelocity(iiLon+1,iiLat+1,iAlt,iNorth_,iBlock)
     Vars(25) = &
          (  rLon)*(  rLat)*IVelocity(iiLon  ,iiLat  ,iAlt,iUp_,iBlock) + &
          (1-rLon)*(  rLat)*IVelocity(iiLon+1,iiLat  ,iAlt,iUp_,iBlock) + &
          (  rLon)*(1-rLat)*IVelocity(iiLon  ,iiLat+1,iAlt,iUp_,iBlock) + &
          (1-rLon)*(1-rLat)*IVelocity(iiLon+1,iiLat+1,iAlt,iUp_,iBlock)
     Vars(26) = &
          (  rLon)*(  rLat)*Potential(iiLon,iiLat,iiAlt,iBlock) + &
          (1-rLon)*(  rLat)*Potential(iiLon+1,iiLat,iiAlt,iBlock) + &
          (  rLon)*(1-rLat)*Potential(iiLon,iiLat+1,iiAlt,iBlock) + &
          (1-rLon)*(1-rLat)*Potential(iiLon+1,iiLat+1,iiAlt,iBlock)


   Vars(27) = &  
          (  rLon)*(  rLat)*JouleHeating(iiLon,iiLat,iiAlt) + &                    
          (1-rLon)*(  rLat)*JouleHeating(iiLon+1,iiLat,iiAlt) + &                  
          (  rLon)*(1-rLat)*JouleHeating(iiLon,iiLat+1,iiAlt) + &                  
          (1-rLon)*(1-rLat)*JouleHeating(iiLon+1,iiLat+1,iiAlt)                    

   Vars(28) = &  
          (  rLon)*(  rLat)*JouleHeating(iiLon,iiLat,iiAlt) + &                    
          (1-rLon)*(  rLat)*JouleHeating(iiLon+1,iiLat,iiAlt) + &                  
          (  rLon)*(1-rLat)*JouleHeating(iiLon,iiLat+1,iiAlt) + &                  
          (1-rLon)*(1-rLat)*JouleHeating(iiLon+1,iiLat+1,iiAlt)                    




     write(iUnit,output_format) Vars

  enddo

end subroutine output_1dall

