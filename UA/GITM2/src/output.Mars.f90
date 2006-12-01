
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

     if (OutputType(iOutputType) == '3DALL')     IsFound = .true.
     if (OutputType(iOutputType) == '1DALL')     IsFound = .true.
     if (OutputType(iOutputType) == '3DNEUTRAL') IsFound = .true.
     if (OutputType(iOutputType) == '3DION')     IsFound = .true.
     if (OutputType(iOutputType) == '3DELECTRO') IsFound = .true.
     if (OutputType(iOutputType) == '2D')        IsFound = .true.

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
       write(*,*) "Writing Output files at iStep : ",iStep, tSimulation

  call calc_physics(iBlock)
  call chapman_integrals(iBlock)
  call calc_rates(iBlock)
!  call calc_efield(iBlock)

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

!!  case ('3DSRC')
!!
!!     !! open file
!!     open(unit=iOutputUnit_, &
!!          file=dir//"/"//blk_str//"_"//cTime//"."//OutputType(iOutputType),&
!!          status="unknown")
!!
!!     call write_head_blocks
!!     call write_head_time
!!     call write_head_version
!!     call output_3dsrc
!!
!!     !! close file
!!     close(unit=iOutputUnit_)
!!
!!  case ('3DAUR')
!!
!!     !! open file
!!     open(unit=iOutputUnit_, &
!!          file=dir//"/"//blk_str//"_"//cTime//"."//OutputType(iOutputType),&
!!          status="unknown")
!!
!!     call write_head_blocks
!!     call write_head_time
!!     call write_head_version
!!     call output_3daur
!!
!!     !! close file
!!     close(unit=iOutputUnit_)
!!
!!  case ('2DGEO')
!!
!!     !! open file
!!     open(unit=iOutputUnit_, &
!!          file=dir//"/"//blk_str//"_"//cTime//"."//OutputType(iOutputType),&
!!          status="unknown")
!!
!!     call write_head_blocks
!!     call write_head_time
!!     call write_head_version
!!     call output_2d
!!
!!     !! close file
!!     close(unit=iOutputUnit_)
!!
!!  case ('2DMAG')
!!
!!     if (iBLK == 1) then
!!
!!        !! open file
!!        open(unit=iOutputUnit_, &
!!             file=dir//"/"//blk_str//"_"//cTime//"."//OutputType(iOutputType),&
!!             status="unknown")
!!
!!        call write_head_time
!!        call write_head_version
!!        call output_2delectro
!!
!!        !! close file
!!        close(unit=iOutputUnit_)
!!
!!     endif

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

!!  case ('3DION')
!!
!!     !! open file
!!     open(unit=iOutputUnit_, &
!!          file=dir//"/"//blk_str//"_"//cTime//"."//OutputType(iOutputType),&
!!          status="unknown")
!!
!!     call write_head_blocks
!!     call write_head_time
!!     call write_head_version
!!     call output_3dion
!!
!!     !! close file
!!     close(unit=iOutputUnit_)

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

    nvars_to_write = 15
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
    write(iOutputUnit_,"(I7,A1,a)")  7, " ", "Field Line Length"
    write(iOutputUnit_,"(I7,A1,a)")  8, " ", "Sigma PP"
    write(iOutputUnit_,"(I7,A1,a)")  9, " ", "Sigma LL"
    write(iOutputUnit_,"(I7,A1,a)") 10, " ", "Sigma H"
    write(iOutputUnit_,"(I7,A1,a)") 11, " ", "Sigma C"
    write(iOutputUnit_,"(I7,A1,a)") 12, " ", "Sigma PL"
    write(iOutputUnit_,"(I7,A1,a)") 13, " ", "Sigma LP"
    write(iOutputUnit_,"(I7,A1,a)") 14, " ", "K^D_{m\phi}"
    write(iOutputUnit_,"(I7,A1,a)") 15, " ", "K^D_{m\lamda}"
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
               LengthMC(iLon,iLat), &
               SigmaHHMC(iLon,iLat), &
               SigmaLLMC(iLon,iLat), &
               SigmaHHMC(iLon,iLat), &
               SigmaCCMC(iLon,iLat), &
               SigmaLPMC(iLon,iLat), &
               SigmaPLMC(iLon,iLat), &
               KDmpMC(iLon,iLat), &
               KdmlMC(iLon,iLat)
       enddo
    enddo

  end subroutine output_2delectro

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine output_3dall

    nvars_to_write = 18
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
    write(iOutputUnit_,"(I7,A1,a)")  4, " ", "LocalTime"
    write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Rho"
    write(iOutputUnit_,"(I7,A1,a)")  6, " ", "Temperature"
    write(iOutputUnit_,"(I7,A1,a)")  7, " ", "[CO]"
    write(iOutputUnit_,"(I7,A1,a)")  8, " ", "[CO2]"
    write(iOutputUnit_,"(I7,A1,a)")  9, " ", "[O]"
    write(iOutputUnit_,"(I7,A1,a)") 10, " ", "[N2]"
    write(iOutputUnit_,"(I7,A1,a)") 11, " ", "Vn (east)"
    write(iOutputUnit_,"(I7,A1,a)") 12, " ", "Vn (north)"
    write(iOutputUnit_,"(I7,A1,a)") 13, " ", "Vn (up)"
    write(iOutputUnit_,"(I7,A1,a)") 14, " ", "[O2+]"
    write(iOutputUnit_,"(I7,A1,a)") 15, " ", "[CO2+]"
    write(iOutputUnit_,"(I7,A1,a)") 16, " ", "[NO+]"
    write(iOutputUnit_,"(I7,A1,a)") 17, " ", "[O+]"
    write(iOutputUnit_,"(I7,A1,a)") 18, " ", "[e-]"
    write(iOutputUnit_,*) ""

    write(iOutputUnit_,*) "BEGIN"

    do iAlt=-1,nAlts+2
       iiAlt = max(min(iAlt,nAlts),1)
       do iLat=1-nGCs,nLats+nGCs
          iiLat = min(max(iLat,1),nLats)
          do iLon=1-nGCs,nLons+nGCs
             iiLon = min(max(iLon,0),nLons+1)
             write(iOutputUnit_,output_format)       &
                  Longitude(iLon,iBlock), &
                  Latitude(iLat,iBlock), &
                  altitude(iAlt),&
                  LocalTime(iiLon), &
                  Rho(iLon,iLat,iAlt,iBlock),&
                  Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt),&
                  NDensityS(iLon,iLat,iAlt,iCO_,iBlock), &
                  NDensityS(iLon,iLat,iAlt,iCO2_,iBlock), &
                  NDensityS(iLon,iLat,iAlt,iO_,iBlock), &
                  NDensityS(iLon,iLat,iAlt,iN2_,iBlock), &
                  velocity(iLon,iLat,iAlt,iEast_,iBlock), &
                  velocity(iLon,iLat,iAlt,iNorth_,iBlock), &
                  velocity(iLon,iLat,iAlt,iUp_,iBlock), &
                  IDensityS(iLon,iLat,iAlt,iO2_,iBlock), &
                  IDensityS(iLon,iLat,iAlt,iCO2P_,iBlock), &
                  IDensityS(iLon,iLat,iAlt,iNOP_,iBlock), &
                  IDensityS(iLon,iLat,iAlt,iOP_,iBlock), &
                  IDensityS(iLon,iLat,iAlt,ie_,iBlock)
          enddo
       enddo
    enddo

  end subroutine output_3dall

end subroutine output

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

subroutine output_1dall(iiLon, iiLat, iBlock, rLon, rLat, iUnit)

  use ModGITM
  use ModPlanet

  implicit none

  integer, intent(in) :: iiLat, iiLon, iBlock, iUnit
  real, intent(in)    :: rLon, rLat

  integer, parameter :: nVars = 12
  real :: Vars(nVars)
  character (len=100) :: output_format
  integer :: iAlt, iiAlt

  write(output_format,"('(1p,',I2,'E11.3)')") nVars

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
  write(iUnit,"(I7,A1,a)")  6, " ", "[CO]"
  write(iUnit,"(I7,A1,a)")  7, " ", "[CO2]"
  write(iUnit,"(I7,A1,a)")  8, " ", "[O]"
  write(iUnit,"(I7,A1,a)")  9, " ", "[N2]"
  write(iUnit,"(I7,A1,a)") 10, " ", "Vn (east)"
  write(iUnit,"(I7,A1,a)") 11, " ", "Vn (north)"
  write(iUnit,"(I7,A1,a)") 12, " ", "Vn (up)"
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
          (  rLon)*(  rLat)*Rho(iiLon,iiLat,iAlt,iBlock) + &
          (1-rLon)*(  rLat)*Rho(iiLon+1,iiLat,iAlt,iBlock) + &
          (  rLon)*(1-rLat)*Rho(iiLon,iiLat+1,iAlt,iBlock) + &
          (1-rLon)*(1-rLat)*Rho(iiLon+1,iiLat+1,iAlt,iBlock)
     Vars(5) = TempUnit(iiLon,iiLat,iAlt) * ( &
          (  rLon)*(  rLat)*Temperature(iiLon,iiLat,iAlt,iBlock) + &
          (1-rLon)*(  rLat)*Temperature(iiLon+1,iiLat,iAlt,iBlock) + &
          (  rLon)*(1-rLat)*Temperature(iiLon,iiLat+1,iAlt,iBlock) + &
          (1-rLon)*(1-rLat)*Temperature(iiLon+1,iiLat+1,iAlt,iBlock))
     Vars(6) = &
          (  rLon)*(  rLat)*NDensityS(iiLon  ,iiLat  ,iAlt,iCO_,iBlock) + &
          (1-rLon)*(  rLat)*NDensityS(iiLon+1,iiLat  ,iAlt,iCO_,iBlock) + &
          (  rLon)*(1-rLat)*NDensityS(iiLon  ,iiLat+1,iAlt,iCO_,iBlock) + &
          (1-rLon)*(1-rLat)*NDensityS(iiLon+1,iiLat+1,iAlt,iCO_,iBlock)
     Vars(7) = &
          (  rLon)*(  rLat)*NDensityS(iiLon  ,iiLat  ,iAlt,iCO2_,iBlock) + &
          (1-rLon)*(  rLat)*NDensityS(iiLon+1,iiLat  ,iAlt,iCO2_,iBlock) + &
          (  rLon)*(1-rLat)*NDensityS(iiLon  ,iiLat+1,iAlt,iCO2_,iBlock) + &
          (1-rLon)*(1-rLat)*NDensityS(iiLon+1,iiLat+1,iAlt,iCO2_,iBlock)
     Vars(8) = &
          (  rLon)*(  rLat)*NDensityS(iiLon  ,iiLat  ,iAlt,iO_,iBlock) + &
          (1-rLon)*(  rLat)*NDensityS(iiLon+1,iiLat  ,iAlt,iO_,iBlock) + &
          (  rLon)*(1-rLat)*NDensityS(iiLon  ,iiLat+1,iAlt,iO_,iBlock) + &
          (1-rLon)*(1-rLat)*NDensityS(iiLon+1,iiLat+1,iAlt,iO_,iBlock)
     Vars(9) = &
          (  rLon)*(  rLat)*NDensityS(iiLon  ,iiLat  ,iAlt,iN2_,iBlock) + &
          (1-rLon)*(  rLat)*NDensityS(iiLon+1,iiLat  ,iAlt,iN2_,iBlock) + &
          (  rLon)*(1-rLat)*NDensityS(iiLon  ,iiLat+1,iAlt,iN2_,iBlock) + &
          (1-rLon)*(1-rLat)*NDensityS(iiLon+1,iiLat+1,iAlt,iN2_,iBlock)
     Vars(10) = &
          (  rLon)*(  rLat)*Velocity(iiLon  ,iiLat  ,iAlt,iNorth_,iBlock) + &
          (1-rLon)*(  rLat)*Velocity(iiLon+1,iiLat  ,iAlt,iNorth_,iBlock) + &
          (  rLon)*(1-rLat)*Velocity(iiLon  ,iiLat+1,iAlt,iNorth_,iBlock) + &
          (1-rLon)*(1-rLat)*Velocity(iiLon+1,iiLat+1,iAlt,iNorth_,iBlock)
     Vars(11) = &
          (  rLon)*(  rLat)*Velocity(iiLon  ,iiLat  ,iAlt,iEast_,iBlock) + &
          (1-rLon)*(  rLat)*Velocity(iiLon+1,iiLat  ,iAlt,iEast_,iBlock) + &
          (  rLon)*(1-rLat)*Velocity(iiLon  ,iiLat+1,iAlt,iEast_,iBlock) + &
          (1-rLon)*(1-rLat)*Velocity(iiLon+1,iiLat+1,iAlt,iEast_,iBlock)
     Vars(12) = &
          (  rLon)*(  rLat)*Velocity(iiLon  ,iiLat  ,iAlt,iUp_,iBlock) + &
          (1-rLon)*(  rLat)*Velocity(iiLon+1,iiLat  ,iAlt,iUp_,iBlock) + &
          (  rLon)*(1-rLat)*Velocity(iiLon  ,iiLat+1,iAlt,iUp_,iBlock) + &
          (1-rLon)*(1-rLat)*Velocity(iiLon+1,iiLat+1,iAlt,iUp_,iBlock)

     write(iUnit,output_format) Vars

  enddo

end subroutine output_1dall

