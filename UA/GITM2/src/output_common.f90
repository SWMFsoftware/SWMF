integer function bad_outputtype()

  use ModInputs, only : OutputType, nOutputTypes

  implicit none

  integer :: iOutputType
  logical :: IsFound

  do iOutputType = 1, nOutputTypes

     IsFound = .false.

     if (OutputType(iOutputType) == '3DALL')     IsFound = .true.
     if (OutputType(iOutputType) == '3DNEU')     IsFound = .true.
     if (OutputType(iOutputType) == '3DION')     IsFound = .true.
     if (OutputType(iOutputType) == '3DTHM')     IsFound = .true.
     if (OutputType(iOutputType) == '3DCHM')     IsFound = .true.
     if (OutputType(iOutputType) == '3DUSR')     IsFound = .true.
     if (OutputType(iOutputType) == '3DGLO')     IsFound = .true.

     if (OutputType(iOutputType) == '2DGEL')     IsFound = .true.
     if (OutputType(iOutputType) == '2DMEL')     IsFound = .true.
     if (OutputType(iOutputType) == '2DUSR')     IsFound = .true.

     if (OutputType(iOutputType) == '1DALL')     IsFound = .true.
     if (OutputType(iOutputType) == '1DGLO')     IsFound = .true.
     if (OutputType(iOutputType) == '1DTHM')     IsFound = .true.
     if (OutputType(iOutputType) == '1DNEW')     IsFound = .true.

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

subroutine output(dir, iBlock, iOutputType)

  use ModSphereInterface, only:iStartBlk
  use ModSatellites, only : CurrentSatellitePosition, CurrentSatelliteName
  use ModGITM
  use ModEUV
  use ModTime
  use ModInputs
  use ModSources
  use ModUserGITM, only: nVarsUser2d, nVarsUser3d

  implicit none

  character (len=*), intent(in) :: dir
  integer, intent(in) :: iBlock
  integer, intent(in) :: iOutputType

  character (len=5) :: proc_str,cBlock, cType
  character (len=14) :: cTime
  integer :: iiLat, iiLon, nGCs
  integer :: iLon,iLat,iAlt, nVars_to_Write, nlines, iBLK,iSpecies
  logical :: done

  real :: LatFind, LonFind
  real :: rLon, rLat

  character (len=2) :: cYear, cMonth, cDay, cHour, cMinute, cSecond

  !! construct naming strings

    if (iOutputType == -1) then
     cType = "1DALL"
  else
     cType = OutputType(iOutputType)
     if (cType(1:2) == "3D" .and. Is1D) then 
        cType(1:2) = "1D"
     endif
  endif

  if (Is1D) then
     iiLat = 1
     iiLon = 1
     rLon = 1.0
     rLat = 1.0
  endif

  if (iOutputType == -1) then

     LatFind = CurrentSatellitePosition(iNorth_)
     LonFind = CurrentSatellitePosition(iEast_)

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

  if ((iProc == 0.and.iBlock == 1).and.(iOutputType /= -1)) &
       write(*,'(a,i7,i5,5i3)') &
       "Writing Output files ("//cType//") at iStep : ",&
       iStep, iTimeArray(1:6)

  if (iOutputType == -1) &
       write(*,'(a,i7,i5,5i3)') &
       "Writing satellite file ("//CurrentSatelliteName//") at iStep : ",&
       iStep, iTimeArray(1:6)

  call calc_physics(iBlock)
  call calc_rates(iBlock)
  call calc_collisions(iBlock)
  call chapman_integrals(iBlock)
  if (.not. Is1D) call calc_efield(iBlock)

  iBLK = iStartBLK + iBlock

  write(cBlock,'(a1,i4.4)') "b",iBLK

  call i2s(mod(iTimeArray(1),100), cYear, 2)
  call i2s(iTimeArray(2), cMonth, 2)
  call i2s(iTimeArray(3), cDay, 2)
  call i2s(iTimeArray(4), cHour, 2)
  call i2s(iTimeArray(5), cMinute, 2)
  call i2s(iTimeArray(6), cSecond, 2)

  cTime = "t"//cYear//cMonth//cDay//"_"//cHour//cMinute//"00"

  !! ---------------------------------------------
  !! Write the binary data files
  !! ---------------------------------------------

  if (iOutputType == -1) then
     open(unit=iOutputUnit_, form="unformatted", &
          file=dir//"/"//CurrentSatelliteName//"_"//cTime//".sat",&
          status="unknown")
  else
     if (cType /= '2DMEL' .or. iBLK == 1) then
        open(unit=iOutputUnit_, form="unformatted", &
             file=dir//"/"//cType//"_"//cTime//"."//cBlock,&
             status="unknown")
     endif
  endif

  nGCs = 2

  select case (cType)

  case ('3DALL')

     nvars_to_write = 13+nSpeciesTotal+nSpecies+nIons
     call output_3dall(iBlock)

  case ('3DNEU')

     nvars_to_write = 8+nSpeciesTotal+nSpecies
     call output_3dneu(iBlock)

  case ('3DION')

     nvars_to_write = 8+nIons+6
     call output_3dion(iBlock)

  case ('3DTHM')

     nvars_to_write = 14
     call output_3dthm(iBlock)

  case ('3DCHM')

     nvars_to_write = 30
     call output_3dchm(iBlock)

  case ('3DGLO')

     nvars_to_write = 3 + 3
     call output_3dglo(iBlock)

  case ('3DUSR')

     if (iBlock == 1) call set_nVarsUser3d
     nvars_to_write = nVarsUser3d
     call output_3duser(iBlock, iOutputUnit_)

  case ('2DGEL')

     nvars_to_write = 13
     call output_2dgel(iBlock)

  case ('2DMEL')

     nvars_to_write = 25
     if (iBLK == 1) call output_2dmel(iBlock)

  case ('2DUSR')

     if (iBlock == 1) call set_nVarsUser2d
     nvars_to_write = nVarsUser2d
     call output_2duser(iBlock, iOutputUnit_)

  case ('1DALL')

     nGCs = 0
     nvars_to_write = 13+nSpeciesTotal+nSpecies+nIons+nSpecies+5
     call output_1dall(iiLon, iiLat, iBlock, rLon, rLat, iOutputUnit_)

  case ('1DGLO')

     nGCs = 0
     nvars_to_write = 6
     call output_1dglo

 case ('1DTHM')
     
     nGCs = 0
     nvars_to_write = 14 + (nspeciestotal*2)
     call output_1dthm

  case ('1DNEW')
     nGCs = 0
     nvars_to_write = 15 + nSpeciesTotal + nSpecies + nIons + nSpecies
     call output_1dnew(iiLon, iiLat, iBlock, rLon, rLat, iOutputUnit_)

  end select

  close(unit=iOutputUnit_)

  !! Now write the header file

  if ((iProc == 0 .and. iBlock == nBlocks) .or. iOutputType == -1) then 

     if (iOutputType == -1) then
        open(unit=iOutputUnit_, &
             file=dir//"/"//CurrentSatelliteName//"_"//cTime//".header",&
             status="unknown")
     else
        open(unit=iOutputUnit_, &
             file=dir//"/"//cType//"_"//cTime//".header",&
             status="unknown")
     endif

     call write_head_blocks
     call write_head_time
     call write_head_version

     if (cType(3:5) == 'USR') then
        call output_header_user(cType, iOutputUnit_)
     elseif (cType(3:5) == 'NEW') then
        call output_header_new
     else
        call output_header
     endif

     close(unit=iOutputUnit_)

  endif

contains

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine output_header

    use ModElectrodynamics, only : nMagLats, nMagLons

    integer :: iOff, iSpecies, iIon

    write(iOutputUnit_,*) "NUMERICAL VALUES"

    write(iOutputUnit_,"(I7,6A)") nvars_to_write, " nvars"
    if (cType(1:2) /= "2D") then 
       write(iOutputUnit_,"(I7,7A)") nAlts+4, " nAltitudes"
    else
       write(iOutputUnit_,"(I7,7A)") 1, " nAltitudes"
    endif
    if (cType(1:2) == "1D") then 
       write(iOutputUnit_,"(I7,7A)") 1, " nLatitudes"
       write(iOutputUnit_,"(I7,7A)") 1, " nLongitudes"
    else
       if (cType(3:5) =="MEL") then
          write(iOutputUnit_,"(I7,A)") nMagLats, " nLatitude"
          write(iOutputUnit_,"(I7,A)") nMagLons+1, " nLongitudes"
          write(iOutputUnit_,*) " "
          write(iOutputUnit_,*) "NO GHOSTCELLS"
       elseif (cType(3:5) =="GEL") then
          write(iOutputUnit_,"(I7,A)") nLats, " nLatitude"
          write(iOutputUnit_,"(I7,A)") nLons, " nLongitudes"
          write(iOutputUnit_,*) " "
          write(iOutputUnit_,*) "NO GHOSTCELLS"
       else
          write(iOutputUnit_,"(I7,7A)") nLats+nGCs*2, " nLatitudes"
          write(iOutputUnit_,"(I7,7A)") nLons+nGCs*2, " nLongitudes"
       endif
    endif
    write(iOutputUnit_,*) ""

    write(iOutputUnit_,*) "VARIABLE LIST"
    write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
    write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
    write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"

    if (cType(3:5) == "GEL") then

       write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Potential"
       write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Pedersen Conductance"
       write(iOutputUnit_,"(I7,A1,a)")  6, " ", "Hall Conductance"
       write(iOutputUnit_,"(I7,A1,a)")  7, " ", "Electron_Average_Energy"
       write(iOutputUnit_,"(I7,A1,a)")  8, " ", "Electron_Energy_Flux"
       write(iOutputUnit_,"(I7,A1,a)")  9, " ", "DivJuAlt"
       write(iOutputUnit_,"(I7,A1,a)") 10, " ", "Pedersen FL Conductance"
       write(iOutputUnit_,"(I7,A1,a)") 11, " ", "Hall FL Conductance"
       write(iOutputUnit_,"(I7,A1,a)") 12, " ", "DivJu FL"
       write(iOutputUnit_,"(I7,A1,a)") 13, " ", "FL Length"

    endif
    
    if (cType(3:5) == "THM") then

       write(iOutputUnit_,"(I7,A1,a)")  4, " ", "EUV Heating"
       write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Conduction"
       write(iOutputUnit_,"(I7,A1,a)")  6, " ", "Molecular Conduction"
       write(iOutputUnit_,"(I7,A1,a)")  7, " ", "Eddy Conduction"
       write(iOutputUnit_,"(I7,A1,a)")  8, " ", "Eddy Adiabatic Conduction"
       write(iOutputUnit_,"(I7,A1,a)")  9, " ", "Chemical Heating"
       write(iOutputUnit_,"(I7,A1,a)")  10, " ", "Auroral Heating"
       write(iOutputUnit_,"(I7,A1,a)")  11, " ", "Joule Heating"
       write(iOutputUnit_,"(I7,A1,a)")  12, " ", "NO Cooling"
       write(iOutputUnit_,"(I7,A1,a)")  13, " ", "O Cooling"
       write(iOutputUnit_,"(I7,A1,a)")  14, " ", "Total Abs EUV"
       if (cType(1:2) == "1D") then
          do iSpecies = 1, nSpeciesTotal
             write(iOutputUnit_,"(I7,A1,a,a)") 11 + iSpecies, " ", &
                  "Production Rate ",cSpecies(iSpecies)
          enddo
          do iSpecies = 1, nSpeciesTotal
             write(iOutputUnit_,"(I7,A1,a,a)") 11 + nSpeciesTotal + iSpecies, " ", &
                  "Loss Rate ",cSpecies(iSpecies)
             
          enddo
       endif
       
    endif

    if (cType(3:5) == "CHM") then

       write(iOutputUnit_,"(I7,A1,a)") 4, " ", "N!D2!U+!N + e"
       write(iOutputUnit_,"(I7,A1,a)") 5, " ", "O!D2!U+!N + e"
       write(iOutputUnit_,"(I7,A1,a)") 6, " ", "N!D2!U+!N + O"
       write(iOutputUnit_,"(I7,A1,a)") 7, " ", "NO!U+!N + e"
       write(iOutputUnit_,"(I7,A1,a)") 8, " ", "N!U+!N + O!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 9, " ", "NO + N"
       write(iOutputUnit_,"(I7,A1,a)") 10, " ","O!U+!N + O!D2!N" 
       write(iOutputUnit_,"(I7,A1,a)") 11, " ", "N + O!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 12, " ", "O!D2!U+!N + N"
       write(iOutputUnit_,"(I7,A1,a)") 13, " ", "O!D2!U+!N + NO"
       write(iOutputUnit_,"(I7,A1,a)") 14, " ", "O!D2!U+!N + N2"
       write(iOutputUnit_,"(I7,A1,a)") 15, " ", "N!D2!U+!N + O!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 16, " ", "N!U+!N + O"
       write(iOutputUnit_,"(I7,A1,a)") 17, " ", "O!+!N + N!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 18, " ", "O(1D) + N!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 19, " ", "O(1D) + O!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 20, " ", "O(1D) + O"
       write(iOutputUnit_,"(I7,A1,a)") 21, " ", "O(1D) + e"
       write(iOutputUnit_,"(I7,A1,a)") 22, " ", "N(2D) + O!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 23, " ", "O!U+!N(2D)+e"
       write(iOutputUnit_,"(I7,A1,a)") 24, " ", "N(2D) + O"
       write(iOutputUnit_,"(I7,A1,a)") 25, " ", "N(2D) + e"
       write(iOutputUnit_,"(I7,A1,a)") 26, " ", "O!U+!N(2D + N!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 27, " ", "O!U+!N(2P) + e"
       write(iOutputUnit_,"(I7,A1,a)") 28, " ", "O!U+!N(2P) + O"
       write(iOutputUnit_,"(I7,A1,a)") 29, " ", "O!U+!N(2P) + N!D2!N"
       write(iOutputUnit_,"(I7,A1,a)") 30, " ", "Chemical Heating Rate"
       
       
    endif

    if (cType(3:5) == "GLO") then

       write(iOutputUnit_,"(I7,A1,a)")  4, " ", "6300 A Emission"
       write(iOutputUnit_,"(I7,A1,a)")  5, " ", "PhotoElectronUp"
       write(iOutputUnit_,"(I7,A1,a)")  6, " ", "PhotoElectronDown"

    endif
       
    if (cType(3:5) == "MEL") then

       write(iOutputUnit_,"(I7,A1,a)")  4, " ", "MLT"
       write(iOutputUnit_,"(I7,A1,a)")  5, " ", "GeoLat"
       write(iOutputUnit_,"(I7,A1,a)")  6, " ", "GeoLon"
       write(iOutputUnit_,"(I7,A1,a)")  7, " ", "Pedersen Conductance"
       write(iOutputUnit_,"(I7,A1,a)")  8, " ", "Hall Conductance"
       write(iOutputUnit_,"(I7,A1,a)")  9, " ", "DivJuAlt"
       write(iOutputUnit_,"(I7,A1,a)") 10, " ", "Field Line Length"
       write(iOutputUnit_,"(I7,A1,a)") 11, " ", "Sigma PP"
       write(iOutputUnit_,"(I7,A1,a)") 12, " ", "Sigma LL"
       write(iOutputUnit_,"(I7,A1,a)") 13, " ", "Sigma H"
       write(iOutputUnit_,"(I7,A1,a)") 14, " ", "Sigma C"
       write(iOutputUnit_,"(I7,A1,a)") 15, " ", "Sigma PL"
       write(iOutputUnit_,"(I7,A1,a)") 16, " ", "Sigma LP"
       write(iOutputUnit_,"(I7,A1,a)") 17, " ", "K^D_{m\phi}"
       write(iOutputUnit_,"(I7,A1,a)") 18, " ", "K^D_{m\lamda}"
       write(iOutputUnit_,"(I7,A1,a)") 19, " ", "Solver A"
       write(iOutputUnit_,"(I7,A1,a)") 20, " ", "Solver B"
       write(iOutputUnit_,"(I7,A1,a)") 21, " ", "Solver C"
       write(iOutputUnit_,"(I7,A1,a)") 22, " ", "Solver D"
       write(iOutputUnit_,"(I7,A1,a)") 23, " ", "Solver E"
       write(iOutputUnit_,"(I7,A1,a)") 24, " ", "Solver S"
       write(iOutputUnit_,"(I7,A1,a)") 25, " ", "DynamoPotential"

    endif

    if (cType(3:5) == "ALL" .or. cType(3:5) == "NEU") then

       write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Rho"
   
       iOff = 4
       do iSpecies = 1, nSpeciesTotal
          write(iOutputUnit_,"(I7,A1,a)")  iOff+iSpecies, " ", &
               "["//cSpecies(iSpecies)//"]"
       enddo
    
       iOff = 4+nSpeciesTotal
       write(iOutputUnit_,"(I7,A1,a)")  iOff+1, " ", "Temperature"
       write(iOutputUnit_,"(I7,A1,a)")  iOff+2, " ", "V!Dn!N (east)"
       write(iOutputUnit_,"(I7,A1,a)")  iOff+3, " ", "V!Dn!N (north)"
       write(iOutputUnit_,"(I7,A1,a)")  iOff+4, " ", "V!Dn!N (up)"

       iOff = 8+nSpeciesTotal
       do iSpecies = 1, nSpecies
          write(iOutputUnit_,"(I7,A1,a)")  iOff+iSpecies, " ",&
               "V!Dn!N (up,"//cSpecies(iSpecies)//")"
       enddo

    endif

    if (cType(3:5) == "ALL" .or. cType(3:5) == "ION") then

       iOff = 3
       if (cType(3:5) == "ALL") iOff = 8+nSpeciesTotal+nSpecies
       do iIon = 1, nIons
          write(iOutputUnit_,"(I7,A1,a)") iOff+iIon, " ", "["//cIons(iIon)//"]"
       enddo

       iOff = iOff+nIons

       write(iOutputUnit_,"(I7,A1,a)") iOff+1, " ", "eTemperature"
       write(iOutputUnit_,"(I7,A1,a)") iOff+2, " ", "iTemperature"
       write(iOutputUnit_,"(I7,A1,a)") iOff+3, " ", "V!Di!N (east)"
       write(iOutputUnit_,"(I7,A1,a)") iOff+4, " ", "V!Di!N (north)"
       write(iOutputUnit_,"(I7,A1,a)") iOff+5, " ", "V!Di!N (up)"

       iOff = iOff + 5

       if (cType(3:5) == "ALL") then

          write(iOutputUnit_,"(I7,A1,a)") iOff+1, " ", "N2 Mixing Ratio"
          write(iOutputUnit_,"(I7,A1,a)") iOff+2, " ", "CH4 Mixing Ratio"
          write(iOutputUnit_,"(I7,A1,a)") iOff+3, " ", "Ar Mixing Ratio"
          write(iOutputUnit_,"(I7,A1,a)") iOff+4, " ", "HCN Mixing Ratio"
          write(iOutputUnit_,"(I7,A1,a)") iOff+5, " ", "H2 Mixing Ratio"

!       write(iOutputUnit_,"(I7,A1,a)") iOff+6, " ", "15N2 Mixing Ratio"
!       write(iOutputUnit_,"(I7,A1,a)") iOff+7, " ", "13CH4 Mixing Ratio"

          iOff = iOff + nSpecies
          write(iOutputUnit_,"(I7,A1,a)") iOff+1, " ", "RadCooling"
          write(iOutputUnit_,"(I7,A1,a)") iOff+2, " ", "EuvHeating"
          write(iOutputUnit_,"(I7,A1,a)") iOff+3, " ", "Conduction"
          write(iOutputUnit_,"(I7,A1,a)") iOff+4, " ", "Heat Balance Total"
          write(iOutputUnit_,"(I7,A1,a)") iOff+5, " ", "Heaing Efficiency"

       else

          write(iOutputUnit_,"(I7,A1,a)") iOff+1, " ", "B0 East"
          write(iOutputUnit_,"(I7,A1,a)") iOff+2, " ", "B0 North"
          write(iOutputUnit_,"(I7,A1,a)") iOff+3, " ", "B0 Vertical"
          write(iOutputUnit_,"(I7,A1,a)") iOff+4, " ", "B0 Magnitude"
          write(iOutputUnit_,"(I7,A1,a)") iOff+5, " ", "Magnetic Latitude"
          write(iOutputUnit_,"(I7,A1,a)") iOff+6, " ", "Magnetic Longitude"

       endif

    endif

    write(iOutputUnit_,*) ""

  end subroutine output_header


 subroutine output_header_new

    use ModElectrodynamics, only : nMagLats, nMagLons

    integer :: iOff, iSpecies, iIon

    write(iOutputUnit_,*) "NUMERICAL VALUES"

    write(iOutputUnit_,"(I7,6A)") nvars_to_write, " nvars"
    write(iOutputUnit_,"(I7,7A)") nAlts, " nAltitudes"
    write(iOutputUnit_,"(I7,7A)") 1, " nLatitudes"
    write(iOutputUnit_,"(I7,7A)") 1, " nLongitudes"
    write(iOutputUnit_,*) ""

    write(iOutputUnit_,*) "VARIABLE LIST"
    write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
    write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Local Time"
    write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Latitude"
    write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Altitude"
    write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Solar Zenith Angle"
    write(iOutputUnit_,"(I7,A1,a)")  6, " ", "Rho"
    iOff = 6
    do iSpecies = 1, nSpeciesTotal
       write(iOutputUnit_,"(I7,A1,a)")  iOff+iSpecies, " ", &
            "["//cSpecies(iSpecies)//"]"
    enddo
    
    iOff = iOff + nSpeciesTotal
    write(iOutputUnit_,"(I7,A1,a)")  iOff+1, " ", "Temperature"
    write(iOutputUnit_,"(I7,A1,a)")  iOff+2, " ", "V!Dn!N (east)"
    write(iOutputUnit_,"(I7,A1,a)")  iOff+3, " ", "V!Dn!N (north)"
    write(iOutputUnit_,"(I7,A1,a)")  iOff+4, " ", "V!Dn!N (up)"

    iOff = iOff + 4
    do iSpecies = 1, nSpecies
       write(iOutputUnit_,"(I7,A1,a)")  iOff+iSpecies, " ",&
            "V!Dn!N (up,"//cSpecies(iSpecies)//")"
    enddo

    iOff = iOff + nSpecies
    do iIon = 1, nIons
       write(iOutputUnit_,"(I7,A1,a)") iOff+iIon, " ", "["//cIons(iIon)//"]"
    enddo

    iOff = iOff+nIons
    write(iOutputUnit_,"(I7,A1,a)") iOff+1, " ", "eTemperature"
    write(iOutputUnit_,"(I7,A1,a)") iOff+2, " ", "iTemperature"
    write(iOutputUnit_,"(I7,A1,a)") iOff+3, " ", "V!Di!N (east)"
    write(iOutputUnit_,"(I7,A1,a)") iOff+4, " ", "V!Di!N (north)"
    write(iOutputUnit_,"(I7,A1,a)") iOff+5, " ", "V!Di!N (up)"

    iOff = iOff + 5
    do iSpecies = 1, nSpecies
       write(iOutputUnit_,"(I7,A1,a)")  iOff+iSpecies, " ", &
            " "//cSpecies(iSpecies)//"Mixing Ratio"
    enddo


    write(iOutputUnit_,*) ""

  end subroutine output_header_new

  !----------------------------------------------------------------

  subroutine write_head_blocks

    if (cType(1:2) == "1D") return

    write(iOutputUnit_,*) "BLOCKS"
    write(iOutputUnit_,"(I7,A)") 1, " nBlocksAlt"
    if (cType /= "2DMEL") then
       write(iOutputUnit_,"(I7,A)") nBlocksLat, " nBlocksLat"
       write(iOutputUnit_,"(I7,A)") nBlocksLon, " nBlocksLon"
    else
       write(iOutputUnit_,"(I7,A)") 1, " nBlocksLat"
       write(iOutputUnit_,"(I7,A)") 1, " nBlocksLon"
    endif
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
    write(iOutputUnit_,"(I7,A)") iTimeArray(7), " Millisecond"
    write(iOutputUnit_,*) ""

  end subroutine write_head_time

  !----------------------------------------------------------------
  !
  !----------------------------------------------------------------

  subroutine write_head_version

    write(iOutputUnit_,*) "VERSION"
    write(iOutputUnit_,*) 2.4+PlanetNum
    write(iOutputUnit_,*) ""

  end subroutine write_head_version

end subroutine output

!!  !----------------------------------------------------------------
!!  !
!!  !----------------------------------------------------------------
!!  
!!  subroutine output_1d(dir, cName, iBlock, Position)
!!  
!!    use ModGITM
!!    use ModEUV
!!    use ModTime
!!    use ModInputs
!!    use ModSources
!!    use ModConstants
!!  
!!    implicit none
!!  
!!    character (len=*), intent(in) :: dir
!!    character (len=*), intent(in) :: cName
!!    integer, intent(in)           :: iBlock
!!    real, intent(in)              :: Position(3)
!!  
!!    character (len=14) :: cTime
!!    integer :: iLon,iLat,iAlt, nvars_to_write, nlines, iBLK, i
!!    integer :: iiLon, iiLat, iiAlt
!!    logical :: done
!!  
!!    real :: LatFind, LonFind
!!    real :: rLon, rLat
!!  
!!    character (len=2) :: cYear, cMonth, cDay, cHour, cMinute, cSecond
!!  
!!    LatFind = Position(iNorth_)
!!    LonFind = Position(iEast_)
!!  
!!    if ((Longitude(0,iBlock)+Longitude(1,iBlock))/2 <=LonFind .and. &
!!         (Longitude(nLons,iBlock)+Longitude(nLons+1,iBlock))/2 >LonFind) then
!!       if ((Latitude(0,iBlock)+Latitude(1,iBlock))/2 <=LatFind .and. &
!!            (Latitude(nLats,iBlock)+Latitude(nLats+1,iBlock))/2 >LatFind) then
!!  
!!          iiLat = -1
!!          iiLon = -1
!!          do iLon = 0,nLons
!!             if (Longitude(iLon,iBlock) <= LonFind .and. &
!!                  Longitude(iLon+1,iBlock) > LonFind) then
!!                iiLon = iLon
!!                rLon = 1.0 - (LonFind - Longitude(iLon,iBlock)) / &
!!                     (Longitude(iLon+1,iBlock) - Longitude(iLon,iBlock))
!!             endif
!!          enddo
!!  
!!          do iLat = 0,nLats
!!             if (Latitude(iLat,iBlock) <= LatFind .and. &
!!                  Latitude(iLat+1,iBlock) > LatFind) then
!!                iiLat = iLat
!!                rLat = 1.0 - (LatFind - Latitude(iLat,iBlock)) / &
!!                     (Latitude(iLat+1,iBlock) - Latitude(iLat,iBlock))
!!             endif
!!          enddo
!!  
!!       else
!!          return
!!       endif
!!    else 
!!       return
!!    endif
!!  
!!    if (iProc == 0 .and. iBlock == 1) &
!!         write(*,'(a,i7,i5,5i3)') &
!!         "Writing Output files at iStep : ",iStep, iTimeArray(1:6)
!!  
!!    call calc_physics(iBlock)
!!    call chapman_integrals(iBlock)
!!    call calc_rates(iBlock)
!!    if (.not. Is1D) call calc_efield(iBlock)
!!  
!!    !! construct naming strings
!!  
!!    call i2s(mod(iTimeArray(1),100), cYear, 2)
!!    call i2s(iTimeArray(2), cMonth, 2)
!!    call i2s(iTimeArray(3), cDay, 2)
!!    call i2s(iTimeArray(4), cHour, 2)
!!    call i2s(iTimeArray(5), cMinute, 2)
!!    call i2s(iTimeArray(6), cSecond, 2)
!!  
!!    cTime = "t"//cYear//cMonth//cDay//"_"//cHour//cMinute//cSecond
!!  
!!    !! open file
!!    open(unit=iOutputUnit_, &
!!         file=dir//"/"//cName//"_"//cTime//"."//"dat",&
!!         status="unknown")
!!  
!!    write(iOutputUnit_,*) ""
!!    write(iOutputUnit_,*) "TIME"
!!    do i=1,7
!!       write(iOutputUnit_,*) iTimeArray(i)
!!    enddo
!!    write(iOutputUnit_,*) ""
!!  
!!    call output_header(.true., nVars_to_Write)
!!  
!!    call output_1dall(iiLon, iiLat, iBlock, rLon, rLat, iOutputUnit_)
!!  
!!    close(unit=iOutputUnit_)
!!  
!!  end subroutine output_1d

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_3dall(iBlock)

  use ModGITM
  use ModInputs

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon, iiAlt, iiLat, iiLon, i

  do iAlt=-1,nAlts+2
     !!! Why ???
     iiAlt = max(min(iAlt,nAlts),1)
     do iLat=-1,nLats+2
        !!! Why ???
        iiLat = min(max(iLat,1),nLats)
        do iLon=-1,nLons+2
           !!! Why ???
           iiLon = min(max(iLon,1),nLons)
           write(iOutputUnit_)       &
                Longitude(iLon,iBlock), &
                Latitude(iLat,iBlock), &
                Altitude_GB(iLon,iLat,iAlt,iBlock),&
                Rho(iLon,iLat,iAlt,iBlock),&
                (NDensityS(iLon,iLat,iAlt,i,iBlock),i=1,nSpeciesTotal), &
                Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt),&
                (Velocity(iLon,iLat,iAlt,i,iBlock),i=1,3), &
                (VerticalVelocity(iLon,iLat,iAlt,i,iBlock),i=1,nSpecies), &
                (IDensityS(iLon,iLat,iAlt,i,iBlock),i=1,nIons), &
                eTemperature(iLon,iLat,iAlt,iBlock)  ,&
                ITemperature(iLon,iLat,iAlt,iBlock)  ,&
                (Ivelocity(iLon,iLat,iAlt,i,iBlock),i=1,3)
        enddo
     enddo
  enddo

end subroutine output_3dall

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_3dneu(iBlock)

  use ModGITM
  use ModInputs

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon, iiAlt, iiLat, iiLon

  do iAlt=-1,nAlts+2
     iiAlt = max(min(iAlt,nAlts),1)
     do iLat=-1,nLats+2
        iiLat = min(max(iLat,1),nLats)
        do iLon=-1,nLons+2
           iiLon = min(max(iLon,1),nLons)
           write(iOutputUnit_)       &
                Longitude(iLon,iBlock), &
                Latitude(iLat,iBlock), &
                Altitude_GB(iLon,iLat,iAlt,iBlock),&
                Rho(iLon,iLat,iAlt,iBlock),&
                NDensityS(iLon,iLat,iAlt,:,iBlock), &
                Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt),&
                velocity(iLon,iLat,iAlt,:,iBlock) , &
                VerticalVelocity(iLon,iLat,iAlt,:,iBlock)
        enddo
     enddo
  enddo

end subroutine output_3dneu

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_3dion(iBlock)

  use ModGITM
  use ModInputs

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon, iiAlt, iiLat, iiLon

  do iAlt=-1,nAlts+2
     do iLat=-1,nLats+2
        do iLon=-1,nLons+2
           write(iOutputUnit_)          &
                Longitude(iLon,iBlock),               &
                Latitude(iLat,iBlock),                &
                Altitude_GB(iLon,iLat,iAlt,iBlock),   &
                IDensityS(iLon,iLat,iAlt,:,iBlock),   &
                eTemperature(iLon,iLat,iAlt,iBlock),  &
                ITemperature(iLon,iLat,iAlt,iBlock),  &
                Ivelocity(iLon,iLat,iAlt,:,iBlock),   &
                B0(iLon,iLat,iAlt,:,iBlock), &
                mLatitude(iLon,iLat,iAlt,iBlock), &
                mLongitude(iLon,iLat,iAlt,iBlock)
        enddo
     enddo
  enddo

end subroutine output_3dion

!----------------------------------------------------------------
!
!----------------------------------------------------------------
subroutine output_3dthm(iBlock)

  use ModGITM
  use ModInputs
  use ModSources
  use ModEuv, only : EuvTotal
  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon, iiAlt, iiLat, iiLon


  do iAlt=-1,nAlts+2
     iiAlt = max(min(iAlt,nAlts),1)
     do iLat=-1,nLats+2
        iiLat = min(max(iLat,1),nLats)
        do iLon=-1,nLons+2
           iiLon = min(max(iLon,1),nLons)

           write(iOutputUnit_)          &
                Longitude(iLon,iBlock),               &
                Latitude(iLat,iBlock),                &
                Altitude_GB(iLon,iLat,iAlt,iBlock),   &
                EuvHeating(iiLon,iiLat,iiAlt,iBlock)*dt*TempUnit(iiLon,iiLat,iiAlt),    &
                Conduction(iiLon,iiLat,iiAlt)*TempUnit(iiLon,iiLat,iiAlt),              &
                MoleConduction(iiLon,iiLat,iiAlt),                             &
                EddyCond(iiLon,iiLat,iiAlt),                                   &
                EddyCondAdia(iiLon,iiLat,iiAlt),                               &
                ChemicalHeatingRate(iiLon,iiLat,iiAlt)*TempUnit(iiLon,iiLat,iiAlt),     &
                AuroralHeating(iiLon,iiLat,iiAlt)*dt*TempUnit(iiLon,iiLat,iiAlt),       &
                JouleHeating(iiLon,iiLat,iiAlt)*dt*TempUnit(iiLon,iiLat,iiAlt),         &
                -NOCooling(iiLon,iiLat,iiAlt)*dt*TempUnit(iiLon,iiLat,iiAlt),           &
                -OCooling(iiLon,iiLat,iiAlt)*dt*TempUnit(iiLon,iiLat,iiAlt),            &
                EuvTotal(iiLon,iiLat,iiAlt,iBlock) * dt
           
        enddo
     enddo
  enddo
     
end subroutine output_3dthm

!----------------------------------------------------------------
!
!----------------------------------------------------------------


subroutine output_1dthm

  use ModGITM
  use ModInputs
  use ModSources
  use ModEuv, only : EuvTotal
  implicit none

  integer :: iAlt, iLat, iLon, iiAlt,iSpecies
  real    :: varsS(nSpeciesTotal),varsL(nSpeciesTotal)

  
  do iAlt=-1,nAlts+2
     iiAlt = max(min(iAlt,nAlts),1)
 
     do iSpecies = 1, nSpeciesTotal 
        varsS(iSpecies) = NeutralSourcesTotal(iSpecies,iiAlt)
        varsL(iSpecies) = NeutralLossesTotal(iSpecies,iiAlt)
     enddo

     write(iOutputUnit_) &
          Longitude(1,1),               &
          Latitude(1,1),                &
          Altitude_GB(1,1,iAlt,1),   &
          EuvHeating(1,1,iiAlt,1)*dt*TempUnit(1,1,iiAlt),    &
          Conduction(1,1,iiAlt)*TempUnit(1,1,iiAlt),              &
          MoleConduction(1,1,iiAlt),                             &
          EddyCond(1,1,iiAlt),                                   &
          EddyCondAdia(1,1,iiAlt),                               &
          ChemicalHeatingRate(1,1,iiAlt)*TempUnit(1,1,iiAlt),     &
          AuroralHeating(1,1,iiAlt)*dt*TempUnit(1,1,iiAlt),       &
          JouleHeating(1,1,iiAlt)*dt*TempUnit(1,1,iiAlt),         &
          -RadCooling(1,1,iiAlt,1)*dt*TempUnit(1,1,iiAlt),           &
          -OCooling(1,1,iiAlt)*dt*TempUnit(1,1,iiAlt),            &
          EuvTotal(1,1,iiAlt,1) * dt,                             &
          varsS, varsL
                   
  enddo

end subroutine output_1dthm

!----------------------------------------------------------------
!
!----------------------------------------------------------------
subroutine output_3dchm(iBlock)

  use ModGITM
  use ModInputs
  use ModSources
  use ModConstants
  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon, iiAlt, iiLat, iiLon, iReact
  real :: vars(nReactions)

  do iAlt=-1,nAlts+2
     iiAlt = max(min(iAlt,nAlts),1)
     do iLat=-1,nLats+2
        iiLat = min(max(iLat,1),nLats)
        do iLon=-1,nLons+2
           iiLon = min(max(iLon,1),nLons)
           do iReact = 1, nReactions
              
              vars(iReact) = ChemicalHeatingSpecies(iiLon,iiLat,iiAlt,iReact) / &
                   Element_Charge
              
              enddo
              
              write(iOutputUnit_) &
                   Longitude(iLon,iBlock),               &
                   Latitude(iLat,iBlock),                &
                   Altitude_GB(iLon,iLat,iAlt,iBlock),   &
                   Vars, &
                   ChemicalHeatingRate(iiLon,iiLat,iiAlt) * &
                   cp(iilon,iiLat,iiAlt,iBlock) *   &
                   Rho(iilon,iiLat,iiAlt,iBlock)*TempUnit(iilon,iiLat,iiAlt) / &
                   Element_Charge
        enddo
     enddo
  enddo

end subroutine output_3dchm

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_3dglo(iBlock)

  use ModGITM
  use ModInputs

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon, iiAlt, iiLat, iiLon

  do iAlt=-1,nAlts+2
     iiAlt = max(min(iAlt,nAlts),1)
     do iLat=-1,nLats+2
        iiLat = min(max(iLat,1),nLats)
        do iLon=-1,nLons+2
           iiLon = min(max(iLon,1),nLons)
              
              write(iOutputUnit_) &
                   Longitude(iLon,iBlock),               &
                   Latitude(iLat,iBlock),                &
                   Altitude_GB(iLon,iLat,iAlt,iBlock),   &
                   vEmissionRate(iiLon,iiLat,iiAlt,i6300_,iBlock), &
                   PhotoEFluxTotal(iiLon,iiLat,iiAlt,iBlock,1),    &
                   PhotoEFluxTotal(iiLon,iiLat,iiAlt,iBlock,2)
                   
        enddo
     enddo
  enddo

end subroutine output_3dglo

!----------------------------------------------------------------
!
!----------------------------------------------------------------
subroutine output_1dglo

  use ModGITM
  use ModInputs

  implicit none

  integer :: iAlt, iLat, iLon, iiAlt

  do iAlt=-1,nAlts+2
     iiAlt = max(min(iAlt,nAlts),1)
 
     write(iOutputUnit_) &
          Longitude(1,1),               &
          Latitude(1,1),                &
          Altitude_GB(1,1,iAlt,1),   &
          vEmissionRate(1,1,iiAlt,i6300_,1), &
          PhotoEFluxTotal(1,1,iiAlt,1,1),    &
          PhotoEFluxTotal(1,1,iiAlt,1,2)
                   
  enddo

end subroutine output_1dglo

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_2dgel(iBlock)

  use ModElectrodynamics
  use ModGITM
  use ModInputs

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon, iiAlt, iiLat, iiLon

  iAlt = 1
  do iLat=1,nLats
     do iLon=1,nLons
        write(iOutputUnit_)       &
             Longitude(iLon,iBlock), &
             Latitude(iLat,iBlock),&
             Altitude_GB(iLon,iLat,iAlt,iBlock), &
             Potential(iLon,iLat,iAlt,iBlock), &
             PedersenConductance(iLon,iLat,iBlock), &
             HallConductance(iLon,iLat,iBlock), &
             ElectronAverageEnergy(iLon,iLat), &
             ElectronEnergyFlux(iLon,iLat), &
             DivJuAlt(iLon,iLat), &
             PedersenFieldLine(iLon, iLat), &
             HallFieldLine(iLon, iLat), &
             DivJuFieldLine(iLon, iLat), &
             LengthFieldLine(iLon, iLat)
     enddo
  enddo

end subroutine output_2dgel

!----------------------------------------------------------------

subroutine output_2dmel(iBlock)
  
  use ModElectrodynamics
  use ModConstants, only:Pi
  use ModGITM
  use ModInputs

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon
  !--------------------------------------------------------------------------

  iAlt = 1
  do iLat=1,nMagLats
     do iLon=1,nMagLons+1
        write(iOutputUnit_)       &
             MagLonMC(iLon,iLat)*Pi/180.0, &
             MagLatMC(iLon,iLat)*Pi/180.0, &
             Altitude_GB(iLon,iLat,iAlt,iBlock), &
             MagLocTimeMC(iLon,iLat)*Pi/180.0, &
             GeoLatMC(iLon,iLat), &
             GeoLonMC(iLon,iLat), &
             SigmaPedersenMC(iLon,iLat), &
             SigmaHallMC(iLon,iLat), &
             DivJuAltMC(iLon,iLat), &
             LengthMC(iLon,iLat), &
             SigmaPPMC(iLon,iLat), &
             SigmaLLMC(iLon,iLat), &
             SigmaHHMC(iLon,iLat), &
             SigmaCCMC(iLon,iLat), &
             SigmaPLMC(iLon,iLat), &
             SigmaLPMC(iLon,iLat), &
             KDpmMC(iLon,iLat), &
             KdlmMC(iLon,iLat), &
             solver_a_mc(iLon,iLat), &
             solver_b_mc(iLon,iLat), &
             solver_c_mc(iLon,iLat), &
             solver_d_mc(iLon,iLat), &
             solver_e_mc(iLon,iLat), &
             solver_s_mc(iLon,iLat), &
             DynamoPotentialMC(iLon,iLat)
     enddo
  enddo

end subroutine output_2dmel



!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_1dall(iiLon, iiLat, iBlock, rLon, rLat, iUnit)

  use ModGITM
  use ModEUV, only: HeatingEfficiency_CB
  use ModSources, only: JouleHeating, RadCooling, EuvHeating, Conduction
                        
  use ModInputs, only: iOutputUnit_
  implicit none

  integer, intent(in) :: iiLat, iiLon, iBlock, iUnit
  real, intent(in)    :: rLon, rLat

  integer, parameter :: nVars = 13+nSpeciesTotal+nSpecies+nIons+nSpecies+5
  real :: Vars(nVars)
  real :: Tmp(0:nLons+1,0:nLats+1)
  integer :: iAlt, iiAlt, iOff, iIon, iSpecies, iDir

  do iAlt=-1,nAlts+2

     iiAlt = max(min(iAlt,nAlts),1)


     Vars(1) = &
          rLon*Longitude(iiLon,iBlock)+(1-rLon)*Longitude(iiLon+1,iBlock)
     Vars(2) = &
          rLat*Latitude(iiLat,iBlock)+(1-rLat)*Latitude(iiLat+1,iBlock)

     Vars(3) = Altitude_GB(iiLon, iiLat, iAlt, iBlock)

     Tmp = Rho(0:nLons+1,0:nLats+1,iAlt,iBlock)
     Vars(4) = inter(Tmp,iiLon,iiLat,rlon,rlat)

     iOff = 4
     do iSpecies = 1, nSpeciesTotal
        Tmp = NDensityS(0:nLons+1,0:nLats+1,iAlt,iSpecies,iBlock)
!        Tmp = NDensityS(0:nLons+1,0:nLats+1,iAlt,iSpecies,iBlock)/NDensity(0:nLons+1,0:nLats+1,iAlt,iBlock)
        Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     Tmp = Temperature(0:nLons+1,0:nLats+1,iAlt,iBlock) * &
          TempUnit(0:nLons+1,0:nLats+1,iAlt)
     iOff = 5+nSpeciesTotal
     Vars(iOff) = inter(Tmp,iiLon,iiLat,rlon,rlat)

     do iDir = 1, 3
        Tmp = Velocity(0:nLons+1,0:nLats+1,iAlt,iDir,iBlock)
        Vars(iOff+iDir) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     iOff = 8+nSpeciesTotal
     do iSpecies = 1, nSpecies
        Tmp = VerticalVelocity(0:nLons+1,0:nLats+1,iAlt,iSpecies,iBlock)
        Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     iOff = 8+nSpeciesTotal+nSpecies
     do iIon = 1, nIons
        Tmp = IDensityS(0:nLons+1,0:nLats+1,iAlt,iIon,iBlock)
        Vars(iOff+iIon) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     iOff = 8+nSpeciesTotal+nSpecies+nIons+1
     Tmp = eTemperature(0:nLons+1,0:nLats+1,iAlt,iBlock)
     Vars(iOff)   = inter(Tmp,iiLon,iiLat,rlon,rlat)

     iOff = iOff+1
     Tmp = iTemperature(0:nLons+1,0:nLats+1,iAlt,iBlock)
     Vars(iOff) = inter(Tmp,iiLon,iiLat,rlon,rlat)

!     do iDir = 1, 3
!        Tmp = IVelocity(0:nLons+1,0:nLats+1,iAlt,iDir,iBlock)
!        Vars(iOff+iDir) = inter(Tmp,iiLon,iiLat,rlon,rlat)
!     enddo

        iOff = iOff + 1
        Tmp = IVelocity(0:nLons+1,0:nLats+1,iAlt,iEast_,iBlock)
        Vars(iOff) = inter(Tmp,iiLon,iiLat,rlon,rlat)

        iOff = iOff + 1
        Tmp = IVelocity(0:nLons+1,0:nLats+1,iAlt,iNorth_,iBlock)
        Vars(iOff) = inter(Tmp,iiLon,iiLat,rlon,rlat)

        iOff = iOff + 1
        Tmp = IVelocity(0:nLons+1,0:nLats+1,iAlt,iUp_,iBlock)
        Vars(iOff) = inter(Tmp,iiLon,iiLat,rlon,rlat)

        do iSpecies = 1, nSpecies
           Tmp = NDensityS(0:nLons+1,0:nLats+1,iAlt,iSpecies,iBlock)/NDensity(0:nLons+1,0:nLats+1,iAlt,iBlock)
           Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,rlon,rlat)
        enddo

!        iOff = iOff + 1
!        Tmp = NDensityS(0:nLons+1,0:nLats+1,iAlt,iSpecies,iBlock)/NDensity(0:nLons+1,0:nLats+1,iAlt,iBlock)
!        Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,rlon,rlat)

        iOff = iOff + nSpecies
        Vars(iOff+1) = Dt*RadCooling(1,1,iiAlt,iBlock)*TempUnit(1,1,iiAlt)

        Vars(iOff+2) = Dt*EuvHeating(1,1,iiAlt,iBlock)*TempUnit(1,1,iiAlt)

        Vars(iOff+3) = Conduction(1,1,iiAlt)*TempUnit(1,1,iiAlt)

        Vars(iOff+4) = Dt*EuvHeating(1,1,iiAlt,iBlock)*TempUnit(1,1,iiAlt) - &
         Dt*RadCooling(1,1,iiAlt,iBlock)*TempUnit(1,1,iiAlt) + &
         Conduction(1,1,iiAlt)*TempUnit(1,1,iiAlt)

        Vars(iOff+5) = HeatingEfficiency_CB(1,1,iiAlt,iBlock)

     write(iOutputUnit_) Vars

  enddo

contains

  real function inter(variable, iiLon, iiLat, rLon, rLat) result(PointValue)

    implicit none

    real :: variable(0:nLons+1, 0:nLats+1), rLon, rLat
    integer :: iiLon, iiLat

    PointValue = &  
          (  rLon)*(  rLat)*Variable(iiLon  ,iiLat  ) + &
          (1-rLon)*(  rLat)*Variable(iiLon+1,iiLat  ) + &
          (  rLon)*(1-rLat)*Variable(iiLon  ,iiLat+1) + &
          (1-rLon)*(1-rLat)*Variable(iiLon+1,iiLat+1)

  end function inter

end subroutine output_1dall

subroutine output_1dnew(iiLon, iiLat, iBlock, rLon, rLat, iUnit)

  use ModGITM
  use ModEUV, only : Sza, HeatingEfficiency_CB
  use ModSources, only: JouleHeating
  use ModInputs, only: iOutputUnit_
  implicit none

  integer, intent(in) :: iiLat, iiLon, iBlock, iUnit
  real, intent(in)    :: rLon, rLat

  integer, parameter :: nVars = 13+nSpeciesTotal+nSpecies+nIons+nSpecies+2
  real :: Vars(nVars)
  real :: Tmp(0:nLons+1,0:nLats+1)
  integer :: iAlt, iiAlt, iOff, iIon, iSpecies, iDir

  do iAlt=-1,nAlts+2

     iiAlt = max(min(iAlt,nAlts),1)

     Vars(1) = &
          rLon*Longitude(iiLon,iBlock)+(1-rLon)*Longitude(iiLon+1,iBlock)
     Vars(2) = &
          rLon*LocalTime(iiLon)+(1-rLon)*LocalTime(iiLon+1)
     Vars(3) = &
          rLat*Latitude(iiLat,iBlock)+(1-rLat)*Latitude(iiLat+1,iBlock)

     Vars(4) = Sza(iiLon,iiLat,iBlock)

     Vars(5) = Altitude_GB(iiLon, iiLat, iAlt, iBlock)

     Tmp = Rho(0:nLons+1,0:nLats+1,iAlt,iBlock)
     Vars(6) = inter(Tmp,iiLon,iiLat,rlon,rlat)

     iOff = 6
     do iSpecies = 1, nSpeciesTotal
        Tmp = NDensityS(0:nLons+1,0:nLats+1,iAlt,iSpecies,iBlock)
        Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     Tmp = Temperature(0:nLons+1,0:nLats+1,iAlt,iBlock) * &
          TempUnit(0:nLons+1,0:nLats+1,iAlt)
     iOff = 7+nSpeciesTotal
     Vars(iOff) = inter(Tmp,iiLon,iiLat,rlon,rlat)

     do iDir = 1, 3
        Tmp = Velocity(0:nLons+1,0:nLats+1,iAlt,iDir,iBlock)
        Vars(iOff+iDir) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     iOff = 10+nSpeciesTotal
     do iSpecies = 1, nSpecies
        Tmp = VerticalVelocity(0:nLons+1,0:nLats+1,iAlt,iSpecies,iBlock)
        Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     iOff = 10+nSpeciesTotal+nSpecies
     do iIon = 1, nIons
        Tmp = IDensityS(0:nLons+1,0:nLats+1,iAlt,iIon,iBlock)
        Vars(iOff+iIon) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     iOff = iOff+1
     Tmp = eTemperature(0:nLons+1,0:nLats+1,iAlt,iBlock)
     Vars(iOff)   = inter(Tmp,iiLon,iiLat,rlon,rlat)

     iOff = iOff+1
     Tmp = iTemperature(0:nLons+1,0:nLats+1,iAlt,iBlock)
     Vars(iOff) = inter(Tmp,iiLon,iiLat,rlon,rlat)

     iOff = iOff
     do iDir = 1, 3
        Tmp = IVelocity(0:nLons+1,0:nLats+1,iAlt,iDir,iBlock)
        Vars(iOff+iDir) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     iOff = iOff
     do iSpecies = 1, nSpecies
        Tmp = NDensityS(0:nLons+1,0:nLats+1,iAlt,iSpecies,iBlock)/NDensity(0:nLons+1,0:nLats+1,iAlt,iBlock)
        Vars(iOff+iSpecies) = inter(Tmp,iiLon,iiLat,rlon,rlat)
     enddo

     write(iOutputUnit_) Vars

  enddo

contains

  real function inter(variable, iiLon, iiLat, rLon, rLat) result(PointValue)

    implicit none

    real :: variable(0:nLons+1, 0:nLats+1), rLon, rLat
    integer :: iiLon, iiLat

    PointValue = &  
          (  rLon)*(  rLat)*Variable(iiLon  ,iiLat  ) + &
          (1-rLon)*(  rLat)*Variable(iiLon+1,iiLat  ) + &
          (  rLon)*(1-rLat)*Variable(iiLon  ,iiLat+1) + &
          (1-rLon)*(1-rLat)*Variable(iiLon+1,iiLat+1)

  end function inter

end subroutine output_1dnew



