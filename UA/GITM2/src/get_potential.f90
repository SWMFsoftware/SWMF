subroutine get_potential(iBlock)

  use ModGITM
  use ModTime
  use ModIndicesInterfaces
  use ModInputs
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock

  integer :: iError, iLat, iLon, iAlt
  logical :: IsFirstTime = .true.
  logical :: IsFirstPotential(nBlocksMax) = .true.
  logical :: IsFirstAurora(nBlocksMax) = .true.

  logical :: UseIMF = .false.
  logical :: UseHPI = .false.

  real, dimension(-1:nLons+2,-1:nLats+2) :: TempPotential, Grid

  character (len=100), dimension(100) :: Lines
  character (len=100) :: TimeLine

  real :: temp

  call start_timing("Get AMIE Potential")
  call report("get_potential",2)

  iError = 0

  if (index(cPlanet,"Earth") == 0) then 

     potential = 0.0
     ElectronAverageEnergy = 0.1
     ElectronEnergyFlux = 0.0001
     return

  endif

  if (IsFirstTime) then

     if (index(cAMIEFileNorth,"none") > 0) then

        Lines(1) = "#BACKGROUND"
        Lines(2) = "UA/DataIn/"

        UseHPI = .true.
        call get_IMF_Bz(CurrentTime, temp, iError)
        call IO_SetIMFBz(temp)
        if (iError /= 0) then
           write(*,*) "Can not find IMF Bz."
           write(*,*) "Setting potential to Millstone HPI."
           Lines(3) = "millstone_hpi"    ! Change to "zero" if you want
        else
           write(*,*) "Setting potential to Weimer [1996]."
           Lines(3) = "weimer96"    ! Change to "zero" if you want
           UseIMF = .true.
        endif
        Lines(4) = "ihp"
        Lines(5) = "idontknow"
        Lines(6) = ""

     else

        if (index(cAMIEFileNorth,"mhd") > 0) then

           Lines(1) = "#MHDFILE"
           Lines(2) = cAMIEFileSouth
           write(TimeLine,'(i4)') iTimeArray(1)
           Lines(3) = TimeLine
           write(TimeLine,'(i2)') iTimeArray(2)
           Lines(4) = TimeLine
           write(TimeLine,'(i2)') iTimeArray(3)
           Lines(5) = TimeLine
           Lines(6) = ""

           UseIMF = .false.

        else

           Lines(1) = "#AMIEFILES"
           Lines(2) = cAMIEFileNorth
           Lines(3) = cAMIEFileSouth
           Lines(4) = ""
           Lines(5) = ""
           Lines(6) = ""

           UseIMF = .false.

        endif

     endif

     Lines(7) = "#DEBUG"
     Lines(8) = "0"
     Lines(9) = "0"
     Lines(10) = ""
     Lines(11) = "#END"

     call IE_set_inputs(Lines)

     call IE_Initialize(iError)

     if (iError /= 0) then
        write(*,*) &
             "Code Error in IE_Initialize called from get_potential.f90"
        write(*,*) "Error : ",iError
        call stop_gitm("Stopping in get_potential")
     endif

     call IO_SetnMLTs(nLons+4)
     call IO_SetnLats(nLats+4)
      
     IsFirstTime = .false.

  endif

  call IO_SetTime(CurrentTime)

  call IO_SetNorth

  if (UseIMF) then

     call get_IMF_Bz(CurrentTime, temp, iError)
     call IO_SetIMFBz(temp)

     if (iError /= 0) then
        write(*,*) &
             "Code Error in get_IMF_Bz called from get_potential.f90"
        write(*,*) "Code : ",iError
        call stop_gitm("Stopping in get_potential")
     endif

     if (iDebugLevel > 1) write(*,*) "==> Bz : ",temp

     call get_IMF_By(CurrentTime, temp, iError)
     call IO_SetIMFBy(temp)

     if (iError /= 0) then
        write(*,*) &
             "Code Error in get_IMF_By called from get_potential.f90"
        call stop_gitm("Stopping in get_potential")
     endif

     if (iDebugLevel > 1) write(*,*) "==> IMF By : ",temp

     call get_SW_V(CurrentTime, temp, iError)
     call IO_SetSWV(temp)

     if (iError /= 0) then
        write(*,*) "Code Error in get_sw_v called from get_potential.f90"
        call stop_gitm("Stopping in get_potential")
     endif

     if (iDebugLevel > 1) write(*,*) "==> Solar Wind Velocity : ",temp

     call get_kp(CurrentTime, temp, iError)
     call IO_Setkp(temp)

     if (iError /= 0) then
        write(*,*) "Code Error in get_kp called from get_potential.f90"
        call stop_gitm("Stopping in get_potential")
     endif

  endif

  if (UseHPI) then

     call get_HPI(CurrentTime, temp, iError)
     call IO_SetHPI(temp)

     if (iError /= 0) then
        write(*,*) "Code Error in get_hpi called from get_potential.f90"
        call stop_gitm("Stopping in get_potential")
     endif

  endif

  if (index(cAMIEFileNorth,"none") <= 0 .and. &
     index(cAMIEFileNorth,"SPS") <= 0 .and. iBlock == 1) then 
     if (iDebugLevel > 1) &
          write(*,*) "==> Reading AMIE values for time :",CurrentTime
     call get_AMIE_values(CurrentTime)
  endif

  if (iDebugLevel > 1) write(*,*) "==> Setting up IE Grid"

  if (floor((tSimulation-dt)/DtPotential) /= &
       floor((tsimulation)/DtPotential) .or. IsFirstPotential(iBlock)) then

     if (iDebugLevel > 1) write(*,*) "==> Getting Potential"

     Potential(:,:,:,iBlock) = 0.0

     do iAlt=-1,nAlts+2

        call start_timing("setgrid")
        call IO_SetGrid(                    &
             MLT(-1:nLons+2,-1:nLats+2,iAlt), &
             MLatitude(-1:nLons+2,-1:nLats+2,iAlt,iBlock), iError)
        call end_timing("setgrid")

        if (iError /= 0) then
           write(*,*) "Error in routine get_potential (IO_SetGrid):"
           write(*,*) iError
           call stop_gitm("Stopping in get_potential")
        endif

        if (iDebugLevel > 1 .and. iAlt == 1) &
             write(*,*) "==> Getting IE potential"

        call start_timing("getpotential")
        call IO_GetPotential(TempPotential, iError)
        call end_timing("getpotential")

        if (iError /= 0) then
           write(*,*) "Error in get_potential (IO_GetPotential):"
           write(*,*) iError
           call stop_gitm("Stopping in get_potential")
        endif

        Potential(:,:,iAlt,iBlock) = TempPotential

        !----------------------------------------------
        ! Another example of user output

        if (iAlt == 1) then 
           UserData2d(1:nLons,1:nLats,1,1,iBlock) = &
                TempPotential(1:nLons,1:nLats)/1000.0
        endif

     enddo

     IsFirstPotential(iBlock) = .false.

  endif

  if (floor((tSimulation-dt)/DtAurora) /= &
       floor((tsimulation)/DtAurora) .or. IsFirstAurora(iBlock)) then

     if (iDebugLevel > 1) write(*,*) "==> Getting Aurora"

     iAlt = nAlts+1

     call start_timing("setgrid")
     call IO_SetGrid(                    &
          MLT(-1:nLons+2,-1:nLats+2,iAlt), &
          MLatitude(-1:nLons+2,-1:nLats+2,iAlt,iBlock), iError)
     call end_timing("setgrid")

     if (iError /= 0) then
        write(*,*) "Error in routine get_potential (IO_SetGrid):"
        write(*,*) iError
        call stop_gitm("Stopping in get_potential")
     endif

     call IO_GetAveE(ElectronAverageEnergy, iError)
     if (iError /= 0) then
        write(*,*) "Error in get_potential (IO_GetAveE):"
        write(*,*) iError
        call stop_gitm("Stopping in get_potential")
     endif

     do iLat=-1,nLats+2
        do iLon=-1,nLons+2
           if (ElectronAverageEnergy(iLon,iLat) < 0.0) then
              ElectronAverageEnergy(iLon,iLat) = 0.1
              write(*,*) "i,j Negative : ",iLon,iLat,&
                   ElectronAverageEnergy(iLon,iLat)
           endif
           if (ElectronAverageEnergy(iLon,iLat) > 100.0) then
              write(*,*) "i,j Positive : ",iLon,iLat,&
                   ElectronAverageEnergy(iLon,iLat)
              ElectronAverageEnergy(iLon,iLat) = 0.1
           endif
        enddo
     enddo

     call IO_GetEFlux(ElectronEnergyFlux, iError)
     if (iError /= 0) then
        write(*,*) "Error in get_potential (IO_GetEFlux):"
        write(*,*) iError
        call stop_gitm("Stopping in get_potential")
     endif

     if (iDebugLevel > 2) &
          write(*,*) "==> Max, electron_ave_ene : ", &
          maxval(ElectronAverageEnergy), &
          maxval(ElectronEnergyFlux)

     IsFirstAurora(iBlock) = .false.

  endif

  if (iDebugLevel > 1) &
       write(*,*) "==> Min, Max, CPC Potential : ", &
       minval(Potential(:,:,:,iBlock))/1000.0, &
       maxval(Potential(:,:,:,iBlock))/1000.0, &
       (maxval(Potential(:,:,:,iBlock))-minval(Potential(:,:,:,iBlock)))/1000.0

  call end_timing("Get AMIE Potential")

end subroutine get_potential
