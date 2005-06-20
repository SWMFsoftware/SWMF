
subroutine advance

  use ModConstants
  use ModGITM
  use ModTime
  use ModEUV
  use ModInputs

  implicit none

  integer, external :: jday

  integer :: iBlock, iAlt, iLat, iLon,ispecies
  real*8 :: DTime

  call report("advance",1)
  call start_timing("advance")

  ! Godunov Splitting reverses the hor and ver each time step

  if (.not. UseStatisticalModelsOnly) then

        call add_sources
        if (.not. Is1D) call advance_horizontal_all
        call advance_vertical_all

  else

     Dt = DtStatisticalModels

  endif

  if (iDebugLevel > 0) write(*,*) "=> MaxTemp : ",maxval(temperature)*TempUnit(1,1,nalts)

  tSimulation = tSimulation + dt
  CurrentTime = StartTime + tSimulation

  call time_real_to_int(CurrentTime, iTimeArray)

  DTime = CurrentTime - VernalTime
  do while (DTime > SecondsPerYear)
     VernalTime = VernalTime+int(DaysPerYear)*Rotation_Period
     DTime = CurrentTime - VernalTime
  enddo
  iDay  = DTime / Rotation_Period
  uTime = (DTime / Rotation_Period - iDay) * Rotation_Period

  iJulianDay = jday(iTimeArray(1), iTimeArray(2), iTimeArray(3)) 
!  utime = iTimeArray(4)*3600.0 + iTimeArray(5)*60.0 + &
!       iTimeArray(6) + iTimeArray(7)/1000.0

  if (UseStatisticalModelsOnly) then
     call init_msis
     call init_iri
     call init_b0
  endif

  call end_timing("advance")

contains

  !==========================================================================
  subroutine advance_vertical_all

    call report("advance_vertical_all",1)
    call start_timing("vertical_all")

    do iBlock = 1, nBlocks

       do iLon = 1, nLons ; do iLat = 1, nLats
          call advance_vertical(iLon,iLat,iBlock)
       end do; end do

    end do

    call end_timing("vertical_all")

  end subroutine advance_vertical_all

  !==========================================================================
  subroutine advance_horizontal_all

    call report("advance_horizontal_all",1)
    call start_timing("horizontal_all")

    call exchange_messages_sphere

    do iBlock = 1, nBlocks

       call calc_physics(iBlock)
       call calc_rates(iBlock)

       call advance_horizontal(iBlock)

    end do

!    Two stage
!
!     do iBlock = 1, nBlocks
!        call calc_physics(iBlock)
!        ! call calc_rates(iBlock)
!        call advance_horizontal(iBlock)
!     end do
!
!     call exchange_messages_sphere

    call end_timing("horizontal_all")

  end subroutine advance_horizontal_all

end subroutine advance
