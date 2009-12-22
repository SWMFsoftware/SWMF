
subroutine get_log_info(GlobalMinTemp, GlobalMaxTemp, &
     GlobalMinVertVel, GlobalMaxVertVel, AverageTemp, AverageVertVel, &
     TotalVolume)

  use ModGITM

  real, intent(out) :: GlobalMinTemp, GlobalMaxTemp
  real, intent(out) :: GlobalMinVertVel, GlobalMaxVertVel
  real, intent(out) :: AverageTemp, AverageVertVel
  real, intent(out) :: TotalVolume

  integer :: iBlock, iSpecies
  !--------------------------------------------------------------------------

  GlobalMaxTemp    = 0.0
  GlobalMinTemp    = 1.0e32
  GlobalMaxVertVel = 0.0
  GlobalMinVertVel = 1.0e32
  AverageTemp      = 0.0
  AverageVertVel   = 0.0
  TotalVolume    = 0.0

  do iBlock = 1, nBlocks

     call calc_rates(iBlock)

     GlobalMaxTemp = max(GlobalMaxTemp, &
          maxval(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)* &
          TempUnit(1:nLons,1:nLats,1:nAlts)))
     GlobalMinTemp = min(GlobalMinTemp, &
          minval(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)* &
          TempUnit(1:nLons,1:nLats,1:nAlts)))
     AverageTemp = AverageTemp + &
          sum(Temperature(1:nLons,1:nLats,1:nAlts,iBlock) &
          *   TempUnit(1:nLons,1:nLats,1:nAlts) &
          *   CellVolume(1:nLons,1:nLats,1:nAlts,iBlock))

     GlobalMaxVertVel = max(GlobalMaxVertVel, &
          maxval(VerticalVelocity(1:nLons,1:nLats,1:nAlts,:,iBlock)))
     GlobalMinVertVel = min(GlobalMinVertVel, &
          minval(VerticalVelocity(1:nLons,1:nLats,1:nAlts,:,iBlock)))
     do iSpecies = 1, nSpecies
        AverageVertVel = AverageVertVel + &
             sum(VerticalVelocity(1:nLons,1:nLats,1:nAlts,iSpecies,iBlock) &
             *   CellVolume(1:nLons,1:nLats,1:nAlts,iBlock))
     enddo

     TotalVolume = TotalVolume + &
          sum(CellVolume(1:nLons,1:nLats,1:nAlts,iBlock))

  enddo

  AverageTemp    = AverageTemp
  AverageVertVel = AverageVertVel / nSpecies

  

end subroutine get_log_info

!==============================================================================

subroutine logfile(dir)

  use ModGITM
  use ModInputs
  use ModTime
  use ModMpi
  use ModIndices
  use ModIndicesInterfaces
  use ModUtilities, ONLY: flush_unit

  implicit none

  character (len=*), intent(in) :: dir
  character (len=8) :: cIter

  real    :: minTemp, maxTemp, localVar, minVertVel, maxVertVel
  real    :: AverageTemp, AverageVertVel, TotalVolume, Bx, By, Bz, Vx, Hpi
  integer :: iError

  if (.not. IsOpenLogFile .and. iProc == 0) then

     IsOpenLogFile = .true.
     call CON_io_unit_new(iLogFileUnit_)

     write(cIter,"(i8.8)") iStep

     open(unit=iLogFileUnit_, &
          file=dir//"/log"//cIter//".dat",status="replace")

     write(iLogFileUnit_,'(a)') "GITM2 log file"
     write(iLogFileUnit_,'(a,L2)') "## Inputs from UAM.in" 
      write(iLogFileUnit_,'(a,L2)') "# Resart=", dorestart
     write(iLogFileUnit_,'(4(a,f9.3))') "# Eddy coef: ", EddyDiffusionCoef, &
          " Eddy P0: ",EddyDiffusionPressure0,&
          " Eddy P1: ",EddyDiffusionPressure1,&
          " Eddy Scaling: ",EddyScaling
     write(iLogFileUnit_,'(2(a,L2))') "# Statistical Models Only: ",usestatisticalmodelsonly,&
          " Apex: ",useApex
     if (useEUVdata) then
        write(iLogFileUnit_,'(a,L2,a)') "# EUV Data: ",useEUVdata, "File: ", cEUVFile
     else
        write(iLogFileUnit_,'(a,L2)') "# EUV Data: ",useEUVdata
     endif
      write(iLogFileUnit_,'(a,a15)') "# AMIE: ", cAmieFileNorth,cAmieFileSouth
      write(iLogFileUnit_,'(3(a,L2))') "# Solar Heating: ",useSolarHeating ,&
          " Joule Heating: ",useJouleHeating ,&
          " Auroral Heating: ", useAuroralHeating
       write(iLogFileUnit_,'(2(a,L2))') "# NO Cooling: ", useNOCooling, &
          " O Cooling: ", useOCooling
       write(iLogFileUnit_,'(3(a,L2))') "# Conduction: ",useConduction ,&
          " Turbulent Conduction: ", useTurbulentCond,&
          " Updated Turbulent Conduction: ",useUpdatedTurbulentCond
       write(iLogFileUnit_,'(3(a,L2))') "# Pressure Grad: ",usePressureGradient ,&
          " Ion Drag: ", useIonDrag,&
          " Neutral Drag: ", useNeutralDrag
       write(iLogFileUnit_,'(3(a,L2))') "# Viscosity: ", useViscosity,&
          " Coriolis: ", useCoriolis,&
          " Gravity: ",useGravity 
       write(iLogFileUnit_,'(3(a,L2))') "# Ion Chemistry: ",useIonChemistry ,&
          " Ion Advection: ",useIonAdvection ,&
          " Neutral Chemistry: ", useNeutralChemistry
          
       write(iLogFileUnit_,'(a)') " "
       write(iLogFileUnit_,'(a)') "#START"
       write(iLogFileUnit_,'(a)') &
            "   iStep yyyy mm dd hh mm ss  ms      dt "// &
            "min(T) max(T) mean(T) min(VV) max(VV) mean(VV) F107 F107A "// &
            "By Bz Vx HPI"
          

  endif

  call get_log_info(MinTemp, MaxTemp, MinVertVel, MaxVertVel, &
       AverageTemp, AverageVertVel, TotalVolume)

  localVar = TotalVolume
  call MPI_REDUCE(localVar, TotalVolume, 1, MPI_REAL, MPI_SUM, &
       0, iCommGITM, iError)

  localVar = MinTemp
  call MPI_REDUCE(localVar, minTemp, 1, MPI_REAL, MPI_MIN, &
       0, iCommGITM, iError)

  localVar = MaxTemp
  call MPI_REDUCE(localVar, maxTemp, 1, MPI_REAL, MPI_MAX, &
       0, iCommGITM, iError)

  localVar = MinVertVel
  call MPI_REDUCE(localVar, minVertVel, 1, MPI_REAL, MPI_MIN, &
       0, iCommGITM, iError)

  localVar = MaxVertVel
  call MPI_REDUCE(localVar, maxVertVel, 1, MPI_REAL, MPI_MAX, &
       0, iCommGITM, iError)

  LocalVar = AverageTemp
  call MPI_REDUCE(LocalVar, AverageTemp, 1, MPI_REAL, MPI_SUM, &
       0, iCommGITM, iError)

  LocalVar = AverageVertVel
  call MPI_REDUCE(LocalVar, AverageVertVel, 1, MPI_REAL, MPI_SUM, &
       0, iCommGITM, iError)

  if (iProc == 0) then

     AverageTemp = AverageTemp / TotalVolume
     AverageVertVel = AverageVertVel / TotalVolume

     call get_f107(CurrentTime, f107, iError)
     call get_f107A(CurrentTime, f107A, iError)
      call get_IMF_By(CurrentTime, by, iError)
      call get_IMF_Bz(CurrentTime, bz, iError)
      call get_sw_v(CurrentTime, Vx, iError)
      call get_hpi(CurrentTime,Hpi,iError)

     write(iLogFileUnit_,"(i8,i5,5i3,i4,f8.4,6f13.5,6f9.1)") &
          iStep, iTimeArray, dt, minTemp, maxTemp, AverageTemp, &
          minVertVel, maxVertVel, AverageVertVel,&
         f107,f107A,By,Bz,Vx, Hpi

     call flush_unit(iLogFileUnit_)
  endif

end subroutine logfile
