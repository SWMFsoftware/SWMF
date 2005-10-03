
subroutine write_output

  use ModTime
  use ModInputs
  use ModGITM

  implicit none

  real, external :: get_timing
  integer :: i, nMLTsTmp,nLatsTmp, iBlock
  logical :: IsDone
  real    :: t
  character (len=4) :: sTimeUnit

  if (floor((tSimulation-dt)/DtReport) /= &
       floor((tsimulation)/DtReport) .and. iDebugLevel >= 0) then
     if (IsFramework) then
        if(iProc==0)write(*,"(a,i6,a,3i2.2)") "UA:GITM2 iStep ", iStep, &
             ", Time : ",iTimeArray(4:6)
     else
        t = get_timing("GITM")
        if (t < 120.0) then
           sTimeUnit = " sec"
        else 
           if (t < 7200.0) then
              t = t/60.0
              sTimeUnit = " min"
           else
              t = t/3600.0
              sTimeUnit = " hrs"
           endif
        endif
        write(*,"(a,i6,a,3i2.2,a,f10.2,a)") "iStep ", iStep, &
             ", Time : ",iTimeArray(4:6), &
             ", RealTime : ",t,sTimeUnit
     endif
  endif

  IsDone = .false.
  do i = 1, nOutputTypes
     if (floor((tSimulation-dt)/DtPlot(i)) /= &
          floor((tsimulation)/DtPlot(i)) .or. tSimulation == 0.0) then
        if (.not. IsDone) then
           if (.not.UseApex .and. .not.Is1D) call UA_calc_electrodynamics(nMLTsTmp, nLatsTmp)
           IsDone = .true.
        endif
        do iBlock = 1, nBlocks
           call output("UA/data/",iBlock, i)
        enddo
     endif
  enddo

  call move_satellites

  if (floor((tSimulation-dt)/DtRestart) /= &
       floor((tsimulation)/DtRestart)) then
     call write_restart("UA/restartOUT/")
  endif

end subroutine write_output

