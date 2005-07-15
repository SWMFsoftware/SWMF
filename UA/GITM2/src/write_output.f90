
subroutine write_output

  use ModTime
  use ModInputs
  use ModGITM

  implicit none

  real, external :: get_timing
  integer :: i, nMLTsTmp,nLatsTmp, iBlock
  logical :: IsDone

  if (floor((tSimulation-dt)/DtReport) /= &
       floor((tsimulation)/DtReport) .and. iDebugLevel >= 0) then
     if (IsFramework) then
        if(iProc==0)write(*,"(a,i6,a,3i2.2)") "UA:GITM2 iStep ", iStep, &
             ", Time : ",iTimeArray(4:6)
     else
        write(*,"(a,i6,a,3i2.2,a,f10.2,a)") "iStep ", iStep, &
             ", Time : ",iTimeArray(4:6), &
             ", RealTime : ",get_timing("GITM")/60.0," min"
     endif
  endif

  IsDone = .false.
  do i = 1, nOutputTypes
     if (floor((tSimulation-dt)/DtPlot(i)) /= &
          floor((tsimulation)/DtPlot(i)) .or. tSimulation == 0.0) then
        if (.not. IsDone) then
           if (.not.UseApex) call UA_calc_electrodynamics(nMLTsTmp, nLatsTmp)
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

