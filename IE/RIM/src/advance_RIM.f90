
subroutine advance_RIM

  use ModRIM
  use ModParamRIM
  use ModIoRIM

  implicit none

  integer :: iLon, iLat, iError, iFile, iT1, iT2

  iError = 0

  call IO_SetTime(CurrentTime)
  call IO_SetNorth

  if (maxval(Latitude) > HighLatBoundary) call set_imf

  call distribute

  if (nSolve == 0 .and. maxval(Jr) == 0.0) then
     potential = 0.0
     call gather
     return
  endif

  if (maxval(Latitude) > HighLatBoundary) &
       call IO_GetPotential(EmpiricalPotential, iError)

  if (iError /= 0) then
     write(*,*) "Error : ", iError
     call stop_RIM("Stopping in advance_RIM, call to IO_GetPotential")
  endif

  do iLon=0,nLons+1 
     nEmpiricalLats = 0
     do iLat=1,nLats
        if (abs(Latitude(iLon,iLat)) > HighLatBoundary) then
           nEmpiricalLats = nEmpiricalLats + 1
           potential(iLon,iLat) = EmpiricalPotential(iLon,nEmpiricalLats)
        endif
     enddo
  enddo

  call get_conductance
  call conductance_gradients
  call solve
  call gather

  if(nSolve>0)then
     do iFile=1,nFile
        if (dt_output(iFile) > 0) then
           iT1 = floor((CurrentTime-StartTime)/dt_output(iFile))
           iT2 = floor((OldTime-StartTime)/dt_output(iFile))
        else
           iT1 = 0
           iT2 = 0
        endif
        if ((dn_output(iFile).gt.0.and.mod(nSolve,dn_output(iFile)).eq.0).or.&
             (iT1 > iT2)) then
           call write_output_RIM(iFile)
        endif
     end do
  end if

  nSolve=nSolve+1

end subroutine advance_RIM
