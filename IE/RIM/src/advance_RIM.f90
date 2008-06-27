
subroutine advance_RIM

  use ModRIM
  use ModParamRIM
  use ModIoRIM

  implicit none

  integer :: iLon, iLat, iError, iFile

  call IO_SetTime(CurrentTime)
  call IO_SetNorth

  call set_imf

  call distribute

  call get_conductance

  call IO_GetPotential(EmpiricalPotential, iError)

  if (iError /= 0) then
     write(*,*) "Error : ", iError
     call stop_RIM("Stopping in advance_RIM, call to IO_GetPotential")
  endif

  do iLon=0,nLons+1 
     nEmpiricalLats = 1
     do iLat=1,nLats
        if (abs(Latitude(iLon,iLat)) > HighLatBoundary) then
           potential(iLon,iLat) = EmpiricalPotential(iLon,nEmpiricalLats)
           if (iLon == 20) write(*,*) "pot : ",iLat,potential(iLon,iLat)
           nEmpiricalLats = nEmpiricalLats + 1
        endif
     enddo
  enddo

  nSolve=nSolve+1

  call conductance_gradients
  call solve
  call gather

  if(nSolve>0)then
     do iFile=1,nFile
        if (dn_output(iFile).gt.0.and.mod(nSolve,dn_output(iFile)).eq.0) then
           call write_output_RIM(iFile)
        endif
     end do
  end if

  write(*,*) "nSolve : ", nSolve

end subroutine advance_RIM
