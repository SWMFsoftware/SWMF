
subroutine gather

  use ModRIM
  use ModMpi
  use ModProcIE
  use ModParamRIM, only: iDebugLevel

  implicit none

  integer :: iError, iSize

  ! Let's try message passing some of this stuff to Processor 0

  PotentialAll = -1.0e32

  ! Fill Arrays
  PotentialAll(iProc*nLons+1:iProc*nLons+nLons, :) = Potential

  if (iProc == 0) then
     PotentialAll(0, :) = Potential(0,:)
  endif

  if (iProc == nProc-1) then
     PotentialAll((iProc+1)*nLons+1, :) = Latitude(nLons+1,:)
  endif

  iSize = nLats * (nProc*nLons+2)

  localVar = PotentialAll
  call MPI_REDUCE(localVar, PotentialAll, iSize, MPI_REAL, MPI_MAX, &
       0, iComm, iError)

  if (iProc == 0) then
     cpcpn = maxval(PotentialAll(:,0:nLats/2)) - &
             minval(PotentialAll(:,0:nLats/2))
     cpcps = maxval(PotentialAll(:,0:nLats/2)) - &
             minval(PotentialAll(:,0:nLats/2))
  endif

  if (iDebugLevel > 1) then 
     write(*,*) "CPCP of total solution : ", &
          (maxval(PotentialAll) - minval(PotentialAll))/1000.0, " kV"
     write(*,*) "CPCP South/North : ", cpcps/1000.0, cpcpn/1000.0, " kV"
  endif

end subroutine gather
