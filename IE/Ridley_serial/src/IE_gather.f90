
subroutine IE_gather

  use ModMpi
  use ModProcIE
  use ModIonosphere

  implicit none

  real :: localVar(2*IONO_nTheta-1, IONO_nPsi)
  integer :: iError, iSize
  !--------------------------------------------------------------------------
  IONO_phi = -1.0e32

  if (iProc == 0) then
     IONO_phi(1:IONO_nTheta,:) = IONO_NORTH_Phi
  endif

  if (iProc == nProc-1) then
     IONO_phi(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_Phi
  endif

  if (nProc > 1) then
     localVar = IONO_Phi
     iSize = (2*IONO_nTheta-1) * IONO_nPsi
     iError = 0
     call MPI_REDUCE(localVar, IONO_Phi, iSize, MPI_REAL, MPI_MAX, &
          0, iComm, iError)
  endif

end subroutine IE_gather
