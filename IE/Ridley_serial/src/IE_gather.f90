
subroutine IE_gather

  use ModMpi
  use ModProcIE
  use ModIonosphere

  implicit none

  real, dimension(2*IONO_nTheta-1, IONO_nPsi):: localVar, localVar1,localVar2
  integer :: iError, iSize
  !--------------------------------------------------------------------------
  IONO_phi = -1.0e32

  if (iProc == 0) then
     IONO_phi(1:IONO_nTheta,:) = IONO_NORTH_Phi
     IONO_IonNumFlux(1:IONO_nTheta,:) = IONO_NORTH_IonNumFlux
     IONO_Joule(1:IONO_nTheta,:) = IONO_NORTH_Joule
  endif

  if (iProc == nProc-1) then
     IONO_phi(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_Phi
     IONO_IonNumFlux(1:IONO_nTheta,:) = IONO_SOUTH_IonNumFlux
     IONO_Joule(1:IONO_nTheta,:) = IONO_SOUTH_Joule
  endif

  if (nProc > 1) then
     localVar = IONO_Phi
     localVar1 = IONO_IonNumFlux
     localVar2 = IONO_Joule
     iSize = (2*IONO_nTheta-1) * IONO_nPsi
     iError = 0
     call MPI_REDUCE(localVar, IONO_Phi, iSize, MPI_REAL, MPI_MAX, &
          0, iComm, iError)
     call MPI_REDUCE(localVar1, IONO_IonNumFlux, iSize, MPI_REAL, MPI_MAX, &
          0, iComm, iError)
     call MPI_REDUCE(localVar2, IONO_Joule, iSize, MPI_REAL, MPI_MAX, &
          0, iComm, iError)
     
  endif

end subroutine IE_gather
