
subroutine IE_gather

  use ModMpi
  use ModProcIE
  use ModIonosphere

  implicit none

  real :: localVar(2*IONO_nTheta-1, IONO_nPsi)
  integer :: iError, iSize
  !--------------------------------------------------------------------------
  IONO_phi = -1.0e32
  IONO_IonNumFlux = -1.0e32
  IONO_Joule = -1.0e32
  IONO_Jr = -1.0e32

  if (iProc == 0) then
     IONO_phi(1:IONO_nTheta,:) = IONO_NORTH_Phi
     IONO_IonNumFlux(1:IONO_nTheta,:) = IONO_NORTH_IonNumFlux
     IONO_Joule(1:IONO_nTheta,:) = IONO_NORTH_Joule
     IONO_Jr(1:IONO_nTheta,:) = IONO_NORTH_Jr
  endif

  if (iProc == nProc-1) then
     IONO_phi(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_Phi
     IONO_IonNumFlux(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_IonNumFlux
     IONO_Joule(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_Joule
     IONO_Jr(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_Jr
  endif

  if (nProc > 1) then

     iSize = (2*IONO_nTheta-1) * IONO_nPsi

     localVar = IONO_Phi
     iError = 0
     call MPI_REDUCE(localVar, IONO_Phi, iSize, MPI_REAL, MPI_MAX, &
          0, iComm, iError)

     localVar = IONO_IonNumFlux
     iError = 0
     call MPI_REDUCE(localVar, IONO_IonNumFlux, iSize, MPI_REAL, MPI_MAX, &
          0, iComm, iError)

     localVar = IONO_Joule
     iError = 0
     call MPI_REDUCE(localVar, IONO_Joule, iSize, MPI_REAL, MPI_MAX, &
          0, iComm, iError)
     
     localVar = IONO_Jr
     iError = 0
     call MPI_REDUCE(localVar, IONO_Jr, iSize, MPI_REAL, MPI_MAX, &
          0, iComm, iError)
     
  endif

end subroutine IE_gather
