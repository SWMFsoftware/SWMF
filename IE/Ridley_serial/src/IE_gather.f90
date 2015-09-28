!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

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
  IONO_Ave_E = -1.0e32
  IONO_Eflux = -1.0e32
  IONO_SigmaP = -1.0e32
  IONO_SigmaH = -1.0e32

  if (iProc == 0) then
     IONO_phi(1:IONO_nTheta,:) = IONO_NORTH_Phi
     IONO_IonNumFlux(1:IONO_nTheta,:) = IONO_NORTH_IonNumFlux
     IONO_Joule(1:IONO_nTheta,:) = IONO_NORTH_Joule
     IONO_Jr(1:IONO_nTheta,:) = IONO_NORTH_Jr
     IONO_Ave_E(1:IONO_nTheta,:) = IONO_NORTH_Ave_E
     IONO_Eflux(1:IONO_nTheta,:) = IONO_NORTH_EFlux
     IONO_SigmaP(1:IONO_nTheta,:) = IONO_NORTH_SigmaP
     IONO_SigmaH(1:IONO_nTheta,:) = IONO_NORTH_SigmaH
  endif

  if (iProc == nProc-1) then
     IONO_phi(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_Phi
     IONO_IonNumFlux(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_IonNumFlux
     IONO_Joule(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_Joule
     IONO_Jr(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_Jr
     IONO_Ave_E(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_Ave_E
     IONO_Eflux(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_EFlux
     IONO_SigmaP(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_SigmaP
     IONO_SigmaH(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_SigmaH
  endif

  if (nProc > 1) then

     iSize = (2*IONO_nTheta-1) * IONO_nPsi

     localVar = IONO_Phi
     call MPI_REDUCE(localVar, IONO_Phi, iSize, MPI_REAL, MPI_MAX, &
          0, iComm, iError)

     localVar = IONO_IonNumFlux
     call MPI_REDUCE(localVar, IONO_IonNumFlux, iSize, MPI_REAL, MPI_MAX, &
          0, iComm, iError)

     localVar = IONO_Joule
     call MPI_REDUCE(localVar, IONO_Joule, iSize, MPI_REAL, MPI_MAX, &
          0, iComm, iError)
     
     localVar = IONO_Jr
     call MPI_REDUCE(localVar, IONO_Jr, iSize, MPI_REAL, MPI_MAX, &
          0, iComm, iError)

     localVar = IONO_Ave_E
     call MPI_REDUCE(localVar, IONO_Ave_E, iSize, MPI_REAL, MPI_MAX, &
          0, iComm, iError)

     localVar = IONO_Eflux
     call MPI_REDUCE(localVar, IONO_Eflux, iSize, MPI_REAL, MPI_MAX, &
          0, iComm, iError)     

     localVar = IONO_SigmaP
     call MPI_REDUCE(localVar, IONO_SigmaP, iSize, MPI_REAL, MPI_MAX, &
          0, iComm, iError)     

     localVar = IONO_SigmaH
     call MPI_REDUCE(localVar, IONO_SigmaH, iSize, MPI_REAL, MPI_MAX, &
          0, iComm, iError)     

  endif

end subroutine IE_gather
