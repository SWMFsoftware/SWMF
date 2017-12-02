!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine IE_gather

  use ModMpi
  use ModProcIE
  use ModIonosphere

  implicit none

  integer :: iError, n
  !--------------------------------------------------------------------------
  IONO_phi        = -Huge(1.0)
  IONO_IonNumFlux = -Huge(1.0)
  IONO_Joule      = -Huge(1.0)
  IONO_Jr         = -Huge(1.0)
  IONO_Ave_E      = -Huge(1.0)
  IONO_Eflux      = -Huge(1.0)
  IONO_SigmaP     = -Huge(1.0)
  IONO_SigmaH     = -Huge(1.0)

  if (iProc == 0) then
     IONO_phi(1:IONO_nTheta,:)        = IONO_NORTH_Phi
     IONO_IonNumFlux(1:IONO_nTheta,:) = IONO_NORTH_IonNumFlux
     IONO_Joule(1:IONO_nTheta,:)      = IONO_NORTH_Joule
     IONO_Jr(1:IONO_nTheta,:)         = IONO_NORTH_Jr
     IONO_Ave_E(1:IONO_nTheta,:)      = IONO_NORTH_Ave_E
     IONO_Eflux(1:IONO_nTheta,:)      = IONO_NORTH_EFlux
     IONO_SigmaP(1:IONO_nTheta,:)     = IONO_NORTH_SigmaP
     IONO_SigmaH(1:IONO_nTheta,:)     = IONO_NORTH_SigmaH
  endif

  if (iProc == nProc-1) then
     IONO_phi(IONO_nTheta:2*IONO_nTheta-1,:)        = IONO_SOUTH_Phi
     IONO_IonNumFlux(IONO_nTheta:2*IONO_nTheta-1,:) = IONO_SOUTH_IonNumFlux
     IONO_Joule(IONO_nTheta:2*IONO_nTheta-1,:)      = IONO_SOUTH_Joule
     IONO_Jr(IONO_nTheta:2*IONO_nTheta-1,:)         = IONO_SOUTH_Jr
     IONO_Ave_E(IONO_nTheta:2*IONO_nTheta-1,:)      = IONO_SOUTH_Ave_E
     IONO_Eflux(IONO_nTheta:2*IONO_nTheta-1,:)      = IONO_SOUTH_EFlux
     IONO_SigmaP(IONO_nTheta:2*IONO_nTheta-1,:)     = IONO_SOUTH_SigmaP
     IONO_SigmaH(IONO_nTheta:2*IONO_nTheta-1,:)     = IONO_SOUTH_SigmaH
  endif

  if (nProc > 1) then
     n = size(IONO_Phi)
     call MPI_reduce_real_array(IONO_Phi,        n, MPI_MAX, 0, iComm, iError)
     call MPI_reduce_real_array(IONO_IonNumFlux, n, MPI_MAX, 0, iComm, iError)
     call MPI_reduce_real_array(IONO_Joule,      n, MPI_MAX, 0, iComm, iError)
     call MPI_reduce_real_array(IONO_Jr,         n, MPI_MAX, 0, iComm, iError)
     call MPI_reduce_real_array(IONO_Ave_E,      n, MPI_MAX, 0, iComm, iError)
     call MPI_reduce_real_array(IONO_Eflux,      n, MPI_MAX, 0, iComm, iError)
     call MPI_reduce_real_array(IONO_SigmaP,     n, MPI_MAX, 0, iComm, iError)
     call MPI_reduce_real_array(IONO_SigmaH,     n, MPI_MAX, 0, iComm, iError)

  endif

end subroutine IE_gather
