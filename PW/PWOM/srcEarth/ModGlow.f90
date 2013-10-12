!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModGlow
  ! Module for Glowex and connected subroutines


  INTEGER JMAX, NBINS, LMAX, NMAJ, NEX, NW, NC, NST, NEI, NF
  PARAMETER (JMAX=92)
  PARAMETER (NBINS=84)
  PARAMETER (LMAX=59)
  PARAMETER (NMAJ=3)
  PARAMETER (NEX=20)
  PARAMETER (NW=20)
  PARAMETER (NC=10)
  PARAMETER (NST=6)
  PARAMETER (NEI=10)
  PARAMETER (NF=4)

  REAL  PHOTOTF(601),phototp(601)

  REAL PLG(9,4),CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC, &
       DAY,DF,DFA,APD,APDF,APT(4)
  INTEGER IYR
  INTEGER IFACTOR
  REAL EFLUX(NF), EZERO(NF), &
       SZA, DIP,  EFRAC, &
       ZO(JMAX), ZN2(JMAX), ZO2(JMAX), ZNO(JMAX), &
       ZNS(JMAX), ZND(JMAX), ZRHO(JMAX), ZE(JMAX), &
       ZCOL(NMAJ,JMAX),ZTN(JMAX),ZMAJ(NMAJ,JMAX),ZZ(JMAX), &
       ZTI(JMAX), ZTE(JMAX), &
       WAVE1(LMAX), WAVE2(LMAX), SFLUX(LMAX), &
       ENER(NBINS), DEL(NBINS), PHITOP(NBINS), &
       PESPEC(NBINS,JMAX), SESPEC(NBINS,JMAX), &
       PHOTOI(NST,NMAJ,JMAX),PHOTOD(NST,NMAJ,JMAX),PHONO(NST,JMAX), &
       QTI(JMAX),AURI(NMAJ,JMAX),PIA(NMAJ,JMAX),SION(NMAJ,JMAX), &
       UFLX(NBINS,JMAX), DFLX(NBINS,JMAX), AGLW(NEI,NMAJ,JMAX), &
       EHEAT(JMAX), TEZ(JMAX), ECALC(JMAX), &
       ZXDEN(NEX,JMAX), ZETA(NW,JMAX), ZCETA(NC,NW,JMAX),VCB(NW), &
       NNN(NMAJ)
  INTEGER ISCALE,JLOCAL,IERR
  
  real    :: DrGlow = 2.0e6
  integer :: nCellGlow = 390
  
contains
  subroutine get_ionization(nCell, AltD_C, IonRateO_C)
    implicit none
    integer, intent(in) :: nCell
    real,    intent(in) :: AltD_C(nCell)
    real,    intent(out):: IonRateO_C(nCell) 
    integer :: iCell, iMin,iMax
    real    :: x, Dx1, Dx2
    !------------------------------------------------------------------------
    call glowex
    do iCell = 1,nCell
       x = ( AltD_C(iCell) - AltD_C(1) ) / DrGlow + 1 
       iMin = floor  (x)
       iMax = ceiling(x)
       if (iMax >= nCellGlow) then
          IonRateO_C(iCell) = &
               Phototf(nCellGlow-2)
          cycle
       endif
       ! Set Interpolation weights
       Dx1 = x   - iMin ; Dx2 = 1.0 - Dx1
       
       ! Interpolate
       IonRateO_C(iCell) = &
            Dx2 * Phototf(iMin) + Dx1 * Phototf(iMax)

    enddo
        
  end subroutine get_ionization
  
end Module ModGlow
