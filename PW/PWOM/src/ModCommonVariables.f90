
module ModCommonVariables

  use ModCommonPlanet
  
  use ModParameters
  
  integer iUnitInput, iUnitOutput, iUnitGraphics, iUnitRestart, &
       iUnitCollision, iUnitSourceGraphics,iLine,nLine
  
  
  !    MSIS86 parameters
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
  
  Real wHorizontal
  CHARACTER(7) TypeSolver
  Logical IsImplicit,DoLog

  Real AR12top(2),AR23top(2),CellVolumeTop(2)

  REAL cMax_O(MaxGrid), cMax_H(MaxGrid),cMax_e(MaxGrid),&
       cfl(MaxGrid), cfl_O(MaxGrid), cfl_H(MaxGrid),cfl_e(MaxGrid),&
       MaxCfl,cfl_dt
  
  REAL CELLNW(nVar,MaxGrid)
  REAL AR12(MaxGrid),AR23(MaxGrid),DAREA(MaxGrid), CellVolume_C(MaxGrid)
  
  REAL ALTD(MaxGrid),RAD(MaxGrid),RBOUND(MaxGrid)
  ! plasma parameters, U=velocity, D=density, P=pressure, T=electrons
  ! the parameter is followed by the quantity e.g. TOXY=temp of oxygen
  REAL UOXYG(MaxGrid),DOXYG(MaxGrid),POXYG(MaxGrid),&
       TOXYG(MaxGrid),UHYD(MaxGrid),&
       DHYD(MaxGrid),PHYD(MaxGrid),THYD(MaxGrid),UHEL(MaxGrid),&
       DHEL(MaxGrid),PHEL(MaxGrid),&
       THEL(MaxGrid),UELECT(MaxGrid),DELECT(MaxGrid),PELECT(MaxGrid),&
       TELECT(MaxGrid)
     !Changed for Saturn      
     
     
  REAL EFIELD(MaxGrid),GRAVTY(MaxGrid),CURR(MaxGrid),EfieldConstant
  REAL Centrifugal(MaxGrid)
  REAL FCLSNE(MaxGrid),FCLSHE(MaxGrid),FCLSNH(MaxGrid),&
       FCLSNO(MaxGrid),&
       ECLSNE(MaxGrid),ECLSHE(MaxGrid),ECLSNH(MaxGrid),ECLSNO(MaxGrid),&
       ADMSO(MaxGrid),ADMSH(MaxGrid),ADMSHE(MaxGrid),ADMSE(MaxGrid),&
       YYH(MaxGrid)
      
  REAL TCONO(MaxGrid),TCONH(MaxGrid),TCONHE(MaxGrid),&
       TCONE(MaxGrid)
  REAL QOXYG(MaxGrid),QHEL(MaxGrid),QHYD(MaxGrid),&
       QELECT(MaxGrid)
  REAL ELFXIN
  !      REAL PHOTOTF(MaxGrid)

!GABOR These variables are never used !!!
  CHARACTER*4 NAME(2),ISDATE(3),ISTIME(2)

!    IsRestart = .true. indicates restart from a file, 
!    otherwise no file is read
  logical :: IsRestart
  logical IsVariableDt
  


     
  REAL  TCSFO,TCSFHE,TCSFH,TCSFE,TCBGO,TCBGHE,TCBGH,TCBGE
  REAL  CZHN2,CZHO2,CZHO,CZHOX,CZHEN2,CZHEO2,CZHEHE,&
       CZHEO,CZHEH,CZHEOX,CZHEHD,XTNMAX

  REAL  ETOP,CURTIM0,CURTIM,CURRMN,CURRMX
  REAL  USURFO,PSURFO,DSURFO,TSURFO,WSURFO,USURFH,PSURFH,&
       DSURFH,TSURFH,WSURFH,USURHE,PSURHE,DSURHE,TSURHE,WSURHE,USURFE,&
       PSURFE,DSURFE,TSURFE,WSURFE
  REAL  UBGNDO,PBGNDO,DBGNDO,TBGNDO,WBGNDO,UBGNDH,PBGNDH,&
       DBGNDH,TBGNDH,WBGNDH,UBGNHE,PBGNHE,DBGNHE,TBGNHE,WBGNHE,UBGNDE,&
       PBGNDE,DBGNDE,TBGNDE,WBGNDE,&
       UBGNDO2,PBGNDO2,DBGNDO2,TBGNDO2,&
       UBGNDH2,PBGNDH2,DBGNDH2,TBGNDH2,&
       UBGNHE2,PBGNHE2,DBGNHE2,TBGNHE2

  REAL  ALTMIN,ALTMAX,DRBND,DTR1,DTR2

  REAL  GAMMA,GMIN1,GMIN2,GPL1,GPL2,GM12,GRAR,GREC
  REAL  RGASO,RGASH,RGASHE,RGASE,CPO,CPH,CPHE,CPE,CVO,CVH,CVHE,CVE
  REAL  XAMU,XMSO,XMSH,XMSHE,XMSE,RTOXEL,RTHDEL,RTHEEL
  INTEGER NDIM,NDIM1,NSTEP,NPRINT,NSTPMX,NDIM2,NDIMM
  REAL  DT,DTMX,TIME,TMAX,DTX1,DTX2
  INTEGER NPT1,NPT2,NPT3,NPT4,NPT5,NPT6,NCL,NTS
  REAL  H0,H1E1,H1O1,H1H1,H1E2,H1O2,H1H2,H2,H3,H4 &
       ,HLPE,HLPE0,HLPO,HLPH,HLPHE
  INTEGER IYD,IART
  REAL  UT,SEC,GLAT,GLONG,STL,F107A,F107,AP(7),GMLAT,GMLONG
  REAL  PHOTOTF(MaxGrid),phototp(MaxGrid)
  REAL  TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0, &
       G0,RL,DD,DB14
  REAL PTM(8),PDM(8,7)
  INTEGER ISW
  REAL SW(25),SWC(25)
  REAL TINFG,GB,ROUT,TT(15)
  REAL GSURF,RE,Omega
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




       
end module ModCommonVariables
