
module ModCommonVariables

  use ModCommonPlanet
  
  use ModParameters
  
  integer iUnitInput, iUnitOutput, iUnitGraphics, iUnitRestart, &
       iUnitCollision, iUnitSourceGraphics
  
  

  
  real :: wHorizontal
  character(len=7) :: TypeSolver
  logical :: DoLog
  real :: uJoule2
  real :: AR12top(2),AR23top(2),CellVolumeTop(2)

  real :: cMax_O(MaxGrid), cMax_H(MaxGrid),cMax_e(MaxGrid),&
       cfl(MaxGrid), cfl_O(MaxGrid), cfl_H(MaxGrid),cfl_e(MaxGrid),&
       MaxCfl,cfl_dt
  
  REAL CELLNW(nVar,MaxGrid)
  REAL AR12(MaxGrid),AR23(MaxGrid),DAREA(MaxGrid), CellVolume_C(MaxGrid)
  integer NEXP
  
  REAL ALTD(MaxGrid),RAD(MaxGrid),RBOUND(MaxGrid)
  ! plasma parameters, U=velocity, D=density, P=pressure, T=electrons
  ! the parameter is followed by the quantity e.g. TOXY=temp of oxygen

  real :: State_GV(-1:MaxGrid,nVar),StateOld_GV(-1:MaxGrid,nVar),&
          SoundSpeed_GI(0:MaxGrid,nIon),Source_CV(MaxGrid,nVar)
  real :: HeatCon_GI(0:maxGrid,nIon)
     
  REAL EFIELD(MaxGrid),GRAVTY(MaxGrid),CURR(MaxGrid),EfieldConstant
  REAL Centrifugal(MaxGrid)
  real :: AveIE, EfluxIE
  !scaling reference for peak topside heat flux  in aurora during quiet times
  real, parameter  :: EtopAurora = 2.0e-3, EfluxRef = 7.0  
  real, parameter  :: EtopMin = 5e-4   !minimum heat flux (due to polar rain)
  real, parameter  :: EtopPhotoElectrons = 1e-3   !Peak Heat Flux from photoelectrons
  Logical :: UsePhotoElectronHeatFlux = .true., UseAuroralHeatFlux = .true., &
             UseCuspHeatFlux = .true.

  REAL QOXYG(MaxGrid),QHEL(MaxGrid),QHYD(MaxGrid),&
       QELECT(MaxGrid)
  REAL ELFXIN

  REAL  CZHN2,CZHO2,CZHO,CZHOX,CZHEN2,CZHEO2,CZHEHE,&
       CZHEO,CZHEH,CZHEOX,CZHEHD,XTNMAX

  REAL  ETOP,CURTIM0,CURTIM,CURRMN,CURRMX

  REAL  UBGNDO,PBGNDO,DBGNDO,TBGNDO,WBGNDO,UBGNDH,PBGNDH,&
       DBGNDH,TBGNDH,WBGNDH,UBGNHE,PBGNHE,DBGNHE,TBGNHE,WBGNHE,UBGNDE,&
       PBGNDE,DBGNDE,TBGNDE,WBGNDE,&
       UBGNDO2,PBGNDO2,DBGNDO2,TBGNDO2,&
       UBGNDH2,PBGNDH2,DBGNDH2,TBGNDH2,&
       UBGNHE2,PBGNHE2,DBGNHE2,TBGNHE2

  REAL  ALTMIN,ALTMAX,DTR1,DTR2
  real :: DrBnd =2.0e6

  REAL  GAMMA,GMIN1,GMIN2,GPL1,GPL2,GM12,GRAR,GREC
  REAL  CPO,CPH,CPHE,CPE,CVO,CVH,CVHE,CVE
  REAL  XAMU
  real :: Mass_I(nIon),MassElecIon_I(nIon-1),RGAS_I(nIon)
  INTEGER NDIM,NDIM1,NSTEP,NPRINT,NSTPMX,NDIM2,NDIMM
  REAL  DT,DTMX,TIME,TMAX,DTX1,DTX2
  INTEGER NPT1,NPT2,NPT3,NPT4,NPT5,NPT6,NCL,NTS
  REAL  H0,H1E1,H1O1,H1H1,H1E2,H1O2,H1H2,H2,H3,H4 &
       ,HLPE,HLPE0,HLPO,HLPH,HLPHE
  Logical :: UseStaticAtmosphere=.false.
  INTEGER :: IYD=76183,IART=1
  REAL    ::UT,SEC,GLAT,GLONG,STL,F107A=60.,F107=60.,GMLAT,GMLONG
  real    :: AP(7)=(/4.,4.,4.,4.,4.,4.,4./)
  REAL  IonRateO_C(MaxGrid)
!  REAL  TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0, &
!       G0,RL,DD,DB14
!  REAL PTM(8),PDM(8,7)
!  INTEGER ISW
!  REAL SW(25),SWC(25)
!  REAL TINFG,GB,ROUT,TT(15)
  REAL GSURF,RE,Omega
  
  real :: CollisionFreq_IIC(nIon, nSpecies,MaxGrid)
  real :: HeatFlowCoef_II(nIon, nSpecies)
  real :: FricHeatCoef_II(nIon, nSpecies)
  real :: MassFracCoef_II(nIon, nSpecies)


       
end module ModCommonVariables
