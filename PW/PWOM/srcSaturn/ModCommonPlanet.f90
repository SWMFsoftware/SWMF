!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModCommonPlanet
  use ModParameters
!  character(6) NamePlanet
  character(len=100), parameter :: NamePlanet = 'SATURN'
  
  integer,parameter :: nVar=16
  real,   parameter :: rLowerBoundary  = 6.1668e9
  real,   parameter :: rPlanet = 60268.E5

  REAL XH2(MaxGrid),XH(MaxGrid),XH2O(MaxGrid),XCH4(MaxGrid), &
       XTN(MAXGRID)
  
  REAL CLHpH3p(MaxGrid),CLELHp(MaxGrid),CLELH3p(MaxGrid), &
       CLHpH(MaxGrid),CLELH(MaxGrid),CLH2pH3p(MaxGrid), &
       CLH2pHp(MaxGrid),CLELH2p(MaxGrid)
  REAL FFHpp1(MaxGrid),FFHpp3(MaxGrid),FFHpp4(MaxGrid), &
       FFHpc2(MaxGrid),FFHpc3(MaxGrid), FFHpc9(MaxGrid),&
       FFHpc8(MaxGrid),FFHpr1(MaxGrid),FFH3pc1(MaxGrid),FFH3pc2(MaxGrid), &
       FFH3pc6(MaxGrid),FFH3pc7(MaxGrid),FFH3pr2(MaxGrid), &
       FFH2pp2(MaxGrid),FFH2pc9(MaxGrid),FFH2pc1(MaxGrid)


  REAL HLPE,HLPE0,HLPion1,HLPion2,HLPion3,HLPHE
  
  REAL jp1,jp2,jp3,jp4,kc1,kc2,kc3,kc6,kc7,kc8,kr1,kr2
  REAL kc9(MaxGrid)
  
  integer, parameter :: nIon = 4, nSpecies=8
  integer, parameter :: Ion1_ = 1, &    !Saturn:H3+   
       Ion2_ = 2, &    !Saturn:H+
       Ion3_ = 3, &    !Saturn:H2+
       Neutral1_= 5, & !Saturn:H2
       Neutral2_= 6, & !Saturn:H
       Neutral3_= 7, & !Saturn:H2O
       Neutral4_= 8    !Saturn:CH4
  
  ! named state variables
  integer, parameter :: RhoH3_=1, uH3_=2,  pH3_=3,  Th3_=4, &
       RhoH_=5, uH_=6,  pH_=7,  Th_=8, &
       RhoH2_=9,uH2_=10,pH2_=11,Th2_=12, &
       RhoE_=13,uE_=14, pE_=15, Te_=16
  

  integer,parameter :: iRho_I(nIon) = (/1,5,9,13/),&
       iU_I(nIon)   = (/2,6,10,14/),&
       iP_I(nIon)   = (/3,7,11,15/),&
       iT_I(nIon)   = (/4,8,12,16/)
  
  integer,parameter :: iRho(3) = (/1,4,7/),&
                       iU(3)   = (/2,5,8/),&
                       iP(3)   = (/3,6,9/),iTe_=10

  ! Dummy values for neutral parameters needed for ModAurora
!  integer, parameter :: O_=1, O2_=2, N2_=3, H_=4, He_=5 ! orig
  integer, parameter :: O_=1, O2_=2, N2_=3, He_=4 ! new
  ! Actual values for neutral parameters needed for plotting
  integer, parameter :: nNeutral = 4, H2_=1, H_=2, H2O_=3, CH4_=4
  real :: NDensity_CI(MaxGrid,nNeutral),NeutralPressure_C(MaxGrid)
  
  ! For Plotting
  integer, parameter :: nPlotVar = 20
  character(len=100), parameter :: NamePlotVar= &
 'r Lat Lon uH3 uH uH2 ue lgnH3 lgnH lgnH2 lgne TH3 TH TH2 Te MH3 MH MH2 Me Ef Pe g'
  integer, parameter :: nPlotVarNeutral = 6
  character(len=79), parameter :: NamePlotVarNeutral= &
 'r Lat Lon [H2] [H] [H2O] [CH4] g'


end Module ModCommonPlanet

