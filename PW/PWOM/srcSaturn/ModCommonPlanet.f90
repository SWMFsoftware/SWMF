Module ModCommonPlanet
  use ModParameters
  character(6) NamePlanet
  parameter (NamePlanet = 'Saturn')
  
  integer,parameter :: nVar=12
  real,   parameter :: rLowerBoundary  = 6.1668e9
  real,   parameter :: rPlanet = 60268.E5

  REAL XH2(MaxGrid),XH(MaxGrid),XH2O(MaxGrid),XCH4(MaxGrid), &
       XTN(MAXGRID)
  
  REAL CLHpH3p(MaxGrid),CLELHp(MaxGrid),CLELH3p(MaxGrid), &
       CLHpH(MaxGrid),CLELH(MaxGrid)
  REAL FFHpp1(MaxGrid),FFHpp3(MaxGrid),FFHpp4(MaxGrid), &
       FFHpc2(MaxGrid),FFHpc3(MaxGrid), FFHpc9(MaxGrid),&
       FFHpc8(MaxGrid),FFHpr1(MaxGrid),FFH3pc1(MaxGrid),FFH3pc2(MaxGrid), &
       FFH3pc6(MaxGrid),FFH3pc7(MaxGrid),FFH3pr2(MaxGrid)
  
  
  REAL jp1,jp2,jp3,jp4,kc1,kc2,kc3,kc6,kc7,kc8,kr1,kr2
  REAL kc9(MaxGrid)
  
  integer, parameter :: nIon = 3, nSpecies=6
  integer, parameter :: Ion1_ = 1, &    !Saturn:H3+   
       Ion2_ = 2, &    !Saturn:H+
       Neutral1_= 3, & !Saturn:H2
       Neutral2_= 4, & !Saturn:H
       Neutral3_= 5, & !Saturn:H2O
       Neutral4_= 6    !Saturn:CH4
  
  ! named state variables
  integer, parameter :: RhoO_=1, uO_=2,  pO_=3,  To_=4, &
       RhoH_=5, uH_=6,  pH_=7,  Th_=8, &
       RhoE_=9,uE_=10, pE_=11, Te_=12
  

  integer,parameter :: iRho_I(nIon) = (/1,5,9/),&
       iU_I(nIon)   = (/2,6,10/),&
       iP_I(nIon)   = (/3,7,11/),&
       iT_I(nIon)   = (/4,8,12/)
  
  integer,parameter :: iRho(2) = (/1,4/),&
                       iU(2)   = (/2,5/),&
                       iP(2)   = (/3,6/),iTe_=7

  ! Neutral parameters,densities,pressures (needed for ModAurora) 
  integer, parameter :: nNeutral = 5, O_=1, O2_=2, N2_=3, H_=4, He_=5
  real :: NDensity_CI(MaxGrid,nNeutral),NeutralPressure_C(MaxGrid)
  
  ! For Plotting
  integer, parameter :: nPlotVar = 16
  character(len=79), parameter :: NamePlotVar= &
 'r Lat Lon uH3 uH ue lgnH3 lgnH lgne TH3 TH Te MH3 MH Me Ef Pe g'
  


end Module ModCommonPlanet

