Module ModCommonPlanet
  use ModParameters
  character(6) NamePlanet
  parameter (NamePlanet = 'Earth ')

  integer,parameter :: nVar=16
  real,   parameter :: rLowerBoundary  = 6.55677E8
  real,   parameter :: rPlanet = 6356.77E5
  
  REAL XN2(MaxGrid),XO2(MaxGrid),XO(MaxGrid),XH(MaxGrid),&
       XHE(MaxGrid),XTN(MaxGrid)
  

  REAL CLOXN2(MaxGrid),CLOXO2(MaxGrid),CLOXO(MaxGrid),&
       CLOXHE(MaxGrid),&
       CLOXH(MaxGrid),CLOXHD(MaxGrid),CLOXHL(MaxGrid),CLOXEL(MaxGrid),&
       CLHN2(MaxGrid),&
       CLHO2(MaxGrid),CLHHE(MaxGrid),CLHO(MaxGrid),CLHH(MaxGrid),&
       CLHOX(MaxGrid),&
       CLHHL(MaxGrid),CLHEL(MaxGrid),CLHEN2(MaxGrid),CLHEO2(MaxGrid),&
       CLHEHE(MaxGrid),&
       CLHEO(MaxGrid),CLHEH(MaxGrid),CLHEOX(MaxGrid),CLHEHD(MaxGrid),&
       CLHEEL(MaxGrid),&
       CLELN2(MaxGrid),CLELO2(MaxGrid),CLELHE(MaxGrid),CLELO(MaxGrid),&
       CLELH(MaxGrid),&
       CLELOX(MaxGrid),CLELHL(MaxGrid),CLELHD(MaxGrid)
  REAL FFOX1(MaxGrid),FFOX2(MaxGrid),FFOX3(MaxGrid),&
       FFOX4(MaxGrid),FFOX5(MaxGrid),&
       FFOX6(MaxGrid),FFHYD1(MaxGrid),FFHYD2(MaxGrid),FFHE1(MaxGrid),&
       FFHE2(MaxGrid)
  
  REAL  CTOXN2,CTOXO2,CTOXO,CTOXHE,CTOXH,CTOXHD,CTOXHL,&
       CTOXEL,CTHN2,CTHO2,CTHHE,CTHO,CTHH,CTHOX,CTHHL,CTHEL,CTHEN2,&
       CTHEO2,CTHEHE,CTHEO,CTHEH,CTHEOX,CTHEHD,CTHEEL,CTELN2,CTELO2,&
       CTELHE,CTELO,CTELH,CTELOX,CTELHL,CTELHD
  REAL  CMOXN2,CMOXO2,CMOXO,CMOXHE,CMOXH,CMOXHD,CMOXHL,&
       CMOXEL,CMHN2,CMHO2,CMHHE,CMHO,CMHH,CMHOX,CMHHL,CMHEL,CMHEN2,&
       CMHEO2,CMHEHE,CMHEO,CMHEH,CMHEOX,CMHEHD,CMHEEL,CMELN2,CMELO2,&
       CMELHE,CMELO,CMELH,CMELOX,CMELHL,CMELHD
  
  ! Neutral parameters,densities,pressures
  integer, parameter :: nNeutral = 5, O_=1, O2_=2, N2_=3, H_=4, He_=5
  real :: NDensity_CI(MaxGrid,nNeutral),NeutralPressure_C(MaxGrid)

  !collision terms
  
  integer, parameter :: nIon = 4, nSpecies=9
  integer, parameter :: Ion1_ = 1, &    !Earth:O+ ,Saturn:H3+
                        Ion2_ = 2, &    !Earth:H+ ,Saturn:H+
                        Ion3_ = 3, &    !Earth:He+,Saturn:none
                        Ion4_ = 4, &    !Earth:e  ,Saturn:e
                        Neutral1_= 5, & !Earth:O, Saturn:H2
                        Neutral2_= 6, & !Earth:H, Saturn:H
                        Neutral3_= 7, & !Earth:O2, Saturn:H2O
                        Neutral4_= 8, & !Earth:N2, Saturn:CH4
                        Neutral5_= 9    !Earth:He, Saturn:none

  ! named state variables
  integer, parameter :: RhoO_=1, uO_=2,  pO_=3,  To_=4, &
                        RhoH_=5, uH_=6,  pH_=7,  Th_=8, &
                        RhoHe_=9,uHe_=10,pHe_=11,The_=12, &
                        RhoE_=13,uE_=14, pE_=15, Te_=16

  integer,parameter :: iRho_I(4) = (/1,5,9,13/),&
                       iU_I(4)   = (/2,6,10,14/),&
                       iP_I(4)   = (/3,7,11,15/),&
                       iT_I(4)   = (/4,8,12,16/)

  ! for reduced state used in ModPwImplicit
  integer,parameter :: iRho(3) = (/1,4,7/),&
                       iU(3)   = (/2,5,8/),&
                       iP(3)   = (/3,6,9/),iTe_=10
                       
  
  integer, parameter :: nPlotVar = 20
  character(len=79), parameter :: NamePlotVar= &
 'r Lat Lon uO uH uHe ue lgnO lgnH lgnHe lgne TO TH THe Te MO MH MHe Me Ef Pe g'

end Module ModCommonPlanet
    
