Module ModCommonPlanet
  use ModParameters
  character(6) NamePlanet
  parameter (NamePlanet = 'Earth ')
  
  REAL XN2(MaxGrid),XO2(MaxGrid),XO(MaxGrid),XH(MaxGrid),&
       XHE(MaxGrid),XTN(MaxGrid)
  
  REAL CFOXN2(MaxGrid),CFOXO2(MaxGrid),CFOXO(MaxGrid),&
       CFOXHE(MaxGrid),&
       CFOXH(MaxGrid),CFOXHD(MaxGrid),CFOXHL(MaxGrid),CFOXEL(MaxGrid),&
       CFHN2(MaxGrid),&
       CFHO2(MaxGrid),CFHHE(MaxGrid),CFHO(MaxGrid),CFHH(MaxGrid),&
       CFHOX(MaxGrid),&
       CFHHL(MaxGrid),CFHEL(MaxGrid),CFHEN2(MaxGrid),CFHEO2(MaxGrid),&
       CFHEHE(MaxGrid),&
       CFHEO(MaxGrid),CFHEH(MaxGrid),CFHEOX(MaxGrid),CFHEHD(MaxGrid),&
       CFHEEL(MaxGrid),&
       CFELN2(MaxGrid),CFELO2(MaxGrid),CFELHE(MaxGrid),CFELO(MaxGrid),&
       CFELH(MaxGrid),&
       CFELOX(MaxGrid),CFELHL(MaxGrid),CFELHD(MaxGrid)
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


!!! These are saturn parameters that are put here because they are called in
!!! a Saturn case in solver_rusanov.f90. I will make this more general later
  REAL CFHpH(MaxGrid), CFHPH2(MaxGrid), CTHpH, CTHpH2, &
       CFHpH3(MaxGrid), CTHpH3p, &
       CFHpEL(MaxGrid), CTHpEL, &
       CFH3pH(MaxGrid), CFH3pH2(MaxGrid), CTH3pH2, &
       CFH3pHp(MaxGrid), CTH3pHp, &
       CFH3pEL(MaxGrid), CTH3pEL, CFHpH3p(MaxGrid), CTH3pH
 
      
    end Module ModCommonPlanet
    
