Module ModCommonPlanet
  use ModParameters
  character(6) NamePlanet
  parameter (NamePlanet = 'Saturn')
  
  REAL XH2(MaxGrid),XH(MaxGrid),XH2O(MaxGrid),XCH4(MaxGrid), &
       XTN(MAXGRID)
  REAL CFHpH3p(MaxGrid),CFELHp(MaxGrid),CFELH3p(MaxGrid), &
       CFHpH(MaxGrid),CFELH2(MaxGrid),CFH3pHp(MaxGrid),CFHpEL(MaxGrid), &
       CFH3pEL(MaxGrid), &
       CFHpH2(MaxGrid),CFH3pH(MaxGrid),CFH3pH2(MaxGrid),CFELH(MaxGrid) 
     
  REAL CLHpH3p(MaxGrid),CLELHp(MaxGrid),CLELH3p(MaxGrid), &
       CLHpH(MaxGrid),CLELH(MaxGrid)
  REAL FFHpp1(MaxGrid),FFHpp3(MaxGrid),FFHpp4(MaxGrid), &
       FFHpc2(MaxGrid),FFHpc3(MaxGrid), &
       FFHpc8(MaxGrid),FFHpr1(MaxGrid),FFH3pc1(MaxGrid),FFH3pc2(MaxGrid), &
       FFH3pc6(MaxGrid),FFH3pc7(MaxGrid),FFH3pr2(MaxGrid)
  
  REAL  CTHpH,CTHpH2,CTHpH3p,CTHpEL,CTH3pH,CTH3pH2,CTH3pHp, &
       CTH3pEL,CTELH,CTELH2,CTELHp,CTELH3p
  REAL  CMHpH,CMHpH2,CMHpH3p,CMHpEL,CMH3pH,CMH3pH2,CMH3pHp, &
       CMH3pEL,CMELH,CMELH2,CMELHp,CMELH3p
  
  !      REAL jp1,jp2,jp3,jp4,kc1,kc2,kc3,kc6,kc7,kc8,kr1,k2r
  

!!! These are earth parameters that are put here because they are called in
!!! an earth case in solver_rusanov.f90. I will make this more general later

  REAL  CTOXN2,CTOXO2,CTOXO,CTOXHE,CTOXH,CTOXHD,CTOXHL,&
       CTOXEL,CTHN2,CTHO2,CTHHE,CTHO,CTHH,CTHOX,CTHHL,CTHEL,CTHEN2,&
       CTHEO2,CTHEHE,CTHEO,CTHEH,CTHEOX,CTHEHD,CTHEEL,CTELN2,CTELO2,&
       CTELHE,CTELO,CTELOX,CTELHL,CTELHD
  REAL  CMOXN2,CMOXO2,CMOXO,CMOXHE,CMOXH,CMOXHD,CMOXHL,&
       CMOXEL,CMHN2,CMHO2,CMHHE,CMHO,CMHH,CMHOX,CMHHL,CMHEL,CMHEN2,&
       CMHEO2,CMHEHE,CMHEO,CMHEH,CMHEOX,CMHEHD,CMHEEL,CMELN2,CMELO2,&
       CMELHE,CMELO,CMELOX,CMELHL,CMELHD

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
       CFELOX(MaxGrid),CFELHL(MaxGrid),CFELHD(MaxGrid)
end Module ModCommonPlanet

