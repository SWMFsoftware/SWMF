
CALEX This subroutine calculates the collision frequencies and
CALEX then calculates the momentum and energy collision terms
      SUBROUTINE COLLIS(N)
      use ModCommonVariables
C     
C
C
C
      DO 250 I=1,N
      TOX2=TOXYG(I)*TOXYG(I)
      TOX3=TOX2*TOXYG(I)
      TOX4=TOX2*TOX2
      AAH=SQRT(THYD(I)+XTN(I)/16.+1.2E-8*UHYD(I)*UHYD(I))
      BBH=SQRT(XTN(I)+TOXYG(I)/16.+1.2E-8*UOXYG(I)*UOXYG(I))
      CCH=1.-9.149E-4*TOXYG(I)+4.228E-7*TOX2-6.790E-11*TOX3+
     ;4.225E-15*TOX4
      IF (TOXYG(I).GT.1700.) GO TO 201
      FFOXYG=FFOX4(I)*(1.-1.290E-3*TOXYG(I)+6.233E-7*TOX2)
      GO TO 205
201   FFOXYG=FFOX5(I)*(1.-1.410E-3*TOXYG(I)+6.036E-7*TOX2)
205   ADMSO(I)=FFOX1(I)+FFOX2(I)*AAH*DHYD(I)+(FFOX3(I)*BBH+
     ;FFOXYG+FFOX6(I)*CCH)*DOXYG(I)


      ADMSH(I)=FFHYD1(I)*BBH*DOXYG(I)+FFHYD2(I)*AAH*DHYD(I)
      ADMSHE(I)=FFHE1(I)+FFHE2(I)*DHEL(I)
      ADMSE(I)=RTHDEL*ADMSH(I)+RTOXEL*ADMSO(I)+RTHEEL*ADMSHE(I)
C
      TROX=0.5*(XTN(I)+TOXYG(I))
      TRHYD=0.5*(XTN(I)+THYD(I))
      TRHEL=0.5*(XTN(I)+THEL(I))
      T1OX=SQRT(XTN(I)+TOXYG(I)/16.)

CALEX These are (reduced temperatures) * (m1+m2) raised to the 1.5
CALEX as shown on page 86 Nagy. This is for use in collision freqs
CALEX of coulomb collisions below.       
      T1OXH=(TOXYG(I)+16.*THYD(I))**1.5
      T1OXHE=(TOXYG(I)+4.*THEL(I))**1.5
      T1HEH=(THEL(I)+4.*THYD(I))**1.5
      TE32=TELECT(I)**1.5
      DTE32=DELECT(I)/TE32
CALEX calculate collision frequencies. 

CALEX cf of O+ and O
      CFOXO(I)=CLOXO(I)*SQRT(TROX)*(1.-0.064*ALOG10(TROX))**2
CALEX this looks like cf of O+ and H but see comments for cloxh
      CFOXH(I)=CLOXH(I)*T1OX
CALEX cf of O+ and H+. coulomb collisions as eq 4.142 p.95 in Nagy
      CFOXHD(I)=CLOXHD(I)*DHYD(I)/T1OXH
CALEX cf of O+ and He+. coulomb collisions as eq 4.142 p.95 in Nagy      
      CFOXHL(I)=CLOXHL(I)*DHEL(I)/T1OXHE
      
      CFOXEL(I)=CLOXEL(I)*DTE32
CALEX  cf of H+ and O      
      CFHO(I)=CLHO(I)*SQRT(THYD(I))*(1.-0.047*ALOG10(THYD(I)))**2
CALEX cf of H+ and H      
      CFHH(I)=CLHH(I)*SQRT(TRHYD)*(1.-0.083*ALOG10(TRHYD))**2
CALEX cf of H+ and O+. coulomb collisions as eq 4.142 p.95 in Nagy      
      CFHOX(I)=CLHOX(I)*DOXYG(I)/T1OXH
CALEX cf of H+ and He+. coulomb collisions as eq 4.142 p.95 in Nagy      
      CFHHL(I)=CLHHL(I)*DHEL(I)/T1HEH
      
      CFHEL(I)=CLHEL(I)*DTE32
CALEX cf of He+ and HE      
      CFHEHE(I)=CLHEHE(I)*SQRT(TRHEL)*(1.-0.093*ALOG10(TRHEL))**2
CALEX cf of He+ and O+      
      CFHEOX(I)=CLHEOX(I)*DOXYG(I)/T1OXHE
CALEX cf of He+ and H+      
      CFHEHD(I)=CLHEHD(I)*DHYD(I)/T1HEH
CALEX cfheel=7.43E-3/M_e * n_e / T_e^1.5, this seems like it should be
CALEX the cf between He and electrons but the formula does not match 4.144
      CFHEEL(I)=CLHEEL(I)*DTE32
CALEX cf for electrons and O, same as 4.144      
      CFELOX(I)=CLELOX(I)*DOXYG(I)/TE32
CALEX cf for electrons and He      
      CFELHL(I)=CLELHL(I)*DHEL(I)/TE32
CALEX cf for el and H      
      CFELHD(I)=CLELHD(I)*DHYD(I)/TE32
C
CALEX velocity difference needed for source terms
      UHDOX=UHYD(I)-UOXYG(I)
      UHDHL=UHYD(I)-UHEL(I)
      UHDEL=UHYD(I)-UELECT(I)
      UHEOX=UHEL(I)-UOXYG(I)
      UHEEL=UHEL(I)-UELECT(I)
      UOXEL=UOXYG(I)-UELECT(I)

CALEX This calculates collision source terms: 
CALEX fclsn1=n*((u2-u1)*cf12+(u3-u1)*cf13+...)
      FCLSNO(I)=DOXYG(I)*(UHDOX*CFOXHD(I)+UHEOX*CFOXHL(I)-
     $UOXEL*CFOXEL(I)-UOXYG(I)*(CFOXN2(I)+CFOXO2(I)+CFOXO(I)
     $+CFOXH(I)+CFOXHE(I)))
      FCLSNH(I)=DHYD(I)*(-UHDOX*CFHOX(I)-UHDHL*CFHHL(I)-
     $UHDEL*CFHEL(I)-UHYD(I)*(CFHN2(I)+CFHO2(I)+CFHO(I)
     $+CFHH(I)+CFHHE(I)))
      FCLSHE(I)=DHEL(I)*(UHDHL*CFHEHD(I)-UHEOX*CFHEOX(I)
     $-UHEEL*CFHEEL(I)-UHEL(I)*(CFHEN2(I)+CFHEO2(I)+CFHEO(I)+
     $CFHEHE(I)+CFHEHE(I)))

CALEX UHLEL is not defined anywhere in this program!!! based on
CALEX UOXEL I assume it is the relative velocity between helium
CALEX and electrionsand so I define it as such below
      UHLEL=UHEEL
      

CALEX     $UHDEL,CFELHD(I),UELECT(I),CFELN2(I),CFELO2(I),CFELO(I),
CALEX     $CFELHE(I),CFELH(I)',
CALEX     $ DELECT(I),UOXEL,CFELOX(I)
     

      FCLSNE(I)=DELECT(I)*(UOXEL*CFELOX(I)+UHLEL*CFELHL(I)+
     $UHDEL*CFELHD(I)-UELECT(I)*(CFELN2(I)+CFELO2(I)+CFELO(I)+
     $CFELHE(I)+CFELH(I)))
C
CALEX calculate square of velocity differences for use in the
CALEX energy collision term
      UHDOX=UHDOX*UHDOX
      UHDHL=UHDHL*UHDHL
      UHDEL=UHDEL*UHDEL
      UHEOX=UHEOX*UHEOX
      UHEEL=UHEEL*UHEEL
      UOXEL=UOXEL*UOXEL
      UOXN=UOXYG(I)**2
      UHEN=UHEL(I)**2
      UHDN=UHYD(I)**2
      UELN=UELECT(I)**2
CALEX these are temperature differences needed in order to calculate
CALEX the energy collision term 
      THDOX=THYD(I)-TOXYG(I)
      THDHL=THYD(I)-THEL(I)
      THDEL=THYD(I)-TELECT(I)
      THDN=THYD(I)-XTN(I)
      THEOX=THEL(I)-TOXYG(I)
      THEEL=THEL(I)-TELECT(I)
CALEX Can you define a variable THEN?!
      THEN=THEL(I)-XTN(I)
      TOXEL=TOXYG(I)-TELECT(I)
      TOXN=TOXYG(I)-XTN(I)
      TELN=TELECT(I)-XTN(I)

CALEX These are the energy collision terms as seen in eq 4.86 in Nagy
      ECLSNO(I)=DOXYG(I)*(-TOXN*(CTOXN2*CFOXN2(I)+CTOXO2*CFOXO2(I)+
     $CTOXO*CFOXO(I)+CTOXHE*CFOXHE(I)+CTOXH*CFOXH(I))+
     $THDOX*CTOXHD*CFOXHD(I)+THEOX*CTOXHL*CFOXHL(I)-
     $TOXEL*CTOXEL*CFOXEL(I)+UOXN*(CMOXN2*CFOXN2(I)+CMOXO2*CFOXO2(I)+
     $CMOXO*CFOXO(I)+CMOXHE*CFOXHE(I)+CMOXH*CFOXH(I))+
     $UHDOX*CMOXHD*CFOXHD(I)+UHEOX*CMOXHL*CFOXHL(I)+
     $UOXEL*CMOXEL*CFOXEL(I))
      ECLSNH(I)=DHYD(I)*(-THDN*(CTHN2*CFHN2(I)+CTHO2*CFHO2(I)+
     $CTHO*CFHO(I)+CTHHE*CFHHE(I)+CTHH*CFHH(I))-
     $THDOX*CTHOX*CFHOX(I)-THDHL*CTHHL*CFHHL(I)-
     $THDEL*CTHEL*CFHEL(I)+UHDN*(CMHN2*CFHN2(I)+CMHO2*CFHO2(I)+
     $CMHO*CFHO(I)+CMHHE*CFHHE(I)+CMHH*CFHH(I))+
     $UHDOX*CMHOX*CFHOX(I)+UHDHL*CMHHL*CFHHL(I)+
     $UHDEL*CMHEL*CFHEL(I))
      ECLSHE(I)=DHEL(I)*(-THEN*(CTHEN2*CFHEN2(I)+CTHEO2*CFHEO2(I)+
     $CTHEO*CFHEO(I)+CTHEHE*CFHEHE(I)+CTHEH*CFHEH(I))+
     $THDHL*CTHEHD*CFHEHD(I)-THEOX*CTHEOX*CFHEOX(I)-
     $THEEL*CTHEEL*CFHEEL(I)+UHEN*(CMHEN2*CFHEN2(I)+CMHEO2*CFHEO2(I)+
     $CMHEO*CFHEO(I)+CMHEHE*CFHEHE(I)+CMHEH*CFHEH(I))+
     $UHEOX*CMHEOX*CFHEOX(I)+UHDHL*CMHEHD*CFHEHD(I)+
     $UHEEL*CMHEEL*CFHEEL(I))
      ECLSNE(I)=DELECT(I)*(-TELN*(CTELN2*CFELN2(I)+CTELO2*CFELO2(I)+
     $CTELO*CFELO(I)+CTELHE*CFELHE(I)+CTELH*CFELH(I))+
     $TOXEL*CTELOX*CFELOX(I)+THEEL*CTELHL*CFELHL(I)+
     $THDEL*CTELHD*CFELHD(I)+UELN*(CMELN2*CFELN2(I)+CMELO2*CFELO2(I)+
     $CMELO*CFELO(I)+CMELHE*CFELHE(I)+CMELH*CFELH(I))+
     $UOXEL*CMELOX*CFELOX(I)+UHEEL*CMELHL*CFELHL(I)+
     $UHDEL*CMELHD*CFELHD(I))
C
CALEX calculate heat conductivities
      TCONO(I)=HLPO*(DOXYG(I)/DELECT(I))*TOXYG(I)**2.5
      TCONE(I)=HLPE*TELECT(I)**2.5
      TCONH(I)=HLPH*(DHYD(I)/DELECT(I))*THYD(I)**2.5
      TCONH(I)=TCONH(I)/(1.+(0.7692*(CFHN2(I)+CFHO2(I))+
     $1.0962*CFHO(I))/CFHOX(I))
      TCONHE(I)=HLPHE*(DHEL(I)/XMSHE)*THEL(I)
      TCONHE(I)=TCONHE(I)/(0.99*CFHEN2(I)+0.99*CFHEO2(I)+
     $1.02*CFHEO(I)+1.48*CFHEHE(I)+2.22*CFHEH(I)+
     $1.21*CFHEOX(I)+2.23*CFHEHD(I))      
250   CONTINUE
      RETURN
      END
