
CALEX This subroutine calculates the collision frequencies and
CALEX then calculates the momentum and energy collision terms
      SUBROUTINE COLLIS(N)
      use ModCommonVariables
C     
C
C
C
      DO 250 I=1,N
         
C**********************************************************************
C Determine the mass sources from chemistry
C**********************************************************************

      ADMSH(I)=(FFHpp1(I)+FFHpp3(I)+FFHpp4(I)+FFHpc2(I)*DHYD(I)/XMSH
     ;+FFHpc3(I)*DHYD(I)/XMSH+FFHpc8(I)*DHYD(I)/XMSH
     ;+FFHpr1(I)*DHYD(I)*DELECT(I)/XMSH/XMSE)*XMSH

      ADMSO(I)=(FFH3pc1(I)+FFH3pc2(I)*DHYD(I)/XMSH+FFH3pc6(I)*DOXYG(I)/XMSO
     ;+FFH3pc7(I)*DOXYG(I)/XMSO
     ;+FFH3pr2(I)*(TELECT(I)**(-.65))*DOXYG(I)*DELECT(I)/XMSO/XMSE)*XMSO
      
      ADMSHE(I)=0.
      ADMSE(I)=RTHDEL*ADMSH(I)+RTOXEL*ADMSO(I)




C**********************************************************************
C Calculate collision frequencies. 
C**********************************************************************

      TROX=0.5*(XTN(I)+TOXYG(I))
      TRHYD=0.5*(XTN(I)+THYD(I))
      TRHEL=0.5*(XTN(I)+THEL(I))
      T1OX=SQRT(XTN(I)+TOXYG(I)/16.)

C These are (reduced temperatures) * (m1+m2) raised to the 1.5
C as shown on page 86 Nagy. This is for use in collision freqs
C of coulomb collisions below.       
      T1HpH3p=(TOXYG(I)+3.*THYD(I))**1.5
      T1OXHE=(TOXYG(I)+4.*THEL(I))**1.5
      T1HEH=(THEL(I)+4.*THYD(I))**1.5
c      write(*,*) TELECT(I)
      TE32=TELECT(I)**1.5
      DTE32=DELECT(I)/TE32


C H+ and H3+         
         CFHpH3p(I)=CLHpH3p(I)*DOXYG(I)/T1HpH3p
C electron H+ and electron H3+
         CFELHp(I)  = CLELHp(I)*DHYD(I)/TE32
         CFELH3p(I) = CLELH3p(I)*DOXYG(I)/TE32
                
C  ion neutrals
         CFHpH(I)=CLHpH(I)*SQRT(TRHYD)*(1.-.083*ALOG10(TRHYD))**2.

C electron H, e H2 done in collis
         CFELH(I)= CLELH(I)*(1.-1.35E-4*TELECT(I))*SQRT(TELECT(I))
         CFELH2(I)=getcfeh2(TELECT(I),XH2(I),XMSE,UELECT(I))

!------------------------------------------------------------------------------
! Turn off collision frequency terms to find the instability in mom sources

!         CFHpH(I)  =0.
!         CFHpH2(I) =0.
!         CFHpH3p(I)=0.
!         CFELHp(I) =0.
!         CFELH3p(I)=0.
!         CFH3pH(I) =0.
!         CFH3pH2(I)=0.

         
!------------------------------------------------------------------------------

C Now get the inverse collision freq
         CFH3pHp(I)=DHYD(I)/DOXYG(I)*CFHpH3p(I)
         CFHpEL(I) = DELECT(I)/DHYD(I)*CFELHp(I)
         CFH3pEL(I)=DELECT(I)/DOXYG(I)*CFELH3p(I)




C**********************************************************************
C Determin the momentum source terms
C**********************************************************************

C Velocity difference needed for source terms
      UHDOX=UHYD(I)-UOXYG(I)
      UHDHL=UHYD(I)-UHEL(I)
      UHDEL=UHYD(I)-UELECT(I)
c      UHEOX=UHEL(I)-UOXYG(I)
c      UHEEL=UHEL(I)-UELECT(I)
      UOXEL=UOXYG(I)-UELECT(I)

C This calculates collision source terms: 
C fclsn1=n*((u2-u1)*cf12+(u3-u1)*cf13+...)
      FCLSNO(I)=DOXYG(I)*(UHDOX*CFH3pHp(I)-
     $UOXEL*CFH3pEL(I)-UOXYG(I)*(CFH3pH(I)+CFH3pH2(I)))
      
      FCLSNH(I)=DHYD(I)*(-UHDOX*CFHpH3p(I)-
     $UHDEL*CFHpEL(I)-UHYD(I)*(CFHpH(I)+CFHpH2(I)))

! if implicit then take out the contribution to the term that is being made
! implicit
!      If (IsImplicit) then
!         FCLSNH(I)=FCLSNH(I)+uHyd(I)*CFHpH(I)*dHyd(I)
!         
!      endif

      FCLSHE(I)=0.
 
      FCLSNE(I)=DELECT(I)*(UOXEL*CFELH3p(I)+
     $UHDEL*CFELHp(I)-UELECT(I)*(CFELH(I)+CFELH2(I)))



C**********************************************************************
C Determine the energy source terms
C**********************************************************************
C
C Calculate square of velocity differences for use in the
C  energy collision term
      UHDOX=UHDOX*UHDOX
C      UHDHL=UHDHL*UHDHL
      UHDEL=UHDEL*UHDEL
C      UHEOX=UHEOX*UHEOX
C      UHEEL=UHEEL*UHEEL
      UOXEL=UOXEL*UOXEL
      UOXN=UOXYG(I)**2
C      UHEN=UHEL(I)**2
      UHDN=UHYD(I)**2
      UELN=UELECT(I)**2
CALEX these are temperature differences needed in order to calculate
CALEX the energy collision term 
      THDOX=THYD(I)-TOXYG(I)
C      THDHL=THYD(I)-THEL(I)
      THDEL=THYD(I)-TELECT(I)
      THDN=THYD(I)-XTN(I)
      THEOX=THEL(I)-TOXYG(I)
C      THEEL=THEL(I)-TELECT(I)
CALEX Can you define a variable THEN?!
CALEX      THEN=THEL(I)-XTN(I)
      TOXEL=TOXYG(I)-TELECT(I)
      TOXN=TOXYG(I)-XTN(I)
      TELN=TELECT(I)-XTN(I)
      
!      if (IsImplicit) then
!         TOXN=-XTN(I)
!         THDN=-XTN(I)
!      endif

C These are the energy collision terms as seen in eq 4.86 in Nagy
      ECLSNO(I)=DOXYG(I)*(-TOXN*(CTH3pH*CFH3pH(I)+CTH3pH2*CFH3pH2(I))+
     $THDOX*CTH3pHp*CFH3pHp(I)-
     $TOXEL*CTH3pEL*CFH3pEL(I)+UOXN*(CMH3pH*CFH3pH(I)+CMH3pH2*CFH3pH2(I)
     $)+UHDOX*CMH3pH*CFH3pH(I)+
     $UOXEL*CMH3pEL*CFH3pEL(I))
      ECLSNH(I)=DHYD(I)*(-THDN*(CTHpH*CFHpH(I)+CTHpH2*CFHpH2(I))-
     $THDOX*CTHpH3p*CFHpH3p(I)-
     $THDEL*CTHpEL*CFHpEL(I)+UHDN*(CMHpH*CFHpH(I)+CMHpH2*CFHpH2(I))+
     $UHDOX*CMHpH3p*CFHpH3p(I)+
     $UHDEL*CMHpEL*CFHpEL(I))
      ECLSHE(I)=0.
C      ECLSHE(I)=DHEL(I)*(-THEN*(CTHEN2*CFHEN2(I)+CTHEO2*CFHEO2(I)+
C     $CTHEO*CFHEO(I)+CTHEHE*CFHEHE(I)+CTHEH*CFHEH(I))+
C     $THDHL*CTHEHD*CFHEHD(I)-THEOX*CTHEOX*CFHEOX(I)-
C     $THEEL*CTHEEL*CFHEEL(I)+UHEN*(CMHEN2*CFHEN2(I)+CMHEO2*CFHEO2(I)+
C     $CMHEO*CFHEO(I)+CMHEHE*CFHEHE(I)+CMHEH*CFHEH(I))+
C     $UHEOX*CMHEOX*CFHEOX(I)+UHDHL*CMHEHD*CFHEHD(I)+
C     $UHEEL*CMHEEL*CFHEEL(I))
      ECLSNE(I)=DELECT(I)*(-TELN*(CTELH*CFELH(I)+CTELH2*CFELH2(I))+
     $TOXEL*CTELH3p*CFELH3p(I)+
     $THDEL*CTELHp*CFELHp(I)+UELN*(CMELH*CFELH(I)+CMELH2*CFELH2(I))+
     $UOXEL*CMELH3p*CFELH3p(I)+UHDEL*CMELHp*CFELHp(I))
C

C No sources 
c      if (I .gt. 465) then
!         ADMSH(I)=0.
!         ADMSO(I)=0.
!         ADMSHE(I)=0.
!         ADMSE(I)=0.
!c
!         FCLSNO(I)=0.
!         FCLSNH(I)=0.
!         FCLSHE(I)=0.
!         FCLSNE(I)=0.
!c
!         ECLSNO(I)=0.
!         ECLSNH(I)=0.
!         ECLSHE(I)=0.
!         ECLSNE(I)=0.
c      endif

C**********************************************************************
C Calculate heat conductivities
C**********************************************************************

      TCONO(I)=HLPO*(DOXYG(NDIM)/DELECT(NDIM))*TOXYG(I)**2.5
      TCONE(I)=HLPE*TELECT(I)**2.5
      TCONH(I)=HLPH*(DHYD(NDIM)/DELECT(NDIM))*THYD(I)**2.5
      TCONHE(I)=0.
C      TCONHE(I)=HLPHE*THEL(I)**2.5
            
250   CONTINUE
CALEX      WRITE(*,*) DOXYG(NDIM),DHYD(NDIM)
CALEX      CLOSE(unit=37)
      RETURN
      END
