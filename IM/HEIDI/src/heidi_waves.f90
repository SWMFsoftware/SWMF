! File name: heidi_waves_010.f90
!
! Contains: wave-particle interaction routines for HEIDI
!	WAPARA
!	ANISCH
!
! Last Modified: March 2006, Mike Liemohn
!
! **********************************************************************
!                              WAPARA 
!       Routine calculates normalized waves pa diffusion coeff
!	NOTE: This subroutine has not been upgraded for mram03.
!***********************************************************************
SUBROUTINE WAPARA 

  use ModHeidiSize
  use ModHeidiMain
  use ModHeidiWaves
  use ModIoUnit, ONLY : io_unit_new
  
  IMPLICIT NONE

  CHARACTER*80 HEADER
  CHARACTER*5 ST1
  CHARACTER*2 ST2
  REAL :: keVerg
  INTEGER :: KN,L,KK
  DATA keVerg/1.602E-9/

  integer :: iUnit16! = 22
  integer :: iUnit18! = 23
  integer :: iUnit20! = 24

  IF (T.EQ.0.0) 	   ST1='fc_12'	! initial
  IF (T.EQ.1800.)    ST1='fc_12'	! in 0.5 hour
  IF (T/3600.EQ.1.)  ST1='fc_14'	! in 1 hours
  IF (T/3600.EQ.2.)  ST1='fc_16'	! in 2 hours
  IF (T/3600.EQ.3.)  ST1='fc_18'	! in 3 hours
  IF (T/3600.GE.4.)  ST1='fc_20'	! in 4 hours
  IF(S.EQ.2)ST2='_h'
  IF(S.EQ.3)ST2='he'
  IF(S.EQ.4)ST2='_o'   

  !	OPEN(UNIT=22,NAME='wmd'//ST1//ST2//'.in',STATUS='OLD') ! multiion
  !        WRITE(22,7) T/3600
  
  iUnit16 =io_unit_new()
  
  OPEN(UNIT=iUnit16,FILE='wmdfc_16'//ST2//'.in',STATUS='OLD')
  READ(iUnit16,20) HEADER
  
  iUnit18 =io_unit_new()
  
  OPEN(UNIT=iUnit18,FILE='wmdfc_18'//ST2//'.in',STATUS='OLD')
  READ(iUnit18,20) HEADER
  
  iUnit20 =io_unit_new()

  OPEN(UNIT=iUnit20,FILE='wmdfc_20'//ST2//'.in',STATUS='OLD')
  READ(iUnit20,20) HEADER
7 FORMAT(2X,19HNormalized Daa, T =,F8.0)

  !.......make a look up table for normalized diff coeff
  !        degrad=pi/180.
  !        etn0min=0.1
  !        E1=ALOG10(etn0min)
  !        etn0max=1e4
  !        E2=ALOG10(etn0max)
  !        DEW=(E2-E1)/(ENG-1)
  !        DO KN=1,ENG
  !         EX=E1+DEW*(KN-1)
  !         ENOR(KN)=10**EX
  !        ENDDO
  DO KN=1,ENG
     !         write(22,17)kn,enor(kn)
     read(iUnit16,17)kk,enor(kn,1) 
     read(iUnit18,17)kk,enor(kn,2)
     read(iUnit20,17)kk,enor(kn,3)
     !	 CALL WPIDC(ENOR(KN),S,0.15,0.1,0.2,500.,4.)
     !         DO L=1,NPA-1
     !          NDAAJ(KN,L)=NDAA3(L)          ! normalized diff coeff Daa
     !         ENDDO
     !         write(22,27)(ndaaj(kn,l),l=1,LO)
     read(iUnit16,27) (ndaaj(kn,l,1),l=1,LO)
     read(iUnit18,27) (ndaaj(kn,l,2),l=1,LO)
     read(iUnit20,27) (ndaaj(kn,l,3),l=1,LO)
  enddo

17 FORMAT(I6,E13.4)
20 FORMAT(A80)
27 FORMAT(8(1PE12.3))
  CLOSE(iUnit16)
  CLOSE(iUnit18)
  CLOSE(iUnit20)

  RETURN
END SUBROUTINE WAPARA
!
! End of subroutine WAPARA
!

  !**********************************************************************
  !			 	ANISCH
  !     Routine calculates the temperature anisotropy
  !            defining the strong growth region
  !     NOTE: This subroutine has not been updated for mram03!!
  !**********************************************************************
SUBROUTINE ANISCH

  use ModHeidiSize
  use ModHeidiCurrents
  use ModHeidiIO
  use ModHeidiWaves
  use ModHeidiMain

  implicit none

  REAL :: DWAVE(NPA),CMRA(SLEN),AVDAA(NPA),TAVDAA(NPA),cubL,   &
       DAA(NE,NPA,Slen),DN(2,25),EPO(2,25),AII(2,25),XFR(NR,NT),   &
       Bline,RLD,PHID,CAML,degrad,rlmas,pmas,cv,gausgam,   &
       WGG,GRP,VGP,rat,bf,ZE,QFAC,omega,y,res1,res2,fac6,band,esu,   &
       bwave2,fnorm,ETN,tang,delf,pr1,pr2,prm1,prm2,keVerg,   &
       tfac11,tfac21,daabar,delml,muboun,taudaa,afir,cjper,pero,   &
       paro,wuc,wlc,delw,asec,erf,eden1!,funt
  INTEGER :: MINP(NT),MAXP(NT),kprint,I,J,K,L,j_2,IML,ical,KN,iwn,   &
       IER, kres
  CHARACTER*5 ST1
  CHARACTER*2 ST2
  CHARACTER*5 ST4
  CHARACTER*80 HEADER
  external :: erf!,funt

  integer :: iUnitOWans=11
  integer :: iUnitWgr=9
  integer :: iUnitWans=8
  integer :: iUnitWdc=22
  integer :: iUnitJper=3

  !real    :: funt(nPa,nR,nT),funi(nPa,nR,nT)


  DATA pmas/1.672E-24/, keVerg/1.602E-9/,    &
       cv/3E10/, esu/4.803E-10/, gausgam/1.E-5/

  !call get_IntegralH(funt)
  !call get_IntegralI(funi)


  kprint=	25		! 27 (E=51), 30 (E=100), 33 (E=200)
  IF(S.EQ.2)ST2='_h'
  IF(S.EQ.3)ST2='he'
  IF(S.EQ.4)ST2='_o'

  IF (T.EQ.0.) ST1='wpap0'
  IF (T/3600.GE.1.AND.T/3600.LE.24) ST1='wpap1'	 
  IF (T/3600.GT.24.AND.T/3600.LE.48.) ST1='wpap2'
  IF (T/3600.GT.48.AND.T/3600.LE.60.) ST1='wpap3'	  
  IF (T/3600.GT.60.AND.T/3600.LE.72.) ST1='wpap4'	  
  IF (T/3600.GT.72.AND.T/3600.LE.84.) ST1='wpap5'
  IF (T/3600.GT.84.AND.T/3600.LE.90.) ST1='wpap6'
  IF (T/3600.GT.90.AND.T/3600.LE.92.) ST1='wpap7'
  IF (T/3600.GT.92.AND.T/3600.LE.94.) ST1='wpap8'	  
  IF (T/3600.GT.94.AND.T/3600.LE.96.) ST1='wpa10'
  IF (T/3600.EQ.97) ST1='wpa11'
  IF (T/3600.EQ.98) ST1='wpa12'
  IF (T/3600.EQ.99) ST1='wpa13'
  IF (T/3600.EQ.100) ST1='wpa14'
  IF (T/3600.EQ.101) ST1='wpa15'
  IF (T/3600.EQ.102) ST1='wpa16'
  IF (T/3600.EQ.103) ST1='wpa17'
  IF (T/3600.EQ.104) ST1='wpa18'
  IF (T/3600.EQ.105) ST1='wpa19'
  IF (T/3600.EQ.106) ST1='wpa20'
  IF (T/3600.EQ.107) ST1='wpa21'
  IF (T/3600.EQ.108) ST1='wpa22'
  IF (T/3600.GT.108) ST1='wpa23'
  !	IF (T/3600.GT.96.AND.T/3600.LE.98.) ST1='wpa10'	  
  !        IF (T/3600.GT.98.AND.T/3600.LE.100.) ST1='wpa11'	  
  !        IF (T/3600.GT.100.AND.T/3600.LE.102.) ST1='wpa12'	  
  !        IF (T/3600.GT.102.AND.T/3600.LE.104.) ST1='wpa13'	  
  !        IF (T/3600.GT.104.AND.T/3600.LE.108.) ST1='wpa14'	  

  !        IF (AMOD(T,7200).NE.0.) GO TO	100

  IF (T/3600.LE.96.) ST4='rna10'
  IF (T/3600.EQ.97.) ST4='rna11'
  IF (T/3600.EQ.98.) ST4='rna12'
  IF (T/3600.EQ.99.) ST4='rna13'
  IF (T/3600.EQ.100) ST4='rna14'
  IF (T/3600.EQ.101) ST4='rna15'
  IF (T/3600.EQ.102) ST4='rna16'
  IF (T/3600.EQ.103) ST4='rna17'
  IF (T/3600.EQ.104) ST4='rna18'
  IF (T/3600.EQ.105) ST4='rna19'
  IF (T/3600.EQ.106) ST4='rna20'
  IF (T/3600.GE.107) ST4='rna21'
  !        IF (T/3600.EQ.108) ST4='rna22'


  !	OPEN(UNIT=10,FILE=ST4//'he.ans',STATUS='OLD')
  !	READ(10,20) HEADER
  OPEN(UNIT=iUnitOWans,FILE=ST4//'_o.wans',STATUS='OLD')
  READ(iUnitOWans,20) HEADER
  OPEN(UNIT=iUnitWgr,FILE=ST1//ST2//'.wgr',STATUS='UNKNOWN')
  WRITE(iUnitWgr,55) T/3600,KP
55 FORMAT(2HT=,F6.1,2X,3HKp=,F3.1,X,   &
       21HL  PHI  GAIN[dB]  XFR)
  OPEN(UNIT=iUnitWans,FILE=ST1//ST2//'.wans',STATUS='UNKNOWN')
  WRITE(iUnitWans,56) T/3600,KP     
56 FORMAT(2HT=,F6.1,2X,3HKp=,F3.1,X,   &
       'L   PHI     ANIS   EDEN[keV/cm3]   RNHT[1/cm3]',3X,   &
       'PPER[keV/cm3] PPAR[keV/cm3]')
  OPEN(UNIT=iUnitWdc,FILE=ST1//ST2//'.wdc',STATUS='UNKNOWN')
  WRITE(iUnitWdc,7) T/3600,KP     
7 FORMAT(2X,3HT =,F6.1,2X,4HKp =,F6.2) 
  !	OPEN(UNIT=9,FILE='rcc10.in',STATUS='OLD')
  !	READ(9,20) HEADER

  rlmas=MAS(S)*1E3	! mass of reson part in [g]
  degrad=pi/180.
  Bline=0.
  CMRA(1)=0.
  DO IML=2,14
     !	DO IML=2,ISO
     CMRA(IML)=AMLA(IML)*DEGRAD
     CAML=SQRT(1+3*SIN(CMRA(IML))*SIN(CMRA(IML)))*COS(CMRA(IML))
     Bline=Bline+CAML*(CMRA(IML)-CMRA(IML-1))
  ENDDO

  !.......make a look up table for normalized diff coeff
  cubL=0.32*esu/(cv*MAS(2)*1E3*2*pi)
  !	 wmL=0.058
  !	 delwL=0.0145
  !	 wucL=0.06
  !	 wlcL=0.04
  !	 delfL=delwL*cubL/(5.)**3
  !	IF (T/3600.EQ.2.) THEN 
  !	 wm=0.16
  !	 delw=0.04
  !	 wuc=0.18
  !	 wlc=0.14
  !	 delf=delw*cubL/(4.25)**3 
  !	ENDIF
  !	IF (T/3600.EQ.3.) THEN 
  !	 wm=0.18
  !	 delw=0.045
  !	 wuc=0.22
  !	 wlc=0.15
  !	 delf=delw*cubL/(3.75)**3
  !	ENDIF
  !        IF (T/3600.GE.4.) THEN  
  !	 wm=0.20
  !	 delw=0.05
  !	 wuc=0.22
  !	 wlc=0.16
  !	 delf=delw*cubL/(3.5)**3
  !	ENDIF

  CALL PRESSURES
  DWAVE(1)=0.
  DO I=2,IO
     DO J=1,JO
        !	  READ(10,551) RLD,PHID,AII(1,j),EPO(1,j),DN(1,j)  ! moder He+
        READ(iUnitOWans,551) RLD,PHID,AII(2,j),eden1,DN(2,j),pero,paro
        EPO(2,j)=paro/DN(2,j)	!O+
     ENDDO
     DO J=1,JO
        if (JO.eq.24) j_2=j
        if (JO.eq.48) j_2=j/2+1 
        DO K=2,KO
           DO L=1,LO
              ATAW(I,J,K,L)=0.
              GTAW(I,J,K,L)=0.
              DO IML=1,ISO
                 DAA(K,L,IML)=0.
              ENDDO
           ENDDO
        ENDDO
        !.......write anisotropy, E parallel, RC dens 
        !	WRITE(iUnitWans,551)LZ(I),PHI(J),ANIS(I,J,S),EPAR(I,J,S),RNHT(I,J,S)
        WRITE(iUnitWans,551)LZ(I),PHI(J),ANIS(I,J,S),EDEN(I,J,S),RNHT(I,J,S),   &
             PPER(I,J,S),PPAR(I,J,S)

        !	go to 15				! no waves

        !	READ(9,10) rli,ppi,wcd(i,j),wm
        !          wuc=wm+0.02
        !          wlc=wm-0.03
        !          delw=wm/4.
        !          delf=delw*cubL/(3.75)**3
        IF (S.EQ.2) &
             CALL WGRDAR(LZ(I),XNE(I,J),RNHT(I,J,S),EPAR(I,J,S),ANIS(I,J,S),&
             !      	DN(1,j_2),EPO(1,j_2),AII(1,j_2),   &
             2e-15,0.06,0.02,	   &			! He+
             DN(2,j_2),EPO(2,j_2),AII(2,j_2),WGG,xfr(i,j),GRP,VGP)	! O+
        !     	WCD(I,J)=WGG*LZ(I)*Bline*4*RE*100.		! Bline in [cm]
        WCD(I,J)=WGG*LZ(I)*Bline*2*RE*100.*8.6859	! +-5 deg, in dB
     enddo
  enddo

  !	go to 100					! no waves

  DO J=1,JO
     ical=1
     DO I=2,IO-1
        rat=XNE(i,j)/XNE(i+1,j)
        if (rat.ge.10.) then
           !	   print*,'pp grad, l=',lz(i)
           if (ical.eq.1) minp(j)=i
           !	   print*,'minp(j)=',minp(j)
           maxp(j)=i
           wcd(i,j)=wcd(i,j)*2.
           ical=ical+1
           !	  print*,'maxp(j)=',maxp(j)
        endif
     ENDDO
     wcd(minp(j)-1,j)=wcd(minp(j)-1,j)*1.5
     !	 wcd(minp(j)-3,j)=wcd(minp(j)-3,j)*0.75
     !	 wcd(minp(j)-4,j)=wcd(minp(j)-4,j)*0.75
     !	 wcd(minp(j)-5,j)=wcd(minp(j)-5,j)*0.75
     !	 wcd(minp(j)-6,j)=wcd(minp(j)-6,j)*0.75
     wcd(maxp(j)+1,j)=wcd(maxp(j)+1,j)*1.5
  ENDDO

  DO I=2,IO
     DO J=1,JO
        IF(WCD(I,J).GT.70.0) THEN
           BWAVE(I,J)=10.			! Wave B field [nT]
        ELSE
           !	   BWAVE(I,J)=EXP((2*WCD(I,J)-40.)/8.6859)
           BWAVE(I,J)=10*10**((WCD(I,J)-70.)/20)
        ENDIF
        WRITE(iUnitWgr,10)LZ(I),PHI(J),WCD(I,J),xfr(i,j),BWAVE(I,J)
        !	BWAVE(I,J)=1.0		 		! Wave B field [nT^2/Hz]
     ENDDO
  ENDDO

  !	go to 100					! no waves

  DO I=2,IO
     DO J=1,JO
        if (xfr(i,j).le.0.16) iwn=1
        if (xfr(i,j).eq.0.18) iwn=2
        if (xfr(i,j).ge.0.20) iwn=3
        wuc=xfr(i,j)+0.02
        wlc=xfr(i,j)-0.03
        delw=xfr(i,j)/4.
        delf=delw*cubL/(3.75)**3
        !	BWAVE(I,J)=1.0
        IF(WCD(I,J).GT.30.0) THEN	! W gr & coll damping check
           DO IML=1,14
              !	 DO 2 IML=1,ISO
              bf=BE(I,IML)
              ZE=XNE(I,J)*BE(I,IML)/BE(I,1)
              QFAC=keVerg*8*PI*ZE/bf/bf
              omega=esu*bf/(rlmas*cv)
              y=(wuc-xfr(i,j))/delw
              res1=erf(y,IER)
              y=(xfr(i,j)-wlc)/delw
              res2=erf(y,IER)
              fac6=res1+res2
              band=(SQRT(pi)/2.)*(2*pi*delf)*fac6
              !          bwave2=bwave(i,j)*band*gausgam**2
              !	  bwave2=0.5*0.5*gausgam**2			! nT
              bwave2=bwave(i,j)*bwave(i,j)*gausgam**2	! nT
              fnorm=omega*bwave2/bf/bf
              DO K=2,KO
                 ETN=EKEV(K)*QFAC
                 DO KN=1,ENG-1
                    IF(ENOR(KN,IWN).LE.ETN.AND.ENOR(KN+1,IWN).GT.ETN) THEN
                       KRES=KN
                       IF(ABS(ENOR(KN,IWN)-ETN).GT.ABS(ENOR(KN+1,IWN)-ETN))    &
                            KRES=KN+1
                       DO L=1,LO-1
                          DAA(K,L,IML)=NDAAJ(KRES,L,IWN)*fnorm	! denorm [1/s]
                       enddo
                    ENDIF
                 enddo
              enddo

           enddo
           !.....zero out bounce-averaged quantities and bounce-average
           tang=asin(1.0)
           do k=2,KO
              if (lz(i).eq.(4).and.k.eq.kprint)then
                 write(iUnitWdc,550)
                 write(iUnitWdc,*)'L=',LZ(I),' MLT=',MLT(J),' Y=',xfr(i,j),   &
                      ' fnorm=',fnorm,' bwave2=',bwave2,' band=',band,   &
                      ' delf=',delf
              endif
              do l=2,LO-1
                 avdaa(L)=0.0
                 tavdaa(L)=0.0
                 do 524 iml=1,13
                    !      do 524 iml=1,ISO-1
                    pr1=ZRpabn(i,L,iml)
                    pr2=ZRpabn(i,L,iml+1)
                    if(pr1.eq.-1.0)pr1=asin(1.0)
                    if(pr2.eq.-1.0)pr2=asin(1.0)
                    if((pr1.ge.tang).or.(pr2.ge.tang))go to 524
                    prm1=amla(iml)*degrad
                    prm2=amla(iml+1)*degrad
                    tfac11=cos(prm1)**7.*cos(pr1)
                    tfac21=cos(prm2)**7.*cos(pr2)
                    daabar=(tfac11*daa(k,L,iml)+tfac21*daa(k,L,iml+1))/2.
                    delml=prm2-prm1
                    avdaa(L)=avdaa(L)+daabar*delml
524              continue
                 !.....write out bounce-averaged results	
                 MUBOUN=MU(L)+0.5*WMU(L)
                 DWAVE(L)=AVDAA(L)*DT*(1.-MUBOUN*MUBOUN)/MUBOUN 
                 if (lz(i).eq.(4).and.k.eq.kprint)then
                    !taudaa=dwave(l)/DT/MUBOUN/FUNT(MUBOUN)
                    taudaa=dwave(l)/DT/MUBOUN/(0.5*(FUNT(l+1,i,j)-FUNT(l,i,j)))
                    if(taudaa.lt.1e-20) taudaa=1e-20
                    write(iUnitWdc,552) ekev(k),PAbn(L),taudaa,1/taudaa
                 endif
                 AFIR=DWAVE(L)/FACMU(L,i,j)/DMU(L)/WMU(L)
                 ATAW(I,J,K,L)=AFIR
                 ASEC=DWAVE(L-1)/FACMU(L,i,j)/DMU(L-1)/WMU(L)
                 GTAW(I,J,K,L)=ASEC
              enddo
           enddo
        ENDIF
     enddo
  enddo

  !.......Calculation of perpendicular current [nA/m2]
  OPEN(UNIT=iUnitJper,FILE=ST1//ST2//'.jper',STATUS='UNKNOWN')
  WRITE(iUnitJper,156) T/3600,KP     
156 FORMAT(2HT=,F4.1,2X,3HKp=,F3.1,X,   &
         'L   PHI     JPER[nA/m2]') 
  DO I=2,IO-1
     DO J=1,JO
        CJPER=-1.6E3*((PPER(I+1,J,S)-PPER(I,J,S))/DL1+3*(PPER(I,J,S)-   &
             PPAR(I,J,S))/LZ(I))/BE(I,1)
        WRITE(iUnitJper,551) LZ(I),PHI(J),CJPER
     ENDDO
  ENDDO

100 CONTINUE	
  CLOSE(iUnitJper)
  CLOSE(iUnitWans)
  CLOSE(iUnitWgr)
!  CLOSE(10)
  CLOSE(iUnitOWans)
  CLOSE(iUnitWdc)
10 FORMAT(F5.2,F10.6,5E13.4)
20 FORMAT(A80)
25 FORMAT(1X,E10.3,/10(1PE12.3))
550 format(1x,36hbounce-aver diff coeff & time scales,/,   &
         1x,4heqpa,4x,9hdaa [1/s],4x,8htdaa [s])
551 FORMAT(F5.2,F10.6,6(1PE12.3))
552 format(2F9.3,3(1PE12.3))

  RETURN
END SUBROUTINE ANISCH
!
! End of subroutine ANISCH
!
