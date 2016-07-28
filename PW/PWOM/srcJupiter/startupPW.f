 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ALEX(10/11/04): I THINK THAT THIS SUBROUTINE INITIALIZES THE GRID AND
C     AND SETS CONSTANTS. IT MUST BE INITIALIZED EVERY TIME THE CODE RUNS
C     EVEN IF RESTARTING FROM A FILE.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE STRT
C
C
      use ModCommonVariables
      use ModCommonPlanet,ONLY: HLPion1,HLPion2,HLPE,HLPE0,NDensity_CI
      use ModNumConst, ONLY:cTwoPi
      use ModPhotoElectron
      use ModPWOM  ,ONLY: UseAurora,UseIndicies, UseIE, iLine
      use ModCouplePWOMtoSE, only: get_se_for_pwom
C
      NPT1=14
      NPT2=16
      NPT3=30
      NPT4=35
      NPT5=60
      NPT6=70
1     FORMAT(5X,I5)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     DEFINE THE GAS SPECIFIC HEAT AT CONSTANT PRESSURE (CP),          C
C           THE AVERAGE MOLECULAR MASS (AVMASS), THE ACTUAL            C
C           GAS CONSTANT (RGAS) AND THE SPECIFIC HEAT RATIO (GAMMA)    C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C

CALEX not I use O for H3 and he for H2
C Gas constant = k_Boltzmann/AMU
      RGAS=8.314E7
C Adiabatic index
      GAMMA=5./3.
C AMU in gramms
      XAMU=1.6606655E-24
C Mass of atomic H3 in gramms
      Mass_I(Ion1_)=3.0237*XAMU
C Mass of atomic H in gramms
      Mass_I(Ion2_)=1.00797*XAMU
C Mass of H2 in gramms
      Mass_I(Ion3_)=2.0159*XAMU
C Mass of electron in gramms
      Mass_I(nIon)=9.109534E-28
C Relative mass of H3 to electron
      MassElecIon_I(Ion1_)=Mass_I(nIon)/Mass_I(Ion1_)
C Relative mass of atomic H to electron
      MassElecIon_I(Ion2_)=Mass_I(nIon)/Mass_I(Ion2_)
C Relative mass of H2 to electron
      MassElecIon_I(Ion3_)=Mass_I(nIon)/Mass_I(Ion3_)
C kB/m_H3
      RGAS_I(Ion1_)=RGAS*XAMU/Mass_I(Ion1_)
C kB/m_H
      RGAS_I(Ion2_)=RGAS*XAMU/Mass_I(Ion2_)
C kB/m_H2
      RGAS_I(Ion3_)=RGAS*XAMU/Mass_I(Ion3_)
C kB/m_e
      RGAS_I(nIon)=RGAS*XAMU/Mass_I(nIon)
      GMIN1=GAMMA-1.
      GMIN2=GMIN1/2.
      GPL1=GAMMA+1.
      GPL2=GPL1/2.
      GM12=GMIN1/GAMMA/2.
      GRAR=GAMMA/GMIN2
      GREC=1./GAMMA
      CPion1=GAMMA*RGAS_I(Ion1_)/GMIN1
      CPion2=GAMMA*RGAS_I(Ion2_)/GMIN1
      CPion3=GAMMA*RGAS_I(Ion3_)/GMIN1
      CPE=GAMMA*RGAS_I(nIon)/GMIN1
      CVion1=RGAS_I(Ion1_)/GMIN1
      CVion2=RGAS_I(Ion2_)/GMIN1
      CVion3=RGAS_I(Ion3_)/GMIN1
      CVE=RGAS_I(nIon)/GMIN1


c      Omega=1./37800.

      Omega=1.0/35280.0 * cTwoPi
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     READ FIELD ALIGNED CURRENT DENSITY AT LOWEST GRID POINT          C
C        (UNIT=AMPERE/M**2)                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C      READ (5,2) CURR(1)

      call set_vertical_grid

! Jupiter-specific
      do I=1,nDim
         GRAVTY(I)=-1.2657786e23/RAD(I)**2 ! accel in cm/s^2, w/ Rad in cm
         Centrifugal(I)=RAD(I)*((sin((90.-GLAT)*3.14159/180.))**2)*Omega**2
      enddo

      CURR(1)=2.998E2*CURR(1)
!      CURTIM=150.
!      CURTIM0=500.
      CURRMN=CURR(1)*(RAD(1)/(RAD(1)-DRBND))**NEXP
      CURRMX=CURR(1)*(RAD(1)/(RAD(NDIM)+DRBND))**NEXP
      
      do k=2,nDim
         CURR(k)=CURR(1)*(RAD(1)/(RAD(k)))**NEXP
      enddo



C      SGN1=1.
C      IF (CURR(1).LT.0.) SGN1=-1.
C      IF (ABS(CURR(1)).LT.1.E-4) SGN1=0.
C      CURRMX=SGN1*0.2998*(RAD(1)/(RAD(NDIM)+DRBND))**NEXP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     DEFINE THE NEUTRAL ATMOSPHERE MODEL                              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C IF GEOMAGNETIC COORDINATES ARE SPECIFIED, SET IART TO 1 !!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     DEFINE DATE                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SOLAR MAXIMUM WINTER
C      IYD=80360
C      F107A=180.
C      F107=180.
C SOLAR MAXIMUM SUMMER
C      IYD=80172
C      F107A=180.
C      F107=180.
C SOLAR MINIMUM WINTER
C      IYD=84360
C      F107A=60.
C      F107=60.
C SOLAR MINIMUM SUMMER
C     IYD=84172
C     F107A=60.
C     F107=60.
C      SEC=43200.
C      STL=12.
C      GMLAT=80.
C      GMLONG=0.
C      IART=1
C      GLAT=80.
C      GLONG=0.
C      IART=0
C FEB 20, 1990
CALEX IYD=year_day of year

c      IYD=90051
c      SEC=20.75*3600.
c      F107A=180.
c      F107=189.5
c      IART=0
c      GLONG=325.13
c      GLAT=70.47
C END
CALEX GGM determines geomagnetic lat. from geographic lat and lon
CALEX since the dipole on Saturn is aligned with the rotation axis
CALEX I have set GMLONG=GLONG
CALEX      CALL GGM(IART,GLONG,GLAT,GMLONG,GMLAT)

      DO 49 I=1,7
      AP(I)=50.
49    CONTINUE 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C                                                                      C
      !assum GMlat and GMlong is same in glat and glon and smlat smlon
      gmLat=SmLat
      gmLon=SmLon
      gLat=SmLat
      gLon=SmLon
      gLat2=-SmLat
      gLon2=SmLon
      
      CALL JupiterAtmos(XH2,XH,XH2O,XCH4,XTN)
      NDensity_CI(1:nDim,H2_) = XH2(1:NDIM)
      NDensity_CI(1:nDim,H_)  = XH(1:NDIM)
      NDensity_CI(1:nDim,H2O_)= XH2O(1:NDIM)
      NDensity_CI(1:nDim,CH4_)= XCH4(1:NDIM)


      !get the SE fluxes from SE first (call here to get Ionization rate
      if (.not.allocated(SeDens_C)) allocate(SeDens_C(nDim))
      if (.not.allocated(SeFlux_C)) allocate(SeFlux_C(nDim))
      if (.not.allocated(SeHeat_C)) allocate(SeHeat_C(nDim))
      if (.not.allocated(PhotoIonRate_IC)) 
     &     allocate(PhotoIonRate_IC(nIon-1,nDim))
      if (.not.allocated(SecIonRate_IC)) 
     &     allocate(SecIonRate_IC(nIon-1,nDim))
      if (.not.allocated(TotalIonRate_IC)) 
     &     allocate(TotalIonRate_IC(nIon-1,nDim))
      
      if ((floor((Time+1.0e-5)/DtGetSe)/=floor((Time+1.0e-5-DT)/DtGetSe))
     &     .and.DoCoupleSE) then 
         call get_se_for_pwom(Time,UTsec,iLine,(/GmLat,GmLon/),
     &        (/GLAT,GLONG/),(/GLAT2,GLONG2/),
     &        State_GV(1:nDim,RhoE_)/Mass_I(nIon),State_GV(1:nDim,Te_),
     &        Efield(1:nDim),Ap,F107,F107A,IYD,SeDens_C, SeFlux_C, SeHeat_C,
     &        PhotoIonRatePW_IC=PhotoIonRate_IC,
     &        SecIonRatePW_IC=SecIonRate_IC)

         TotalIonRate_IC(:,:) = PhotoIonRate_IC(:,:) + SecIonRate_IC(:,:)
         ! Divide the Ionization rate from SE by oxygen density to get
         ! production in units of ions/cc/s rather than ions/s
         !IonRateO_C(1:nDim)=IonRateO_C(1:nDim)/XO(1:nDim)
      endif
!      write(*,*) 'PhotoIonRate_IC'
!      write(*,*) PhotoIonRate_IC(1:3,1)
!      write(*,*) PhotoIonRate_IC(1:3,2)
      write(*,*) PhotoIonRate_IC(1:3,3),SecIonRate_IC(1:3,3)
      write(*,*) TotalIonRate_IC(1:3,3)
!      stop
      if((.not.DoCoupleSE) .or. (.not.UseFeedbackFromSE)) then
         SeDens_C(:)=0.0
         SeFlux_C(:)=0.0
         SeHeat_C(:)=0.0
      endif
C
      CALL STRT1
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     DEFINE TOPSIDE ELECTRON HEAT FLUX AND PARAMETRIC HEAT SOURCES    C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C      READ (5,2) ETOP,ELFXIN

CALEX I don't know how to fix this heat input for Saturn.
CALEX In Dee's thesis she says you need electon heat flux to be a minimum
CALEX of 20E-3 ergs cm^2 /s.      
C      ETOP=1.0E-3
C      ETOP=20.0E-3   ! original value
C      ETOP=1.0E-3    ! previously used
      ETOP=0.0        ! set to zero for SE coupling
      ELFXIN=0.
C
C      ELFXIN=9.
C
2     FORMAT(6X,1PE15.4)
     
C      READ (5,2) HEATI1,HEATI2,ELHEAT
      HEATI1=0.
      HEATI2=0.
C
      HEATI1=1.0E-7
calex origionally heatI1 was 0 and heat I2 was not 0, but i have a suspicion that
Calex this should not be true. I think the two refers to helium which we are not 
calex looking at.
c      HEATI2=2.5E-11
C
      ELHEAT=0.
      
      HEATA1=3.5E7
      HEATA2=2.0E8
      HEATA3=1.5E8
      HEATS1=2.*1.0E7**2
      HEATS2=2.*2.5E7**2
      HEATS3=2.*1.0E7**2
      DO 53 K=1,NDIM
      HEATX1=EXP(-(ALTD(K)-HEATA1)**2/HEATS1)
      HEATX2=EXP(-(ALTD(K)-HEATA2)**2/HEATS2)
      HEATX3=EXP(-(ALTD(K)-HEATA3)**2/HEATS3)
C KGS What's all this? delete?
      QOXYG(K)=(HEATI1*HEATX1+HEATI2*HEATX2)/
     #         (State_GV(K,RhoH3_)+State_GV(K,RhoH_))
calex origionally it was qhyd=qoxy/16 but I think 16 should be 3 since
calex I am letting oxy stand in for H3+
      QHYD(K)=QOXYG(K)/3.
      QHEL(K)=QOXYG(K)/4.
C
      QOXYG(K)=0.
C
      QOXYG(K)=0.
      QHYD(K)=0.
      QHEL(K)=0.

      QELECT(K)=ELHEAT*HEATX3
53    CONTINUE
      DO 54 K=1,NDIM,10
      
52    FORMAT(5(1PE15.4))
54    CONTINUE
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     READ STARTING TIME, TERMINATION TIME AND BASIC TIME STEP         C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C      READ (5,4) TIME,TMAX
CCC      TIME=0.
CCC      TMAX=1.E6
      
      
4     FORMAT(6X,2(1PE15.4))
C      READ (5,2) DT   1/20.0 is a good value
C     Read in the time step
      write(*,*) dt
c      DT=1./1.
      DTX1=DT
      DTR1=DTX1/DRBND
      DTX2=DT*NTS
      DTR2=DTX2/DRBND
    
      H0=0.5/DRBND
      H1E1=1./DTX1
      H1O1=1./DTX1
      H1H1=1./DTX1
      H1E2=1./DTX2
      H1O2=1./DTX2
      H1H2=1./DTX2
      H2=1./DRBND/DRBND
      H3=0.5*H2
      H4=0.5*H0
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     DEFINE THE HEAT CONDUCTION PARAMETERS                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CALEX terms with "surf" in them refer to surface values      
      HLPE0=GMIN1/RGAS_I(nIon)
      HLPE=1.23E-6*GMIN1/RGAS_I(nIon)
      HLPion1=2.86E-8*(Mass_I(nIon)/Mass_I(Ion1_))*GMIN1/RGAS_I(Ion1_)
      HLPion2=7.37E-8*(Mass_I(nIon)/Mass_I(Ion2_))*GMIN1/RGAS_I(Ion2_)
! KGS fix this
      HLPion3=7.37E-8*(Mass_I(nIon)/Mass_I(Ion3_))*GMIN1/RGAS_I(Ion3_)

 
      
CALEX These are the heat conductivities at the lower boundary. Note:
CALEX that no allowance is made to take into account the effect of
CALEX neutrals on the heat conduction as was done at earth.     
      HeatCon_GI(0,Ion1_)=HLPion1*(State_GV(0,RhoH3_)/State_GV(0,RhoE_))*State_GV(0,Th3_)**2.5
      HeatCon_GI(0,nIon)=HLPE*State_GV(0,Te_)**2.5
      HeatCon_GI(0,Ion2_)=HLPion2*(State_GV(0,RhoH_)/State_GV(0,RhoE_))*State_GV(0,Th_)**2.5
      HeatCon_GI(0,Ion3_)=HLPion3*(State_GV(0,RhoH2_)/State_GV(0,RhoE_))*State_GV(0,Th2_)**2.5

      
C!      HeatCon_GI(0,Ion1_)=HLPion1*State_GV(0,Th3_)**2.5
C!      HeatCon_GI(0,nIon)=HLPE*State_GV(0,Te_)**2.5
C!      HeatCon_GI(0,Ion2_)=HLPion2*State_GV(0,Th_)**2.5
C!      HeatCon_GI(0,Ion2_)E=HLPHE*State_GV(0,The_)**2.5
!*remove      CALL JupiterAtmos(1,ALTMAX,XNH2,XNH,XNH2O,XNCH4,TEMP)
      XTNMAX=XTN(NDIM+1)

      ETOP=ETOP*DRBND/1.23E-6
      CALL PW_set_upper_bc
3     FORMAT(4X,I6)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     PRINT INITIAL AND BOUNDARY PARAMETERS                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
      DO 60 I=1,NDIM
      ALTD(I)=ALTD(I)/1.E5
60    CONTINUE
      ALTMIN=ALTMIN/1.E5
      ALTMAX=ALTMAX/1.E5
      ETOP1=ETOP*1.23E-6/DRBND
      CALL COLLIS(NDIM,State_GV(-1:nDim+2,:))
      CALL PW_calc_efield(nDim,State_GV(-1:nDim+2,:))

      RETURN
      END



      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ALEX(10/11/04): I THINK THIS IS A COLD START ROUTINE FOR DEALING WITH
C     A LACK OF RESTART FILE. BASICALLY THE VELOCITY IS SET TO ZERO AND 
C     THE OTHER GRID VALUES ARE SET TO VALUES THAT WON'T CRASH.
C     (11/9/04) I also see that source terms and collision terms are set
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      
      SUBROUTINE STRT1
      use ModCommonVariables
      use ModPWOM, ONLY:IsRestart
!     REAL jp1,jp2,jp3,jp4,kc1,kc2,kc3,kc6,kc7,kc8,kr1,kr2
      real DensityHp,DensityH3p,DensityH2p
C     ALEX define the reaction rates, label by reaction number
C     ALEX j is for photochemistry, k is for regular chemistry
      
      !multiply all saturn photoproduction rates by 3.389 (Rsatfromsun/Rjupfromsun)**2
c     ! H2+hnu     --> H+ + H + e
!!!      jp1=1.9E-11*3.389
!      jp1=9.5E-11
      !             --> H2+ + e
!!!      jp2=9.9E-10*3.389
      !jp2=5.4E-10
      ! H+hnu      --> H+
!!!      jp3=1.0E-9*3.389
      !jp3=7.3E-10
      ! H2O+hnu    --> H+ + OH +e
!!!      jp4=4.2E-10*3.389
      !jp4=1.3E-10

      ! H2+ + H2    --> H3+ +H
      kc1=2.E-9
      ! H+ + H2 + M --> H3+ + M
      kc2=3.2E-29
      !kc2=0.0
      
      ! H+ + CH4    --> CH3+ + H2  
      !             --> CH4+ + H
      kc3=4.5E-9
      !kc3=4.15E-9
      !kc3=0.0
      ! H3+ + CH4   --> CH5+ + H2 
      kc6=2.4E-9
      !kc6=0.0
      ! H3+ + H2O   --> H3O+ + H2
      kc7=5.3E-9
      !kc7=0.0
      ! H+ + H2O    --> H2O+ + H
      kc8=8.2E-9
      !kc8=0.0
      ! H+ + e      --> H + hnu
      kr1=1.91E-10
      !kr1=2.0E-12
      !kr1=0.0
      ! H3+ + e     --> H2 + H
      !             --> 3H
      kr2=1.73E-6
      !kr2=4.6E-6
!note kr2 includes H3+ + e --> H2 + H and --> 3H

      ! H+ + H2(nu>3) --> H2+ + H
      do i=1,nDim 
         call get_rate(ALTD(i)/100000.0,kc9(i))
      enddo
      
      
      !kc9(:)=0.0*kc9(:)
C     
C     
C     
C     
C     PHIOX=7.00E-7
C     
C     DEFINE THE HE PHOTOIONIZATION RATE
C     
      PHIHE=1.30E-7
C     
C     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     C
C     DEFINE THE GAS PARAMETERS AT THE LOWER BOUNDARY                  C
C     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     C
      State_GV(-1:0,uH3_)=0.
      State_GV(-1:0,uH_)=0.
      State_GV(-1:0,uH2_)=0.
      State_GV(-1:0,uE_)=0.
! *remove      CALL JupiterAtmos(1,ALTMIN,XNH2,XNH,XNH2O,XNCH4,TEMP)

CALEX I pretend that for plasma parameters, O is H3 and HE is
CALEX chemical equilibrium value for H2+ this allow me to just
CALEX change the chemistry but leave the rest of the code the same  
      State_GV(-1:0,Th3_)=XTN(0)
      State_GV(-1:0,Th_)=XTN(0)
      State_GV(-1:0,Th2_)=XTN(0)
      State_GV(-1:0,Te_)=XTN(0)

C KGS this subroutine needs to be modified
      call calc_chemical_equilibrium(DensityHp,DensityH3p,DensityH2p)
      State_GV(-1:0,RhoH3_)=Mass_I(Ion1_)*DensityH3p
      State_GV(-1:0,RhoH_)=Mass_I(Ion2_)*DensityHp
      State_GV(-1:0,RhoH2_)=Mass_I(Ion3_)*DensityH2p
      write(*,*) 'H+(1400km)=',DensityHp,', H3+(1400km)=',DensityH3p,
     &     ', H2+(1400km)=',DensityH2p


C I have used numerically calculated chemical equilibrium
C solution for T=800k.       
!      State_GV(0,RhoH3_)=Mass_I(Ion1_)*6489.69
!c      State_GV(0,RhoHe_)=Mass_I(Ion3_)*jp2/kc1
!      State_GV(0,RhoH_)E=0.
!      State_GV(0,RhoH_)=Mass_I(Ion2_)*1343.64
C I have used numerically calculated chemical equilibrium
C solution for T=1000k.       
c      State_GV(0,RhoH3_)=Mass_I(Ion1_)*4725.0
c      State_GV(0,RhoH_)E=0.
c      State_GV(0,RhoH_)=Mass_I(Ion2_)*368.0

C I have used numerically calculated chemical equilibrium
C solution for T=1500k.       
c      State_GV(0,RhoH3_)=Mass_I(Ion1_)*560.0
c      State_GV(0,RhoH_)E=0.
c      State_GV(0,RhoH_)=Mass_I(Ion2_)*30.0


C I have used numerically calculated chemical equilibrium
C solution for T=1500k. with enhanced water and decreased CH4      
c      State_GV(0,RhoH3_)=Mass_I(Ion1_)*11435.41
c      State_GV(0,RhoH_)E=0.
c      State_GV(0,RhoH_)=Mass_I(Ion2_)*1463.48

C I have used numerically calculated chemical equilibrium
C solution for T=100k. with reduced CH4 enhanced h2o      
C      State_GV(0,RhoH3_)=Mass_I(Ion1_)*5509.0
C      State_GV(0,RhoH_)E=0.
C      State_GV(0,RhoH_)=Mass_I(Ion2_)*1124.0


      State_GV(-1:0,RhoE_)=MassElecIon_I(Ion3_)*State_GV(-1:0,RhoH2_) +
     &     MassElecIon_I(Ion2_)*State_GV(-1:0,RhoH_) + 
     &     MassElecIon_I(Ion1_)*State_GV(-1:0,RhoH3_)
      State_GV(-1:0,pH3_)=RGAS_I(Ion1_)*State_GV(-1:0,Th3_)*State_GV(-1:0,RhoH3_)
      State_GV(-1:0,pH_)=RGAS_I(Ion2_)*State_GV(-1:0,Th_)*State_GV(-1:0,RhoH_)
      State_GV(-1:0,pH2_)=RGAS_I(Ion3_)*State_GV(-1:0,Th2_)*State_GV(-1:0,RhoH2_)
      State_GV(-1:0,pE_)=RGAS_I(nIon)*State_GV(-1:0,Te_)*State_GV(-1:0,RhoE_)
      SoundSpeed_GI(0,Ion1_)=SQRT(GAMMA*RGAS_I(Ion1_)*State_GV(0,Th3_))
      SoundSpeed_GI(0,Ion2_)=SQRT(GAMMA*RGAS_I(Ion2_)*State_GV(0,Th_))
      SoundSpeed_GI(0,Ion3_)=SQRT(GAMMA*RGAS_I(Ion3_)*State_GV(0,Th2_))
      SoundSpeed_GI(0,nIon)=SQRT(GAMMA*RGAS_I(nIon)*State_GV(0,Te_))

C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     

! Start collision freq at zero
      CollisionFreq_IIC(:,:,:)= 0.0
      HeatFlowCoef_II(:,:)    = 0.0
      FricHeatCoef_II(:,:)    = 0.0


      DO 20 I=1,NDIM
CALEX SOURCE COEF?
C jp = photochemical reaction rate
C kc = collisional reaction rate      
         FFHpp1(I)=jp1*XH2(I)
         FFHpp3(I)=jp3*XH(I)
         FFHpp4(I)=jp4*XH2O(I)
         FFHpc2(I)=-kc2*XH2(I)*XH2(I)
         FFHpc3(I)=-kc3*XCH4(I)
         FFHpc8(I)=-kc8*XH2O(I)
         FFHpc9(I)=-kc9(I)*XH2(I)
         FFHpr1(I)=-kr1
CALEX write out source coeff
CALEX         write(26,*) FFHpp1(I),FFHpp3(I),FFHpp4(I),FFHpc2(I),FFHpc3(I),FFHpc8(I),FFHpr1(I)
        

!         FFH3pc1(I)=kc1*(jp2/kc1)*XH2(I)   ! KGS H2+ folded in here
         FFH3pc1(I)=kc1*XH2(I)              ! KGS Make sure this is multiplied by H2+ now
         FFH3pc2(I)=kc2*XH2(I)*XH2(I)
         FFH3pc6(I)=-kc6*XCH4(I)
         FFH3pc7(I)=-kc7*XH2O(I)
         FFH3pr2(I)=-kr2

CALEX write out source coeff
CALEX         write(27,*) FFH3pc1(I),FFH3pc2(I),FFH3pc6(I),FFH3pc7(I),FFH3pr2(I)

         FFH2pp2(I)=jp2*XH2(I)
         FFH2pc9(I)=kc9(I)*XH2(I)
         FFH2pc1(I)=-kc1*XH2(I)

CALEX CL=COLLISION COEF, CF=collision freq ?         
C KGS CF(1,2) = CL(1,2) * n(2)
C KGS n(2)m(2)CF(2,1) = n(1)m(1)CF(1,2)
CALEX the coulomb collisions
CAlex H+ and H3+
c KGS and H2+   from eq 4.143?
         CLHpH3p(I)=1.905*4.**1.5/Mass_I(Ion1_) ! KGS orig
         CLHpH3p(I)=1.100*4.**1.5/Mass_I(Ion1_) ! KGS correct value?
         CLH2pH3p(I)=3.17*6.**0.5/Mass_I(Ion1_) ! KGS double-check
         CLH2pHp(I)=3.81*2.**0.5/Mass_I(Ion2_)  ! KGS double-check
CALEX electron H+ and electron H3+
c KGS and electron H2+
         CLELHp(I)=54.5/Mass_I(Ion2_)
         CLELH3p(I)=54.5/Mass_I(Ion1_)
         CLELH2p(I)=54.5/Mass_I(Ion3_)

CALEX  ion neutrals
C KGS looks similar to Nagy 4.88 Maxwell Molecule Collisions
C KGS but some terms seem to be missing
         ! H+ - H2
         CLHpH(I)=2.65E-10*XH(I)
!         CFHpH2(I)=2.6E-9*XH2(I)*(.82/.667)**.5
         CollisionFreq_IIC(Ion2_,Neutral1_,I)=
     &        2.6E-9*XH2(I)*(.82/.667)**.5

! KGS 2-2 Coulomb Collision is in collisionPW.f

         ! H3+ - H
!         CFH3pH(I)=2.6E-9*XH(I)*(.667/.75)**.5
         CollisionFreq_IIC(Ion1_,Neutral2_,I)=
     &        2.6E-9*XH(I)*(.667/.75)**.5

         ! H3+ - H2
!         CFH3pH2(I)=2.6E-9*XH2(I)*(.82/1.2)**.5
         CollisionFreq_IIC(Ion1_,Neutral1_,I)=
     &        2.6E-9*XH2(I)*(.82/1.2)**.5

C KGS make sure I should be using Maxwell Molecule Collisions here
         ! H2+ - H2
         CollisionFreq_IIC(Ion3_,Neutral1_,I)=
     &        2.6E-9*XH2(I)*(.82/1.0)**.5    ! KGS check these
         ! H2+ - H
         CollisionFreq_IIC(Ion3_,Neutral2_,I)=
     &        2.6E-9*XH(I)*(.667/.667)**.5

CALEX electron H, e H2 done in collis
         CLELH(I)=4.5E-9*XH(I)

         
 20   CONTINUE
c      GRAVTY(nDim)=GRAVTY(nDim)*0.0
cAlex
      NEXP=3
      WRITE(*,*) CURR(1),RAD(1),NEXP
cendAlex      
      CRRX=CURR(1)*RAD(1)**NEXP
      DO 25 I=1,NDIM
         CURR(I)=CRRX/RAD(I)**NEXP
 25   CONTINUE
C     SGN1=1.
C     IF (CURR(1).LT.0.) SGN1=-1.
C     IF (ABS(CURR(1)).LT.1.E-4) SGN1=0.
C     DO 30 I=NDIM-250,NDIM
C     CURR(I)=SGN1*0.2998*(RAD(1)/RAD(I))**NEXP
C     30    CONTINUE
      
CALEX CT & CM = energy collision term coef?
CALEX Based on Nagy p. 83, I believe that CT is the coeff
CALEX of the term that is due to temperature difference or
CALEX heat flow between species, and CM is the term due to 
CALEX frictional heating between species moving through each other

CALEX CTOXN2 = 3*R_o*M_o/(M_o+M_{N2}) see nagy p.83
      !CTHpH=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion2_))
      !CTHpH2=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+2.*Mass_I(Ion2_))
      !CTHpH3p=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion1_))
      !CTHpEL=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(nIon))
      
      HeatFlowCoef_II(Ion3_,Neutral2_)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion2_))
      HeatFlowCoef_II(Ion3_,Neutral1_)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion3_))
      HeatFlowCoef_II(Ion3_,Ion1_)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion1_))
      HeatFlowCoef_II(Ion3_,Ion2_)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion2_))
      HeatFlowCoef_II(Ion3_,nIon)=3.*RGAS_I(Ion3_)*Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(nIon))

      MassFracCoef_II(Ion3_,:) = HeatFlowCoef_II(Ion3_,:) / (3.0*RGAS_I(Ion3_))

      HeatFlowCoef_II(Ion2_,Neutral2_)=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion2_))
      HeatFlowCoef_II(Ion2_,Neutral1_)=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+2.*Mass_I(Ion2_))
      HeatFlowCoef_II(Ion2_,Ion1_)=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion1_))
      HeatFlowCoef_II(Ion2_,Ion3_)=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion3_))
      HeatFlowCoef_II(Ion2_,nIon)=3.*RGAS_I(Ion2_)*Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(nIon))

      MassFracCoef_II(Ion2_,:) = HeatFlowCoef_II(Ion2_,:) / (3.0*RGAS_I(Ion2_))
!      CTH3pH=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
!      CTH3pH2=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+2.*Mass_I(Ion2_))
!      CTH3pHp=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
!      CTH3pEL=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(nIon))

      HeatFlowCoef_II(Ion1_,Neutral2_)=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
      HeatFlowCoef_II(Ion1_,Neutral1_)=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+2.*Mass_I(Ion2_))
      HeatFlowCoef_II(Ion1_,Ion2_)=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
      HeatFlowCoef_II(Ion1_,Ion3_)=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(Ion3_))
      HeatFlowCoef_II(Ion1_,nIon)=3.*RGAS_I(Ion1_)*Mass_I(Ion1_)/(Mass_I(Ion1_)+Mass_I(nIon))

      MassFracCoef_II(Ion1_,:) = HeatFlowCoef_II(Ion1_,:) / (3.0*RGAS_I(Ion1_))
!      CTELH=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+Mass_I(Ion2_))
!      CTELH2=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+2.*Mass_I(Ion2_))
!      CTELHp=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+Mass_I(Ion2_))
!      CTELH3p=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+Mass_I(Ion1_))
      
      HeatFlowCoef_II(nIon,Neutral2_)=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+Mass_I(Ion2_))
      HeatFlowCoef_II(nIon,Neutral1_)=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+2.*Mass_I(Ion2_))
      HeatFlowCoef_II(nIon,Ion3_)=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+Mass_I(Ion3_))
      HeatFlowCoef_II(nIon,Ion2_)=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+Mass_I(Ion2_))
      HeatFlowCoef_II(nIon,Ion1_)=3.*RGAS_I(nIon)*Mass_I(nIon)/(Mass_I(nIon)+Mass_I(Ion1_))

      MassFracCoef_II(nIon,:) = HeatFlowCoef_II(nIon,:) / (3.0*RGAS_I(nIon))

CALEX CMOXN2 = M_{N2}/(M_o+M_{N2}) see nagy p.83
!      CMHpH=Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion2_))
!      CMHpH2=2.*XAMU/(Mass_I(Ion2_)+2.*XAMU)
!      CMHpH3p=Mass_I(Ion1_)/(Mass_I(Ion2_)+Mass_I(Ion1_))
!      CMHpEL=Mass_I(nIon)/(Mass_I(Ion2_)+Mass_I(nIon))

      FricHeatCoef_II(Ion3_,Neutral2_)=Mass_I(Ion3_)/(Mass_I(Ion3_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion3_,Neutral1_)=2.*XAMU/(Mass_I(Ion3_)+2.*XAMU)
      FricHeatCoef_II(Ion3_,Ion2_)=Mass_I(Ion2_)/(Mass_I(Ion3_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion3_,Ion1_)=Mass_I(Ion1_)/(Mass_I(Ion3_)+Mass_I(Ion1_))
      FricHeatCoef_II(Ion3_,nIon)=Mass_I(nIon)/(Mass_I(Ion3_)+Mass_I(nIon))

      FricHeatCoef_II(Ion2_,Neutral2_)=Mass_I(Ion2_)/(Mass_I(Ion2_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion2_,Neutral1_)=2.*XAMU/(Mass_I(Ion2_)+2.*XAMU)
      FricHeatCoef_II(Ion2_,Ion3_)=Mass_I(Ion3_)/(Mass_I(Ion2_)+Mass_I(Ion3_))
      FricHeatCoef_II(Ion2_,Ion1_)=Mass_I(Ion1_)/(Mass_I(Ion2_)+Mass_I(Ion1_))
      FricHeatCoef_II(Ion2_,nIon)=Mass_I(nIon)/(Mass_I(Ion2_)+Mass_I(nIon))
       
!      CMH3pH=Mass_I(Ion2_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
!      CMH3pH2=2.*XAMU/(Mass_I(Ion1_)+2.*XAMU)
!      CMH3pHp=Mass_I(Ion2_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
!      CMH3pEL=Mass_I(nIon)/(Mass_I(Ion1_)+Mass_I(nIon))
       
      FricHeatCoef_II(Ion1_,Neutral2_)=Mass_I(Ion2_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion1_,Neutral1_)=2.*XAMU/(Mass_I(Ion1_)+2.*XAMU)
      FricHeatCoef_II(Ion1_,Ion3_)=Mass_I(Ion3_)/(Mass_I(Ion1_)+Mass_I(Ion3_))
      FricHeatCoef_II(Ion1_,Ion2_)=Mass_I(Ion2_)/(Mass_I(Ion1_)+Mass_I(Ion2_))
      FricHeatCoef_II(Ion1_,nIon)=Mass_I(nIon)/(Mass_I(Ion1_)+Mass_I(nIon))

!      CMELH=Mass_I(Ion2_)/(Mass_I(nIon)+Mass_I(Ion2_))
!      CMELH2=2.*XAMU/(Mass_I(nIon)+2.*XAMU)
!      CMELHp=Mass_I(Ion2_)/(Mass_I(nIon)+Mass_I(Ion2_))
!      CMELH3p=Mass_I(Ion1_)/(Mass_I(nIon)+Mass_I(Ion1_))
      
      FricHeatCoef_II(nIon,Neutral2_)=Mass_I(Ion2_)/(Mass_I(nIon)+Mass_I(Ion2_))
      FricHeatCoef_II(nIon,Neutral1_)=2.*XAMU/(Mass_I(nIon)+2.*XAMU)
      FricHeatCoef_II(nIon,Ion3_)=Mass_I(Ion3_)/(Mass_I(nIon)+Mass_I(Ion3_))
      FricHeatCoef_II(nIon,Ion2_)=Mass_I(Ion2_)/(Mass_I(nIon)+Mass_I(Ion2_))
      FricHeatCoef_II(nIon,Ion1_)=Mass_I(Ion1_)/(Mass_I(nIon)+Mass_I(Ion1_))
      
C     ALEX(10/11/04): 
C     TRY SETTING THE PLASMA PARAMETERS HERE TO THE SURFACE VALUES

      print *,'Restart: ',IsRestart

      if(IsRestart) RETURN
      IsRestart = .true.
      ! this needs to be changed
      do K=1,NDIM
         State_GV(K,uH3_)=0
         State_GV(K,pH3_)=State_GV(0,pH3_)*exp(-(ALTD(k)-1400.E5)/5000.E5)
         State_GV(K,RhoH3_)=State_GV(0,RhoH3_)*exp(-(ALTD(k)-1400.E5)/5000.E5)
         State_GV(K,Th3_)=State_GV(0,Th3_)
         State_GV(K,uH_)=0
         State_GV(K,pH_)=State_GV(0,pH_)*exp(-(ALTD(k)-1400.E5)/5000.E5)
         State_GV(K,RhoH_)=State_GV(0,RhoH_)*exp(-(ALTD(k)-1400.E5)/5000.E5)
         State_GV(K,Th_)=State_GV(0,Th_)
         State_GV(K,uH2_)=0
         State_GV(K,pH2_)=State_GV(0,pH2_)*exp(-(ALTD(k)-1400.E5)/5000.E5)
         State_GV(K,RhoH2_)=State_GV(0,RhoH2_)*exp(-(ALTD(k)-1400.E5)/5000.E5)
         State_GV(K,Th2_)=State_GV(0,Th2_)
         State_GV(K,RhoE_)=State_GV(0,RhoE_)*exp(-(ALTD(k)-1400.E5)/5000.E5)
         State_GV(K,uE_)=0
         State_GV(K,pE_)=State_GV(0,pE_)*exp(-(ALTD(k)-1400.E5)/5000.E5)
         State_GV(K,Te_)=State_GV(0,Te_)
         
         
      enddo
      
      RETURN
      END
