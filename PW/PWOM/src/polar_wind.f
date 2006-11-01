CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     THIS PROGRAM SOLVES THE SIMPLE HYDRODYNAMIC EQUATIONS WITH       C
C     HEAT CONDUCTION AND FRICTION.                                    C
C                                                                      C
C     LIST OF LOGICAL I/O UNITS:                                       C
C                                                                      C
C          5  -  INPUT CONTROL CARDS                                   C
C          16 -  DETAILED LINE PRINTER OUTPUT                          C
C          7  -  RESTART OUTPUT FILE                                   C
C          8  -  RESTART INPUT FILE                                    C
C          9  -  COMPRESSED OUTPUT FILE FOR COLOR GRAPHICS             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCiUnitRestartCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C photoelectron impact ionization modified for no conjugate hemisphere
C P. E. sources. See Subroutine Glowex, variable 'rat' D. Cannata 9May90
C 
C ADDED SUBROUTINE CALL FOR IMPACT IONIZATION "call precip" PRESENTLY
C HARD WIRED FOR POLAR RAIN FOR INVARIENT LATITUDES GREATER THAN 75 DEG.
C CHANGES FOR CUSP OR AURORAL PRECIPITATION CAN BE MADE BY CHANGING
C ALF1(x) AND FLUX1(x) ARRAY, x=1->4. MODIFIFIED BY DICK AND STEVE 6-90.
C
C 2 time step acceleration implemented. Subroutines with names ending in 
C 1 or 2 refer to lower and upper regions of flux tube; those ending in
C W treat the whole flux tube. NTS = ratio of larger time step to
C smaller time step. NCL = number of cells in lower region.
C MSTEP = counter used to keep track  of time steps; whole field line
C done for MSTEP = NTS.  S. Guiter 14May90

C Variable time step option implemented 8/05 A. Glocer

      Subroutine POLAR_WIND

      Use ModPWOM, only: DTpolarwind,nLine
      use ModIoUnit, ONLY: UnitTmp_
      use ModCommonVariables
      use ModFieldLine
      INTEGER NOTP(100)

C     define the output files and attaching units
      character*100 :: NameRestart
      Logical IsFirstCall
      Data IsFirstCall / .true./
      Real Jr
      !-----------------------------------------------------------------------

      Call get_field_line(dOxyg(:), uOxyg(:), pOxyg(:), TOxyg(:),     
     &                 dHel(:), uHel(:), pHel(:), THel(:),         
     &                 dHyd(:), uHyd(:), pHyd(:), THyd(:),         
     &                 dElect(:), uElect(:), pElect(:), TElect(:), 
     &                 GMLAT,GMLONG,Jr,wHorizontal,
     &                 iUnitOutput=iUnitOutput, iUnitGraphics=iUnitGraphics,   
     &                 NameRestart=NameRestart,    
     &                 iLine=iLine, Time=Time,MaxLineTime=Tmax,
     &                 TypeSolver=TypeSolver,IsVariableDT=IsVariableDT,
     &                 IsRestart=IsRestart,DToutput=DToutput)
      
      DT=DTpolarwind
     

      CURR(1)      = Jr
      !NameInput    = 'pw.input'
      !NameOutput   = 'log.out'
      !NameGraphics = 'plots.out'
      !NameSourceGraphics='plot_sources.out'
      !NameRestart  = 'restart.out'
      !NameCollision= 'plots_collision.out'


C
      NTS = 1
      NCL = 60


21    FORMAT(2X,I8)

      KSTEP=1
      NSTEP=1
      MSTEP=1
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     READ THE NUMBER OF CELLS                                         C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C      READ (5,21) NDIM
C     nDim is the actual gridsize used
      if (NamePlanet == 'Earth ') nDim=390
C Use nDim=550 for top boundary at 2 Re in the earth case
c      nDim=550
c      nDim=400
      if (NamePlanet == 'Saturn') nDim=800
c      if (NamePlanet == 'Saturn') nDim=930

c      nDim=680
c      nDim=1500
c      nDim=1180
c      nDim=750
      WRITE (iUnitOutput,21) NDIM
      NDIM2=NDIM-1
      NDIM1=NDIM+1
      NDIMM=NDIM+2
      
      
!!!!      read(iUnitInput,*) IsRestart
!!!!      if (IsRestart) then
!!!!         write(*,*) 'Is Restart', IsRestart
!!!!      endif

!!!!      read(iUnitInput,*) IsVariableDt
!!!!      if (IsVariableDt) then
!!!!         write(*,*) 'IsVariableDT', IsVariableDt
!!!!      endif
      

C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     DEFINE THE INITIAL CONDITIONS IN THE CELLS AND AT CELL           C
C           BOUNDARIES                                                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C

!!!!      TIME=0
      
!!!!      if(IsRestart)then
!!!!         OPEN(UNIT=8,FILE='restart.in',STATUS='OLD')
!!!!         READ (8,2001) TIME,DDT1,NS
!!!!
!!!!
!!!!         READ (8,2002) (XXX,UOXYG(K),POXYG(K),DOXYG(K),TOXYG(K),K=1,NDIM)
!!!!         READ (8,2002) (XXX,UHEL(K),PHEL(K),DHEL(K),THEL(K),K=1,NDIM)
!!!!         READ (8,2002) (XXX,UHYD(K),PHYD(K),DHYD(K),THYD(K),K=1,NDIM)
!!!!         READ (8,2002) (XXX,UELECT(K),PELECT(K),DELECT(K),TELECT(K),
!!!!     &        K=1,NDIM)
!!!!         CLOSE(UNIT=8)
 2001    FORMAT(2(1PE16.6),I10)
!!!!! 2002    FORMAT(5(1PE16.6))
 2002    FORMAT(5(1PE25.16))
!!!!      endif


111   CALL STRT



      write(iUnitOutput,*) 'time:',time

      if (IsFirstCall) then
         CALL PRNTSM
         !CALL prntCollision
         !CALL PRNT_Sources
         IsFirstCall = .false.
      endif
!******************************************************************************
! This is the start of a while loop that carries out the main steps of
! the simulation
!******************************************************************************



1     CONTINUE
      IF (MSTEP .NE. NTS) THEN
         
         if (TypeSolver == 'Godunov') then
            CALL Solver('O1',UOXYG,DOXYG,POXYG,USURFO,DSURFO
     ;           ,PSURFO,DBGNDO,PBGNDO,UBGNDO,ADMSO,ECLSNO,FCLSNO,QOXYG,RGASO,
     ;           CELLNW(3,:),CELLNW(1,:),CELLNW(2,:),CELLNW(4,:))

            CALL Solver('H1',UHYD,DHYD,PHYD,USURFH,DSURFH
     ;           ,PSURFH,DBGNDH,PBGNDH,UBGNDH,ADMSH,ECLSNH,FCLSNH,QHYD,RGASH,
     ;           CELLNW(7,:),CELLNW(5,:),CELLNW(6,:),CELLNW(8,:))

            if (NamePlanet == 'Earth ') then
               CALL Solver('He1',uHel,dHel,pHel,uSurHe,dSurHe
     ;              ,pSurHe,DBGNHe,PBGNHe,UBGNHe,ADMSHe,ECLSHe,FCLSHe,qHel,RGASHe,
     ;              CELLNW(11,:),CELLNW(9,:),CELLNW(10,:),CELLNW(12,:))
            endif
            if (NamePlanet == 'Saturn') then
               call clfmhe
            endif

         endif
         
         if (TypeSolver == 'Rusanov') then
            
            call rusanov_solver('O1',RgasO,
     &           dOxyg, uOxyg, pOxyg, tOxyg, 
     &           dSurfO, uSurfO, pSurfO, 
     &           dBgndO, uBgndO, pBgndO, 
     &           dBgndO2, uBgndO2, pBgndO2, 
     &           AdmsO, FCLSNO, ECLSNO, 
     &           CELLNW(3,:), CELLNW(1,:), CELLNW(2,:),CELLNW(4,:))

            call rusanov_solver('H1',RgasH,
     &           dHYD, uHYD, pHYD,tHyd, 
     &           dSurfH, uSurfH, pSurfH, 
     &           dBgndH, uBgndH, pBgndH, 
     &           dBgndH2, uBgndH2, pBgndH2, 
     &           AdmsH, FCLSNH, ECLSNH, 
     &           CELLNW(7,:), CELLNW(5,:), CELLNW(6,:),CELLNW(8,:))

            if (NamePlanet == 'Earth ') then
               call rusanov_solver('He1',RgasHe,
     &              dHel, uHel, pHel,tHel, 
     &              dSurHe, uSurHe, pSurHe, 
     &              dBgnHe, uBgnHe, pBgnHe, 
     &              dBgnHe2, uBgnHe2, pBgnHe2, 
     &              AdmsHe, FCLSHe, ECLSHe, 
     &              CELLNW(11,:), CELLNW(9,:), CELLNW(10,:),CELLNW(12,:))
            endif
            if (NamePlanet == 'Saturn') then
               call clfmhe
            endif

         endif

!         pcHYD = maxval(abs(PHYD-CELLNW(6,:))/CELLNW(6,:))
!         pcOXY = maxval(abs(POXYG-CELLNW(2,:))/CELLNW(2,:))
!         pc = max (pcHYD,pcOXY)

C         CURHLP=EXP(-0.5*(TIME-CURTIM0)**2/CURTIM**2)
         DO 6 K=1,NCL
            UOXYG(K)=CELLNW(1,K)
            POXYG(K)=CELLNW(2,K)
            DOXYG(K)=CELLNW(3,K)
            TOXYG(K)=CELLNW(4,K)
            UHYD(K)=CELLNW(5,K)
            PHYD(K)=CELLNW(6,K)
            DHYD(K)=CELLNW(7,K)
            THYD(K)=CELLNW(8,K)
            UHEL(K)=CELLNW(9,K)
            PHEL(K)=CELLNW(10,K)
            DHEL(K)=CELLNW(11,K)
            THEL(K)=CELLNW(12,K)
c            write(*,*)'RTHDEL*DHYD(K),RTOXEL*DOXYG(K),RTHEEL*DHEL(K)',  RTHDEL*DHYD(K),RTOXEL*DOXYG(K),RTHEEL*DHEL(K)
            DELECT(K)=RTHDEL*DHYD(K)+RTOXEL*DOXYG(K)+RTHEEL*DHEL(K)
            UELECT(K)=(RTHDEL*DHYD(K)*UHYD(K)+RTOXEL*DOXYG(K)*
     $           UOXYG(K)+RTHEEL*DHEL(K)*UHEL(K)
     $           -1.8965E-18*CURR(K))/DELECT(K)
 6       CONTINUE
!         call calc_cfl
!         if (TypeSolver == 'Godunov') then
!            call diffusion(DHYD,DSURFH,DBGNDH,cMax_H)
!            call diffusion(DOXYG,DSURFO,DBGNDO,cMax_O)
!            
!            call diffusion(UHYD,USURFH,UBGNDH,cMax_H)
!            call diffusion(UOXYG,USURFO,UBGNDO,cMax_O)
!            
!            call diffusion(PHYD,PSURFH,PBGNDH,cMax_H)
!            call diffusion(POXYG,PSURFO,PBGNDO,cMax_O)
!         endif

         
         do k=1,nDim
            DELECT(K)=RTHDEL*DHYD(K)+RTOXEL*DOXYG(K)+RTHEEL*DHEL(K)
            UELECT(K)=(RTHDEL*DHYD(K)*UHYD(K)+RTOXEL*DOXYG(K)*
     $           UOXYG(K)+RTHEEL*DHEL(K)*UHEL(K)
     $           -1.8965E-18*CURR(K))/DELECT(K)
         enddo

         
         CALL NEWBGD
         CALL ELFLDW
      ELSE
         if (TypeSolver == 'Godunov') then
            CALL Solver('O2',UOXYG,DOXYG,POXYG,USURFO,DSURFO
     ;           ,PSURFO,DBGNDO,PBGNDO,UBGNDO,ADMSO,ECLSNO,FCLSNO,QOXYG,RGASO,
     ;           CELLNW(3,:),CELLNW(1,:),CELLNW(2,:),CELLNW(4,:))
            
            CALL Solver('H2',UHYD,DHYD,PHYD,USURFH,DSURFH
     ;           ,PSURFH,DBGNDH,PBGNDH,UBGNDH,ADMSH,ECLSNH,FCLSNH,QHYD,RGASH,
     ;           CELLNW(7,:),CELLNW(5,:),CELLNW(6,:),CELLNW(8,:))
C     ALEX Check to see if the planet is earth or saturn, then do the update for He+ or H2+
            if (NamePlanet .eq. 'Earth ') then
               CALL Solver('He2',uHel,dHel,pHel,uSurHe,dSurHe
     ;              ,pSurHe,DBGNHe,PBGNHe,UBGNHe,ADMSHe,ECLSHe,FCLSHe,qHel,RGASHe,
     ;              CELLNW(11,:),CELLNW(9,:),CELLNW(10,:),CELLNW(12,:))
            endif
            if (NamePlanet .eq. 'Saturn') then
C     ALEX            do K=1,NDIM
C     ALEX               CELLNW(11,K)=XMSHE*jp2/kc1
C     ALEX               CELLNW(9,K)=0.
C     ALEX               CELLNW(10,K)=RGASHE*CELLNW(11,K)*800.
C     ALEX               CELLNW(12,K)=800.
C     ALEX            enddo
               call clfmhe
            endif

         endif

         if (TypeSolver == 'Rusanov') then
            
            call rusanov_solver('O2',RgasO,
     &           dOxyg, uOxyg, pOxyg,tOxyg, 
     &           dSurfO, uSurfO, pSurfO, 
     &           dBgndO, uBgndO, pBgndO, 
     &           dBgndO2, uBgndO2, pBgndO2, 
     &           AdmsO, FCLSNO, ECLSNO, 
     &           CELLNW(3,:), CELLNW(1,:), CELLNW(2,:),CELLNW(4,:))

            call rusanov_solver('H2',RgasH,
     &           dHYD, uHYD, pHYD, tHyd,
     &           dSurfH, uSurfH, pSurfH, 
     &           dBgndH, uBgndH, pBgndH, 
     &           dBgndH2, uBgndH2, pBgndH2, 
     &           AdmsH, FCLSNH, ECLSNH, 
     &           CELLNW(7,:), CELLNW(5,:), CELLNW(6,:),CELLNW(8,:))

            if (NamePlanet == 'Earth ') then
               call rusanov_solver('He2',RgasHe,
     &              dHel, uHel, pHel,tHel, 
     &              dSurHe, uSurHe, pSurHe, 
     &              dBgnHe, uBgnHe, pBgnHe, 
     &              dBgnHe2, uBgnHe2, pBgnHe2, 
     &              AdmsHe, FCLSHe, ECLSHe, 
     &              CELLNW(11,:), CELLNW(9,:), CELLNW(10,:),CELLNW(12,:))
            endif
            if (NamePlanet == 'Saturn') then
               call clfmhe
            endif
         endif
          

!         pcHYD = maxval(abs(PHYD-CELLNW(6,:))/CELLNW(6,:))

!         pcOXY = maxval(abs(POXYG-CELLNW(2,:))/CELLNW(2,:))

!         pc = max (pcHYD,pcOXY)


!         CURHLP=EXP(-0.5*(TIME-CURTIM0)**2/CURTIM**2)
!         CURHLP=EXP(-0.5*(TIME-CURTIM0)**2/CURTIM**2)
         DO 2 K=1,NDIM
            UOXYG(K)=CELLNW(1,K)
            POXYG(K)=CELLNW(2,K)
            DOXYG(K)=CELLNW(3,K)
            TOXYG(K)=CELLNW(4,K)
            UHYD(K)=CELLNW(5,K)
            PHYD(K)=CELLNW(6,K)
            DHYD(K)=CELLNW(7,K)
            THYD(K)=CELLNW(8,K)
            UHEL(K)=CELLNW(9,K)
            PHEL(K)=CELLNW(10,K)
            DHEL(K)=CELLNW(11,K)
            THEL(K)=CELLNW(12,K)

            DELECT(K)=RTHDEL*DHYD(K)+RTOXEL*DOXYG(K)+RTHEEL*DHEL(K)
            UELECT(K)=(RTHDEL*DHYD(K)*UHYD(K)+RTOXEL*DOXYG(K)*
     $           UOXYG(K)+RTHEEL*DHEL(K)*UHEL(K)
     $           -1.8965E-18*CURR(K))/DELECT(K)
 2       CONTINUE
!         call calc_cfl
!         if (TypeSolver == 'Godunov') then
!            call diffusion(DHYD,DSURFH,DBGNDH,cMax_H)
!            call diffusion(DOXYG,DSURFO,DBGNDO,cMax_O)
!            
!            call diffusion(UHYD,USURFH,UBGNDH,cMax_H)
!            call diffusion(UOXYG,USURFO,UBGNDO,cMax_O)
!            
!            call diffusion(PHYD,PSURFH,PBGNDH,cMax_H)
!            call diffusion(POXYG,PSURFO,PBGNDO,cMax_O)
!         endif            

         do k=1,nDim
            DELECT(K)=RTHDEL*DHYD(K)+RTOXEL*DOXYG(K)+RTHEEL*DHEL(K)
            UELECT(K)=(RTHDEL*DHYD(K)*UHYD(K)+RTOXEL*DOXYG(K)*
     $           UOXYG(K)+RTHEEL*DHEL(K)*UHEL(K)
     $           -1.8965E-18*CURR(K))/DELECT(K)

         enddo
         

        

         CALL CLTEMPW
         CALL CLFME2W
         CALL NEWBGD
         CALL COLLIS(NDIM)
         CALL ELFLDW
         
      ENDIF
      TIME=TIME+DT
!      call calc_cfl
c      call calcdt(pc)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC               
C     This determines whether or not to output a restart file.
C     There are two methods. Either output on iteration for a 
C     Fixed timestep, or output on frequency for a variable DT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (.not. IsVariableDt) then
         IF (NSTEP.EQ.NOTP(KSTEP)) Then
            KSTEP=KSTEP+1
c            CALL PRNTOP
c            CALL PRNTHE
c            CALL PRNTHP
c            CALL PRNTEL
c            CALL PRNTEF
            write(iUnitOutput,*) 'time:',time
            CALL PRNTSM
            !CALL prntCollision
            !CALL PRNT_Sources
            open(UnitTmp_, FILE=NameRestart)

            WRITE (UnitTmp_,2001) TIME,DT,NSTEP
            Write (UnitTmp_,*) GMLAT, GMLONG

            WRITE (UnitTmp_,2002)(RAD(K),UOXYG(K),POXYG(K),DOXYG(K),TOXYG(K),K=1,NDIM)
            WRITE (UnitTmp_,2002)(RAD(K),UHEL(K),PHEL(K),DHEL(K),THEL(K),K=1,NDIM)
            WRITE (UnitTmp_,2002)(RAD(K),UHYD(K),PHYD(K),DHYD(K),THYD(K),K=1,NDIM)
            WRITE (UnitTmp_,2002)(RAD(K),UELECT(K),PELECT(K),DELECT(K),TELECT(K),
     ;           K=1,NDIM)
            close(UnitTmp_)
         endif
      endif


      IF (IsVariableDt) then
         IF (floor(Time/DToutput).NE.floor((Time-DT)/DToutput) ) Then
            KSTEP=KSTEP+1
c            CALL PRNTOP
c            CALL PRNTHE
c            CALL PRNTHP
c            CALL PRNTEL
c            CALL PRNTEF
            write(iUnitOutput,*) 'time:',time
            CALL PRNTSM
            !CALL prntCollision
            !CALL PRNT_Sources


            open(UnitTmp_, FILE=NameRestart)
            WRITE (UnitTmp_,2001) TIME,DT,NSTEP
            Write (UnitTmp_,*) GMLAT, GMLONG
            WRITE (UnitTmp_,2002)(RAD(K),UOXYG(K),POXYG(K),DOXYG(K),TOXYG(K),K=1,NDIM)
            WRITE (UnitTmp_,2002)(RAD(K),UHEL(K),PHEL(K),DHEL(K),THEL(K),K=1,NDIM)
            WRITE (UnitTmp_,2002)(RAD(K),UHYD(K),PHYD(K),DHYD(K),THYD(K),K=1,NDIM)
            WRITE (UnitTmp_,2002)(RAD(K),UELECT(K),PELECT(K),DELECT(K),TELECT(K),
     ;           K=1,NDIM)
            close(UnitTmp_)
         endif
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!      IF (TIME.GT.TMAX) Then
!         Call Put_Mod_PW(dOxyg, uOxyg, pOxyg, TOxyg,     
!     &                   dHel, uHel, pHel, THel,         
!     &                   dHyd, uHyd, pHyd, THyd,         
!     &                   dElect, uElect, pElect, TElect, 
!     &                   IsRestart,IsVariableDt,Time,DT,DToutput,DT,
!     &                   TypeSolver,GMLAT,GMLONG) 
!
!         RETURN
!      endif
      NSTEP=NSTEP+1
      IF (MSTEP .EQ. NTS) MSTEP=0 
      MSTEP=MSTEP+1
      if (.not. IsVariableDT) then
         IF (KSTEP.GT.NPT) RETURN
      endif
      IF (MSTEP .NE. NTS) THEN
         if (TypeSolver == 'Godunov') then
            CALL Solver('O3',UOXYG,DOXYG,POXYG,USURFO,DSURFO
     ;           ,PSURFO,DBGNDO,PBGNDO,UBGNDO,ADMSO,ECLSNO,FCLSNO,QOXYG,RGASO,
     ;           CELLNW(3,:),CELLNW(1,:),CELLNW(2,:),CELLNW(4,:))
            
            CALL Solver('H3',UHYD,DHYD,PHYD,USURFH,DSURFH
     ;           ,PSURFH,DBGNDH,PBGNDH,UBGNDH,ADMSH,ECLSNH,FCLSNH,QHYD,RGASH,
     ;           CELLNW(7,:),CELLNW(5,:),CELLNW(6,:),CELLNW(8,:))
            
C     ALEX Check to see if the planet is earth or saturn, then do the update for He+ or H2+
            if (NamePlanet .eq. 'Earth ') then
               CALL Solver('He3',uHel,dHel,pHel,uSurHe,dSurHe
     ;              ,pSurHe,DBGNHe,PBGNHe,UBGNHe,ADMSHe,ECLSHe,FCLSHe,qHel,RGASHe,
     ;              CELLNW(11,:),CELLNW(9,:),CELLNW(10,:),CELLNW(12,:))
            endif
            if (NamePlanet .eq. 'Saturn') then
CALEX            do K=1,NDIM
CALEX               CELLNW(11,K)=XMSHE*jp2/kc1
CALEX               CELLNW(9,K)=0.
CALEX               CELLNW(10,K)=RGASHE*CELLNW(11,K)*800.
CALEX               CELLNW(12,K)=800.
CALEX            enddo
               call clfmhe
            endif
         endif
         
         if (TypeSolver == 'Rusanov') then
            
            call rusanov_solver('O3',RgasO,
     &           dOxyg, uOxyg, pOxyg,tOxyg, 
     &           dSurfO, uSurfO, pSurfO, 
     &           dBgndO, uBgndO, pBgndO, 
     &           dBgndO2, uBgndO2, pBgndO2, 
     &           AdmsO, FCLSNO, ECLSNO, 
     &           CELLNW(3,:), CELLNW(1,:), CELLNW(2,:),CELLNW(4,:))

            call rusanov_solver('H3',RgasH,
     &           dHYD, uHYD, pHYD, tHyd,
     &           dSurfH, uSurfH, pSurfH, 
     &           dBgndH, uBgndH, pBgndH, 
     &           dBgndH2, uBgndH2, pBgndH2, 
     &           AdmsH, FCLSNH, ECLSNH, 
     &           CELLNW(7,:), CELLNW(5,:), CELLNW(6,:),CELLNW(8,:))

            if (NamePlanet == 'Earth ') then
               call rusanov_solver('He3',RgasHe,
     &              dHel, uHel, pHel,tHel, 
     &              dSurHe, uSurHe, pSurHe, 
     &              dBgnHe, uBgnHe, pBgnHe, 
     &              dBgnHe2, uBgnHe2, pBgnHe2, 
     &              AdmsHe, FCLSHe, ECLSHe, 
     &              CELLNW(11,:), CELLNW(9,:), CELLNW(10,:),CELLNW(12,:))
            endif
            if (NamePlanet == 'Saturn') then
               call clfmhe
            endif
         endif

!         pcHYD = maxval(abs(PHYD-CELLNW(6,:))/CELLNW(6,:))
!         pcOXY = maxval(abs(POXYG-CELLNW(2,:))/CELLNW(2,:))
!         pc = max (pcHYD,pcOXY)


!         CURHLP=EXP(-0.5*(TIME-CURTIM0)**2/CURTIM**2)
         DO 7 K=1,NCL
            UOXYG(K)=CELLNW(1,K)
            POXYG(K)=CELLNW(2,K)
            DOXYG(K)=CELLNW(3,K)
            TOXYG(K)=CELLNW(4,K)
            UHYD(K)=CELLNW(5,K)
            PHYD(K)=CELLNW(6,K)
            DHYD(K)=CELLNW(7,K)
            THYD(K)=CELLNW(8,K)
            UHEL(K)=CELLNW(9,K)
            PHEL(K)=CELLNW(10,K)
            DHEL(K)=CELLNW(11,K)
            THEL(K)=CELLNW(12,K)
            DELECT(K)=RTHDEL*DHYD(K)+RTOXEL*DOXYG(K)+RTHEEL*DHEL(K)
            UELECT(K)=(RTHDEL*DHYD(K)*UHYD(K)+RTOXEL*DOXYG(K)*
     $           UOXYG(K)+RTHEEL*DHEL(K)*UHEL(K)
     $           -1.8965E-18*CURR(K))/DELECT(K)
 7       CONTINUE
!         call calc_cfl
!         if (TypeSolver == 'Godunov') then
!            call diffusion(DHYD,DSURFH,DBGNDH,cMax_H)
!            call diffusion(DOXYG,DSURFO,DBGNDO,cMax_O)
!            
!            call diffusion(UHYD,USURFH,UBGNDH,cMax_H)
!            call diffusion(UOXYG,USURFO,UBGNDO,cMax_O)
!            
!            call diffusion(PHYD,PSURFH,PBGNDH,cMax_H)
!            call diffusion(POXYG,PSURFO,PBGNDO,cMax_O)
!         endif
         do k=1,nDim
            DELECT(K)=RTHDEL*DHYD(K)+RTOXEL*DOXYG(K)+RTHEEL*DHEL(K)
            UELECT(K)=(RTHDEL*DHYD(K)*UHYD(K)+RTOXEL*DOXYG(K)*
     $           UOXYG(K)+RTHEEL*DHEL(K)*UHEL(K)
     $           -1.8965E-18*CURR(K))/DELECT(K)
         enddo

          

         CALL NEWBGD
         CALL ELFLDW


      ELSE
         CALL CLTEMPW
         CALL CLFME2W
CALEX call the solver
         if (TypeSolver == 'Godunov') then
            CALL Solver('O4',UOXYG,DOXYG,POXYG,USURFO,DSURFO
     ;           ,PSURFO,DBGNDO,PBGNDO,UBGNDO,ADMSO,ECLSNO,FCLSNO,QOXYG,RGASO,
     ;           CELLNW(3,:),CELLNW(1,:),CELLNW(2,:),CELLNW(4,:))
            
            CALL Solver('H4',UHYD,DHYD,PHYD,USURFH,DSURFH
     ;           ,PSURFH,DBGNDH,PBGNDH,UBGNDH,ADMSH,ECLSNH,FCLSNH,QHYD,RGASH,
     ;           CELLNW(7,:),CELLNW(5,:),CELLNW(6,:),CELLNW(8,:))
!            CURHLP=EXP(-0.5*(TIME-CURTIM0)**2/CURTIM**2)
            
C     ALEX Check to see if the planet is earth or saturn, then do the update for He+ or H2+
         
            if (NamePlanet .eq. 'Earth ') then
               CALL Solver('He4',uHel,dHel,pHel,uSurHe,dSurHe
     ;              ,pSurHe,DBGNHe,PBGNHe,UBGNHe,ADMSHe,ECLSHe,FCLSHe,qHel,RGASHe,
     ;              CELLNW(11,:),CELLNW(9,:),CELLNW(10,:),CELLNW(12,:))
            endif
            
            if (NamePlanet .eq. 'Saturn') then
CALEX            do K=1,NDIM
CALEX               CELLNW(11,K)=XMSHE*jp2/kc1
CALEX               CELLNW(9,K)=0.
CALEX               CELLNW(10,K)=RGASHE*CELLNW(11,K)*800.
CALEX               CELLNW(12,K)=800.
CALEX            enddo
               call clfmhe
            endif
         endif

         if (TypeSolver == 'Rusanov') then
            
            call rusanov_solver('O4',RgasO,
     &           dOxyg, uOxyg, pOxyg,tOxyg, 
     &           dSurfO, uSurfO, pSurfO, 
     &           dBgndO, uBgndO, pBgndO, 
     &           dBgndO2, uBgndO2, pBgndO2, 
     &           AdmsO, FCLSNO, ECLSNO, 
     &           CELLNW(3,:), CELLNW(1,:), CELLNW(2,:),CELLNW(4,:))

            call rusanov_solver('H4',RgasH,
     &           dHYD, uHYD, pHYD, tHyd,
     &           dSurfH, uSurfH, pSurfH, 
     &           dBgndH, uBgndH, pBgndH, 
     &           dBgndH2, uBgndH2, pBgndH2, 
     &           AdmsH, FCLSNH, ECLSNH, 
     &           CELLNW(7,:), CELLNW(5,:), CELLNW(6,:),CELLNW(8,:))

            if (NamePlanet == 'Earth ') then
               call rusanov_solver('He4',RgasHe,
     &              dHel, uHel, pHel, tHel,
     &              dSurHe, uSurHe, pSurHe, 
     &              dBgnHe, uBgnHe, pBgnHe, 
     &              dBgnHe2, uBgnHe2, pBgnHe2, 
     &              AdmsHe, FCLSHe, ECLSHe, 
     &              CELLNW(11,:), CELLNW(9,:), CELLNW(10,:),CELLNW(12,:))
            endif
            if (NamePlanet == 'Saturn') then
               call clfmhe
            endif
         endif
!         pcHYD = maxval(abs(PHYD-CELLNW(6,:))/CELLNW(6,:))
!         pcOXY = maxval(abs(POXYG-CELLNW(2,:))/CELLNW(2,:))
!         pc = max (pcHYD,pcOXY)

         DO 4 K=1,NDIM
            UOXYG(K)=CELLNW(1,K)
            POXYG(K)=CELLNW(2,K)
            DOXYG(K)=CELLNW(3,K)
            TOXYG(K)=CELLNW(4,K)
            UHYD(K)=CELLNW(5,K)
            PHYD(K)=CELLNW(6,K)
            DHYD(K)=CELLNW(7,K)
            THYD(K)=CELLNW(8,K)
            UHEL(K)=CELLNW(9,K)
            PHEL(K)=CELLNW(10,K)
            DHEL(K)=CELLNW(11,K)
            THEL(K)=CELLNW(12,K)
            DELECT(K)=RTHDEL*DHYD(K)+RTOXEL*DOXYG(K)+RTHEEL*DHEL(K)
            UELECT(K)=(RTHDEL*DHYD(K)*UHYD(K)+RTOXEL*DOXYG(K)*
     $           UOXYG(K)+RTHEEL*DHEL(K)*UHEL(K)
     $           -1.8965E-18*CURR(K))/DELECT(K)
 4       CONTINUE
!         call calc_cfl
!         if (TypeSolver == 'Godunov') then
!            call diffusion(DHYD,DSURFH,DBGNDH,cMax_H)
!            call diffusion(DOXYG,DSURFO,DBGNDO,cMax_O)
!            
!            call diffusion(UHYD,USURFH,UBGNDH,cMax_H)
!            call diffusion(UOXYG,USURFO,UBGNDO,cMax_O)
!            
!            call diffusion(PHYD,PSURFH,PBGNDH,cMax_H)
!            call diffusion(POXYG,PSURFO,PBGNDO,cMax_O)
!         endif
         do k=1,nDim
            DELECT(K)=RTHDEL*DHYD(K)+RTOXEL*DOXYG(K)+RTHEEL*DHEL(K)
            UELECT(K)=(RTHDEL*DHYD(K)*UHYD(K)+RTOXEL*DOXYG(K)*
     $           UOXYG(K)+RTHEEL*DHEL(K)*UHEL(K)
     $           -1.8965E-18*CURR(K))/DELECT(K)
         enddo

         CALL NEWBGD
         CALL COLLIS(NDIM)
         CALL ELFLDW

      ENDIF
      TIME=TIME+DT
!      call calc_cfl
c      call calcdt(pc)


c      write(*,*) DT

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC               
C     This determines whether or not to output a restart file.
C     There are two methods. Either output on iteration for a 
C     Fixed timestep, or output on frequency for a variable DT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF (.not. IsVariableDt) then
         IF (NSTEP.EQ.NOTP(KSTEP)) Then
            KSTEP=KSTEP+1
c            CALL PRNTOP
c            CALL PRNTHE
c            CALL PRNTHP
c            CALL PRNTEL
c            CALL PRNTEF
            write(iUnitOutput,*) 'time:',time
            CALL PRNTSM
            !CALL prntCollision
            !CALL PRNT_Sources
            
            open(UnitTmp_, FILE=NameRestart)
            WRITE (UnitTmp_,2001) TIME,DT,NSTEP
            Write (UnitTmp_,*) GMLAT, GMLONG

            WRITE (UnitTmp_,2002)(RAD(K),UOXYG(K),POXYG(K),DOXYG(K),TOXYG(K),K=1,NDIM)
            WRITE (UnitTmp_,2002)(RAD(K),UHEL(K),PHEL(K),DHEL(K),THEL(K),K=1,NDIM)
            WRITE (UnitTmp_,2002)(RAD(K),UHYD(K),PHYD(K),DHYD(K),THYD(K),K=1,NDIM)
            WRITE (UnitTmp_,2002)(RAD(K),UELECT(K),PELECT(K),DELECT(K),TELECT(K),
     ;           K=1,NDIM)
            close(UnitTmp_)
         endif
      endif
      IF (IsVariableDt) then
         IF (floor(Time/DToutput).NE.floor((Time-DT)/DToutput) ) Then
            KSTEP=KSTEP+1
c            CALL PRNTOP
c            CALL PRNTHE
c            CALL PRNTHP
c            CALL PRNTEL
c            CALL PRNTEF
            write(iUnitOutput,*) 'time:',time
            CALL PRNTSM
            !CALL prntCollision
            !CALL PRNT_Sources
            
            open(UnitTmp_, FILE=NameRestart)
            WRITE (UnitTmp_,2001) TIME,DT,NSTEP
            Write (UnitTmp_,*) GMLAT, GMLONG

            WRITE (UnitTmp_,2002)(RAD(K),UOXYG(K),POXYG(K),DOXYG(K),TOXYG(K),K=1,NDIM)
            WRITE (UnitTmp_,2002)(RAD(K),UHEL(K),PHEL(K),DHEL(K),THEL(K),K=1,NDIM)
            WRITE (UnitTmp_,2002)(RAD(K),UHYD(K),PHYD(K),DHYD(K),THYD(K),K=1,NDIM)
            WRITE (UnitTmp_,2002)(RAD(K),UELECT(K),PELECT(K),DELECT(K),TELECT(K),
     ;           K=1,NDIM)
            close(UnitTmp_)
         endif
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF (TIME.GT.TMAX) Then 
         Call put_field_line(dOxyg, uOxyg, pOxyg, TOxyg,     
     &                       dHel, uHel, pHel, THel,         
     &                       dHyd, uHyd, pHyd, THyd,         
     &                       dElect, uElect, pElect, TElect, 
     &                       GMLAT,GMLONG,Jr,wHorizontal,
     &                       Time=Time ) 

         RETURN
      endif
      
      NSTEP=NSTEP+1
      IF (MSTEP .EQ. NTS) MSTEP=0 
      MSTEP=MSTEP+1

      if (.not.IsVariableDT) then
         IF (KSTEP.GT.NPT) RETURN
      endif

      GO TO 1
      
      END
