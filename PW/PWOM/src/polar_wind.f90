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

      Use ModPWOM, only: DtVertical,nLine,IsStandAlone,DoSavePlot
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
     &                 IsRestart=IsRestart,DToutput=DToutput,nAlt=nDim,DoLog=DoLog,
     &                 nStep=nStep,IsImplicit=IsImplicit,IsImplicitAll=IsImplicitAll)
      
      DT=DtVertical

      


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
!      NSTEP=1
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
!      if (NamePlanet == 'Earth ') nDim=390
C Use nDim=550 for top boundary at 2 Re in the earth case
c      nDim=550
c      nDim=400
!      if (NamePlanet == 'Saturn') nDim=800
c      if (NamePlanet == 'Saturn') nDim=930

c      nDim=680
c      nDim=1500
c      nDim=1180
c      nDim=750

      
      NDIM2=NDIM-1
      NDIM1=NDIM+1
      NDIMM=NDIM+2
      
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     DEFINE THE INITIAL CONDITIONS IN THE CELLS AND AT CELL           C
C           BOUNDARIES                                                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C

         

      CALL STRT

      



      if (IsFirstCall .and. DoSavePlot) then
         CALL PW_print_plot
         IsFirstCall = .false.
      endif

!******************************************************************************
! This is the start of a while loop that carries out the main steps of
! the simulation
!******************************************************************************

      TIMELOOP: DO
C      If (UseFullImplicit) then
C         call implicit_solver
C      else
C        advect and add collisional and chemistry sources
C        calculate heatflux for ions and update ion temp.
C        calculate heatflux for electrons, update electron temp.
C        update boundaries
C        calculate collisional and chemistry source
C        calculate electric field with new electron temp.

         call advect
         CALL CLTEMPW
         CALL CLFME2W
         CALL NEWBGD
         CALL COLLIS(NDIM)
         CALL ELFLDW         

         TIME=TIME+DT
         NSTEP=NSTEP+1

         IF (floor((Time+1.0e-5)/DToutput)/=floor((Time+1.0e-5-DT)/DToutput) ) 
     &        CALL PW_print_plot

C        Reverse order of advection and implicit temperature update

         CALL CLTEMPW
         CALL CLFME2W
         CALL advect
C      endif

C     finish update by calculating boundaries, collision source and electric field
C     these will be used in the next time step
      CALL NEWBGD
      CALL COLLIS(NDIM)
      CALL ELFLDW

      TIME=TIME+DT
      NSTEP=NSTEP+1
      IF (floor((Time+1.0e-5)/DToutput)/=floor((Time+1.0e-5-DT)/DToutput) ) 
     &     CALL PW_print_plot

      IF (IsStandAlone .and. 
     &     floor((Time+1.0e-5)/DToutput)/=floor((Time+1.0e-5-2.0*DT)/DToutput) )
     &     call PW_write_restart(
     &     nDim,RAD,GmLat,GmLong,Time,DT,nStep,NameRestart,     
     &     dOxyg, uOxyg, pOxyg, TOxyg,         
     &     dHel , uHel , pHel , THel ,         
     &     dHyd , uHyd , pHyd , THyd ,         
     &     dElect, uElect, pElect, TElect)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      

      IF (TIME+1.0e-5 >= TMAX) Then 
         Call put_field_line(dOxyg, uOxyg, pOxyg, TOxyg,     
     &                       dHel, uHel, pHel, THel,         
     &                       dHyd, uHyd, pHyd, THyd,         
     &                       dElect, uElect, pElect, TElect, 
     &                       GMLAT,GMLONG,Jr,wHorizontal,
     &                       Time=Time,nStep=nStep,r_C=RAD ) 

         RETURN
      endif
      
      enddo TIMELOOP

      contains      
C==============================================================================

      subroutine advect

      if (TypeSolver == 'Godunov') then
         CALL Solver('O1',UOXYG,DOXYG,POXYG,USURFO,DSURFO,
     ;        PSURFO,DBGNDO,PBGNDO,UBGNDO,ADMSO,ECLSNO,FCLSNO,QOXYG,RGASO,
     ;        CELLNW(3,:),CELLNW(1,:),CELLNW(2,:),CELLNW(4,:))

         CALL Solver('H1',UHYD,DHYD,PHYD,USURFH,DSURFH
     ;        ,PSURFH,DBGNDH,PBGNDH,UBGNDH,ADMSH,ECLSNH,FCLSNH,QHYD,RGASH,
     ;        CELLNW(7,:),CELLNW(5,:),CELLNW(6,:),CELLNW(8,:))

         if (NamePlanet == 'Earth ')
     &        CALL Solver('He1',uHel,dHel,pHel,uSurHe,dSurHe
     ;        ,pSurHe,DBGNHe,PBGNHe,UBGNHe,ADMSHe,ECLSHe,FCLSHe,qHel,RGASHe,
     ;        CELLNW(11,:),CELLNW(9,:),CELLNW(10,:),CELLNW(12,:))

      else if (TypeSolver == 'Rusanov') then
            
         call rusanov_solver(Ion1_,RgasO,
     &        dOxyg, uOxyg, pOxyg, tOxyg, 
     &        dSurfO, uSurfO, pSurfO, 
     &        dBgndO, uBgndO, pBgndO, 
     &        dBgndO2, uBgndO2, pBgndO2, 
     &        AdmsO, FCLSNO, ECLSNO, 
     &        CELLNW(3,:), CELLNW(1,:), CELLNW(2,:),CELLNW(4,:))

         call rusanov_solver(Ion2_,RgasH,
     &        dHYD, uHYD, pHYD,tHyd, 
     &        dSurfH, uSurfH, pSurfH, 
     &        dBgndH, uBgndH, pBgndH, 
     &        dBgndH2, uBgndH2, pBgndH2, 
     &        AdmsH, FCLSNH, ECLSNH, 
     &        CELLNW(7,:), CELLNW(5,:), CELLNW(6,:),CELLNW(8,:))
         
         if (NamePlanet == 'Earth ')call rusanov_solver(Ion3_,RgasHe,
     &        dHel, uHel, pHel,tHel, 
     &        dSurHe, uSurHe, pSurHe, 
     &        dBgnHe, uBgnHe, pBgnHe, 
     &        dBgnHe2, uBgnHe2, pBgnHe2, 
     &        AdmsHe, FCLSHe, ECLSHe, 
     &        CELLNW(11,:), CELLNW(9,:), CELLNW(10,:),CELLNW(12,:))
      endif

      if (NamePlanet == 'Saturn') call clfmhe
      
      DO K=1,NDIM
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
c     write(*,*)'RTHDEL*DHYD(K),RTOXEL*DOXYG(K),RTHEEL*DHEL(K)',  RTHDEL*DHYD(K),RTOXEL*DOXYG(K),RTHEEL*DHEL(K)
      enddo

      do k=1,nDim
         DELECT(K)=RTHDEL*DHYD(K)+RTOXEL*DOXYG(K)+RTHEEL*DHEL(K)
         UELECT(K)=(RTHDEL*DHYD(K)*UHYD(K)+RTOXEL*DOXYG(K)*
     $        UOXYG(K)+RTHEEL*DHEL(K)*UHEL(K)
     $        -1.8965E-18*CURR(K))/DELECT(K)
      enddo

      end subroutine advect

      END

