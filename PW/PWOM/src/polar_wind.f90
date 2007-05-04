
Subroutine POLAR_WIND
  ! Discussion:
  ! This subroutine solves the polar wind solution along a given field line
  ! called by advance_line. A. Glocer
  !
  ! photoelectron impact ionization modified for no conjugate hemisphere
  ! P. E. sources. See Subroutine Glowex, variable 'rat' D. Cannata 9May90
  ! 
  ! ADDED SUBROUTINE CALL FOR IMPACT IONIZATION "call precip" PRESENTLY
  ! HARD WIRED FOR POLAR RAIN FOR INVARIENT LATITUDES GREATER THAN 75 DEG.
  ! CHANGES FOR CUSP OR AURORAL PRECIPITATION CAN BE MADE BY CHANGING
  ! ALF1(x) AND FLUX1(x) ARRAY, x=1->4. MODIFIFIED BY DICK AND STEVE 6-90.
  !
  
  Use ModPWOM, only: DtVertical,nLine,IsStandAlone,DoSavePlot,iLine
  use ModIoUnit, ONLY: UnitTmp_
  use ModCommonVariables
  use ModFieldLine
  INTEGER NOTP(100)
  
  !     define the output files and attaching units
  character*100 :: NameRestart
  Logical,save  :: IsFirstCall=.true.
  Real Jr
  !-----------------------------------------------------------------------
  
  Call get_field_line(dOxyg(:), uOxyg(:), pOxyg(:), TOxyg(:),      &    
       dHel(:), uHel(:), pHel(:), THel(:),                         &
       dHyd(:), uHyd(:), pHyd(:), THyd(:),                         &
       dElect(:), uElect(:), pElect(:), TElect(:),                 &
       GMLAT,GMLONG,Jr,wHorizontal,                                &
       iUnitOutput=iUnitOutput, iUnitGraphics=iUnitGraphics,       &
       NameRestart=NameRestart,                                    &
       iLine=iLine, Time=Time,MaxLineTime=Tmax,                    &
       TypeSolver=TypeSolver,IsVariableDT=IsVariableDT,            &
       IsRestart=IsRestart,DToutput=DToutput,nAlt=nDim,DoLog=DoLog,&
       nStep=nStep,IsImplicit=IsImplicit,IsImplicitAll=IsImplicitAll)
  
  DT=DtVertical
  
  
  
  
  CURR(1)      = Jr
  
  NTS = 1
  NCL = 60
  
  
21 FORMAT(2X,I8)
  
  KSTEP=1
  MSTEP=1
  
  NDIM2=NDIM-1
  NDIM1=NDIM+1
  NDIMM=NDIM+2
  
  CALL STRT
  
  if (IsFirstCall .and. DoSavePlot) then
     CALL PW_print_plot
     if(iLine == nLine) IsFirstCall = .false.
  endif
  
  !******************************************************************************
  ! This is the start of a while loop that carries out the main steps of
  ! the simulation
  !******************************************************************************
  
  TIMELOOP: DO
     !     If (UseFullImplicit) then
     !        call implicit_solver
     !     else
     !       advect and add collisional and chemistry sources
     !       calculate heatflux for ions and update ion temp.
     !       calculate heatflux for electrons, update electron temp.
     !       update boundaries
     !       calculate collisional and chemistry source
     !       calculate electric field with new electron temp.
     
     call advect
     CALL CLTEMPW
     CALL CLFME2W
     CALL PW_set_upper_bc
     CALL COLLIS(NDIM)
     CALL ELFLDW         
     
     TIME=TIME+DT
     NSTEP=NSTEP+1
     
     IF (floor((Time+1.0e-5)/DToutput)/=floor((Time+1.0e-5-DT)/DToutput) )& 
          CALL PW_print_plot
     
     !       Reverse order of advection and implicit temperature update
     
     CALL CLTEMPW
     CALL CLFME2W
     CALL advect
     !     endif
     
     !    finish update by calculating boundaries, collision source and electric field
     !    these will be used in the next time step
     CALL PW_set_upper_bc
     CALL COLLIS(NDIM)
     CALL ELFLDW
     
     TIME=TIME+DT
     NSTEP=NSTEP+1
     IF (floor((Time+1.0e-5)/DToutput)/=floor((Time+1.0e-5-DT)/DToutput) ) &
          CALL PW_print_plot
     
     IF (IsStandAlone .and. &
          floor((Time+1.0e-5)/DToutput)/=floor((Time+1.0e-5-2.0*DT)/DToutput) )&
          call PW_write_restart(&
          nDim,RAD,GmLat,GmLong,Time,DT,nStep,NameRestart, &    
          dOxyg, uOxyg, pOxyg, TOxyg,                      &
          dHel , uHel , pHel , THel ,                      &
          dHyd , uHyd , pHyd , THyd ,                      &
          dElect, uElect, pElect, TElect)                  
     
     
     IF (TIME+1.0e-5 >= TMAX) Then 
        Call put_field_line(dOxyg, uOxyg, pOxyg, TOxyg,     &
             dHel, uHel, pHel, THel,                        & 
             dHyd, uHyd, pHyd, THyd,                        &
             dElect, uElect, pElect, TElect,                &
             GMLAT,GMLONG,Jr,wHorizontal,                   &
             Time=Time,nStep=nStep,r_C=RAD   ) 
        
        RETURN
     endif
     
  enddo TIMELOOP
  
contains      
  !==============================================================================
  
  subroutine advect
    
    if (TypeSolver == 'Godunov') then
       CALL Solver('O1',UOXYG,DOXYG,POXYG,USURFO,DSURFO,                 &
             PSURFO,DBGNDO,PBGNDO,UBGNDO,ADMSO,ECLSNO,FCLSNO,QOXYG,RGASO,&
             CELLNW(3,:),CELLNW(1,:),CELLNW(2,:),CELLNW(4,:))

       CALL Solver('H1',UHYD,DHYD,PHYD,USURFH,DSURFH,                  &
            PSURFH,DBGNDH,PBGNDH,UBGNDH,ADMSH,ECLSNH,FCLSNH,QHYD,RGASH,&
            CELLNW(7,:),CELLNW(5,:),CELLNW(6,:),CELLNW(8,:))

         if (NamePlanet == 'Earth ')&
              CALL Solver('He1',uHel,dHel,pHel,uSurHe,dSurHe,              &
              pSurHe,DBGNHe,PBGNHe,UBGNHe,ADMSHe,ECLSHe,FCLSHe,qHel,RGASHe,&
              CELLNW(11,:),CELLNW(9,:),CELLNW(10,:),CELLNW(12,:))

      else if (TypeSolver == 'Rusanov') then
            
         call rusanov_solver(Ion1_,RgasO,&
              dOxyg, uOxyg, pOxyg, tOxyg,& 
              dSurfO, uSurfO, pSurfO,    &
              dBgndO, uBgndO, pBgndO,    &
              dBgndO2, uBgndO2, pBgndO2, &
              AdmsO, FCLSNO, ECLSNO,     &
              CELLNW(3,:), CELLNW(1,:), CELLNW(2,:),CELLNW(4,:))
         
         call rusanov_solver(Ion2_,RgasH,&
             dHYD, uHYD, pHYD,tHyd,      &
             dSurfH, uSurfH, pSurfH,     &
             dBgndH, uBgndH, pBgndH,     &
             dBgndH2, uBgndH2, pBgndH2,  &
             AdmsH, FCLSNH, ECLSNH,      & 
             CELLNW(7,:), CELLNW(5,:), CELLNW(6,:),CELLNW(8,:))
         
         if (NamePlanet == 'Earth ')           &
              call rusanov_solver(Ion3_,RgasHe,&
              dHel, uHel, pHel,tHel,           &
              dSurHe, uSurHe, pSurHe,          &
              dBgnHe, uBgnHe, pBgnHe,          &
              dBgnHe2, uBgnHe2, pBgnHe2,       & 
              AdmsHe, FCLSHe, ECLSHe,          &  
              CELLNW(11,:), CELLNW(9,:), CELLNW(10,:),CELLNW(12,:))
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
      enddo

      do k=1,nDim
         DELECT(K)=RTHDEL*DHYD(K)+RTOXEL*DOXYG(K)+RTHEEL*DHEL(K)
         UELECT(K)=(RTHDEL*DHYD(K)*UHYD(K)+RTOXEL*DOXYG(K)* &
              UOXYG(K)+RTHEEL*DHEL(K)*UHEL(K) &
              -1.8965E-18*CURR(K))/DELECT(K)
      enddo
      
    end subroutine advect
    
  end Subroutine POLAR_WIND

