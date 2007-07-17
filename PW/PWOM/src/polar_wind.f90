subroutine polar_wind

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
  
  use ModPWOM, only: DtVertical,nLine,IsStandAlone,DoSavePlot,iLine,IsFullyImplicit
  use ModIoUnit, ONLY: UnitTmp_
  use ModCommonVariables
  use ModFieldLine
  INTEGER NOTP(100)
  
  !     define the output files and attaching units
  character*100 :: NameRestart
  logical,save  :: IsFirstCall=.true.
  real    :: Jr
  real    :: NewState_GV(-1:maxGrid,nVar)

  !-----------------------------------------------------------------------
  
  call get_field_line(State_GV(1:maxGrid,:),                       &
       GMLAT,GMLONG,Jr,wHorizontal,                                &
       iUnitOutput=iUnitOutput, iUnitGraphics=iUnitGraphics,       &
       NameRestart=NameRestart,                                    &
       iLine=iLine, Time=Time,MaxLineTime=Tmax,                    &
       TypeSolver=TypeSolver,IsVariableDT=IsVariableDT,            &
       IsRestart=IsRestart,DToutput=DToutput,nAlt=nDim,DoLog=DoLog,&
       nStep=nStep)
  
  DT=DtVertical
  DtImpl = Dt
  DtExpl = Dt*0.1

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
     If (IsFullyImplicit) then
        call PW_implicit_update
     else
        !       advect and add collisional and chemistry sources
        !       calculate heatflux for ions and update ion temp.
        !       calculate heatflux for electrons, update electron temp.
        !       update boundaries
        !       calculate collisional and chemistry source
        !       calculate electric field with new electron temp.
        
        
        call advect
     endif
     
     CALL PW_iheat_flux
     
     CALL PW_eheat_flux
     
     CALL PW_set_upper_bc
     CALL COLLIS(NDIM)
     CALL PW_calc_efield         
     TIME=TIME+DT


     

     NSTEP=NSTEP+1
     
     IF (floor((Time+1.0e-5)/DToutput)/=floor((Time+1.0e-5-DT)/DToutput) )& 
          CALL PW_print_plot
     
     !       Reverse order of advection and implicit temperature update
     
     CALL PW_iheat_flux

     CALL PW_eheat_flux
     If (IsFullyImplicit) then
        call PW_implicit_update
     else
        CALL advect
     endif

     !    finish update by calculating boundaries, collision source 
     !    and electric field
     !    these will be used in the next time step

     CALL PW_set_upper_bc
     CALL COLLIS(NDIM)
     CALL PW_calc_efield
     
     TIME=TIME+DT
     NSTEP=NSTEP+1


     IF (floor((Time+1.0e-5)/DToutput)/=floor((Time+1.0e-5-DT)/DToutput) ) &
          CALL PW_print_plot
     
     IF (IsStandAlone .and. &
          floor((Time+1.0e-5)/DToutput)/=floor((Time+1.0e-5-2.0*DT)/DToutput) )&
          call PW_write_restart(&
          nDim,RAD,GmLat,GmLong,Time,DT,nStep,NameRestart, &    
          State_GV)                  
       
     IF (TIME+1.0e-5 >= TMAX) Then 
        Call put_field_line(State_GV(1:maxGrid,:),    &
             GMLAT,GMLONG,Jr,wHorizontal,&
             Time=Time,nStep=nStep,r_C=RAD   ) 
        
        RETURN
     endif
     
  enddo TIMELOOP
  
contains      
  !============================================================================
  
  subroutine advect
    
    NewState_GV = State_GV
    if (TypeSolver == 'Godunov') then
       do iIon=1,nIon-1
          CALL Solver(iIon,Dt,State_GV(:,iRho_I(iIon):iT_I(iIon)),&
               Source_CV(:,iRho_I(iIon)),Source_CV(:,iP_I(iIon)),&
               Source_CV(:,iU_I(iIon)),&
               RGAS_I(iIon),NewState_GV(:,iRho_I(iIon):iT_I(iIon)))
       enddo

      else if (TypeSolver == 'Rusanov') then
         do iIon=1,nIon-1
            call rusanov_solver(iIon,RGAS_I(iIon),dt,   &
                 State_GV(:,iRho_I(iIon):iT_I(iIon)),&
                 Source_CV(:,iRho_I(iIon)), Source_CV(:,iU_I(iIon)),&
                 Source_CV(:,iP_I(iIon)),  &
                 NewState_GV(:,iRho_I(iIon):iT_I(iIon)))
         enddo
      endif
            
      State_GV(:,:) = NewState_GV(:,:)

      ! Set electron density and velocity
      State_GV(1:nDim,RhoE_)=0.0
      State_GV(1:nDim,uE_)  =0.0
      do k=1,nDim
         do iIon=1,nIon-1
            State_GV(K,RhoE_) = &
                 State_GV(k,RhoE_)+MassElecIon_I(iIon)*State_GV(K,iRho_I(iIon))
            State_GV(K,uE_)= &
                 State_GV(k,uE_)+ &
                 (MassElecIon_I(iIon)*State_GV(K,iRho_I(iIon))&
                  *State_GV(K,iU_I(iIon)))
         enddo
         
         State_GV(K,uE_)=(State_GV(K,uE_) -1.8965E-18*CURR(K))/State_GV(K,RhoE_)
      enddo
    end subroutine advect
    
  end Subroutine POLAR_WIND

