!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine PW_initialize
  
  use ModPlanetConst, ONLY: Planet_, NamePlanet_I
  use ModNumConst, ONLY: cDegToRad
  use ModMpi
  use ModIoUnit, ONLY: io_unit_new,UnitTmp_
  use ModPwom
  use ModCommonPlanet,ONLY: nIon,iRho_I,iU_I,iP_I,iT_I
  use ModCommonVariables, ONLY:IYD,ALTD,Mass_I,DRBND,MassElecIon_I,Rgas_I,XAMU,&
       RGAS,Ion1_,Ion2_,Ion3_,nIon
  use ModTimeConvert, ONLY: time_int_to_real
  use ModPwTime
  use ModAurora, ONLY: init_aurora
  use ModPwWaves,ONLY: wave_init
  use ModCouplePWOMtoSE, ONLY: init_pwom_se_coupling
  use ModPhotoelectron, ONLY: PrecipEnergyMin, PrecipEnergyMax, &
       PrecipEnergyMean, PrecipEnergyFlux, UseFixedPrecip, DoCoupleSE,&
       PolarRainEMin, PolarRainEMax, &
       PolarRainEMean, PolarRainEFlux, UsePolarRain, IsVerboseSE
  use ModOvation, ONLY: UseOvation, StartTimeOvation=>StartTime, &
       OvationEmin,OvationEmax
  use ModParticle, ONLY: init_particle, put_to_particles, bury_line,&
       nLineParticle=>nLine,iLineGlobalParticle_I=>iLineGlobal_I, &
       read_restart_particle
  use CON_axes,         ONLY: init_axes
  implicit none

  ! Temporary variables
  real:: ddt1, xxx
  integer:: ns, iPe, iError,iIon
  integer:: iYear, iDOY

  integer, external :: julianday

  ! AMU in grams
  real, parameter :: AMUinGrams=1.6606655E-24
  
  !for particles
  real,allocatable :: Density_IC(:,:),Velocity_IC(:,:),Temperature_IC(:,:)
  
  !---------------------------------------------------------------------------
  !***************************************************************************
  !  Set the number of fieldlines that each processor solves for
  !***************************************************************************
  
  if (nTotalLine < nProc) &
       call con_stop(&
       "PW ERROR:nTotalLine<nProc. Reduce number of procs for PW in LAYOUT.in")

  if (iProc < mod(nTotalLine,nProc)) then
     nLine= (nTotalLine+nProc-1)/nProc
  else
     nLine= nTotalLine/nProc
  endif

  allocate(nLine_P(0:nProc-1), nLineBefore_P(0:nProc-1))
  call MPI_allgather(nLine,1,MPI_INTEGER, nLine_P,1,MPI_INTEGER,iComm,iError)
  nLineBefore_P(0) = 0
  do iPe = 1, nProc - 1
     nLineBefore_P(iPe) = sum(nLine_P(0:iPe-1))
  end do

  !\
  ! Set the Time parameters
  !/
  iYear = iStartTime(1)
  iDOY  = julianday(iStartTime(1),iStartTime(2),iStartTime(3)) 
  IYD=mod(iYear,100)*1000+iDOY
  call time_int_to_real(iStartTime,CurrentTime)
  StartTime=CurrentTime
  
  !make ovation starttime match simulation start time
  if (UseOvation) StartTimeOvation=StartTime

  !\
  ! Set axes for coord transform when in standalone mode
  !/
  if (IsStandAlone) call init_axes(StartTime)

  !\
  ! Allocate arrays for simulation
  !/
  if (.not.allocated(r_C)) allocate(r_C(nAlt),&
       State_CVI(nAlt,nVar,nLine),&
       GeoMagLat_I(nLine),GeoMagLon_I(nLine),          &
       ThetaLine_I(nLine), PhiLine_I(nLine),           &
       ThetaLineOld_I(nLine), PhiLineOld_I(nLine),           &
       xLine_I(nLine),yLine_I(nLine),zLine_I(nLine),          &
       xLineOld_I(nLine),yLineOld_I(nLine),zLineOld_I(nLine), &
       UthetaLine_I(nLine),UphiLine_I(nLine),          &
       UxLine_I(nLine),UyLine_I(nLine),UzLine_I(nLine),       &
       OmegaLine_I(nLine),                      &
       JrLine_I(nLine),EfluxLine_I(nLine),AvELine_I(nLine),&
       iThetaLine_I(nLine),iPhiLine_I(nLine), &
       NameRestartIn(nLine), NameRestart(nLine), NameGraphics(nLine),&
       NameOutput(nLine),  iUnitRestart(nLine),iUnitRestartIn(nLine),&
       iUnitGraphics(nLine),iUnitOutput(nLine), iLineGlobal(nLine), &
       Dt_I(nLine))
  
  call wave_init(nAlt)
  
  !**************************************************************************
  !  Define file names and unit numbers, and open for reading and writing.
  !***************************************************************************
  NameSourceGraphics = 'PW/plot_sources.out'
  NameCollision      = 'PW/plots_collision.out'
  NamePhiNorth       = 'PW/North.dat'
  NamePhiSouth       = 'PW/South.dat'

  do iLine=1,nLine
     if (iproc .lt. mod(nTotalLine,nProc)) then
        iLineGlobal(iLine)=&
             iproc*ceiling(real(nTotalLine)/real(nProc))+iLine
     else
        iLineGlobal(iLine)=&
             (mod(NTotalLine,nProc))*ceiling(real(nTotalLine)/real(nProc)) &
             + ((iproc)-mod(nTotalLine,nProc))                        &
             *floor(real(nTotalLine)/real(nProc))+iLine
     endif
     write(NameRestartIn(iLine),"(a,i4.4,a)") &
          'PW/restartIN/restart_iline',iLineGlobal(iLine),'.dat'
     write(NameRestart(iLine),"(a,i4.4,a)") &
          'PW/restartOUT/restart_iline',iLineGlobal(iLine),'.dat'
     
     !Setup log files
     if (nLog == -1) then
        write(NameOutput(iLine),"(a,i4.4,a)") &
             'PW/log_iline',iLineGlobal(iLine),'.out'   
        iUnitOutput(iLine)  = io_unit_new()
        open(iUnitOutput(iLine),FILE=NameOutput(iLine))    
     elseif(nLog ==0) then
        !do nothing in this case
     elseif(nLog==iLineGlobal(iLine)) then
        write(NameOutput(iLine),"(a,i4.4,a)") &
             'PW/log_iline',iLineGlobal(iLine),'.out'      
        iUnitOutput(iLine)  = io_unit_new()
        open(iUnitOutput(iLine),FILE=NameOutput(iLine))     
     else
     end if

     
     
  enddo

!******************************************************************************
!  Read the restart file
!******************************************************************************

  if(IsRestart)then
     call PW_read_restart
  else
     do iLine = 1, nLine
        ThetaLine_I (iLine) = 10.0 * cDegToRad
        PhiLine_I   (iLine) = 0.0
     end do
     Time=0.0
  endif
  
  CurrentTime=StartTime+Time
  
  ! Set the output plot files
  do iLine =1, nLine
     if(ThetaLine_I(iLine) <= 90 * cDegToRad)then
        write(NameGraphics(iLine),"(a,i4.4,a)") &
             'PW/plots/north_plots_iline',iLineGlobal(iLine),'.out'
     else
        write(NameGraphics(iLine),"(a,i4.4,a)") &
             'PW/plots/south_plots_iline',iLineGlobal(iLine),'.out'
     endif
     
     open(UnitTmp_,FILE=NameGraphics(iLine),STATUS='replace')
     close(UnitTmp_)
  enddo
  !****************************************************************************
  ! Set vertical field-line grid
  !****************************************************************************
  call set_vertical_grid

  !****************************************************************************
  ! Use Get_GITM to bring in neutral atmosphere from GITM
  !****************************************************************************
  !call GetNeutralData

  !****************************************************************************
  !  Set parameters for reading in potential and time of simulation
  !****************************************************************************
 
  TimeMax = Tmax

  !****************************************************************************
  ! Read information from IE file, and get the velocities
  !****************************************************************************
  
  if(.not.UseIe)then
     call PW_get_electrodynamics

     !initialize field line locations
     call initial_line_location
  end if


  if (UseAurora) call init_aurora

  ! initialize the SE model (note cannot use fix precip and ovation 
  ! simultaneously)
  if(DoCoupleSE) then
     if(UseFixedPrecip .and. UsePolarRain .and. .not.UseOvation) then
        call init_pwom_se_coupling(IsVerboseSE,nAlt,nLine,iLineGlobal,ALTD,&
             PrecipEminPwIn=PrecipEnergyMin,PrecipEmaxPwIn=PrecipEnergyMax, &
             PrecipEmeanPwIn=PrecipEnergyMean,PrecipEfluxPwIn=PrecipEnergyFlux,&
             PolarRainEminPwIn=PolarRainEMin,PolarRainEmaxPwIn=PolarRainEMax, &
             PolarRainEmeanPwIn=PolarRainEMean,&
             PolarRainEfluxPwIn=PolarRainEFlux)
     elseif(.not.UseFixedPrecip .and. UsePolarRain .and. UseOvation) then
        call init_pwom_se_coupling(IsVerboseSE,nAlt,nLine,iLineGlobal,ALTD,&
             PolarRainEminPwIn=PolarRainEMin,PolarRainEmaxPwIn=PolarRainEMax, &
             PolarRainEmeanPwIn=PolarRainEMean,&
             PolarRainEfluxPwIn=PolarRainEFlux,&
             OvationEminPwIn=OvationEmin,      &
             OvationEmaxPwIn=OvationEmax)
     elseif(UseFixedPrecip) then
        call init_pwom_se_coupling(IsVerboseSE,nAlt,nLine,iLineGlobal,ALTD,&
             PrecipEminPwIn=PrecipEnergyMin,PrecipEmaxPwIn=PrecipEnergyMax, &
             PrecipEmeanPwIn=PrecipEnergyMean,PrecipEfluxPwIn=PrecipEnergyFlux)
     elseif(UsePolarRain) then
        call init_pwom_se_coupling(IsVerboseSE,nAlt,nLine,iLineGlobal,ALTD,&
             PolarRainEminPwIn=PolarRainEMin,PolarRainEmaxPwIn=PolarRainEMax, &
             PolarRainEmeanPwIn=PolarRainEMean,&
             PolarRainEfluxPwIn=PolarRainEFlux)
     elseif(UseOvation) then
        call init_pwom_se_coupling(IsVerboseSE,nAlt,nLine,iLineGlobal,ALTD,&
             OvationEminPwIn=OvationEmin,      &
             OvationEmaxPwIn=OvationEmax)
     else
        call init_pwom_se_coupling(IsVerboseSE,nAlt,nLine,iLineGlobal,ALTD)
     end if
  end if
  

  !Set mass of species for each fluid/particle population
  select case(NamePlanet_I(Planet_))
  case('EARTH')
     ! Mass of atomic O in grams
     Mass_I(Ion1_)=15.994*AMUinGrams
     ! Mass of atomic H in grams
     Mass_I(Ion2_)=1.00797*AMUinGrams
     ! Mass of atomic He in grams
     Mass_I(Ion3_)=4.0026*AMUinGrams
     ! Mass of electron in grams
     Mass_I(nIon)=9.109534E-28

     ! Relative mass of atomic O to electron
     MassElecIon_I(Ion1_)=Mass_I(nIon)/Mass_I(Ion1_)
     ! Relative mass of atomic H to electron
     MassElecIon_I(Ion2_)=Mass_I(nIon)/Mass_I(Ion2_)
     ! Relative mass of atomic He to electron
     MassElecIon_I(Ion3_)=Mass_I(nIon)/Mass_I(Ion3_)
     ! kB/m_O
     RGAS_I(Ion1_)=RGAS*XAMU/Mass_I(Ion1_)
     ! kB/m_H
     RGAS_I(Ion2_)=RGAS*XAMU/Mass_I(Ion2_)
     ! kB/m_He
     RGAS_I(Ion3_)=RGAS*XAMU/Mass_I(Ion3_)
     ! kB/m_e
     RGAS_I(nIon)=RGAS*XAMU/Mass_I(nIon)
  case('JUPITER')
     ! Mass of atomic H3 in grams
      Mass_I(Ion1_)=3.0237*AMUinGrams
      ! Mass of atomic H in grams
      Mass_I(Ion2_)=1.00797*AMUinGrams
      ! Mass of H2 in grams
      Mass_I(Ion3_)=2.0159*AMUinGrams
      ! Mass of electron in grams
      Mass_I(nIon)=9.109534E-28
      ! Relative mass of H3 to electron
      MassElecIon_I(Ion1_)=Mass_I(nIon)/Mass_I(Ion1_)
      ! Relative mass of atomic H to electron
      MassElecIon_I(Ion2_)=Mass_I(nIon)/Mass_I(Ion2_)
      ! Relative mass of H2 to electron
      MassElecIon_I(Ion3_)=Mass_I(nIon)/Mass_I(Ion3_)
      ! kB/m_H3
      RGAS_I(Ion1_)=RGAS*XAMU/Mass_I(Ion1_)
      ! kB/m_H
      RGAS_I(Ion2_)=RGAS*XAMU/Mass_I(Ion2_)
      ! kB/m_H2
      RGAS_I(Ion3_)=RGAS*XAMU/Mass_I(Ion3_)
      ! kB/m_e
      RGAS_I(nIon)=RGAS*XAMU/Mass_I(nIon)
   case('SATURN')
      ! Mass of atomic H3 in grams
      Mass_I(Ion1_)=3.0237*AMUinGrams
      ! Mass of atomic H in grams
      Mass_I(Ion2_)=1.00797*AMUinGrams
      ! Mass of H2 in grams
      Mass_I(Ion3_)=2.0159*AMUinGrams
      ! Mass of electron in grams
      Mass_I(nIon)=9.109534E-28
      ! Relative mass of H3 to electron
      MassElecIon_I(Ion1_)=Mass_I(nIon)/Mass_I(Ion1_)
      ! Relative mass of atomic H to electron
      MassElecIon_I(Ion2_)=Mass_I(nIon)/Mass_I(Ion2_)
      ! Relative mass of H2 to electron
      MassElecIon_I(Ion3_)=Mass_I(nIon)/Mass_I(Ion3_)
      ! kB/m_H3
      RGAS_I(Ion1_)=RGAS*XAMU/Mass_I(Ion1_)
      ! kB/m_H
      RGAS_I(Ion2_)=RGAS*XAMU/Mass_I(Ion2_)
      ! kB/m_H2
      RGAS_I(Ion3_)=RGAS*XAMU/Mass_I(Ion3_)
      ! kB/m_e
      RGAS_I(nIon)=RGAS*XAMU/Mass_I(nIon)
   end select
      
  if (UseParticles) then
     
     call init_particle(nAltParticles,AltMinParticles,AltMaxParticles,&
          TypeParticleGrid)
     nLineParticle=nLine
     allocate(iLineGlobalParticle_I(nLine))
     iLineGlobalParticle_I=iLineGlobal

     !set altitude index for fluid to particle transition
     iAltParticle=ceiling((AltMinParticles-ALTD(1))/DRBND)

     !for each line restart or 
     !fill the state variables, intially sample the particles 
     !and then bury the line. Each lines particles will be disintered before 
     !advancing the calculation
     if(IsRestart .and. .not.DoInitAltParticles) then
        do iLine=1,nLine
           call read_restart_particle(iLine)
           call bury_line(iLine)
        enddo
     else
        
        allocate(Density_IC(nIon-1,nAlt),Velocity_IC(nIon-1,nAlt),Temperature_IC(nIon-1,nAlt))
        do iLine=1,nLine
           do iIon=1,nIon-1
              Density_IC(iIon,:)=State_CVI(1:nAlt,iRho_I(iIon),iLine)/Mass_I(iIon)
              Velocity_IC(iIon,:)=State_CVI(1:nAlt,iU_I(iIon),iLine)
              Temperature_IC(iIon,:)=State_CVI(1:nAlt,iT_I(iIon),iLine)
           enddo
           call put_to_particles(nAlt,nIon-1,ALTD(1:nAlt), &
                DoInitAltParticles,Density_IC,Velocity_IC,Temperature_IC)
           call bury_line(iLine)
           
        enddo
        deallocate(Density_IC,Velocity_IC,Temperature_IC)
     endif
  endif
end subroutine PW_initialize

!=============================================================================
integer function julianday(year, mon, day) result(Julian_Day)
  
  implicit none
  
  integer :: i
  integer, dimension(1:12) :: dayofmon
  integer :: year, mon, day
  
  dayofmon(1) = 31
  dayofmon(2) = 28
  dayofmon(3) = 31
  dayofmon(4) = 30
  dayofmon(5) = 31
  dayofmon(6) = 30
  dayofmon(7) = 31
  dayofmon(8) = 31
  dayofmon(9) = 30
  dayofmon(10) = 31
  dayofmon(11) = 30
  dayofmon(12) = 31
  
  if (mod(year,4).eq.0) dayofmon(2) = dayofmon(1) + 1
  Julian_Day = 0
  do i = 1, mon-1
     Julian_Day = Julian_Day + dayofmon(i)
  enddo
  Julian_Day = Julian_Day + day
  
end function julianday
