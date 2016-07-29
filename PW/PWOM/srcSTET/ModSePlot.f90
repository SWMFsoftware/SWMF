Module ModSePlot
  implicit none
  
  private !except
  
  real, public :: DtSavePlot=300.0
  
  public :: plot_state
  public :: plot_state_pot
  public :: plot_omni_iono
  public :: plot_omni_iono_pot
  public :: plot_omni_line
  public :: plot_omni_pot
  public :: plot_along_field
contains
  !============================================================================
  ! save state plot for verification
  subroutine plot_state(iLine,nStep,time,iphiup,iphidn,phiup,phidn)
    use ModSeGrid,     ONLY: FieldLineGrid_IC,EqAngleGrid_IG,Bfield_IC,&
                             nIono,nPlas, nAngle, nEnergy, nLine, nPoint,&
                             nThetaAlt_II,rPlanetCM
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi

    integer, intent(in) :: iLine, nStep
    real,    intent(in) :: time
    real,    intent(in) :: iphiup(nLine,0:nAngle,0:2*nIono+1,nEnergy+1),&
                           iphidn(nLine,0:nAngle,0:2*nIono+1,nEnergy+1), &
                           phiup(nLine,0:nAngle,0:nPlas+1,nEnergy+1), &
                           phidn(nLine,0:nAngle,0:nPlas+1,nEnergy+1)
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)
    !grid parameters
    integer, parameter :: nDim =2, nVar=7, S_=1, PA_=2,B_=1
    integer, parameter :: E1_=3, E2_=4, E3_=5, E4_=6, E5_=7
    !set the corresponding energy channels for E1-E5
    integer, parameter :: EE1_=5, EE2_=10, EE3_=20, EE4_=40, EE5_=80
    character(len=100),parameter :: NamePlotVar='S PA B PA E1 E2 E3 E4 E5 g r'
    character(len=13) :: NamePlot='StatePlot.out'
    character(len=*),parameter :: NameHeader='SE output'
    character(len=5) :: TypePlot='ascii'
    integer :: iAngle,iAngleDn,iPoint,iIono,iPlas
    !--------------------------------------------------------------------------
    allocate(Coord_DII(nDim,nPoint,2*nAngle),PlotState_IIV(nPoint,2*nAngle,nVar))

!    do iLine=1,nLine
       PlotState_IIV = 0.0
       Coord_DII     = 0.0
       
       !Set Coordinates along field line and PA
       do iPoint=1,nPoint
          do iAngle=1,2*nAngle
             !set coord based on up or down region
             if(iAngle<nAngle) then
                Coord_DII(S_,iPoint,iAngle) = FieldLineGrid_IC(iLine,iPoint)&
                     /rPlanetCM
                Coord_DII(PA_,iPoint,iAngle)= EqAngleGrid_IG(iLine,iAngle)
             else
                iAngleDn = 2*nAngle-iAngle
                Coord_DII(S_,iPoint,iAngle) = FieldLineGrid_IC(iLine,iPoint)&
                     /rPlanetCM
                Coord_DII(PA_,iPoint,iAngle)= &
                     cPi-EqAngleGrid_IG(iLine,iAngleDn)
             endif
             
             !set plotstate based on up or down region  
             if (iAngle <= nThetaAlt_II(iLine,iPoint))then
                PlotState_IIV(iPoint,iAngle,B_)  = Bfield_IC(iLine,iPoint)
                PlotState_IIV(iPoint,iAngle,PA_) = EqAngleGrid_IG(iLine,iAngle)
                !choose phiup or iphiup by spatial region
                if (iPoint <= nIono)then
                   iIono=iPoint
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        iphiup(iLine,iAngle,iIono,EE1_)
                   PlotState_IIV(iPoint,iAngle,E2_) = &
                        iphiup(iLine,iAngle,iIono,EE2_)
                   PlotState_IIV(iPoint,iAngle,E3_) = &
                        iphiup(iLine,iAngle,iIono,EE3_)
                   PlotState_IIV(iPoint,iAngle,E4_) = &
                        iphiup(iLine,iAngle,iIono,EE4_)
                   PlotState_IIV(iPoint,iAngle,E5_) = &
                        iphiup(iLine,iAngle,iIono,EE5_)
                elseif(iPoint >nIono .and. iPoint <=nPoint-nIono) then
                   iPlas=iPoint-nIono
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        phiup(iLine,iAngle,iPlas,EE1_)
                   PlotState_IIV(iPoint,iAngle,E2_) = &
                        phiup(iLine,iAngle,iPlas,EE2_)
                   PlotState_IIV(iPoint,iAngle,E3_) = &
                        phiup(iLine,iAngle,iPlas,EE3_)
                   PlotState_IIV(iPoint,iAngle,E4_) = &
                        phiup(iLine,iAngle,iPlas,EE4_)
                   PlotState_IIV(iPoint,iAngle,E5_) = &
                        phiup(iLine,iAngle,iPlas,EE5_)
                else
                   iIono=iPoint-nPlas
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        iphiup(iLine,iAngle,iIono,EE1_)
                   PlotState_IIV(iPoint,iAngle,E2_) = &
                        iphiup(iLine,iAngle,iIono,EE2_)
                   PlotState_IIV(iPoint,iAngle,E3_) = &
                        iphiup(iLine,iAngle,iIono,EE3_)
                   PlotState_IIV(iPoint,iAngle,E4_) = &
                        iphiup(iLine,iAngle,iIono,EE4_)
                   PlotState_IIV(iPoint,iAngle,E5_) = &
                        iphiup(iLine,iAngle,iIono,EE5_)
                endif
             elseif(iAngle > 2*nAngle-nThetaAlt_II(iLine,iPoint)) then
                iAngleDn = 2*nAngle-iAngle
                PlotState_IIV(iPoint,iAngle,B_)  = Bfield_IC(iLine,iPoint)
                PlotState_IIV(iPoint,iAngle,PA_) = &
                     cPi-EqAngleGrid_IG(iLine,iAngleDn)
                !choose phidn or iphidn by spatial region
                if (iPoint <= nIono)then
                   iIono=iPoint
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        iphidn(iLine,iAngleDn,iIono,EE1_)
                   PlotState_IIV(iPoint,iAngle,E2_) = &
                        iphidn(iLine,iAngleDn,iIono,EE2_)
                   PlotState_IIV(iPoint,iAngle,E3_) = &
                        iphidn(iLine,iAngleDn,iIono,EE3_)
                   PlotState_IIV(iPoint,iAngle,E4_) = &
                        iphidn(iLine,iAngleDn,iIono,EE4_)
                   PlotState_IIV(iPoint,iAngle,E5_) = &
                        iphidn(iLine,iAngleDn,iIono,EE5_)
                elseif(iPoint >nIono .and. iPoint <=nPoint-nIono) then
                   iPlas=iPoint-nIono
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        phidn(iLine,iAngleDn,iPlas,EE1_)
                   PlotState_IIV(iPoint,iAngle,E2_) = &
                        phidn(iLine,iAngleDn,iPlas,EE2_)
                   PlotState_IIV(iPoint,iAngle,E3_) = &
                        phidn(iLine,iAngleDn,iPlas,EE3_)
                   PlotState_IIV(iPoint,iAngle,E4_) = &
                        phidn(iLine,iAngleDn,iPlas,EE4_)
                   PlotState_IIV(iPoint,iAngle,E5_) = &
                        phidn(iLine,iAngleDn,iPlas,EE5_)
                else
                   iIono=iPoint-nPlas
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        iphidn(iLine,iAngleDn,iIono,EE1_)
                   PlotState_IIV(iPoint,iAngle,E2_) = &
                        iphidn(iLine,iAngleDn,iIono,EE2_)
                   PlotState_IIV(iPoint,iAngle,E3_) = &
                        iphidn(iLine,iAngleDn,iIono,EE3_)
                   PlotState_IIV(iPoint,iAngle,E4_) = &
                        iphidn(iLine,iAngleDn,iIono,EE4_)
                   PlotState_IIV(iPoint,iAngle,E5_) = &
                        iphidn(iLine,iAngleDn,iIono,EE5_)
                endif
             else
                PlotState_IIV(iPoint,iAngle,:)=0.0
             endif
          enddo
       enddo
       
       !Plot grid for given line
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
 !   end do
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_state
  !============================================================================
  ! save state plot for verification in case with potential. Only plot 
  ! at a particular energy, iEnergyIn, since now grid is energy dependent
  subroutine plot_state_pot(iLine,iEnergyIn,nStep,time,iphiup,iphidn,phiup,phidn)
    use ModSeGrid,     ONLY: FieldLineGrid_IC,EqAngleGrid_IG,Bfield_IC,&
                             nIono,nPlas, nAngle, nEnergy, nLine, nPoint,&
                             nThetaAlt_IIC,EnergyGrid_I,MaxAlt_IC,rPlanetCM
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi

    integer, intent(in) :: iLine, iEnergyIn,nStep
    real,    intent(in) :: time
    real,    intent(in) :: iphiup(nLine,0:nAngle,0:2*nIono+1,nEnergy+1),&
                           iphidn(nLine,0:nAngle,0:2*nIono+1,nEnergy+1), &
                           phiup(nLine,0:nAngle,0:nPlas+1,nEnergy+1), &
                           phidn(nLine,0:nAngle,0:nPlas+1,nEnergy+1)
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)

    !grid parameters
    integer, parameter :: nDim =2, nVar=3, S_=1, PA_=2,B_=1
    integer, parameter :: E1_=3
    character(len=100),parameter :: NamePlotVar='S PA B PA Phi g r'
    character(len=100) :: NamePlot
    character(len=*),parameter :: NameHeader='SE output'
    character(len=5) :: TypePlot='ascii'
    integer :: iAngle,iAngleDn,iPoint,iIono,iPlas
        logical :: IsFirstCall=.true.
    !--------------------------------------------------------------------------
    allocate(Coord_DII(nDim,nPoint,2*nAngle),PlotState_IIV(nPoint,2*nAngle,nVar))

!    do iLine=1,nLine
       PlotState_IIV = 0.0
       Coord_DII     = 0.0
       
       !Set Coordinates along field line and PA
       ALONG_LINE: do iPoint=1,nPoint
          ANGLE: do iAngle=1,2*nAngle
             !set coord based on up or down region
             if(iAngle<nAngle) then
                Coord_DII(S_,iPoint,iAngle) = FieldLineGrid_IC(iLine,iPoint)&
                     /rPlanetCM
                Coord_DII(PA_,iPoint,iAngle)= EqAngleGrid_IG(iLine,iAngle)
             else
                iAngleDn = 2*nAngle-iAngle
                Coord_DII(S_,iPoint,iAngle) = FieldLineGrid_IC(iLine,iPoint)&
                     /rPlanetCM
                Coord_DII(PA_,iPoint,iAngle)= &
                     cPi-EqAngleGrid_IG(iLine,iAngleDn)
             endif
             
             !if above maxalt for energy then plot state is zero and cycle
             if(iPoint>MaxAlt_IC(iLine,iEnergyIn) &
                  .and. iPoint<nPoint-MaxAlt_IC(iLine,iEnergyIn)) then
                PlotState_IIV(iPoint,iAngle,:)=0.0
                cycle ANGLE
             endif

             !set plotstate based on up or down region  
             if (iAngle <= nThetaAlt_IIC(iLine,iEnergyIn,iPoint))then
                PlotState_IIV(iPoint,iAngle,B_)  = Bfield_IC(iLine,iPoint)
                PlotState_IIV(iPoint,iAngle,PA_) = EqAngleGrid_IG(iLine,iAngle)
                !choose phiup or iphiup by spatial region
                if (iPoint <= nIono)then
                   iIono=iPoint
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        iphiup(iLine,iAngle,iIono,iEnergyIn)
                 elseif(iPoint >nIono .and. iPoint <=nPoint-nIono) then
                   iPlas=iPoint-nIono
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        phiup(iLine,iAngle,iPlas,iEnergyIn)
                else
                   iIono=iPoint-nPlas
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        iphiup(iLine,iAngle,iIono,iEnergyIn)
                endif
             elseif(iAngle > 2*nAngle-nThetaAlt_IIC(iLine,iEnergyIn,iPoint)) then
                iAngleDn = 2*nAngle-iAngle
                PlotState_IIV(iPoint,iAngle,B_)  = Bfield_IC(iLine,iPoint)
                PlotState_IIV(iPoint,iAngle,PA_) = &
                     cPi-EqAngleGrid_IG(iLine,iAngleDn)
                !choose phidn or iphidn by spatial region
                if (iPoint <= nIono)then
                   iIono=iPoint
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        iphidn(iLine,iAngleDn,iIono,iEnergyIn)
                elseif(iPoint >nIono .and. iPoint <=nPoint-nIono) then
                   iPlas=iPoint-nIono
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        phidn(iLine,iAngleDn,iPlas,iEnergyIn)
                else
                   iIono=iPoint-nPlas
                   PlotState_IIV(iPoint,iAngle,E1_) = &
                        iphidn(iLine,iAngleDn,iIono,iEnergyIn)
                endif
             else
                PlotState_IIV(iPoint,iAngle,:)=0.0
             endif
          enddo ANGLE
       enddo ALONG_LINE
  
       ! set name for plotfile
       write(NamePlot,"(a,i4.4,a)") 'State_',floor(EnergyGrid_I(iEnergyIn)+.001),'.out'
       
       if(IsFirstCall) then
          !Plot grid for given line
          call save_plot_file(NamePlot, TypePositionIn='rewind', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
               nDimIn=nDim,CoordIn_DII=Coord_DII,                &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
          IsFirstCall = .false.
       else
          call save_plot_file(NamePlot, TypePositionIn='append', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
               nDimIn=nDim,CoordIn_DII=Coord_DII,                &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
       endif
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_state_pot
  
  !============================================================================
  ! plot omnidirectional flux in the ionosphere
  subroutine plot_omni_iono(iLine,nStep,time,specup,specdn,IsIono1)
    use ModSeGrid,     ONLY: FieldLineGrid_IC, nIono, nEnergy, nLine, &
         nPoint,EnergyGrid_I,iLineGlobal_I,rPlanetCM
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi
    
    integer, intent(in) :: iLine, nStep
    real,    intent(in) :: time
    real,    intent(in) :: specup(nLine,nEnergy,nPoint),&
         specdn(nLine,nEnergy,nPoint)
    logical, intent(in) :: IsIono1
    
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)

    !grid parameters
    integer, parameter :: nDim =2, nVar=3, S_=2, E_=1
    !variable parameters
    integer, parameter :: Flux_=1, FluxUp_=2, FluxDn_=3
    
    character(len=100),parameter :: NamePlotVar='E[eV] Alt[km] Omni[cm-3eV-1s-1] OmniUp[cm-3eV-1s-1] OmniDn[cm-3eV-1s-1] g r'
    character(len=100) :: NamePlot1='OmniIono1.out',NamePlot2='OmniIono.out'
    character(len=100) :: NamePlot
    
    character(len=*),parameter :: NameHeader='SE output iono'
    character(len=5) :: TypePlot='ascii'
    integer :: iAngle,iAngleDn,iIono,iEnergy, iPoint

    logical :: IsFirstCall1=.true.,IsFirstCall2=.true.
    !--------------------------------------------------------------------------
    
    allocate(Coord_DII(nDim,nEnergy,nIono),PlotState_IIV(nEnergy,nIono,nVar))
    
    !    do iLine=1,nLine
    PlotState_IIV = 0.0
    Coord_DII     = 0.0
    
    !Set values
    if (IsIono1) then
       do iEnergy=1,nEnergy
          do iIono=1,nIono
             !relate iPoint and iIono
             iPoint = iIono
             
             Coord_DII(E_,iEnergy,iIono) = EnergyGrid_I(iEnergy)             
             Coord_DII(S_,iEnergy,iIono) = FieldLineGrid_IC(iLine,iPoint)/1e5
             ! set plot state
             PlotState_IIV(iEnergy,iIono,Flux_)  = &
                  specup(iLine,iEnergy,iPoint)-specdn(iLine,iEnergy,iPoint)
             PlotState_IIV(iEnergy,iIono,FluxUp_)  = &
                  specup(iLine,iEnergy,iPoint)
             PlotState_IIV(iEnergy,iIono,FluxDn_)  = &
                  specdn(iLine,iEnergy,iPoint)
          enddo
       enddo
       
       ! set name for plotfile
       write(NamePlot,"(a,i4.4,a)") 'OmniIono1_',iLineGlobal_I(iLine),'.out'
  
       !Plot grid for given line
       if(IsFirstCall1) then
          call save_plot_file(NamePlot, TypePositionIn='rewind', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
               nDimIn=nDim,CoordIn_DII=Coord_DII,                &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
          IsFirstCall1 = .false.
       else
          call save_plot_file(NamePlot, TypePositionIn='append', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
               nDimIn=nDim,CoordIn_DII=Coord_DII,                &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
       endif
  
    else
       do iEnergy=1,nEnergy
          do iIono=1,nIono
             !relate iPoint and iIono
             iPoint = nPoint-iIono + 1
             Coord_DII(E_,iEnergy,iIono) = EnergyGrid_I(iEnergy)             
             Coord_DII(S_,iEnergy,iIono) = &
                  (FieldLineGrid_IC(iLine,nPoint) &
                  - FieldLineGrid_IC(iLine,iPoint))/1e5
             !set plot state
             PlotState_IIV(iEnergy,iIono,Flux_)  = &
                  specup(iLine,iEnergy,iPoint)-specdn(iLine,iEnergy,iPoint)
             PlotState_IIV(iEnergy,iIono,FluxUp_)  = &
                  specup(iLine,iEnergy,iPoint)
             PlotState_IIV(iEnergy,iIono,FluxDn_)  = &
                  specdn(iLine,iEnergy,iPoint)
          enddo
       enddo
       ! set name for plotfile
       write(NamePlot,"(a,i4.4,a)") 'OmniIono2_',iLineGlobal_I(iLine),'.out'
       
       !Plot grid for given line
       if(IsFirstCall2) then
          call save_plot_file(NamePlot, TypePositionIn='rewind', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
               nDimIn=nDim,CoordIn_DII=Coord_DII,                &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
          IsFirstCall2 = .false.
       else
          call save_plot_file(NamePlot, TypePositionIn='append', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
               nDimIn=nDim,CoordIn_DII=Coord_DII,                &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
       endif
  
    endif
    
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_omni_iono
     
  !============================================================================
  ! plot omnidirectional flux in the ionosphere
  subroutine plot_omni_iono_pot(iLine,nStep,time,specup,specdn,IsIono1)
    use ModSeGrid,     ONLY: FieldLineGrid_IC, nIono, nEnergy, nLine, &
         nPoint,KineticEnergy_IIC,iLineGlobal_I
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi
    
    integer, intent(in) :: iLine, nStep
    real,    intent(in) :: time
    real,    intent(in) :: specup(nLine,nEnergy,nPoint),&
         specdn(nLine,nEnergy,nPoint)
    logical, intent(in) :: IsIono1
    
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)

    !grid parameters
    integer, parameter :: nDim =2, nVar=3, S_=2, E_=1
    !variable parameters
    integer, parameter :: Flux_=1, FluxUp_=2, FluxDn_=3
    
    character(len=100),parameter :: NamePlotVar='E[eV] Alt[km] Omni[cm-3eV-1s-1] OmniUp[cm-3eV-1s-1] OmniDn[cm-3eV-1s-1] g r'
    character(len=100) :: NamePlot1='OmniIono1.out',NamePlot2='OmniIono.out'
    character(len=100) :: NamePlot
    
    character(len=*),parameter :: NameHeader='SE output iono'
    character(len=5) :: TypePlot='ascii'
    integer :: iAngle,iAngleDn,iIono,iEnergy, iPoint

    logical :: IsFirstCall1=.true.,IsFirstCall2=.true.
    !--------------------------------------------------------------------------
    
    allocate(Coord_DII(nDim,nEnergy,nIono),PlotState_IIV(nEnergy,nIono,nVar))
    
    !    do iLine=1,nLine
    PlotState_IIV = 0.0
    Coord_DII     = 0.0
    
    !Set values
    if (IsIono1) then
       do iEnergy=1,nEnergy
          do iIono=1,nIono
             !relate iPoint and iIono
             iPoint = iIono
             
             Coord_DII(E_,iEnergy,iIono) = KineticEnergy_IIC(iLine,iEnergy,iPoint)             
             Coord_DII(S_,iEnergy,iIono) = FieldLineGrid_IC(iLine,iPoint)/1e5
             ! set plot state
             PlotState_IIV(iEnergy,iIono,Flux_)  = &
                  specup(iLine,iEnergy,iPoint)-specdn(iLine,iEnergy,iPoint)
             PlotState_IIV(iEnergy,iIono,FluxUp_)  = &
                  specup(iLine,iEnergy,iPoint)
             PlotState_IIV(iEnergy,iIono,FluxDn_)  = &
                  specdn(iLine,iEnergy,iPoint)
          enddo
       enddo
       
       ! set name for plotfile
       write(NamePlot,"(a,i4.4,a)") 'OmniIono1_',iLineGlobal_I(iLine),'.out'
  
       !Plot grid for given line
       if(IsFirstCall1) then
          call save_plot_file(NamePlot, TypePositionIn='rewind', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
               nDimIn=nDim,CoordIn_DII=Coord_DII,                &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
          IsFirstCall1 = .false.
       else
          call save_plot_file(NamePlot, TypePositionIn='append', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
               nDimIn=nDim,CoordIn_DII=Coord_DII,                &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
       endif
  
    else
       do iEnergy=1,nEnergy
          do iIono=1,nIono
             !relate iPoint and iIono
             iPoint = nPoint-iIono + 1
             Coord_DII(E_,iEnergy,iIono) = KineticEnergy_IIC(iLine,iEnergy,iPoint)             
             Coord_DII(S_,iEnergy,iIono) = &
                  (FieldLineGrid_IC(iLine,nPoint) &
                  - FieldLineGrid_IC(iLine,iPoint))/1e5
             !set plot state
             PlotState_IIV(iEnergy,iIono,Flux_)  = &
                  specup(iLine,iEnergy,iPoint)-specdn(iLine,iEnergy,iPoint)
             PlotState_IIV(iEnergy,iIono,FluxUp_)  = &
                  specup(iLine,iEnergy,iPoint)
             PlotState_IIV(iEnergy,iIono,FluxDn_)  = &
                  specdn(iLine,iEnergy,iPoint)
          enddo
       enddo
       ! set name for plotfile
       write(NamePlot,"(a,i4.4,a)") 'OmniIono2_',iLine,'.out'
       
       !Plot grid for given line
       if(IsFirstCall2) then
          call save_plot_file(NamePlot, TypePositionIn='rewind', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
               nDimIn=nDim,CoordIn_DII=Coord_DII,                &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
          IsFirstCall2 = .false.
       else
          call save_plot_file(NamePlot, TypePositionIn='append', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
               nDimIn=nDim,CoordIn_DII=Coord_DII,                &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
       endif
  
    endif
    
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_omni_iono_pot
     
  !============================================================================
  ! plot omnidirectional flux along the entire line
  subroutine plot_omni_line(iLine,nStep,time,specup,specdn)
    use ModSeGrid,     ONLY: FieldLineGrid_IC, nEnergy, nLine, &
         nPoint,EnergyGrid_I,iLineGlobal_I
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi
    
    integer, intent(in) :: iLine, nStep
    real,    intent(in) :: time
    real,    intent(in) :: specup(nLine,nEnergy,nPoint),&
         specdn(nLine,nEnergy,nPoint)
    
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)

    !grid parameters
    integer, parameter :: nDim =2, nVar=3, S_=2, E_=1
    !variable parameters
    integer, parameter :: Flux_=1, FluxUp_=2, FluxDn_=3
    
    character(len=100),parameter :: NamePlotVar='E[eV] Alt[km] Omni[cm-3eV-1s-1] OmniUp[cm-3eV-1s-1] OmniDn[cm-3eV-1s-1] g r'
    character(len=100) :: NamePlot1='OmniIono1.out',NamePlot2='OmniIono.out'
    character(len=100) :: NamePlot
    
    character(len=*),parameter :: NameHeader='SE output Fluxes'
    character(len=5) :: TypePlot='ascii'
    integer :: iAngle,iAngleDn,iIono,iEnergy, iPoint

    logical :: IsFirstCall=.true.
    !--------------------------------------------------------------------------
    
    allocate(Coord_DII(nDim,nEnergy,nPoint),PlotState_IIV(nEnergy,nPoint,nVar))
    
    !    do iLine=1,nLine
    PlotState_IIV = 0.0
    Coord_DII     = 0.0
    
    !Set values
    do iEnergy=1,nEnergy
       do iPoint=1,nPoint
          Coord_DII(E_,iEnergy,iPoint) = EnergyGrid_I(iEnergy)             
          Coord_DII(S_,iEnergy,iPoint) = FieldLineGrid_IC(iLine,iPoint)/1e5
          ! set plot state
          PlotState_IIV(iEnergy,iPoint,Flux_)  = &
               specup(iLine,iEnergy,iPoint)-specdn(iLine,iEnergy,iPoint)
          PlotState_IIV(iEnergy,iPoint,FluxUp_)  = &
               specup(iLine,iEnergy,iPoint)
          PlotState_IIV(iEnergy,iPoint,FluxDn_)  = &
               specdn(iLine,iEnergy,iPoint)
       enddo
    enddo
    
    ! set name for plotfile
    write(NamePlot,"(a,i4.4,a)") 'Fluxes_',iLineGlobal_I(iLine),'.out'
    
    !Plot grid for given line
    if(IsFirstCall) then
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
       IsFirstCall = .false.
    else
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
    endif
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_omni_line
  
  
  !===========================================================================  
  !============================================================================
  ! plot omnidirectional flux along the entire line in case of potential. Note 
  ! that the energy axis is now non-uniform since we are looking at kinetic 
  ! energy output.
  subroutine plot_omni_pot(iLine,nStep,time,specup,specdn)
    use ModSeGrid,     ONLY: FieldLineGrid_IC, nEnergy, nLine, &
         nPoint,KineticEnergy_IIC,iLineGlobal_I, UsePwRegion, nPwRegion, nIono
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi
    
    integer, intent(in) :: iLine, nStep
    real,    intent(in) :: time
    real,    intent(in) :: specup(nLine,nEnergy,nPoint),&
         specdn(nLine,nEnergy,nPoint)
    
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)
    !grid parameters
    integer, parameter :: nDim =2, nVar=3, S_=2, E_=1
    !variable parameters
    integer, parameter :: Flux_=1, FluxUp_=2, FluxDn_=3
    
    character(len=100),parameter :: NamePlotVar='E[eV] Alt[km] Omni[cm-3eV-1s-1] OmniUp[cm-3eV-1s-1] OmniDn[cm-3eV-1s-1] g r'
    character(len=100) :: NamePlot1='OmniIono1.out',NamePlot2='OmniIono.out'
    character(len=100) :: NamePlot
    
    character(len=*),parameter :: NameHeader='SE output Fluxes'
    character(len=5) :: TypePlot='ascii'
    integer :: iAngle,iAngleDn,iIono,iEnergy, iPoint, nPointUsed
    logical :: IsFirstCall=.true.
    !--------------------------------------------------------------------------

    nPointUsed = nPoint
    if(UsePwRegion) nPointUsed=nIono+nPwRegion
    
    allocate(Coord_DII(nDim,nEnergy,nPointUsed),PlotState_IIV(nEnergy,nPointUsed,nVar))
    
    !    do iLine=1,nLine
    PlotState_IIV = 0.0
    Coord_DII     = 0.0
    
    !Set values
    do iEnergy=1,nEnergy
       do iPoint=1,nPointUsed
          Coord_DII(E_,iEnergy,iPoint) = KineticEnergy_IIC(iLine,iEnergy,iPoint)
          Coord_DII(S_,iEnergy,iPoint) = FieldLineGrid_IC(iLine,iPoint)/1e5
          ! set plot state
          PlotState_IIV(iEnergy,iPoint,Flux_)  = &
               specup(iLine,iEnergy,iPoint)-specdn(iLine,iEnergy,iPoint)
          PlotState_IIV(iEnergy,iPoint,FluxUp_)  = &
               specup(iLine,iEnergy,iPoint)
          PlotState_IIV(iEnergy,iPoint,FluxDn_)  = &
               specdn(iLine,iEnergy,iPoint)
       enddo
    enddo
    
    ! set name for plotfile
    write(NamePlot,"(a,i4.4,a)") 'Fluxes_',iLineGlobal_I(iLine),'.out'
    
    !Plot grid for given line
    if(IsFirstCall) then
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
       IsFirstCall = .false.
    else
       write(*,*)'!!! nPoint, nPointUsed=', nPoint, nPointUsed
       write(*,*)'!!! max(specup)=', maxval(specup(:,:,1:nPointUsed))
       write(*,*)'!!! max(specdn)=', maxval(specdn(:,:,1:nPointUsed))
       write(*,*)'!!! max(KineticEnergy_IIC)=',maxval(KineticEnergy_IIC(:,:,1:nPointUsed))
       write(*,*)'!!! max( FieldLineGrid_IC)=',maxval(FieldLineGrid_IC(:,1:nPointUsed))
       write(*,*)'!!! max(Var)=', maxval(PlotState_IIV)
       write(*,*)'!!! max(coord)=', maxval(Coord_DII)


       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
    endif
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_omni_pot
  
  
  !============================================================================ 


  !============================================================================
  ! 1D output plots of integrated quantities along the field 
  ! (potential, heating rate, SE number density, SE number flux)
  subroutine plot_along_field(iLine,time,HeatingRate_IC,NumberDens_IC,&
       NumberFlux_IC,TotalIonizationRate_IC)
    use ModSeGrid,     ONLY: FieldLineGrid_IC, DeltaPot_IC, nLine, nPoint, &
                             MinEnergy_IC, EnergyGrid_I,Efield_IC,iLineGlobal_I,&
                             rPlanetCM

    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg
    integer, intent(in) :: iLine
    real, intent(in) :: time,HeatingRate_IC(nLine,nPoint), &
         NumberDens_IC(nLine,nPoint),NumberFlux_IC(nLine,nPoint),&
         TotalIonizationRate_IC(nLine,nPoint)

    real, allocatable   :: Coord_I(:), PlotState_IV(:,:)
    integer, parameter :: nDim =1, nVar=6, Pot_=1, Qe_=2,Nse_=3,Fse_=4,E_=5,&
                          rate_=6
    character(len=100),parameter :: NamePlotVar=&
         'S Pot[eV] Qe[eV/cm3/s] Nse[/cc] Fse[/cm2/s] E[V/m] IonRate[/cm3/s] g r'
    character(len=100) :: NamePlot
    character(len=*),parameter :: NameHeader='Integrated output'
    character(len=5) :: TypePlot='ascii'
    integer :: iPoint
    logical,save :: IsFirstCall =.true.
    !--------------------------------------------------------------------------
    allocate(Coord_I(nPoint), PlotState_IV(nPoint,nVar))
    
    PlotState_IV = 0.0
    Coord_I     = 0.0
    
    !Set Coordinates along field line and PA
    do iPoint=1,nPoint
       Coord_I(iPoint) = FieldLineGrid_IC(iLine,iPoint)/rPlanetCM
       PlotState_IV(iPoint,Pot_) = DeltaPot_IC(iLine,iPoint)
       PlotState_IV(iPoint,Qe_)  = HeatingRate_IC(iLine,iPoint)
       PlotState_IV(iPoint,Nse_) = NumberDens_IC(iLine,iPoint)
       PlotState_IV(iPoint,Fse_) = NumberFlux_IC(iLine,iPoint)
       PlotState_IV(iPoint,E_)   = Efield_IC(iLine,iPoint)
       PlotState_IV(iPoint,rate_)= TotalIonizationRate_IC(iLine,iPoint)
    enddo
    
    ! set name for plotfile
    write(NamePlot,"(a,i4.4,a)") 'STET_1D_iLine',iLineGlobal_I(iLine),'.out'
    
    !Plot grid for given line. Overwrite old results on firstcall
    if(IsFirstCall) then
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn= 1,TimeIn=Time,     &
            nDimIn=nDim,CoordIn_I=Coord_I,                &
            VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/))
       if (iLine==nLine) IsFirstCall = .false.
    else
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn= 1,TimeIn=Time,     &
            nDimIn=nDim,CoordIn_I=Coord_I,                &
            VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/))
    end if
    
    
    deallocate(Coord_I, PlotState_IV)
  end subroutine plot_along_field
  

  !============================================================================

end Module ModSePlot
