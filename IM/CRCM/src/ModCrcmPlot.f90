Module ModCrcmPlot
  implicit none
  private ! except
  public :: Crcm_plot
  character(len=5),  public    :: TypePlot   = 'ascii'
  logical,           public    :: DoSavePlot = .false.
  real,              public    :: DtOutput = 10.0
  character(len=17), parameter :: NameHeader = 'CRCM output_var11'  
contains
  !=============================================================================
  subroutine Crcm_plot(nLat,nLon, X_C,Y_C,Pressure_C,PressureHot_C,Rho_C,Beq_C,&
       Volume_C,Potential_C,Time,Dt)
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModCrcmInitialize, ONLY: IsRestart
    integer, intent(in) :: nLat, nLon
    real,    intent(in) :: X_C(nLat,nLon), Y_C(nLat,nLon), Time, Dt
    real,    intent(in) :: Pressure_C(nLat,nLon),Rho_C(nLat,nLon), & 
                           Beq_C(nLat,nLon),Volume_C(nLat,nLon),   &
                           Potential_C(nLat,nLon),PressureHot_C(nLat,nLon)
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)
    integer             :: nVar=6,iLat,iLon
    integer, parameter  :: x_=1, y_=2, nDim=2, P_=1,Phot_=2,Rho_=3,Beq_=4,Vol_=5, Pot_=6
    real                :: theta, phi
    character(len=14)   :: NamePlotEq = 'IM/CRCMeq.outs'
    character(len=79), parameter :: NamePlotVar='x y P[nP] Phot[nP] n[/m3] Beq[T] Vol[m3/Wb] Pot[Volts] g'
    real, parameter              :: gamma = 1.66666666666667
    logical, save                :: IsFirstCall = .true.
    !--------------------------------------------------------------------------
    
    allocate(Coord_DII(nDim,nLat,nLon+1), &
             PlotState_IIV(nLat,nLon+1,nVar))

    !Set Coords
    Coord_DII(x_,:, 1:nLon) = X_C(:,1:nLon)
    Coord_DII(y_,:, 1:nLon) = Y_C(:,1:nLon)
    
    !fill ghost cells of Coords
    Coord_DII(x_,:, nLon+1) = X_C(:,1)
    Coord_DII(y_,:, nLon+1) = Y_C(:,1)
    
    !Set plot data
    PlotState_IIV(:,1:nLon,P_)   = Pressure_C(:,1:nLon)    
    PlotState_IIV(:,1:nLon,Phot_)= PressureHot_C(:,1:nLon)    
    PlotState_IIV(:,1:nLon,Rho_) = Rho_C   (:,1:nLon)    
    PlotState_IIV(:,1:nLon,Beq_) = Beq_C   (:,1:nLon)    
    PlotState_IIV(:,1:nLon,Vol_) = Volume_C(:,1:nLon)    
    PlotState_IIV(:,1:nLon,Pot_) = Potential_C(:,1:nLon)    

    !fill ghost cells of plot data
    PlotState_IIV(:,nLon+1,P_) = Pressure_C(:,1)
    PlotState_IIV(:,nLon+1,Phot_) = PressureHot_C(:,1)
    PlotState_IIV(:,nLon+1,Rho_) = Rho_C   (:,1)    
    PlotState_IIV(:,nLon+1,Beq_) = Beq_C   (:,1)    
    PlotState_IIV(:,nLon+1,Vol_) = Volume_C(:,1)    
    PlotState_IIV(:,nLon+1,Pot_) = Potential_C(:,1)    


    !write plot
    if (IsRestart) then
       call save_plot_file(NamePlotEq, TypePositionIn='append',              &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,                 &
            NameVarIn = NamePlotVar, nStepIn= nint(Time/Dt),TimeIn=Time,     &
            nDimIn=2,CoordIn_DII=Coord_DII,                                  &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/gamma/))
    else
       if(IsFirstCall) then
          call save_plot_file(NamePlotEq, TypePositionIn='rewind',           &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,              &
               NameVarIn = NamePlotVar, nStepIn= nint(Time/Dt),TimeIn=Time,  &
               nDimIn=2,CoordIn_DII=Coord_DII,                               &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/gamma/))
          IsFirstCall = .false.
       else
          call save_plot_file(NamePlotEq, TypePositionIn='append',           &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,              &
               NameVarIn = NamePlotVar, nStepIn= nint(Time/Dt),TimeIn=Time,  &
               nDimIn=2,CoordIn_DII=Coord_DII,                               &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/gamma/))
          end if
    end if
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine Crcm_plot
end Module ModCrcmPlot
