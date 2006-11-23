module ModFieldLine

  use ModParameters
  implicit none

  private
  public put_field_line, get_field_line

  real, Dimension(maxGrid)  ::  dOxygPW,uOxygPW,pOxygPW,TOxygPW, &
       dHelPW,uHelPW,pHelPW,THelPW,     &
       dHydPW,uHydPW,pHydPW,THydPW,     &
       dElectPW,uElectPW,pElectPW,TElectPW

  logical :: IsRestartPW, IsVariableDtPW 
  real    :: TimePW,MaxLineTimePW,DToutputPW, DTpolarwindPW,GeoMagLatPW,&
       GeoMagLonPW,JrPW,nDimPW
  real    :: wHorizontalPW
  integer :: iUnitInputPW,      &
       iUnitOutputPW,iUnitGraphicsPW,                 &
       iUnitSourceGraphicsPW,iUnitRestartPW,          &
       iUnitCollisionPW,iUnitRestartInPW,iLinePW,nLinePW
  CHARACTER(7) :: TypeSolverPW
  character*100 :: NameRestartPW

contains

  !***************************************************************************
  !  Put polarwind variables into Mod_PW for passing 
  !***************************************************************************

  subroutine put_field_line(dOxyg, uOxyg, pOxyg, TOxyg,     &
       dHel, uHel, pHel, THel,                              &
       dHyd, uHyd, pHyd, THyd,                              &
       dElect, uElect, pElect, TElect,                      &
       GeoMagLat,GeoMagLon,Jr,wHorizontal,                  &
       iUnitOutput,iUnitGraphics, NameRestart,iLine,Time,   &
       MaxLineTime,TypeSolver,IsVariableDt,IsRestart,DToutput,nDim )

    use ModParameters

    real, intent(in),dimension(maxGrid):: dOxyg, uOxyg, pOxyg, TOxyg,     &
         dHel, uHel, pHel, THel,         &
         dHyd, uHyd, pHyd, THyd,         &
         dElect, uElect, pElect, TElect
    real, optional, intent(in) :: Time,MaxLineTime,DToutput

    real,    intent(in)     :: GeoMagLat,GeoMagLon,Jr,wHorizontal              
    integer, optional,intent(in)     :: iUnitOutput,iUnitGraphics,iLine,nDim
    character*100,optional,intent(in):: NameRestart
    character(7),optional,intent(in)::TypeSolver
    logical,optional,intent(in) :: IsVariableDt,IsRestart
    !-------------------------------------------------------------------------

    dOxygPW (:) = dOxyg(:)
    uOxygPW (:) = uOxyg(:)
    pOxygPW (:) = pOxyg(:)
    TOxygPW (:) = TOxyg(:)
    dHelPW  (:) = dHel(:)
    uHelPW  (:) = uHel(:)
    pHelPW  (:) = pHel(:)
    THelPW  (:) = THel(:)
    dHydPW  (:) = dHyd(:)
    uHydPW  (:) = uHyd(:)
    pHydPW  (:) = pHyd(:)
    THydPW  (:) = THyd(:)
    dElectPW(:) = dElect(:)
    uElectPW(:) = uElect(:)
    pElectPW(:) = pElect(:)
    TElectPW(:) = TElect(:) 
    GeoMagLatPW = GeoMagLat
    GeoMagLonPW = GeoMagLon
    JrPW        = Jr
    wHorizontalPW   = wHorizontal
    

    if (present(nDim))          nDimPW = nDim
    
    if (present(Time))          TimePW = Time
    if (present(MaxLineTime))   MaxLineTimePW = MaxLineTime
    if (present(iUnitGraphics)) iUnitGraphicsPW =iUnitGraphics
    if (present(NameRestart))   NameRestartPW=NameRestart
    if (present(iLine))         iLinePW = iLine
    if (present(iUnitOutput))   iUnitOutputPW=iUnitOutput
    if (present(TypeSolver))    TypeSolverPW =TypeSolver
    if (present(IsVariableDt))  IsVariableDtPW=IsVariableDt
    if (present(IsRestart))     IsRestartPW=IsRestart
    if (present(DToutput))      DToutputPW=DToutput
  end subroutine put_field_line

  !***************************************************************************
  !  Get polarwind variables from Mod_PW 
  !***************************************************************************

  subroutine get_field_line(dOxyg, uOxyg, pOxyg, TOxyg,     &
       dHel, uHel, pHel, THel,                              &
       dHyd, uHyd, pHyd, THyd,                              &
       dElect, uElect, pElect, TElect,                      &
       GeoMagLat,GeoMagLon,Jr,wHorizontal,                  &
       iUnitOutput,iUnitGraphics, NameRestart,iLine,Time,   &
       MaxLineTime,TypeSolver,IsVariableDt,IsRestart,DToutput, nDim)

    use ModParameters

    real, intent(out),dimension(maxGrid):: dOxyg, uOxyg, pOxyg, TOxyg,     &
         dHel, uHel, pHel, THel,         &
         dHyd, uHyd, pHyd, THyd,         &
         dElect, uElect, pElect, TElect


    real,    intent(out)     :: GeoMagLat, GeoMagLon,Jr, wHorizontal           
    
    character*100,optional,intent(out):: NameRestart
    character(7),optional,intent(out):: TypeSolver
    logical,optional,intent(out)      :: IsVariableDt,IsRestart
    real, optional, intent(out)       :: Time,MaxLineTime,DToutput
    integer, optional,intent(out)     :: iUnitOutput,iUnitGraphics,iLine,nDim
    
    dOxyg (:) = dOxygPW(:)
    uOxyg (:) = uOxygPW(:)
    pOxyg (:) = pOxygPW(:)
    TOxyg (:) = TOxygPW(:)
    dHel  (:) = dHelPW(:)
    uHel  (:) = uHelPW(:)
    pHel  (:) = pHelPW(:)
    THel  (:) = THelPW(:)
    dHyd  (:) = dHydPW(:)
    uHyd  (:) = uHydPW(:)
    pHyd  (:) = pHydPW(:)
    THyd  (:) = THydPW(:)
    dElect(:) = dElectPW(:)
    uElect(:) = uElectPW(:)
    pElect(:) = pElectPW(:)
    TElect(:) = TElectPW(:) 
    GeoMagLat = GeoMagLatPW
    GeoMagLon = GeoMagLonPW
    Jr        = JrPW
    wHorizontal = wHorizontalPW
   

    

    if (present(nDim))          nDim = nDimPW
    if (present(Time))          Time = TimePW
    if (present(MaxLineTime))   MaxLineTime = MaxLineTimePW
    if (present(iUnitGraphics)) iUnitGraphics =iUnitGraphicsPW
    if (present(NameRestart))   NameRestart=NameRestartPW
    if (present(iLine))         iLine = iLinePW
    if (present(iUnitOutput))   iUnitOutput=iUnitOutputPW
    if (present(TypeSolver))    TypeSolver =TypeSolverPW
    if (present(IsVariableDt))  IsVariableDt=IsVariableDtPW
    if (present(IsRestart))     IsRestart=IsRestartPW
    if (present(DToutput))      DToutput=DToutputPW

  end subroutine get_field_line

end module ModFieldLine
