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
  real    :: TimePW,TmaxPW,DToutputPW, DTpolarwindPW,GeoMagLatPW,&
       GeoMagLonPW,JrPW
  real    :: wHorizontalPW
  integer :: iUnitInputPW,      &
       iUnitOutputPW,iUnitGraphicsPW,                 &
       iUnitSourceGraphicsPW,iUnitRestartPW,          &
       iUnitCollisionPW,iUnitRestartInPW,iLinePW,nLinePW
  CHARACTER(7) :: TypeSolverPW

contains

  !******************************************************************************
  !  Put polarwind variables into Mod_PW for passing 
  !******************************************************************************

  subroutine put_field_line(dOxyg, uOxyg, pOxyg, TOxyg,     &
       dHel, uHel, pHel, THel,         &
       dHyd, uHyd, pHyd, THyd,         &
       dElect, uElect, pElect, TElect, &
       IsRestart,IsVariableDt,Time,DT,DToutput,DTpolarwind,&
       TypeSolver,GeoMagLat,GeoMagLon,Jr,wHorizontal,      &
       iUnitOutput,iUnitGraphics,NameRestart,iLine,nLine )

    use ModParameters

    real, intent(in),dimension(maxGrid):: dOxyg, uOxyg, pOxyg, TOxyg,     &
         dHel, uHel, pHel, THel,         &
         dHyd, uHyd, pHyd, THyd,         &
         dElect, uElect, pElect, TElect

    real,    intent(in)     :: Time, DT, DToutput, DTpolarwind,GeoMagLat, &
         GeoMagLon,Jr, wHorizontal                  
    integer, intent(in)     :: iUnitInput,iUnitOutput,iUnitGraphics,      &
         iUnitSourceGraphics,iUnitRestart,          &
         iUnitCollision,iUnitRestartIn,iLine,nLine     
    logical, intent(in)     :: IsRestart,IsVariableDt
    CHARACTER(7), intent(in):: TypeSolver

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
    IsRestartPW = IsRestart
    IsVariableDtPW=IsVariableDT
    TimePW      = Time
    TmaxPW      = Time+DT
    DToutputPW  = DToutput
    DTpolarwindPW = DTpolarwind
    TypeSolverPW  = TypeSolver
    GeoMagLatPW = GeoMagLat
    GeoMagLonPW = GeoMagLon
    JrPW        = Jr
    wHorizontalPW   = wHorizontal
    iUnitOutputPW   =iUnitOutput
    iUnitGraphicsPW =iUnitGraphics
    iUnitSourceGraphicsPW=iUnitSourceGraphics
    iUnitCollisionPW=iUnitCollision
    NameRestartPW=NameRestart
    iLinePW = iLine
    nLinePW = nLine

  end subroutine put_field_line

  !******************************************************************************
  !  Get polarwind variables from Mod_PW 
  !*****************************************************************************


  subroutine get_field_line(dOxyg, uOxyg, pOxyg, TOxyg,     &
       dHel, uHel, pHel, THel,         &
       dHyd, uHyd, pHyd, THyd,         &
       dElect, uElect, pElect, TElect, &
       IsRestart,IsVariableDt,Time,Tmax,DToutput,DTpolarwind,&
       TypeSolver,GeoMagLat,GeoMagLon,Jr,wHorizontal,      &
       iUnitOutput,iUnitGraphics,NameRestart,iLine,nLine )

    use ModParameters
    use ModPass

    real, intent(out),dimension(maxGrid):: dOxyg, uOxyg, pOxyg, TOxyg,     &
         dHel, uHel, pHel, THel,         &
         dHyd, uHyd, pHyd, THyd,         &
         dElect, uElect, pElect, TElect

    real,    intent(out)     :: Time,DToutput, DTpolarwind,Tmax,&
         GeoMagLat,GeoMagLon, Jr, wHorizontal

    integer,intent(out)      :: iUnitInput,iUnitOutput,iUnitGraphics,      &
         iUnitSourceGraphics,iUnitRestart,          &
         iUnitCollision,iUnitRestartIn,iLine,nLine
    logical, intent(out)     :: IsRestart,IsVariableDt
    CHARACTER(7),intent(out) :: TypeSolver


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
    IsRestart = IsRestartPW
    IsVariableDT = IsVariableDtPW
    Time      = TimePW
    Tmax      = TmaxPW
    DToutput  = DToutputPW
    DTpolarwind = DTpolarwindPW
    TypeSolver  = TypeSolverPW
    GeoMagLat = GeoMagLatPW
    GeoMagLon = GeoMagLonPW
    Jr        = JrPW
    wHorizontal = wHorizontalPW
    iUnitInput    =iUnitInputPW
    iUnitOutput   =iUnitOutputPW
    iUnitGraphics =iUnitGraphicsPW
    iUnitSourceGraphics=iUnitSourceGraphicsPW
    iUnitRestart  =iUnitRestartPW
    iUnitCollision=iUnitCollisionPW
    iUnitRestartIn=iUnitRestartInPW
    iLine=iLinePW
    nLine=nLinePW
  end subroutine get_field_line

end module ModFieldLine
