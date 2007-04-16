module ModFieldLine

  use ModParameters
  implicit none

  private
  public put_field_line, get_field_line

  real, Dimension(maxGrid)  ::  dOxygPw_C,uOxygPw_C,pOxygPw_C,TOxygPW, &
       dHelPw_C,uHelPw_C,pHelPw_C,THelPW,     &
       dHydPw_C,uHydPw_C,pHydPw_C,THydPW,     &
       dElectPw_C,uElectPw_C,pElectPw_C,TElectPW,rPw_C

  logical :: IsRestartPW, IsVariableDtPW, DoLogPW,IsImplicitPW,IsImplicitAllPW
  real    :: TimePW,MaxLineTimePW,DToutputPW, DTpolarwindPW,GeoMagLatPW,&
       GeoMagLonPW,JrPW,nDimPW
  real    :: wHorizontalPW
  integer :: nStepPw,iUnitInputPW,      &
       iUnitOutputPW,iUnitGraphicsPW,                 &
       iUnitSourceGraphicsPW,iUnitRestartPW,          &
       iUnitCollisionPW,iUnitRestartInPW,iLinePW,nLinePW
  CHARACTER(7) :: TypeSolverPW
  character*100 :: NameRestartPW

contains

  !***************************************************************************
  !  Put polarwind variables into Mod_PW for passing 
  !***************************************************************************

  subroutine put_field_line(dOxyg_CI, uOxyg_CI, pOxyg_CI, TOxyg,     &
       dHel_CI, uHel_CI, pHel_CI, THel,                              &
       dHyd_CI, uHyd_CI, pHyd_CI, THyd,                              &
       dElect_CI, uElect_CI, pElect_CI, TElect,                      &
       GeoMagLat_I,GeoMagLon_I,Jr,wHorizontal,                  &
       iUnitOutput,iUnitGraphics, NameRestart,iLine,Time,   &
       MaxLineTime,TypeSolver,IsVariableDt,IsRestart,DToutput,nAlt,DoLog,&
       nStep,r_C,IsImplicit,IsImplicitAll)

    use ModParameters

    real, intent(in),dimension(maxGrid):: dOxyg_CI, uOxyg_CI, pOxyg_CI, TOxyg,     &
         dHel_CI, uHel_CI, pHel_CI, THel,         &
         dHyd_CI, uHyd_CI, pHyd_CI, THyd,         &
         dElect_CI, uElect_CI, pElect_CI, TElect
    real, optional, intent(in) :: Time,MaxLineTime,DToutput
    real, optional, intent(in) :: r_C(maxGrid)
    real,    intent(in)     :: GeoMagLat_I,GeoMagLon_I,Jr,wHorizontal              
    integer, optional,intent(in)     :: iUnitOutput,iUnitGraphics,iLine,nAlt,nStep
    character*100,optional,intent(in):: NameRestart
    character(7),optional,intent(in)::TypeSolver
    logical,optional,intent(in) :: IsVariableDt,IsRestart,DoLog,IsImplicit,IsImplicitAll
    !-------------------------------------------------------------------------

    dOxygPw_C (:) = dOxyg_CI(:)
    uOxygPw_C (:) = uOxyg_CI(:)
    pOxygPw_C (:) = pOxyg_CI(:)
    TOxygPW (:) = TOxyg(:)
    dHelPw_C  (:) = dHel_CI(:)
    uHelPw_C  (:) = uHel_CI(:)
    pHelPw_C  (:) = pHel_CI(:)
    THelPW  (:) = THel(:)
    dHydPw_C  (:) = dHyd_CI(:)
    uHydPw_C  (:) = uHyd_CI(:)
    pHydPw_C  (:) = pHyd_CI(:)
    THydPW  (:) = THyd(:)
    dElectPw_C(:) = dElect_CI(:)
    uElectPw_C(:) = uElect_CI(:)
    pElectPw_C(:) = pElect_CI(:)
    TElectPW(:) = TElect(:) 
    GeoMagLatPW = GeoMagLat_I
    GeoMagLonPW = GeoMagLon_I
    JrPW        = Jr
    wHorizontalPW   = wHorizontal
    

    if (present(nAlt))          nDimPW = nAlt
    
    if (present(Time))          TimePW = Time
    if (present(MaxLineTime))   MaxLineTimePW  = MaxLineTime
    if (present(iUnitGraphics)) iUnitGraphicsPW= iUnitGraphics
    if (present(NameRestart))   NameRestartPW  = NameRestart
    if (present(iLine))         iLinePW        = iLine
    if (present(iUnitOutput))   iUnitOutputPW  = iUnitOutput
    if (present(TypeSolver))    TypeSolverPW   = TypeSolver
    if (present(IsVariableDt))  IsVariableDtPW = IsVariableDt
    if (present(IsRestart))     IsRestartPW    = IsRestart
    if (present(DToutput))      DToutputPW     = DToutput
    if (present(DoLog))         DoLogPW        = DoLog
    if (present(nStep))         nStepPW        = nStep
    if (present(r_c))           rPw_C(:)       = r_C(:)
    if (present(IsImplicit))    IsImplicitPW   = IsImplicit
    if (present(IsImplicitAll)) IsImplicitAllPW   = IsImplicitAll


  end subroutine put_field_line

  !***************************************************************************
  !  Get polarwind variables from Mod_PW 
  !***************************************************************************

  subroutine get_field_line(dOxyg_CI, uOxyg_CI, pOxyg_CI, TOxyg,     &
       dHel_CI, uHel_CI, pHel_CI, THel,                              &
       dHyd_CI, uHyd_CI, pHyd_CI, THyd,                              &
       dElect_CI, uElect_CI, pElect_CI, TElect,                      &
       GeoMagLat_I,GeoMagLon_I,Jr,wHorizontal,                  &
       iUnitOutput,iUnitGraphics, NameRestart,iLine,Time,   &
       MaxLineTime,TypeSolver,IsVariableDt,IsRestart,DToutput, nAlt,DoLog,&
       nStep,r_C,IsImplicit,IsImplicitAll)

    use ModParameters

    real, intent(out),dimension(maxGrid):: dOxyg_CI, uOxyg_CI, pOxyg_CI, TOxyg,     &
         dHel_CI, uHel_CI, pHel_CI, THel,         &
         dHyd_CI, uHyd_CI, pHyd_CI, THyd,         &
         dElect_CI, uElect_CI, pElect_CI, TElect


    real,    intent(out)     :: GeoMagLat_I, GeoMagLon_I,Jr, wHorizontal           
    
    character*100,optional,intent(out):: NameRestart
    character(7),optional,intent(out):: TypeSolver
    logical,optional,intent(out)      :: IsVariableDt,IsRestart,DoLog,IsImplicit,IsImplicitAll
    real, optional, intent(out)       :: Time,MaxLineTime,DToutput
    real, optional, intent(out)        :: r_C(maxGrid)
    integer, optional,intent(out)     :: iUnitOutput,iUnitGraphics,iLine,nAlt,nStep
    !--------------------------------------------------------------------------
    
    dOxyg_CI (:) = dOxygPw_C(:)
    uOxyg_CI (:) = uOxygPw_C(:)
    pOxyg_CI (:) = pOxygPw_C(:)
    TOxyg (:) = TOxygPW(:)
    dHel_CI  (:) = dHelPw_C(:)
    uHel_CI  (:) = uHelPw_C(:)
    pHel_CI  (:) = pHelPw_C(:)
    THel  (:) = THelPW(:)
    dHyd_CI  (:) = dHydPw_C(:)
    uHyd_CI  (:) = uHydPw_C(:)
    pHyd_CI  (:) = pHydPw_C(:)
    THyd  (:) = THydPW(:)
    dElect_CI(:) = dElectPw_C(:)
    uElect_CI(:) = uElectPw_C(:)
    pElect_CI(:) = pElectPw_C(:)
    TElect(:) = TElectPW(:) 
    GeoMagLat_I = GeoMagLatPW
    GeoMagLon_I = GeoMagLonPW
    Jr        = JrPW
    wHorizontal = wHorizontalPW
   

    

    if (present(nAlt))          nAlt = nDimPW
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
    if (present(DoLog))         DoLog=DoLogPW
    if (present(nStep))         nStep        = nStepPW
    if (present(r_C))           r_C(:)       = rPw_C(:)
    if (present(IsImplicit))    IsImplicit   = IsImplicitPW
    if (present(IsImplicitAll)) IsImplicitAll= IsImplicitAllPW
  end subroutine get_field_line

end module ModFieldLine
