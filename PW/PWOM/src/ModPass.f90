Module ModPass
  use ModParameters
  
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


end Module ModPass
