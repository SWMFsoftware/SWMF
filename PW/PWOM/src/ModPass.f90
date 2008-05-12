module ModFieldLine

  use ModParameters
  use ModCommonPlanet,ONLY: nVar
  implicit none

  private
  public put_field_line, get_field_line

  real,allocatable    ::  PassState_CV(:,:)

  logical :: DoLogPW
  real    :: TimePW,MaxLineTimePW,DToutputPW, DtPw,GeoMagLatPW,&
       GeoMagLonPW,JrPW,nDimPW,uJoule2PW
  real    :: wHorizontalPW
  integer :: nStepPw,iUnitInputPW,      &
       iUnitOutputPW,                 &
       iUnitSourceGraphicsPW,iUnitRestartPW,          &
       iUnitCollisionPW,iUnitRestartInPW,iLinePW,nLinePW
  CHARACTER(7) :: TypeSolverPW
  character*100 :: NameRestartPW
  real,allocatable  :: rPW_C(:)
contains

  !***************************************************************************
  !  Put polarwind variables into Mod_PW for passing 
  !***************************************************************************

  subroutine put_field_line(nAlt,State_CV, &
       GeoMagLat_I,GeoMagLon_I,Jr,wHorizontal,uJoule2,                  &
       iUnitOutput, NameRestart,iLine,Time,   &
       MaxLineTime,TypeSolver,DToutput,DoLog,&
       nStep,r_C,Dt)

    use ModParameters
    use ModCommonPlanet,ONLY: nVar
    integer, intent(in) :: nAlt
    real, intent(in) :: State_CV(nAlt,nVar)
    real, optional, intent(in) :: Time,MaxLineTime,DToutput,uJoule2
    real, optional, intent(in) :: r_C(nAlt)
    real,    intent(in)     :: GeoMagLat_I,GeoMagLon_I,Jr,wHorizontal  
    integer, optional,intent(in)  :: iUnitOutput,iLine,nStep
    character*100,optional,intent(in):: NameRestart
    character(7),optional,intent(in)::TypeSolver
    logical,optional,intent(in) :: DoLog
    real,optional,intent(in) :: Dt
    !-------------------------------------------------------------------------
    if (.not. allocated(PassState_CV)) allocate(PassState_CV(nAlt,nVar))
    if (.not. allocated(rPW_C))        allocate(rPW_C(nAlt))
    PassState_CV(:,:) = State_CV(:,:)
    GeoMagLatPW = GeoMagLat_I
    GeoMagLonPW = GeoMagLon_I
    JrPW        = Jr
    wHorizontalPW   = wHorizontal
    

    if (present(Time))          TimePW = Time
    if (present(MaxLineTime))   MaxLineTimePW  = MaxLineTime
    if (present(NameRestart))   NameRestartPW  = NameRestart
    if (present(iLine))         iLinePW        = iLine
    if (present(iUnitOutput))   iUnitOutputPW  = iUnitOutput
    if (present(TypeSolver))    TypeSolverPW   = TypeSolver
    if (present(DToutput))      DToutputPW     = DToutput
    if (present(DoLog))         DoLogPW        = DoLog
    if (present(nStep))         nStepPW        = nStep
    if (present(r_c))           rPw_C(:)       = r_C(:)
    if (present(uJoule2))       uJoule2PW      = uJoule2
    if (present(Dt))            DtPw           = Dt

  end subroutine put_field_line

  !***************************************************************************
  !  Get polarwind variables from Mod_PW 
  !***************************************************************************

  subroutine get_field_line(nAlt,State_CV,& 
       GeoMagLat_I,GeoMagLon_I,Jr,wHorizontal,uJoule2,                  &
       iUnitOutput, NameRestart,iLine,Time,   &
       MaxLineTime,TypeSolver,DToutput,DoLog,&
       nStep,r_C,Dt)

    use ModParameters
    use ModCommonPlanet, ONLY: nVar
    integer,intent(in) :: nAlt
    real, intent(out):: State_CV(nAlt,nVar)


    real,    intent(out)     :: GeoMagLat_I, GeoMagLon_I,Jr, wHorizontal           
    
    character*100,optional,intent(out):: NameRestart
    character(7),optional,intent(out):: TypeSolver
    logical,optional,intent(out)      :: DoLog
    real, optional, intent(out)       :: Time,MaxLineTime,DToutput,uJoule2
    real, optional, intent(out)        :: r_C(nAlt),Dt
    integer, optional,intent(out)     :: iUnitOutput,iLine,nStep
    !--------------------------------------------------------------------------
    if (.not. allocated(PassState_CV)) allocate(PassState_CV(nAlt,nVar))
    if (.not. allocated(rPW_C)) allocate(rPW_C(nAlt))
    State_CV(:,:)= PassState_CV(:,:)
    GeoMagLat_I = GeoMagLatPW
    GeoMagLon_I = GeoMagLonPW
    Jr        = JrPW
    wHorizontal = wHorizontalPW
   

    

    if (present(Time))          Time = TimePW
    if (present(MaxLineTime))   MaxLineTime = MaxLineTimePW
    if (present(NameRestart))   NameRestart=NameRestartPW
    if (present(iLine))         iLine = iLinePW
    if (present(iUnitOutput))   iUnitOutput=iUnitOutputPW
    if (present(TypeSolver))    TypeSolver =TypeSolverPW
    if (present(DToutput))      DToutput=DToutputPW
    if (present(DoLog))         DoLog=DoLogPW
    if (present(nStep))         nStep        = nStepPW
    if (present(r_C))           r_C(:)       = rPw_C(:)
    if (present(uJoule2))       uJoule2      = uJoule2PW
    if (present(Dt))            Dt           = DtPw
  end subroutine get_field_line

end module ModFieldLine
