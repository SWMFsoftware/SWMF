!^CFG COPYRIGHT UM
!==============================================================================
subroutine advect_all_points

  use ModIO, ONLY: Plot_Type, nFile, Plot_, &
       NameLine_I, nLine_I, XyzStartLine_DII
  implicit none
  integer :: iFile, iPlotFile, nLine
  character(len=*), parameter :: NameSub='advect_all_points'
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  call set_oktest(NameSub, DoTest, DoTestMe)
  do iFile = Plot_, nFile
     if (index(Plot_Type(iFile),'lin') /= 1) CYCLE
     iPlotFile = iFile - Plot_
     if (NameLine_I(iPlotFile) /= 'A') CYCLE
     nLine = nLine_I(iPlotFile)
     if(DoTestMe)write(*,*) NameSub,' advect',nLine,' points for file',iFile
     call advect_points(nLine, XyzStartLine_DII(:,1:nLine,iPlotFile))
  end do

end subroutine advect_all_points
!==============================================================================
subroutine advect_points(nPoint, Xyz_DI)

  use ModMain,    ONLY: time_accurate, nStage, Dt
  use ModAdvance, ONLY: nVar, Rho_, RhoUx_, RhoUz_, Ux_, Uz_
  use ModProcMH,  ONLY: nProc, iComm
  use ModMpi

  implicit none

  integer, intent(in)    :: nPoint
  real,    intent(inout) :: Xyz_DI(3,nPoint)

  integer, parameter :: Weight_= 0
  integer, parameter :: nState = 4 ! Rho, RhoUx, RhoUy, RhoUz

  ! automatic arrays
  real, dimension(Weight_:nState, nPoint) :: State_VI, StateAll_VI

  ! Temporary variables
  real    :: StateMhd_V(nVar)
  real    :: Weight
  integer :: iPoint, iError

  character(len=*), parameter :: NameSub='advect_points'
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  if(.not.time_accurate) RETURN
  call set_oktest(NameSub, DoTest, DoTestMe)

  ! Get weight, density and momentum on local PE for all points
  do iPoint = 1, nPoint
     call get_satellite_data(Xyz_DI(:,iPoint),.false.,StateMhd_V,Weight)
     State_VI(Weight_ , iPoint) = Weight
     State_VI(1:nState, iPoint) = StateMhd_V(Rho_:RhoUz_)
  end do

  ! Add up contributions from all the processors
  if(nProc > 1)then
     call MPI_allreduce(State_VI, StateAll_VI, (1+nState)*nPoint, &
          MPI_REAL, MPI_SUM, iComm, iError)
     State_VI = StateAll_VI
  end if

  ! Get velocities
  do iPoint = 1, nPoint

     ! Check total weight and divide by it if necessary
     Weight = State_VI(Weight_, iPoint)

     if(Weight <= 0.0)then
        write(*,*)'track_lines WARNING: end of line is not found'
     elseif(abs(Weight - 1) > 1.e-5)then
        write(*,*)'track_lines WARNING: first order accurate, iPoint, Xyz=',&
             iPoint, Xyz_DI(:,iPoint)
        State_VI(1:nState, iPoint) = State_VI(1:nState, iPoint) / Weight
     end if

     ! Convert momenta to velocity
     State_VI(Ux_:Uz_, iPoint) = &
          State_VI(RhoUx_:RhoUz_, iPoint)/State_VI(Rho_, iPoint)
  end do

  ! Move points with the velocity
  Xyz_DI = Xyz_DI + Dt * State_VI(Ux_:Uz_,:)

!!!
  if(DoTestMe)write(*,*)NameSub,' new Xyz_DI=',Xyz_DI

  ! Move back along field lines
!!! to be implemented

end subroutine advect_points
