!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModAdvectPoints

  use BATL_lib, ONLY: &
       test_start, test_stop

  implicit none

  private ! except

  public:: advect_all_points ! advect all points defined by plots

contains
  !============================================================================

  subroutine advect_all_points

    use ModIO, ONLY: Plot_Type, nFile, Plot_, &
         NameLine_I, nLine_I, XyzStartLine_DII

    integer :: iFile, iPlotFile, nLine

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'advect_all_points'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    call timing_start(NameSub)
    do iFile = Plot_, nFile
       if (index(Plot_Type(iFile),'lin') /= 1) CYCLE
       iPlotFile = iFile - Plot_
       if (NameLine_I(iPlotFile) /= 'A') CYCLE
       nLine = nLine_I(iPlotFile)
       if(DoTest)write(*,*) NameSub,' advect',nLine,' points for file',iFile
       call advect_points(nLine, XyzStartLine_DII(:,1:nLine,iPlotFile))
    end do
    call timing_stop(NameSub)

    call test_stop(NameSub, DoTest)
  end subroutine advect_all_points
  !============================================================================

  subroutine advect_points(nPoint, Xyz_DI)

    use ModMain, ONLY: time_accurate, Dt, nStage

    integer, intent(in)    :: nPoint
    real,    intent(inout) :: Xyz_DI(3,nPoint)

    ! use to be an automatic array
    real, dimension(:,:), allocatable :: XyzOld_DI

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'advect_points'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(.not.time_accurate) RETURN
    call timing_start(NameSub)
    if(nStage == 1)then
       ! Full step uses StateOld
       call advect_points1(0.0, Dt, nPoint, Xyz_DI, Xyz_DI)
    else
       allocate(XyzOld_DI(3,nPoint));
       XyzOld_DI = Xyz_DI
       ! Half step uses StateOld
       call advect_points1(1.0, Dt/2, nPoint, XyzOld_DI, Xyz_DI)
       ! Full step uses (State+StateOld)/2
       call advect_points1(0.5, Dt, nPoint, XyzOld_DI, Xyz_DI)
       deallocate(XyzOld_DI)
    end if
    call timing_stop(NameSub)

    call test_stop(NameSub, DoTest)
  end subroutine advect_points
  !============================================================================

  subroutine advect_points1(WeightOldState, Dt, nPoint, XyzOld_DI, Xyz_DI)

    use ModAdvance, ONLY: Rho_, RhoUz_, Ux_, Uz_, &
         iRho_I, iRhoUx_I, iRhoUy_I, iRhoUz_I
    use ModCurrent, ONLY: get_point_data
    use BATL_lib,   ONLY: nBlock, nProc, iComm
    use ModMpi

    real,    intent(in)    :: WeightOldState      ! Weight of old state
    real,    intent(in)    :: Dt                  ! Time step
    integer, intent(in)    :: nPoint              ! Number of points
    real,    intent(inout) :: XyzOld_DI(3,nPoint) ! Original position
    real,    intent(inout) :: Xyz_DI(3,nPoint)    ! Current position

    integer, parameter :: Weight_= 0
    integer, parameter :: nState = 4 ! Rho, RhoUx, RhoUy, RhoUz

    ! automatic arrays
    real, dimension(:,:), allocatable :: State_VI

    ! Temporary variables
    real    :: Weight
    integer :: iPoint, iError

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'advect_points1'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(DoTest)write(*,*)NameSub,' nPoint=',nPoint
    if(DoTest)write(*,*)NameSub,' old Xyz_DI=',Xyz_DI

    ! Allocate arrays that used to be automatic
    allocate(State_VI(Weight_:nState,nPoint))

    ! Get weight, density and momentum on local PE for all points
    do iPoint = 1, nPoint
       call get_point_data(WeightOldState, Xyz_DI(:,iPoint), 1, nBlock, &
            Rho_, RhoUz_, State_VI(:, iPoint))
    end do

    ! Add up contributions from all the processors
    if(nProc > 1)call MPI_allreduce(MPI_IN_PLACE, State_VI, size(State_VI), &
         MPI_REAL, MPI_SUM, iComm, iError)

    ! Get velocities
    do iPoint = 1, nPoint

       ! Check total weight and divide by it if necessary
       Weight = State_VI(Weight_, iPoint)

       if(Weight <= 0.0)then
          write(*,*)NameSub,' WARNING: iPoint, Weight=',iPoint, Weight
       elseif(abs(Weight - 1) > 1.e-5)then
          write(*,*)NameSub,' WARNING: iPoint, Weight, Xyz=',&
               iPoint, Weight, Xyz_DI(:,iPoint)
          State_VI(1:nState, iPoint) = State_VI(1:nState, iPoint) / Weight
       end if

       ! Convert momenta to velocity
       State_VI(iRhoUx_I,iPoint) = State_VI(iRhoUx_I,iPoint)/State_VI(iRho_I,iPoint)
       State_VI(iRhoUy_I,iPoint) = State_VI(iRhoUy_I,iPoint)/State_VI(iRho_I,iPoint)
       State_VI(iRhoUz_I,iPoint) = State_VI(iRhoUz_I,iPoint)/State_VI(iRho_I,iPoint)

    end do

    ! Move points with the velocity
    Xyz_DI = XyzOld_DI + Dt * State_VI(Ux_:Uz_,:)

!!!
    if(DoTest)then
       write(*,*)NameSub,' Dt        =',Dt
       write(*,*)NameSub,' Velocity  =',State_VI(Ux_:Uz_,:)
       write(*,*)NameSub,' new Xyz_DI=',Xyz_DI(1:3,1)
    end if
    ! Move back along field lines
!!! to be implemented

    ! Deallocate arrays
    deallocate(State_VI)

    call test_stop(NameSub, DoTest)
  end subroutine advect_points1
  !============================================================================

  subroutine advect_test

    ! Create a non-trivial flow field and advect a number of points

    use ModProcMH,   ONLY: iProc
    use ModMain,     ONLY: nI, nJ, nK, nBlock, Dt
    use ModAdvance,  ONLY: StateOld_VGB, State_VGB, Rho_, RhoUx_, RhoUy_, RhoUz_
    use ModGeometry, ONLY: x2
    use BATL_lib,    ONLY: Xyz_DGB, x_, z_
    use ModNumConst, ONLY: cTwoPi
    use ModIoUnit,   ONLY: UnitTmp_
    use ModUtilities, ONLY: open_file, close_file

    integer, parameter :: nStep=100, nPoint = 1000
    real,    parameter :: Ux = 1.0/cTwoPi, Uz = -2.0/cTwoPi
    integer :: iPoint, iStep
    real :: Xyz_DI(3,nPoint)
    real :: Time

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'advect_test'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    call timing_start(NameSub)

    ! Initial positions
    Xyz_DI = 0.0
    do iPoint=1, nPoint
       Xyz_DI(1,iPoint) = (iPoint-0.5)/nPoint*0.9*x2
    end do

    ! Rotation with Omega=1 around the Z axis
    State_VGB(Rho_,:,:,:,1:nBlock)   = 1.0
    State_VGB(RhoUy_,:,:,:,1:nBlock) = 0.0
    State_VGB(RhoUx_,:,:,:,1:nBlock) = -Xyz_DGB(z_,:,:,:,1:nBlock)
    State_VGB(RhoUz_,:,:,:,1:nBlock) = +Xyz_DGB(x_,:,:,:,1:nBlock)

    ! Initial Time
    Time = 0.0

    ! Time step
    Dt = cTwoPi/nStep

    if(iProc==0)then
       call open_file(FILE='advect_test.log')
       write(UnitTmp_,'(a)')'Shifted circle'
       write(UnitTmp_,'(a)')'x1 y1 z1 xM yM zM xN yN zN'
    end if

    do iStep = 1, nStep

       call timing_step(iStep)

       StateOld_VGB(Rho_:RhoUz_,1:nI,1:nJ,1:nK,1:nBlock) = &
            State_VGB(Rho_:RhoUz_,1:nI,1:nJ,1:nK,1:nBlock)

       State_VGB(RhoUx_,:,:,:,1:nBlock) = -(Xyz_DGB(z_,:,:,:,1:nBlock) - Time*Uz)
       State_VGB(RhoUz_,:,:,:,1:nBlock) = +(Xyz_DGB(x_,:,:,:,1:nBlock) - Time*Ux)

       ! Rotation
       call advect_points(nPoint, Xyz_DI)

       Time = Time + Dt

       if(iProc==0)write(UnitTmp_,'(9f15.10)') &
            Xyz_DI(:,1),Xyz_DI(:,nPoint/2), Xyz_DI(:,nPoint)
    end do

    if(iProc==0)call close_file

    call timing_stop(NameSub)

    call timing_report

    call test_stop(NameSub, DoTest)
  end subroutine advect_test
  !============================================================================

end module ModAdvectPoints
!==============================================================================
