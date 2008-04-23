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
  call timing_start(NameSub)
  do iFile = Plot_, nFile
     if (index(Plot_Type(iFile),'lin') /= 1) CYCLE
     iPlotFile = iFile - Plot_
     if (NameLine_I(iPlotFile) /= 'A') CYCLE
     nLine = nLine_I(iPlotFile)
     if(DoTestMe)write(*,*) NameSub,' advect',nLine,' points for file',iFile
     call advect_points(nLine, XyzStartLine_DII(:,1:nLine,iPlotFile))
  end do
  call timing_stop(NameSub)

end subroutine advect_all_points

!==============================================================================
subroutine advect_points(nPoint, Xyz_DI)

  use ModMain, ONLY: time_accurate, Dt, nStage
  
  implicit none

  integer, intent(in)    :: nPoint
  real,    intent(inout) :: Xyz_DI(3,nPoint)

  ! use to be an automatic array
  real, dimension(:,:), allocatable :: XyzOld_DI
  integer :: iError

  character(len=*), parameter :: NameSub = 'advect_points'
  !--------------------------------------------------------------------------
  if(.not.time_accurate) RETURN
  call timing_start(NameSub)
  if(nStage == 1)then
     ! Full step uses StateOld
     call advect_points1(0.0, Dt, nPoint, Xyz_DI, Xyz_DI)
  else
     allocate(XyzOld_DI(3,nPoint), stat=iError); 
     call alloc_check(iError,"XyzOld_DI")
     XyzOld_DI = Xyz_DI
     ! Half step uses StateOld
     call advect_points1(1.0, Dt/2, nPoint, XyzOld_DI, Xyz_DI)
     ! Full step uses (State+StateOld)/2
     call advect_points1(0.5, Dt, nPoint, XyzOld_DI, Xyz_DI)
     deallocate(XyzOld_DI)
  end if
  call timing_stop(NameSub)

end subroutine advect_points

!==============================================================================
subroutine advect_points1(WeightOldState, Dt, nPoint, XyzOld_DI, Xyz_DI)

  use ModMain,    ONLY: nBlock
  use ModAdvance, ONLY: nVar, Rho_, RhoUx_, RhoUz_, Ux_, Uz_
  use ModProcMH,  ONLY: nProc, iComm
  use ModMpi

  implicit none

  real,    intent(in)    :: WeightOldState      ! Weight of old state
  real,    intent(in)    :: Dt                  ! Time step
  integer, intent(in)    :: nPoint              ! Number of points
  real,    intent(inout) :: XyzOld_DI(3,nPoint) ! Original position
  real,    intent(inout) :: Xyz_DI(3,nPoint)    ! Current position

  integer, parameter :: Weight_= 0
  integer, parameter :: nState = 4 ! Rho, RhoUx, RhoUy, RhoUz

  ! automatic arrays
!!$  real, dimension(Weight_:nState, nPoint) :: State_VI, StateAll_VI
  real, dimension(:,:), allocatable :: State_VI, StateAll_VI

  ! Temporary variables
  real    :: Weight
  integer :: iPoint, iError

  character(len=*), parameter :: NameSub='advect_points1'
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  call set_oktest(NameSub, DoTest, DoTestMe)

  if(DoTestMe)write(*,*)NameSub,' nPoint=',nPoint
  if(DoTestMe)write(*,*)NameSub,' old Xyz_DI=',Xyz_DI

  ! Allocate arrays that used to be automatic
  allocate(State_VI(Weight_:nState,nPoint), &
       StateAll_VI(Weight_:nState,nPoint), stat=iError)
  call alloc_check(iError,"advect_points1 arrays")

  ! Get weight, density and momentum on local PE for all points
  do iPoint = 1, nPoint
     call get_point_data(WeightOldState, Xyz_DI(:,iPoint), 1, nBlock, &
          Rho_, RhoUz_, State_VI(:, iPoint))
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
        write(*,*)NameSub,' WARNING: iPoint, Weight=',iPoint, Weight
     elseif(abs(Weight - 1) > 1.e-5)then
        write(*,*)NameSub,' WARNING: iPoint, Weight, Xyz=',&
             iPoint, Weight, Xyz_DI(:,iPoint)
        State_VI(1:nState, iPoint) = State_VI(1:nState, iPoint) / Weight
     end if

     ! Convert momenta to velocity
     State_VI(Ux_:Uz_, iPoint) = &
          State_VI(RhoUx_:RhoUz_, iPoint)/State_VI(Rho_, iPoint)

  end do

  ! Move points with the velocity
  Xyz_DI = XyzOld_DI + Dt * State_VI(Ux_:Uz_,:)

!!!
  if(DoTestMe)then
     write(*,*)NameSub,' Dt        =',Dt
     write(*,*)NameSub,' Velocity  =',State_VI(Ux_:Uz_,:)
     write(*,*)NameSub,' new Xyz_DI=',Xyz_DI(1:3,1)
  end if
  ! Move back along field lines
!!! to be implemented

  ! Deallocate arrays that used to be utomatic
  deallocate(State_VI, StateAll_VI)

end subroutine advect_points1

!==============================================================================

subroutine get_point_data(WeightOldState, XyzIn_D, iBlockMin, iBlockMax, &
     iVarMin, iVarMax, StateCurrent_V)

  ! Interpolate the (new and/or old) state vector from iVarMin to iVarMax and 
  ! the current (if iVarMax=nVar+3) for input position 
  ! XyzIn_D given in Cartesian coordinates. The interpolated state 
  ! is second order accurate everywhere except where there is a 
  ! resolution change in more than one direction for the cell centers 
  ! surrounding the given position. In these exceptional cases the 
  ! interpolated state is first order accurate. The interpolation algorithm
  ! is based on trilinear interpolation, but it is generalized for
  ! trapezoidal hexahedrons.

  use ModNumConst
  use ModVarIndexes, ONLY : nVar, Bx_, By_, Bz_
  use ModProcMH
  use ModMain, ONLY : nI, nJ, nK, nCells, nBlock, unusedBLK, BlkTest
  use ModAdvance, ONLY : State_VGB, StateOld_VCB
  use ModGeometry, ONLY : XyzStart_BLK, dx_BLK, dy_BLK, dz_BLK
  use ModGeometry, ONLY : UseCovariant     
  use ModParallel, ONLY : NeiLev

  implicit none

  ! Weight of the old state in the interpolation
  real, intent(in)  :: WeightOldState

  ! Input position is in generalized coordinates
  real, intent(in)  :: XyzIn_D(3)

  ! Block index range (typically 1:nBlock or iBlock:iBlock)
  integer, intent(in) :: iBlockMin, iBlockMax

  ! Do we need to calculate currents
  integer, intent(in) :: iVarMin, iVarMax

  ! Weight and interpolated state at the input position
  real, intent(out) :: StateCurrent_V(0:iVarMax-iVarMin+1)

  ! Local variables

  ! Maximum index for state variable and number of state variables
  integer :: iStateMax, nState

  ! Shall we calculate the current
  logical :: DoCurrent

  ! Position in generalized coordinates
  real :: Xyz_D(3)

  ! Cell size and buffer size for current block
  real,    dimension(3) :: Dxyz_D, DxyzInv_D, DxyzLo_D, DxyzHi_D

  ! Position of cell center to the lower index direction
  integer, dimension(3) :: IjkLo_D 

  ! Position of point and current cell center
  real :: x, y, z, xI, yJ, zK

  ! Bilinear weights
  real    :: WeightX,WeightY,WeightZ,WeightXyz

  ! Dimension, cell, block index and MPI error code
  integer :: iDim,i,j,k,iLo,jLo,kLo,iHi,jHi,kHi,iBlock,iError

  ! Testing
  logical :: DoTest,DoTestMe
  !----------------------------------------------------------------------------
  ! Calculate maximum index and the number of state variables
  iStateMax = min(iVarMax, nVar)
  nState    = iStateMax - iVarMin + 1

  ! Convert to generalized coordinates if necessary

  if(UseCovariant)then                
     call xyz_to_gen(XyzIn_D,Xyz_D)
  else                                
     Xyz_D = XyzIn_D
  end if                              
 

  ! Set state and weight to zero, so MPI_reduce will add it up right
  StateCurrent_V = cZero

  ! Loop through all blocks
  BLOCK: do iBlock = iBlockMin, iBlockMax
     if(unusedBLK(iBlock)) CYCLE

!     if(iBlock == BlkTest)then
!        call set_oktest('get_point_data', DoTest, DoTestMe)
!     else
        DoTest = .false.; DoTestMe = .false.
!     end if

     if(DoTestMe)write(*,*)'get_point_data called with XyzIn_D=',XyzIn_D

     ! Put cell size of current block into an array
     Dxyz_D(1)=Dx_BLK(iBlock)
     Dxyz_D(2)=Dy_BLK(iBlock)
     Dxyz_D(3)=Dz_BLK(iBlock)

     ! Set buffer zone according to relative size of neighboring block
     do iDim = 1, 3
        ! Block at the lower index side
        select case(NeiLev(2*iDim-1,iBlock))
        case(1)
           DxyzLo_D(iDim) = 1.5*Dxyz_D(iDim)
        case(-1)
           DxyzLo_D(iDim) = 0.75*Dxyz_D(iDim)
        case default
           DxyzLo_D(iDim) = Dxyz_D(iDim)
        end select
        ! Check if point is inside the buffer zone on the lower side
        if(Xyz_D(iDim)<XyzStart_BLK(iDim,iBlock) - DxyzLo_D(iDim)) CYCLE BLOCK

        ! Block at the upper index side
        select case(NeiLev(2*iDim,iBlock))
        case(1)
           DxyzHi_D(iDim) = 1.5*Dxyz_D(iDim)
        case(-1)
           DxyzHi_D(iDim) = 0.75*Dxyz_D(iDim)
        case default
           DxyzHi_D(iDim) = Dxyz_D(iDim)
        end select
        ! Check if point is inside the buffer zone on the upper side
        if(Xyz_D(iDim) > XyzStart_BLK(iDim,iBlock) + &
             (nCells(iDim)-1)*Dxyz_D(iDim) + DxyzHi_D(iDim)) CYCLE BLOCK
     end do

     ! Find closest cell center indexes towards the lower index direction
     IjkLo_D = floor((Xyz_D - XyzStart_BLK(:,iBlock))/Dxyz_D)+1

     ! Set the size of the box for bilinear interpolation

     ! At resolution change the box size is the sum
     ! average of the cell size of the neighboring blocks

     ! Also make sure that IjkLo_D is not out of bounds
     do iDim = 1,3
        if(IjkLo_D(iDim) < 1)then
           IjkLo_D(iDim)   = 0
           DxyzInv_D(iDim) = 1/DxyzLo_D(iDim)
        elseif(IjkLo_D(iDim) >= nCells(iDim))then
           IjkLo_D(iDim)   = nCells(iDim)
           DxyzInv_D(iDim) = 1/DxyzHi_D(iDim)
        else
           DxyzInv_D(iDim) = 1/Dxyz_D(iDim)
        end if
     end do

     if(DoTest)then
        write(*,*)'Point found at iProc,iBlock,iLo,jLo,kLo=',&
             iProc,iBlock,IjkLo_D
        write(*,*)'iProc, XyzStart_BLK,Dx_BLK=',iProc, &
             XyzStart_BLK(:,iBlock),dx_BLK(iBlock)
     end if

     ! Set the index range for the physical cells
     iLo = max(IjkLo_D(1),1)
     jLo = max(IjkLo_D(2),1)
     kLo = max(IjkLo_D(3),1)
     iHi = min(IjkLo_D(1)+1,nI)
     jHi = min(IjkLo_D(2)+1,nJ)
     kHi = min(IjkLo_D(3)+1,nK)

     ! Put the point position into scalars
     x = Xyz_D(1); y = Xyz_D(2); z = Xyz_D(3)

     ! Loop through the physical cells to add up their contribution
     do k = kLo, kHi
        zK = XyzStart_BLK(3,iBlock) + (k-1)*Dxyz_D(3)
        WeightZ = 1 - DxyzInv_D(3)*abs(z-zK)
        do j = jLo, jHi
           yJ = XyzStart_BLK(2,iBlock) + (j-1)*Dxyz_D(2)
           WeightY = 1 - DxyzInv_D(2)*abs(y-yJ)
           do i = iLo, iHi
              xI = XyzStart_BLK(1,iBlock) + (i-1)*Dxyz_D(1)
              WeightX = 1 - DxyzInv_D(1)*abs(x-xI)

              WeightXyz = WeightX*WeightY*WeightZ

              if(WeightXyz>0.0)then
                 ! Add up the weight
                 StateCurrent_V(0) = StateCurrent_V(0) + WeightXyz

                 ! Add contibutions from the old state if required
                 if(WeightOldState > 0.0) &
                      StateCurrent_V(1:nState) = StateCurrent_V(1:nState) + &
                      WeightXyz * WeightOldState * &
                      StateOld_VCB(iVarMin:iStateMax,i,j,k,iBlock)

                 ! Add contibutions from the current state if required
                 if(WeightOldState < 1.0) &
                      StateCurrent_V(1:nState) = StateCurrent_V(1:nState) + &
                      WeightXyz * (1 - WeightOldState) * &
                      State_VGB(iVarMin:iStateMax,i,j,k,iBlock)

                 ! The current is always based on the new state
                 if(iVarMax == nVar + 3)call add_current

                 if(DoTest)write(*,*)'Contribution iProc,i,j,k,WeightXyz=',&
                      iProc,i,j,k,WeightXyz
              end if
           end do
        end do
     end do
  end do BLOCK

contains

  !============================================================================
  subroutine add_current

    ! Add current to the current part of StateCurrent
    use ModMain, ONLY: prolong_order
    use ModGeometry,ONLY:body_BLK
    real :: Current_D(3)
    real :: DxInv, DyInv, DzInv
    integer :: iLo,iHi,jLo,jHi,kLo,kHi

    iLo=i-1; jLo=j-1; kLo=k-1; iHi=i+1; jHi=j+1; kHi=k+1

    if(prolong_order==1)then
       ! Avoid the ghost cells at resolution changes
       if(i==1 .and.DxyzLo_D(1)/=Dxyz_D(1))iLo=1
       if(i==nI.and.DxyzHi_D(1)/=Dxyz_D(1))iHi=nI
       if(j==1 .and.DxyzLo_D(2)/=Dxyz_D(2))jLo=1
       if(j==nJ.and.DxyzHi_D(2)/=Dxyz_D(2))jHi=nJ
       if(k==1 .and.DxyzLo_D(3)/=Dxyz_D(3))kLo=1
       if(k==nK.and.DxyzHi_D(3)/=Dxyz_D(3))kHi=nK
    end if

    DxInv = 1/((iHi-iLo)*Dxyz_D(1))
    DyInv = 1/((jHi-jLo)*Dxyz_D(2))
    DzInv = 1/((kHi-kLo)*Dxyz_D(3))

    if(UseCovariant)then                          
       call covariant_curlb(i,j,k,iBlock,Current_D,.not.body_BLK(iBlock))
    else                                          
       Current_D(1) = &
            (State_VGB(Bz_,i,jHi,k,iBlock)-State_VGB(Bz_,i,jLo,k,iBlock))*DyInv- &
            (State_VGB(By_,i,j,kHi,iBlock)-State_VGB(By_,i,j,kLo,iBlock))*DzInv

       Current_D(2) = &
            (State_VGB(Bx_,i,j,kHi,iBlock)-State_VGB(Bx_,i,j,kLo,iBlock))*DzInv- &
            (State_VGB(Bz_,iHi,j,k,iBlock)-State_VGB(Bz_,iLo,j,k,iBlock))*DxInv

       Current_D(3) = &
            (State_VGB(By_,iHi,j,k,iBlock)-State_VGB(By_,iLo,j,k,iBlock))*DxInv- &
            (State_VGB(Bx_,i,jHi,k,iBlock)-State_VGB(Bx_,i,jLo,k,iBlock))*DyInv

    end if                                       
    

    StateCurrent_V(nState+1:nState+3) = StateCurrent_V(nState+1:nState+3) &
         + WeightXyz * Current_D

  end subroutine add_current

end subroutine get_point_data

!==============================================================================
subroutine advect_test

  ! Create a non-trivial flow field and advect a number of points

  use ModProcMH,   ONLY: iProc
  use ModMain,     ONLY: nI, nJ, nK, nBlock, Dt
  use ModAdvance,  ONLY: StateOld_VCB, State_VGB, Rho_, RhoUx_, RhoUy_, RhoUz_
  use ModGeometry, ONLY: x2, x_BLK, y_BLK, z_BLK
  use ModNumConst, ONLY: cTwoPi
  use ModIoUnit,   ONLY: UnitTmp_
  implicit none
  integer, parameter :: nStep=100, nPoint = 1000
  real,    parameter :: Ux = 1.0/cTwoPi, Uz = -2.0/cTwoPi
  integer :: iPoint, iStep
  real :: Xyz_DI(3,nPoint)
  real :: Time

  character(len=*), parameter :: NameSub = 'advect_test'
  !-------------------------------------------------------------------------

  call timing_start(NameSub)

  ! Initial positions
  Xyz_DI = 0.0
  do iPoint=1, nPoint
     Xyz_DI(1,iPoint) = (iPoint-0.5)/nPoint*0.9*x2
  end do

  ! Rotation with Omega=1 around the Z axis
  State_VGB(Rho_,:,:,:,1:nBlock)   = 1.0
  State_VGB(RhoUy_,:,:,:,1:nBlock) = 0.0
  State_VGB(RhoUx_,:,:,:,1:nBlock) = -z_BLK(:,:,:,1:nBlock)
  State_VGB(RhoUz_,:,:,:,1:nBlock) = +x_BLK(:,:,:,1:nBlock)

  ! Initial Time
  Time = 0.0

  ! Time step
  Dt = cTwoPi/nStep

  if(iProc==0)then
     open(UnitTmp_,file='advect_test.log')
     write(UnitTmp_,'(a)')'Shifted circle'
     write(UnitTmp_,'(a)')'x1 y1 z1 xM yM zM xN yN zN'
  end if

  do iStep = 1, nStep

     call timing_step(iStep)

     StateOld_VCB(Rho_:RhoUz_,:,:,:,1:nBlock) = &
          State_VGB(Rho_:RhoUz_,1:nI,1:nJ,1:nK,1:nBlock)

     State_VGB(RhoUx_,:,:,:,1:nBlock) = -(z_BLK(:,:,:,1:nBlock) - Time*Uz)
     State_VGB(RhoUz_,:,:,:,1:nBlock) = +(x_BLK(:,:,:,1:nBlock) - Time*Ux)

     ! Rotation
     call advect_points(nPoint, Xyz_DI)

     Time = Time + Dt

     if(iProc==0)write(UnitTmp_,'(9f15.10)') &
          Xyz_DI(:,1),Xyz_DI(:,nPoint/2), Xyz_DI(:,nPoint)
  end do

  if(iProc==0)close(UnitTmp_)

  call timing_stop(NameSub)

  call timing_report

end subroutine advect_test
