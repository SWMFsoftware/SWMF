!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program advect

  use BATL_lib, ONLY: MaxDim, nDim, nDimAmr, nI, nJ, nK, &
       MinI, MaxI, MinJ, MaxJ, MinK, MaxK, iProc, barrier_mpi, r_, &
       Xyz_DGB, CellSize_DB, message_pass_cell, &
       StringTest, &
       XyzTest_D, iTest, jTest, kTest, iBlockTest, iProcTest, &
       XyzTest2_D, iTest2, jTest2, kTest2, iBlockTest2, iProcTest2, &
       read_test_param, find_test_cell, test_start, test_stop

  implicit none

  ! Geometry
  character(len=20):: TypeGeometry = 'cartesian'

  ! Square of the radius of the sphere
  real :: BlobRadius2 = 25.0
  real :: BlobCenter_D(nDim) = 0.0

  ! Maximum refinement level
  integer:: MaxLevel = 3

  ! Dynamic adaptation
  logical:: DoAdapt = .true.

  ! Initial refinement radius
  real:: RefineRadius2 = -1.0

  ! Refinement criteria
  real :: RhoCoarsen = 1.4
  real :: RhoRefine  = 1.5

  ! Final simulation time, frequency of plots
  real :: TimeMax = 20.0
  real :: TimeStartLocalStep = 0.49999*20.0
  real :: DtPlot = 1.0

  ! Velocity field
  logical:: UseConstantVelocity = .false.
  real :: Velocity_D(nDim) = 0.0
  real :: RadialVelocity   = 0.0
  real :: AngularVelocity  = 0.0

  ! Linear field. Initially f = sum(Linear_D*Xyz_D)
  real :: Linear_D(nDim) = 1.0

  ! Spatial order of accuracy and beta parameter for the TVD limiter
  integer :: nOrder = 2
  real    :: BetaLimiter = 1.5 ! 1 <= Beta <= 2 for nOrder=2
  logical :: DoLimitVolumeState = .false.

  ! Use fixed or local time stepping
  logical :: UseLocalStep = .false.

  ! Parameters for the explicit scheme
  real :: Cfl = 0.8

  ! Size of the computational domain
  real :: &
       DomainMax_D(3) = (/ +10.0, +5.0, +5.0 /), &
       DomainMin_D(3) = (/ -10.0, -5.0, -5.0 /)

  ! Periodicity
  logical:: IsPeriodicGrid_D(MaxDim) = .true.

  ! Time step counter and simulation time
  integer :: iStep
  real :: Time, TimePlot, Dt

  integer, parameter:: nVar = 2, Rho_ = 1, Lin_=2

  ! Cell centered state
  real, allocatable :: State_VGB(:,:,:,:,:), StateOld_VCB(:,:,:,:,:)

  ! Face centered flux for one block
  real:: Flux_VFD(nVar,1:nI+1,1:nJ+1,1:nK+1,nDim)

  ! Normal velocity times area
  real:: vFaceNormal_DF(nDim,1:nI+1,1:nJ+1,1:nK+1)

  ! Face centered flux for conservation fix
  real, allocatable, dimension(:,:,:,:,:):: Flux_VXB, Flux_VYB, Flux_VZB

  logical:: DoTest
  character(len=*), parameter:: NameSub = 'advect_main'
  !--------------------------------------------------------------------------
  call initialize

  call test_start(NameSub, DoTest, DoTestAll=.true.)
  if(DoTest)call barrier_mpi

  call timing_start('ADVECT')
  do
     if(DoTest)write(*,*)NameSub,' advance iProc, iStep=',iProc,iStep

     call find_test_cell

     if(iProc==iProcTest .and. StringTest /= '')then
        write(*,*)'Test cell iTest, jTest, kTest, iBlockTest, iProcTest=', &
             iTest, jTest, kTest, iBlockTest, iProcTest
        write(*,*)'Test cell Xyz_D=', XyzTest_D(1:nDim)
        write(*,*)'Test cell size =', CellSize_DB(1:nDim,iBlockTest)
     end if

     if(iProc==iProcTest2 .and. StringTest /= '')then
        write(*,*)'2nd Test cell iTest, jTest, kTest, iBlockTest, iProcTest=',&
             iTest2, jTest2, kTest2, iBlockTest2, iProcTest2
        write(*,*)'2nd Test cell Xyz_D=', XyzTest2_D(1:nDim)
        write(*,*)'2nd Test cell size =', CellSize_DB(1:nDim,iBlockTest2)
     end if

     ! Save plot at required frequency
     if( Time >= TimePlot - 1e-10 )then

        if(DoTest)write(*,*)NameSub,' saving plots iProc=',iProc

        call timing_start('save_plot')
        call save_plot
        call timing_stop('save_plot')

        if(DoTest)write(*,*)NameSub,' save logfile iProc=',iProc

        call timing_start('save_log')
        call save_log
        call timing_stop('save_log')

        if(DoTest)write(*,*)NameSub,' finished plots iProc=',iProc

        TimePlot = TimePlot + DtPlot
     end if

     if(DoTest)write(*,*)NameSub,' advance_expl iProc=',iProc
     if(DoTest)call barrier_mpi

     call timing_start('explicit')
     UseLocalStep = (Time > TimeStartLocalStep .and. nDimAmr == nDim)
     if(UseLocalStep)then
        call advance_localstep
     else
        call advance_explicit
     end if
     call timing_stop('explicit')

     if(DoTest)write(*,*)NameSub,' update time iProc=',iProc
     if(DoTest)call barrier_mpi

     ! Update time
     iStep = iStep + 1
     Time  = Time + Dt

     if(DoAdapt)then
        if(DoTest)write(*,*)NameSub,' adapt grid iProc=',iProc
        if(DoTest)call barrier_mpi

        call timing_start('message_pass')
        call message_pass_cell(nVar, State_VGB)
        call timing_stop('message_pass')

        call set_boundary(Time)

        call timing_start('amr')
        call adapt_grid
        call timing_stop('amr')
     end if

     call timing_step(iStep)

     if(Time >= TimeMax - 1e-10) EXIT

  end do
  call timing_stop('ADVECT')
  call timing_report_total

  if(DoTest)write(*,*)NameSub,' finalize iProc=',iProc

  call finalize

  call test_stop(NameSub, DoTest)

contains

  !===========================================================================
  subroutine adapt_grid

    use BATL_lib, ONLY: regrid_batl, nBlock, Unused_B, iNode_B, &
         iStatusNew_A, Refine_, Coarsen_, IsRzGeometry, Xyz_DGB

    real :: Factor
    integer:: iBlock, i, j, k
    !------------------------------------------------------------------------
    Factor = 1.0
    BLOCK: do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       ! Check refinement first
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(IsRzGeometry) Factor = 1./Xyz_DGB(r_,i,j,k,iBlock)
          if(State_VGB(Rho_,i,j,k,iBlock) > Factor*RhoRefine) then
             iStatusNew_A(iNode_B(iBlock)) = Refine_
             CYCLE BLOCK
          end if
       end do; end do; end do
       ! If not refined, check for coarsening
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(IsRzGeometry) Factor = 1./Xyz_DGB(r_,i,j,k,iBlock)
          if(State_VGB(Rho_,i,j,k,iBlock) < Factor*RhoCoarsen)then
             iStatusNew_A(iNode_B(iBlock)) = Coarsen_
             CYCLE BLOCK
          end if
       end do; end do; end do
    end do BLOCK

    call regrid_batl(nVar, State_VGB, DoBalanceEachLevelIn=UseLocalStep)

  end subroutine adapt_grid
  !===========================================================================
  function exact_v(Xyz_D, Time)

    use BATL_lib,    ONLY: IsCartesianGrid, IsRzGeometry, IsPeriodic_D, &
         IsRotatedCartesian, GridRot_DD, DomainSize_D
    use ModNumConst, ONLY: cHalfPi
    use ModCoordTransform, ONLY: rot_matrix_z

    real:: exact_v(nVar)

    real, intent(in):: Xyz_D(nDim)
    real, intent(in):: Time

    ! Return the exact solution at Xyz_D and Time

    ! Square of the radius of the circle/sphere
    real:: XyzShift_D(nDim)

    real :: r1, r2, Rho, Rot_DD(3,3)
    !-------------------------------------------------------------------------
    ! Move position back to initial point
    if(RadialVelocity > 0.0) then
       r1 = sqrt(sum(Xyz_D**2))
       r2 = max(0.0, r1 - Time*RadialVelocity)
       XyzShift_D = (r2/r1)*Xyz_D
    elseif(AngularVelocity /= 0.0)then
       Rot_DD = rot_matrix_z(Time*AngularVelocity)
       XyzShift_D = matmul(Xyz_D, Rot_DD(1:nDim,1:nDim))
    else
       XyzShift_D = Xyz_D - Time*Velocity_D
    end if

    ! Take periodicity into account for Cartesian and RZ geometries only
    if(IsCartesianGrid .or. IsRotatedCartesian) then

       ! Rotate into generalized coordinates to get periodicity right
       if(IsRotatedCartesian) &
            XyzShift_D = matmul(XyzShift_D, GridRot_DD(1:nDim,1:nDim))

       where(IsPeriodic_D(1:nDim))
          XyzShift_D = &
               modulo(XyzShift_D - DomainMin_D(1:nDim), DomainSize_D(1:nDim)) &
               + DomainMin_D(1:nDim)
       end where

       ! Rotate back to XYZ coordinates
       if(IsRotatedCartesian) &
            XyzShift_D = matmul(GridRot_DD(1:nDim,1:nDim), XyzShift_D)
    end if
    r2 = sum((XyzShift_D - BlobCenter_D)**2)

    Rho = 1.0
    if(r2 < BlobRadius2) Rho = Rho + cos(cHalfPi*sqrt(r2/BlobRadius2))**2

    if(RadialVelocity > 0)then
       ! Scale by 1/r1**(D-1) so radial flow is stationary solution
       exact_v(Rho_) = Rho/r1**(nDim-1)
    elseif(IsRzGeometry) then
       exact_v(Rho_) = Rho/Xyz_D(r_)
    else
       exact_v(Rho_) = Rho
    end if

    exact_v(Lin_) = sum(Linear_D*(Xyz_D - Velocity_D*Time))

  end function exact_v

  !===========================================================================
  subroutine initialize

    use BATL_lib, ONLY: init_mpi, init_batl, init_grid_batl, &
         iProc, MaxDim, MaxBlock, nBlock, Unused_B, Xyz_DGB, &
         iTree_IA, MaxLevel_, BetaProlong, IsNodeBasedGrid, IsCartesian, &
         rRound0, rRound1

    use BATL_amr, ONLY: UseSimpleRefinement
    use BATL_geometry, ONLY: GridRot_DD, IsRotatedCartesian

    use ModReadParam, ONLY: read_file, read_init, &
         read_line, read_command, read_var

    use ModNumConst, ONLY: cDegToRad

    character(len=100):: StringLine, NameCommand

    integer:: nRoot_D(MaxDim) = (/4,4,2/)
    logical, allocatable:: DoRefine_B(:)
    real :: BlobRadius, Rgen_I(3)
    integer :: iDim, i, j, k, iBlock, iLevel
    logical:: IsNodeBasedRead = .true.
    logical:: UseUniformAxis
    !------------------------------------------------------------------------

    call init_mpi

    call read_file('PARAM.in')
    call read_init
    READPARAM: do
       if(.not.read_line(StringLine) ) EXIT READPARAM
       if(.not.read_command(NameCommand)) CYCLE READPARAM
       select case(NameCommand)
       case("#GRID")
          do iDim = 1, nDim
             call read_var('nRoot_D', nRoot_D(iDim))
          end do
          do iDim = 1, nDim
             call read_var('DomainMin_D', DomainMin_D(iDim))
             call read_var('DomainMax_D', DomainMax_D(iDim))
          end do
       case("#GRIDGEOMETRY")
          call read_var('TypeGeometry', TypeGeometry)
          if(TypeGeometry == 'roundcube') then
             call read_var('rRound0',rRound0)
             call read_var('rRound1',rRound1)
          end if
       case("#UNIFORMAXIS")
          call read_var('UseUniformaxis', UseUniformAxis)
       case("#GRIDTYPE")
          call read_var('IsNodeBasedGrid', IsNodeBasedRead)
       case("#PERIODIC")
          call read_var('IsPeriodic1', IsPeriodicGrid_D(1))
          if(nDim > 1) call read_var('IsPeriodic2', IsPeriodicGrid_D(2))
          if(nDim > 2) call read_var('IsPeriodic3', IsPeriodicGrid_D(3))
       case("#ADAPT")
          call read_var('DoAdapt', DoAdapt)
       case("#AMR")
          call read_var('MaxLevel', MaxLevel)
       case("#AMRINIT")
          call read_var('RefineRadius', RefineRadius2)
          RefineRadius2 = RefineRadius2**2
       case("#AMRCRITERIA")
          call read_var('RhoCoarsen', RhoCoarsen)
          call read_var('RhoRefine',  RhoRefine)
       case("#BLOB")
          call read_var('BlobRadius', BlobRadius)
          BlobRadius2 = BlobRadius**2
          do iDim = 1, nDim
             call read_var('BlobCenter_D', BlobCenter_D(iDim))
          end do
       case("#VELOCITY")
          do iDim = 1, nDim
             call read_var('Velocity_D', Velocity_D(iDim))
          end do
       case("#RADIALFLOW")
          call read_var('RadialVelocity', RadialVelocity)
       case("#ROTATINGFLOW")
          call read_var('AngularVelocity', AngularVelocity)
          AngularVelocity = AngularVelocity*cDegToRad
       case("#SCHEME")
          call read_var('nOrder', nOrder)
          if(nOrder > 1)then
             call read_var('BetaLimiter', BetaLimiter)
             call read_var('DoLimitVolumeState', DoLimitVolumeState)
          end if
       case("#TIMESTEP")
          call read_var('Cfl', Cfl)
       case("#LOCALSTEP")
          call read_var('TimeStartLocalStep', TimeStartLocalStep)
       case("#SAVEPLOT")
          call read_var('DtPlot', DtPlot)
       case("#STOP")
          call read_var('TimeMax', TimeMax)
       case("#TEST", "#TESTIJK", "#TESTXYZ", "#TEST2IJK", "#TEST2XYZ", &
            "#VERBOSE")
          call read_test_param(NameCommand)
       case default
          call CON_stop(NameSub//' unknown command='//trim(NameCommand))
       end select
    end do READPARAM

    ! Set initial refinement radius if not set
    if(RefineRadius2 < 0.0) RefineRadius2 = BlobRadius2

    UseConstantVelocity = RadialVelocity == 0.0 .and. AngularVelocity == 0.0

    ! Setup generalized radius array
    if(index(TypeGeometry,'genr') > 0) Rgen_I = &
         (/ DomainMin_D(1), (DomainMin_D(1)**3*DomainMax_D(1))**0.25, &
         DomainMax_D(1) /)

    ! Note that the periodicity will be fixed based on TypeGeometry
    call init_batl( &
         MaxBlockIn     = 8000, &          
         CoordMinIn_D   = DomainMin_D,  &
         CoordMaxIn_D   = DomainMax_D,  &
         nRootIn_D      = nRoot_D,      & 
         TypeGeometryIn = TypeGeometry, &
         IsPeriodicIn_D = IsPeriodicGrid_D, &
         RgenIn_I       = Rgen_I, &
         UseUniformAxisIn = UseUniformAxis)

    if(IsRotatedCartesian) then
       ! Rotate velocity vector the same way as the grid
       Velocity_D = matmul(GridRot_DD(1:nDim,1:nDim), Velocity_D)

       ! Rotate the blob center the same way as the grid
       BlobCenter_D = matmul(GridRot_DD(1:nDim,1:nDim), BlobCenter_D)

       ! Rotate the test cell location
       XyzTest_D = matmul(GridRot_DD, XyzTest_D)
    end if

    ! The default value is .not.IsRzGeometry
    UseSimpleRefinement = IsCartesian .or. IsRotatedCartesian

    IsNodeBasedGrid = IsNodeBasedRead

    ! Allow only MaxLevel levels of refinement
    iTree_IA(MaxLevel_,:) = MaxLevel

    allocate(DoRefine_B(MaxBlock))
    do iLevel = 1, MaxLevel
       DoRefine_B = .false.
       LOOPBLOCK: do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(sum((Xyz_DGB(1:nDim,i,j,k,iBlock) - BlobCenter_D)**2) &
                  < RefineRadius2)then
                DoRefine_B(iBlock) = .true.
                CYCLE LOOPBLOCK
             end if
          end do; end do; end do
       end do LOOPBLOCK

       call init_grid_batl(DoRefine_B)
       UseSimpleRefinement = IsCartesian .or. IsRotatedCartesian

    end do
    deallocate(DoRefine_B)

    ! Initial time step and time
    iStep    = 0
    Time     = 0.0
    TimePlot = 0.0

    allocate( &
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), &
         StateOld_VCB(nVar,nI,nJ,nK,MaxBlock), &
         Flux_VXB(nVar,nJ,nK,2,MaxBlock), &
         Flux_VYB(nVar,nI,nK,2,MaxBlock), &
         Flux_VZB(nVar,nI,nJ,2,MaxBlock) )

    State_VGB = 0.0
    Flux_VFD = 0.0
    Flux_VXB = 0.0
    Flux_VYB = 0.0
    Flux_VZB = 0.0

    ! Initialize the state as a sphere with a cos^2 density profile
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          State_VGB(:,i,j,k,iBlock) &
               = exact_v(Xyz_DGB(1:nDim,i,j,k,iBlock), Time)
       end do; end do; end do
    end do

    ! Set BetaProlong to 1.5
    BetaProlong = 1.5

    if(iProc==0)then
       call timing_active(.true.)
       call timing_step(0)
       call timing_depth(-1)
       call timing_report_style('tree')
    end if

  end subroutine initialize

  !===========================================================================
  subroutine save_log

    ! Calculate the totals on processor 0
    use BATL_lib, ONLY: nBlock, Unused_B, CellVolume_GB, Xyz_DGB, &
         iComm, iProc, nProc
    use ModMpi,    ONLY: MPI_reduce, MPI_REAL, MPI_SUM
    use ModIoUnit, ONLY: UnitTmp_

    integer:: iBlock, i, j, k, iError
    integer, parameter:: Volume_=nVar+1, Error_=nVar+1, nTotal=2*nVar+1
    real:: TotalPe_I(nTotal), Total_I(nTotal)

    character(len=100):: NameFile = 'advect.log'
    logical :: DoInitialize = .true.
    !------------------------------------------------------------------------
    Total_I = 0.0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          Total_I(1:nVar) = Total_I(1:nVar) + &
               CellVolume_GB(i,j,k,iBlock)*State_VGB(:,i,j,k,iBlock)

          Total_I(Volume_) = Total_I(Volume_) + CellVolume_GB(i,j,k,iBlock)

          Total_I(Error_+1:Error_+nVar) = Total_I(Error_+1:Error_+nVar) &
               + CellVolume_GB(i,j,k,iBlock) &
               *abs(State_VGB(:,i,j,k,iBlock)  &
               -    exact_v(Xyz_DGB(1:nDim,i,j,k,iBlock), Time))
       end do; end do; end do
    end do

    if(nProc > 1)then
       TotalPe_I = Total_I
       call MPI_reduce(TotalPe_I, Total_I, nTotal, MPI_REAL, MPI_SUM, 0, &
            iComm, iError)
    end if

    if(iProc /= 0) RETURN

    ! Divide total error by total mass to get relative error
    Total_I(Error_+1) = Total_I(Error_+1) / Total_I(1)

    if(DoInitialize)then
       open(UnitTmp_,file=NameFile, status='replace')
       write(UnitTmp_,'(a)')'Advection test for BATL'
       write(UnitTmp_,'(a)')'step time mass linear volume error_rho error_lin'

       DoInitialize = .false.
    else
       open(UnitTmp_,file=NameFile,position='append')
    end if
    write(UnitTmp_,'(i8, 100es15.6)') iStep, Time, Total_I
    close(UnitTmp_)

  end subroutine save_log

  !===========================================================================

  subroutine save_plot

    use BATL_lib, ONLY: MaxDim, nBlock, Unused_B, &
         iComm, nProc, iProc, iNode_B, &
         TypeGeometry, IsCylindrical, IsSpherical, IsRLonLat, IsGenRadius, &
         Phi_, nDimAmr, CoordMin_D, CoordMax_D, nRgen, LogRgen_I, &
         CellVolume_GB, CellSize_DB, Xyz_DGB, CoordMin_DB, CoordMax_DB, &
         rRound0, rRound1, SqrtNDim, nDim

    use ModMpi,    ONLY: MPI_REAL, MPI_INTEGER, MPI_MIN, MPI_SUM, MPI_reduce
    use ModIoUnit, ONLY: UnitTmp_
    use ModKind,   ONLY: nByteReal
    use ModNumConst, ONLY: cPi, cTwoPi, cHalfPi

    character(len=100):: NameSnapshot, NameFile
    real:: CellSizeMin_D(MaxDim), CellSizeMinAll_D(MaxDim), &
         CellSizePlot_D(MaxDim)
    real:: PlotMax_D(MaxDim), PlotMin_D(MaxDim)
    integer :: iDim, iBlock, i, j, k, iError
    integer :: nCell, nCellAll, nPlotDim, iPlot, nPlot, iXyz
    character (len=1) :: c0
    !-----------------------------------------------------------------------

    ! Calculate minimum cell size
    do iDim = 1, MaxDim
       CellSizeMin_D(iDim) = &
            minval(CellSize_DB(iDim,1:nBlock), MASK=.not.Unused_B(1:nBlock))
    end do
    call MPI_reduce(CellSizeMin_D, CellSizeMinAll_D, MaxDim, MPI_REAL, &
         MPI_MIN, 0, iComm, iError)

    ! Swap theta and phi to conform the r-lon-lat grid used in PostIDL
    if(IsSpherical) &
         CellSizeMinAll_D = CellSizeMinAll_D( (/1,3,2/) )

    ! Set plot sizes
    CellSizePlot_D = CellSizeMinAll_D

    ! uncomment "if" to create sturctured grids if there is no AMR
    ! if(MaxLevel > 0)then 
    ! Indicates to PostIDL that there is AMR in first element
    CellSizePlot_D(1) = -1.0
    ! Indicate full AMR by setting all values to -1
    if(nDimAmr == nDim) CellSizePlot_D = -1.0
    ! end if

    nPlot = 1
    if(nDim == 3) nPlot = 2

    do iPlot = 1, nPlot
       if(nDim==1)then
          ! 1D cut along X axis
          write(NameSnapshot,'(a,i7.7)') 'plots/cut_var_1_n',iStep
       elseif(iPlot == 1)then
          ! 2D plot in Z=0 plane
          write(NameSnapshot,'(a,i7.7)') 'plots/z=0_var_1_n',iStep
          iXyz = 3
       else
          ! 2D plot in Y=0 plane
          write(NameSnapshot,'(a,i7.7)') 'plots/y=0_var_2_n',iStep
          iXyz = 2
       end if

       ! write data from all processors into separate files
       write(NameFile,'(a,i4.4,a)') trim(NameSnapshot)//'_pe',iProc,'.idl'
       open(UnitTmp_, file=NameFile, status='replace', form='unformatted')
       nCell = 0
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          if(IsSpherical .and. iPlot == 1)then
             if(CoordMin_DB(2,iBlock) > cHalfPi) CYCLE
             if(CoordMax_DB(2,iBlock) < cHalfPi) CYCLE
          elseif((IsSpherical .or. IsRLonLat) .and. iPlot == 2)then
             if(  CoordMin_DB(Phi_,iBlock) > 0  .and.  &
                  CoordMax_DB(Phi_,iBlock) < cPi       ) CYCLE
             if(  CoordMin_DB(Phi_,iBlock) > cPi .and. &
                  CoordMax_DB(Phi_,iBlock) < cTwoPi    ) CYCLE
          elseif(IsCylindrical .and. iPlot == 2)then
             if(  CoordMin_DB(Phi_,iBlock) > 1e-6          .and. &
                  CoordMax_DB(Phi_,iBlock) < cPi    - 1e-6       ) CYCLE   
             if(  CoordMin_DB(Phi_,iBlock) > cPi    + 1e-6 .and. &
                  CoordMax_DB(Phi_,iBlock) < cTwoPi - 1e-6       ) CYCLE
          elseif(nDim > 1)then
             if(CoordMin_DB(iXyz,iBlock) > 1e-6) CYCLE
             if(CoordMax_DB(iXyz,iBlock) <-1e-6) CYCLE
          end if
          do k = 1, nK 
             if(nDim > 2 .and. (iPlot == 1 .and. .not. IsSpherical &
                  .or.          iPlot == 2 .and.       IsSpherical) )then
                ! Check for sign change of Y or Z along 3rd coordinate
                if(product(Xyz_DGB(iXyz,1,1,k-1:k+1:2,iBlock)) > 0) CYCLE
             end if
             do j = 1, nJ; 
                if(  iPlot == 1 .and.       IsSpherical .or. &
                     iPlot == 2 .and. .not. IsSPherical)then
                   ! Check for sign change of Y or Z along second coordinate
                   if(product(Xyz_DGB(iXyz,1,j-1:j+1:2,1,iBlock)) > 0) CYCLE
                end if
                do i = 1, nI
                   nCell = nCell + 1
                   write(UnitTmp_) CellSize_DB(1,iBlock), &
                        Xyz_DGB(:,i,j,k,iBlock), &
                        CoordMin_DB(:,iBlock) &
                        + ((/i,j,k/)-0.5)*CellSize_DB(:,iBlock), &
                        State_VGB(:,i,j,k,iBlock), &
                        exact_v(Xyz_DGB(:,i,j,k,iBlock), Time), &
                        CellVolume_GB(i,j,k,iBlock), real(iNode_B(iBlock)), &
                        real(iProc), real(iBlock)
                end do
             end do
          end do
       end do
       close(UnitTmp_)

       call MPI_reduce(nCell, nCellAll, 1, MPI_INTEGER, MPI_SUM, 0, &
            iComm, iError)

       ! Write header file
       if(iProc == 0)then

          nPlotDim = min(2,nDim)

          ! set plotting range

          if(IsSpherical)then
             ! Conform with the r-lon-lat grid used in PostIDL
             PlotMin_D(1) = CoordMin_D(1)
             PlotMax_D(1) = CoordMax_D(1)
             PlotMin_D(2) = CoordMin_D(3)
             PlotMax_D(2) = CoordMax_D(3)
             PlotMin_D(3) = cHalfPi - CoordMax_D(2)
             PlotMax_D(3) = cHalfPi - CoordMin_D(2)
          else
             PlotMin_D = CoordMin_D
             PlotMax_D = CoordMax_D
          end if

          if(nDim == 1)then
             PlotMin_D(2:3) = -1e-10
             PlotMax_D(2:3) = +1e-10
          elseif(iPlot == 1)then
             PlotMin_D(3) = -1e-10
             PlotMax_D(3) = +1e-10
          else
             PlotMin_D(2) = -1e-10 
             PlotMax_D(2) = +1e-10
          end if

          NameFile = trim(NameSnapshot)//'.h'
          open(UnitTmp_,file=NameFile,status="replace")
          write(UnitTmp_,'(a)') '#HEADFILE'
          write(UnitTmp_,'(a)') NameFile
          write(UnitTmp_,'(i8,a18)') nProc, 'nProc'        
          write(UnitTmp_,'(l8,a18)') .true.,' DoSaveBinary'
          write(UnitTmp_,'(i8,a18)')nByteReal,' nByteReal'
          write(UnitTmp_,*)

          write(UnitTmp_,'(a)') '#NDIM'
          write(UnitTmp_,'(i8,a18)') nDim, 'nDim'
          write(UnitTmp_,*)
           
          write(UnitTmp_,'(a)') '#NSTEP'
          write(UnitTmp_,'(i8,a18)')iStep, 'nStep'
          write(UnitTmp_,*)

          write(UnitTmp_,'(a)') '#TIMESIMULATION'
          write(UnitTmp_,'(1pe18.10,a18)')Time, 'TimeSimulation'
          write(UnitTmp_,*)        

          write(UnitTmp_,'(a)') '#PLOTRANGE'
          do iDim = 1, nDim
             write(c0,'(i1)') iDim
             write(UnitTmp_,'(1pe18.10,a18)') &
                  PlotMin_D(iDim), 'Coord'//c0//'Min'
             write(UnitTmp_,'(1pe18.10,a18)') &
                  PlotMax_D(iDim), 'Coord'//c0//'Max'
          enddo
          write(UnitTmp_,*)

          write(UnitTmp_,'(a)') '#PLOTRESOLUTION'
          do iDim = 1, nDim
             write(c0,'(i1)') iDim
             write(UnitTmp_,'(1pe18.10,a18)') &
                  CellSizePlot_D(iDim), 'DxSavePlot'//c0
          enddo
          write(UnitTmp_,*)

          write(UnitTmp_,'(a)') '#CELLSIZE'
          do iDim = 1, nDim
             write(c0,'(i1)') iDim
             write(UnitTmp_,'(1pe18.10,a18)') &
                  CellSizeMinAll_D(iDim), 'CellSizeMin'//c0
          enddo
          write(UnitTmp_,*)

          write(UnitTmp_,'(a)') '#NCELL'
          write(UnitTmp_,'(i10,a18)') nCellAll, 'nCellPlot'
          write(UnitTmp_,*)

          write(UnitTmp_,'(a)') '#PLOTVARIABLE'
          write(UnitTmp_,'(i8,a18)') 2*nVar+7, 'nPlotVar'
          write(UnitTmp_,'(a)') &
               'coord1 coord2 coord3 rho lin rhoexact linexact volume '// &
               'node proc block none' 
          write(UnitTmp_,'(a)')   '1 1 1'        ! units
          write(UnitTmp_,*)

          write(UnitTmp_,'(a)') '#GRIDGEOMETRYLIMIT'
          if(IsRLonLat)then
             write(UnitTmp_,'(a,a18)') 'spherical'//TypeGeometry(8:20), &
                  'TypeGeometry'
          else
             write(UnitTmp_,'(a,a18)')TypeGeometry, 'TypeGeometry'
          end if       
          if(TypeGeometry == 'roundcube')then
             write(UnitTmp_,'(es13.5," rRound0")') rRound0
             write(UnitTmp_,'(es13.5," rRound1")') rRound1
             write(UnitTmp_,'(es13.5," SqrtNDim")') SqrtNDim
          end if
          if(IsGenRadius)then
             write(UnitTmp_,'(i8,a)')  nRgen, ' nRgen'
             write(UnitTmp_,'(es13.5," LogRgen")') LogRgen_I
          end if
          write(UnitTmp_,*)

          write(UnitTmp_,'(a)') '#OUTPUTFORMAT'
          write(UnitTmp_,'(a)') 'real4'
          write(UnitTmp_,*)

          close(UnitTmp_)
       end if
    end do

  end subroutine save_plot
  !===========================================================================
  subroutine save_plot_block

    use ModPlotFile, ONLY: save_plot_file
    use BATL_lib, ONLY: iProc, &
         nBlock, Unused_B, CoordMin_DB, CoordMax_DB, CellSize_DB
    
    integer :: iBlock
    character(len=100):: NameFile
    character (len=10) :: TypePosition = 'rewind'
    !---------------------------------------------------------------------

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       write(NameFile,'(a,i3.3,a,i5.5,a)') &
            'advect_pe',iProc,'_blk',iBlock,'.out'
       
       call save_plot_file(NameFile,     &
            TypeFileIn='real4',          &
            TypePositionIn=TypePosition, &
            nStepIn = iStep, &
            TimeIn  = Time, &
            nDimIn  = nDim, &
            CoordMinIn_D = CoordMin_DB(1:nDim,iBlock)        &
            +              0.5*CellSize_DB(1:nDim,iBlock),   &
            CoordMaxIn_D = CoordMax_DB(1:nDim,iBlock)        &
            -              0.5*CellSize_DB(1:nDim,iBlock),   &
            VarIn_VIII = State_VGB(:,1:nI,1:nJ,1:nK,iBlock))
    end do

    TypePosition = 'append'

  end subroutine save_plot_block
  !===========================================================================
  subroutine finalize

    use BATL_lib, ONLY: clean_batl, clean_mpi
    !------------------------------------------------------------------------
    call save_plot
    call save_log

    call clean_batl
    call clean_mpi

  end subroutine finalize
  !===========================================================================
  subroutine set_boundary(TimeBc)

    use BATL_lib, ONLY: nBlock, Unused_B
    real, intent(in):: TimeBc

    ! Set ghost cells outside computational domain for time = TimeBc
    integer:: iBlock
    !------------------------------------------------------------------------
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       call set_boundary_block(iBlock, TimeBc)
    end do
  end subroutine set_boundary
  !===========================================================================
  subroutine set_boundary_block(iBlock, TimeBc)

    use BATL_lib, ONLY: nDim, IsPeriodic_D, Xyz_DGB, iNode_B, get_tree_position

    integer, intent(in):: iBlock
    real, intent(in):: TimeBc

    real:: PositionMin_D(MaxDim), PositionMax_D(MaxDim), State_V(nVar)

    integer:: i, j, k
    !--------------------------------------------------------------------------
    call get_tree_position(iNode_B(iBlock), PositionMin_D, PositionMax_D)

    if(PositionMin_D(1) < 1e-10)then
       ! Min in dimension 1
       if(IsPeriodic_D(1))then
          ! Update Lin including first physical cells (!)
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, 2
             State_V = exact_v(Xyz_DGB(1:nDim,i,j,k,iBlock), TimeBc)
             State_VGB(Lin_,i,j,k,iBlock) = State_V(Lin_)
          end do; end do; end do
       else
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, 0
             State_V = exact_v(Xyz_DGB(1:nDim,i,j,k,iBlock), TimeBc)
             State_VGB(:,i,j,k,iBlock) = State_V
          end do; end do; end do
       end if
    end if

    if(PositionMax_D(1) > 1 - 1e-10)then
       ! Max in dimension 1
       if(IsPeriodic_D(1))then
          ! Update Lin including first physical cell (!)
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = nI-1, MaxI
             State_V = exact_v(Xyz_DGB(1:nDim,i,j,k,iBlock), TimeBc)
             State_VGB(Lin_,i,j,k,iBlock) = State_V(Lin_)
          end do; end do; end do
       else
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = nI+1, MaxI
             State_V = exact_v(Xyz_DGB(1:nDim,i,j,k,iBlock), TimeBc)
             State_VGB(:,i,j,k,iBlock) = State_V
          end do; end do; end do
       end if
    end if

    if(nDim > 1)then

       if(PositionMin_D(2) < 1e-10)then
          ! Min in dimension 2
          if(IsPeriodic_D(2))then
             ! Update Lin including first physical cells (!)
             do k = MinK, MaxK; do j = MinJ, 2; do i = MinI, MaxI
                State_V = exact_v(Xyz_DGB(1:nDim,i,j,k,iBlock), TimeBc)
                State_VGB(Lin_,i,j,k,iBlock) = State_V(Lin_)
             end do; end do; end do
          else
             do k = MinK, MaxK; do j = MinJ, 0; do i = MinI, MaxI
                State_V = exact_v(Xyz_DGB(1:nDim,i,j,k,iBlock), TimeBc)
                State_VGB(:,i,j,k,iBlock) = State_V
             end do; end do; end do
          end if
       end if

       if(PositionMax_D(2) > 1 - 1e-10)then
          ! Max in dimension 2
          if(IsPeriodic_D(2))then
             ! Update Lin including first physical cells (!)
             do k = MinK, MaxK; do j = nJ-1, MaxJ; do i = MinI, MaxI
                State_V = exact_v(Xyz_DGB(1:nDim,i,j,k,iBlock), TimeBc)
                State_VGB(Lin_,i,j,k,iBlock) = State_V(Lin_)
             end do; end do; end do
          else
             do k = MinK, MaxK; do j = nJ+1, MaxJ; do i = MinI, MaxI
                State_V = exact_v(Xyz_DGB(1:nDim,i,j,k,iBlock), TimeBc)
                State_VGB(:,i,j,k,iBlock) = State_V
             end do; end do; end do
          end if
       end if
    end if

    if(nDim > 2)then
       if(PositionMin_D(3) < 1e-10)then
          ! Min in dimension 3
          if(IsPeriodic_D(3))then
             ! Update Lin including first physical cells (!)
             do k = MinK, 2; do j = MinJ, MaxJ; do i = MinI, MaxI
                State_V = exact_v(Xyz_DGB(1:nDim,i,j,k,iBlock), TimeBc)
                State_VGB(Lin_,i,j,k,iBlock) = State_V(Lin_)
             end do; end do; end do
          else
             do k = MinK, 0; do j = MinJ, MaxJ; do i = MinI, MaxI
                State_V = exact_v(Xyz_DGB(1:nDim,i,j,k,iBlock), TimeBc)
                State_VGB(:,i,j,k,iBlock) = State_V
             end do; end do; end do
          end if
       end if

       if(PositionMax_D(3) == 1 - 1e-10)then
          ! Max in dimension 3
          if(IsPeriodic_D(3))then
             ! Update Lin including first physical cells (!)
             do k = nK-1, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
                State_V = exact_v(Xyz_DGB(1:nDim,i,j,k,iBlock), TimeBc)
                State_VGB(Lin_,i,j,k,iBlock) = State_V(Lin_)
             end do; end do; end do
          else
             do k = nK+1, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
                State_V = exact_v(Xyz_DGB(1:nDim,i,j,k,iBlock), TimeBc)
                State_VGB(:,i,j,k,iBlock) = State_V
             end do; end do; end do
          end if
       end if
    end if

  end subroutine set_boundary_block
  !===========================================================================
  subroutine calc_face_flux(iBlock)

    use ModNumConst, ONLY: i_DD
    use BATL_lib, ONLY: CellVolume_GB

    integer, intent(in):: iBlock

    real:: StateLeft_VFD( nVar,1:nI+1,1:nJ+1,1:nK+1,nDim)
    real:: StateRight_VFD(nVar,1:nI+1,1:nJ+1,1:nK+1,nDim)
    real:: Slope_VGD(nVar,0:nI+1,0:nJ+1,0:nK+1,nDim)
    real:: State_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    real:: vInv, vFaceNormal

    integer :: iDim, i, j, k, Di, Dj, Dk

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_face_flux'
    !------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(DoTest)then
       StateLeft_VFD  = -777.7
       StateRight_VFD = -777.7
       vFaceNormal_DF = -777.7
       Flux_VFD       = -777.7
    end if

    if(nOrder==1)then
       do iDim = 1, nDim
          Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
          do k = 1, nK+Dk; do j =1, nJ+Dj; do i = 1, nI+Di
             StateLeft_VFD( :,i,j,k,iDim) = State_VGB(:,i-Di,j-Dj,k-Dk,iBlock)
             StateRight_VFD(:,i,j,k,iDim) = State_VGB(:,i,j,k,iBlock)
          end do; end do; end do
       end do
    else
       State_VG = State_VGB(:,:,:,:,iBlock)
       if(DoLimitVolumeState)then
          do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
             State_VG(:,i,j,k) = State_VG(:,i,j,k)*CellVolume_GB(i,j,k,iBlock)
          end do; end do; end do
       end if
       call limit_slope(iBlock, State_VG, Slope_VGD)

       do iDim = 1, nDim
          Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
          do k = 1, nK+Dk; do j =1, nJ+Dj; do i = 1, nI+Di
             StateLeft_VFD( :,i,j,k,iDim) = State_VG(:,i-Di,j-Dj,k-Dk) &
                  + 0.5*Slope_VGD(:,i-Di,j-Dj,k-Dk,iDim)
             StateRight_VFD(:,i,j,k,iDim) = State_VG(:,i,j,k) &
                  - 0.5*Slope_VGD(:,i,j,k,iDim)
          end do; end do; end do

          if(DoLimitVolumeState)then
             ! Divide by volume averaged over the face
             ! Maybe some more accurate approximation is needed???
             do k = 1, nK+Dk; do j =1, nJ+Dj; do i = 1, nI+Di
                vInv = 2/(abs(CellVolume_GB(i-Di,j-Dj,k-Dk,iBlock)) &
                     +    abs(CellVolume_GB(i,j,k,iBlock)))

                StateLeft_VFD( :,i,j,k,iDim)=vInv*StateLeft_VFD( :,i,j,k,iDim)
                StateRight_VFD(:,i,j,k,iDim)=vInv*StateRight_VFD(:,i,j,k,iDim)
             end do; end do; end do
          end if

       end do
    end if

    ! Calculate fluxes
    call calc_vface_normal(iBlock)

    do iDim = 1, nDim
       Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
       do k = 1, nK+Dk; do j =1, nJ+Dj; do i = 1, nI+Di

          vFaceNormal = vFaceNormal_DF(iDim,i,j,k)
          ! Calculate upwinded flux
          if(vFaceNormal > 0)then
             Flux_VFD(:,i,j,k,iDim)=vFaceNormal*StateLeft_VFD(:,i,j,k,iDim)
          else
             Flux_VFD(:,i,j,k,iDim)=vFaceNormal*StateRight_VFD(:,i,j,k,iDim)
          end if
       end do; end do; end do
    end do

    if(DoTest)then
       write(*,*)NameSub, ': StateLeft_VFD =', &
            StateLeft_VFD(:,iTest,jTest,kTest,:)
       write(*,*)NameSub, ': StateRight_VFD=', &
            StateRight_VFD(:,iTest,jTest,kTest,:)
       write(*,*)NameSub,': vFaceNormal_D=',vFaceNormal_DF(:,iTest,jTest,kTest)
       write(*,*)NameSub, ': Flux_VFD = ', Flux_VFD(:,iTest,jTest,kTest,:)
    end if

    call test_stop(NameSub, DoTest, iBlock)

  end subroutine calc_face_flux
  !===========================================================================
  subroutine limit_slope(iBlock, State_VG, Slope_VGD)

    ! Calculate TVD limited slopes Slope_GD of State_VG

    use ModNumConst, ONLY: i_DD

    integer, intent(in) :: iBlock  ! Passed for testing only
    real,    intent(in) :: State_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real,    intent(out):: Slope_VGD(nVar,0:nI+1,0:nJ+1,0:nK+1,nDim)

    real    :: SlopeLeft_V(nVar), SlopeRight_V(nVar)
    integer :: iDim, i, j, k, Di, Dj, Dk

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'limit_slope'
    !----------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(DoTest)Slope_VGD = -777.7

    do iDim = 1, nDim
       Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
       do k = 1-Dk, nK+Dk; do j = 1-Dj, nJ+Dj; do i = 1-Di, nI+Di
          SlopeLeft_V  = State_VG(:,i,j,k)  - State_VG(:,i-Di,j-Dj,k-Dk)
          SlopeRight_V = State_VG(:,i+Di,j+Dj,k+Dk) - State_VG(:,i,j,k)
          Slope_VGD(:,i,j,k,iDim) = &
               (sign(0.5, SlopeLeft_V) + sign(0.5, SlopeRight_V))*min( &
               BetaLimiter*abs(SlopeLeft_V), &
               BetaLimiter*abs(SlopeRight_V), &
               0.5*abs(SlopeLeft_V + SlopeRight_V))
       end do; end do; end do
    end do

    if(DoTest)then
       write(*,*)            NameSub, ': State_VG(:,iTest-1:iTest+1)=', &
            State_VG(:,iTest-1:iTest+1,jTest,kTest)
       if(nDim>1) write(*,*) NameSub, ': State_VG(:,jTest-1:jTest+1)=', &
            State_VG(:,iTest,jTest-1:jTest+1,kTest)
       if(nDim>2) write(*,*) NameSub, ': State_VG(:,kTest-1:kTest+1)=', &
            State_VG(:,iTest,jTest,kTest-1:kTest+1)
       write(*,*) NameSub, ': Slope1=', Slope_VGD(:,iTest,jTest,kTest,:)
    end if

    call test_stop(NameSub, DoTest, iBlock)

  end subroutine limit_slope
  !===========================================================================
  subroutine calc_vface_normal(iBlock)

    use ModNumConst, ONLY: i_DD
    use BATL_lib, ONLY: IsCartesian, IsRzGeometry, Xyz_DGB, &
         CellFace_DB, CellFace_DFB, FaceNormal_DDFB, Dim1_, Dim2_

    integer, intent(in):: iBlock

    ! Calculate either Flux_DF fluxes or vFaceNormal_DF for CFL condition

    real:: XyzFace_D(nDim), vFaceNormal

    integer:: iDim, Di, Dj, Dk, i, j, k
    !-----------------------------------------------------------------------

    ! Calculate upwinded fluxes 
    do iDim = 1, nDim
       Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
       do k = 1, nK+Dk; do j =1, nJ+Dj; do i = 1, nI+Di

          ! Calculate velocity
          if(RadialVelocity > 0. .or. AngularVelocity > 0.)then
             XyzFace_D = 0.5*(Xyz_DGB(1:nDim,i,j,k,iBlock) &
                  +           Xyz_DGB(1:nDim,i+Di,j+Dj,k+Dk,iBlock))
             if(RadialVelocity > 0.)then
                Velocity_D = XyzFace_D/sqrt(sum(XyzFace_D**2))*RadialVelocity
             else
                Velocity_D(Dim1_) = -AngularVelocity*XyzFace_D(Dim1_)
                Velocity_D(Dim2_) =  AngularVelocity*XyzFace_D(Dim2_)
             end if
          end if

          ! Calculate normal velocity times face area
          if(IsCartesian)then
             vFaceNormal = CellFace_DB(iDim,iBlock)*Velocity_D(iDim)
          elseif(IsRzGeometry)then
             vFaceNormal = CellFace_DFB(iDim,i,j,k,iBlock)*Velocity_D(iDim)
          else
             vFaceNormal = sum(FaceNormal_DDFB(:,iDim,i,j,k,iBlock)*Velocity_D)
          end if

          vFaceNormal_DF(iDim,i,j,k) = vFaceNormal
       end do; end do; end do
    end do

  end subroutine calc_vface_normal

  !============================================================================

  subroutine advance_explicit

    use BATL_lib, ONLY: message_pass_cell, nBlock, Unused_B, &
         IsCartesianGrid, IsCartesian, &
         CellSize_DB, CellVolume_B, CellVolume_GB, &
         iComm, nProc, store_face_flux, apply_flux_correction
    use ModNumConst, ONLY: i_DD
    use ModMpi

    integer:: iStage, iDim, iBlock, i, j, k, Di, Dj, Dk, iError
    real:: DtInv, DtPe, DtStage, InvVolume, Flux

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'advance_explicit'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    DtInv = 0.0

    ! Calculate time step limit for each block
    if((IsCartesianGrid) .and. UseConstantVelocity)then
       ! Velocity is constant and positive so simply check cellsize per block
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          DtInv = max(DtInv, sum(abs(Velocity_D)/CellSize_DB(1:nDim,iBlock)) )
       end do
    else
       ! Cell width is estimated as Volume/FaceArea
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          ! Calculate 
          call calc_vface_normal(iBlock)

          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             Flux = 0.0
             do iDim = 1, nDim
                Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)

                Flux = Flux + max(abs(vFaceNormal_DF(iDim,i,j,k)), &
                     abs(vFaceNormal_DF(iDim,i+Di,j+Dj,k+Dk)))
             end do
             DtInv = max(DtInv, Flux/CellVolume_GB(i,j,k,iBlock))
          end do; end do; end do
       end do
    end if
    Dt = min( Cfl/DtInv, TimeMax - Time)

    if(nProc > 1)then
       DtPe = Dt
       call MPI_allreduce(DtPe, Dt, 1, MPI_REAL, MPI_MIN, iComm, iError)
    end if

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       StateOld_VCB(:,:,:,:,iBlock) = State_VGB(:,1:nI,1:nJ,1:nK,iBlock)
    end do

    do iStage = 1, nOrder

       if(DoTest)write(*,*)NameSub,': iProc, iStage, State=', &
            iProc, iStage, State_VGB(:,iTest,jTest,kTest,iBlockTest)

       DtStage = (Dt*iStage)/nOrder

       call set_boundary(Time + 0.5*Dt*(iStage-1))

       if(DoTest)write(*,*)NameSub,' after set_boundary State=', &
            State_VGB(:,iTest,jTest,kTest,iBlockTest)

       call timing_start('message_pass')
       call message_pass_cell(nVar, State_VGB)
       call timing_stop('message_pass')

       if(DoTest)write(*,*)NameSub,' after message_pass State=', &
            State_VGB(:,iTest,jTest,kTest,iBlockTest)

       ! This is only needed to fix Lin if periodic BCs were applied
       call set_boundary(Time + 0.5*Dt*(iStage-1))

       if(DoTest)write(*,*)NameSub,' after set_boundary State=', &
            State_VGB(:,iTest,jTest,kTest,iBlockTest)

       call timing_start('update')
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          call calc_face_flux(iBlock)

          call store_face_flux(iBlock, nVar, Flux_VFD, &
               Flux_VXB, Flux_VYB, Flux_VZB, &
               DtIn = DtStage, &
               DoStoreCoarseFluxIn = .true.)

          State_VGB(:,1:nI,1:nJ,1:nK,iBlock) = StateOld_VCB(:,:,:,:,iBlock)

          if(IsCartesian) InvVolume = 1/CellVolume_B(iBlock)
          do iDim = 1, nDim
             Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                if(.not.IsCartesian) InvVolume = 1/CellVolume_GB(i,j,k,iBlock)
                State_VGB(:,i,j,k,iBlock) = State_VGB(:,i,j,k,iBlock) &
                     - DtStage*InvVolume* &
                     (Flux_VFD(:,i+Di,j+Dj,k+Dk,iDim) - Flux_VFD(:,i,j,k,iDim))
             end do; end do; end do

          end do
       end do

       if(DoTest)write(*,*)NameSub,': before flux corr. State=', &
            State_VGB(:,iTest,jTest,kTest,iBlockTest)

       call apply_flux_correction(nVar, 0, State_VGB, &
            Flux_VXB=Flux_VXB, Flux_VYB=Flux_VYB, Flux_VZB=Flux_VZB)

       if(DoTest)write(*,*)NameSub,': after  flux corr. State=', &
            State_VGB(:,iTest,jTest,kTest,iBlockTest)

       call timing_stop('update')

    end do

    call test_stop(NameSub, DoTest)

  end subroutine advance_explicit
  !===========================================================================
  subroutine advance_localstep

    use BATL_lib, ONLY: message_pass_cell, MaxBlock, nBlock, Unused_B, &
         IsCartesian, CellSize_DB, CellVolume_B, CellVolume_GB, &
         iComm, nProc, store_face_flux, apply_flux_correction, &
         nLevelMin, nLevelMax, min_tree_level, &
         iTree_IA, Level_, iNode_B, CellFace_DB, DiLevelNei_IIIB
    use ModNumConst, ONLY: i_DD
    use ModMpi

    integer:: nStage, iStage, iStageBlock, iDim, iBlock
    integer:: i, j, k, Di, Dj, Dk, iError
    real :: DtMin, DtMax, DtMinPe, DtMaxPe
    real :: TimeStage, TimeBlock, DtLocal, InvVolume
    real, save, allocatable:: Dt_B(:), Time_B(:), TimeOld_B(:)
    integer, save, allocatable:: iStage_B(:)

    integer:: iLevelMin

    logical:: DoTest, DoTestBlock
    character(len=*), parameter:: NameSub = 'advance_localstep'
    !----------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(DoTest)then
       write(*,*) NameSub,' starting with nLevelMin, nLevelMax=', &
            nLevelMin, nLevelMax
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          if(is_incorrect_block(iBlock, Time)) write(*,*) '!!! just starting'
       end do
    end if

    if(.not.allocated(Dt_B))then
       allocate(Dt_B(MaxBlock), Time_B(MaxBlock), TimeOld_B(MaxBlock), &
            iStage_B(MaxBlock))
       Dt_B      = 0.0
    end if

    ! For newly created blocks this is needed
    iStage_B(1:nBlock)  = 1
    Time_B(1:nBlock)    = Time
    TimeOld_B(1:nBlock) = Time

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       Dt_B(iBlock) = Cfl / sum( Velocity_D/CellSize_DB(1:nDim,iBlock) )
       Time_B(iBlock) = Time
    end do
    DtMin = minval(Dt_B(1:nBlock), MASK = .not. Unused_B(1:nBlock))
    DtMax = maxval(Dt_B(1:nBlock), MASK = .not. Unused_B(1:nBlock))

    if(nProc > 1)then
       DtMinPe = DtMin
       call MPI_allreduce(DtMinPe, DtMin, 1, MPI_REAL, MPI_MIN, iComm, iError)
       DtMaxPe = DtMax
       call MPI_allreduce(DtMaxPe, DtMax, 1, MPI_REAL, MPI_MAX, iComm, iError)
    end if

    DtMax = floor( (DtMax + 1e-10)/DtMin )*DtMin
    Dt_B(1:nBlock) = floor( (Dt_B(1:nBlock) + 1e-10)/DtMin )*DtMin

    ! Exploit the fact that time step is proportional to cell size !!!
    Dt = DtMax
    nStage = nint(DtMax / DtMin)
    TimeStage = Time

    do iStage = 1, 2*nStage

       ! The blocks between iLevelMin
       iLevelMin = min_tree_level(iStage)
       if(DoTest)write(*,*) NameSub,' starting iStage, iLevelMin=', &
            iStage, iLevelMin

       ! Update the boundary conditions for these blocks
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          if(iTree_IA(Level_,iNode_B(iBlock)) < iLevelMin-1) CYCLE
          call set_boundary_block(iBlock, Time_B(iBlock))
       end do

       ! Messagepass for blocks between iLevelMin-1 and up.
       call timing_start('message_pass')
       call message_pass_cell(nVar, State_VGB, &
            TimeOld_B=TimeOld_B, Time_B=Time_B, iLevelMin=iLevelMin-1)
       call timing_stop('message_pass')

       ! Update the boundary conditions for Linear variable
       ! for periodic boundaries. Normally this is not really needed.
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          if(iTree_IA(Level_,iNode_B(iBlock)) < iLevelMin-1) CYCLE
          call set_boundary_block(iBlock, Time_B(iBlock))
       end do

       call timing_start('update')
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          TimeBlock = Time_B(iBlock)
          if(TimeBlock > TimeStage*(1+1e-10)) CYCLE

          DoTestBlock = DoTest .and. iBlock == iBlockTest

          iStageBlock = iStage_B(iBlock)

          if(DoTestBlock) &
               write(*,*) NameSub,' iStage, iStageBlock, TimeOld_B, Time_B=', &
               iStage, iStageBlock, TimeOld_B(iBlock), TimeBlock

          if(DoTestBlock)write(*,*) NameSub,' starting stage with state=', &
               State_VGB(:,iTest,jTest,kTest,iBlockTest), &
               exact_v(Xyz_DGB(1:nDim,iTest,jTest,kTest,iBlockTest), TimeBlock)

          if(DoTest)then
             if(is_incorrect_block(iBlock, Time_B(iBlock), UseGhost=.true.)) &
                  then
                write(*,*) '!!! at start of iStage=', iStage
                write(*,*) '!!! iLevelMin, iLevelB=', &
                     iLevelMin, iTree_IA(Level_,iNode_B(iBlock))
                write(*,*) '!!! DiLevelNei_IIIB(:,0,0)=', &
                     DiLevelNei_IIIB(:,0,0,iBlock)
                write(*,*) '!!! DiLevelNei_IIIB(0,:,0)=', &
                     DiLevelNei_IIIB(0,:,0,iBlock)
             end if
          end if

          call calc_face_flux(iBlock)

          DtLocal = Dt_B(iBlock) * iStageBlock/2.0

          if(iStageBlock==2) call store_face_flux(iBlock, nVar, Flux_VFD, &
               Flux_VXB, Flux_VYB, Flux_VZB, &
               DtIn = DtLocal, &
               DoStoreCoarseFluxIn = .true.)

          if(iStageBlock == 1)then
             StateOld_VCB(:,:,:,:,iBlock) = State_VGB(:,1:nI,1:nJ,1:nK,iBlock)
          else
             State_VGB(:,1:nI,1:nJ,1:nK,iBlock) = StateOld_VCB(:,:,:,:,iBlock)

             if(DoTestBlock)write(*,*) NameSub,' after setting to StateOld=', &
                  State_VGB(:,iTest,jTest,kTest,iBlockTest)
          end if

          if(IsCartesian) InvVolume = 1.0/CellVolume_B(iBlock)
          do iDim = 1, nDim
             Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                if(.not.IsCartesian) InvVolume = 1/CellVolume_GB(i,j,k,iBlock)
                State_VGB(:,i,j,k,iBlock) = State_VGB(:,i,j,k,iBlock) &
                     - DtLocal*InvVolume* &
                     (Flux_VFD(:,i+Di,j+Dj,k+Dk,iDim) - Flux_VFD(:,i,j,k,iDim))
             end do; end do; end do
          end do

          ! Update local time and stage for the block
          TimeOld_B(iBlock) = Time_B(iBlock)
          Time_B(iBlock)    = Time_B(iBlock) + Dt_B(iBlock)/2
          iStage_B(iBlock)  = 3 - iStageBlock

          if(DoTest)then
             if(is_incorrect_block(iBlock, Time_B(iBlock))) &
                  write(*,*) NameSub, 'after update in iStage=', iStage
          end if
          if(DoTestBlock)then
             ! Note old time
             write(*,*) NameSub,' DtLocal, InvVolume=', DtLocal, InvVolume
             write(*,*) NameSub,' flux x left, exact=', &
                  Flux_VFD(:,iTest,jTest,kTest,1), &
                  CellFace_DB(1,iBlock)*Velocity_D(1)* &
                  ( exact_v(Xyz_DGB(1:nDim,iTest-1,jTest,kTest,iBlockTest),&
                  TimeBlock) &
                  + exact_v(Xyz_DGB(1:nDim,iTest  ,jTest,kTest,iBlockTest),&
                  TimeBlock))/2
             write(*,*) NameSub,' flux x right, exact=', &
                  Flux_VFD(:,iTest+1,jTest,kTest,1), &
                  CellFace_DB(1,iBlock)*Velocity_D(1)* &
                  ( exact_v(Xyz_DGB(1:nDim,iTest  ,jTest,kTest,iBlockTest),&
                  TimeBlock) &
                  + exact_v(Xyz_DGB(1:nDim,iTest+1,jTest,kTest,iBlockTest),&
                  TimeBlock) )/2

             ! Note new time
             write(*,*) NameSub,' after update state, exact=', &
                  State_VGB(:,iTest,jTest,kTest,iBlockTest), &
                  exact_v(Xyz_DGB(1:nDim,iTest,jTest,kTest,iBlockTest),&
                  Time_B(iBlock))
          end if
       end do

       ! Apply flux correction at the end of full steps at level n-1 or below
       if(modulo(iStage, 4) == 0) then
          call apply_flux_correction(nVar, 0, State_VGB, &
               Flux_VXB=Flux_VXB, Flux_VYB=Flux_VYB, Flux_VZB=Flux_VZB, &
               iStageIn = iStage/2)

          if(DoTest)then
             write(*,*) NameSub,' after flux correction state, exact=', &
                  State_VGB(:,iTest,jTest,kTest,iBlockTest), &
                  exact_v(Xyz_DGB(1:nDim,iTest,jTest,kTest,iBlockTest), &
                  Time_B(iBlockTest))

             do iBlock = 1, nBlock
                if(Unused_B(iBlock)) CYCLE
                if(is_incorrect_block(iBlock, Time_B(iBlock))) &
                     write(*,*) NameSub,' after flux correction at iStage=',&
                     iStage
             end do
          end if
       end if

       TimeStage = TimeStage + DtMin/2

       call timing_stop('update')
    end do

    if(  maxval(Time_B, MASK=.not.Unused_B) > &
         minval(Time_B, MASK=.not.Unused_B) + 1e-10)then

       write(*,*)'ERROR: minval(Time_B)=',minval(Time_B, MASK=.not.Unused_B)
       write(*,*)'ERROR: maxval(Time_B)=',maxval(Time_B, MASK=.not.Unused_B)
       write(*,*)'ERROR: minloc(Time_B)=',minloc(Time_B, MASK=.not.Unused_B)
       write(*,*)'ERROR: maxloc(Time_B)=',maxloc(Time_B, MASK=.not.Unused_B)

       call CON_stop('time step problem')

    end if

    if(DoTest)write(*,*) NameSub,' State=', &
         State_VGB(:,iTest,jTest,kTest,iBlockTest), &
         exact_v(Xyz_DGB(1:nDim,iTest,jTest,kTest,iBlockTest), Time+Dt)

    call test_stop(NameSub, DoTest)

  end subroutine advance_localstep
  !===========================================================================
  logical function is_incorrect_block(iBlock, Time, UseGhost)

    ! Check block against the analytic solution
  
    integer, intent(in):: iBlock
    real,    intent(in):: Time
    logical, optional, intent(in):: UseGhost

    real:: State_V(nVar)
    integer:: i, j, k, iMin, iMax, jMin, jMax, kMin, kMax
    
    character(len=*), parameter:: NameSub = 'is_incorrect_block'
    !------------------------------------------------------------------------
    is_incorrect_block = .false.

    if(present(UseGhost))then
       iMin = MinI; iMax = MaxI 
       jMin = MinJ; jMax = MaxJ
       kMin = MinK; kMax = MaxK
    else
       iMin = 1; iMax = nI
       jMin = 1; jMax = nJ 
       kMin = 1; kMax = nK
    end if
       
    do k=kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
      State_V = exact_v(Xyz_DGB(:,i,j,k,iBlock), Time)
      if(abs(State_V(Lin_) - State_VGB(Lin_,i,j,k,iBlock)) > 1e-6)then

         is_incorrect_block = .true.

         write(*,*)'ERROR at iStep=', iStep
         write(*,*)'i, j, k, iBlock, iProc=', i, j, k, iBlock, iProc
         write(*,*)'Xyz, Time  =', Xyz_DGB(:,i,j,k,iBlock), Time
         write(*,*)'State_VGB  =', State_VGB(Lin_,i,j,k,iBlock)
         write(*,*)'Exact value=', State_V(Lin_)
      end if
    end do; end do; end do

  end function is_incorrect_block

end program advect

include 'external_routines.f90'
