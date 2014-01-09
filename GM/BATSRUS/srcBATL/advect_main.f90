!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program advect

  use BATL_lib, ONLY: nDim, nDimAmr, nI, nJ, nK, &
       MinI, MaxI, MinJ, MaxJ, MinK, MaxK, iProc, barrier_mpi, r_

  implicit none

  ! Geometry
  character(len=20):: TypeGeometry = 'cartesian'

  ! Square of the radius of the sphere
  real :: BlobRadius2 = 25.0
  real :: BlobCenter_D(nDim) = 0.0

  ! Maximum refinement level
  integer:: MaxLevel = 3

  ! Final simulation time, frequency of plots
  real :: TimeMax = 20.0
  real :: TimeStartLocalStep = 0.49999*20.0
  real :: DtPlot = 1.0

  ! Velocity field
  logical:: UseConstantVelocity = .false.
  real :: Velocity_D(nDim) = 0.0
  real :: RadialVelocity   = 0.0
  real :: AngularVelocity  = 0.0

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

  ! Time step counter and simulation time
  integer :: iStep
  real :: Time, TimePlot, Dt

  integer, parameter:: nVar = 1, Rho_ = 1

  ! Cell centered state
  real, allocatable :: State_VGB(:,:,:,:,:), StateOld_VCB(:,:,:,:,:)

  ! Face centered flux for one block
  real:: Flux_VFD(nVar,1:nI+1,1:nJ+1,1:nK+1,nDim)

  ! Normal velocity times area
  real:: vFaceNormal_DF(nDim,1:nI+1,1:nJ+1,1:nK+1)

  ! Face centered flux for conservation fix
  real, allocatable, dimension(:,:,:,:,:):: Flux_VXB, Flux_VYB, Flux_VZB

  logical, parameter :: DoTest = .false.
  character(len=*), parameter:: NameSub = 'advect_main'
  !--------------------------------------------------------------------------
  call initialize
  
  if(DoTest)write(*,*)NameSub,' starting iProc=',iProc
  if(DoTest)call barrier_mpi

  call timing_start('ADVECT')
  do
     if(DoTest)write(*,*)NameSub,' advance iProc, iStep=',iProc,iStep

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

     if(DoTest)write(*,*)NameSub,' adapt grid iProc=',iProc
     if(DoTest)call barrier_mpi

     call timing_start('amr')
     call adapt_grid
     call timing_stop('amr')

     call timing_step(iStep)

     if(Time >= TimeMax - 1e-10) EXIT

  end do
  call timing_stop('ADVECT')
  call timing_report_total

  if(DoTest)write(*,*)NameSub,' finalize iProc=',iProc

  call finalize

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
          if(State_VGB(Rho_,i,j,k,iBlock) > Factor*1.5) then
             iStatusNew_A(iNode_B(iBlock)) = Refine_
             CYCLE BLOCK
          end if
       end do; end do; end do
       ! If not refined, check for coarsening
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(IsRzGeometry) Factor = 1./Xyz_DGB(r_,i,j,k,iBlock)
          if(State_VGB(Rho_,i,j,k,iBlock) < Factor*1.4)then
             iStatusNew_A(iNode_B(iBlock)) = Coarsen_
             CYCLE BLOCK
          end if
       end do; end do; end do
    end do BLOCK

    call regrid_batl(nVar,State_VGB,DoBalanceEachLevelIn=UseLocalStep)

  end subroutine adapt_grid
  !===========================================================================
  real function exact_density(Xyz_D)

    use BATL_lib,    ONLY: IsCartesianGrid, IsRzGeometry, IsPeriodic_D
    use ModNumConst, ONLY: cHalfPi
    use ModCoordTransform, ONLY: rot_matrix_z

    real, intent(in):: Xyz_D(nDim)

    ! Square of the radius of the circle/sphere
    real:: XyzShift_D(nDim)

    real :: r1, r2, Rho, rot_DD(3,3)

    real:: DomainSize_D(nDim)
    !-------------------------------------------------------------------------
    ! Move position back to initial point
    if(RadialVelocity > 0.0) then
       r1 = sqrt(sum(Xyz_D**2))
       r2 = max(0.0, r1 - Time*RadialVelocity)
       XyzShift_D = (r2/r1)*Xyz_D
    elseif(AngularVelocity /= 0.0)then
       rot_DD = rot_matrix_z(Time*AngularVelocity)
       XyzShift_D = matmul(Xyz_D, rot_DD(1:nDim,1:nDim))
    else
       XyzShift_D = Xyz_D - Time*Velocity_D
    end if

    ! Take periodicity into account for Cartesian and RZ geometries only
    if(IsCartesianGrid) then
       DomainSize_D = DomainMax_D(1:nDim) - DomainMin_D(1:nDim)
       where(IsPeriodic_D(1:nDim))
          XyzShift_D = modulo(XyzShift_D - DomainMin_D(1:nDim), DomainSize_D) &
               + DomainMin_D(1:nDim)
       end where
    end if
    r2 = sum((XyzShift_D - BlobCenter_D)**2)

    Rho = 1.0
    if(r2 < BlobRadius2) Rho = Rho + cos(cHalfPi*sqrt(r2/BlobRadius2))**2

    if(RadialVelocity > 0)then
       ! Scale by 1/r1**(D-1) so radial flow is stationary solution
       exact_density = Rho/r1**(nDim-1)
    elseif(IsRzGeometry) then
       exact_density = Rho/Xyz_D(r_)
    else
       exact_density = Rho
    end if

  end function exact_density

  !===========================================================================
  subroutine initialize

    use BATL_lib, ONLY: init_mpi, init_batl, init_grid_batl, &
         iProc, MaxDim, MaxBlock, nBlock, Unused_B, Xyz_DGB, &
         iTree_IA, MaxLevel_, BetaProlong, IsNodeBasedGrid, IsCartesian

    use BATL_amr, ONLY: UseSimpleRefinement

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
       case("#UNIFORMAXIS")
          call read_var('UseUniformaxis', UseUniformAxis)
       case("#GRIDTYPE")
          call read_var('IsNodeBasedGrid', IsNodeBasedRead)
       case("#AMR")
          call read_var('MaxLevel', MaxLevel)
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
       case default
          call CON_stop(NameSub//' unknown command='//trim(NameCommand))
       end select
    end do READPARAM

    UseConstantVelocity = RadialVelocity == 0.0 .and. AngularVelocity == 0.0

    ! Setup generalized radus array
    if(index(TypeGeometry,'genr') > 0) Rgen_I = &
         (/ DomainMin_D(1), (DomainMin_D(1)**3*DomainMax_D(1))**0.25, DomainMax_D(1) /)

    ! Note that the periodicity will be fixed based on TypeGeometry
    call init_batl( &
         MaxBlockIn     = 8000, &          
         CoordMinIn_D   = DomainMin_D,  &
         CoordMaxIn_D   = DomainMax_D,  &
         nRootIn_D      = nRoot_D,      & 
         TypeGeometryIn = TypeGeometry, &
         IsPeriodicIn_D = (/.true., .true., .true./), &
         RgenIn_I       = Rgen_I, &
         UseUniformAxisIn = UseUniformAxis)

    ! The default value is .not.IsRzGeometry
    UseSimpleRefinement = IsCartesian

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
                  < BlobRadius2)then
                DoRefine_B(iBlock) = .true.
                CYCLE LOOPBLOCK
             end if
          end do; end do; end do
       end do LOOPBLOCK

       call init_grid_batl(DoRefine_B)
       UseSimpleRefinement = IsCartesian

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
               = exact_density(Xyz_DGB(1:nDim,i,j,k,iBlock))
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
    integer, parameter:: Volume_=nVar+1, Error_=nVar+2
    real:: TotalPe_I(Error_), Total_I(Error_)

    character(len=100):: NameFile = 'advect.log'
    logical :: DoInitialize = .true.
    !------------------------------------------------------------------------
    Total_I = 0.0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          Total_I(1:nVar) = Total_I(1:nVar) + &
               CellVolume_GB(i,j,k,iBlock)*State_VGB(1:nVar,i,j,k,iBlock)

          Total_I(Volume_) = Total_I(Volume_) + CellVolume_GB(i,j,k,iBlock)

          Total_I(Error_) = Total_I(Error_) &
               + CellVolume_GB(i,j,k,iBlock) &
               *abs(State_VGB(Rho_,i,j,k,iBlock)  &
               -    exact_density(Xyz_DGB(1:nDim,i,j,k,iBlock)))
       end do; end do; end do
    end do

    if(nProc > 1)then
       TotalPe_I = Total_I
       call MPI_reduce(TotalPe_I, Total_I, Error_, MPI_REAL, MPI_SUM, 0, &
            iComm, iError)
    end if

    if(iProc /= 0) RETURN

    ! Divide total error by total mass to get relative error
    Total_I(Error_) = Total_I(Error_) / Total_I(1)

    if(DoInitialize)then
       open(UnitTmp_,file=NameFile, status='replace')
       write(UnitTmp_,'(a)')'Advection test for BATL'
       write(UnitTmp_,'(a)')'step time mass volume error'

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
         CellVolume_GB, CellSize_DB, Xyz_DGB, CoordMin_DB, CoordMax_DB

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
                        exact_density(Xyz_DGB(:,i,j,k,iBlock)), &
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
          write(UnitTmp_,'(a)')         NameFile
          write(UnitTmp_,'(i8,a)')      nProc,        ' nProc'
          write(UnitTmp_,'(i8,a)')      iStep,        ' n_step'
          write(UnitTmp_,'(1pe13.5,a)') Time,         ' t'
          write(UnitTmp_,'(6(1pe18.10),a)') &
               (PlotMin_D(iDim),PlotMax_D(iDim),iDim=1,MaxDim),' plot_range'
          write(UnitTmp_,'(6(1pe18.10),i10,a)') &
               CellSizePlot_D, &
               CellSizeMinAll_D, nCellAll,            ' plot_dx, dxmin, ncell'
          write(UnitTmp_,'(i8,a)')     nVar+8,        ' nplotvar'
          write(UnitTmp_,'(i8,a)')     1,             ' neqpar'
          write(UnitTmp_,'(10es13.5)') 0.0            ! eqpar
          write(UnitTmp_,'(a)')        &
               'coord1 coord2 coord3 rho rhoexact volume node proc block none' 
          write(UnitTmp_,'(a)')        '1 1 1'        ! units
          write(UnitTmp_,'(l8,a)')     .true.,        ' IsBinary' 
          write(UnitTmp_,'(i8,a)')     nByteReal,     ' nByteReal'
          if(IsRLonLat)then
             write(UnitTmp_,'(a)')     'spherical'//TypeGeometry(8:20)
          else
             write(UnitTmp_,'(a)')     TypeGeometry
          end if
          if(IsGenRadius)then
             write(UnitTmp_,'(i8,a)')   nRgen,        ' nRgen'
             write(UnitTmp_,'(es13.5," LogRgen")') LogRgen_I
          end if
          write(UnitTmp_,'(a)')        'real4'        ! type of .out file
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
  subroutine set_boundary

    use BATL_lib, ONLY: nDim, IsPeriodic_D, nBlock, Unused_B, &
         DiLevelNei_IIIB, Unset_, Xyz_DGB
    ! Set ghost cells outside computational domain

    integer:: iBlock, i, j, k
    !------------------------------------------------------------------------
    if(all(IsPeriodic_D(1:nDim))) RETURN

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       if(.not.IsPeriodic_D(1))then

          ! Min in dimension 1
          if(DiLevelNei_IIIB(-1,0,0,iBlock) == Unset_)then

             do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, 0
                State_VGB(1,i,j,k,iBlock) &
                     = exact_density(Xyz_DGB(1:nDim,i,j,k,iBlock))
             end do; end do; end do

          end if
          if(DiLevelNei_IIIB(+1,0,0,iBlock) == Unset_)then
             do k = MinK, MaxK; do j = MinJ, MaxJ; do i = nI+1, MaxI
                State_VGB(1,i,j,k,iBlock) &
                     = exact_density(Xyz_DGB(1:nDim,i,j,k,iBlock))
             end do; end do; end do
          end if
       end if

       if(.not.IsPeriodic_D(2) .and. nDim > 1)then
          ! Min in dimension 2
          if(DiLevelNei_IIIB(0,-1,0,iBlock) == Unset_)then
             do k = MinK, MaxK; do j = MinJ, 0; do i = MinI, MaxI
                State_VGB(1,i,j,k,iBlock) &
                     = exact_density(Xyz_DGB(1:nDim,i,j,k,iBlock))
             end do; end do; end do
          end if
          if(DiLevelNei_IIIB(0,+1,0,iBlock) == Unset_)then
             do k = MinK, MaxK; do j = nJ+1, MaxJ; do i = MinI, MaxI
                State_VGB(1,i,j,k,iBlock) &
                     = exact_density(Xyz_DGB(1:nDim,i,j,k,iBlock))
             end do; end do; end do
          end if
       end if

       if(.not.IsPeriodic_D(3) .and. nDim > 2)then
          ! Min in dimension 3
          if(DiLevelNei_IIIB(0,0,-1,iBlock) == Unset_)then
             do k = MinK, 0; do j = MinJ, MaxJ; do i = MinI, MaxI
                State_VGB(1,i,j,k,iBlock) &
                     = exact_density(Xyz_DGB(1:nDim,i,j,k,iBlock))
             end do; end do; end do
          end if
          if(DiLevelNei_IIIB(0,0,+1,iBlock) == Unset_)then
             do k = nK+1, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
                State_VGB(1,i,j,k,iBlock) &
                     = exact_density(Xyz_DGB(1:nDim,i,j,k,iBlock))
             end do; end do; end do
          end if
       end if
    end do

  end subroutine set_boundary
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
    !------------------------------------------------------------------------

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
       call limit_slope(State_VG, Slope_VGD)

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

  end subroutine calc_face_flux
  !===========================================================================
  subroutine limit_slope(State_VG, Slope_VGD)

    ! Calculate TVD limited slopes Slope_GD of State_VG

    use ModNumConst, ONLY: i_DD

    real, intent(in) :: State_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(out):: Slope_VGD(nVar,0:nI+1,0:nJ+1,0:nK+1,nDim)

    real    :: SlopeLeft_V(nVar), SlopeRight_V(nVar)
    integer :: iDim, i, j, k, Di, Dj, Dk
    !----------------------------------------------------------------------
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

  end subroutine limit_slope
  !===========================================================================
  subroutine calc_vface_normal(iBlock)

    use ModNumConst, ONLY: i_DD
    use BATL_lib, ONLY: IsCartesian, IsRzGeometry, x_, y_, Xyz_DGB, &
         CellFace_DB, CellFace_DFB, FaceNormal_DDFB

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
                Velocity_D(x_) = -AngularVelocity*XyzFace_D(y_)
                Velocity_D(y_) =  AngularVelocity*XyzFace_D(x_)
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

    logical, parameter:: DoTest = .false.
    character(len=*), parameter:: NameSub = 'advance_explicit'
    !--------------------------------------------------------------------------
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

       if(DoTest)write(*,*)NameSub,' advance_expl iProc, iStage=',iProc,iStage

       call set_boundary

       call timing_start('message_pass')
       call message_pass_cell(nVar, State_VGB, &
            DoSendCornerIn=.true., nProlongOrderIn=2)
       call timing_stop('message_pass')

       if(DoTest)write(*,*)NameSub,' finished message_pass iProc=',iProc

       call timing_start('update')
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          call calc_face_flux(iBlock)

          DtStage = (Dt*iStage)/nOrder

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

       call apply_flux_correction(nVar, 0, State_VGB, &
            Flux_VXB=Flux_VXB, Flux_VYB=Flux_VYB, Flux_VZB=Flux_VZB)

       call timing_stop('update')

    end do

  end subroutine advance_explicit
  !===========================================================================
  subroutine advance_localstep

    use BATL_lib, ONLY: message_pass_cell, MaxBlock, nBlock, Unused_B, &
         IsCartesian, CellSize_DB, CellVolume_B, CellVolume_GB, &
         iComm, nProc, store_face_flux, apply_flux_correction
    use ModNumConst, ONLY: i_DD
    use ModMpi

    integer:: nStage, iStage, iStageBlock, iDim, iBlock
    integer:: i, j, k, Di, Dj, Dk, iError
    real :: DtMin, DtMax, DtMinPe, DtMaxPe
    real :: TimeStage, DtLocal, InvVolume
    real, save, allocatable:: Dt_B(:), Time_B(:), TimeOld_B(:)
    integer, save, allocatable:: iStage_B(:)

    logical, parameter:: DoTest = .false.
    character(len=*), parameter:: NameSub = 'advance_localstep'
    !-----------------------------------------------------------------------
    if(.not.allocated(Dt_B))then
       allocate(Dt_B(MaxBlock), Time_B(MaxBlock), TimeOld_B(MaxBlock), &
            iStage_B(MaxBlock))
       iStage_B  = 1
       Time_B    = Time
       TimeOld_B = Time
       Dt_B      = 0.0
    end if

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

       call timing_start('message_pass')
       call message_pass_cell(nVar, State_VGB, &
            DoSendCornerIn=.true., nProlongOrderIn=2, &
            TimeOld_B=TimeOld_B, Time_B=Time_B)

       call timing_stop('message_pass')

       call set_boundary

       call timing_start('update')
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          if(Time_B(iBlock) > TimeStage*(1+1e-10)) CYCLE

          iStageBlock = iStage_B(iBlock)

          call calc_face_flux(iBlock)

          DtLocal =  Dt_B(iBlock) * iStageBlock/2.0

          if(iStageBlock==2) call store_face_flux(iBlock, nVar, Flux_VFD, &
               Flux_VXB, Flux_VYB, Flux_VZB, &
               DtIn = DtLocal, &
               DoStoreCoarseFluxIn = .true.)

          if(iStageBlock == 1)then
             StateOld_VCB(:,:,:,:,iBlock) = State_VGB(:,1:nI,1:nJ,1:nK,iBlock)
          else
             State_VGB(:,1:nI,1:nJ,1:nK,iBlock) = StateOld_VCB(:,:,:,:,iBlock)
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

       end do

       ! Apply flux correction at the end of full steps at level n-1 or below
       if(modulo(iStage,4) == 0) call apply_flux_correction(nVar,0,State_VGB, &
            Flux_VXB=Flux_VXB, Flux_VYB=Flux_VYB, Flux_VZB=Flux_VZB, &
            iStageIn = iStage/2)

       TimeStage = TimeStage + DtMin/2

       call timing_stop('update')
    end do

    if( maxval(Time_B, MASK=.not.Unused_B) > &
         minval(Time_B, MASK=.not.Unused_B) + 1e-10)then

       write(*,*)'ERROR: minval(Time_B)=',minval(Time_B, MASK=.not.Unused_B)
       write(*,*)'ERROR: maxval(Time_B)=',maxval(Time_B, MASK=.not.Unused_B)
       write(*,*)'ERROR: minloc(Time_B)=',minloc(Time_B, MASK=.not.Unused_B)
       write(*,*)'ERROR: maxloc(Time_B)=',maxloc(Time_B, MASK=.not.Unused_B)

       call CON_stop('time step problem')

    end if

  end subroutine advance_localstep

end program advect

!=============================================================================
subroutine CON_stop(String)
  use BATL_lib, ONLY: iProc
  use ModMpi, ONLY: MPI_abort, MPI_COMM_WORLD
  implicit none
  integer:: iError, nError
  character (len=*), intent(in) :: String
  !--------------------------------------------------------------------------
  write(*,*)'CON_stop called on processor ',iProc,' with String='
  write(*,*) String

  call MPI_abort(MPI_COMM_WORLD, nError, iError)
  stop
end subroutine CON_stop
