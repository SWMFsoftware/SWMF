!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
module ModSatelliteFile

  use ModUtilities, ONLY: open_file, close_file

  implicit none
  save
  private !Except

  public:: read_satellite_parameters  ! read satellite file input parameters
  public:: set_satellite_file_status  ! open, open to append or close the file
  public:: read_satellite_input_files ! read satellite trajectories
  public:: set_satellite_flags        ! find the processor/block for satellite
  public:: get_satellite_ray ! map field line from satellite ^CFG IF RAYTRACE
  public:: gm_trace_sat      ! map field line from satellite ^CFG IF RAYTRACE
  
  logical, public :: DoSaveSatelliteData = .false. ! save satellite data?
  integer, public :: nSatellite = 0                ! number of satellites

  integer, parameter :: MaxSatellite=300

  real, public :: TimeSatStart_I(MaxSatellite) = 0.
  real, public :: TimeSatEnd_I(MaxSatellite) = 0.

  ! These variables are public for write_logfile only !!! Should be improved
  ! Names and unit numbers for satellite files
  character(len=50), public:: NameSat_I(MaxSatellite)
  integer, public:: iUnitSat_I(MaxSatellite) = -1
  logical, public:: IsFirstWriteSat_I(MaxSatellite) = .true.

  ! current positions
  real, public:: XyzSat_DI(3,MaxSatellite) 

  ! variables to control time output format 
  character(len=100), public :: TimeSat_I(MaxSatellite)

  ! variables to write to the satellite files
  character(len=500), public :: StringSatVar_I(MaxSatellite)

  logical, public:: DoTrackSatellite_I(MaxSatellite) = .false.

  integer,public:: iCurrent_satellite_position(MaxSatellite)=1

  ! Local variables
  character(len=100) :: NameFile_I(MaxSatellite)
  logical:: IsOpen_I(MaxSatellite) = .false.
  logical:: UseSatFile_I(MaxSatellite)   = .true.
  integer:: nPointTraj_I(MaxSatellite)
  integer:: iProcSat_I(MaxSatellite)
  integer:: iBlockSat_I(MaxSatellite)
  real, allocatable   :: XyzSat_DII(:,:,:)
  real, allocatable   :: TimeSat_II(:, :)

  character(len=3)  :: TypeSatCoord_I(MaxSatellite)

  ! Time limits (in seconds) for satellite trajectory cut 
  ! for .not. time_accurate session.
  ! If a steady-state simulation is run for a specific moment of time
  ! (set in  StartTime), the TimeSatStart_I determines the starting point of 
  ! the satellite trajectory, while TimeSatEnd_I determines the trajectory 
  ! ending point.
  ! Both determine the considered trajectory cut.
  ! Unlike in time_accurate sessions, after each dn_output simulation 
  ! steps the satellite variables for ALL the trajectory cut are 
  ! saved in file.


contains
  !============================================================================
  subroutine read_satellite_parameters(NameCommand)

    use ModProcMH,    ONLY: iProc
    use ModIO,        ONLY: nFile, MaxFile, Satellite_, plot_dimensional, &
         Dn_Output, Dt_Output, plot_type, NamePlotDir, TypeCoordPlot_I
    use ModUtilities, ONLY: check_dir
    use ModReadParam, ONLY: read_var
    use ModIO,        ONLY: NamePrimitiveVarPlot
    
    character(len=*), intent(in) :: NameCommand

    integer :: iSat, iFile
    character(len=100):: StringSatellite
    character (len=3) :: satellite_var
    character(len=*), parameter :: NameSub = 'read_satellite_parameters'
    !------------------------------------------------------------------------
    select case(NameCommand)
    case("#SATELLITE")
       call read_var('nSatellite', nSatellite)
       if(nSatellite <= 0) RETURN
       if(iProc==0) call check_dir(NamePlotDir)
       nFile = max(nFile, Satellite_ + nSatellite)
       if (nFile > MaxFile .or. nSatellite > MaxSatellite)&
            call stop_mpi(&
            'The number of output files is too large in #SATELLITE:'&
            //' nFile > MaxFile .or. nSatellite > MaxSatellite')

       do iSat = 1, nSatellite
          iFile = Satellite_ + iSat
          call read_var('StringSatellite', StringSatellite)
          ! Satellite output frequency
          ! Note that we broke with tradition here so that the
          ! dt_output will always we read!  This may be changed
          ! in later distributions
          call read_var('DnOutput', dn_output(iFile))
          call read_var('DtOutput', dt_output(iFile))

          ! Satellite inputfile name or the satellite name
          call read_var('NameTrajectoryFile',&
               NameSat_I(iSat))
          if(index(StringSatellite,'eqn')>0 &
               .or. index(StringSatellite,'Eqn')>0 .or. &
               index(StringSatellite,'EQN')>0 ) then
             UseSatFile_I(iSat) = .false.
          else
             UseSatFile_I(iSat) = .true.
          end if

          ! Satellite variables
          if(index(StringSatellite,'VAR')>0 .or. &
               index(StringSatellite,'var')>0 )then
             satellite_var='var'
             plot_dimensional(iFile) = index(StringSatellite,'VAR')>0
             TimeSat_I(iSat) = 'step date'
             call read_var('NameSatelliteVars',StringSatVar_I(iSat))
          elseif(index(StringSatellite,'MHD')>0 .or. &
               index(StringSatellite,'mhd')>0)then
             satellite_var='mhd'
             plot_dimensional(iFile)= index(StringSatellite,'MHD')>0
             TimeSat_I(iSat) = 'step date'
             StringSatVar_I(iSat)=NamePrimitiveVarPlot//' jx jy jz'
          elseif(index(StringSatellite,'FUL')>0 .or. &
               index(StringSatellite,'ful')>0)then
             satellite_var='ful'
             plot_dimensional(ifile)= index(StringSatellite,'FUL')>0
             TimeSat_I(iSat) = 'step date'
             StringSatVar_I(iSat)=&
                  NamePrimitiveVarPlot//' b1x b1y b1z e jx jy jz'
          else
             call stop_mpi(&
                  'Variable definition (mhd,ful,var) missing' &
                  //' from StringSatellite='//StringSatellite)
          end if

          !Change by DTW, July 2007
          !Add ray-tracing variables if 'ray' is present.
          if (index(StringSatellite,'ray')>0 .or. &
               index(StringSatellite,'RAY')>0) then
             StringSatVar_I(iSat) = trim(StringSatVar_I(iSat)) // &
                  ' theta1 phi1 status theta2 phi2'
          endif

          plot_type(iFile) = "satellite"

          ! Determine the time output format to use in the 
          ! satellite files.  This is loaded by default above, 
          ! but can be input in the log_string line.
          if(index(StringSatellite,'none')>0) then
             TimeSat_I(iSat) = 'none'
          elseif((index(StringSatellite,'step')>0) .or. &
               (index(StringSatellite,'date')>0) .or. &
               (index(StringSatellite,'time')>0)) then
             TimeSat_I(iSat) = ''
             if(index(StringSatellite,'step')>0) &
                  TimeSat_I(iSat) = 'step'
             if(index(StringSatellite,'date')>0) &
                  TimeSat_I(iSat) = trim(TimeSat_I(iSat))//' date'
             if(index(StringSatellite,'time')>0) &
                  TimeSat_I(iSat) = trim(TimeSat_I(iSat))//' time'
          end if

          ! Recognize coordinate system name if present
          if (index(StringSatellite,'GEO') > 0) TypeCoordPlot_I(iFile) = 'GEO'
          if (index(StringSatellite,'GSE') > 0) TypeCoordPlot_I(iFile) = 'GSE'
          if (index(StringSatellite,'GSM') > 0) TypeCoordPlot_I(iFile) = 'GSM'
          if (index(StringSatellite,'MAG') > 0) TypeCoordPlot_I(iFile) = 'MAG'
          if (index(StringSatellite,'SMG') > 0) TypeCoordPlot_I(iFile) = 'SMG'
          if (index(StringSatellite,'HGR') > 0) TypeCoordPlot_I(iFile) = 'HGR'
          if (index(StringSatellite,'HGI') > 0) TypeCoordPlot_I(iFile) = 'HGI'
          if (index(StringSatellite,'HGC') > 0) TypeCoordPlot_I(iFile) = 'HGC'

       end do
    case('#STEADYSTATESATELLITE')
       do iSat = 1, nSatellite
          call read_var('SatelliteTimeStart', TimeSatStart_I(iSat))
          call read_var('SatelliteTimeEnd',   TimeSatEnd_I(iSat))
       end do
    case default
       call stop_mpi(NameSub//' unknown command='//NameCommand)
    end select
  end subroutine read_satellite_parameters
  !============================================================================
  subroutine set_satellite_file_status(iSat,TypeStatus)

    use ModMain,   ONLY: n_step
    use ModIoUnit, ONLY: io_unit_new
    use ModIO,     ONLY: NamePlotDir

    integer, intent(in) :: iSat
    character(LEN=*),intent(in) :: TypeStatus

    integer :: l1, l2
    logical :: DoTest, DoTestMe

    character(len=*), parameter:: NameSub='set_satellite_file_status'
    !------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    select case(TypeStatus)
    case('open')
       l1 = index(NameSat_I(iSat), '/', back=.true.) + 1
       l2 = index(NameSat_I(iSat), '.') - 1
       if (l1-1 <= 0) l1 = 1
       if (l2+1 <= 0) l2 = len_trim(NameSat_I(iSat))

       if(n_step < 1000000)then
          write(NameFile_I(iSat),'(a,i6.6,a)')trim(NamePlotDir)//&
               'sat_'//NameSat_I(iSat)(l1:l2)//'_n',n_step,'.sat'
       else
          write(NameFile_I(iSat),'(a,i8.8,a)')trim(NamePlotDir)//&
               'sat_'//NameSat_I(iSat)(l1:l2)//'_n',n_step,'.sat'
       end if
       if(DoTestMe) then
          write(*,*) NameSub,': satellitename:', &
               NameSat_I(iSat), 'status =', TypeStatus
          write(*,*) 'iSat,l1,l2: ', iSat, l1, l2
          write(*,*) NameSub,': NameFile_I(iSat):', trim(NameFile_I(iSat))
       end if

       iUnitSat_I(iSat) = io_unit_new()
       call open_file(iUnitSat_I(iSat), file=NameFile_I(iSat))

       IsOpen_I(iSat) = .true.
    case('append')
       if(.not.IsOpen_I(iSat))then
          iUnitSat_I(iSat) = io_unit_new()
          call open_file(iUnitSat_I(iSat), FILE=trim(NameFile_I(iSat)), &
               STATUS='old', POSITION='append')
          IsOpen_I(iSat) = .true.
       end if
    case('close')
       if (IsOpen_I(iSat)) call close_file(iUnitSat_I(iSat))
       IsOpen_I(iSat) = .false.
    case default
       call stop_mpi(NameSub//': unknown TypeStatus='//TypeStatus)
    end select

  end subroutine set_satellite_file_status
  !============================================================================

  subroutine read_satellite_input_files

    use ModProcMH,      ONLY: iProc, iComm
    use ModMain,        ONLY: MaxDim, lVerbose, TypeCoordSystem, StartTime
    use CON_axes,       ONLY: transform_matrix
    use ModTimeConvert, ONLY: time_int_to_real
    use ModIoUnit,      ONLY: UnitTmp_
    use ModIo,          ONLY: iUnitOut, write_prefix
    use ModKind,        ONLY: Real8_
    use ModMpi

    integer :: iError, i, iSat , nPoint

    ! One line of input
    character (len=100) :: line

    integer      :: iTime_I(7)
    real         :: Xyz_D(MaxDim)
    real(Real8_) :: DateTime
    integer      :: MaxPoint
    real, allocatable:: Time_I(:), Xyz_DI(:,:)
    character(len=100):: NameFile

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'read_satellite_input_files'
    !--------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    ! Count maximum number of points by reading all satellite files
    MaxPoint = 0
    if(iProc == 0)then
       SATELLITES1: do iSat=1, nSatellite
          if(.not.UseSatFile_I(iSat)) CYCLE SATELLITES1
          NameFile = NameSat_I(iSat)
          call open_file(file=NameFile, status="old")
          nPoint = 0

          TypeSatCoord_I(iSat) = TypeCoordSystem
          READFILE1: do
             read(UnitTmp_,'(a)', iostat = iError ) line
             if (iError /= 0) EXIT READFILE1
             if(index(line,'#START')>0)then
                READPOINTS1: do
                   read(UnitTmp_,*, iostat=iError) iTime_I, Xyz_D
                   if (iError /= 0) EXIT READFILE1
                   ! Add new point
                   nPoint = nPoint + 1
                end do READPOINTS1
             end if
          end do READFILE1
          call close_file
          MaxPoint = max(MaxPoint, nPoint)
       end do SATELLITES1

    end if

    ! Tell all processors the maximum number of points
    call MPI_Bcast(MaxPoint, 1, MPI_INTEGER, 0, iComm, iError)

    ! allocate arrays depending on number of points
    allocate(Time_I(MaxPoint), Xyz_DI(MaxDim, MaxPoint))
    allocate(XyzSat_DII(3, nSatellite, MaxPoint))
    allocate(TimeSat_II(nSatellite, MaxPoint))

    ! Read the trajectories
    SATELLITES: do iSat=1, nSatellite

       if(.not.UseSatFile_I(iSat)) CYCLE SATELLITES

       ! Read file on the root processor
       if (iProc == 0) then

          NameFile = NameSat_I(iSat)

          if(lVerbose>0)then
             call write_prefix; write(iUnitOut,*) NameSub, &
                  " reading: ",trim(NameFile)
          end if

          call open_file(file=NameFile, status="old")
          nPoint = 0

          ! Read the file: read #COOR TypeCoord, #START and points
          ! Default coordinate system is the one used by BATSRUS (or GSM?)
          TypeSatCoord_I(iSat) = TypeCoordSystem
          READFILE: do

             read(UnitTmp_,'(a)', iostat = iError ) line

             if (iError /= 0) EXIT READFILE

             if(index(line,'#COOR')>0) &
                  read(UnitTmp_,'(a)') TypeSatCoord_I(iSat)

             if(index(line,'#START')>0)then

                READPOINTS: do

                   read(UnitTmp_,*, iostat=iError) iTime_I, Xyz_D

                   if (iError /= 0) EXIT READFILE

                   ! Add new point
                   nPoint = nPoint + 1

                   ! Store coordinates
                   Xyz_DI(:,nPoint) = Xyz_D

                   ! Convert integer date/time to simulation time
                   call time_int_to_real(iTime_I, DateTime)
                   Time_I(nPoint) = DateTime - StartTime

                enddo READPOINTS

             endif

          enddo READFILE

          call close_file

          if(DoTest)write(*,*) NameSub,': nPoint=',nPoint

          ! Convert the coordinates if necessary
          if(TypeSatCoord_I(iSat) /= TypeCoordSystem)then
             do i = 1, nPoint
                Xyz_DI(:,i) = matmul( &
                     transform_matrix( Time_I(i), &
                     TypeSatCoord_I(iSat), TypeCoordSystem), Xyz_DI(:,i) )
             end do
          end if

       end if

       ! Tell the number of points to the other processors
       call MPI_Bcast(nPoint, 1, MPI_INTEGER, 0, iComm, iError)
       nPointTraj_I(iSat) = nPoint

       ! Tell the other processors the satellite time
       call MPI_Bcast(Time_I, nPoint, MPI_REAL, 0, iComm, iError)

       ! Tell the other processors the coordinates
       call MPI_Bcast(Xyz_DI, MaxDim*nPoint, MPI_REAL, 0, iComm, iError)

       ! Store time and positions for satellite iSat on all PE-s

       TimeSat_II(iSat, 1:nPoint) = Time_I(1:nPoint)
       do i = 1, nPoint
          XyzSat_DII(:, iSat, i) = Xyz_DI(:, i)
       end do

       if(DoTest)then
          nPoint = min(10,nPoint)
          write(*,*) NameSub,': tSat=', TimeSat_II(iSat,1:nPoint)
          write(*,*) NameSub,': xSat=', XyzSat_DII(1,iSat,1:nPoint)
          write(*,*) NameSub,': ySat=', XyzSat_DII(2,iSat,1:nPoint)
          write(*,*) NameSub,': zSat=', XyzSat_DII(3,iSat,1:nPoint)
       end if

    end do SATELLITES

    deallocate(Time_I, Xyz_DI)

  end subroutine read_satellite_input_files

  !==========================================================================

  subroutine set_satellite_flags(iSat)

    use BATL_lib, ONLY: find_grid_block

    integer, intent(in) :: iSat
    logical :: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'set_satellite_flags'
    !--------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    if (DoTestMe) write(*,*) NameSub,' starting for iSat=', iSat

    call set_satellite_positions(iSat)

    if(DoTestMe)write(*,*) NameSub,' DoTrackSatellite=', &
         DoTrackSatellite_I(iSat)
    if(.not.DoTrackSatellite_I(iSat)) RETURN !Position is not defined

    call find_grid_block(XyzSat_DI(:,iSat), &
         iProcSat_I(iSat), iBlockSat_I(iSat))

    if (iProcSat_I(iSat) < 0) DoTrackSatellite_I(iSat) = .false.

    if (DoTestMe) write(*,*)'set_satellite_flags iPE,iBLK,TrackSatellite=', &
         iProcSat_I(iSat), iBlockSat_I(iSat), DoTrackSatellite_I(iSat) 

  end subroutine set_satellite_flags

  !============================================================================

  subroutine set_satellite_positions(iSat)
    use ModProcMH
    use ModMain, ONLY :time_simulation
    use ModNumConst

    integer, intent(in) :: iSat
    integer :: i
    real    :: dtime

    logical :: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'set_satellite_positions'
    !--------------------------------------------------------------------------
    call set_oktest('set_satellite_positions',DoTest, DoTestMe)

    if (UseSatFile_I(iSat)) then

       if (nPointTraj_I(iSat) > 0) then

          i = icurrent_satellite_position(iSat)

          if(DoTestMe)write(*,*) NameSub,' nPoint, iPoint, TimeSim, TimeSat=',&
               nPointTraj_I(iSat), i, Time_Simulation, TimeSat_II(iSat,i)

          do while ((i < nPointTraj_I(iSat)) .and.   &
               (TimeSat_II(iSat,i) <= Time_Simulation))
             i = i + 1
          enddo

          icurrent_satellite_position(iSat) = i

          if(DoTestMe)write(*,*) NameSub,' final iPoint=', i

          if (i == nPointTraj_I(iSat) .and. &
               TimeSat_II(iSat,i) <= Time_Simulation .or. i==1) then 

             DoTrackSatellite_I(iSat) = .false.
             XyzSat_DI(:,iSat) = 0.0

             if(DoTestMe)write(*,*) NameSub,' DoTrackSat=.false.'

          else

             DoTrackSatellite_I(iSat) = .true.

             dTime = 1.0 - (TimeSat_II(iSat,i) - Time_Simulation) / &
                  (TimeSat_II(iSat,i) - TimeSat_II(iSat,i-1) + 1.0e-6)

             XyzSat_DI(:,iSat) = dTime * XyzSat_DII(:,iSat,i) + &
                  (1.0 - dTime) * XyzSat_DII(:,iSat,i-1) 

             if(DoTestMe)write(*,*) NameSub,' XyzSat=', XyzSat_DI(:,iSat)

          endif

       endif

    else 

       call satellite_trajectory_formula(iSat)

    end if

  end subroutine set_satellite_positions

  !============================================================================

  subroutine satellite_trajectory_formula(iSat)

    integer, intent(in) :: iSat
    character (len=100) :: name_string
    real :: Xvect(3)
    !-------------------------------------------------------------------------
    name_string = trim(NameSat_I(iSat))
    Xvect(:) = XyzSat_DI(:,iSat)

    ! Case should be for a specific satellite.  The trajectories can depend
    ! on the 'real' time so that the satellite knows where it is at.  For
    ! example, Cassini could be if'd on the date so that the code knows 
    ! whether the run is for a time near Earth, Jupiter or Saturn. 

    ! This routine should always set the TrackSatellite(iSat) flag. When
    ! the satellite is at a useless position the time should return a
    ! do not track flag (DoTrackSatellite_I(iSat) = .false.).

    select case(name_string)
    case ('earth')
       Xvect(1) = 5.0
       Xvect(2) = 5.0
       Xvect(3) = 5.0
       DoTrackSatellite_I(iSat) = .true.
    case ('cassini')
       Xvect(1) = 5.0
       Xvect(2) = 5.0
       Xvect(3) = 5.0
       DoTrackSatellite_I(iSat) = .true.
    case default 
       Xvect(1) = 1.0
       Xvect(2) = 1.0
       Xvect(3) = 1.0
       DoTrackSatellite_I(iSat) = .false.
    end select

  end subroutine satellite_trajectory_formula

  !============================================================================

  subroutine get_satellite_ray(iSatIn, SatRayVar_I)

    use ModRaytrace
    use ModMpi
    use BATL_size,ONLY: nI, nJ, nK
    use BATL_lib, ONLY: IsCartesianGrid, CellSize_DB, CoordMin_DB, xyz_to_coord

    integer, intent(in) :: iSatIn
    real,    intent(out):: SatRayVar_I(5)

    character(len=*), parameter :: NameSub = 'get_satellite_ray'
    integer  :: iDir, iBLK, iDim
    real     :: Coord_D(3), RayVars(3,2,nI,nJ,nK)
    real     :: Dx1, Dx2, Dy1, Dy2, Dz1, Dz2
    integer  :: i1, i2, j1, j2, k1, k2, iNear, jNear, kNear
    integer  :: i, j, k

    !--------------------------------------------------------------------------
    
    ! Initialize to zero
    SatRayVar_I = 0.0

    ! Only use this if we're on the correct node.
    if (iProc /= iProcSat_I(iSatIn)) RETURN

    iBLK = iBlockSat_I(iSatIn)
    if (iBLK == 0) RETURN

    if (IsCartesianGrid) then 
       Coord_D = XyzSat_DI(:,iSatIn)
    else
       call xyz_to_coord(XyzSat_DI(:,iSatIn), Coord_D)
    end if

    ! Normalize coordinates to the cell center indexes
    Coord_D = (Coord_D - CoordMin_DB(:,iBLK)) / CellSize_DB(:,iBlk) + 0.5

    ! Set location assuming point is inside block.
    i1 = floor(Coord_D(1))
    j1 = floor(Coord_D(2))  
    k1 = floor(Coord_D(3))
    i2 = ceiling(Coord_D(1))
    j2 = ceiling(Coord_D(2))
    k2 = ceiling(Coord_D(3))

    ! If Coord_D is outside of block, change i,j,k in order to extrapolate.
    if(any( Coord_D < 1) .or. any(Coord_D > (/nI, nJ, nK/))) then
       i1 = min(nI-1, max(1, i1));   i2 = i1 + 1
       j1 = min(nJ-1, max(1, j1));   j2 = j1 + 1
       k1 = min(nK-1, max(1, k1));   k2 = k1 + 1
    endif

    ! Set interpolation weights
    Dx1 = Coord_D(1) - i1; Dx2 = 1.0 - Dx1
    Dy1 = Coord_D(2) - j1; Dy2 = 1.0 - Dy1
    Dz1 = Coord_D(3) - k1; Dz2 = 1.0 - Dz1

    ! Calculate the nearest point.
    iNear = min( nI, max(nint(Coord_D(1)),1) )
    jNear = min( nJ, max(nint(Coord_D(2)),1) )
    kNear = min( nK, max(nint(Coord_D(3)),1) )

    ! Copy ray tracing values to new array so allow changing of values.
    RayVars = ray(1:3,1:2,1:nI,1:nJ,1:nK,iBLK)

    ! Use the ray status of the nearest point to the satellite.
    SatRayVar_I(3) = RayVars(3, 1, iNear,jNear,kNear)

    ! For each direction along the ray, determine if all lines surrounding 
    ! point are open or closed, then set SatRayVar_I accordingly.
    do iDir=1,2

       if ( any(RayVars(3,1,i1:i2,j1:j2,k1:k2) < 1) .or. &
            any(RayVars(3,1,i1:i2,j1:j2,k1:k2) == iDir) ) then
          ! One or more lines is open in direction iDir, must use nearest point
          do iDim=1,2
             SatRayVar_I(iDim + 3*(iDir-1)) = &
                  RayVars(iDim,iDir,iNear,jNear,kNear)
          end do

       else   ! All lines closed in direction iDir, interpolate.

          ! If the satellite is near the 0/360 degree boundary in longitude,
          ! the result of the interpolation will be incorrect.  Adjust
          ! longitudes accordingly.
          if (any(RayVars(2,iDir,i1:i2,j1:j2,k1:k2)>330.0) .AND. &
               any(RayVars(2,iDir,i1:i2,j1:j2,k1:k2)<30.0)) then

             do i=i1,i2
                do j=j1,j2
                   do k=k1,k2
                      if (RayVars(2,iDir,i,j,k) < 30.0) &
                           RayVars(2,iDir,i,j,k) = RayVars(2,iDir,i,j,k) + 360
                   enddo
                enddo
             enddo
             !forall(i=i1:i2,j=j1:j2,k=k1:k2,RayVars(2,iDir,i,j,k)<30.0)
             !   RayVars(2,iDir,i,j,k) = RayVars(2,iDir,i,j,k) + 360.0
             !end forall
          endif

          do iDim=1,2
             SatRayVar_I(iDim + 3*(iDir-1)) =                         &
                  Dz2*(   Dy2*(   Dx2*RayVars(iDim,iDir,i1,j1,k1)   &
                  +                Dx1*RayVars(iDim,iDir,i2,j1,k1))  &
                  +        Dy1*(   Dx2*RayVars(iDim,iDir,i1,j2,k1)   &
                  +                Dx1*RayVars(iDim,iDir,i2,j2,k1))) &
                  +Dz1*(   Dy2*(   Dx2*RayVars(iDim,iDir,i1,j1,k2)   &
                  +                Dx1*RayVars(iDim,iDir,i2,j1,k2))  &
                  +        Dy1*(   Dx2*RayVars(iDim,iDir,i1,j2,k2)   &
                  +                Dx1*RayVars(iDim,iDir,i2,j2,k2)))
          end do
       endif

    end do

    ! Ensure that longitude wraps around 0/360 degree boundary correctly.
    if(SatRayVar_I(2) <   0.0) SatRayVar_I(2) = SatRayVar_I(2) + 360.0
    if(SatRayVar_I(5) <   0.0) SatRayVar_I(5) = SatRayVar_I(5) + 360.0
    if(SatRayVar_I(2) > 360.0) SatRayVar_I(2) = SatRayVar_I(2) - 360.0
    if(SatRayVar_I(5) > 360.0) SatRayVar_I(5) = SatRayVar_I(5) - 360.0

  end subroutine get_satellite_ray
  !============================================================================
  
  subroutine GM_trace_sat(SatXyz_D, SatRay_D)
    
    use ModProcMH,    ONLY: iComm, iProc
    use ModRayTrace,  ONLY: DoExtractState, DoExtractUnitSi, &
         rIonosphere
    use ModVarIndexes,ONLY: nVar
    use ModMain,      ONLY:time_simulation,TypeCoordSystem
    use CON_line_extract, ONLY: line_init, line_collect, line_get, line_clean
    use ModNumConst,  ONLY: cRadToDeg
    use CON_axes,     ONLY: transform_matrix
    use CON_planet_field, ONLY: map_planet_field
    use ModPhysics,   ONLY: rBody
    implicit none
    
    real, intent(in) :: SatXyz_D(3) ! Satellite Position
    real, intent(out)   :: SatRay_D(3)
    real :: SatXyzIono_D(3),SatXyzEnd_D(3),SatXyzEnd2_D(3),B2
    
    integer            :: nStateVar
    integer            :: nPoint
    
    real, pointer :: PlotVar_VI(:,:)
    
    integer :: nLine, nVarOut, iHemisphere
    
    logical :: IsParallel = .true., IsOpen=.true.
    
    character(len=*), parameter :: NameSub = 'GM_trace_sat'
    logical :: DoTest, DoTestMe
    real    :: L, Rsat, Rxy2
    ! Conversion matrix between SM and GM coordinates 
    ! (to be safe initialized to unit matrix)
    real :: GmSm_DD(3,3) = reshape( (/ &
         1.,0.,0., &
         0.,1.,0., &
         0.,0.,1. /), (/3,3/) )
    !-------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)
    
    DoExtractState = .true.
    DoExtractUnitSi= .false.
    
    ! Set the number lines and variables to be extracted
    nLine     = 1
    nStateVar = 4
    if(DoExtractState) nStateVar = nStateVar + nVar
    
    if (sum(SatXyz_D(1:3)**2) > rBody**2) then
       ! Initialize CON_line_extract
       call line_init(nStateVar)
       
       ! Obtain the line data
       call ray_lines(nLine, (/IsParallel/), SatXyz_D)
       
       ! Collect lines from all PE-s to Proc 0
       call line_collect(iComm,0)
       
       if(iProc==0)then
          call line_get(nVarOut, nPoint)
          if(nVarOut /= nStateVar)call stop_mpi(NameSub//': nVarOut error')
          allocate(PlotVar_VI(0:nVarOut, nPoint))
          call line_get(nVarOut, nPoint, PlotVar_VI, DoSort=.true.)
       end if
    
       call line_clean
       
       
       ! Only iProc 0 stores result 
       if(iProc == 0) then
          SatXyzEnd_D = PlotVar_VI(2:4,nPoint)
          deallocate(PlotVar_VI)
       endif
       
       !Now Trace in opposite direction to make sure line is closed
       ! Initialize CON_line_extract
       call line_init(nStateVar)
       
       ! Obtain the line data
       call ray_lines(nLine, (/.not.IsParallel/), SatXyz_D)
       
       ! Collect lines from all PE-s to Proc 0
       call line_collect(iComm,0)
       
       if(iProc==0)then
          call line_get(nVarOut, nPoint)
          if(nVarOut /= nStateVar)call stop_mpi(NameSub//': nVarOut error')
          allocate(PlotVar_VI(0:nVarOut, nPoint))
          call line_get(nVarOut, nPoint, PlotVar_VI, DoSort=.true.)
       end if
       
       call line_clean
       
       ! Only iProc 0 stores result 
       if(iProc /= 0) RETURN
       
       SatXyzEnd2_D = PlotVar_VI(2:4,nPoint)
       
       B2=sum(PlotVar_VI(5:7,1)**2)
       
       deallocate(PlotVar_VI)
       
       ! Only iProc 0 works for returning line info
       !  if(iProc /= 0) RETURN
       
       ! Check that line is closed
       if (sum(SatXyzEnd_D(:)**2.0) > 8.0*rIonosphere**2.0 .or. &
            sum(SatXyzEnd2_D(:)**2.0) > 8.0*rIonosphere**2.0) then
          IsOpen=.true.
       else
          IsOpen=.false.
       endif
       
       if (.not. IsOpen) then
          call map_planet_field(Time_Simulation, SatXyzEnd_D, &
               TypeCoordSystem//' NORM', rIonosphere, SatXyzIono_D, iHemisphere)
          
          ! Transformation matrix between the SM(G) and GM coordinates
          GmSm_DD = transform_matrix(time_simulation,'SMG',TypeCoordSystem)
          ! Convert GM position into RB position
          SatXyzIono_D = matmul(SatXyzIono_D, GmSm_DD)
          
          ! Convert XYZ to Lat-Lon-IsOpen
          ! Calculate  -90 < latitude = asin(z)  <  90
          SatRay_D(1) = cRadToDeg * asin(SatXyzIono_D(3)/rIonosphere)
          ! Calculate 0 < longitude = atan2(y,x) < 360
          SatRay_D(2) =  &
               modulo(cRadToDeg *atan2(SatXyzIono_D(2),SatXyzIono_D(1)),360.0)

          ! set closed flag
          SatRay_D(3)=3.0
       else
          SatRay_D(1)=-100.0
          SatRay_D(2)=-200.0
          SatRay_D(3)=0.0
          
       endif
    else
       !When planet is inside rBody use dipole assumption
       Rsat=sqrt(sum(SatXyz_D(1:3)**2))
       Rxy2=sum(SatXyz_D(1:2)**2)
       if (Rxy2>0.0) then
          L=Rsat**3.0/Rxy2
          SatRay_D(1)=acos(1/sqrt(L))*cRadToDeg
       else
          SatRay_D(1)=90.0
       endif
       
       SatRay_D(2) =  &
            modulo(cRadToDeg *atan2(SatXyz_D(2),SatXyz_D(1)),360.0)
       !set closed flag
       SatRay_D(3)=3.0              
    endif
  end subroutine GM_trace_sat
  !============================================================================

  
end module ModSatelliteFile
