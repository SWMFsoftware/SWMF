!^CFG COPYRIGHT UM
!==============================================================================
module ModSatelliteFile
  use ModSize,   ONLY: MaxBlock
  use ModIO,     ONLY: MaxFile   !!! may not be needed

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

  integer, parameter :: MaxSatelliteFile=300

  real, public :: TimeSatStart_I(MaxSatelliteFile) = 0.
  real, public :: TimeSatEnd_I(MaxSatelliteFile) = 0.

  ! These variables are public for write_logfile only !!! Should be improved
  ! Names and unit numbers for satellite files
  character(len=50), public:: Satellite_name(MaxSatelliteFile)
  integer, public:: iUnitSat_I(MaxSatelliteFile) = -1
  logical, public:: Satellite_first_write(MaxSatelliteFile) = .true.

  ! current positions
  real, public:: XSatellite(MaxSatelliteFile, 3) 

  ! variables to control time output format 
  character(len=100), public :: sat_time(MaxFile) !!! should be indexed by iSat

  ! variables to write to the satellite files
  character(len=100), public :: satellite_vars(MaxFile) !!! should be indexed by iSat

  logical, public:: DoTrackSatellite_I(MaxSatelliteFile) = .false.

  ! Local variables
  character(len=100) :: NameFile_I(MaxSatelliteFile)
  logical:: IsOpen_I(MaxSatelliteFile) = .false.
  logical:: SatelliteInBLK(MaxSatelliteFile, MaxBlock)= .false.
  logical:: UseSatelliteFile(MaxSatelliteFile)   = .true.
  integer:: Satellite_Npts(MaxSatelliteFile)
  integer,public:: iCurrent_satellite_position(MaxSatelliteFile)=1
  integer:: iPEsatellite(MaxSatelliteFile)
  integer:: iBLKsatellite(MaxSatelliteFile)
  real, allocatable   :: XSatellite_traj(:, :, :)
  real, allocatable   :: Satellite_Time(:, :)

  character(len=3)  :: TypeSatCoord_I(MaxSatelliteFile)

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
         Dn_Output, Dt_Output, plot_type, NamePlotDir
    use ModUtilities, ONLY: check_dir
    use ModVarIndexes,ONLY: NamePrimitiveVar
    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand

    integer :: iFile
    character(len=100):: StringSatellite
    character (len=3) :: satellite_var, satellite_form
    character(len=*), parameter :: NameSub = 'read_satellite_parameters'
    !------------------------------------------------------------------------
    select case(NameCommand)
    case("#SATELLITE")
       call read_var('nSatellite', nSatellite)
       if(nSatellite <= 0) RETURN
       if(iProc==0) call check_dir(NamePlotDir)
       nFile = max(nFile, Satellite_ + nSatellite)
       if (nFile > MaxFile .or. nSatellite > MaxSatelliteFile)&
            call stop_mpi(&
            'The number of output files is too large in #SATELLITE:'&
            //' nFile > MaxFile .or. nSatellite > MaxSatelliteFile')

       do iFile = Satellite_+1, satellite_ + nSatellite
          call read_var('StringSatellite', StringSatellite)
          ! Satellite output frequency
          ! Note that we broke with tradition here so that the
          ! dt_output will always we read!  This may be changed
          ! in later distributions
          call read_var('DnOutput',dn_output(ifile))
          call read_var('DtOutput',dt_output(ifile))


          ! Satellite inputfile name or the satellite name
          call read_var('NameTrajectoryFile',&
               Satellite_name(ifile-satellite_))
          if(index(StringSatellite,'eqn')>0 &
               .or. index(StringSatellite,'Eqn')>0 .or. &
               index(StringSatellite,'EQN')>0 ) then
             UseSatelliteFile(ifile-satellite_) = .false.
          else
             UseSatelliteFile(ifile-satellite_) = .true.
          end if

          ! Satellite variables
          if(index(StringSatellite,'VAR')>0 .or. &
               index(StringSatellite,'var')>0 )then
             satellite_var='var'
             plot_dimensional(ifile)= index(StringSatellite,'VAR')>0
             sat_time(ifile) = 'step date'
             call read_var('NameSatelliteVars',satellite_vars(ifile))
          elseif(index(StringSatellite,'MHD')>0 .or. &
               index(StringSatellite,'mhd')>0)then
             satellite_var='mhd'
             plot_dimensional(ifile)= index(StringSatellite,'MHD')>0
             sat_time(ifile) = 'step date'
             satellite_vars(ifile)=NamePrimitiveVar//' jx jy jz'
          elseif(index(StringSatellite,'FUL')>0 .or. &
               index(StringSatellite,'ful')>0)then
             satellite_var='ful'
             plot_dimensional(ifile)= index(StringSatellite,'FUL')>0
             sat_time(ifile) = 'step date'
             satellite_vars(ifile)=&
                  NamePrimitiveVar//' b1x b1y b1z e jx jy jz'
          else
             call stop_mpi(&
                  'Variable definition (mhd,ful,var) missing' &
                  //' from StringSatellite='//StringSatellite)
          end if

          !Change by DTW, July 2007
          !Add ray-tracing variables if 'ray' is present.
          if (index(StringSatellite,'ray')>0 .or. &
               index(StringSatellite,'RAY')>0) then
             satellite_vars(ifile) = trim(satellite_vars(ifile)) // &
                  ' theta1 phi1 status theta2 phi2'
          endif

          plot_type(ifile) = "satellite"

          ! Determine the time output format to use in the 
          ! satellite files.  This is loaded by default above, 
          ! but can be input in the log_string line.
          if(index(StringSatellite,'none')>0) then
             sat_time(ifile) = 'none'
          elseif((index(StringSatellite,'step')>0) .or. &
               (index(StringSatellite,'date')>0) .or. &
               (index(StringSatellite,'time')>0)) then
             sat_time(ifile) = ''
             if(index(StringSatellite,'step')>0) &
                  sat_time(ifile) = 'step'
             if(index(StringSatellite,'date')>0) &
                  write(sat_time(ifile),'(a)') &
                  sat_time(ifile)(1:len_trim(sat_time(ifile)))&
                  //' date'
             if(index(StringSatellite,'time')>0) &
                  write(sat_time(ifile),'(a)') &
                  sat_time(ifile)(1:len_trim(sat_time(ifile)))&
                  //' time'
          end if

       end do
    case('#STEADYSTATESATELLITE')
       do iFile = 1, nSatellite
          call read_var('SatelliteTimeStart', TimeSatStart_I(ifile))
          call read_var('SatelliteTimeEnd',   TimeSatEnd_I(ifile))
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
       l1 = index(Satellite_name(iSat), '/', back=.true.) + 1
       l2 = index(Satellite_name(iSat), '.') - 1
       if (l1-1 <= 0) l1 = 1
       if (l2+1 <= 0) l2 = len_trim(Satellite_name(iSat))

       if(n_step < 1000000)then
          write(NameFile_I(iSat),'(a,i6.6,a)')trim(NamePlotDir)//&
               'sat_'//Satellite_Name(iSat)(l1:l2)//'_n',n_step,'.sat'
       else
          write(NameFile_I(iSat),'(a,i8.8,a)')trim(NamePlotDir)//&
               'sat_'//Satellite_Name(iSat)(l1:l2)//'_n',n_step,'.sat'
       end if
       if(DoTestMe) then
          write(*,*) NameSub,': satellitename:', &
               Satellite_name(iSat), 'status =', TypeStatus
          write(*,*) 'iSat,l1,l2: ', iSat, l1, l2
          write(*,*) NameSub,': NameFile_I(iSat):', NameFile_I(iSat)
       end if

       iUnitSat_I(iSat) = io_unit_new()
       open(iUnitSat_I(iSat), file=trim(NameFile_I(iSat)), status='replace')

       IsOpen_I(iSat) = .true.
    case('append')
       if(.not.IsOpen_I(iSat))then
          iUnitSat_I(iSat) = io_unit_new()
          open(iUnitSat_I(iSat), file=trim(NameFile_I(iSat)), status='old', &
               POSITION='append')
          IsOpen_I(iSat) = .true.
       end if
    case('close')
       if (IsOpen_I(iSat)) close(iUnitSat_I(iSat))
       IsOpen_I(iSat) = .false.
    case default
       call stop_mpi(NameSub//': unknown TypeStatus='//TypeStatus)
    end select

  end subroutine set_satellite_file_status
  !============================================================================

  subroutine read_satellite_input_files

    use ModProcMH,      ONLY: iProc, iComm
    use ModMain,        ONLY: nDim, lVerbose, TypeCoordSystem, StartTime
    use CON_axes,       ONLY: transform_matrix
    use ModTimeConvert, ONLY: time_int_to_real
    use ModIo,          ONLY: iUnitOut, write_prefix
    use ModIoUnit,      ONLY: UnitTmp_
    use ModKind,        ONLY: Real8_
    use ModMpi

    integer :: iError, i, iSat , nPoint

    ! One line of input
    character (len=100) :: line

    integer      :: iTime_I(7)
    real         :: Xyz_D(nDim)
    real(Real8_) :: DateTime
    integer      :: MaxPoint
    real, allocatable:: Time_I(:), Xyz_DI(:,:)
    character(len=100):: NameFile

    character(len=*), parameter :: NameSub = 'read_satellite_input_files'

    logical :: DoTest, DoTestMe

    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    ! Count maximum number of points by reading all satellite files
    MaxPoint = 0
    if(iProc == 0)then
       SATELLITES1: do iSat=1, nSatellite
          if(.not.UseSatelliteFile(iSat)) CYCLE SATELLITES1
          NameFile = Satellite_Name(iSat)
          open(UnitTmp_, file=NameFile, status="old", iostat = iError)
          if (iError /= 0) call stop_mpi(NameSub // &
               ' ERROR1: unable to open file ' // trim(NameFile))
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
          close(UnitTmp_)
          MaxPoint = max(MaxPoint, nPoint)
       end do SATELLITES1

    end if

    ! Tell all processors the maximum number of points
    call MPI_Bcast(MaxPoint, 1, MPI_INTEGER, 0, iComm, iError)

    ! allocate arrays depending on number of points
    allocate(Time_I(MaxPoint), Xyz_DI(nDim, MaxPoint))
    allocate(XSatellite_traj(nSatellite, MaxPoint, 3))
    allocate(Satellite_Time(nSatellite, MaxPoint))

    ! Read the trajectories
    SATELLITES: do iSat=1, nSatellite

       if(.not.UseSatelliteFile(iSat)) CYCLE SATELLITES

       ! Read file on the root processor
       if (iProc == 0) then

          NameFile = Satellite_Name(iSat)

          if(lVerbose>0)then
             call write_prefix; write(iUnitOut,*) NameSub, &
                  " reading: ",trim(NameFile)
          end if

          open(UnitTmp_, file=NameFile, status="old", iostat = iError)
          if (iError /= 0) call stop_mpi(NameSub // &
               ' ERROR2: unable to open file ' // trim(NameFile))
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

          close(UnitTmp_)

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
       Satellite_Npts(iSat) = nPoint

       ! Tell the other processors the satellite time
       call MPI_Bcast(Time_I, nPoint, MPI_REAL, 0, iComm, iError)

       ! Tell the other processors the coordinates
       call MPI_Bcast(Xyz_DI, nDim*nPoint, MPI_REAL, 0, iComm, iError)

       ! Store time and positions for satellite iSat on all PE-s

       Satellite_Time(iSat, 1:nPoint) = Time_I(1:nPoint)
       do i = 1, nPoint
          xSatellite_traj(iSat, i, :) = Xyz_DI(:, i)
       end do

       if(DoTest)then
          nPoint = min(10,nPoint)
          write(*,*) NameSub,': tSat=', Satellite_Time( iSat, 1:nPoint)
          write(*,*) NameSub,': xSat=', xSatellite_traj(iSat, 1:nPoint,1)
          write(*,*) NameSub,': ySat=', xSatellite_traj(iSat, 1:nPoint,2)
          write(*,*) NameSub,': zSat=', xSatellite_traj(iSat, 1:nPoint,3)
       end if

    end do SATELLITES

    deallocate(Time_I, Xyz_DI)

  end subroutine read_satellite_input_files

  !==========================================================================

  subroutine set_satellite_flags(iSat)
    use ModProcMH
    use ModMain, ONLY : nDim,nI,nJ,nK,nBlockMax,PROCtest,unusedBLK
    use ModGeometry, ONLY : XyzStart_BLK,dx_BLK,dy_BLK,dz_BLK
    use ModGeometry, ONLY : UseCovariant               
    use ModNumConst
    use ModMpi

    integer, intent(in) :: iSat
    integer :: iPE,iPEtmp, iBLK, iBLKtemp
    real    :: xSat,ySat,zSat
    integer :: i,j,k, iError
    real,dimension(nDim)::GenOut_D
    logical :: DoTest, DoTestMe

    !--------------------------------------------------------------------------
    call set_oktest('set_satellite_flags',DoTest, DoTestMe)

    if (DoTestMe) &
         write(*,*)'Starting set_satellite_flags',&
         nSatellite,', call set_satellite_positions'

    call set_satellite_positions(iSat)
    if(.not.DoTrackSatellite_I(iSat)) RETURN !Position is not defined

    if(UseCovariant)then                  
       call xyz_to_gen(XSatellite(iSat,:),GenOut_D)
       xSat=GenOut_D(1)
       ySat=GenOut_D(2)
       zSat=GenOut_D(3)
    else                                   
       xSat=XSatellite(iSat,1)
       ySat=XSatellite(iSat,2)
       zSat=XSatellite(iSat,3)
    end if

    iPE = -1
    iBLKtemp = -1

    do iBLK = 1, nBlockMax

       SatelliteInBLK(isat,iBLK) =.not.unusedBLK(iBLK).and.&
            xSat >  XyzStart_BLK(1,iBLK) - cHalf*dx_BLK(iBLK) .and. &
            xSat <= XyzStart_BLK(1,iBLK) + (nI-cHalf)*dx_BLK(iBLK) .and. &
            ySat >  XyzStart_BLK(2,iBLK) - cHalf*dy_BLK(iBLK) .and. &
            ySat <= XyzStart_BLK(2,iBLK) + (nJ-cHalf)*dy_BLK(iBLK) .and. &
            zSat >  XyzStart_BLK(3,iBLK) - cHalf*dz_BLK(iBLK) .and. &
            zSat <= XyzStart_BLK(3,iBLK) + (nK-cHalf)*dz_BLK(iBLK)

       if(SatelliteInBLK(isat,iBLK))then 
          iPE = iProc
          iBLKtemp = iBLK
       end if

    end do

    call MPI_ALLREDUCE(iPE, iPEtmp, 1, MPI_INTEGER, MPI_MAX, iComm, iError)
    iPEsatellite(isat) = iPEtmp

    if (iPEsatellite(isat) >= 0) call MPI_Bcast(iBLKtemp, &
         1, MPI_INTEGER, iPEsatellite(isat), iComm, iError)
    iBLKsatellite(isat) = iBLKtemp

    if (iPEsatellite(isat) == -1) DoTrackSatellite_I(isat) = .false.

    if (DoTestMe) write(*,*)'set_satellite_flags (Proc',PROCtest,')(isat=', &
         isat,'): iPE,iBLK,TrackSatellite:', &
         iPEsatellite(isat),iBLKsatellite(isat),DoTrackSatellite_I(isat) 

  end subroutine set_satellite_flags

  !============================================================================

  subroutine set_satellite_positions(iSat)
    use ModProcMH
    use ModMain, ONLY : nI,nJ,nK,n_step,nBlockMax,StartTime,time_simulation, &
         time_accurate
    use ModGeometry, ONLY : x1,x2,y1,y2,z1,z2,XyzStart_BLK,dx_BLK,dy_BLK,dz_BLK
    use ModGeometry, ONLY : TypeGeometry               
    use ModNumConst

    integer, intent(in) :: iSat
    integer :: i, iBLK
    real    :: xSat,ySat,zSat
    real    :: dtime

    logical :: DoTest, DoTestMe

    !--------------------------------------------------------------------------
    if (iProc==0) &
         call set_oktest('set_satellite_positions',DoTest, DoTestMe)

    if (UseSatelliteFile(iSat)) then

       if (Satellite_Npts(iSat) > 0) then

          i = icurrent_satellite_position(iSat)

          do while ((i < Satellite_Npts(iSat)) .and.   &
               (Satellite_Time(iSat,i) < Time_Simulation))
             i = i + 1
          enddo

          icurrent_satellite_position(iSat) = i

          if ((i == Satellite_Npts(iSat).and. &
               Satellite_Time(iSat,i) <= Time_Simulation ).or.(i==1)) then 

             DoTrackSatellite_I(iSat) = .false.

          else

             DoTrackSatellite_I(iSat) = .true.

             dTime = 1.0 - (Satellite_Time(iSat,i) - Time_Simulation) / &
                  (Satellite_Time(iSat,i) - Satellite_Time(iSat,i-1) + 1.0e-6)

             xSatellite(iSat,:) = dTime * xSatellite_traj(iSat,i,:) + &
                  (1.0 - dTime) * xSatellite_traj(iSat,i-1,:) 

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

    name_string = trim(Satellite_name(iSat))
    Xvect(:) = XSatellite(iSat,:)

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
  !^CFG IF RAYTRACE BEGIN
  subroutine get_satellite_ray(iSatIn, SatRayVar_I)

    use ModRaytrace
    use ModMpi
    use ModMain,       ONLY: nI, nJ, nK, nBlock
    use ModGeometry,   ONLY: XyzStart_BLK, dx_BLK, dy_BLK, dz_BLK, UseCovariant

    integer, intent(in) :: iSatIn
    real,    intent(out):: SatRayVar_I(5)

    character(len=*), parameter :: NameSub = 'get_satellite_ray'
    integer  :: iDir, iBLK, iDim
    real     :: Xyz_D(3),  Dxyz_D(3), RayVars(3,2,nI,nJ,nK)
    real     :: Dx1, Dx2, Dy1, Dy2, Dz1, Dz2
    integer  :: i1, i2, j1, j2, k1, k2, iNear, jNear, kNear
    integer  :: i, j, k  ! Used in the FORALL statement

    !--------------------------------------------------------------------------

    ! Only use this if we're on the correct node.
    if (iProc /= iPEsatellite(iSatIn)) then
       do iDim=1,5
          SatRayVar_I(iDim) = 0.0
       enddo
       RETURN
    endif

    iBLK = iBLKSatellite(iSatIn)
    if (iBLK == 0) RETURN

    Dxyz_D(1) = dx_BLK(iBLK)
    Dxyz_D(2) = dy_BLK(iBLK)
    Dxyz_D(3) = dz_BLK(iBLK)

    if (UseCovariant) then 
       call xyz_to_gen(Xsatellite(iSatIn,:),Xyz_D)
    else
       Xyz_D = Xsatellite(iSatIn,:)
    endif

    ! Normalize coordinates
    do iDim=1,3
       Xyz_D(iDim) = (Xyz_D(iDim) - XyzStart_BLK(iDim,iBLK)) / Dxyz_D(iDim)
    end do

    ! Set location assuming point is inside block.
    i1 = floor(Xyz_D(1))
    j1 = floor(Xyz_D(2))  
    k1 = floor(Xyz_D(3))
    i2 = ceiling(Xyz_D(1))
    j2 = ceiling(Xyz_D(2))
    k2 = ceiling(Xyz_D(3))

    ! If Xyz_D is outside of block, change i,j,k in order to extrapolate.
    if(any( Xyz_D < 1) .or. any(Xyz_D > (/nI, nJ, nK/))) then
       i1 = min(nI-1, max(1, i1));   i2 = i1 + 1
       j1 = min(nJ-1, max(1, j1));   j2 = j1 + 1
       k1 = min(nK-1, max(1, k1));   k2 = k1 + 1
    endif

    ! Set interpolation weights
    Dx1 = Xyz_D(1) - i1; Dx2 = 1.0 - Dx1
    Dy1 = Xyz_D(2) - j1; Dy2 = 1.0 - Dy1
    Dz1 = Xyz_D(3) - k1; Dz2 = 1.0 - Dz1

    ! Calculate the nearest point.
    iNear = min( nI, max(nint(Xyz_D(1)),1) )
    jNear = min( nJ, max(nint(Xyz_D(2)),1) )
    kNear = min( nK, max(nint(Xyz_D(3)),1) )

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
  
  subroutine GM_trace_sat(SatXyz_D,SatRay_D)
    
    use ModProcMH,    ONLY: iComm, iProc
    use ModRayTrace,  ONLY: NameVectorField, DoExtractState, DoExtractUnitSi, &
         rIonosphere
    use ModVarIndexes,ONLY: nVar
    use ModMain,      ONLY: n_step,time_accurate,time_simulation,TypeCoordSystem
    use CON_line_extract, ONLY: line_init, line_collect, line_get, line_clean
    use ModNumConst,  ONLY: cTiny, cRadToDeg
    use CON_axes,     ONLY: transform_matrix
    use CON_planet_field, ONLY: map_planet_field
    use ModPhysics,   ONLY: rBody
    implicit none
    
    real, intent(in) :: SatXyz_D(3) ! Satellite Position
    real, intent(out)   :: SatRay_D(3)
    real :: SatXyzIono_D(3),SatXyzEnd_D(3),SatXyzEnd2_D(3),B2
    
    character(len=100) :: NameFile, NameStart, NameVar, StringTitle
    integer            :: nLineFile, nStateVar, nPlotVar
    integer            :: iPoint, nPoint, iPointNext, nPoint1
    
    real, pointer :: PlotVar_VI(:,:)
    
    integer :: iPlotFile, iLine, nLine, nVarOut, iHemisphere
    
    logical :: IsParallel = .true., IsOpen=.true.
    
    character(len=*), parameter :: NameSub = 'GM_trace_sat'
    logical :: DoTest, DoTestMe
    
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
          ! Calculate -180 < longitude = atan2(y,x) < 180
          SatRay_D(2) = cRadToDeg * atan2(SatXyzIono_D(2),SatXyzIono_D(1))
          if (SatXyzIono_D(2) < cTiny .and. SatXyzIono_D(1)<0.0) then 
             SatRay_D(2)=180.0
          endif
          if (SatXyzIono_D(1) < cTiny .and. SatXyzIono_D(2)<0.0) then 
             SatRay_D(2)=-90.0
          endif
          
          ! set closed flag
          SatRay_D(3)=3.0
       else
          SatRay_D(1)=-100.0
          SatRay_D(2)=-200.0
          SatRay_D(3)=0.0
          
       endif
    else
       !When planet is inside rBody treat as if on closed line
       SatRay_D(1)=-100.0
       SatRay_D(2)=-200.0
       SatRay_D(3)=0.0              
    endif
  end subroutine GM_trace_sat
  !============================================================================
  !^CFG END RAYTRACE
  
end module ModSatelliteFile
