!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_wrapper
  use ModConst, ONLY: cDegToRad
  use CON_coupler, ONLY: OH_,IH_,PT_,SC_, Couple_CC, Grid_C, &
       iCompSourceCouple, set_coord_system
  use CON_time
  implicit none
  SAVE

  private ! except

  ! Coupling with CON
  public:: PT_set_param
  public:: PT_init_session
  public:: PT_run
  public:: PT_save_restart
  public:: PT_finalize

  ! Point coupling
  public:: PT_get_grid_info
  public:: PT_find_points

  ! GM coupling
  public:: PT_put_from_gm

  ! OH coupling
  public:: PT_put_from_oh
  public:: PT_get_for_oh

  ! IH coupling
  public:: PT_put_from_ih

  ! SC coupling
  public:: PT_put_from_sc

  ! codes describing status of coupling with the SWMF components (OH, Ih, Sc)
  integer:: IhCouplingCode
  integer:: OhCouplingCode
  integer:: ScCouplingCode

  ! coupling operation counter (need for debugging)
  integer:: nRecvFromOH=0
  integer:: nSentToOH=0

  !----------------------------Coupling with field lines ----------------------
  ! Coupling via field line grid  with MHD components
  public:: PT_do_extract_lines
  public:: PT_put_coupling_param
  public:: PT_adjust_lines

  ! Parameters for coupling to MHD via moving lagrangian grid
  real             :: DataInputTime = 0.0, PTTime = 0.0
  ! MHD data array MHData_VIB(LagrID_:nMHData, 1:nVertexMax, 1:nLine)
  real,    pointer :: MHData_VIB(:, :, :)
  ! Number of actally used grid vertexes per each line, nVertex_B(1:nLine)
  integer, pointer :: nVertex_B(:)
  ! Grid:
  ! Mxx point number on the magnetic field line
  integer          :: nVertexMax=20000
  ! Dimensions of the grid formed by the line intersections with the spherical
  ! "origin" surface nLon*nLat, uniform in latitude grid:
  integer          :: nLon = 4, nLat = 4
  ! The radius of said origin surface, in UnitX as used in PT_set_param:
  real             :: ROrigin = 2.50
  ! Size of angular grid, in latitude and longitude, at origin
  ! surface R=ROrigin, in radians
  real             :: LonMin = -10.0*cDegToRad
  real             :: LonMax =  10.0*cDegToRad
  real             :: LatMin =  25.0*cDegToRad
  real             :: LatMax =  90.0*cDegToRad
  logical          :: DoCheck = .true., DoInit = .true.
contains
  !============================================================================
  subroutine PT_set_param(CompInfo, TypeAction)

    use CON_comp_info
    use ModReadParam
    use CON_bline,   ONLY: BL_set_grid, UseBLine_C
    use ModConst,    ONLY: rSun
    use CON_physics, ONLY: get_time
    
    
    ! Arguments
    type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp
    character (len=*), intent(in)     :: TypeAction ! What to do

    ! Contains the PARAM.in segment
    character(len=lStringLine), allocatable :: StringLineF_I(:)

    integer :: iComm, iProc, nProc, nThread, iTrue, nVar
    character(len=20):: NameVar
    
    character(len=*), parameter:: NameSub = 'PT_set_param'
    !--------------------------------------------------------------------------
    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,&
            Use        =.true., &
            NameVersion='PARMISAN', &
            Version    =1.0)

    case('MPI')
       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
    case('CHECK')
       ! AMPS could check now the input parameters for consistency
       if(UseBLine_C(PT_).and.DoCheck)then
          DoCheck = .false.   !To do this only once
          call get_time(tSimulationOut = PTTime)
          DataInputTime = PTTime
       end if
    case('READ')
       ! get section of PARAM.in that contains the PT module
    case('STDOUT')
    case('FILEOUT')
    case('GRID')
       ! Grid info depends on BATSRUS
       if(UseBLine_C(PT_))then
          !
          ! Test version;  UnitX=1.0. If the choice of UnitX, which is the
          ! unit of coordinates in the BL coupler, is hardwired, UnitX
          ! should be set to rSun from
          !
          ! use ModConst, ONLY: rSun
          !
          ! Otherwise, it should be set to some value, to be read and
          ! provided by PT/AMPS
          !
          call BL_set_grid(TypeCoordSystem='HGR', UnitX=rSun)
       end if
    case default
       call CON_stop(NameSub//': PT_ERROR: unknown TypeAction='//TypeAction)
    end select

  end subroutine PT_set_param
  !============================================================================

  subroutine PT_init_session(iSession, TimeSimulation)
    use CON_bline,  ONLY: BL_init, UseBLine_C, BL_get_origin_points
    !INPUT PARAMETERS:
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    integer::code

    character(len=*), parameter:: NameSub = 'PT_init_session'
    !--------------------------------------------------------------------------
    if(UseBLine_C(PT_).and.DoInit)then
       DoInit = .false.   ! Do this only once
       nullify(MHData_VIB); nullify(nVertex_B)
       call BL_init(nVertexMax, nLon, nLat,  &
            MHData_VIB, nVertex_B)
       call BL_get_origin_points(ROrigin, LonMin, LonMax, LatMin, LatMax)
    end if
  end subroutine PT_init_session
  !============================================================================

  subroutine PT_finalize(TimeSimulation)
    use CON_bline, ONLY: save_mhd
    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter:: NameSub = 'PT_finalize'
    !--------------------------------------------------------------------------
    if(DataInputTime > PTTime)then
       PTTime = DataInputTime
       call save_mhd(PTTime)
    end if
  end subroutine PT_finalize
  !============================================================================

  subroutine PT_save_restart(TimeSimulation)
    use CON_bline, ONLY: save_mhd
    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time
    character(len=*), parameter:: NameSub = 'PT_save_restart'
    !--------------------------------------------------------------------------
    if(DataInputTime > PTTime)then
       PTTime = DataInputTime
       call save_mhd(PTTime)
    end if

  end subroutine PT_save_restart
  !============================================================================

  subroutine PT_run(TimeSimulation, TimeSimulationLimit)
    use CON_bline, ONLY:  nLine, NameVar_V, nMHData, UseBLine_C, BL_update_r
    use CON_bline, ONLY: save_mhd
    !INPUT/OUTPUT ARGUMENTS:
    real, intent(inout):: TimeSimulation   ! current time of component

    !INPUT ARGUMENTS:
    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    real:: xEarth(3)
    character(len=*), parameter:: NameSub = 'PT_run'
    integer::ForceReachingSimulationTimeLimit=0 
    !--------------------------------------------------------------------------

    ! update the location of the Earth in the coupled on the AMPS side when
    ! coupling with IH as active
    if (IhCouplingCode==1) then
       call GetEarthLocation(xEarth)
       call set_earth_locaton_hgi(xEarth)
    end if
    !
    ! if UseBLine_C(PT_), available: DataInputTime, MHData_VIB, nVertex_B
    ! stub  for amps_get_bline is put to the end of the file
    if (UseBLine_C(PT_)) then
       call BL_update_r

       if (DataInputTime > PTTime) then
          PTTime = DataInputTime
          call save_mhd(PTTime)
       end if


!      TimeSimulation = TimeSimulationLimit  !For test purpose only
!      RETURN                                !For test purpose only
    end if

  end subroutine PT_run
  !============================================================================

  subroutine PT_get_grid_info(nDimOut, iGridOut, iDecompOut)

    ! Provide information about AMPS grid. Set number of ion fluids.

    integer, intent(out):: nDimOut    ! grid dimensionality
    integer, intent(out):: iGridOut   ! grid index (increases with AMR)
    integer, intent(out):: iDecompOut ! decomposition index

    character(len=*), parameter:: NameSub = 'PT_get_grid_info'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//' is not implemented in this version')
  end subroutine PT_get_grid_info
  !============================================================================

  subroutine PT_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

    integer, intent(in) :: nDimIn                ! dimension of position vector
    integer, intent(in) :: nPoint                ! number of positions
    real,    intent(in) :: Xyz_DI(nDimIn,nPoint) ! positions
    integer, intent(out):: iProc_I(nPoint)       ! processor owning position
    character(len=*), parameter:: NameSub = 'PT_find_points'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//' is not implemented in this version')

  end subroutine PT_find_points
  !============================================================================

  subroutine PT_put_from_gm( &
       NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)

    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI
    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data

    real, intent(out), optional, allocatable:: Pos_DI(:,:) ! Position vectors

    character(len=*), parameter:: NameSub = 'PT_put_from_gm'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//' is not implemented in this version')
  end subroutine PT_put_from_gm
  !============================================================================

  subroutine PT_put_from_oh( &
       NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)

    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI
    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data
    
    real, intent(out), optional, allocatable:: Pos_DI(:,:) ! Position vectors

    character(len=*), parameter:: NameSub = 'PT_put_from_oh'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//' is not implemented in this version')
  end subroutine PT_put_from_oh
  !============================================================================

  subroutine PT_put_from_ih( &
       NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)

    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI
    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data

    real, intent(out), optional, allocatable:: Pos_DI(:,:) ! Position vectors

    character(len=*), parameter:: NameSub = 'PT_put_from_ih'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//' is not implemented in this version')
  end subroutine PT_put_from_ih
  !============================================================================

  subroutine PT_put_from_sc(NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)
    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI
    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data

    real, intent(out), optional, allocatable:: Pos_DI(:,:) ! Position vectors

    character(len=*), parameter:: NameSub = 'PT_put_from_sc'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//' is not implemented in this version')
  end subroutine PT_put_from_sc
  !============================================================================
  subroutine PT_get_for_oh(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
       Data_VI)

    ! Get data from PT to OH

    logical,          intent(in):: IsNew   ! true for new point array
    character(len=*), intent(in):: NameVar ! List of variables
    integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
    integer,          intent(in):: nDimIn  ! Dimensionality of positions
    integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

    real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
    real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array

    integer::i,j

    character(len=*), parameter:: NameSub = 'PT_get_for_oh'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//' is not implemented in this version')
  end subroutine PT_get_for_oh
  !============================================================================
  subroutine PT_do_extract_lines(DoExtract)
    ! Interface routine to be called from super-structure on all PEs
    use CON_coupler, ONLY: i_proc0, i_comm, is_proc0
    use ModMpi
    logical, intent(out):: DoExtract

    integer :: iError

    ! when restarting, line data is available, i.e. ready to couple with mh;
    ! get value at SP root and broadcast to all SWMF processors

    ! The logical to control if trace the field lines originally
    ! (DoExtract=.true.)
    ! or not should be shaped on the root PE of the PT model and broadcast
    ! over all PEs of the SWMF

    character(len=*), parameter:: NameSub = 'PT_do_extract_lines'
    !--------------------------------------------------------------------------
    if(is_proc0(PT_)) DoExtract = .true.
    call MPI_Bcast(DoExtract, 1, MPI_LOGICAL, i_proc0(PT_), i_comm(), iError)
  end subroutine PT_do_extract_lines
  !============================================================================
  subroutine PT_put_coupling_param(Source_, TimeIn)
    use CON_bline,  ONLY: Lower_
    integer,        intent(in) :: Source_
    real,           intent(in) :: TimeIn
    !--------------------------------------------------------------------------
    if(DataInputTime >= TimeIn)RETURN
    ! New coupling time, get it and save old state
    DataInputTime = TimeIn
    if(Source_==Lower_)then
       ! Do what is needed with the MHD data about to be gone
       ! call do_something_with_MHData_VIB_array
       MHData_VIB(1:, :, :) = 0.0
    else
       call CON_stop("Time in IH-PT coupling differs from that in SC-PT")
    end if
  end subroutine PT_put_coupling_param
  !============================================================================
  ! Called from coupler after the updated grid point lo<cation are
  ! received from the other component (SC, IH). Determines whether some
  ! grid points should be added/deleted
  subroutine PT_adjust_lines(Source_)
    use CON_bline,          ONLY: &
         iOffset_B, BL_adjust_lines,  Lower_, Upper_, nLine
    integer, intent(in) :: Source_

    integer :: iLine  ! Loop variable

    character(len=*), parameter:: NameSub = 'PT_adjust_lines'
    !--------------------------------------------------------------------------
    call BL_adjust_lines(Source_)
    if(Source_ == Lower_)then
       do iLine = 1, nLine
          ! Offset the array, allocated at the lagrangian grid, if needed
          ! call offset(iBlock, iOffset=iOffset_B(iBlock))
       end do
    end if
    ! Called after the grid points are received from the
    ! component, nullify offset.
    if(Source_ == Upper_)iOffset_B(1:nLine) = 0
  end subroutine PT_adjust_lines
  !============================================================================
end module PT_wrapper
!=========Interface =================
subroutine parmisan_get_bline(&
         DataInputTime, & ! real, intent in
         nVertexMax   , & ! integer, intent in
         nLine        , & ! integer, intent in
         nVertex_B    , & ! integer, dimension(1:nLine), intent in
         nMHData      , & ! integer, intent it
         NameVar_V    , & ! character(len=10), dimension(0:nMHData) intent in
         MHData_VIB)! real,intent in,dimension(0:nMHData,1:nVertexMax,1:nLine)
  use CON_coupler, ONLY:  PT_, is_proc0
  implicit none
  real,              intent(in) :: DataInputTime
  integer,           intent(in) :: nVertexMax
  integer,           intent(in) :: nLine
  integer,           intent(in) :: nVertex_B(1:nLine)
  integer,           intent(in) :: nMHData
  character(len=10), intent(in) :: NameVar_V(0:nMHData)
  real,              intent(in) :: MHData_VIB(0:nMHData, 1:nVertexMax, 1:nLine)
  !----------------------------------------------------------------------------
end subroutine parmisan_get_bline
