!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE SP
!^CMP FILE IH
module CON_couple_mh_sp
  use CON_axes
  use CON_coupler
  use IH_wrapper, ONLY: &
       IH_check_ready_for_sp    ,&  ! If returns .false., the run stops.
       IH_synchronize_refinement,&  ! Updates the info about IH grid
       IH_extract_line          ,&  ! Makes mf lines for given set of oriigins
       IH_put_particles         ,&  ! Put "points" (coords+IDs) from SP
       IH_n_particle            ,&  ! Number of "points" on a given Proc in IH
       IH_get_particle_coords   ,&  ! Coords of a "point" for a given # in array
       IH_get_particle_indexes  ,&  ! Stored IDs of a "point" for a given #
       IH_xyz_to_coord          ,&  ! Conversion Cartesian <-> generalized
       IH_coord_to_xyz          ,&  ! coords in the IH.
       IH_get_for_mh            ,&  ! Gets MHD information from a given cell
       IH_Grid, IH_LineGrid, IH_LocalLineGrid  ! Grid descriptors
  !^CMP IF SC BEGIN
  use SC_wrapper, ONLY: SC_check_ready_for_sp, SC_synchronize_refinement,    &
       SC_extract_line, SC_put_particles, SC_n_particle,                     &
       SC_get_particle_coords, SC_get_particle_indexes, SC_xyz_to_coord,     &
       SC_get_for_mh, SC_Grid, SC_LineGrid, SC_LocalLineGrid
  !^CMP END SC
  !^CMP IF OH BEGIN
  use OH_wrapper, ONLY: OH_check_ready_for_sp, OH_synchronize_refinement,    &
       OH_extract_line, OH_put_particles, OH_n_particle,                     &
       OH_get_particle_coords, OH_get_particle_indexes, OH_xyz_to_coord,     &
       OH_coord_to_xyz, OH_get_for_mh, OH_Grid, OH_LineGrid, OH_LocalLineGrid
  !^CMP END OH
  use SP_wrapper, ONLY: &
       SP_do_extract_lines      ,&  ! If returns .true., extract the mf lines
       SP_put_coupling_param    ,&  ! Set time and interaface bounds
       SP_adjust_lines              ! Process if needed the updated mf lines
  !^CMP IF PT BEGIN
  use PT_wrapper, ONLY: &
       PT_do_extract_lines      ,&  ! If returns .true., extract the mf lines
       PT_put_coupling_param    ,&  ! Set time and interaface bounds
       PT_adjust_lines              ! Process if needed the updated mf lines
  !^CMP END PT
  use CON_bline, ONLY:  BL_     ,&  ! SP_ or PT_
       RScMin, RScMax           ,&  !^CMP IF SC
       RIhMin, RIhMax           ,&
       ROhMin, ROhMax           ,&  !^CMP IF OH
       BlMh_DD                  ,&  ! transformation matrix
       TimeBl, TypeCoordMh      ,&  ! local time; coord system of source
       IsSource4BL_C            ,&  ! list of used MH models
       BL_init_foot_points      ,&  ! Initialize footpoint array
       BL_get_bounds            ,&  ! Provides RScMin/Max and/or RIhMin/Max
       BL_n_particle            ,&  ! Number of "points" in a given line in SP
       BL_put_from_mh           ,&  ! Put MHD info from SC or IH to SP
       BL_is_interface_block    ,&  ! Skips unused block
       BL_interface_point_coords,&  ! Check if the point is within interface
       BL_put_line              ,&  ! Put particle Xyz from SC/IH to SP
       BL_Grid, BL_LocalGrid        ! Grid descriptors (global and local)

  implicit none
  private ! Except
  public::couple_mh_sp_init
  public::couple_ih_sp
  public::couple_sc_sp              !^CMP IF SC
  public::couple_oh_sp              !^CMP IF OH

  type(RouterType),save,private::RouterIhBl        ! IH (MHD data) => SP
  type(RouterType),save,private::RouterLineIhBl    ! IH (Particle coords)=>SP
  !^CMP IF SC BEGIN
  type(RouterType),save,private::RouterScBl        ! SC MHD data => BL
  type(RouterType),save,private::RouterLineScBl    ! SC (Particle coords)=>BL
  !^CMP END SC
  !^CMP IF OH BEGIN
  type(RouterType),save,private::RouterOhBl        ! OH MHD data => BL
  type(RouterType),save,private::RouterLineOhBl    ! OH (Particle coords)=>BL
  !^CMP END OH

  ! Rescaling coefficients
  real,    public    :: UnitBl2UnitMh, UnitMh2UnitBl

  ! Three-dimensional coordinate vector
  integer, parameter :: nDim = 3
  ! We send: (1) the global index of the magnetic field line
  ! and (2) the particle number along the magnetic field line
  integer, parameter :: nAux = 2
  ! Misc
  logical :: DoTest, DoTestMe
  character(LEN=*), parameter :: NameSub='couple_mh_sp'
  logical ::DoInit = .true., DoExtract = .false.
  real, save :: DataInputTimeLast = -1.0
contains
  !============================================================================
  subroutine BL_put_coupling_param(Source_, TimeIn)
    integer, intent(in) :: Source_
    real,    intent(in) :: TimeIn
    !--------------------------------------------------------------------------
    select case(BL_)
       case(SP_)
          call SP_put_coupling_param(Source_, TimeIn)
       case(PT_)
          call PT_put_coupling_param(Source_, TimeIn)
       case default
          call CON_stop('The target model it BL is not allowed')
       end select
     end subroutine BL_put_coupling_param
  !============================================================================
  subroutine BL_adjust_lines(Source_)
    integer, intent(in) :: Source_
    !--------------------------------------------------------------------------
    select case(BL_)
    case(SP_)
       call SP_adjust_lines(Source_)
    case(PT_)
       call PT_adjust_lines(Source_)
    case default
       call CON_stop('The target model it BL is not allowed')
    end select
  end subroutine BL_adjust_lines
  !============================================================================
  subroutine transform_matrix_and_coef(tNow, MH_)
    real,    intent(in) :: tNow
    integer, intent(in) ::  MH_
    ! Transformation matrix and coefficients
    !--------------------------------------------------------------------------
    TimeBl = tNow
    UnitBl2UnitMh = Grid_C(BL_)%UnitX/Grid_C(MH_)%UnitX
    UnitMh2UnitBl = 1/UnitBl2UnitMh
    TypeCoordMh   = Grid_C(MH_)%TypeCoord
    BlMh_DD       = transform_matrix(TimeSim = TimeBl, &
         TypeCoordIn = TypeCoordMh, TypeCoordOut = Grid_C(BL_)%TypeCoord)
  end subroutine transform_matrix_and_coef
  !============================================================================
  subroutine i_particle(nDim, Xyz_D, nIndex, iIndex_I, IsInterfacePoint)
    ! Just converts the number of point in the grid "lines" to real one,
    ! to be used as an input parameter for the mapping procedure
    integer,intent(in)    :: nDim, nIndex
    logical,intent(out)   :: IsInterfacePoint
    real,   intent(inout) :: Xyz_D(nDim)
    integer,intent(inout) :: iIndex_I(nIndex)
    !--------------------------------------------------------------------------
    IsInterfacePoint = .true.; Xyz_D(1)= real(iIndex_I(1))
  end subroutine i_particle
  !============================================================================
  subroutine couple_mh_sp_init
    use CON_physics, ONLY: get_time
    use ModConst
    ! whether need to extract new line (equal to .not.DoRestart in BL
    real :: tNow
    !--------------------------------------------------------------------------
    if(.not.DoInit)RETURN    ! The initialization can be done only once
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    call get_time(tSimulationOut = tNow)
    ! determine whether BL is ready for coupling with MH
    select case(BL_)
    case(SP_)
       call SP_do_extract_lines(DoExtract)
    case(PT_)
       call PT_do_extract_lines(DoExtract)
    case default
       call CON_stop('The BL model is not allowed')
    end select
    if(IsSource4BL_C(SC_))call couple_sc_sp_init(tNow)  !^CMP IF SC
    if(IsSource4BL_C(IH_))call couple_ih_sp_init(tNow)
    if(IsSource4BL_C(OH_))call couple_oh_sp_init(tNow)  !^CMP IF OH
    if(DoExtract.and.is_proc(BL_))call BL_init_foot_points
    ! After couple_mh_sp_init
    ! (a) set of the Lagrangian point coordinates is in BL
    ! (b) logical DoInit is .true. to prevent the coordinate
    ! info on BL to be overritten with that from MH models
  end subroutine couple_mh_sp_init
  !============================================================================
  subroutine couple_sc_sp_init(tNow)               !^CMP IF SC BEGIN
    real,   intent(in) :: tNow
    logical            :: IsReady
    integer, parameter :: Lower_ = 0
    integer            :: nLength
    !--------------------------------------------------------------------------
    call SC_check_ready_for_sp(IsReady)
    if(.not.IsReady) call CON_stop(&
         "SC component not ready for "//NameSub//" correct PARAM.in")
    call set_couple_var_info(SC_, BL_)
    ! Initialize coupler from SC (source )to BL (target) for sendding data
    ! from SC to BL. The points at which the data shuld be provided are sent
    ! from BL to SC, at the stage of the router construction. Two BL grid
    ! indexes are also sent from BL to SC (nMappedPointIndex = nAux = 2).
    call init_router(                             &
         GridSource        = SC_Grid             ,&
         GridTarget        = BL_Grid             ,&
         nMappedPointIndex = nAux                ,&
         Router            = RouterScBl          )
    ! Router to send advected point coordinates from SC to BL
    call init_router(SC_LineGrid, BL_Grid, RouterLineScBl)
    if(.not.RouterScBl%IsProc)RETURN
    if(.not.DoExtract)RETURN
    ! Now, in BL there are only origin points
    if(is_proc(BL_))then
       call BL_put_coupling_param(SC_, tNow)
       call BL_get_bounds(RScMin, RScMax)
    end if
    call SC_synchronize_refinement(RouterScBl%iProc0Source,&
         RouterScBl%iCommUnion)
    call  transform_matrix_and_coef(tNow, SC_)
    ! use router intended for sending MHD data      (SC => BL)
    ! to send the the origin points of the MF lines (BL => SC).
    call set_router(&
         GridSource            = SC_Grid                  ,&
         GridTarget            = BL_LocalGrid             ,&
         Router                = RouterScBl               ,&
         n_interface_point_in_block = BL_n_particle       ,&
         interface_point_coords= BL_interface_point_coords,&! Xyz_D in BL
         is_interface_block    = BL_is_interface_block    ,&! used line
         mapping               = mapping_sp_to_sc         ,&!=>Coord_D in SC
         interpolate           = interpolation_amr_gc)
    if(is_proc(SC_))then
       nLength = nlength_buffer_source(RouterScBl)
       ! Trace the MF lines in SC, down to rScMin, up to rScMax
       call SC_extract_line(&
            Xyz_DI     = RouterScBl%BufferSource_II(1:nDim,&! Coord_D in SC
            1:nLength)                                    ,&
            iTraceMode = Lower_                           ,&
            iIndex_II  = nint(RouterScBl% BufferSource_II( &! Global block
            nDim+1:nDim+nAux,1:nLength))                  ,&! and ID in BL
            RSoftBoundary = RScMax*UnitBl2UnitMh)           ! Limit line in SC
       ! Now in SC are the parts of lines from the RScMn to RScMax
       ! + 1 Point above RScMax per line
    end if
    call BL_put_lines_from_sc ! Now these lines are in BL
  end subroutine couple_sc_sp_init
  !============================================================================
  subroutine mapping_sp_to_sc(nDimIn, XyzBl_D, nDimOut, CoordSc_D, &
       IsInterfacePoint)
    integer, intent(in) :: nDimIn, nDimOut
    real,    intent(in) :: XyzBl_D(nDimIn)
    real,    intent(out):: CoordSc_D(nDimOut)
    logical, intent(out):: IsInterfacePoint
    real                :: XyzSc_D(nDim)
    !--------------------------------------------------------------------------
    IsInterfacePoint = .true.; XyzSc_D = matmul(XyzBl_D, BlMh_DD)*UnitBl2UnitMh
    call SC_xyz_to_coord(XyzSc_D, CoordSc_D)
  end subroutine mapping_sp_to_sc
  !============================================================================
  subroutine BL_put_lines_from_sc
    !--------------------------------------------------------------------------
    call construct_router_from_source(                      &
         GridSource                 = SC_LocalLineGrid     ,&
         GridTarget                 = BL_Grid              ,&
         Router                     = RouterLineScBl       ,&
         n_interface_point_in_block = SC_n_particle        ,&
         interface_point_coords     = i_particle           ,&! iParticle=>x_D(1)
         mapping                    = mapping_line_sc_to_sp) ! x_D(1)=>BLindexes
    call couple_comp(                 RouterLineScBl       ,&
         nVar                       = nDim                 ,&
         fill_buffer = SC_get_coord_for_sp_and_transform   ,&! Xyz_D for BL
         apply_buffer               = BL_put_line)           ! Put in place
    ! Only points below RMinBl are removed.
  end subroutine BL_put_lines_from_sc
  !============================================================================
  subroutine mapping_line_sc_to_sp(nDimIn, XyzIn_D, nDimOut, CoordBl_D, &
       IsInterfacePoint)
    ! See comments in mapping_line_ih_to_sp
    integer, intent(in) :: nDimIn, nDimOut
    real,    intent(in) :: XyzIn_D(nDimIn)
    real,    intent(out):: CoordBl_D(nDimOut)
    logical, intent(out):: IsInterfacePoint
    integer:: iIndex_I(2), iParticle
    !--------------------------------------------------------------------------
    IsInterfacePoint = .true.; iParticle = nint(XyzIn_D(1))
    call SC_get_particle_indexes(iParticle, iIndex_I)
    CoordBl_D = coord_grid_d(BL_Grid, iIndex_I(1),[iIndex_I(2),1,1])
  end subroutine mapping_line_sc_to_sp
  !============================================================================
  subroutine SC_get_coord_for_sp_and_transform(&
       nPartial, iGetStart, Get, w, XyzBl_D, nDimOut)
    integer,            intent(in):: nPartial, iGetStart, nDimOut
    type(IndexPtrType), intent(in):: Get
    type(WeightPtrType),intent(in):: w
    real,              intent(out):: XyzBl_D(nDimOut)
    real :: XyzSc_D(nDim)
    integer :: iParticle
    !--------------------------------------------------------------------------
    iParticle = Get%iCB_II(1,iGetStart)
    call SC_get_particle_coords(iParticle, XyzSc_D)
    ! coord transformation to BL
    XyzBl_D = matmul(BlMh_DD, XyzSc_D)*UnitMh2UnitBl
  end subroutine SC_get_coord_for_sp_and_transform
  !============================================================================
  subroutine couple_sc_sp(DataInputTime)
    ! use CON_global_message_pass
    real,intent(in)::DataInputTime
    integer:: nLength
    !--------------------------------------------------------------------------
    if(.not.RouterScBl%IsProc)RETURN
    if(DoTest.and.is_proc0(SC_))&
         write(*,'(a,es12.5)')NameSub//': couple to SC,time=', DataInputTime
    if(is_proc(BL_))then
       call BL_put_coupling_param(Source_ = SC_, TimeIn  = DataInputTime)
       if(IsSource4BL_C(IH_))then
          call BL_get_bounds(rMinIn = RScMin, rMaxIn = RScMax,&
               rBufferUpIn = RIhMin)
       else
          call BL_get_bounds(rMinIn = RScMin, rMaxIn = RScMax)
       end if
    end if
    call  transform_matrix_and_coef(DataInputTime, SC_)
    call SC_synchronize_refinement(RouterScBl%iProc0Source,&
         RouterScBl%iCommUnion)
    if(DataInputTime > DataInputTimeLast)then
       ! The advected particles from SC are taken
       if(.not.DoInit) call BL_put_lines_from_sc
       ! Some particles may be lost, recover line intergity
       if(is_proc(BL_))call BL_adjust_lines(SC_)
    end if
    ! Set router SC=> BL to  receive MHD data
    call set_router(&
         GridSource             = SC_Grid                  ,&
         GridTarget             = BL_LocalGrid             ,&
         Router                 = RouterScBl               ,&
         n_interface_point_in_block = BL_n_particle        ,&
         interface_point_coords = BL_interface_point_coords,&
         is_interface_block     = BL_is_interface_block    ,&
         mapping                = mapping_sp_to_sc         ,&
         interpolate            = interpolation_amr_gc)
    call couple_comp(RouterScBl                            ,&
         nVar                   = nVarBuffer               ,&
         fill_buffer            = SC_get_for_mh            ,&
         apply_buffer           = BL_put_from_mh)
    !
    ! The MHD data within the heliocentric distances RScMin<R<RScMax
    ! By the way get particle coords from the router buffer to SC
    if(is_proc(SC_))then
       nLength = nlength_buffer_source(RouterScBl)
       call SC_put_particles(Xyz_DI =                       &
            RouterScBl%BufferSource_II(1:nDim, 1:nLength)  ,&
            iIndex_II = nint(RouterScBl%BufferSource_II(    &
            nDim+1:nDim+nAux, 1:nLength)) )
    end if
  end subroutine couple_sc_sp                !^CMP END SC
  !============================================================================
  subroutine couple_ih_sp_init(tNow)
    real,   intent(in) :: tNow
    logical            :: IsReady
    integer, parameter :: Upper_ = 1
    integer            :: nLength
    ! Check whether IH is ready for coupling with BL;
    !--------------------------------------------------------------------------
    call IH_check_ready_for_sp(IsReady)
    if(.not.IsReady) call CON_stop(&
         "IH component not ready for "//NameSub//" correct PARAM.in")
    call set_couple_var_info(IH_, BL_)
    ! Initialize coupler from IH (source )to BL (target)
    ! Data will be copied from IH to BL, however, the points
    ! in which the data shuld be provided will be sent from BL to
    ! IH, at the stage of the router construction. To further benefit
    ! from this opportunity, two BL grid indexes are also sent from
    ! BL to IH, at this stage, therefore nMappedPointIndex=nAux = 2.
    call init_router(                             &
         GridSource           = IH_Grid          ,&
         GridTarget           = BL_Grid          ,&
         nMappedPointIndex    = nAux             ,&
         Router               = RouterIhBl       )
    ! Set local GD on the Particle_I structure
    call init_router(IH_LineGrid, BL_Grid, RouterLineIhBl)
    if(.not.RouterIhBl%IsProc)RETURN
    ! Router to send particles is initialized.
    if(.not.DoExtract)RETURN
    if(is_proc(BL_))then
       call BL_put_coupling_param(Source_ = IH_, TimeIn = tNow)
       if(IsSource4BL_C(SC_))then              !^CMP IF SC BEGIN
          call BL_get_bounds(RScMax, RIhMax)
       else                                    !^CMP END SC
          call BL_get_bounds(RIhMin, RIhMax)
       end if                                  !^CMP IF SC
    end if
    call IH_synchronize_refinement(RouterIhBl%iProc0Source,&
         RouterIhBl%iCommUnion)
    call  transform_matrix_and_coef(tNow, IH_)
    ! Put in place the origin points of the MF lines.
    call set_router(                                            &
         GridSource                 = IH_Grid                  ,&
         GridTarget                 = BL_LocalGrid             ,&
         Router                     = RouterIhBl               ,&
         n_interface_point_in_block = BL_n_particle            ,&
         interface_point_coords     = BL_interface_point_coords,&
         is_interface_block         = BL_is_interface_block    ,&
         mapping                    = mapping_sp_to_ih         ,&
         interpolate                = interpolation_amr_gc)
    ! Now in IH are 1 point above RScMax per each line
    if(is_proc(IH_))then
       nLength = nlength_buffer_source(RouterIhBl)
       call IH_extract_line(Xyz_DI  = RouterIhBl%               &
            BufferSource_II(1:nDim,1:nLength)                  ,&
            iTraceMode              = Upper_                   ,&
            iIndex_II               = nint(RouterIhBl%          &
            BufferSource_II(nDim+1:nDim+nAux,1:nLength))       ,&
            RSoftBoundary           = RIhMax*UnitBl2UnitMh)
       ! Now in IH are the parts of lines from RScMax + 1 point above
       ! to RIhMax + 1 point above
    end if
    call BL_put_lines_from_ih
  end subroutine couple_ih_sp_init
  !============================================================================
  subroutine mapping_sp_to_ih(nDimIn, XyzBl_D, nDimOut, CoordIh_D, &
       IsInterfacePoint)
    integer, intent(in) :: nDimIn, nDimOut
    real,    intent(in) :: XyzBl_D(nDimIn)
    real,    intent(out):: CoordIh_D(nDimOut)
    logical, intent(out):: IsInterfacePoint
    real                :: XyzIh_D(nDim)
    !--------------------------------------------------------------------------
    IsInterfacePoint = .true.
    XyzIh_D = matmul(XyzBl_D, BlMh_DD)*UnitBl2UnitMh
    call IH_xyz_to_coord(XyzIh_D, CoordIh_D)
  end subroutine mapping_sp_to_ih
  !============================================================================
  subroutine BL_put_lines_from_ih
    !--------------------------------------------------------------------------
    call construct_router_from_source(&
         GridSource                 = IH_LocalLineGrid    ,&
         GridTarget                 = BL_Grid             ,&
         Router                     = RouterLineIhBl      ,&
         n_interface_point_in_block = IH_n_particle       ,&
         interface_point_coords     = i_particle          ,&
         mapping                    = mapping_line_ih_to_sp)
    call couple_comp(RouterLineIhBl                       ,&
         nVar = 3                                         ,&
         fill_buffer = IH_get_coord_for_sp_and_transform  ,&
         apply_buffer= BL_put_line)
  end subroutine BL_put_lines_from_ih
  !============================================================================
  subroutine mapping_line_ih_to_sp(nDimIn, XyzIn_D, nDimOut, CoordOut_D, &
       IsInterfacePoint)
    integer, intent(in) :: nDimIn, nDimOut
    real,    intent(in) :: XyzIn_D(nDimIn)
    real,    intent(out):: CoordOut_D(nDimOut)
    logical, intent(out):: IsInterfacePoint
    integer:: iIndex_I(2), iParticle
    !--------------------------------------------------------------------------
    ! Convert back the above converted particle number
    IsInterfacePoint = .true.; iParticle = nint(XyzIn_D(1))
    call IH_get_particle_indexes(iParticle, iIndex_I)
    ! Convert particle indexes to BL gen coords, so that the nearest_cell
    ! program in the router coud figure out where it should be put in BL
    CoordOut_D = coord_grid_d(BL_Grid, iIndex_I(1), [iIndex_I(2),1,1])
  end subroutine mapping_line_ih_to_sp
  !============================================================================
  subroutine IH_get_coord_for_sp_and_transform(&
       nPartial, iGetStart, Get, w, XyzBl_D, nDimOut)
    integer,            intent(in):: nPartial, iGetStart, nDimOut
    type(IndexPtrType), intent(in):: Get
    type(WeightPtrType),intent(in):: w
    real,              intent(out):: XyzBl_D(nDimOut)
    real :: XyzIh_D(nDim)
    !--------------------------------------------------------------------------
    call IH_get_particle_coords(Get%iCB_II(1, iGetStart), XyzIh_D)
    ! perform transformation before returning
    XyzBl_D = matmul(BlMh_DD, XyzIh_D)*UnitMh2UnitBl
  end subroutine IH_get_coord_for_sp_and_transform
  !============================================================================
  subroutine couple_ih_sp(DataInputTime)
    real,   intent(in) :: DataInputTime
    integer            :: nLength
    !--------------------------------------------------------------------------
    if(.not.RouterIhBl%IsProc)then
       ! If OH is used, DoInit is reset to .false. there
       if(DoInit)DoInit = IsSource4BL_C(OH_)
       if(.not.IsSource4BL_C(OH_))DataInputTimeLast = DataInputTime
       RETURN
    end if
    if(DoTest.and.is_proc0(IH_))&
         write(*,'(a,es12.5)')NameSub//': couple to IH,time=', DataInputTime
    if(is_proc(BL_))then
       call BL_put_coupling_param(Source_ = IH_, TimeIn  = DataInputTime)

       call BL_get_bounds(rMinIn = RIhMin, rMaxIn = RIhMax, &
            rBufferLoIn = RScMax)
    end if
    call  transform_matrix_and_coef(DataInputTime, IH_)
    call IH_synchronize_refinement(RouterIhBl%iProc0Source,&
         RouterIhBl%iCommUnion)
    if(DataInputTime > DataInputTimeLast)then
       if(.not.DoInit) call BL_put_lines_from_ih
       ! Some points may be lost lost in IH
       if(is_proc(BL_))call BL_adjust_lines(IH_)
    end if
    ! In points at RIhMin < R < RIhMax get the MHD data
    ! Set router IH=> BL to  receive MHD data
    call set_router(                                            &
         GridSource                 = IH_Grid                  ,&
         GridTarget                 = BL_LocalGrid             ,&
         Router                     = RouterIHBl               ,&
         n_interface_point_in_block = BL_n_particle            ,&
         interface_point_coords     = BL_interface_point_coords,&
         is_interface_block         = BL_is_interface_block    ,&
         mapping                    = mapping_sp_to_ih         ,&
         interpolate                = interpolation_amr_gc)
    call couple_comp(                 RouterIhBl               ,&
         nVar                       = nVarBuffer               ,&
         fill_buffer                = IH_get_for_mh,            &
         apply_buffer               = BL_put_from_mh)
    ! By the way get coords for particles in IH from the router buffer
    if(is_proc(IH_))then
       nLength = nlength_buffer_source(RouterIhBl)
       if(IsSource4BL_C(SC_))call sort_out_sc_particles     !^CMP IF SC
       call IH_put_particles(Xyz_DI =  RouterIhBl%              &
            BufferSource_II(1:nDim, 1:nLength)                 ,&
            iIndex_II               = nint(RouterIhBl%          &
            BufferSource_II(nDim+1:nDim+nAux, 1:nLength)))
    end if
    ! If OH is used, DoInit and DataInputTimeLast are reset in OH
    if(DoInit)DoInit = IsSource4BL_C(OH_)
    if(.not.IsSource4BL_C(OH_))DataInputTimeLast = DataInputTime
  contains                                       !^CMP IF SC BEGIN
    !==========================================================================
    subroutine sort_out_sc_particles
      ! Sort out particles with R < RScMax to be advected by the SC
      integer:: iParticle, iParticleNew
      real   :: Xyz_D(nDim), Coord_D(nDim)
      !------------------------------------------------------------------------
      iParticleNew = 0
      do iParticle = 1, nLength
         Coord_D = RouterIhBl%BufferSource_II(1:nDim, iParticle)
         call IH_coord_to_xyz(Coord_D, Xyz_D)
         if(norm2(Xyz_D) < RScMax*UnitBl2UnitMh)CYCLE
         iParticleNew  = iParticleNew + 1
         RouterIhBl%BufferSource_II(:, iParticleNew) = &
              RouterIhBl%BufferSource_II(:, iParticle)
      end do
      nLength = iParticleNew
    end subroutine sort_out_sc_particles                 !^CMP END SC
    !==========================================================================
  end subroutine couple_ih_sp
  !============================================================================
  !^CMP IF OH BEGIN
  subroutine couple_oh_sp_init(tNow)
    real,   intent(in) :: tNow
    logical            :: IsReady
    integer, parameter :: Upper_ = 1
    integer            :: nLength
    ! Check whether OH is ready for coupling with BL;
    !--------------------------------------------------------------------------
    call OH_check_ready_for_sp(IsReady)
    if(.not.IsReady) call CON_stop(&
         "OH component not ready for "//NameSub//" correct PARAM.in")
    call set_couple_var_info(OH_, BL_)
    ! Initialize coupler from OH (source )to BL (target)
    ! Data will be copied from OH to BL, however, the points
    ! in which the data shuld be provided will be sent from BL to
    ! IH, at the stage of the router construction. To further benefit
    ! from this opportunity, two BL grid indexes are also sent from
    ! BL to OH, at this stage, therefore nMappedPointIndex = 2.
    call init_router(                             &
         GridSource           = OH_Grid          ,&
         GridTarget           = BL_Grid          ,&
         nMappedPointIndex    = nAux             ,&
         Router               = RouterOhBl        )
    ! Set local GD on the Particle_I structure
    call init_router(OH_LineGrid, BL_Grid, RouterLineOhBl)
    if(.not.RouterOhBl%IsProc)RETURN
    ! Router to send particles is initialized.
    if(.not.DoExtract)RETURN
    if(is_proc(BL_))then
       call BL_put_coupling_param(Source_ = OH_, TimeIn = tNow)
       if(IsSource4BL_C(IH_))then
          call BL_get_bounds(RIhMax, ROhMax)
       else
          call BL_get_bounds(ROhMin, ROhMax)
       end if
    end if
    call OH_synchronize_refinement(RouterOhBl%iProc0Source,&
         RouterOhBl%iCommUnion)
    call  transform_matrix_and_coef(tNow, OH_)
    ! Put in place the origin points of the MF lines.
    call set_router(                                            &
         GridSource                 = OH_Grid                  ,&
         GridTarget                 = BL_LocalGrid             ,&
         Router                     = RouterOhBl               ,&
         n_interface_point_in_block = BL_n_particle            ,&
         interface_point_coords     = BL_interface_point_coords,&
         is_interface_block         = BL_is_interface_block    ,&
         mapping                    = mapping_sp_to_OH         ,&
         interpolate                = interpolation_amr_gc)
    ! Now in OH-semi-router is 1 point above RRIhMax per each line
    if(is_proc(OH_))then
       nLength = nlength_buffer_source(RouterOhBl)
       call OH_extract_line(Xyz_DI  = RouterOhBl%               &
            BufferSource_II(1:nDim,1:nLength)                  ,&
            iTraceMode              = Upper_                   ,&
            iIndex_II               = nint(RouterOhBl%          &
            BufferSource_II(nDim+1:nDim+nAux,1:nLength))       ,&
            RSoftBoundary           = ROhMax*UnitBl2UnitMh)
       ! Now in OH are the parts of lines from RIHMax + 1 point above
       ! to ROhMax + 1 Point above
    end if
    call BL_put_lines_from_oh
  end subroutine couple_oh_sp_init
  !============================================================================
  subroutine mapping_sp_to_oh(nDimIn, XyzBl_D, nDimOut, CoordOh_D, &
       IsInterfacePoint)
    integer, intent(in) :: nDimIn, nDimOut
    real,    intent(in) :: XyzBl_D(nDimIn)
    real,    intent(out):: CoordOh_D(nDimOut)
    logical, intent(out):: IsInterfacePoint
    real                :: XyzOh_D(nDim)
    !--------------------------------------------------------------------------
    IsInterfacePoint = .true.
    XyzOh_D = matmul(XyzBl_D, BlMh_DD)*UnitBl2UnitMh
    call OH_xyz_to_coord(XyzOh_D, CoordOh_D)
  end subroutine mapping_sp_to_oh
  !============================================================================
  subroutine BL_put_lines_from_oh
    !--------------------------------------------------------------------------
    call construct_router_from_source(&
         GridSource                 = OH_LocalLineGrid    ,&
         GridTarget                 = BL_Grid             ,&
         Router                     = RouterLineOhBl      ,&
         n_interface_point_in_block = OH_n_particle       ,&
         interface_point_coords     = i_particle          ,&
         mapping                    = mapping_line_oh_to_sp)
    call couple_comp(RouterLineOhBl                       ,&
         nVar = 3                                         ,&
         fill_buffer = OH_get_coord_for_sp_and_transform  ,&
         apply_buffer= BL_put_line)
  end subroutine BL_put_lines_from_oh
  !============================================================================
  subroutine mapping_line_oh_to_sp(nDimIn, XyzIn_D, &
       nDimOut, CoordOut_D, IsInterfacePoint)
    integer, intent(in) :: nDimIn, nDimOut
    real,    intent(in) :: XyzIn_D(nDimIn)
    real,    intent(out):: CoordOut_D(nDimOut)
    logical, intent(out):: IsInterfacePoint
    integer:: iIndex_I(2), iParticle
    !--------------------------------------------------------------------------
    ! Convert back the above converted particle number
    IsInterfacePoint = .true.; iParticle = nint(XyzIn_D(1))
    call OH_get_particle_indexes(iParticle, iIndex_I)
    ! Convert particle indexes to BL gen coords, so that the nearest_cell
    ! program in the router coud figure out where it should be put in BL
    CoordOut_D = coord_grid_d(BL_Grid, iIndex_I(1), [iIndex_I(2),1,1])
  end subroutine mapping_line_oh_to_sp
  !============================================================================
  subroutine OH_get_coord_for_sp_and_transform(&
       nPartial, iGetStart, Get, w, XyzBl_D, nDimOut)
    integer,            intent(in):: nPartial, iGetStart, nDimOut
    type(IndexPtrType), intent(in):: Get
    type(WeightPtrType),intent(in):: w
    real,              intent(out):: XyzBl_D(nDimOut)
    real :: XyzOh_D(nDim)
    !--------------------------------------------------------------------------
    call OH_get_particle_coords(Get%iCB_II(1, iGetStart), XyzOh_D)
    ! perform transformation before returning
    XyzBl_D = matmul(BlMh_DD, XyzOh_D)*UnitMh2UnitBl
  end subroutine OH_get_coord_for_sp_and_transform
  !============================================================================
  subroutine couple_oh_sp(DataInputTime)
    real,   intent(in) :: DataInputTime
    integer            :: nLength
    !--------------------------------------------------------------------------
    if(.not.RouterOhBl%IsProc)then
       DoInit = .false.
       DataInputTimeLast = DataInputTime
       RETURN
    end if
    if(DoTest.and.is_proc0(OH_))&
         write(*,'(a,es12.5)')NameSub//': couple to OH,time=', DataInputTime
    if(is_proc(BL_))then
       call BL_put_coupling_param(Source_ = OH_, TimeIn = DataInputTime)
       call BL_get_bounds(rMinIn = ROhMin, rMaxIn = ROhMax, &
            rBufferLoIn = RIhMax)
    end if
    call  transform_matrix_and_coef(DataInputTime, OH_)
    call OH_synchronize_refinement(RouterOhBl%iProc0Source,&
         RouterOhBl%iCommUnion)
    if(DataInputTime > DataInputTimeLast)then
       if(.not.DoInit) call BL_put_lines_from_oh
       ! Now the full magnetic line is available, including probably
       ! some points lost in IH
       if(is_proc(BL_))call BL_adjust_lines(OH_)
    end if
    ! In points at RIhMin < R < RIhMax get the MHD data
    ! Set router IH=> BL to  receive MHD data
    call set_router(                                            &
         GridSource                 = OH_Grid                  ,&
         GridTarget                 = BL_LocalGrid             ,&
         Router                     = RouterOhBl               ,&
         n_interface_point_in_block = BL_n_particle            ,&
         interface_point_coords     = BL_interface_point_coords,&
         is_interface_block         = BL_is_interface_block    ,&
         mapping                    = mapping_sp_to_oh         ,&
         interpolate                = interpolation_amr_gc)
    call couple_comp(                 RouterOhBl               ,&
         nVar                       = nVarBuffer               ,&
         fill_buffer                = OH_get_for_mh            ,&
         apply_buffer               = BL_put_from_mh)
    ! By the way get coords for particles in OH from the router buffer
    if(is_proc(OH_))then
       nLength = nlength_buffer_source(RouterOhBl)
       if(IsSource4Bl_C(IH_))call sort_out_ih_particles
       call OH_put_particles(Xyz_DI =  RouterOhBl%              &
            BufferSource_II(1:nDim, 1:nLength)                 ,&
            iIndex_II               = nint(RouterOhBl%          &
            BufferSource_II(nDim+1:nDim+nAux, 1:nLength)))
    end if
    DoInit = .false.
    DataInputTimeLast = DataInputTime
  contains
    !==========================================================================
    subroutine sort_out_ih_particles
      ! Sort out particles to be advected by the IH, which have
      ! R in BL coordinates less than RIhMax**2
      integer:: iParticle, iParticleNew
      real   :: Xyz_D(3), Coord_D(3)
      !------------------------------------------------------------------------
      iParticleNew = 0
      do iParticle = 1, nLength
         Coord_D = RouterOhBl%BufferSource_II(1:nDim, iParticle)
         call OH_coord_to_xyz(Coord_D, Xyz_D)
         if(norm2(Xyz_D) < RIhMax*UnitBl2UnitMh)CYCLE
         iParticleNew  = iParticleNew + 1
         RouterOhBl%BufferSource_II(:, iParticleNew) = &
              RouterOhBl%BufferSource_II(:, iParticle)
      end do
      nLength = iParticleNew
    end subroutine sort_out_ih_particles
    !==========================================================================
  end subroutine couple_oh_sp         !^CMP END OH
  !============================================================================
end module CON_couple_mh_sp
!==============================================================================
