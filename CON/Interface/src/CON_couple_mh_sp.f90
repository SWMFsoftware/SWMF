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
       IH_LineDD      ! 1D "grid": Coord = real(iParticle + iProc*nParticleMax)

  !^CMP IF SC BEGIN
  use SC_wrapper, ONLY: SC_check_ready_for_sp, SC_synchronize_refinement,    &
       SC_extract_line, SC_put_particles, SC_n_particle,                     &
       SC_get_particle_coords, SC_get_particle_indexes, SC_xyz_to_coord,     &
       SC_get_for_mh, SC_LineDD
  !^CMP END SC
  use SP_wrapper, ONLY: &
       SP_check_ready_for_mh    ,&  ! If returns .false., extract the mf lines
       SP_get_bounds_comp       ,&  ! Provides RScMin/Max and/or RIhMin/Max
       SP_put_coupling_param    ,&  ! Set time and interaface bounds
       SP_adjust_lines              ! Process if needed the updated mf lines
  use CON_mflampa, ONLY: &
       MF_n_particle            ,&  ! Number of "points" in a given line in SP
       MF_put_from_mh           ,&  ! Put MHD info from SC or IH to SP
       MF_interface_point_coords,&  ! Check if the point is within interface
       MF_put_line                  ! Put particle Xyz from SC/IH to SP
       !get_bounds                   ! Provides RScMin/Max and/or RIhMin/Max
  implicit none

  private ! Except
  public::couple_mh_sp_init
  public::couple_ih_sp
  public::couple_sc_sp              !^CMP IF SC

  type(GridType)     , save::SP_Grid           ! Target (Particle coords)
  type(LocalGridType), save::SP_LocalGrid      ! Target (MHD data)
  type(GridType)     , save::IH_Grid           ! Source (MHD data)
  type(RouterType),save,private::RouterIhSp        ! IH (MHD data) => SP
  type(GridType)     , save::IH_LineGrid       ! Misc
  type(LocalGridType), save::IH_LocalLineGrid  ! Source (MHD data)
  type(RouterType),save,private::RouterLineIhSp    ! IH (Particle coords)=>SP
  !^CMP IF SC BEGIN
  type(GridType)     , save::SC_Grid           ! Source (MHD data)
  type(RouterType),save,private::RouterScSp        ! SC MHD data => SP
  type(GridType)     , save::SC_LineGrid       ! Misc
  type(LocalGridType), save::SC_LocalLineGrid  ! Source (Particle Coords)
  type(RouterType),save,private::RouterLineScSp    ! SC (Particle coords)=>SP
  !^CMP END SC

  ! Three-dimensional grids in SC and IH,
  ! three-component coordinate vector
  integer, parameter :: nDim = 3

  ! We send: (1) the global index of the magnetic field line
  ! and (2) the particle number along the magnetic field line
  integer, parameter :: nAux = 2

  ! Misc
  integer :: nLength, iError

  ! Transformation matrices and the unit of length ratios
  real    :: ScToSp_DD(3,3), UnitSp2UnitSc, UnitSc2UnitSp !^CMP IF SC
  real    :: IhToSp_DD(3,3), UnitSp2UnitIh, UnitIh2UnitSp

  logical::DoTest,DoTestMe
  character(LEN=*),parameter::NameSub='couple_mh_sp'
  real :: tNow

  integer, parameter :: Lower_ = 0
  !^CMP IF SC BEGIN
  ! Short ID of MHD components. If there is only one component,
  ! it has an identifier Lower_, otherwise, SC has identifier Lower_
  ! while IH is Upper
  integer, parameter :: Upper_ = 1
  ! solar corona and lower and upper boundaries
  real:: RScMin, RScMax             !^CMP END SC
  ! inner heliosphere lower and upper boundaries
  real:: RIhMin, RIhMax
  logical ::DoInit=.true.
contains
  !============================================================================
  subroutine couple_mh_sp_init
    use CON_physics, ONLY: get_time
    use ModConst
    ! whether a component is ready for coupling
    logical :: IsReady
    ! whether need to extract new line (equal to .not.DoRestart in SP
    logical :: DoExtract
    !--------------------------------------------------------------------------
    if(.not.DoInit)RETURN
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    ! The initialization can be done only once

    call get_time(tSimulationOut=tNow)

    ! determine whether SP is ready for coupling with MH
    call SP_check_ready_for_mh(IsReady)
    ! if not ready => need to extract new field lines
    DoExtract = .not.IsReady

    if(use_comp(SC_))call couple_sc_sp_init  !^CMP IF SC
    if(use_comp(IH_))call couple_ih_sp_init

    ! After couple_mh_sp_init
    ! (a) set of low and upper boundaries in the coupler
    ! (b) set of the Lagrangian point coordinates in SP
    ! (c) logical DoInit is .true. to prevent the coordinate
    ! info on SP to be overritten with that from SC or IH
  contains
    !==========================================================================
    subroutine couple_sc_sp_init             !^CMP IF SC BEGIN
      ! get the value of SC boundaries as set in SP
      !------------------------------------------------------------------------
      call SP_get_bounds_comp(Lower_, RScMin, RScMax)
      ! Check whether SC is ready for coupling with SP;
      call SC_check_ready_for_sp(IsReady)
      if(.not.IsReady) call CON_stop(&
           "SC component not ready for "//NameSub//" correct PARAM.in")
      call set_couple_var_info(SC_, SP_)
      UnitSp2UnitSc = Grid_C(SP_)%UnitX/Grid_C(SC_)%UnitX
      UnitSc2UnitSp = 1/UnitSp2UnitSc

      ! Initialize coupler from SC (source )to SP (target)
      ! Data will be copied from SC to SP, however, the points
      ! in which the data shuld be provided will be sent from SP to
      ! SC, at the stage of the router construction. To further benefit
      ! from this opportunity, two SP grid indexes are also sent from
      ! SP to SC,at this stage, therefore nMappedPointIndex = 2.
      call init_coupler(iCompSource = SC_          ,&
           iCompTarget       = SP_                 ,&
           GridSource        = SC_Grid             ,&
           GridTarget        = SP_Grid             ,&
           nMappedPointIndex = nAux                ,&
           Router            = RouterScSp          ,&
           LocalGridTarget   = SP_LocalGrid)

      ! Set local GD on the Particle_I structure
      call set_standard_grid_descriptor(SC_LineDD,Grid=SC_LineGrid)
      call init_router(SC_LineGrid, SP_Grid, RouterLineScSp)
      if(.not.RouterScSp%IsProc)RETURN
      if(is_proc(SC_))call set_local_gd(iProc = i_proc(),   &
           Grid = SC_LineGrid, LocalGrid = SC_LocalLineGrid)
      ! Router to send particles is initialized.
      if(.not.DoExtract)RETURN
      !
      ! If DoExtract, then
      !
      if(is_proc(SP_))&
           call SP_put_coupling_param(Lower_, RScMin, RScMax, tNow)
      call SC_synchronize_refinement(RouterScSp%iProc0Source,&
           RouterScSp%iCommUnion)
      ScToSp_DD=transform_matrix(tNow,&
           Grid_C(SC_)%TypeCoord, Grid_C(SP_)%TypeCoord)
      ! use router intended for sending MHD data  (SC => SP)
      ! to send the particle coordinates backward (SP=>SC).

      ! Put in place the origin points of the MF lines
      call set_router(&
           GridSource            = SC_Grid                  ,&
           GridTarget            = SP_LocalGrid             ,&
           Router                = RouterScSp               ,&
           n_interface_point_in_block = MF_n_particle       ,&
           interface_point_coords= MF_interface_point_coords,&
           mapping               = mapping_sp_to_sc         ,&
           interpolate           = interpolation_amr_gc)
      if(is_proc(SC_))then
         nLength = nlength_buffer_source(RouterScSp)
         call SC_extract_line(&
              Xyz_DI     = RouterScSp%BufferSource_II(1:nDim,&
              1:nLength)                                    ,&
              iTraceMode = Lower_                           ,&
              iIndex_II  = nint(RouterScSp% BufferSource_II( &
              nDim+1:nDim+nAux,1:nLength))                  ,&
              RSoftBoundary = RScMax*UnitSp2UnitSc)
      end if
      call SP_get_lines_from_sc
    end subroutine couple_sc_sp_init                 !^CMP END SC
    !==========================================================================
    subroutine couple_ih_sp_init
      ! get the value of IH boundaries as set in SP
      !------------------------------------------------------------------------
      if(use_comp(SC_))then  !^CMP IF SC BEGIN
         call SP_get_bounds_comp(Upper_, RIhMin, RIhMax)
      else                   !^CMP END SC
         call SP_get_bounds_comp(Lower_, RIhMin, RIhMax)
      end if                 !^CMP IF SC

      ! Check whether IH is ready for coupling with SP;
      call IH_check_ready_for_sp(IsReady)
      if(.not.IsReady) call CON_stop(&
           "IH component not ready for "//NameSub//" correct PARAM.in")
      call set_couple_var_info(IH_, SP_)
      UnitSp2UnitIh = Grid_C(SP_)%UnitX/Grid_C(IH_)%UnitX
      UnitIh2UnitSp = 1/UnitSp2UnitIh

      ! Initialize coupler from IH (source )to SP (target)
      ! Data will be copied from IH to SP, however, the points
      ! in which the data shuld be provided will be sent from SP to
      ! IH, at the stage of the router construction. To further benefit
      ! from this opportunity, two SP grid indexes are also sent from
      ! SP to IH, at this stage, therefore nMappedPointIndex = 2.
      if(use_comp(SC_))then          !^CMP IF SC BEGIN
         ! Do not reset already installed SP_LocalGrid
         call init_coupler(iCompSource = IH_          ,&
              iCompTarget          = SP_              ,&
              GridSource           = IH_Grid          ,&
              GridTarget           = SP_Grid          ,&
              nMappedPointIndex    = nAux             ,&
              Router               = RouterIhSp)
      else                           !^CMP END SC
         call init_coupler(iCompSource = IH_          ,&
              iCompTarget          = SP_              ,&
              GridSource           = IH_Grid          ,&
              GridTarget           = SP_Grid          ,&
              nMappedPointIndex    = nAux             ,&
              Router               = RouterIhSp       ,&
              LocalGridTarget      = SP_LocalGrid)
      end if                         !^CMP IF SC
      ! Set local GD on the Particle_I structure
      call set_standard_grid_descriptor(IH_LineDD,Grid=&
           IH_LineGrid)
      call init_router(IH_LineGrid, SP_Grid, RouterLineIhSp)
      if(.not.RouterIhSp%IsProc)RETURN
      if(is_proc(IH_))call set_local_gd(iProc = i_proc(),&
           Grid      = IH_LineGrid                      ,&
           LocalGrid = IH_LocalLineGrid)
      ! Router to send particles is initialized.
      if(.not.DoExtract)RETURN
      if(is_proc(SP_))then
         if(use_comp(SC_))then             !^CMP IF SC BEGIN
            call SP_put_coupling_param(&
                 iModelIn    = Upper_, &
                 rMinIn      = RScMax, &
                 rMaxIn      = RIhMax, &
                 TimeIn      = tNow, &
                 rBufferUpIn = RIhMin)
         else                              !^CMP END SC
            call SP_put_coupling_param(&
                 iModelIn    = Lower_, &
                 rMinIn      = RIhMin, &
                 rMaxIn      = RIhMax,&
                 TimeIn      = tNow,&
                 rBufferLoIn = RScMax)
         end if                            !^CMP IF SC
      end if
      call IH_synchronize_refinement(RouterIhSp%iProc0Source,&
           RouterIhSp%iCommUnion)
      IhToSp_DD=transform_matrix(tNow,&
           Grid_C(IH_)%TypeCoord, Grid_C(SP_)%TypeCoord)
      ! Put in place the origin points of the MF lines.
      call set_router(                                            &
           GridSource                 = IH_Grid                  ,&
           GridTarget                 = SP_LocalGrid             ,&
           Router                     = RouterIHSp               ,&
           n_interface_point_in_block = MF_n_particle            ,&
           interface_point_coords     = MF_interface_point_coords,&
           mapping                    = mapping_sp_to_IH         ,&
           interpolate                = interpolation_amr_gc)
      if(is_proc(IH_))then
         nLength = nlength_buffer_source(RouterIhSp)
         call IH_extract_line(Xyz_DI  = RouterIhSp%               &
              BufferSource_II(1:nDim,1:nLength)                  ,&
              iTraceMode              = Upper_                   ,&
              iIndex_II               = nint(RouterIhSp%          &
              BufferSource_II(nDim+1:nDim+nAux,1:nLength))       ,&
              RSoftBoundary           = RIhMax*UnitSp2UnitIh)
      end if
      call SP_get_lines_from_ih
    end subroutine couple_ih_sp_init
    !==========================================================================
  end subroutine couple_mh_sp_init
  !============================================================================
  !^CMP IF SC BEGIN
  subroutine couple_sc_sp(DataInputTime)
    use CON_global_message_pass

    real,intent(in)::DataInputTime
    integer:: nLength
    !--------------------------------------------------------------------------
    if(.not.RouterScSp%IsProc)RETURN

    tNow=DataInputTime
    if(is_proc(SP_))call SP_put_coupling_param(&
         iModelIn    = Lower_, &
         rMinIn      = RScMin,&
         rMaxIn      = RScMax,&
         TimeIn      = DataInputTime,&
         rBufferUpIn = RIhMin)
    ScToSp_DD=transform_matrix(tNow,&
         Grid_C(SC_)%TypeCoord, Grid_C(SP_)%TypeCoord)
    call SC_synchronize_refinement(RouterScSp%iProc0Source,&
         RouterScSp%iCommUnion)
    if(.not.DoInit) call SP_get_lines_from_sc
    if(is_proc(SP_))call SP_adjust_lines(DoInit)
    ! Set router SC=> SP to  receive MHD data
    call set_router(&
         GridSource             = SC_Grid                  ,&
         GridTarget             = SP_LocalGrid             ,&
         Router                 = RouterScSp               ,&
         n_interface_point_in_block = MF_n_particle        ,&
         interface_point_coords = MF_interface_point_coords,&
         mapping                = mapping_sp_to_sc         ,&
         interpolate            = interpolation_amr_gc)
    call couple_comp(RouterScSp                            ,&
         nVar                   = nVarBuffer               ,&
         fill_buffer = SC_get_for_sp_and_transform         ,&
         apply_buffer           = MF_put_from_mh)
    ! By the way get particle coords from the router buffer
    if(is_proc(SC_))then
       nLength = nlength_buffer_source(RouterScSp)
       call SC_put_particles(Xyz_DI =                       &
            RouterScSp%BufferSource_II(1:nDim, 1:nLength)  ,&
            iIndex_II = nint(RouterScSp%BufferSource_II(    &
            nDim+1:nDim+nAux, 1:nLength)))
    end if
  end subroutine couple_sc_sp
  !============================================================================
  subroutine SP_get_lines_from_sc
    !--------------------------------------------------------------------------
    call construct_router_from_source(                      &
         GridSource                 = SC_LocalLineGrid     ,&
         GridTarget                 = SP_Grid              ,&
         Router                     = RouterLineScSp       ,&
         n_interface_point_in_block = SC_n_particle        ,&
         interface_point_coords     = i_particle           ,&
         mapping                    = mapping_line_sc_to_sp)
    call couple_comp(                 RouterLineScSp       ,&
         nVar                       = nDim                 ,&
         fill_buffer = SC_get_coord_for_sp_and_transform   ,&
         apply_buffer               = MF_put_line)
  end subroutine sp_get_lines_from_sc
  !============================================================================
  subroutine mapping_sp_to_sc(nDimIn, XyzIn_D, nDimOut,     &
       CoordOut_D, IsInterfacePoint)
    integer, intent(in) :: nDimIn, nDimOut
    real,    intent(in) :: XyzIn_D(nDimIn)
    real,    intent(out):: CoordOut_D(nDimOut)
    logical, intent(out):: IsInterfacePoint

    real                :: XyzTemp_D(nDim)
    !--------------------------------------------------------------------------
    IsInterfacePoint = .true.
    XyzTemp_D = matmul(XyzIn_D, ScToSp_DD)*UnitSp2UnitSc
    call SC_xyz_to_coord(XyzTemp_D, CoordOut_D)
  end subroutine mapping_sp_to_sc
  !============================================================================
  subroutine mapping_line_sc_to_sp(nDimIn, XyzIn_D, nDimOut, &
       CoordOut_D, IsInterfacePoint)
    ! See comments in mapping_line_ih_to_sp
    integer, intent(in) :: nDimIn, nDimOut
    real,    intent(in) :: XyzIn_D(nDimIn)
    real,    intent(out):: CoordOut_D(nDimOut)
    logical, intent(out):: IsInterfacePoint

    integer:: iIndex_I(2), iParticle
    !--------------------------------------------------------------------------
    IsInterfacePoint = .true.; iParticle = nint(XyzIn_D(1))
    call SC_get_particle_indexes(iParticle, iIndex_I)
    CoordOut_D = coord_grid_d(SP_Grid, iIndex_I(1),[iIndex_I(2),1,1])
  end subroutine mapping_line_sc_to_sp
  !============================================================================
  subroutine SC_get_for_sp_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V

    integer:: iVarBx, iVarBz, iVarMx, iVarMz
    !--------------------------------------------------------------------------
    ! get buffer with variables
    call SC_get_for_mh(nPartial,iGetStart,Get,w,State_V,nVar)
    ! indices of variables
    iVarBx = iVar_V(BxCouple_);   iVarBz = iVar_V(BzCouple_)
    iVarMx = iVar_V(RhoUxCouple_);iVarMz = iVar_V(RhoUzCouple_)
    ! perform transformation before returning
    State_V(iVarBx:iVarBz) = matmul(ScToSp_DD,State_V(iVarBx:iVarBz))
    State_V(iVarMx:iVarMz) = matmul(ScToSp_DD,State_V(iVarMx:iVarMz))
  end subroutine SC_get_for_sp_and_transform
  !============================================================================
  subroutine SC_get_coord_for_sp_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V

    !--------------------------------------------------------------------------
    call SC_get_particle_coords(Get%iCB_II(1,iGetStart),State_V)
    ! coord transformation
    State_V = matmul(ScToSp_DD, State_V)*UnitSc2UnitSp
  end subroutine SC_get_coord_for_sp_and_transform
  !============================================================================
  !=============================================== !^CMP END SC
  subroutine couple_ih_sp(DataInputTime)
    real,intent(in)::DataInputTime
    !--------------------------------------------------------------------------
    if(.not.RouterIhSp%IsProc)RETURN
    tNow = DataInputTime
    if(is_proc(SP_))call SP_put_coupling_param(&
         iModelIn    = Upper_, &
         rMinIn      = RIhMin, &
         rMaxIn      = RIhMax, &
         TimeIn      = DataInputTime, &
         rBufferLoIn = RScMax)
    IhToSp_DD=transform_matrix(tNow,&
         Grid_C(IH_)%TypeCoord, Grid_C(SP_)%TypeCoord)
    call IH_synchronize_refinement(RouterIhSp%iProc0Source,&
         RouterIhSp%iCommUnion)
    if(.not.DoInit) call SP_get_lines_from_ih
    if(is_proc(SP_))call SP_adjust_lines(DoInit)

    ! Set router IH=> SP to  receive MHD data
    call set_router(                                            &
         GridSource                 = IH_Grid                  ,&
         GridTarget                 = SP_LocalGrid             ,&
         Router                     = RouterIHSp               ,&
         n_interface_point_in_block = MF_n_particle            ,&
         interface_point_coords     = MF_interface_point_coords,&
         mapping                    = mapping_sp_to_ih         ,&
         interpolate                = interpolation_amr_gc)
    call couple_comp(                 RouterIhSp               ,&
         nVar                       = nVarBuffer               ,&
         fill_buffer              = IH_get_for_sp_and_transform,&
         apply_buffer               = MF_put_from_mh)
    ! By the way get particle coords from the router buffer
    if(is_proc(IH_))then
       nLength = nlength_buffer_source(RouterIhSp)
       if(use_comp(SC_))call sort_out_sc_particles     !^CMP IF SC
       call IH_put_particles(Xyz_DI =  RouterIhSp%              &
            BufferSource_II(1:nDim, 1:nLength)                 ,&
            iIndex_II               = nint(RouterIhSp%          &
            BufferSource_II(nDim+1:nDim+nAux, 1:nLength)))
    end if
    DoInit=.false.
  contains                                       !^CMP IF SC BEGIN
    !==========================================================================
    subroutine sort_out_sc_particles
      ! Sort out particles to be advected by the SC, which have
      ! R2 in SP coordinates less than RScMax**2
      integer:: iParticle, iParticleNew
      real   :: Xyz_D(3), Coord_D(3), RScMax2

      !------------------------------------------------------------------------
      iParticleNew = 0; RScMax2 = (RScMax*UnitSp2UnitIh)**2
      do iParticle = 1, nLength
         Coord_D = RouterIhSp%BufferSource_II(1:nDim, iParticle)
         call IH_coord_to_xyz(Coord_D, Xyz_D)
         if(sum(Xyz_D**2) < RScMax2)CYCLE
         iParticleNew  = iParticleNew + 1
         RouterIhSp%BufferSource_II(:, iParticleNew) = &
              RouterIhSp%BufferSource_II(:, iParticle)
      end do
      nLength = iParticleNew
    end subroutine sort_out_sc_particles                 !^CMP END SC
    !==========================================================================
  end subroutine couple_ih_sp
  !============================================================================
  subroutine SP_get_lines_from_ih
    !--------------------------------------------------------------------------
    call construct_router_from_source(&
         GridSource                 = IH_LocalLineGrid    ,&
         GridTarget                 = SP_Grid             ,&
         Router                     = RouterLineIhSp      ,&
         n_interface_point_in_block = IH_n_particle       ,&
         interface_point_coords     = i_particle          ,&
         mapping                    = mapping_line_ih_to_sp)
    call couple_comp(RouterLineIhSp                       ,&
         nVar = 3                                         ,&
         fill_buffer = IH_get_coord_for_sp_and_transform  ,&
         apply_buffer= MF_put_line)
  end subroutine SP_get_lines_from_ih
  !============================================================================
  subroutine mapping_sp_to_ih(nDimIn, XyzIn_D, nDimOut, CoordOut_D, &
       IsInterfacePoint)
    integer, intent(in) :: nDimIn, nDimOut
    real,    intent(in) :: XyzIn_D(nDimIn)
    real,    intent(out):: CoordOut_D(nDimOut)
    logical, intent(out):: IsInterfacePoint
    real                :: XyzTemp_D(nDim)
    !--------------------------------------------------------------------------
    IsInterfacePoint = .true.
    XyzTemp_D = matmul(XyzIn_D, IhToSp_DD)*UnitSp2UnitIh
    call IH_xyz_to_coord(XyzTemp_D, CoordOut_D)
  end subroutine mapping_sp_to_ih
  !============================================================================
  subroutine i_particle(nDim, Xyz_D, nIndex, iIndex_I,IsInterfacePoint)
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
  subroutine mapping_line_ih_to_sp(nDimIn, XyzIn_D, &
       nDimOut, CoordOut_D, IsInterfacePoint)
    integer, intent(in) :: nDimIn, nDimOut
    real,    intent(in) :: XyzIn_D(nDimIn)
    real,    intent(out):: CoordOut_D(nDimOut)
    logical, intent(out):: IsInterfacePoint

    integer:: iIndex_I(2), iParticle
    !--------------------------------------------------------------------------
    ! Convert back the above converted particle number
    IsInterfacePoint = .true.; iParticle = nint(XyzIn_D(1))
    call IH_get_particle_indexes(iParticle, iIndex_I)
    ! Convert particle indexes to SP gen coords, so that the nearest_cell
    ! program in the router coud figure out where it should be put in SP
    CoordOut_D = coord_grid_d(SP_Grid, iIndex_I(1), [iIndex_I(2),1,1])
  end subroutine mapping_line_ih_to_sp
  !============================================================================
  subroutine IH_get_for_sp_and_transform(&
       nPartial, iGetStart, Get, w, State_V, nVar)
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V

    integer:: iVarBx, iVarBz, iVarMx, iVarMz
    !--------------------------------------------------------------------------
    call IH_get_for_mh(nPartial,iGetStart,Get,w,State_V,nVar)
    ! indices of variables
    iVarBx = iVar_V(BxCouple_);   iVarBz = iVar_V(BzCouple_)
    iVarMx = iVar_V(RhoUxCouple_);iVarMz = iVar_V(RhoUzCouple_)
    ! perform transformation before returning
    State_V(iVarBx:iVarBz) = matmul(IhToSp_DD,State_V(iVarBx:iVarBz))
    State_V(iVarMx:iVarMz) = matmul(IhToSp_DD,State_V(iVarMx:iVarMz))
  end subroutine IH_get_for_sp_and_transform
  !============================================================================
  subroutine IH_get_coord_for_sp_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V

    !--------------------------------------------------------------------------
    call IH_get_particle_coords(Get%iCB_II(1, iGetStart), State_V)
    ! perform transformation before returning
    State_V = matmul(IhToSp_DD, State_V)*UnitIh2UnitSp
  end subroutine IH_get_coord_for_sp_and_transform
  !============================================================================
end module CON_couple_mh_sp
