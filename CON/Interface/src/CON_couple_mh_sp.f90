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
       SP_do_extract_lines      ,&  ! If returns .true., extract the mf lines  
       SP_put_coupling_param    ,&  ! Set time and interaface bounds
       SP_adjust_lines              ! Process if needed the updated mf lines
  use PT_wrapper, ONLY: &
       PT_do_extract_lines      ,&  ! If returns .true., extract the mf lines  
       PT_put_coupling_param    ,&  ! Set time and interaface bounds
       PT_adjust_lines              ! Process if needed the updated mf lines
  use CON_bline, ONLY:  BL_,     &
       RScMin, RScMax,           &  !^CMP IF SC  
       RIhMin, RIhMax,           &
       BL_init_foot_points,      &  ! Initialize footpoint array 
       BL_get_bounds,            &  ! Provides RScMin/Max and/or RIhMin/Max
       BL_n_particle            ,&  ! Number of "points" in a given line in SP
       BL_put_from_mh           ,&  ! Put MHD info from SC or IH to SP
       BL_interface_point_coords,&  ! Check if the point is within interface
       BL_put_line                  ! Put particle Xyz from SC/IH to SP
                 
  implicit none

  private ! Except
  public::couple_mh_sp_init
  public::couple_ih_sp
  public::couple_sc_sp              !^CMP IF SC

  type(GridType)     , save::BL_Grid           ! Target (Particle coords)
  type(LocalGridType), save::BL_LocalGrid      ! Target (MHD data)
  type(GridType)     , save::IH_Grid           ! Source (MHD data)
  type(RouterType),save,private::RouterIhBl    ! IH (MHD data) => SP
  type(GridType)     , save::IH_LineGrid       ! Misc
  type(LocalGridType), save::IH_LocalLineGrid  ! Source (MHD data)
  type(RouterType),save,private::RouterLineIhBl    ! IH (Particle coords)=>SP
  !^CMP IF SC BEGIN
  type(GridType)     , save::SC_Grid           ! Source (MHD data)
  type(RouterType),save,private::RouterScBl        ! SC MHD data => BL
  type(GridType)     , save::SC_LineGrid       ! Misc
  type(LocalGridType), save::SC_LocalLineGrid  ! Source (Particle Coords)
  type(RouterType),save,private::RouterLineScBl    ! SC (Particle coords)=>BL
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
  real    :: ScToBl_DD(3,3), UnitBl2UnitSc, UnitSc2UnitBl !^CMP IF SC
  real    :: IhToBl_DD(3,3), UnitBl2UnitIh, UnitIh2UnitBl

  logical::DoTest,DoTestMe
  character(LEN=*),parameter::NameSub='couple_mh_sp'
  real :: tNow

  logical ::DoInit=.true.
contains
  !============================================================================
  subroutine BL_do_extract_lines(DoExtract)
    logical, intent(out) ::  DoExtract
    !--------------------------------------------------------------------------
    select case(BL_)
       case(SP_)
          call SP_do_extract_lines(DoExtract)
       case(PT_)
          call PT_do_extract_lines(DoExtract)
       case default
          call CON_stop('The target model it BL is not allowed')
       end select
  end subroutine BL_do_extract_lines
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
  subroutine couple_mh_sp_init
    use CON_physics, ONLY: get_time
    use ModConst
    ! whether need to extract new line (equal to .not.DoRestart in BL
    logical :: DoExtract
    !--------------------------------------------------------------------------
    if(.not.DoInit)RETURN
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    ! The initialization can be done only once

    call get_time(tSimulationOut=tNow)

    ! determine whether BL is ready for coupling with MH
    call BL_do_extract_lines(DoExtract)
    if(use_comp(SC_))call couple_sc_sp_init  !^CMP IF SC
    if(use_comp(IH_))call couple_ih_sp_init
    if(DoExtract.and.is_proc(BL_))call BL_init_foot_points

    ! After couple_mh_sp_init
    ! (a) set of low and upper boundaries in the coupler
    ! (b) set of the Lagrangian point coordinates in BL
    ! (c) logical DoInit is .true. to prevent the coordinate
    ! info on BL to be overritten with that from SC or IH
  contains
    !==========================================================================
    subroutine couple_sc_sp_init             !^CMP IF SC BEGIN
      logical :: IsReady
      integer, parameter :: Lower_=0 
      !------------------------------------------------------------------------
      call SC_check_ready_for_sp(IsReady)
      if(.not.IsReady) call CON_stop(&
           "SC component not ready for "//NameSub//" correct PARAM.in")
      call set_couple_var_info(SC_, BL_)
      UnitBl2UnitSc = Grid_C(BL_)%UnitX/Grid_C(SC_)%UnitX
      UnitSc2UnitBl = 1/UnitBl2UnitSc

      ! Initialize coupler from SC (source )to BL (target)
      ! Data will be copied from SC to BL, however, the points
      ! in which the data shuld be provided will be sent from BL to
      ! SC, at the stage of the router construction. To further benefit
      ! from this opportunity, two BL grid indexes are also sent from
      ! BL to SC,at this stage, therefore nMappedPointIndex = 2.
      call init_coupler(iCompSource = SC_          ,&
           iCompTarget       = BL_                 ,&
           GridSource        = SC_Grid             ,&
           GridTarget        = BL_Grid             ,&
           nMappedPointIndex = nAux                ,&
           Router            = RouterScBl          ,&
           LocalGridTarget   = BL_LocalGrid)

      ! Set local GD on the Particle_I structure
      call set_standard_grid_descriptor(SC_LineDD,Grid=SC_LineGrid)
      call init_router(SC_LineGrid, BL_Grid, RouterLineScBl)
      if(.not.RouterScBl%IsProc)RETURN
      if(is_proc(SC_))call set_local_gd(iProc = i_proc(),   &
           Grid = SC_LineGrid, LocalGrid = SC_LocalLineGrid)
      ! Router to send particles is initialized.
      if(.not.DoExtract)RETURN
      !
      ! If DoExtract, then
      !
      if(is_proc(BL_))then
         call BL_put_coupling_param(SC_, tNow)
         call BL_get_bounds(RScMin, RScMax)
      end if
      call SC_synchronize_refinement(RouterScBl%iProc0Source,&
           RouterScBl%iCommUnion)
      ScToBl_DD=transform_matrix(tNow,&
           Grid_C(SC_)%TypeCoord, Grid_C(BL_)%TypeCoord)
      ! use router intended for sending MHD data  (SC => BL)
      ! to send the particle coordinates backward (BL=>SC).

      ! Put in place the origin points of the MF lines
      call set_router(&
           GridSource            = SC_Grid                  ,&
           GridTarget            = BL_LocalGrid             ,&
           Router                = RouterScBl               ,&
           n_interface_point_in_block = BL_n_particle       ,&
           interface_point_coords= BL_interface_point_coords,&
           mapping               = mapping_sp_to_sc         ,&
           interpolate           = interpolation_amr_gc)
      if(is_proc(SC_))then
         nLength = nlength_buffer_source(RouterScBl)
         call SC_extract_line(&
              Xyz_DI     = RouterScBl%BufferSource_II(1:nDim,&
              1:nLength)                                    ,&
              iTraceMode = Lower_                           ,&
              iIndex_II  = nint(RouterScBl% BufferSource_II( &
              nDim+1:nDim+nAux,1:nLength))                  ,&
              RSoftBoundary = RScMax*UnitBl2UnitSc)
         !
         ! Now in SC are the parts of lines from the inner
         ! boundary of SC to RScMax + 1 Point above
         !
      end if
      call BL_put_lines_from_sc
      !
      ! Now in BL are the parts of lines from RScMin 
      ! to RScMax + 1 point above
      !
    end subroutine couple_sc_sp_init                 !^CMP END SC
    !==========================================================================
    subroutine couple_ih_sp_init
      logical :: IsReady
      integer, parameter:: Upper_ = 1
      !------------------------------------------------------------------------
      ! Check whether IH is ready for coupling with BL;
      call IH_check_ready_for_sp(IsReady)
      if(.not.IsReady) call CON_stop(&
           "IH component not ready for "//NameSub//" correct PARAM.in")
      call set_couple_var_info(IH_, BL_)
      UnitBl2UnitIh = Grid_C(BL_)%UnitX/Grid_C(IH_)%UnitX
      UnitIh2UnitBl = 1/UnitBl2UnitIh

      ! Initialize coupler from IH (source )to BL (target)
      ! Data will be copied from IH to BL, however, the points
      ! in which the data shuld be provided will be sent from BL to
      ! IH, at the stage of the router construction. To further benefit
      ! from this opportunity, two BL grid indexes are also sent from
      ! BL to IH, at this stage, therefore nMappedPointIndex = 2.
      if(use_comp(SC_))then          !^CMP IF SC BEGIN
         ! Do not reset already installed BL_LocalGrid
         call init_coupler(iCompSource = IH_          ,&
              iCompTarget          = BL_              ,&
              GridSource           = IH_Grid          ,&
              GridTarget           = BL_Grid          ,&
              nMappedPointIndex    = nAux             ,&
              Router               = RouterIhBl)
      else                           !^CMP END SC
         call init_coupler(iCompSource = IH_          ,&
              iCompTarget          = BL_              ,&
              GridSource           = IH_Grid          ,&
              GridTarget           = BL_Grid          ,&
              nMappedPointIndex    = nAux             ,&
              Router               = RouterIhBl       ,&
              LocalGridTarget      = BL_LocalGrid)
      end if                         !^CMP IF SC
      ! Set local GD on the Particle_I structure
      call set_standard_grid_descriptor(IH_LineDD,Grid=&
           IH_LineGrid)
      call init_router(IH_LineGrid, BL_Grid, RouterLineIhBl)
      if(.not.RouterIhBl%IsProc)RETURN
      if(is_proc(IH_))call set_local_gd(iProc = i_proc(),&
           Grid      = IH_LineGrid                      ,&
           LocalGrid = IH_LocalLineGrid)
      ! Router to send particles is initialized.
      if(.not.DoExtract)RETURN
      if(is_proc(BL_))then
         call BL_put_coupling_param(&
              Source_     = IH_,    &
              TimeIn      = tNow)
         if(use_comp(SC_))then                   !^CMP IF SC BEGIN
            call BL_get_bounds(RScMax, RIhMax)
         else                                    !^CMP END SC
            call BL_get_bounds(RIhMin, RIhMax)
         end if                                  !^CMP IF SC                          
      end if
      call IH_synchronize_refinement(RouterIhBl%iProc0Source,&
           RouterIhBl%iCommUnion)
      IhToBl_DD=transform_matrix(tNow,&
           Grid_C(IH_)%TypeCoord, Grid_C(BL_)%TypeCoord)
      ! Put in place the origin points of the MF lines.
      call set_router(                                            &
           GridSource                 = IH_Grid                  ,&
           GridTarget                 = BL_LocalGrid             ,&
           Router                     = RouterIHBl               ,&
           n_interface_point_in_block = BL_n_particle            ,&
           interface_point_coords     = BL_interface_point_coords,&
           mapping                    = mapping_sp_to_IH         ,&
           interpolate                = interpolation_amr_gc)
      !
      ! Now in IH are 1 point above RScMax per each line
      !
      if(is_proc(IH_))then
         nLength = nlength_buffer_source(RouterIhBl)
         call IH_extract_line(Xyz_DI  = RouterIhBl%               &
              BufferSource_II(1:nDim,1:nLength)                  ,&
              iTraceMode              = Upper_                   ,&
              iIndex_II               = nint(RouterIhBl%          &
              BufferSource_II(nDim+1:nDim+nAux,1:nLength))       ,&
              RSoftBoundary           = RIhMax*UnitBl2UnitIh)
         !
         ! Now in IH are the parts of lines from RScMax + 1 point above
         ! to RIhMax + 1 Point above
         !
      end if
      call BL_put_lines_from_ih
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
    if(.not.RouterScBl%IsProc)RETURN
    if(DoTest.and.is_proc0(SC_))&
         write(*,'(a,es12.5)')NameSub//': couple to SC,time=', DataInputTime
    tNow=DataInputTime
    if(is_proc(BL_))then
       call BL_put_coupling_param(&
            Source_ = SC_,        &
            TimeIn  = DataInputTime)
       call BL_get_bounds(&
            rMinIn      = RScMin,&
            rMaxIn      = RScMax,&
            rBufferUpIn = RIhMin)
    end if
    !
    ! Coordinates for all particles are nullified. Those which are
    ! not sent (or lost) keep to be zeroed and may be thus found.
    !
    ScToBl_DD=transform_matrix(tNow,&
         Grid_C(SC_)%TypeCoord, Grid_C(BL_)%TypeCoord)
    call SC_synchronize_refinement(RouterScBl%iProc0Source,&
         RouterScBl%iCommUnion)
    if(.not.DoInit) call BL_put_lines_from_sc
    !
    ! The advected particles from SC are taken  
    !
    if(is_proc(BL_))call BL_adjust_lines(SC_)
    !
    ! Some particles may be lost, after adjustment the line intergity
    ! is recovered
    ! 
    ! Set router SC=> BL to  receive MHD data
    call set_router(&
         GridSource             = SC_Grid                  ,&
         GridTarget             = BL_LocalGrid             ,&
         Router                 = RouterScBl               ,&
         n_interface_point_in_block = BL_n_particle        ,&
         interface_point_coords = BL_interface_point_coords,&
         mapping                = mapping_sp_to_sc         ,&
         interpolate            = interpolation_amr_gc)
    call couple_comp(RouterScBl                            ,&
         nVar                   = nVarBuffer               ,&
         fill_buffer = SC_get_for_sp_and_transform         ,&
         apply_buffer           = BL_put_from_mh)
    !
    ! The MHD data within the heliocentric distances RScMin<R<RScMax 
    ! By the way get particle coords from the router buffer to SC
    if(is_proc(SC_))then
       nLength = nlength_buffer_source(RouterScBl)
       call SC_put_particles(Xyz_DI =                       &
            RouterScBl%BufferSource_II(1:nDim, 1:nLength)  ,&
            iIndex_II = nint(RouterScBl%BufferSource_II(    &
            nDim+1:nDim+nAux, 1:nLength)))
    end if
  end subroutine couple_sc_sp
  !============================================================================
  subroutine BL_put_lines_from_sc
    !--------------------------------------------------------------------------
    call construct_router_from_source(                      &
         GridSource                 = SC_LocalLineGrid     ,&
         GridTarget                 = BL_Grid              ,&
         Router                     = RouterLineScBl       ,&
         n_interface_point_in_block = SC_n_particle        ,&
         interface_point_coords     = i_particle           ,&
         mapping                    = mapping_line_sc_to_sp)
    call couple_comp(                 RouterLineScBl       ,&
         nVar                       = nDim                 ,&
         fill_buffer = SC_get_coord_for_sp_and_transform   ,&
         apply_buffer               = BL_put_line)
  end subroutine BL_put_lines_from_sc
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
    XyzTemp_D = matmul(XyzIn_D, ScToBl_DD)*UnitBl2UnitSc
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
    CoordOut_D = coord_grid_d(BL_Grid, iIndex_I(1),[iIndex_I(2),1,1])
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
    iVarBx = iVar_V(BxCouple_)   ; iVarBz = iVar_V(BzCouple_)
    iVarMx = iVar_V(RhoUxCouple_); iVarMz = iVar_V(RhoUzCouple_)
    ! perform transformation before returning
    State_V(iVarBx:iVarBz) = matmul(ScToBl_DD,State_V(iVarBx:iVarBz))
    State_V(iVarMx:iVarMz) = matmul(ScToBl_DD,State_V(iVarMx:iVarMz))
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
    State_V = matmul(ScToBl_DD, State_V)*UnitSc2UnitBl
  end subroutine SC_get_coord_for_sp_and_transform
  !============================================================================
  !!^CMP END SC
  subroutine couple_ih_sp(DataInputTime)
    real,intent(in)::DataInputTime
    !--------------------------------------------------------------------------
    if(.not.RouterIhBl%IsProc)RETURN
    if(DoTest.and.is_proc0(IH_))&
         write(*,'(a,es12.5)')NameSub//': couple to IH,time=', DataInputTime
    tNow = DataInputTime
    if(is_proc(BL_))then
       call BL_put_coupling_param(&
            Source_ = IH_,        &
            TimeIn  = DataInputTime)
       call BL_get_bounds(        &
            rMinIn      = RIhMin, &
            rMaxIn      = RIhMax, &
            rBufferLoIn = RScMax)
    end if
    IhToBl_DD=transform_matrix(tNow,&
         Grid_C(IH_)%TypeCoord, Grid_C(BL_)%TypeCoord)
    call IH_synchronize_refinement(RouterIhBl%iProc0Source,&
         RouterIhBl%iCommUnion)
    if(.not.DoInit) call BL_put_lines_from_ih
    !
    ! Now the full magnetic line is available, including probably
    ! some points lost in IH
    ! 
    if(is_proc(BL_))call BL_adjust_lines(IH_)
    ! In points at RIhMin < R < RIhMax get the MHD data
    ! Set router IH=> BL to  receive MHD data
    call set_router(                                            &
         GridSource                 = IH_Grid                  ,&
         GridTarget                 = BL_LocalGrid             ,&
         Router                     = RouterIHBl               ,&
         n_interface_point_in_block = BL_n_particle            ,&
         interface_point_coords     = BL_interface_point_coords,&
         mapping                    = mapping_sp_to_ih         ,&
         interpolate                = interpolation_amr_gc)
    call couple_comp(                 RouterIhBl               ,&
         nVar                       = nVarBuffer               ,&
         fill_buffer              = IH_get_for_sp_and_transform,&
         apply_buffer               = BL_put_from_mh)
    ! By the way get coords for particles in IH from the router buffer
    if(is_proc(IH_))then
       nLength = nlength_buffer_source(RouterIhBl)
       if(use_comp(SC_))call sort_out_sc_particles     !^CMP IF SC
       call IH_put_particles(Xyz_DI =  RouterIhBl%              &
            BufferSource_II(1:nDim, 1:nLength)                 ,&
            iIndex_II               = nint(RouterIhBl%          &
            BufferSource_II(nDim+1:nDim+nAux, 1:nLength)))
    end if
    DoInit=.false.
  contains                                       !^CMP IF SC BEGIN
    !==========================================================================
    subroutine sort_out_sc_particles
      ! Sort out particles to be advected by the SC, which have
      ! R2 in BL coordinates less than RScMax**2
      integer:: iParticle, iParticleNew
      real   :: Xyz_D(3), Coord_D(3), RScMax2

      !------------------------------------------------------------------------
      iParticleNew = 0; RScMax2 = (RScMax*UnitBl2UnitIh)**2
      do iParticle = 1, nLength
         Coord_D = RouterIhBl%BufferSource_II(1:nDim, iParticle)
         call IH_coord_to_xyz(Coord_D, Xyz_D)
         if(sum(Xyz_D**2) < RScMax2)CYCLE
         iParticleNew  = iParticleNew + 1
         RouterIhBl%BufferSource_II(:, iParticleNew) = &
              RouterIhBl%BufferSource_II(:, iParticle)
      end do
      nLength = iParticleNew
    end subroutine sort_out_sc_particles                 !^CMP END SC
    !==========================================================================
  end subroutine couple_ih_sp
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
  subroutine mapping_sp_to_ih(nDimIn, XyzIn_D, nDimOut, CoordOut_D, &
       IsInterfacePoint)
    integer, intent(in) :: nDimIn, nDimOut
    real,    intent(in) :: XyzIn_D(nDimIn)
    real,    intent(out):: CoordOut_D(nDimOut)
    logical, intent(out):: IsInterfacePoint
    real                :: XyzTemp_D(nDim)
    !--------------------------------------------------------------------------
    IsInterfacePoint = .true.
    XyzTemp_D = matmul(XyzIn_D, IhToBl_DD)*UnitBl2UnitIh
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
    ! Convert particle indexes to BL gen coords, so that the nearest_cell
    ! program in the router coud figure out where it should be put in BL
    CoordOut_D = coord_grid_d(BL_Grid, iIndex_I(1), [iIndex_I(2),1,1])
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
    State_V(iVarBx:iVarBz) = matmul(IhToBl_DD,State_V(iVarBx:iVarBz))
    State_V(iVarMx:iVarMz) = matmul(IhToBl_DD,State_V(iVarMx:iVarMz))
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
    State_V = matmul(IhToBl_DD, State_V)*UnitIh2UnitBl
  end subroutine IH_get_coord_for_sp_and_transform
  !============================================================================
end module CON_couple_mh_sp
