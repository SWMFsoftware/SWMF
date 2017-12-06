!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE SP
!^CMP FILE IH
module CON_couple_mh_sp
  use CON_axes
  use CON_coupler
  use IH_wrapper, ONLY: IH_synchronize_refinement, IH_extract_line,           &
       IH_get_for_sp, IH_get_a_line_point, IH_add_to_line, IH_n_particle,     &
       IH_line_interface_point, IH_get_particle_indexes, IH_check_particles,  &
        IH_LineDD, IH_get_particle_coords, IH_xyz_to_coord, IH_coord_to_xyz   
  !^CMP IF SC BEGIN 
  use SC_wrapper, ONLY: SC_synchronize_refinement, SC_extract_line,           &
       SC_get_for_sp, SC_get_a_line_point, SC_add_to_line, SC_n_particle,     &
       SC_line_interface_point, SC_get_particle_indexes, SC_check_particles,  & 
       SC_LineDD, SC_get_particle_coords, SC_xyz_to_coord    
  !^CMP END SC
  use SP_wrapper, ONLY: &
       SP_put_from_sc, SP_put_from_ih, SP_put_input_time, SP_put_line,        &
       SP_n_particle, SP_check_if_do_extract, SP_get_bounds_comp,             &
       SP_interface_point_coords_for_ih, SP_interface_point_coords_for_sc,    &
       SP_interface_point_coords_for_ih_extract, SP_set_line_foot,            &
       SP_copy_old_state, SP_adjust_lines, SP_assign_lagrangian_coords
       
  implicit none
  
  private !Except
  public::couple_mh_sp_init
  public::couple_ih_sp              
  public::couple_sc_sp              !^CMP IF SC

  type(GridType),save::SP_Grid           !Target (Particle coords)
  type(LocalGridType),       save::SP_LocalGrid      !Target (MHD data)
  type(GridType),save::IH_Grid           !Source (MHD data)
  type(RouterType),save,private::RouterIhSp        !IH (MHD data) => SP 
  type(GridType),save::IH_LineGrid       !Misc
  type(LocalGridType),       save::IH_LocalLineGrid  !Source (MHD data) 
  type(RouterType),save,private::RouterLineIhSp    !IH (Particle coords)=>SP    
  !^CMP IF SC BEGIN
  type(GridType),save::SC_Grid !Source (MHD data)  
  type(RouterType),save,private::RouterScSp        !SC MHD data => SP
  type(GridType),save::SC_LineGrid       !Misc
  type(LocalGridType),       save::SC_LocalLineGrid  !Source (Particle Coords) 
  type(RouterType),save,private::RouterLineScSp    !SC (Particle coords)=>SP    
  !^CMP END SC
  !\
  !Three-dimensional grids in SC and IH, 
  !three-component coordinate vector
  !/
  integer, parameter :: nDim = 3
  !\
  ! We send: (1) the global index of the magnetic field line
  ! and (2) the particle number along the magnetic field line
  !/ 
  integer, parameter :: nAux = 2
  !\
  ! Misc
  integer :: nLength, iError
  !\
  ! Transformation matrices
  real :: ScToSp_DD(3,3) !^CMP IF SC
  real :: IhToSp_DD(3,3)
 
  logical::DoTest,DoTestMe
  character(LEN=*),parameter::NameSub='couple_mh_sp'
  real :: tNow
 
  integer, parameter :: Lower_ = 1 
  !^CMP IF SC BEGIN
  ! Short ID of MHD components. If there is only one component,
  ! it has an identifier Lower_, otherwise, SC has identifier Lower_
  ! while IH is Upper
  integer, parameter :: Upper_ = 2 
  ! solar corona and lower and upper boundaries
  real:: RScMin, RScMax             !^CMP END SC
  ! inner heliosphere lower and upper boundaries
  real:: RIhMin, RIhMax
contains
  !==================================================================
  subroutine couple_mh_sp_init
    use CON_physics, ONLY: get_time
    use ModConst
    logical,save::DoInit=.true. 
    integer, parameter:: &
         iInterfaceOrigin = 0, iInterfaceEnd = 1
    ! whether need to extract new line (equal to .not.DoRestart
    logical:: DoExtract
    !----------------------------------------------------------------------
    if(.not.DoInit)return
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    !The initialization can be done only once

    call get_time(tSimulationOut=tNow)
    if(is_proc(SP_))call SP_put_input_time(tNow)

    ! determine whether need to extract new field lines
    call SP_check_if_do_extract(DoExtract)

    if(use_comp(SC_))call couple_sc_sp_init  !^CMP IF SC         
    if(use_comp(IH_))call couple_ih_sp_init
    if(DoExtract.and.is_proc(SP_))then
       call SP_set_line_foot
       call SP_assign_lagrangian_coords
    end if
    DoInit=.false.
  contains
    !===============================
    subroutine couple_sc_sp_init             !^CMP IF SC BEGIN
      ! get the value of SC boundaries as set in SP
      call SP_get_bounds_comp(Lower_, RScMin, RScMax)
      ! account for coord unit difference
      RScMin = RScMin * Grid_C(SP_)%UnitX / Grid_C(SC_)%UnitX
      RScMax = RScMax * Grid_C(SP_)%UnitX / Grid_C(SC_)%UnitX
      call SC_check_particles()
      call set_couple_var_info(SC_, SP_)
      !\
      !Initialize coupler from SC (source )to SP (target)
      !Data will be copied from SC to SP, however, the points
      !in which the data shuld be provided will be sent from SP to
      !SC, at the stage of the router construction. To further benefit
      !from this opportunity, two SP grid indexes are also sent from
      !SP to SC,at this stage, therefore nMappedPointIndex = 2. 
      call init_coupler(iCompSource = SC_,          &
           iCompTarget = SP_,                       &
           GridSource = SC_Grid,&
           GridTarget = SP_Grid,&
           nMappedPointIndex    = nAux,             &
           Router          = RouterScSp        ,    &
           LocalGridTarget   = SP_LocalGrid)
      !/
      !\
      ! Set local GD on the Particle_I structure
      call set_standard_grid_descriptor(SC_LineDD,Grid=SC_LineGrid)
      call init_router(SC_LineGrid, SP_Grid, RouterLineScSp)
      if(.not.RouterScSp%IsProc)RETURN
      if(is_proc(SC_))call set_local_gd(iProc = i_proc(),   &
           Grid = SC_LineGrid, LocalGrid = SC_LocalLineGrid)
      !Router to send particles is initialized. 
      call SC_synchronize_refinement(RouterScSp%iProc0Source,&
           RouterScSp%iCommUnion)
      ScToSp_DD=transform_matrix(tNow,&                 
           Grid_C(SC_)%TypeCoord, Grid_C(SP_)%TypeCoord)
      if(DoExtract)then
         !use router intended for sending MHD data  (SC => SP) 
         !to send the particle coordinates backward (SP=>SC).
         call set_router(&
              GridSource              = SC_Grid, &
              GridTarget              = SP_LocalGrid, &
              Router                = RouterScSp, &
              n_interface_point_in_block = SP_n_particle,&
              interface_point_coords= SP_interface_point_coords_for_sc, &
              mapping               = mapping_sp_to_sc, &
              interpolate           = interpolation_amr_gc)
         if(is_proc(SC_))then
            !Put in place the origin points of the MF lines.
            nLength = nlength_buffer_source(RouterScSp)
            call SC_extract_line(&
                 nLine             = nLength, &
                 XyzOrigin_DI      = &
                 RouterScSp%BufferSource_II(1:nDim,1:nLength), &
                 iTraceMode        = iInterfaceOrigin, &
                 nIndex            = nAux, &
                 iIndexOrigin_II   = nint(RouterScSp%&
                 BufferSource_II(nDim+1:nDim+nAux,1:nLength)),&
                 RSoftBoundaryIn   =  RScMax, & 
                 UseInputInGenCoord= .true.)
         end if
      end if
      !First coupling with the particle info and data exchange
      call exchange_data_sc_sp(DoInit, DoExtract)
    end subroutine couple_sc_sp_init                 !^CMP END SC
    !==============================
    subroutine couple_ih_sp_init           
      ! get the value of IH boundaries as set in SP
      if(use_comp(SC_))then  !^CMP IF SC BEGIN
         call SP_get_bounds_comp(Upper_, RIhMin, RIhMax)
      else                   !^CMP END SC
         call SP_get_bounds_comp(Lower_, RIhMin, RIhMax)
      end if                 !^CMP IF SC
      ! account for coord unit difference
      RIhMin = RIhMin * Grid_C(SP_)%UnitX / Grid_C(IH_)%UnitX
      RIhMax = RIhMax * Grid_C(SP_)%UnitX / Grid_C(IH_)%UnitX
      call IH_check_particles()
      call set_couple_var_info(IH_, SP_)
      !\
      !Initialize coupler from IH (source )to SP (target)
      !Data will be copied from IH to SP, however, the points
      !in which the data shuld be provided will be sent from SP to
      !IH, at the stage of the router construction. To further benefit
      !from this opportunity, two SP grid indexes are also sent from
      !SP to IH, at this stage, therefore nMappedPointIndex = 2.
      call init_coupler(iCompSource = IH_,          &
           iCompTarget = SP_,                       &
           GridSource = IH_Grid,          &
           GridTarget = SP_Grid,          &
           nMappedPointIndex    = nAux,             &
           Router          = RouterIhSp)
      ! Set local GD on the Particle_I structure
      call set_standard_grid_descriptor(IH_LineDD,Grid=&
           IH_LineGrid)
      call init_router(IH_LineGrid, SP_Grid, RouterLineIhSp)
      if(.not.RouterIhSp%IsProc)RETURN
      if(.not.use_comp(SC_))then !^CMP IF SC BEGIN
         if(is_proc(SP_))call set_local_gd(i_proc(),&
              SP_Grid, SP_LocalGrid)
      end if                     !^CMP END SC
      if(is_proc(IH_))call set_local_gd(&
           iProc = i_proc(), &
           Grid = IH_LineGrid, &
           LocalGrid = IH_LocalLineGrid)
      !Router to send particles is initialized. 
      call IH_synchronize_refinement(RouterIhSp%iProc0Source,&
           RouterIhSp%iCommUnion)
      IhToSp_DD=transform_matrix(tNow,&                 
           Grid_C(IH_)%TypeCoord, Grid_C(SP_)%TypeCoord)
      if(DoExtract)then
         !use the router intended for sending MHD data IH => SP 
         !to send the particle coordinates backward (SP=>IH)
         if(is_proc(SP_))&
              call set_semi_router_from_target(&
              GridSource  = IH_Grid, &
              GridTarget  = SP_LocalGrid, &
              Router                = RouterIHSp, &
              n_interface_point_in_block = SP_n_particle,&
              interface_point_coords= &
              SP_interface_point_coords_for_ih_extract, &
              mapping               = mapping_sp_to_IH, &
              interpolate           = interpolation_amr_gc) 
         call synchronize_router_target_to_source(RouterIhSp)
         if(is_proc(IH_))then
            !Put in place the origin points of the MF lines.
            nLength = nlength_buffer_source(RouterIhSp)
            call IH_extract_line(&
                 nLine             = nLength,                                &
                 XyzOrigin_DI = RouterIhSp%BufferSource_II(1:nDim,1:nLength),&
                 iTraceMode        = iInterfaceEnd,                          &
                 nIndex            = nAux,                                   &
                 iIndexOrigin_II   = nint(RouterIhSp%&
                 BufferSource_II(nDim+1:nDim+nAux,1:nLength)),               &
                 RSoftBoundaryIn   =  RIhMax,                                & 
                 UseInputInGenCoord= .true.)
         end if
      end if
      !First coupling with the particle info and data exchange
      call exchange_data_ih_sp(DoInit, DoExtract)
    end subroutine couple_ih_sp_init
  end subroutine couple_mh_sp_init
  !=========================================================================
  !^CMP IF SC BEGIN
  subroutine couple_sc_sp(DataInputTime)
    use CON_global_message_pass

    real,intent(in)::DataInputTime
    integer:: nLength
    !-------------------------------------------------------
    if(.not.RouterScSp%IsProc)return

    tNow=DataInputTime
    if(is_proc(SP_))call SP_put_input_time(DataInputTime)

    ! IMPORTANT: 
    ! couple_sc_sp is called BEFORE couple_ih_sp; 
    ! save the current state as old,
    ! separate subroutine is used in order to avoid intersections 
    ! as fluid elements transfer from SC to IH
    if(is_proc(SP_))call SP_copy_old_state
    
    ScToSp_DD=transform_matrix(tNow,&
         Grid_C(SC_)%TypeCoord, Grid_C(SP_)%TypeCoord)
    call SC_synchronize_refinement(RouterScSp%iProc0Source,&
         RouterScSp%iCommUnion)
    call exchange_data_sc_sp(DoInit=.false., DoExtract=.false.)
  end subroutine couple_sc_sp
  !=====================
  subroutine exchange_data_sc_sp(DoInit, DoExtract)
    logical, intent(in) :: DoInit, DoExtract
    !----------------------------------------------------------
    if(DoExtract.or..not.DoInit)then
       !\
       !Lagrangian particles are either extracted or 
       !updated in SC
       !/
       ! Send Largrangian particle coords SC => SP
       call construct_router_from_source(&
            GridSource = SC_LocalLineGrid,     &
            GridTarget = SP_Grid,  &
            Router               = RouterLineScSp,     &
            n_interface_point_in_block = SC_n_particle,&
            interface_point_coords=SC_line_interface_point,&
            mapping              = mapping_line_sc_to_sp)
       call couple_comp(RouterLineScSp, &
            nVar = 3, &
            fill_buffer = SC_get_line_for_sp_and_transform, &
            apply_buffer= SP_put_line)
    end if
    if(is_proc(SP_).and..not.DoInit)call SP_adjust_lines(&
         DoAdjustStart=.true.,DoAdjustEnd=.not.use_comp(IH_))
    
    ! Set router SC=> SP to  receive MHD data 
    call set_router(&
         GridSource  = SC_Grid,                       &
         GridTarget  = SP_LocalGrid,                  &
         Router                = RouterScSp,        &
         n_interface_point_in_block = SP_n_particle,&
         interface_point_coords=                    &
         SP_interface_point_coords_for_sc,          &
         mapping               = mapping_sp_to_sc,  &
         interpolate           = interpolation_amr_gc)
    call couple_comp(RouterScSp,                    &
         nVar = nVarBuffer,                         &
         fill_buffer = SC_get_for_sp_and_transform, &
         apply_buffer= SP_put_from_sc)
    !By the way get particle coords from the router buffer 
    if(is_proc(SC_))then
       nLength = nlength_buffer_source(RouterScSp)
       call SC_add_to_line(                                        &
            nParticle = nLength,                                   &
            Xyz_DI = RouterScSp%BufferSource_II(1:nDim, 1:nLength),&
            nIndex    = nAux,                                      &
            iIndex_II = nint(RouterScSp%BufferSource_II(           &
            nDim+1:nDim+nAux, 1:nLength)),                         &      
            UseInputInGenCoord = .true., DoReplace = .true.)
    end if
  end subroutine exchange_data_sc_sp
  !=========================
  subroutine mapping_sp_to_sc(nDimIn, XyzIn_D, nDimOut, CoordOut_D, &
       IsInterfacePoint)
    integer, intent(in) :: nDimIn
    real,    intent(in) :: XyzIn_D(nDimIn)
    integer, intent(in) :: nDimOut
    real,    intent(out):: CoordOut_D(nDimOut)
    logical, intent(out):: IsInterfacePoint

    real                :: XyzTemp_D(nDim)
    !------------------------------------------
    IsInterfacePoint = .true.; XyzTemp_D=matmul(XyzIn_D, ScToSp_DD)
    XyzTemp_D = XyzTemp_D * Grid_C(SP_)%UnitX / Grid_C(SC_)%UnitX
    call SC_xyz_to_coord(XyzTemp_D, CoordOut_D)
  end subroutine mapping_sp_to_sc
  !================================
  subroutine mapping_line_sc_to_sp(nDimIn, XyzIn_D, nDimOut, CoordOut_D, &
       IsInterfacePoint)
    integer, intent(in) :: nDimIn
    real,    intent(in) :: XyzIn_D(nDimIn)
    integer, intent(in) :: nDimOut
    real,    intent(out):: CoordOut_D(nDimOut)
    logical, intent(out):: IsInterfacePoint
    
    integer:: iIndex_I(2), iParticle
    !------------------------------------------
    IsInterfacePoint = .true.; iParticle = nint(XyzIn_D(1))
    call SC_get_particle_indexes(iParticle, iIndex_I)
    CoordOut_D = coord_grid_d(SP_Grid, iIndex_I(1),(/iIndex_I(2),1,1/))
  end subroutine mapping_line_sc_to_sp             
  !===========================================
  subroutine SC_get_for_sp_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V

    integer:: iVarBx, iVarBz, iVarMx, iVarMz
    !-----------------------------------------
    ! get buffer with variables
    call SC_get_for_sp(nPartial,iGetStart,Get,w,State_V,nVar)
    ! indices of variables 
    iVarBx = iVar_V(BxCouple_);   iVarBz = iVar_V(BzCouple_)
    iVarMx = iVar_V(RhoUxCouple_);iVarMz = iVar_V(RhoUzCouple_)
    ! perform transformation before returning
    State_V(iVarBx:iVarBz) = matmul(ScToSp_DD,State_V(iVarBx:iVarBz))
    State_V(iVarMx:iVarMz) = matmul(ScToSp_DD,State_V(iVarMx:iVarMz))
  end subroutine SC_get_for_sp_and_transform
  !==========================================================
  subroutine SC_get_line_for_sp_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V
    !------------------------------------------------------------
    ! get buffer with variables
    call SC_get_particle_coords(Get%iCB_II(1,iGetStart),State_V)
    ! perform transformation before returning
    State_V = matmul(ScToSp_DD, State_V)
    State_V = State_V * Grid_C(SC_)%UnitX / Grid_C(SP_)%UnitX
  end subroutine SC_get_line_for_sp_and_transform
  !=============================================== !^CMP END SC
  subroutine couple_ih_sp(DataInputTime)     
    real,intent(in)::DataInputTime
    !----------------------------------------------------------
    if(.not.RouterIhSp%IsProc)return

    tNow=DataInputTime
    if(is_proc(SP_))call SP_put_input_time(DataInputTime)
    
    IhToSp_DD=transform_matrix(tNow,&
         Grid_C(IH_)%TypeCoord, Grid_C(SP_)%TypeCoord)
    call IH_synchronize_refinement(RouterIhSp%iProc0Source,&
         RouterIhSp%iCommUnion)  
    call exchange_data_ih_sp(DoInit=.false., DoExtract=.false.)
  end subroutine couple_ih_sp
  !==================================================================
  subroutine exchange_data_ih_sp(DoInit, DoExtract)
    logical, intent(in) :: DoInit, DoExtract
    !---------------
    if(DoExtract.or..not.DoInit)then
       !\
       !Lagrangian particles are either extracted or 
       !updated in IH
       !/
       ! Send Largrangian particle coords IH => SP
       call construct_router_from_source(&
            GridSource = IH_LocalLineGrid,         &
            GridTarget = SP_Grid,      &
            Router                     = RouterLineIhSp,   &
            n_interface_point_in_block = IH_n_particle,    &
            interface_point_coords=IH_line_interface_point,&
            mapping                    = mapping_line_ih_to_sp)
       call couple_comp(RouterLineIhSp, &
            nVar = 3, &
            fill_buffer = IH_get_line_for_sp_and_transform, &
            apply_buffer= SP_put_line)
    end if
    if(is_proc(SP_).and..not.DoInit)call SP_adjust_lines(&
         DoAdjustStart = .not.use_comp(SC_), DoAdjustEnd = .true.)
    !
    ! Set router SC=> SP to  receive MHD data 
    call set_router(                         &
         GridSource    = IH_Grid,              &
         GridTarget    = SP_LocalGrid,         &
         Router      = RouterIHSp,           &
         n_interface_point_in_block = SP_n_particle, &
         interface_point_coords=             &
         SP_interface_point_coords_for_ih,   &
         mapping     = mapping_sp_to_ih,     &
         interpolate = interpolation_amr_gc)
    call couple_comp(RouterIhSp, &
         nVar = nVarBuffer, &
         fill_buffer = IH_get_for_sp_and_transform, &
         apply_buffer= SP_put_from_ih)
    !By the way get particle coords from the router buffer 
    if(is_proc(IH_).and..not.(DoInit.and.DoExtract))then
       nLength = nlength_buffer_source(RouterIhSp)
       if(use_comp(SC_))call sort_out_sc_particles   !^CMP IF SC
       call IH_add_to_line(&
            nParticle = nLength,&
            Xyz_DI    =  RouterIhSp%BufferSource_II(&
            1:nDim, 1:nLength), &
            nIndex    = nAux,   &
            iIndex_II = nint(RouterIhSp%BufferSource_II(&
            nDim+1:nDim+nAux, 1:nLength)),&      
            UseInputInGenCoord = .true.,&
            DoReplace = .true.)
    end if
  contains                                         !^CMP IF SC BEGIN
    subroutine sort_out_sc_particles
      ! Sort out particles to be advected by the SC
      integer:: iParticle, iParticleNew
      real   :: Xyz_D(3), Coord_D(3)
      !--------------------------------
      iParticleNew = 0
      do iParticle = 1, nLength
         Coord_D = RouterIhSp%BufferSource_II(1:nDim, iParticle)
         call IH_coord_to_xyz(Coord_D, Xyz_D) 
         if(sum(Xyz_D**2) < RScMax**2)CYCLE
         iParticleNew  = iParticleNew + 1
         RouterIhSp%BufferSource_II(:, iParticleNew) = &
              RouterIhSp%BufferSource_II(:, iParticle)
      end do
      nLength = iParticleNew     
    end subroutine sort_out_sc_particles                 !^CMP END SC
  end subroutine exchange_data_ih_sp
  !==================================================================!
  subroutine mapping_sp_to_ih(nDimIn, XyzIn_D, nDimOut, CoordOut_D, &
       IsInterfacePoint)
    integer, intent(in) :: nDimIn
    real,    intent(in) :: XyzIn_D(nDimIn)
    integer, intent(in) :: nDimOut
    real,    intent(out):: CoordOut_D(nDimOut)
    logical, intent(out):: IsInterfacePoint
    real                :: XyzTemp_D(nDim)
    !------------------------------------------
    IsInterfacePoint = .true.; XyzTemp_D = matmul(XyzIn_D, IhToSp_DD)
    XyzTemp_D = XyzTemp_D * Grid_C(SP_)%UnitX / Grid_C(IH_)%UnitX
    call IH_xyz_to_coord(XyzTemp_D, CoordOut_D)
  end subroutine mapping_sp_to_ih
  !==================================================================!
  subroutine mapping_line_ih_to_sp(nDimIn, XyzIn_D, nDimOut, CoordOut_D, &
       IsInterfacePoint)
    integer, intent(in) :: nDimIn
    real,    intent(in) :: XyzIn_D(nDimIn)
    integer, intent(in) :: nDimOut
    real,    intent(out):: CoordOut_D(nDimOut)
    logical, intent(out):: IsInterfacePoint
    
    integer:: iIndex_I(2), iParticle
    !------------------------------------------
    IsInterfacePoint = .true.; iParticle = nint(XyzIn_D(1))
    call IH_get_particle_indexes(iParticle, iIndex_I)
    CoordOut_D = coord_grid_d(SP_Grid, iIndex_I(1), (/iIndex_I(2),1,1/))
  end subroutine mapping_line_ih_to_sp
  !====================================      
  subroutine IH_get_for_sp_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V

    integer:: iVarBx, iVarBz, iVarMx, iVarMz
    !------------------------------------------------------------
    ! get buffer with variables
    call IH_get_for_sp(nPartial,iGetStart,Get,w,State_V,nVar)
    ! indices of variables 
    iVarBx = iVar_V(BxCouple_);   iVarBz = iVar_V(BzCouple_)
    iVarMx = iVar_V(RhoUxCouple_);iVarMz = iVar_V(RhoUzCouple_)
    ! perform transformation before returning
    State_V(iVarBx:iVarBz) = matmul(IhToSp_DD,State_V(iVarBx:iVarBz))
    State_V(iVarMx:iVarMz) = matmul(IhToSp_DD,State_V(iVarMx:iVarMz))
  end subroutine IH_get_for_sp_and_transform
  !==================================================================!        
  subroutine IH_get_line_for_sp_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V
    !------------------------------------------------------------
    ! get buffer with variables
    call IH_get_particle_coords(Get%iCB_II(1, iGetStart), State_V)
    ! perform transformation before returning
    State_V = matmul(IhToSp_DD, State_V)
    State_V = State_V * Grid_C(IH_)%UnitX / Grid_C(SP_)%UnitX
  end subroutine IH_get_line_for_sp_and_transform
end Module CON_couple_mh_sp
