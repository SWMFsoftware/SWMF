!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE SP
!^CMP FILE IH
module CON_couple_mh_sp

  use CON_coupler
  use IH_wrapper, ONLY: IH_synchronize_refinement,  &        
       IH_extract_line, IH_get_for_sp, IH_get_a_line_point, IH_add_to_line,      &
       IH_n_particle, IH_LineDD, IH_line_interface_point, IH_check_use_particles,&
       IH_get_particle_indexes, IH_get_particle_coords                              
 

  use SC_wrapper, ONLY: SC_synchronize_refinement, &      !^CMP IF SC BEGIN 
       SC_extract_line, SC_get_for_sp, SC_get_a_line_point,SC_add_to_line,       & 
       SC_n_particle, SC_LineDD, SC_line_interface_point, SC_check_use_particles,& 
       SC_get_particle_indexes, SC_get_particle_coords    !^CMP END SC
  use CON_axes

  use SP_wrapper, ONLY: &
       SP_put_from_sc, SP_put_from_ih, SP_put_input_time, &
       SP_put_line, SP_n_particle, SP_synchronize_grid, &
       SP_do_extract, SP_get_domain_boundary, SP_put_r_min, &
       SP_interface_point_coords_for_ih, SP_interface_point_coords_for_sc, &
       SP_interface_point_coords_for_ih_extract, SP_get_cell_index, &
       SP_copy_old_state, SP_adjust_lines, SP_get_particle_index, &
       SP_assign_lagrangian_coords
       
  implicit none
  
  private !Except
  public::couple_mh_sp_init
  public::couple_ih_sp              
  public::couple_sc_sp              !^CMP IF SC

  type(GridDescriptorType),save::SP_GridDescriptor !Target (Particle coords)
  type(LocalGDType),       save::SP_LocalGD        !Target (MHD data)
  type(GridDescriptorType),save::IH_GridDescriptor !Source (MHD data)
  type(RouterType),save,private::RouterIhSp        !IH (MHD data) => SP 
  type(GridDescriptorType),save::IH_LineGridDesc   !Misc
  type(LocalGDType),       save::IH_LocalLineGD    !Source (Particle Coords)
  type(RouterType),save,private::RouterLineIhSp    !IH (Particle coords)=>SP    
  !^CMP IF SC BEGIN
  type(GridDescriptorType),save::SC_GridDescriptor !Source (MHD data)  
  type(RouterType),save,private::RouterScSp        !SC MHD data => SP
  type(GridDescriptorType),save::SC_LineGridDesc   !Misc
  type(LocalGDType),       save::SC_LocalLineGD    !Source (Particle Coords)     
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
  ! solar corona and ih upper boundaries
  real:: RSc             !^CMP IF SC
  real:: RIh
contains
  !==================================================================
  subroutine couple_mh_sp_init
    use CON_physics, ONLY: get_time
    use ModConst
    logical,save::DoInit=.true. 
   
    integer:: iError

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
    if(is_proc0(SP_))&
         call SP_do_extract(DoExtract)
    call MPI_Bcast(DoExtract, 1, MPI_LOGICAL, i_proc0(SP_), i_comm(), iError)

    ! get the value of SC and IH boundary as set in SP
    if(is_proc0(SP_))&
         call SP_get_domain_boundary(RSc, RIh)
    call MPI_Bcast(RSc, 1, MPI_REAL, i_proc0(SP_), i_comm(), iError)
    call MPI_Bcast(RIh, 1, MPI_REAL, i_proc0(SP_), i_comm(), iError) 
    ! Set SP side of the routers:
    !\
    ! Set grid descriptors for components
    ! Initialize routers
    !/ 
    ! Set pair SC-SP
    if(use_comp(SC_)) then   !^CMP IF SC BEGIN
       call couple_sc_sp_init
       ! put the lower boundary of the domain in SC to SP
       if(DoExtract.and.is_proc(SP_))&
            call SP_put_r_min(Grid_C(SC_)%Coord1_I(1))
    end if                   !^CMP END SC
    ! Set pair IH-SP         
    if(use_comp(IH_)) then
       call couple_ih_sp_init
       if(DoExtract.and.is_proc(SP_)  &
            .and..not.(use_comp(SC_)) &     !^CMP IF SC
            )call SP_put_r_min(Grid_C(IH_)%Coord1_I(1))
    end if
    if(DoExtract.and.is_proc(SP_))&
         call SP_assign_lagrangian_coords
    DoInit=.false.
  contains
    !===============================
    subroutine couple_sc_sp_init             !^CMP IF SC BEGIN
      call SC_check_use_particles()
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
           GridDescriptorSource = SC_GridDescriptor,&
           GridDescriptorTarget = SP_GridDescriptor,&
           nMappedPointIndex    = nAux,             &
           Router          = RouterScSp        ,    &
           LocalGDTarget   = SP_LocalGD)
      !/
      !\
      ! Set local GD on the Particle_I structure
      call set_standard_grid_descriptor(SC_LineDD,GD=SC_LineGridDesc)
      call init_router(SC_LineGridDesc, SP_GridDescriptor,  &
           RouterLineScSp, nMappedPointIndex=0)
      if(is_proc(SC_))call set_local_gd(iProc = i_proc(),   &
            GD = SC_LineGridDesc, LocalGD = SC_LocalLineGD)
      !Router to send particles is initialized. 
      !Source SC_LocalLineGD, target SP_GridDescriptor

      if(.not.RouterScSp%IsProc)RETURN
      call SC_synchronize_refinement(RouterScSp%iProc0Source,&
           RouterScSp%iCommUnion)
      ScToSp_DD=transform_matrix(tNow,&                 
           Grid_C(SC_)%TypeCoord, Grid_C(SP_)%TypeCoord)
      if(DoExtract)then
         !Router for sending MHD data from SC to SP may be also
         !used, at the stage of its construction, to send the
         !particle coordinates backward. Use only a second option
         !now, therefore, we cal only two of three subroutines constituting
         !the set_router routine
         if(is_proc(SP_))&
              call set_semi_router_from_target(&
              GDSource              = SC_GridDescriptor, &
              GDTarget              = SP_LocalGD, &
              Router                = RouterScSp, &
              n_interface_point_in_block = SP_n_particle,&
              interface_point_coords= SP_interface_point_coords_for_sc, &
              mapping               = mapping_sp_to_sc, &
              interpolate           = interpolation_amr_gc)
         call synchronize_router_target_to_source(RouterScSp)
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
                 RSoftBoundaryIn   =  RSc, & 
                 UseInputInGenCoord= .true.)
         end if
      end if
      !First coupling with the particle info and data exchange
      call exchange_data_sc_sp(DoInit, DoExtract)
    end subroutine couple_sc_sp_init                 !^CMP END SC
    !==============================
    subroutine couple_ih_sp_init                     
      call IH_check_use_particles()
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
           GridDescriptorSource = IH_GridDescriptor,&
           GridDescriptorTarget = SP_GridDescriptor,&
           nMappedPointIndex    = nAux,             &
           Router          = RouterIhSp)
      if(.not.use_comp(SC_))then !^CMP IF SC BEGIN
         if(is_proc(SP_))call set_local_gd(i_proc(),&
              SP_GridDescriptor, SP_localGD)
      end if                     !^CMP END SC
      !/
      !\
      ! Set local GD on the Particle_I structure
      call set_standard_grid_descriptor(IH_LineDD,GD=&
           IH_LineGridDesc)
      call init_router(IH_LineGridDesc, SP_GridDescriptor, RouterLineIhSp,&
           nMappedPointIndex=0)
      if(is_proc(IH_))call set_local_gd(&
           iProc = i_proc(), &
           GD = IH_LineGridDesc, &
           LocalGD = IH_LocalLineGD)
      !Router to send particles is initialized. 
      !Source IH_LocalLineGD, target SP_GridDescriptor
      if(.not.RouterIhSp%IsProc)RETURN
      call IH_synchronize_refinement(RouterIhSp%iProc0Source,&
           RouterIhSp%iCommUnion)
      IhToSp_DD=transform_matrix(tNow,&                 
           Grid_C(IH_)%TypeCoord, Grid_C(SP_)%TypeCoord)
      if(DoExtract)then
         !Router for sending MHD data from IH to SP may be also
         !used, at the stage of its construction, to send the
         !particle coordinates backward. Use only a second option
         !now, therefore, we cal only two of three subroutines constituting
         !the set_router routine
         if(is_proc(SP_))&
              call set_semi_router_from_target(&
              GDSource  = IH_GridDescriptor, &
              GDTarget  = SP_LocalGD, &
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
                 nLine             = nLength, &
                 XyzOrigin_DI      = &
                 RouterIhSp%BufferSource_II(1:nDim,1:nLength), &
                 iTraceMode        = iInterfaceEnd, &
                 nIndex            = nAux, &
                 iIndexOrigin_II   = nint(RouterIhSp%&
                 BufferSource_II(nDim+1:nDim+nAux,1:nLength)),&
                 RSoftBoundaryIn   =  RIh, & 
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
    call SP_synchronize_grid(RouterScSp%iCommUnion)
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
       !\
       ! Send Largrangian particle coords from SC to SP
       call construct_router_from_source(&
            GDSource = SC_LocalLineGD,     &
            GDTarget = SP_GridDescriptor,  &
            Router               = RouterLineScSp,     &
            n_interface_point_in_block = SC_n_particle,&
            interface_point_coords=SC_line_interface_point,&
            mapping              = mapping_line_sc_to_sp)
       call couple_comp(RouterLineScSp, &
            nVar = 3, &
            fill_buffer = SC_get_line_for_sp_and_transform, &
            apply_buffer= SP_put_line)
       !Lagrangian particle coordinates are sent
       !/
    end if
    if(is_proc(SP_).and..not.DoInit)&
         call SP_adjust_lines(SC_)
    call SP_synchronize_grid(RouterScSp%iCommUnion)
    !\
    ! Send coordinates of the points in SP to receive MHD data from SC
    call set_router(&
         GDSource  = SC_GridDescriptor,             &
         GDTarget  = SP_LocalGD,                    &
         Router                = RouterScSp,        &
         n_interface_point_in_block = SP_n_particle,&
         interface_point_coords=                    &
         SP_interface_point_coords_for_sc,          &
         mapping               = mapping_sp_to_sc,  &
         interpolate           = interpolation_amr_gc,&
         extra_data            = fix_buffer)
    call couple_comp(RouterScSp, &
         nVar = nVarBuffer, &
         fill_buffer = SC_get_for_sp_and_transform, &
         apply_buffer= SP_put_from_sc)
    !MHD Data from SC to SP are sent
    !/
    !Get particles from the semi-router 
    if(is_proc(SC_))then!.and..not.(DoInit.and.DoExtract))then
       nLength = nlength_buffer_source(RouterScSp)
       call SC_add_to_line(&
            nParticle = nLength,&
            Xyz_DI    =  RouterScSp%BufferSource_II(&
            1:nDim, 1:nLength), &
            nIndex    = nAux,   &
            iIndex_II = nint(RouterScSp%BufferSource_II(&
            nDim+1:nDim+nAux, 1:nLength)),&      
            UseInputInGenCoord = .true.,&
            DoReplace = .true.)
    end if
  end subroutine exchange_data_sc_sp
  !^CMP END SC
  !================================
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
    call SP_synchronize_grid(RouterIhSp%iCommUnion)  
    call exchange_data_ih_sp(DoInit=.false., DoExtract=.false.)
  end subroutine couple_ih_sp
  !==================================================================
  subroutine exchange_data_ih_sp(DoInit, DoExtract)
    logical, intent(in) :: DoInit, DoExtract
    integer:: iParticle, iParticleNew
    !---------------
    if(DoExtract.or..not.DoInit)then
       !\
       !Lagrangian particles are either extracted or 
       !updated in IH
       !/
       !\
       ! Send coordinates of the constructed or advected 
       ! Lagrangian points from IH to SP
       call construct_router_from_source(&
            GDSource = IH_LocalLineGD,         &
            GDTarget = SP_GridDescriptor,      &
            Router                     = RouterLineIhSp,   &
            n_interface_point_in_block = IH_n_particle,    &
            interface_point_coords=IH_line_interface_point,&
            mapping                    = mapping_line_ih_to_sp)
       call couple_comp(RouterLineIhSp, &
            nVar = 3, &
            fill_buffer = IH_get_line_for_sp_and_transform, &
            apply_buffer= SP_put_line)
       !Particle coordinates are sent to SP
       !/
    end if
    if(is_proc(SP_).and..not.DoInit)&
         call SP_adjust_lines(IH_)
    call SP_synchronize_grid(RouterIhSp%iCommUnion)
    !\
    ! Send to IH the  coordinates of points in SP and
    ! send back the MHD data in these points
    call set_router(                         &
         GDSource    = IH_GridDescriptor,    &
         GDTarget    = SP_LocalGD,           &
         Router      = RouterIHSp,           &
         n_interface_point_in_block =        &
         SP_n_particle,                      &
         interface_point_coords=             &
         SP_interface_point_coords_for_ih,   &
         mapping     = mapping_sp_to_ih,     &
         interpolate = interpolation_amr_gc, &
         extra_data            = fix_buffer)
    call couple_comp(RouterIhSp, &
         nVar = nVarBuffer, &
         fill_buffer = IH_get_for_sp_and_transform, &
         apply_buffer= SP_put_from_ih)
    !MHD data are sent
    !/
    !Get particles from the semi-router 
    if(is_proc(IH_).and..not.(DoInit.and.DoExtract))then
       nLength = nlength_buffer_source(RouterIhSp)
       !^CMPIF SC BEGIN
       ! Sort out particles advected by the SC
       !/
       iParticleNew = 0
       do iParticle = 1, nLength
          if(sum(RouterIhSp%BufferSource_II(&
            1:nDim, iParticle)**2) < RSc**2)CYCLE
          iParticleNew  = iParticleNew +1
          RouterIhSp%BufferSource_II(:, iParticleNew) = &
               RouterIhSp%BufferSource_II(:, iParticle)
       end do
       nLength = iParticleNew
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
  end subroutine exchange_data_ih_sp
  !==================================================================!
  subroutine xyz_to_coord(XyzIn_D, CoordOut_D, TypeGeometry)
    ! mapping from generalized coordinates in CompIn 
    ! to generalized coordinates in CompOut
    use ModCoordTransform, ONLY: xyz_to_rlonlat, rlonlat_to_xyz
    real,    intent(in) :: XyzIn_D(nDim)
    real,    intent(out):: CoordOut_D(nDim)
    character(len=*), intent(in):: TypeGeometry
    character(len=*), parameter:: NameSub='CON_couple_mh_sp:xyz_to_coord'
    !------------------------------------------------------------
    ! convert to geometry output type
    if( index(TypeGeometry, 'spherical_lnr') > 0 )then
       call xyz_to_rlonlat(XyzIn_D, CoordOut_D)
       ! convert radius to log(radius)
       if(CoordOut_D(1) > 0) CoordOut_D(1) = log(CoordOut_D(1))
    elseif( index(TypeGeometry, 'spherical') > 0 )then
       call xyz_to_rlonlat(XyzIn_D, CoordOut_D)
    elseif(index(TypeGeometry, 'cartesian') > 0 )then
       CoordOut_D = XyzIn_D
    else
       call CON_stop(NameSub//&
            ': unknown type of geometry '//trim(TypeGeometry))
    end if
  end subroutine xyz_to_coord
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
    call xyz_to_coord(XyzIn_D = XyzTemp_D, CoordOut_D = CoordOut_D, &
         TypeGeometry = Grid_C(IH_)%TypeGeometry)
  end subroutine mapping_sp_to_ih
  !==================================================================!
  subroutine mapping_line_ih_to_sp(nDimIn, XyzIn_D, nDimOut, CoordOut_D, &
       IsInterfacePoint)
    use CON_grid_descriptor
    integer, intent(in) :: nDimIn
    real,    intent(in) :: XyzIn_D(nDimIn)
    integer, intent(in) :: nDimOut
    real,    intent(out):: CoordOut_D(nDimOut)
    logical, intent(out):: IsInterfacePoint
    
    integer:: iIndex_I(2), iParticle, iCell
    !------------------------------------------
    IsInterfacePoint = .true.; 
    iParticle = nint(XyzIn_D(1))
    call IH_get_particle_indexes(iParticle, iIndex_I)
    call SP_get_cell_index(iIndex_I(1), iIndex_I(2), iCell)
    CoordOut_D = xyz_grid_d(SP_GridDescriptor,iIndex_I(1),(/iCell,1,1/))
  end subroutine mapping_line_ih_to_sp
  !=======================================!^CMP IF SC BEGIN
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
    call xyz_to_coord(XyzIn_D = XyzTemp_D, CoordOut_D = CoordOut_D, &
         TypeGeometry = Grid_C(SC_)%TypeGeometry)
  end subroutine mapping_sp_to_sc
  !================================
  subroutine mapping_line_sc_to_sp(nDimIn, XyzIn_D, nDimOut, CoordOut_D, &
       IsInterfacePoint)
    use CON_grid_descriptor
    integer, intent(in) :: nDimIn
    real,    intent(in) :: XyzIn_D(nDimIn)
    integer, intent(in) :: nDimOut
    real,    intent(out):: CoordOut_D(nDimOut)
    logical, intent(out):: IsInterfacePoint
    
    integer:: iIndex_I(2), iParticle, iCell
    !------------------------------------------
    IsInterfacePoint = .true.
    iParticle = nint(XyzIn_D(1))
    call SC_get_particle_indexes(iParticle, iIndex_I)
    call SP_get_cell_index(iIndex_I(1), iIndex_I(2), iCell)
    CoordOut_D = xyz_grid_d(SP_GridDescriptor,iIndex_I(1),(/iCell,1,1/))
  end subroutine mapping_line_sc_to_sp
  !===================================!^CMP END SC        
  subroutine IH_get_for_sp_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)

    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V
    integer:: iVarBx, iVarBz
    !------------------------------------------------------------
    ! get buffer with variables
    call IH_get_for_sp(&
         nPartial,iGetStart,Get,w,State_V,nVar)
    ! indices of variables 
    iVarBx = iVar_V(BxCouple_)
    iVarBz = iVar_V(BzCouple_)
    ! perform transformation before returning
    State_V(iVarBx:iVarBz)=matmul(IhToSp_DD,State_V(iVarBx:iVarBz))
  end subroutine IH_get_for_sp_and_transform
  !==================================================================!        
  subroutine IH_get_line_for_sp_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)
    use ModCoordTransform, ONLY: xyz_to_rlonlat
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V
    !------------------------------------------------------------
    ! get buffer with variables
    call IH_get_particle_coords(Get%iCB_II(1, iGetStart), State_V)
    ! perform transformation before returning
    State_V = matmul(IhToSp_DD, State_V)
  end subroutine IH_get_line_for_sp_and_transform
  !==============================================================!
  subroutine fix_buffer(BlockIndex_I,&
            iGridPointInBlock, nAux, Data_V)
    use CON_coupler, ONLY: GlobalTreeNode_, BLK_, &
         GridPointFirst_, GlobalBlock_
    !\
    ! Mapped point in the buffer is characterized
    ! by any of the block indexes: global tree node number,
    ! global block number, local block number or the number of
    ! the first point in the block...
    integer, intent(in) :: BlockIndex_I(1:4)
    ! and the point number in the block
    integer, intent(in) :: iGridPointInBlock
    integer, intent(in) :: nAux !How many reals you may send
    real,    intent(out):: Data_V(1:nAux)
   
    integer:: iIndexIn, iLine, iIndexOut
    !-------------------------------------------------------------
    iLine    = BlockIndex_I(GlobalBlock_)
    iIndexIn = iGridPointInBlock
    call SP_get_particle_index(iLine, iIndexIn, iIndexOut)
    Data_V(1) = real(iLine); Data_V(2) = real(iIndexOut)
  end subroutine fix_buffer
  !^CMP IF SC BEGIN    ===========================================
  subroutine SC_get_for_sp_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V
    integer:: iVarBx, iVarBz
    !------------------------------------------------------------
    ! get buffer with variables
    call SC_get_for_sp(nPartial,iGetStart,Get,w,State_V,nVar)
    ! indices of variables 
    iVarBx = iVar_V(BxCouple_); iVarBz = iVar_V(BzCouple_)
    ! perform transformation before returning
    State_V(iVarBx:iVarBz) = matmul(ScToSp_DD,State_V(iVarBx:iVarBz))
  end subroutine SC_get_for_sp_and_transform
  !=============================================================
  subroutine SC_get_line_for_sp_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)
    use ModCoordTransform, ONLY: xyz_to_rlonlat
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V
    !------------------------------------------------------------
    ! get buffer with variables
    call SC_get_particle_coords(Get%iCB_II(1,iGetStart),State_V)
    ! perform transformation before returning
    State_V = matmul(ScToSp_DD, State_V)
  end subroutine SC_get_line_for_sp_and_transform
  !=============================================== !^CMP END SC
end Module CON_couple_mh_sp
