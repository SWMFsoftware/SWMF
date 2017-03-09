!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE SP

module CON_couple_mh_sp

  ! This coupler employs the following global arrays.
  ! SP_Xyz_DI - is the array of the Lagrangian points.
  ! A part of this array, as well as the mask 'SP_IsInIH' is only availbale
  ! at the PE set, at which either IH or SP run. 
  ! This coordinates are expressed in terms of
  ! the length units and with respect to the frame of reference defined in IH.
  ! Another part of this array, as well as the mask 'SP_IsInSC' 
  ! is only available at the PE set, at which either SC or SP run. 
  ! This coordinates are expressed in terms of
  ! the length units and with respect to the frame of reference defined in SC.
  ! SP_XyzSP - in the array of all Lagrangian points, in units and in the 
  ! frame of reference defined at SP, is available only at the PE set, 
  ! at which SP runs.

  use CON_coupler

  use IH_wrapper, ONLY: IH_synchronize_refinement, &        !^CMP IF IH
       IH_extract_line, IH_get_for_sp, IH_get_a_line_point,&!^CMP IF IH
       IH_get_scatter_line, IH_add_to_line                  !^CMP IF IH

  use SC_wrapper, ONLY: SC_synchronize_refinement, &        !^CMP IF SC
       SC_extract_line, SC_get_for_sp, SC_get_a_line_point,&!^CMP IF SC
       SC_get_scatter_line                                  !^CMP IF SC

  use CON_global_message_pass
  use CON_axes

  use SP_wrapper, ONLY: &
       SP_put_from_mh, SP_put_input_time, &
       SP_put_line, SP_get_request, SP_get_grid_descriptor_param, &
       SP_get_line_all, SP_get_solar_corona_boundary, SP_put_r_min

  implicit none

  private !Except
  public::couple_mh_sp_init
  public::couple_ih_sp              !^CMP IF IH
  public::couple_sc_sp              !^CMP IF SC

  type(GridDescriptorType),save::SP_GridDescriptor !Target

  type(GridDescriptorType),save::IH_GridDescriptor !Source  !^CMP IF IH
  type(RouterType),save,private::RouterIhSp                 !^CMP IF IH

  type(GridDescriptorType),save::SC_GridDescriptor !Source  !^CMP IF SC
  type(RouterType),save,private::RouterScSp                 !^CMP IF SC

  logical,save::DoInit=.true.
  real,   allocatable:: XyzStored_DI(  :,:)
  integer,allocatable:: iAuxStored_II(:,:)

  integer::iError
  real,save::rBoundIh=21.0                !^CMP IF IH
  real,save::rBoundSc=1.20                !^CMP IF SC

  ! available directions of interface between SP and MH 
  integer, parameter:: &
       iInterfaceBegin = -1, iInterfaceOrigin = 0, iInterfaceEnd = 1
  
  logical::DoTest,DoTestMe
  character(LEN=*),parameter::NameSub='couple_mh_sp'
  real,dimension(3,3)::ScToIh_DD,ScToSp_DD,IhToSp_DD,SpToSc_DD
  real :: tNow

contains
  !==================================================================
  subroutine couple_mh_sp_init
    use CON_physics, ONLY: get_time
    use ModConst

    logical::DoneRestart

    integer:: iGridMin_D(3), iGridMax_D(3), ierror
    real:: Disp_D(3)
    real, pointer:: CoordMisc_DI(:,:)

    ! available directions of interface between SP and MH 
    ! (see subroutine exchange_lines below)
    integer, parameter:: &
         iInterfaceBegin = -1, iInterfaceOrigin = 0, iInterfaceEnd = 1

    ! solar corona boundary
    real:: RSc, RMin

    character(len=*), parameter:: NameSub = 'couple_mh_sp_init'
    !----------------------------------------------------------------------
    if(.not.DoInit)return
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    DoInit=.false.
    !The initialization can be done only once

    call get_time(tSimulationOut=tNow)

    !\
    ! Set grid descriptors for components
    ! Initialize routers
    !/
    call SP_get_grid_descriptor_param(iGridMin_D, iGridMax_D, Disp_D)
    call set_grid_descriptor_id(SP_,&
         nDim = 3, &
         iGridPointMin_D = iGridMin_D, &
         iGridPointMax_D = iGridMax_D, &
         Displacement_D  = Disp_D, &
         GridDescriptor  = SP_GridDescriptor)

    ! get the value of solar corona boundary as set in SP
    call SP_get_solar_corona_boundary(RSc)

    if(use_comp(SC_))then  
       ! Set pair SC-SP
       call set_standard_grid_descriptor(SC_,GridDescriptor=&
            SC_GridDescriptor)
       call init_router(SC_GridDescriptor,SP_GridDescriptor,&
            RouterScSp)
       call SC_synchronize_refinement(RouterScSp%iProc0Source,&
            RouterScSp%iCommUnion)
       ScToSp_DD=transform_matrix(tNow,&                 
            Grid_C(SC_)%TypeCoord, Grid_C(SP_)%TypeCoord)
       call set_couple_var_info(SC_, SP_)
       call exchange_lines(SC_)
       ! get the lower boundary of the domain in SC
       ! and put this info to SP
       call SP_put_r_min(Grid_C(SC_)%Coord1_I(1))
    end if

    if(use_comp(IH_))then
       ! Set pair IH-SP
       call set_standard_grid_descriptor(IH_,GridDescriptor=&
            IH_GridDescriptor)
       call init_router(IH_GridDescriptor,SP_GridDescriptor,&
            RouterIhSp)
       call IH_synchronize_refinement(RouterIhSp%iProc0Source,&
            RouterIhSp%iCommUnion)
       IhToSp_DD=transform_matrix(tNow,&                 
            Grid_C(IH_)%TypeCoord, Grid_C(SP_)%TypeCoord)
       call set_couple_var_info(IH_, SP_)
       call exchange_lines(IH_)
    end if

  contains
    !================================================================    
    subroutine exchange_lines(iMHComp)
      ! MH extracts and sends field lines requested by SP;
      !----------------------------------------------------------------
      ! index of MH component
      integer, intent(in):: iMHComp
      !----------------------------------------------------------------

      ! conversion matrix between SP and MH coordinates
      real:: Convert_DD(3,3)
      ! request coordinates, one per line
      real, pointer:: CoordMisc_DI(:,:)
      ! number of particles per line
      integer,allocatable:: nParticleAtLine_I(:)
      ! requested variables
      integer:: nVar
      character(len=100):: NameVar
      ! particle data
      real, allocatable:: Particle_II(:,:)
      ! field line index in particle data
      integer:: iFLIndex
      ! MPI
      integer:: MH_iProcFrom, SP_iProcTo
      integer:: iProcTo_I(1), iProcFrom_I(1)
      integer:: nParticleThisProc
      integer, allocatable:: nParticleRecv_I(:), nParticleSend_I(:)
      integer, allocatable:: iStatus_II(:,:), iRequestS_I(:), iRequestR_I(:)
      integer:: nRequestS, nRequestR
      integer:: iTag = 0
      real, allocatable:: BuffRecv_I(:), BuffSend_I(:)
      ! loop variables
      integer:: iLine, iBuff, iParticle
      !----------------------------------------------------------------
      select case(iMHComp)
      case(SC_)
         call set_router_from_target_2_stage(&
              GridDescriptorSource = SC_GridDescriptor, &
              GridDescriptorTarget = SP_GridDescriptor, &
              Router               = RouterScSp, &
              get_request_target   = SP_get_request, &
              transform            = transform_sp_to_sc, &
              interpolate_source   = interpolation_amr_gc, &
              put_request_source   = SC_put_request)
         if(is_proc(SC_))&
              call SC_extract_line(&
              ubound(XyzStored_DI,2),  XyzStored_DI, iInterfaceOrigin,&
              ubound(iAuxStored_II,1),iAuxStored_II, RSc)
         call set_router_from_source_2_stage(&
              GridDescriptorSource = SC_GridDescriptor, &
              GridDescriptorTarget = SP_GridDescriptor, &
              Router               = RouterScSp, &
              get_scatter_source   = SC_get_scatter_line, &
              transform            = transform_sc_to_sp, &
              interpolate_source   = interpolation_amr_gc, &
              interpolate_target   = interpolate_sp, &
              put_scatter_target   = SP_put_scatter_from_mh)
         call global_message_pass(RouterScSp, &
              nVar = nVarBuffer, &
              fill_buffer = SC_get_for_sp_and_transform, &
              apply_buffer= SP_put_from_mh)
      case(IH_)
         call set_router_from_target_2_stage(&
              GridDescriptorSource = IH_GridDescriptor, &
              GridDescriptorTarget = SP_GridDescriptor, &
              Router               = RouterIhSp, &
              get_request_target   = SP_get_request, &
              transform            = transform_sp_to_ih, &
              interpolate_source   = interpolation_amr_gc, &
              put_request_source   = IH_put_request)
         if(is_proc(IH_))&
              call IH_extract_line(&
              ubound(XyzStored_DI,2),  XyzStored_DI, iInterfaceEnd,&
              ubound(iAuxStored_II,1),iAuxStored_II)
         call set_router_from_source_2_stage(&
              GridDescriptorSource = IH_GridDescriptor, &
              GridDescriptorTarget = SP_GridDescriptor, &
              Router               = RouterIhSp, &
              get_scatter_source   = IH_get_scatter_line, &
              transform            = transform_ih_to_sp, &
              interpolate_source   = interpolation_amr_gc, &
              interpolate_target   = interpolate_sp, &
              put_scatter_target   = SP_put_scatter_from_mh)
         call global_message_pass(RouterIhSp, &
              nVar = nVarBuffer, &
              fill_buffer = IH_get_for_sp_and_transform, &
              apply_buffer= SP_put_from_mh)
      end select
    end subroutine exchange_lines

  end subroutine couple_mh_sp_init
  !==================================================================!

  subroutine transform(iCompIn, iCompOut, nDimIn, CoordIn_D, nDimOut, CoordOut_D)
    ! transform from generalized coordinates in CompIn 
    ! to generalized coordinates in CompOut
    use ModCoordTransform, ONLY: xyz_to_rlonlat, rlonlat_to_xyz
    integer, intent(in) :: iCompIn
    integer, intent(in) :: iCompOut
    integer, intent(in) :: nDimIn
    real,    intent(in) :: CoordIn_D(nDimIn)
    integer, intent(in) :: nDimOut
    real,    intent(out):: CoordOut_D(nDimOut)

    character(len=100):: TypeGeometryIn, TypeGeometryOut
    real:: XyzTemp_D(nDimOut), CoordTemp_D(nDimIn)
    real:: Convert_DD(nDimIn, nDimOut)
    character(len=*), parameter:: NameSub='CON_couple_mh_sp:transform'
    !------------------------------------------------------------
    if(nDimIn /= 3 .or. nDimOut /= 3)&
         call CON_stop(NameSub//': MH or SP component is not 3D')
    ! convert from geometry input type
    TypeGeometryIn = Grid_C(iCompIn)%TypeGeometry
    if(index(TypeGeometryIn, 'spherical_lnr') > 0 )then
       ! convert log(radius) to radius first
       CoordTemp_D = CoordIn_D
       CoordTemp_D(1) = exp(CoordTemp_D(1))
       call rlonlat_to_xyz(CoordTemp_D, XyzTemp_D)
    elseif( index(TypeGeometryIn, 'spherical') > 0 )then
       call rlonlat_to_xyz(CoordIn_D, XyzTemp_D)
    elseif(index(TypeGeometryIn, 'cartesian') > 0 )then
       XyzTemp_D = CoordIn_D
    else
       call CON_stop(NameSub//&
            ': unkown type of geometry '//trim(TypeGeometryIn))
    end if

    ! matrix for cartesian transform  between components
    Convert_DD = transform_matrix(tNow, &
         Grid_C(iCompIn)%TypeCoord,  Grid_C(iCompOut)%TypeCoord)

    ! geometry of iCompOut component
    TypeGeometryOut = Grid_C(iCompOut)%TypeGeometry
    ! rotate cartesian
    XyzTemp_D = matmul(Convert_DD, XyzTemp_D)
    ! convert to geometry output type
    if( index(TypeGeometryOut, 'spherical_lnr') > 0 )then
       call xyz_to_rlonlat(XyzTemp_D, CoordOut_D)
       ! convert radius to log(radius)
       if(CoordOut_D(1) > 0) CoordOut_D(1) = log(CoordOut_D(1))
    elseif( index(TypeGeometryOut, 'spherical') > 0 )then
       call xyz_to_rlonlat(XyzTemp_D, CoordOut_D)
    elseif(index(TypeGeometryOut, 'cartesian') > 0 )then
       CoordOut_D = XyzTemp_D
    else
       call CON_stop(NameSub//&
            ': unkown type of geometry '//trim(TypeGeometryOut))
    end if
  end subroutine transform

  !==================================================================!
  subroutine transform_sp_to_sc(nDimIn, XyzIn_D, nDimOut, CoordOut_D)
    integer, intent(in) :: nDimIn
    real,    intent(in) :: XyzIn_D(nDimIn)
    integer, intent(in) :: nDimOut
    real,    intent(out):: CoordOut_D(nDimOut)
    !------------------------------------------
    call transform(SP_, SC_,nDimIn, XyzIn_D, nDimOut, CoordOut_D)
  end subroutine transform_sp_to_sc
  !==================================================================!
  subroutine transform_sp_to_ih(nDimIn, XyzIn_D, nDimOut, CoordOut_D)
    integer, intent(in) :: nDimIn
    real,    intent(in) :: XyzIn_D(nDimIn)
    integer, intent(in) :: nDimOut
    real,    intent(out):: CoordOut_D(nDimOut)
    !------------------------------------------
    call transform(SP_,IH_,nDimIn, XyzIn_D, nDimOut, CoordOut_D)
  end subroutine transform_sp_to_ih
  !==================================================================!
  subroutine transform_sc_to_sp(nDimIn, XyzIn_D, nDimOut, CoordOut_D)
    integer, intent(in) :: nDimIn
    real,    intent(in) :: XyzIn_D(nDimIn)
    integer, intent(in) :: nDimOut
    real,    intent(out):: CoordOut_D(nDimOut)
    !------------------------------------------
    call transform(SC_, SP_,nDimIn, XyzIn_D, nDimOut, CoordOut_D)
  end subroutine transform_sc_to_sp
  !==================================================================!
  subroutine transform_ih_to_sp(nDimIn, XyzIn_D, nDimOut, CoordOut_D)
    integer, intent(in) :: nDimIn
    real,    intent(in) :: XyzIn_D(nDimIn)
    integer, intent(in) :: nDimOut
    real,    intent(out):: CoordOut_D(nDimOut)
    !------------------------------------------
    call transform(IH_,SP_,nDimIn, XyzIn_D, nDimOut, CoordOut_D)
  end subroutine transform_ih_to_sp
  !==================================================================!
  subroutine put_request(iComp, nData, &
       nDim, Coord_DI, nIndex, iIndex_II, nAux, Aux_VI)
    use ModCoordTransform, ONLY: rlonlat_to_xyz
    integer, intent(in):: iComp
    integer, intent(in):: nData
    integer, intent(in):: nDim
    real,    intent(in):: Coord_DI(nDim, nData)
    integer, intent(in):: nIndex
    integer, intent(in):: iIndex_II(nIndex, nData)
    integer, intent(in):: nAux
    real,    intent(in):: Aux_VI(nAux,nData)
    integer:: iData
    real:: Xyz_D(nDim)
    character(len=100):: TypeGeometry
    character(len=*),parameter::NameSub='CON_couple_mh_sp:put_request'
    !----------------------------------------------------------
    if(allocated(XyzStored_DI)) deallocate(XyzStored_DI)
    allocate(XyzStored_DI(nDim, nData))
    if(allocated(iAuxStored_II)) deallocate(iAuxStored_II)
    allocate(iAuxStored_II(nAux,nData))

    if(nData==0)&
         RETURN

    ! coordinates of requested field lines' origins
    XyzStored_DI = Coord_DI
    ! indices of paritcles along the field lines
    iAuxStored_II = nint(Aux_VI)
    ! perform transformations based on the type of geometry
    TypeGeometry = Grid_C(iComp)%TypeGeometry
    if( index(TypeGeometry, 'spherical_lnr') > 0 )then
       ! convert to radius from log(radius)
       XyzStored_DI(1,:) = exp(XyzStored_DI(1,:))
       do iData = 1, nData
          call rlonlat_to_xyz(XyzStored_DI(:, iData), Xyz_D)
          XyzStored_DI(:, iData) = Xyz_D
       end do
    elseif(index(TypeGeometry, 'cartesian') > 0 )then
       ! do nothing
    else
       call CON_stop(NameSub//': unkown type of geometry '//trim(TypeGeometry))
    end if
  end subroutine put_request
  !==================================================================!
  subroutine SC_put_request(nData, &
       nDim, Coord_DI, nIndex, iIndex_II, nAux, Aux_VI)
    integer, intent(in):: nData
    integer, intent(in):: nDim
    real,    intent(in):: Coord_DI(nDim, nData)
    integer, intent(in):: nIndex
    integer, intent(in):: iIndex_II(nIndex, nData)
    integer, intent(in):: nAux
    real,    intent(in):: Aux_VI(nAux, nData)
    !----------------------------------------------------------
    call put_request(SC_, nData, nDim, Coord_DI, nIndex, iIndex_II, &
         nAux, Aux_VI)
  end subroutine SC_put_request
  !==================================================================!
  subroutine IH_put_request(nData, &
       nDim, Coord_DI, nIndex, iIndex_II, nAux, Aux_VI)
    integer, intent(in):: nData
    integer, intent(in):: nDim
    real,    intent(in):: Coord_DI(nDim, nData)
    integer, intent(in):: nIndex
    integer, intent(in):: iIndex_II(nIndex, nData)
    integer, intent(in):: nAux
    real,    intent(in):: Aux_VI(nAux, nData)
    !----------------------------------------------------------
    call put_request(IH_, nData, nDim, Coord_DI, nIndex, iIndex_II, &
         nAux, Aux_VI)
  end subroutine IH_put_request
  !==================================================================!
  subroutine SP_put_scatter_from_mh(nData, nDim, Coord_DI, nIndex, iIndex_II)
    integer, intent(in):: nDim
    integer, intent(in):: nData
    real,    intent(in):: Coord_DI(nDim, nData)
    integer, intent(in):: nIndex
    integer, intent(in):: iIndex_II(nIndex, nData)
    !--------------------------------------------------------
    call SP_put_line(nData, Coord_DI, iIndex_II)
  end subroutine SP_put_scatter_from_mh
  !==================================================================!
  subroutine interpolate_sp(&
       nCoord, Coord_I, GridDescriptor, &
       nIndex, iIndex_II, nImage, Weight_I)
    use CON_grid_descriptor
    ! number of indices per data entry
    integer, intent(in):: nCoord
    ! data location on Source
    real,    intent(inout):: Coord_I(nCoord)
    ! grid descriptor
    type(GridDescriptorType):: GridDescriptor
    integer, intent(in) :: nIndex
    integer, intent(out):: iIndex_II(0:nIndex,2**GridDescriptor%nDim)
    integer, intent(out):: nImage
    real,    intent(out):: Weight_I(2**GridDescriptor%nDim)
    integer:: iLine
    !--------------------------
    nImage = 1
    Weight_I(1) = 1.0
    iLine = nint(Coord_I(4))

    iIndex_II(0,  1) = GridDescriptor%DD%Ptr%iDecomposition_II(PE_,iLine)
    iIndex_II(1,  1) = nint(Coord_I(5))
    iIndex_II(2:3,1) = (/1,1/)
    iIndex_II(4,  1) = GridDescriptor%DD%Ptr%iDecomposition_II(BLK_,iLine)
  end subroutine interpolate_sp
  !==================================================================
  !^CMP IF IH BEGIN
  subroutine couple_ih_sp(DataInputTime)     

    use CON_global_message_pass

    real,intent(in)::DataInputTime
    real,dimension(3)::Xyz_D
    !-------------------------------------------------------------------------
    if(.not.RouterIhSp%IsProc)return

    tNow=DataInputTime
    IhToSp_DD=transform_matrix(tNow,&
         Grid_C(IH_)%TypeCoord, Grid_C(SP_)%TypeCoord)
    ScToIh_DD=transform_matrix(tNow,&                   !^CMP IF SC
         Grid_C(SC_)%TypeCoord, Grid_C(IH_)%TypeCoord)  !^CMP IF SC


    call IH_synchronize_refinement(RouterIhSp%iProc0Source,&
         RouterIhSp%iCommUnion)
    call set_router_from_target_2_stage(&
         GridDescriptorSource = IH_GridDescriptor, &
         GridDescriptorTarget = SP_GridDescriptor, &
         Router               = RouterIhSp, &
         get_request_target   = SP_get_request, &
         transform            = transform_sp_to_ih, &
         interpolate_source   = interpolation_amr_gc, &
         put_request_source   = IH_put_request)
    if(is_proc(IH_))&
         call IH_add_to_line(&
         ubound(XyzStored_DI,2),   XyzStored_DI, &
         ubound(iAuxStored_II,1), iAuxStored_II)
    call set_router_from_source_2_stage(&
         GridDescriptorSource = IH_GridDescriptor, &
         GridDescriptorTarget = SP_GridDescriptor, &
         Router               = RouterIhSp, &
         get_scatter_source   = IH_get_scatter_line, &
         transform            = transform_ih_to_sp, &
         interpolate_source   = interpolation_amr_gc, &
         interpolate_target   = interpolate_sp, &
         put_scatter_target   = SP_put_scatter_from_mh)
    call global_message_pass(RouterIhSp, &
         nVar = nVarBuffer, &
         fill_buffer = IH_get_for_sp_and_transform, &
         apply_buffer= SP_put_from_mh)
!    if(is_proc(SP_))then
!       call SP_put_input_time(DataInputTime)
!    end if
    !^CMP IF SC BEGIN
    !This coupler is performed after SC-SP coupling, so that 
    !on SP the updated coordinates are available for those
    !points which passed from SC to IH

  end subroutine couple_ih_sp
  !==================================================================
  logical function is_in_ih(Xyz_D)
    real,dimension(:),intent(in)::Xyz_D
    real:: R2
    R2 = dot_product(Xyz_D,Xyz_D)
    is_in_ih=R2>=rBoundIh**2.and.&
         all(Xyz_D < xyz_max_d(IH_)).and.all(Xyz_D >= xyz_min_d(IH_))
  end function is_in_ih
  !==================================================================!        
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
  !^CMP END IH
  !=========================================================================
  !^CMP IF SC BEGIN
  subroutine couple_sc_sp(DataInputTime)
    use CON_global_message_pass

   real,intent(in)::DataInputTime
   !-------------------------------------------------------
    if(.not.RouterScSp%IsProc)return

    tNow=DataInputTime
    ScToSp_DD=transform_matrix(tNow,&
         Grid_C(SC_)%TypeCoord, Grid_C(SP_)%TypeCoord)

    call SC_synchronize_refinement(RouterScSp%iProc0Source,&
         RouterScSp%iCommUnion)
    call set_router_from_source_2_stage(&
         GridDescriptorSource = SC_GridDescriptor, &
         GridDescriptorTarget = SP_GridDescriptor, &
         Router               = RouterScSp, &
         get_scatter_source   = SC_get_scatter_line, &
         transform            = transform_sc_to_sp, &
         interpolate_source   = interpolation_amr_gc, &
         interpolate_target   = interpolate_sp, &
         put_scatter_target   = SP_put_scatter_from_mh)
    call global_message_pass(RouterScSp, &
         nVar = nVarBuffer, &
         fill_buffer = SC_get_for_sp_and_transform, &
         apply_buffer= SP_put_from_mh)
!    if(is_proc(SP_))then
!       call SP_put_input_time(DataInputTime)  
!       call transform_to_sp_from(SC_)
!    end if
  end subroutine couple_sc_sp
  !-------------------------------------------------------------------------
  logical function is_in_sc(Xyz_D)
    real,dimension(:),intent(in)::Xyz_D
    real::R2
    real,save:: RSCMin2 = -1
    real,save:: RSCMax2 = -1
    if(RSCMin2 == -1) RSCMin2 = dot_product(exp(xyz_min_d(SC_)),(/1,0,0/))**2
    if(RSCMax2 == -1) RSCMax2 = dot_product(exp(xyz_max_d(SC_)),(/1,0,0/))**2
    R2=dot_product(Xyz_D,Xyz_D)
    if(use_comp(IH_))then            !^CMP IF IH BEGIN
       is_in_sc=R2>=rBoundSc**2.and.R2<rBoundIh**2.and.&
            R2<RSCMax2.and.R2>=RSCMin2
    else                             !^CMP END IH
       is_in_sc=R2>=rBoundSc**2.and.&
            R2<RSCMax2.and.R2>=RSCMin2
    end if                           !^CMP IF IH
  end function is_in_sc
  !--------------------------------------------------------------------------
  subroutine SC_get_for_sp_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)

    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V
    integer:: iVarBx, iVarBz
    !------------------------------------------------------------
    ! get buffer with variables
    call SC_get_for_sp(&
         nPartial,iGetStart,Get,w,State_V,nVar)
!write(*,*)State_V
    ! indices of variables 
    iVarBx = iVar_V(BxCouple_)
    iVarBz = iVar_V(BzCouple_)
    ! perform transformation before returning
    State_V(iVarBx:iVarBz)=matmul(ScToSp_DD,State_V(iVarBx:iVarBz))
  end subroutine SC_get_for_sp_and_transform
  !=========================================================================
  !^CMP END SC
end Module CON_couple_mh_sp
