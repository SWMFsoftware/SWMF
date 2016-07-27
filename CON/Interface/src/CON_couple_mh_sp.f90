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
       IH_get_scatter_line                                  !^CMP IF IH

  use SC_wrapper, ONLY: SC_synchronize_refinement, &        !^CMP IF SC
       SC_extract_line, SC_get_for_sp, SC_get_a_line_point,&!^CMP IF SC
       SC_get_scatter_line                                  !^CMP IF SC

  use CON_global_message_pass
  use CON_axes

  use SP_wrapper, ONLY: &
       SP_put_from_mh, SP_put_input_time, &
       SP_put_line, SP_request_line, SP_get_grid_descriptor_param, &
       SP_get_line_all

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
  real,   allocatable:: XyzStored_DI(:,:)
  integer,allocatable:: iAuxStored_I(:)

  real,dimension(:,:),pointer ::Xyz_DI
  logical,dimension(:),pointer :: Is_I
  integer,parameter::nPointMax=5000
  integer::nPoint=0
  integer::iPoint
  integer::iError
  real::bDxyz_I(1:6)!The interpolated values of full B and DXyz
  real::DsResolution,XyzLine_D(3)
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
    integer::nLine

    integer:: iGridMin_D(3), iGridMax_D(3), ierror
    real:: Disp_D(3)
    real, pointer:: CoordMisc_DI(:,:)

    ! available directions of interface between SP and MH 
    ! (see subroutine exchange_lines below)
    integer, parameter:: &
         iInterfaceBegin = -1, iInterfaceOrigin = 0, iInterfaceEnd = 1

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
    end if


    !\
    ! Allocate storage for field line requests
    !/
    nLine = SP_GridDescriptor%DD%Ptr%nBlockAll

    !\
    ! Extract and exchange initial data
    !/
    if(use_comp(SC_))&
         call exchange_lines(SC_)
    if(use_comp(IH_))&
         call exchange_lines(IH_)

    ! reserve memeory for SP grid
    call allocate_vector('SP_Xyz_DI', &
                  SP_GridDescriptor%DD%Ptr%nDim, &
                  product(SP_GridDescriptor%DD%Ptr%nCells_D)*nLine)
    call associate_with_global_vector(CoordMisc_DI, 'SP_Xyz_DI')
    if(is_proc(SP_))&
         call SP_get_line_all(CoordMisc_DI)
    nullify(CoordMisc_DI)
    call bcast_global_vector('SP_Xyz_DI',i_proc0(SP_),i_comm())

    if(use_comp(SC_))then
       if(RouterScSp%IsProc)then
          call allocate_mask('SP_IsInSC','SP_Xyz_DI')
          call set_mask('SP_IsInSC','SP_Xyz_DI',is_in_sc)
       end if
    end if

    if(use_comp(IH_))then
       if(RouterIhSp%IsProc)then
          call allocate_mask('SP_IsInIH','SP_Xyz_DI')
          call set_mask('SP_IsInIH','SP_Xyz_DI',is_in_ih)
       end if
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
              get_request_target   = SP_get_request_for_sc, &
              transform            = transform_sp_to_sc, &
              interpolate_source   = interpolation_amr_gc, &
              put_request_source   = SC_put_request)
         if(is_proc(SC_))&
              call SC_extract_line(&
              ubound(XyzStored_DI,2), XyzStored_DI, iInterfaceOrigin)
         call set_router_from_source_2_stage(&
              GridDescriptorSource = SC_GridDescriptor, &
              GridDescriptorTarget = SP_GridDescriptor, &
              Router               = RouterScSp, &
              get_scatter_source   = SC_get_scatter_line, &
              transform            = transform_sc_to_sp, &
              interpolate_source   = interpolation_amr_gc, &
              interpolate_target   = interpolate_sp, &
              put_scatter_target   = SP_put_scatter_from_sc)
         call global_message_pass(RouterScSp, &
              nVar = 11, &
              fill_buffer = SC_get_for_sp_and_transform, &
              apply_buffer= SP_put_from_mh)
      case(IH_)
         call set_router_from_target_2_stage(&
              GridDescriptorSource = IH_GridDescriptor, &
              GridDescriptorTarget = SP_GridDescriptor, &
              Router               = RouterIhSp, &
              get_request_target   = SP_get_request_for_ih, &
              transform            = transform_sp_to_ih, &
              interpolate_source   = interpolation_amr_gc, &
              put_request_source   = IH_put_request)
         if(is_proc(IH_))&
              call IH_extract_line(&
              ubound(XyzStored_DI,2), XyzStored_DI, iInterfaceEnd)
         call set_router_from_source_2_stage(&
              GridDescriptorSource = IH_GridDescriptor, &
              GridDescriptorTarget = SP_GridDescriptor, &
              Router               = RouterIhSp, &
              get_scatter_source   = IH_get_scatter_line, &
              transform            = transform_ih_to_sp, &
              interpolate_source   = interpolation_amr_gc, &
              interpolate_target   = interpolate_sp, &
              put_scatter_target   = SP_put_scatter_from_ih)
         call global_message_pass(RouterIhSp, &
              nVar = 11, &
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
  subroutine SP_get_request_for_sc(nData, &
       nCoord, Coord_II, iIndex_II, nAux, Aux_VI)
    integer,            intent(out):: nData
    integer,            intent(out):: nCoord
    real,   allocatable,intent(out):: Coord_II(:,:)
    integer,allocatable,intent(out):: iIndex_II(:,:)
    integer,            intent(out):: nAux
    real,   allocatable,intent(out):: Aux_VI(:,:)
    !------------------------------------------------------------
    nCoord = SP_GridDescriptor%nDim
    call SP_request_line(iInterfaceOrigin, nData, Coord_II, iIndex_II, &
         nAux, Aux_VI)
  end subroutine SP_get_request_for_sc
  !==================================================================!
  subroutine SP_get_request_for_ih(nData, &
       nCoord, Coord_II, iIndex_II, nAux, Aux_VI)
    integer,            intent(out):: nData
    integer,            intent(out):: nCoord
    real,   allocatable,intent(out):: Coord_II(:,:)
    integer,allocatable,intent(out):: iIndex_II(:,:)
    integer,            intent(out):: nAux
    real,   allocatable,intent(out):: Aux_VI(:,:)
    !------------------------------------------------------------
    nCoord = SP_GridDescriptor%nDim
    call SP_request_line(iInterfaceEnd, nData, Coord_II, iIndex_II, &
         nAux, Aux_VI)
  end subroutine SP_get_request_for_ih
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
    if(allocated(iAuxStored_I)) deallocate(iAuxStored_I)
    allocate(iAuxStored_I(nData))

    if(nData==0)&
         RETURN

    XyzStored_DI = Coord_DI
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
  subroutine SP_put_scatter_from_sc(nData, nDim, Coord_DI, nIndex, iIndex_II)
    integer, intent(in):: nDim
    integer, intent(in):: nData
    real,    intent(in):: Coord_DI(nDim, nData)
    integer, intent(in):: nIndex
    integer, intent(in):: iIndex_II(nIndex, nData)
    !--------------------------------------------------------
    call SP_put_line(nData, Coord_DI, iIndex_II, iInterfaceOrigin)
  end subroutine SP_put_scatter_from_sc
  !==================================================================!
  subroutine SP_put_scatter_from_ih(nData, nDim, Coord_DI, nIndex, iIndex_II)
    integer, intent(in):: nDim
    integer, intent(in):: nData
    real,    intent(in):: Coord_DI(nDim, nData)
    integer, intent(in):: nIndex
    integer, intent(in):: iIndex_II(nIndex, nData)
    !--------------------------------------------------------
    call SP_put_line(nData, Coord_DI, iIndex_II, iInterfaceEnd)
  end subroutine SP_put_scatter_from_ih
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
  !==================================================================!
  subroutine interpolation_amr_gc_old_interface(&
       nDim, Coord_D, GridDescriptor, &
       nIndex, iIndex_II, nImage, Weight_I)
    use CON_grid_descriptor

    integer, intent(in):: nDim
    real,    intent(inout):: Coord_D(nDim)
    type(GridDescriptorType):: GridDescriptor
    integer, intent(in) :: nIndex
    integer, intent(out):: iIndex_II(0:nIndex,2**nDim)
    integer, intent(out):: nImage
    real,    intent(out):: Weight_I(2**nDim)
    !--------------------------
    call interpolation_amr_gc(nDim, Coord_D, GridDescriptor, &
         nIndex, iIndex_II, nImage, Weight_I)
  end subroutine interpolation_amr_gc_old_interface

  !==================================================================!
  subroutine transform_from_cartesian(iComp)
    use ModCoordTransform, ONLY: xyz_to_rlonlat
    integer,intent(in)::iComp
    real,pointer,dimension(:,:)::Coord_DI
    character(len=100) :: TypeGeometry
    real:: Coord_D(3), Rho
    real, parameter:: cTol = 1E-8
    integer::nU_I(2), iParticle
    real:: SpToMh_DD(3,3)
    !------------------------------------------
    nU_I = ubound_vector('SP_Xyz_DI')
    call associate_with_global_vector(Coord_DI, 'SP_Xyz_DI')
    ! convert from cartesian coordinates if necessary                        
    TypeGeometry = Grid_C(iComp)%TypeGeometry
    SpToMh_DD=transform_matrix(tNow,&
         Grid_C(SP_)%TypeCoord,&
         Grid_C(iComp)%TypeCoord)
    if( index(TypeGeometry, 'spherical_lnr') > 0 )then
       do iParticle = 1, nU_I(2)
          ! rotate cartesian
          Coord_D(:) = matmul(SpToMh_DD, Coord_DI(:,iParticle))
          call xyz_to_rlonlat(Coord_D, Coord_DI(1:3,iParticle))
          ! convert radius to log(radius)
          if(Coord_DI(1,iParticle) > 0)&
               Coord_DI(1,iParticle) = log(Coord_DI(1,iParticle))
       end do
    elseif(index(TypeGeometry, 'cartesian') > 0 )then
       do iParticle = 1, nU_I(2)
          ! rotate cartesian
          Coord_DI(:,iParticle) = matmul(SpToMh_DD, Coord_DI(:,iParticle))
       end do
    end if
    nullify(Coord_DI)
  end subroutine transform_from_cartesian

  !==================================================================!
  subroutine transform_to_cartesian(iComp)
    use ModCoordTransform, ONLY: rlonlat_to_xyz
    integer,intent(in)::iComp
    real,pointer,dimension(:,:)::Coord_DI
    character(len=100) :: TypeGeometry
    real:: Coord_D(3), Rho
    real, parameter:: cTol = 1E-8
    integer::nU_I(2), iParticle
    real:: MhToSp_DD(3,3)
    !------------------------------------------
    nU_I = ubound_vector('SP_Xyz_DI')
    call associate_with_global_vector(Coord_DI, 'SP_Xyz_DI')
    ! convert from cartesian coordinates if necessary                        
    TypeGeometry = Grid_C(iComp)%TypeGeometry
    MhToSp_DD=transform_matrix(tNow,&
         Grid_C(iComp)%TypeCoord,&
         Grid_C(SP_)%TypeCoord)
    if( index(TypeGeometry, 'spherical_lnr') > 0 )then
       do iParticle = 1, nU_I(2)
          ! convert log(radius) to radius
          if(Coord_DI(1,iParticle) > 0)&
               Coord_DI(1,iParticle) = exp(Coord_DI(1,iParticle))
          call rlonlat_to_xyz(Coord_DI(1:3,iParticle), Coord_D)
          ! rotate cartesian
          Coord_DI(:,iParticle) = Coord_D
          Coord_DI(:,iParticle) = matmul(MhToSp_DD, Coord_DI(:,iParticle))
       end do
    elseif(index(TypeGeometry, 'cartesian') > 0 )then
       do iParticle = 1, nU_I(2)
          ! rotate cartesian
          Coord_DI(:,iParticle) = matmul(MhToSp_DD, Coord_DI(:,iParticle))
       end do
    end if
    nullify(Coord_DI)
  end subroutine transform_to_cartesian

  !==================================================================!
  subroutine transform_to_sp_from(iComp)
    integer,intent(in)::iComp
    real,pointer,dimension(:,:)::SP_LocalXyz_DI
    logical,pointer,dimension(:)::Is_I
    integer::nU_I(2),i
    real,dimension(3,3)::MhToSp_DD
    character(LEN=2)::NameComp
    real::LengthRatio
    call get_comp_info(iComp,Name=NameComp)
    MhToSp_DD=transform_matrix(tNow,&
         Grid_C(iComp)%TypeCoord,&
         Grid_C(SP_)%TypeCoord)
    if(DoTest)write(*,*)'Transform SP coordinates from '//NameComp
    call associate_with_global_mask(Is_I,'SP_IsIn'//NameComp)
    call associate_with_global_vector(SP_LocalXyz_DI,'SP_Xyz_DI')
    nU_I=ubound(SP_LocalXyz_DI)
    if(DoTest)write(*,*)nU_I
    LengthRatio=Grid_C(iComp)%UnitX/Grid_C(SP_)%UnitX
    do i=1,nU_I(2)
       if(.not.Is_I(i))CYCLE
       SP_LocalXyz_DI(:,i)=matmul(MhToSp_DD,&
            point_state_v('SP_Xyz_DI',3,i))*LengthRatio
    end do
    nullify(SP_LocalXyz_DI,Is_I)
  end subroutine transform_to_sp_from

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
    call bcast_global_vector('SP_Xyz_DI',&
         RouterIhSp%iProc0Source,&
         RouterIhSp%iCommUnion,&
         'SP_IsInIH')
    call transform_from_cartesian(IH_)
    call set_router(& 
         GridDescriptorSource=IH_GridDescriptor,&
         GridDescriptorTarget=SP_GridDescriptor,&
         Router=RouterIhSp,&
         NameMappingVector='SP_Xyz_DI',&
         NameMask='SP_IsInIH',&
         interpolate=interpolation_amr_gc_old_interface)

    call transform_to_cartesian(IH_)
    if(is_proc(SP_))then
!       call SP_put_input_time(DataInputTime)
!       call transform_to_sp_from(IH_)
    end if

    call global_message_pass(RouterIhSp,&
         nVar=11,&
         fill_buffer=IH_get_for_sp_and_transform,&
         apply_buffer=SP_put_from_mh)
    !^CMP IF SC BEGIN
    !This coupler is performed after SC-SP coupling, so that 
    !on SP the updated coordinates are available for those
    !points which passed from SC to IH

!    if(use_comp(SC_))then              
!       if(is_proc0(SP_))then
!          !Check the points which passed from SC to IH:
!          call associate_with_global_vector(Xyz_DI,'SP_Xyz_DI')
!          call associate_with_global_mask(Is_I,'SP_IsInIH')
!
!          do iPoint=1,nPoint
!             if(Is_I(iPoint))CYCLE
!             !Consider a point which before was   not in IH
!             Xyz_D=matmul(ScToIh_DD,Xyz_DI(:,iPoint))*&
!                  (Grid_C(SC_)%UnitX/Grid_C(IH_)%UnitX) 
!             if(is_in_ih(Xyz_D))& ! Now the point is in IH
!                  Xyz_DI(:,iPoint)=Xyz_D
!             !..that is why we convert it to IH coordinates
!          end do
!          nullify(Xyz_DI)
!          nullify(Is_I)
!       end if
!       call bcast_global_vector('SP_Xyz_DI',&
!            RouterIhSp%iProc0Target,&
!            RouterIhSp%iCommUnion)           !^CMP END SC
!    end if
!    call set_mask('SP_IsInIH','SP_Xyz_DI',is_in_ih)
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
    integer, parameter :: Rho_=1, Ux_=2, Uz_=4, Bx_=5, Bz_=7,&
         BuffX_    =9,BuffZ_=11
    !------------------------------------------------------------
    call IH_get_for_sp(&
         nPartial,iGetStart,Get,w,State_V,nVar)

    State_V(Ux_:Uz_)=&
         transform_velocity(tNow,&
         State_V(Ux_:Uz_),&
         State_V(BuffX_:BuffZ_),&
         Grid_C(IH_)%TypeCoord,Grid_C(SP_)%TypeCoord)

    State_V(Bx_:Bz_)=matmul(IhToSp_DD,State_V(Bx_:Bz_))

    ! transfrom coordinates
    State_V(9:11) =matmul(IhToSp_DD,State_V(9:11))

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
    call bcast_global_vector('SP_Xyz_DI',&
         RouterScSp%iProc0Source,&
         RouterScSp%iCommUnion,&
         'SP_IsInSC')
    call transform_from_cartesian(SC_)
    call set_router(& 
         GridDescriptorSource=SC_GridDescriptor,&
         GridDescriptorTarget=SP_GridDescriptor,&
         Router=RouterScSp,&
         NameMappingVector='SP_Xyz_DI',&
         NameMask='SP_IsInSC',&
         interpolate=interpolation_amr_gc_old_interface)
    
    call transform_to_cartesian(SC_)
    if(is_proc(SP_))then
!       call SP_put_input_time(DataInputTime)  
!       call transform_to_sp_from(SC_)
    end if

    call global_message_pass(RouterScSp,&
         nVar=11,&
         fill_buffer=SC_get_for_sp_and_transform,&
         apply_buffer=SP_put_from_mh)
    call set_mask('SP_IsInSC','SP_Xyz_DI',is_in_sc)
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
    integer, parameter :: Rho_=1, Ux_=2, Uz_=4, Bx_=5, Bz_=7,&
         BuffX_    =9,BuffZ_=11
    !------------------------------------------------------------
    call SC_get_for_sp(&
         nPartial,iGetStart,Get,w,State_V,nVar)

    State_V(Ux_:Uz_)=&
         transform_velocity(tNow,&
         State_V(Ux_:Uz_),&
         State_V(BuffX_:BuffZ_),&
         Grid_C(SC_)%TypeCoord,Grid_C(SP_)%TypeCoord)

    State_V(Bx_:Bz_)=matmul(ScToSp_DD,State_V(Bx_:Bz_))

    ! transform coordinates
    State_V(9:11) =matmul(ScToSp_DD,State_V(9:11))
  end subroutine SC_get_for_sp_and_transform
  !=========================================================================
  !^CMP END SC
end Module CON_couple_mh_sp
