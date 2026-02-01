!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_router

  use CON_world
  use ModMpi
  use CON_grid_descriptor, ONLY: &
       GridType, LocalGridType, coord_grid_d, nearest_grid_points, &
       bilinear_interpolation, interpolation_amr_gc, BLK_, GlobalBlock_, &
       GridPointFirst_, CellCentered_, Nodes_, set_local_gd, &
       set_standard_grid_descriptor, clean_gd, global_i_grid_point_to_icb
  use ModUtilities, ONLY: check_allocate, CON_stop

  ! This file presents the class of routers between the grids, each
  ! of them can be either the uniformly spaced or Octree or Quadric
  ! adaptive block grid                   .
  !
  ! The methods include: allocation, initialization, cleaner and
  ! two different constructors.
  !

  implicit none

  SAVE

  logical, parameter  :: UseUnionComm = .true.
  type DoAddPtrType
     logical, pointer :: DoAdd_I(:)
  end type DoAddPtrType
  ! See CON\_grid\_descriptor about iCB index. In the array iCB\_I
  ! the second index enumerates the grid points belonging to some
  ! list, while the first one numerates the position of (0) PE, at
  ! which the point is localized, (1:nDim) grid point indexes in
  ! the block and (nDim+1), if exists, stores the local block
  ! number.
  type IndexPtrType
     integer, pointer :: iCB_II(:,:)
  end type IndexPtrType

  type WeightPtrType
     real, pointer    :: Weight_I(:)
  end type WeightPtrType

  type RouterType
     ! The router can be set between grids of different dimensions
     character(LEN=3)::Name
     integer:: iCompTarget, iCompSource
     ! For the router between LOCAL grids of a component we use the
     ! communicator of the model for sending-receiving the data,
     ! otherwise the global communicator
     !
     logical :: IsLocal, IsProc
     integer :: iProc, nProc, iComm
     ! If the union group is constructed, then for use with broadcast
     ! we need the union communicator and the root PE ranks in this
     ! communicator
     !
     integer :: iCommUnion, iProc0Source, iProc0Target
     !
     ! As the default we use iCB indexes to construct the router,
     ! hence the grid point is characterized by the
     ! Grid%nDim grid point indexes plus one more index for
     ! the block number. Also we allow to use exactly
     ! Grid%nDim indexes, without the block number which
     ! only seems to be of sence for the component which is localized
     ! at one PE only, or which has exactly one block per PE
     !
     integer :: nIndexSource, nIndexTarget
     !
     ! The total amounts of the buffer segments to be sent-received
     ! to/from the PE. The total amounts of the grid points from which
     ! the data should be got or to which the data should be put,some
     ! data points may be counted more than one time
     !
     integer, pointer :: nGet_P(:), nPut_P(:), nRecv_P(:), nSend_P(:)
     ! iCB indexes and the weight coefficients for the points of the  !
     ! target and source grids, which are connected through the router!
     type(IndexPtrType),  pointer :: iGet_P(:)
     type(IndexPtrType),  pointer :: iPut_P(:)
     type(DoAddPtrType),  pointer :: DoAdd_P(:)
     type(WeightPtrType), pointer :: Get_P(:)
     type(WeightPtrType), pointer :: Put_P(:)
     ! Mapped gird points.
     real, pointer :: BufferSource_II(:,:)
     real, pointer :: BufferTarget_II(:,:)
     ! As long as we do not keep the entire mapping vector, we may want to know
     ! a list of global point numbers for the mapped points, i.e.
     ! nMappedPointIndex=1 integer per each mapped point. It may be convenient
     ! to keep more than one index (a global tree node number, global block
     ! number etc, with the only restriction that these indexes are sufficient
     ! to recover the global point number for the mapped point.)

     integer :: nMappedPointIndex

     integer :: nBufferSource, nBufferTarget
     integer :: nVar, nDim, iAuxStart, iAuxEnd
  end type RouterType
  ! The best effeiciency is usually achieved with the MPI pair
  ! RSend+iRecv. Sometimes the fully non-blocking pair iSend+iRecv
  ! is more stable and rarely it is more efficient either. In
  ! the first case use DoRSend=.true. (default), otherwise use
  ! DoRSend=.false.
  logical, parameter :: DoRSend=.false.

  ! aux integer arrays to put data in BufferS_II in the correct order
  integer, allocatable :: iProc_I(:), iOrder_I(:), iTranslated_P(:), iAux_P(:)
  private :: allocate_get_arrays
  private :: allocate_put_arrays
  private :: allocate_buffer_target
  private :: allocate_buffer_source
  private :: check_router_allocation
  private :: check_size
contains
  !============================================================================
  subroutine init_router(GridSource, GridTarget, Router, &
       nIndexSource, nIndexTarget, nMappedPointIndex)
    type(GridType),    intent(in)  :: GridSource, GridTarget
    type (RouterType), intent(out) :: Router
    integer, optional, intent(in)  :: nIndexSource
    integer, optional, intent(in)  :: nIndexTarget
    integer, optional, intent(in)  :: nMappedPointIndex
    integer :: iPE, iError, nProc
    integer :: iProc0Source,iProc0Target, iProcUnion
    integer :: iGroupUnion, iGroupSource, iGroupTarget, iGroup
    !--------------------------------------------------------------------------
    if(i_proc() < 0)then
       Router%IsProc = .false.
       RETURN
    end if
    Router%iCompSource = GridSource%Domain%Ptr%CompID_
    Router%iCompTarget = GridTarget%Domain%Ptr%CompID_

    ! Check if the grids are both local or both global
    !
    if(GridSource%Domain%Ptr%IsLocal.and.GridTarget%Domain%Ptr%IsLocal)then
       Router%IsLocal=.true.
       if( Router%iCompSource /= Router%iCompTarget)&
            call CON_stop('Do not couple Local grids of different components!')

       Router%iProc        = i_proc(Router%iCompTarget)
       Router%nProc        = n_proc(Router%iCompTarget)
       Router%iComm        = i_comm(Router%iCompTarget)
       Router%iCommUnion   = Router%iComm
       Router%iProc0Source = 0
       Router%iProc0Target = 0
       Router%IsProc       = is_proc(Router%iCompTarget)

    elseif((.not.GridSource%Domain%Ptr%IsLocal)&
         .and.(.not.GridTarget%Domain%Ptr%IsLocal))then
       Router%IsLocal      = .false.
       Router%iProc        = i_proc()
       nProc               = n_proc()
       Router%nProc        = nProc
       Router%iComm        = i_comm()
       Router%IsProc       = is_proc()
       iProc0Source        = i_proc0(Router%iCompSource)
       iProc0Target        = i_proc0(Router%iCompTarget)

       if(UseUnionComm)then
          if(.not.allocated(iAux_P))then
             allocate(iAux_P(0:nProc - 1), stat = iError)
             call check_allocate(iError,'iAux_P')
             allocate(iTranslated_P(0:nProc-1), stat = iError)
             call check_allocate(iError,'iTranslated_P')
             do iPE = 0, nProc - 1
                iAux_P(iPE) = iPE
             end do
          end if

          iGroupSource = i_group(Router%iCompSource)
          iGroupTarget = i_group(Router%iCompTarget)
          if(iProc0Target>iProc0Source)then
             call MPI_GROUP_UNION(iGroupSource, iGroupTarget, &
                  iGroupUnion, iError)
             Router%iProc0Source = 0
             call MPI_GROUP_TRANSLATE_RANKS(i_group(), n_proc(), iAux_P(0),&
                  iGroupUnion, iTranslated_P, iError)
             Router%iProc0Target = iTranslated_P(iProc0Target)
          else
             call MPI_GROUP_UNION(iGroupTarget, iGroupSource, &
                  iGroupUnion, iError)
             Router%iProc0Target=0
             call MPI_GROUP_TRANSLATE_RANKS(i_group(), n_proc(), iAux_P(0),&
                  iGroupUnion, iTranslated_P, iError)
             Router%iProc0Source = iTranslated_P(iProc0Source)
          end if
          call MPI_COMM_CREATE(i_comm(), iGroupUnion, Router%iCommUnion, iError)
          call MPI_group_rank(iGroupUnion, iProcUnion, iError)
          Router%IsProc = iProcUnion /= MPI_UNDEFINED
          if(iProcUnion /= Router%iProc0Target.and.i_proc()==iProc0Target)&
               call CON_stop('Wrongly defined Router%iProc0Target')
          if(iProcUnion /= Router%iProc0Source.and.i_proc()==iProc0Source)&
               call CON_stop(&
               'Wrongly defined Router%iProc0Source')
          call MPI_GROUP_FREE(iGroupUnion, iError)
       else  ! Do not use comm union
          Router%iCommUnion   = Router%iComm
          Router%iProc0Source = iProc0Source
          Router%iProc0Target = iProc0Target
       end if
    else
       call CON_stop('Do not couple a Local grid with a global one')
    end if

    Router%nIndexSource = GridSource%nDim + 1
    if(present(nIndexSource))then
       if(nIndexSource >= Router%nIndexSource - 1.or. nIndexSource==1)then
          Router%nIndexSource=nIndexSource
       else
          write(*,*)'IndexMin=', Router%nIndexSource - 1
          call CON_stop('nIndexSource should be at least IndexMin')
       end if
    end if

    Router%nIndexTarget = GridTarget%nDim + 1

    if(present(nIndexTarget))then
       if(nIndexTarget>=Router%nIndexTarget - 1.or.nIndexTarget==1)then
          Router%nIndexTarget=nIndexTarget
       else
          write(*,*)'IndexMin=', Router%nIndexTarget - 1
          call CON_stop('nIndexTarget should be at least IndexMin')
       end if
    end if

    Router%nMappedPointIndex   = 0
    if(present(nMappedPointIndex))then
       Router%nMappedPointIndex = nMappedPointIndex
       select case(nMappedPointIndex)
       case(0,2)
          ! Allowed so far. Needs more work
       case default
          call CON_stop('Illegal number of mapped point indexes')
       end select
    end if
    ! Set dimensionality of the coordinate vectors, to be sent via
    ! the router while constructing it.

    ! This below choice is only valid for the routers constructed
    ! from target (using set_router). The points in the target are
    ! mapped to the source, hence, the sent point should be at the
    ! same dimensionality as the source. Should be reset for
    ! the router constructed from source by calling
    ! set_router_dim(Router, GridTarget%nDim)
    call set_router_dim(Router, GridSource%nDim)

    nProc=Router%nProc
    nullify(Router%BufferSource_II)
    call allocate_buffer_source(Router, nProc)
    Router%nBufferSource = nProc
    nullify(Router%BufferTarget_II)
    call allocate_buffer_target(Router, nProc)
    Router%nBufferTarget = nProc

    ! Allocation:
    allocate(Router%nGet_P(0:nProc-1),  stat=iError)
    call check_allocate(iError,'nGet_P')
    allocate(Router%nPut_P(0:nProc-1),  stat=iError)
    call check_allocate(iError,'nPut_P')
    allocate(Router%nSend_P(0:nProc-1), stat=iError)
    call check_allocate(iError,'nSend_P')
    allocate(Router%nRecv_P(0:nProc-1), stat=iError)
    call check_allocate(iError,'nRecv_P')
    allocate(Router%iGet_P(0:nProc-1),  stat=iError)
    call check_allocate(iError,'iGet_P')
    allocate(Router%Get_P(0:nProc-1),   stat=iError)
    call check_allocate(iError,'Get_P')

    do iPE = 0, nProc - 1
       nullify(Router%iGet_P(iPE)%iCB_II)
       nullify(Router%Get_P(iPE)%Weight_I)
       call allocate_get_arrays(Router, iPE, 1)
    end do
    allocate(Router%iPut_P(0:nProc-1), stat=iError)
    call check_allocate(iError,'iPut_P')
    allocate(Router%Put_P(0:nProc-1),  stat=iError)
    call check_allocate(iError,'Put_P')
    allocate(Router%DoAdd_P(0:nProc-1),stat=iError)
    call check_allocate(iError,'DoAdd_P')

    do iPE = 0, nProc - 1
       nullify(Router%iPut_P(iPE)%iCB_II)
       nullify(Router%Put_P(iPE)%Weight_I)
       nullify(Router%DoAdd_P(iPE)%DoAdd_I)
       call allocate_put_arrays(Router, iPE, 1)
    end do

    Router%nGet_P  = 0
    Router%nPut_P  = 0
    Router%nSend_P = 0
    Router%nRecv_P = 0
    call check_router_allocation(Router)
  end subroutine init_router
  !============================================================================
  subroutine set_router_dim(Router, nDim)
    ! Set dimensionality of the coordinate vectors, to be sent via
    ! the router while constructing it.

    ! For the routers constructed from target (using set_router) nDim should
    ! equal to GridSource%nDim. The points in the target are mapped to the
    ! source, hence, the sent point should be at the same dimensionality
    ! as the source. Similarly, for the router constructed from source
    ! nDim should be equal to GridTarget%nDim
    type(RouterType),intent(inout) :: Router
    integer,         intent(in)    :: nDim
    !--------------------------------------------------------------------------
    Router%nDim      = nDim
    Router%iAuxStart = Router%nDim + 1
    Router%iAuxEnd   = Router%nDim + Router%nMappedPointIndex
    Router%nVar      = Router%iAuxEnd
  end subroutine set_router_dim
  !============================================================================
  subroutine  allocate_get_arrays(Router, iPE, nLength)
    type(RouterType),intent(inout) :: Router
    integer,         intent(in)    :: iPE, nLength
    integer :: iError
    !--------------------------------------------------------------------------
    if(associated(Router%iGet_P(iPE)%iCB_II))&
         deallocate(Router%iGet_P(iPE)%iCB_II)
    if(associated(Router%Get_P(iPE)%Weight_I))&
         deallocate(Router%Get_P(iPE)%Weight_I)
    allocate(Router%iGet_P(iPE)%iCB_II(&
         0:Router%nIndexSource, max(1, nLength)), stat=iError)
    call check_allocate(iError, 'iGet_P%iCB_II')
    Router%iGet_P(iPE)%iCB_II = 0
    allocate(Router%Get_P(iPE)%Weight_I(max(1, nLength)), stat=iError)
    call check_allocate(iError, 'Get_P%Weight_I')
    Router%Get_P(iPE)%Weight_I = 0.0
  end subroutine allocate_get_arrays
  !============================================================================
  subroutine allocate_buffer_source(Router, nLength)
    type(RouterType),intent(inout) :: Router
    integer,         intent(in)    :: nLength
    integer::iError
    !--------------------------------------------------------------------------
    if(associated(Router%BufferSource_II))&
         deallocate(Router%BufferSource_II)
    allocate(Router%BufferSource_II(Router%nVar, nLength), stat=iError)
    call check_allocate(iError, 'BufferSource_II')
    Router%BufferSource_II = 0.0
    Router%nBufferSource = nLength
  end subroutine allocate_buffer_source
  !============================================================================
  subroutine allocate_put_arrays(Router,iPE,nLength)
    type(RouterType),intent(inout)::Router
    integer,intent(in)::iPE,nLength
    integer::iError
    !--------------------------------------------------------------------------
    if(associated(Router%iPut_P(iPE)%iCB_II))&
         deallocate(Router%iPut_P(iPE)%iCB_II)
    if(associated(Router%DoAdd_P(iPE)%DoAdd_I))&
         deallocate(Router%DoAdd_P(iPE)%DoAdd_I)
    allocate(Router%iPut_P(iPE)%iCB_II(&
         0:Router%nIndexTarget, max(1, nLength)), stat=iError)
    call check_allocate(iError,'iPut_P%iCB_II')
    Router%iPut_P(iPE)%iCB_II = 0
    allocate(Router%Put_P(iPE)%Weight_I(max(1, nLength)), stat=iError)
    call check_allocate(iError,'Put_P%Weight_I')
    Router%Put_P(iPE)%Weight_I = 0.0
    allocate(Router%DoAdd_P(iPE)%DoAdd_I(max(1, nLength)), stat=iError)
    call check_allocate(iError,'DoAdd_P%DoAdd_I')
    Router%DoAdd_P(iPE)%DoAdd_I = .false.
  end subroutine allocate_put_arrays
  !============================================================================
  subroutine allocate_buffer_target(Router, nLength)
    type(RouterType),intent(inout) :: Router
    integer,         intent(in)    :: nLength
    integer::iError
    !--------------------------------------------------------------------------
    if(associated(Router%BufferTarget_II))&
         deallocate(Router%BufferTarget_II)
    allocate(Router%BufferTarget_II(Router%nVar, nLength), stat=iError)
    call check_allocate(iError, 'BufferTarget_II')
    Router%BufferTarget_II = 0.0
    Router%nBufferTarget   =  nLength
  end subroutine allocate_buffer_target
  !============================================================================
  subroutine check_router_allocation(Router)
    type(RouterType),intent(inout)::Router
    integer :: iPE, nTotalPut, nTotalGet, UBound_I(2)
    !--------------------------------------------------------------------------
    do iPE = 0, Router%nProc-1
       if(ubound(Router%iPut_P(iPE)%iCB_II,2) < Router%nPut_P(iPE))then
          call allocate_put_arrays(Router, iPE, Router%nPut_P(iPE))
       end if
       if(ubound(Router%iGet_P(iPE)%iCB_II,2) < Router%nGet_P(iPE))then
          call allocate_get_arrays(Router, iPE, Router%nGet_P(iPE))
       end if
    end do
    if(Router%nVar > 0)then
       UBound_I = ubound(Router%BufferTarget_II)
       if(  Ubound_I(1)/=Router%nVar.or.&
            UBound_I(2) < nlength_buffer_target(Router))&
            call allocate_buffer_target(Router, &
            max(Router%nProc, nlength_buffer_target(Router)))
       UBound_I = ubound(Router%BufferSource_II)
       if(  Ubound_I(1)/=Router%nVar.or.&
            UBound_I(2) < nlength_buffer_source(Router))&
            call allocate_buffer_source(Router, &
            max(Router%nProc,nlength_buffer_source(Router)))
    end if
  end subroutine check_router_allocation
  !============================================================================
  integer function nlength_buffer_source(Router)
    type(RouterType),intent(inout)::Router
    !--------------------------------------------------------------------------
    nlength_buffer_source = sum(Router%nSend_P(:))
  end function nlength_buffer_source
  !============================================================================
  integer function nlength_buffer_target(Router)
    type(RouterType),intent(inout)::Router
    !--------------------------------------------------------------------------
    nlength_buffer_target = sum(Router%nRecv_P(:))
  end function nlength_buffer_target
  !============================================================================
  subroutine set_router(GridSource, GridTarget, Router, &
       is_interface_block, n_interface_point_in_block,  &
       interface_point_coords, mapping, interpolate, extra_data)

    ! GridSource is the Grid for the Source component
    ! GridTarget   the Grid for the Target component
    ! Router is  the router to be set
    ! is_interface_block ogical function which allows to skip the block
    ! if there is no interface points in it.
    ! Optional, if not present then all the
    ! blocks are checked for the presence of the interface points

    ! n_interface_point_in_block:
    ! The function may be also used for sparse 1D grids
    ! in which the real information is available only in the first
    ! n grid points, for the given plocks, the other points being
    ! unused. For this purpose, the function should return this
    ! number of used points. It is not applied, if the rurned value
    ! is larger than the total number of grid points, in the given
    ! plock.

    ! interface_point_coords:
    ! The subroutine which defines if the grid point is inside the
    ! interface layer. Optional, if not present, then all the grid
    ! points (at the target grid) are considered as the interface
    ! layer points
    ! mapping: transformation of the location coordinates between components
    ! interpolate: interpolation subroutine for the Source's grid
    ! extra_data: Procedure determining which data are sent from
    ! target to source together with the semi-router
    type(GridType),      intent(in)   :: GridSource
    type(LocalGridType), intent(in)   :: GridTarget
    type(RouterType),    intent(inout):: Router
    interface
       logical function is_interface_block(iBlockLocal)
         implicit none
         integer, intent(in)::iBlockLocal
       end function is_interface_block
       integer function n_interface_point_in_block(iBlockLocal)
         implicit none
         integer,intent(in) :: iBlockLocal
       end function n_interface_point_in_block
       subroutine interface_point_coords(nDim, Xyz_D, nIndex, iIndex_I,&
            IsInterfacePoint)
         implicit none
         integer,          intent(in) :: nDim, nIndex
         logical,         intent(out) :: IsInterfacePoint
         real,          intent(inout) :: Xyz_D(nDim)
         integer,       intent(inout) :: iIndex_I(nIndex)
       end subroutine interface_point_coords
       !-----------------------------------------------------------
       subroutine mapping(nDimIn, CoordIn_D, nDimOut, CoordOut_D, &
            IsInterfacePoint)
         ! this subroutine maps coordinates between components
         integer, intent(in)  :: nDimIn
         real,    intent(in)  :: CoordIn_D(nDimIn)
         integer, intent(in)  :: nDimOut
         real,    intent(out) :: CoordOut_D(nDimOut)
         logical, intent(out) :: IsInterfacePoint
       end subroutine mapping
       !----------------------------------------------------------------------
       subroutine interpolate(nDim, Coord_D, Grid, nIndex, iIndex_II, nImage,&
            Weight_I)
         ! interpolation on Source's grid;
         ! provides PE and indices to access images of
         ! data location and interpolation weights
         use CON_grid_descriptor
         implicit none
         ! number of indices per data entry
         integer, intent(in)      :: nDim
         ! data location on Source
         real,    intent(inout)   :: Coord_D(nDim)
         ! grid descriptor
         type(GridType) :: Grid
         ! indices of images, their number and interpolation weights
         integer, intent(in)      :: nIndex
         integer, intent(out)     :: iIndex_II(0:nIndex,2**nDim)
         integer, intent(out)     :: nImage
         real,    intent(out)     :: Weight_I(2**nDim)
       end subroutine interpolate
       subroutine extra_data(BlockIndex_I,&
            iGridPointInBlock, nAux, Data_V)
         implicit none

         ! Parameters are also available in CON_grid_descriptor
         integer, parameter :: GlobalTreeNode_=1, BLK_ =2
         integer, parameter :: GridPointFirst_=3, GlobalBlock_= 4

         ! Mapped point in the buffer is characterized
         ! by any of the block indexes: global tree node number,
         ! global block number, local block number or the number of
         ! the first point in the block...
         integer, intent(in) :: BlockIndex_I(1:4)
         ! and the point number in the block
         integer, intent(in) :: iGridPointInBlock
         integer, intent(in) :: nAux ! How many reals you may send
         real,    intent(out):: Data_V(1:nAux)
       end subroutine extra_data
    end interface
    optional:: is_interface_block
    optional:: n_interface_point_in_block
    optional:: interface_point_coords
    optional:: mapping
    optional:: interpolate
    optional:: extra_data
    character(len=*), parameter:: NameSub = 'set_router'
    !--------------------------------------------------------------------------
    ! Return if the processor does not belong to the communicator
    if(.not.Router%IsProc)RETURN
    if(is_proc(Router%iCompTarget))&
         call set_semi_router_from_target(&
         ! the Grid for the Source component
         GridSource,      &
         ! the Grid for the Target component
         GridTarget,      &
         ! the router to be set
         Router,                    &
         ! Logical function which allows to skip the block if there is no
         ! interface points in it. Optional, if not present then all the
         ! blocks are checked for the presence of the interface points
         is_interface_block,        &
         ! The function may be also used for sparse 1D grids
         ! in which the real information is available only in the first
         ! n grid points, for the given plocks, the other points being
         ! unused. For this purpose, the function should return this
         ! number of used points. It is not applied, if the rurned value
         ! is larger than the total number of grid points, in the given
         ! plock.
         n_interface_point_in_block,&
         ! The subroutine which defines if the grid point is inside the
         ! interface layer. Optional, if not present, then all the grid
         ! points (at the target grid) are considered as the interface
         ! layer points
         interface_point_coords,    &
         ! transformation of the location coordinates between components
         mapping,                   &
         ! interpolation subroutine for the Source's grid
         interpolate,               &
         ! Procedure determining which data are sent from target to source
         ! together with the semi-router
         extra_data)
    call synchronize_router_from_target(Router)
    if(is_proc(Router%iCompSource))call update_semi_router_at_source(&
         Router, GridSource, interpolate)
  end subroutine set_router
  !============================================================================
  subroutine set_semi_router_from_target(&
       GridSource, GridTarget, Router, is_interface_block, &
       n_interface_point_in_block, interface_point_coords, &
       mapping, interpolate, extra_data)

    ! GridSource: the Grid for the Source component
    ! GridTarget: the Grid for the Target component
    ! Router:  the router to be set
    ! is_interface_block:
    ! Logical function which allows to skip the block if there is no
    ! interface points in it. Optional, if not present then all the
    ! blocks are checked for the presence of the interface points

    ! n_interface_point_in_block:
    ! The function may be also used for sparse 1D grids
    ! in which the real information is available only in the first
    ! n grid points, for the given plocks, the other points being
    ! unused. For this purpose, the function should return this
    ! number of used points. It is not applied, if the rurned value
    ! is larger than the total number of grid points, in the given
    ! plock.

    ! interface_point_coords:
    ! The subroutine which defines if the grid point is inside the
    ! interface layer. Optional, if not present, then all the grid
    ! points (at the target grid) are considered as the interface
    ! layer points

    ! mapping: transformation of the location coordinates between components
    ! interpolate: interpolation subroutine for the Source's grid
    ! extra_data:
    ! Procedure determining which data are sent from target to source
    ! together with the semi-router

    type(GridType),      intent(in)   :: GridSource
    type(LocalGridType), intent(in)   :: GridTarget
    type(RouterType),    intent(inout):: Router
    interface
       logical function is_interface_block(iBlockLocal)
         implicit none
         integer, intent(in)::iBlockLocal
       end function is_interface_block
       integer function n_interface_point_in_block(iBlockLocal)
         implicit none
         integer,intent(in) :: iBlockLocal
       end function n_interface_point_in_block
       subroutine interface_point_coords( &
            nDim, Xyz_D, nIndex, iIndex_I,&
            IsInterfacePoint)
         implicit none
         integer,          intent(in) :: nDim, nIndex
         logical,         intent(out) :: IsInterfacePoint
         real,          intent(inout) :: Xyz_D(nDim)
         integer,       intent(inout) :: iIndex_I(nIndex)
       end subroutine interface_point_coords
       !-----------------------------------------------------------
       subroutine mapping(nDimIn, CoordIn_D, nDimOut, CoordOut_D, &
            IsInterfacePoint)
         ! this subroutine mapss coordinates between components
         integer, intent(in)  :: nDimIn
         real,    intent(in)  :: CoordIn_D(nDimIn)
         integer, intent(in)  :: nDimOut
         real,    intent(out) :: CoordOut_D(nDimOut)
         logical, intent(out) :: IsInterfacePoint
       end subroutine mapping
       !----------------------------------------------------------------------
       subroutine interpolate(&
            nDim, Coord_D, Grid, &
            nIndex, iIndex_II, nImage, Weight_I)
         ! interpolation on Source's grid;
         ! provides PE and indices to access images of
         ! data location and interpolation weights
         use CON_grid_descriptor
         implicit none
         ! number of indices per data entry
         integer, intent(in)      :: nDim
         ! data location on Source
         real,    intent(inout)   :: Coord_D(nDim)
         ! grid descriptor
         type(GridType) :: Grid
         ! indices of images, their number and interpolation weights
         integer, intent(in)      :: nIndex
         integer, intent(out)     :: iIndex_II(0:nIndex,2**nDim)
         integer, intent(out)     :: nImage
         real,    intent(out)     :: Weight_I(2**nDim)
       end subroutine interpolate

       subroutine extra_data(BlockIndex_I, iGridPointInBlock, nAux, Data_V)

         implicit none

         ! Parameters are also available in CON_grid_descriptor
         integer, parameter:: GlobalTreeNode_=1, BLK_ =2
         integer, parameter:: GridPointFirst_=3, GlobalBlock_= 4

         ! Mapped point in the buffer is characterized
         ! by any of the block indexes: global tree node number,
         ! global block number, local block number or the number of
         ! the first point in the block...
         integer, intent(in) :: BlockIndex_I(1:4)
         ! and the point number in the block
         integer, intent(in) :: iGridPointInBlock
         integer, intent(in) :: nAux ! How many reals you may send
         real,    intent(out):: Data_V(1:nAux)
       end subroutine extra_data
    end interface

    optional:: is_interface_block
    optional:: n_interface_point_in_block
    optional:: interface_point_coords
    optional:: mapping
    optional:: interpolate
    optional:: extra_data

    ! The presence of optional parameters
    logical :: UseMappingFunction, DoCheckBlock, DoCheckBlockSize, &
         DoCheckPoint, DoInterpolate
    ! MPI-related variables
    integer :: iProc
    logical :: DoCountOnly
    integer :: iBlockTo, iProcFrom, iLocalGridPoint
    integer :: iLocalPointLast, iPointInBlock, nInterfacePoint
    logical :: IsInterfacePoint
    integer :: nPutUbound_P(0:Router%nProc-1)
    integer :: nRecvUbound

    real      :: CoordTarget_D(GridTarget%nDim)
    real      :: CoordSource_D(GridSource%nDim)
    integer   :: iCell_D(GridTarget%nDim)
    integer   :: iIndexRecv_I(Router%nIndexTarget)
    ! components ids
    integer   :: iCompSource, iProc0Source, iProcLastSource, iProcStrideSource
    ! dimensionality of components
    integer   :: nDimSource, nDimTarget
    ! number of indices (e.g. cell indices) on components
    integer   :: nIndexSource, nIndexTarget
    ! loop variables
    integer   :: nBuffer, iBlockUsed
    ! offset in buffer
    integer   :: nRecvCumSum

    ! Misc
    integer   :: iBlockMisc,  nAux
    ! Return if the processor does not belong to the communicator
    character(len=*), parameter:: NameSub = 'set_semi_router_from_target'
    !--------------------------------------------------------------------------
    if(.not.Router%IsProc) RETURN
    ! Return if the processor does not belong to the target
    if(.not.is_proc(Router%iCompTarget))RETURN
    ! Stage 1:
    ! on Target PE determine the sources of data and
    ! how much data is to be requested from them

    ! For given PE the index in the communicator is:
    iProc       = Router % iProc
    iCompSource = Router%iCompSource
    if(Router%IsLocal)then
       iProc0Source = 0;  iProcStrideSource = 1
       iProcLastSource = Router%nProc - 1
    else
       iProc0Source      = i_proc0(      iCompSource)
       iProcLastSource   = i_proc_last(  iCompSource)
       iProcStrideSource = i_proc_stride(iCompSource)
    end if
    ! determine which optional actions should be taken
    UseMappingFunction = present(mapping)
    DoInterpolate      = present(interpolate)
    DoCheckBlock       = present(is_interface_block)
    DoCheckBlockSize   = present(n_interface_point_in_block)
    DoCheckPoint       = present(interface_point_coords)

    ! introduced for a better readability
    nDimTarget   = GridTarget%nDim
    nDimSource   = GridSource%nDim
    nIndexTarget = Router%nIndexTarget
    nIndexSource = Router%nIndexSource

    if(.not.UseMappingFunction .and. nDimTarget /= nDimSource )call CON_stop(&
         NameSub//': Mapping is needed for Target%nDim/=Source%nDim')

    ! some data will be sent to Source, determine amount:
    ! cell and block indexes are sent
    nAux = Router%nMappedPointIndex

    ! Check dimensions
    DoCountOnly=.true. ! To enter the loop
    do while(DoCountOnly)
       call check_router_allocation(Router)

       ! Store Upper bounds to control if the allocated
       ! index arrays have sufficient size
       do iProcFrom = iProc0Source, iProcLastSource, iProcStrideSource
          nPutUbound_P(iProcFrom) = ubound(Router%iPut_P(iProcFrom)%iCB_II,2)
       end do
       ! and also the buffer size
       nRecvUbound = ubound(Router%BufferTarget_II, 2)

       ! which processor holds a current image
       call check_size(1, [Router%nBufferTarget], iBuffer_I = iProc_I)
       ! correct order of images in the send buffer
       call check_size(1, [Router%nBufferTarget], iBuffer_I = iOrder_I)
       ! reset basic coupling information
       Router%nPut_P  = 0
       Router%nRecv_P = 0
       ! reset index of a current data entry in the buffer
       nBuffer        = 0
       DoCountOnly    = .false.

       ! If the check shows that the allocated array is not
       ! sufficient, then DoCountOnly will be set to true. The loop
       ! then will be repeated for the second time

       ! Loop over global block number (to be revised)
       BLOCKS: do iBlockUsed = 1, GridTarget%nBlock

          ! Convert to the local block number
          iBlockTo = GridTarget%iIndex_IB(BLK_,iBlockUsed)

          ! Skip the block if there is known to be no interface point in it
          if(DoCheckBlock)then
             if(.not.is_interface_block(iBlockTo))CYCLE BLOCKS
          end if

          ! Prepare a GlobalCellNumber Loop, for a given (octree) block

          ! 1. Global number of the last grid point in the previous
          ! global block, with the global block number = nBlockAll-1
          iLocalPointLast = GridTarget%iIndex_IB(GridPointFirst_,iBlockUsed) - 1
          nInterfacePoint = GridTarget%nPointPerBlock
          if( DoCheckBlockSize)then
             nInterfacePoint = min(nInterfacePoint,&
                  n_interface_point_in_block(iBlockTo))
          end if

          ! Now the subloop runs over the iPointBlock within
          ! the current local block
          POINTS:do iPointInBlock  = 1, nInterfacePoint

             ! Recover the global grid point number
             iLocalGridPoint =   iLocalPointLast + iPointInBlock

             ! Get cell index for this grid point in target.
             call global_i_grid_point_to_icb(&
                  GridTarget,&
                  iLocalGridPoint,&
                  iBlockMisc,     & !<=The output is iBlockUsed, which is
                  iCell_D)          ! the loop variable, cannot be affected

             ! Shape an index array for this grid point in target
             iIndexRecv_I(nIndexTarget) = iBlockTo
             iIndexRecv_I(1:nDimTarget) = iCell_D

             ! Generalized coordinates of the target grid point
             CoordTarget_D = coord_grid_d(GridTarget, iBlockUsed, iCell_D)

             ! Using the grid descriptor for the target we may want
             ! 1. To select or deselect this point from the interface
             ! 2. To convert (with no change in the dimensions) the
             ! grid point coordinates (say, from generalized to
             ! Cartesian), that is why CoordTarget_D has intent inout
             if( DoCheckPoint)then
                call interface_point_coords(&
                     nDimTarget,            &
                     CoordTarget_D,         &
                     nIndexTarget,          &
                     iIndexRecv_I,          &
                     IsInterfacePoint)
                if(.not.IsInterfacePoint)CYCLE POINTS

                ! Otherwise, no point is deselected,
                ! the generalized coordinates keep unchanged
             end if
             if(UseMappingFunction)then

                ! This is a mapping of the
                ! target point coordinates to the generalized
                ! coordinates of source
                call mapping(&
                     nDimTarget,&
                     CoordTarget_D,&
                     nDimSource,&
                     CoordSource_D, &
                     IsInterfacePoint)
                if(.not.IsInterfacePoint)CYCLE POINTS
             else
                CoordSource_D=CoordTarget_D
             end if

             ! Now, we put the point indexes to the router
             ! and its coordinates (and, probably,
             ! its global point number) to the buffer
             call put_point_to_semirouter
          end do POINTS   ! iGlobalPoints
       end do BLOCKS      ! iGlobalBlock
    end do
    ! fix the order of Buffer_I so contiguous chunks of data can be sent
    ! to the appropriate processors of Source,
    ! currently iOrder_I contains indices WITHIN these chunks
    nRecvCumSum = Router%nRecv_P(iProc)
    do iProcFrom = iProc0Source, iProcLastSource, iProcStrideSource
       if(iProcFrom == iProc)CYCLE
       where(iProc_I( 1:nBuffer) == iProcFrom)&
            iOrder_I( 1:nBuffer) = iOrder_I( 1:nBuffer) + nRecvCumSum
       nRecvCumSum = nRecvCumSum + Router%nRecv_P(iProcFrom)
    end do

    ! the correct order is found, apply it
    Router%BufferTarget_II(:,iOrder_I(1:nBuffer)) = &
         Router%BufferTarget_II(:,1:nBuffer)
  contains
    !==========================================================================
    subroutine put_point_to_semirouter
      integer :: iIndexGet_II(0:Router%nIndexSource,&
           2**GridSource%nDim)
      real    :: Weight_I(2**GridSource%nDim)
      real    :: CoordPass_D(GridSource%nDim)
      integer :: iImage, nImage, iProcDoNotAdd
      integer :: iProcLookUp_I(2**GridSource%nDim)
      integer :: nProcToGet, iProcToGet
      integer, parameter :: iPointInBlock_ = 2, iBlockAll_ = 1
      integer :: iAux_I(1:2), iBlockAll
      ! error message containers
      character(len=200):: StringErrorFormat, StringErrorMessage
      !------------------------------------------------------------------------
      CoordPass_D = CoordSource_D
      if( DoInterpolate)then
         call interpolate(  &
              nDimSource,   &
              CoordPass_D,    &
              GridSource,     &
              nIndexSource, &
              iIndexGet_II, &
              nImage,       &
              Weight_I)
      else
         call nearest_grid_points(&
              nDim           = nDimSource, &
              Coord_D        = CoordPass_D, &
              Grid = GridSource, &
              nIndex         = nIndexSource, &
              iIndex_II      = iIndexGet_II, &
              nImage         = nImage, &
              Weight_I       = Weight_I)
      end if
      if(nImage < 1)then
         write(StringErrorFormat,'(a,i3,a)') '(a,',nDimSource,'es15.7)'
         write(StringErrorMessage,StringErrorFormat)&
              NameSub//': Interpolation failed at location ', &
              CoordSource_D
         call CON_stop(StringErrorMessage)
      end if

      ! go over the list of images and process result of interpolation
      nProcToGet = 0
      IMAGES1:do iImage = 1, nImage
         iProcFrom = iIndexGet_II(0, iImage)
         if(nProcToGet > 0)then
            if( any(iProcLookUp_I(1:nProcToGet) == iProcFrom)  &
                 )CYCLE IMAGES1
         end if
         nProcToGet = nProcToGet + 1
         iProcLookUp_I(nProcToGet) = iProcFrom
         ! index of current data entry in the buffer
         Router%nPut_P(iProcFrom)  = Router%nPut_P(iProcFrom) + 1
         Router%nRecv_P(iProcFrom) = Router%nRecv_P(iProcFrom) + 1

         DoCountOnly = DoCountOnly.or.&
              Router%nPut_P(iProcFrom) > nPutUbound_P(iProcFrom) .or. &
              sum(Router%nRecv_P(:)) > nRecvUbound

      end do IMAGES1
      if(.not.DoCountOnly)then
         PROCFROM:do iProcToGet = 1, nProcToGet
            iProcFrom = iProcLookUp_I(iProcToGet)
            Router%iPut_P(iProcFrom)%&
                 iCB_II(1:Router%nIndexTarget,&
                 Router%nPut_P(iProcFrom))&
                 = iIndexRecv_I(1:Router%nIndexTarget)
            Router%iPut_P(iProcFrom)%&
                 iCB_II(0,Router%nPut_P(iProcFrom)) = 1
            Router%Put_P(iProcFrom)%&
                 Weight_I(Router%nPut_P(iProcFrom)) = 1.0
            Router%DoAdd_P(iProcFrom)%&
                 DoAdd_I(Router%nRecv_P(iProcFrom)) = .true.

            ! indices of the location where data has to be put
            ! store processor id for later use
            nBuffer = nBuffer + 1
            iProc_I(nBuffer) = iProcFrom

            ! index of the image in the buffer FOR CURRENT
            ! source processor
            iOrder_I(nBuffer) = Router%nRecv_P(iProcFrom)

            ! fill the buffer to be sent
            ! coordinates on Source
            Router%BufferTarget_II(1:Router%nDim, nBuffer) = CoordSource_D
            if(nAux > 0)then
               if(present(extra_data))then
                  call extra_data(&
                       GridTarget%iIndex_IB(:,iBlockUsed), &
                       iPointInBlock,  nAux, Router%BufferTarget_II(&
                       Router%iAuxStart:Router%iAuxEnd,nBuffer))
               else
                  iBlockAll = GridTarget%iIndex_IB(&
                       GlobalBlock_,iBlockUsed)
                  iAux_I(iPointInBlock_) = iPointInBlock
                  iAux_I(iBlockAll_    ) = iBlockAll
                  Router%BufferTarget_II(Router%iAuxStart:Router%iAuxEnd,&
                       nBuffer) = real(iAux_I(1:nAux))
               end if
            end if
         end do PROCFROM

         ! DoAdd should be set to .false. for the same PE or for
         ! the minimal PE
         if(any(iProcLookUp_I(1:nProcToGet) == iProc))then
            iProcDoNotAdd = iProc
         else
            iProcDoNotAdd = minval(iProcLookUp_I(1:nProcToGet))
         end if
         Router%DoAdd_P(iProcDoNotAdd)%&
              DoAdd_I(Router%nRecv_P(iProcDoNotAdd)) = .false.
      end if
    end subroutine put_point_to_semirouter
    !==========================================================================
  end subroutine set_semi_router_from_target
  !============================================================================
  subroutine synchronize_router_from_target(Router)
    type(RouterType),        intent(inout):: Router

    integer :: nRecvCumSum, nSendCumSum, iBuffer

    ! MPI-related variables
    integer :: iStatusS_II(MPI_STATUS_SIZE, Router%nProc)
    integer :: iStatusR_II(MPI_STATUS_SIZE, Router%nProc)
    integer :: iRequestS_I(Router%nProc), iRequestR_I(Router%nProc)
    integer :: nRequestR, nRequestS, iError, iTag=0
    integer :: iProc, nProc, iProcFrom, iProcTo
    integer :: iProc0Source, iProcLastSource, iProcStrideSource
    integer :: iProc0Target, iProcLastTarget, iProcStrideTarget
    !--------------------------------------------------------------------------
    ! Return if the processor does not belong to the communicator
    if(.not.Router%IsProc) RETURN
    ! identify components
    ! For given PE the index in the communicator is:
    iProc = Router % iProc
    ! total number of processors and on components
    nProc = Router%nProc
    if(Router%IsLocal)then
       iProc0Source = 0;  iProcStrideSource = 1
       iProcLastSource = nProc-1
       iProc0Target = 0;  iProcStrideTarget = 1
       iProcLastTarget = nProc-1
    else
       iProc0Target      = i_proc0(Router%iCompTarget)
       iProcLastTarget   = i_proc_last(Router%iCompTarget)
       iProcStrideTarget = i_proc_stride(Router%iCompTarget)
       iProc0Source      = i_proc0(Router%iCompSource)
       iProcLastSource   = i_proc_last(Router%iCompSource)
       iProcStrideSource = i_proc_stride(Router%iCompSource)
    end if

    ! Stage 2:
    ! send the router info to Source:
    ! first, send the amount of data to be received,
    ! then the info itself (stored in buffer on Target)

    ! post recvs
    nRequestR = 0
    if(is_proc(Router%iCompSource))then
       do iProcFrom = iProc0Target, iProcLastTarget, iProcStrideTarget
          ! Do not wait for the message from self
          if(iProc==iProcFrom)CYCLE
          nRequestR = nRequestR + 1
          iTag = iProcFrom
          call MPI_Irecv(Router % nSend_P(iProcFrom), 1, MPI_INTEGER,&
               iProcFrom, iTag, Router%iComm, iRequestR_I(nRequestR), iError)
       end do
    end if
    ! Make sure all recv's are posted before using an rsend
    if(DoRSend)call MPI_BARRIER(Router%iCommUnion,  iError)

    ! post sends
    nRequestS = 0
    if(is_proc(Router%iCompTarget))then
       do iProcTo =  iProc0Source, iProcLastSource, iProcStrideSource
          ! Copy, if needed
          if(iProc==iProcTo)then
             Router % nSend_P(iProc) = Router % nRecv_P(iProc)
             CYCLE
          end if
          nRequestS = nRequestS + 1
          iTag = iProc
          if(DoRSend)then
             call MPI_Rsend(Router % nRecv_P(iProcTo), 1, MPI_INTEGER,&
                  iProcTo, iTag, Router%iComm, iError)
          else
             call MPI_Isend(Router % nRecv_P(iProcTo), 1, MPI_INTEGER,&
                  iProcTo, iTag, Router%iComm, iRequestS_I(nRequestS), iError)
          end if
       end do
    end if
    ! Finalize transfer
    call MPI_waitall(nRequestR, iRequestR_I, iStatusR_II, iError)
    if(.not.DoRSend)&
         call MPI_waitall(nRequestS, iRequestS_I, iStatusS_II, iError)

    ! send the actual router info
    ! post recvs
    nRequestR = 0
    if(is_proc(Router%iCompSource))then
       ! Allocate source buffer
       call check_router_allocation(Router)
       nRecvCumSum = Router%nSend_P(iProc)
       do iProcFrom = iProc0Target, iProcLastTarget, iProcStrideTarget
          if(Router%nSend_P(iProcFrom) == 0)CYCLE
          ! Do not wait for the message from self
          if(iProc==iProcFrom)CYCLE
          nRequestR = nRequestR + 1
          iTag = iProcFrom
          call MPI_Irecv(Router%BufferSource_II(1:Router%nVar,&
               1+nRecvCumSum:Router%nSend_P(iProcFrom)+nRecvCumSum), &
               Router%nSend_P(iProcFrom)*Router%nVar, MPI_REAL,&
               iProcFrom, iTag, Router%iComm, iRequestR_I(nRequestR), iError)
          nRecvCumSum = nRecvCumSum + Router%nSend_P(iProcFrom)
       end do
    end if
    ! Make sure all recv's are posted before using an rsend
    if(DoRSend)call MPI_BARRIER(Router%iCommUnion,  iError)
    ! post sends
    nRequestS = 0
    if(is_proc(Router%iCompTarget))then
       ! Copy for the data exchange on the same PE
       if(Router%nRecv_P(iProc) > 0)&
            Router%BufferSource_II(:,1:Router%nSend_P(iProc)) = &
            Router%BufferTarget_II(:,1:Router%nRecv_P(iProc))
       nSendCumSum = Router % nRecv_P(iProc)
       do iProcTo = iProc0Source, iProcLastSource, iProcStrideSource
          if(Router % nRecv_P(iProcTo) == 0) CYCLE
          if(iProcTo==iProc)CYCLE
          nRequestS = nRequestS + 1
          iTag = iProc
          if(DoRSend)then
             call MPI_Rsend(Router%BufferTarget_II(1:Router%nVar,&
                  1+nSendCumSum:Router%nRecv_P(iProcTo)+nSendCumSum), &
                  Router%nRecv_P(iProcTo)*Router%nVar, MPI_REAL,&
                  iProcTo, iTag, Router%iComm, iError)
          else
             call MPI_Isend(Router%BufferTarget_II(1:Router%nVar,&
                  1+nSendCumSum:Router%nRecv_P(iProcTo)+nSendCumSum), &
                  Router%nRecv_P(iProcTo)*Router%nVar, MPI_REAL,&
                  iProcTo, iTag, Router%iComm, iRequestS_I(nRequestS), iError)
          end if
          nSendCumSum = nSendCumSum + Router%nRecv_P(iProcTo)
       end do
    end if
    ! Finalize transfer
    call MPI_waitall(nRequestR, iRequestR_I, iStatusR_II, iError)
    if(.not.DoRSend)&
         call MPI_waitall(nRequestS, iRequestS_I, iStatusS_II, iError)
  end subroutine synchronize_router_from_target
  !============================================================================
  subroutine update_semi_router_at_source(Router, GridSource, interpolate)
    integer :: iCompTarget
    type(RouterType),  intent(inout):: Router
    type(GridType),    intent(in)   :: GridSource
    interface
       subroutine interpolate(nDim, Coord_D, Grid, &
            nIndex, iIndex_II, nImage, Weight_I)
         ! interpolation on Source's grid;
         ! provides PE and indices to access images of
         ! data location and interpolation weights
         use CON_grid_descriptor
         implicit none
         ! number of indices per data entry
         integer, intent(in):: nDim
         ! data location on Source
         real,    intent(inout):: Coord_D(nDim)
         ! grid descriptor
         type(GridType):: Grid
         ! indices of images, their number and interpolation weights
         integer, intent(in) :: nIndex
         integer, intent(out):: iIndex_II(0:nIndex,2**nDim)
         integer, intent(out):: nImage
         real,    intent(out):: Weight_I(2**nDim)
       end subroutine interpolate
    end interface

    optional:: interpolate
    integer :: nRecvCumSum
    integer :: iProcTo, iBuffer
    integer :: iProc
    integer :: nDimSource, nIndexSource, nImageMax
    logical :: DoInterpolate
    integer :: iProc0Target, iProcLastTarget, iProcStrideTarget
    !--------------------------------------------------------------------------
    ! Return if the processor does not belong to the communicator
    if(.not.Router%IsProc) RETURN
    if(.not.is_proc(Router%iCompSource))RETURN
    ! For given PE the index in the communicator is:
    iProc = Router % iProc
    ! identify components
    iCompTarget = Router%iCompTarget
    if(Router%IsLocal)then
       iProc0Target = 0;  iProcStrideTarget = 1
       iProcLastTarget = Router%nProc-1
    else
       iProc0Target      = i_proc0(iCompTarget)
       iProcLastTarget   = i_proc_last(iCompTarget)
       iProcStrideTarget = i_proc_stride(iCompTarget)
    end if
    DoInterpolate = present(interpolate)

    ! Stage 3 set semi-router for source
    ! process the data that has been received
    nDimSource  = GridSource%nDim
    nIndexSource= Router%nIndexSource

    nImageMax = 2**nDimSource
    ! prepare containers for router information of Source side
    do iProcTo = iProc0Target, iProcLastTarget, iProcStrideTarget
       Router%nGet_P(iProcTo) = Router%nSend_P(iProcTo)*nImageMax
    end do
    call check_router_allocation(Router)
    Router%nGet_P = 0
    ! fill these containers

    ! First process the buffered info from this PE
    iProcTo = iProc
    do iBuffer = 1, Router%nSend_P(iProc)
       call put_point_to_semirouter
    end do
    nRecvCumSum = Router%nSend_P(iProc)
    do iProcTo =  iProc0Target, iProcLastTarget, iProcStrideTarget
       if(iProcTo == iProc)CYCLE
       do iBuffer = nRecvCumSum + 1, nRecvCumSum + Router%nSend_P(iProcTo)
          call put_point_to_semirouter
       end do
       ! increment the offset
       nRecvCumSum = nRecvCumSum + Router%nSend_P(iProcTo)
    end do
  contains
    !==========================================================================
    subroutine put_point_to_semirouter
      ! interpolation-related variables
      integer :: iIndexGet_II(0:GridSource%nDim+1, 2**GridSource%nDim)
      real    :: CoordSource_D(GridSource%nDim)
      real    :: CoordPass_D(GridSource%nDim)
      real    :: Weight_I(2**GridSource%nDim)
      integer :: iImage, nImage, nImagePart, iToGet
      !------------------------------------------------------------------------
      CoordSource_D = Router%BufferSource_II(1:Router%nDim, iBuffer)
      CoordPass_D = CoordSource_D
      if( DoInterpolate)then
         call interpolate(  &
              nDimSource,   &
              CoordPass_D,    &
              GridSource,     &
              nIndexSource, &
              iIndexGet_II, &
              nImage,       &
              Weight_I)
      else
         call nearest_grid_points(&
              nDim           = nDimSource, &
              Coord_D          = CoordPass_D, &
              Grid             = GridSource, &
              nIndex         = nIndexSource, &
              iIndex_II      = iIndexGet_II, &
              nImage         = nImage, &
              Weight_I       = Weight_I)
      end if
      if(nImage<1)then
         Write(*,*)'Interpolation in CON_router failed on Source PE=',iProc
         write(*,*)'Point Coords=', CoordSource_D, ' nImage=',nImage
         call CON_stop('Failure!!!')
      end if
      nImagePart = count(iIndexGet_II(0,1:nImage)==iProc)
      Router%nGet_P(iProcTo) = Router%nGet_P(iProcTo) + nImagePart
      if(nImagePart==0)then
         Write(*,*)'Interpolation in CON_router failed on Source PE=', iProc
         write(*,*)'Point Coords=', CoordSource_D,  ' nImage=',nImage
         write(*,*)' iProc  iBlock iCell'
         do iImage = 1, nImage
            write(*,*) iIndexGet_II(:, iImage)
         end do
         call CON_stop('No image on the requested PE')
      end if
      ! indices
      do iImage=1,nImage
         if(iProc==iIndexGet_II(0,iImage))then
            iToGet = Router%nGet_P(iProcTo)+1-nImagePart
            Router%iGet_P(iProcTo)%iCB_II(:,iToGet) = iIndexGet_II(:,iImage)
            Router%iGet_P(iProcTo)%iCB_II(0,iToGet) = nImagePart
            Router%Get_P(iProcTo)%Weight_I(iToGet)  = Weight_I(iImage)
            nImagePart = nImagePart - 1
         end if
      end do
    end subroutine put_point_to_semirouter
    !==========================================================================
  end subroutine update_semi_router_at_source
  !============================================================================
  subroutine construct_router_from_source(GridSource, GridTarget, Router,     &
       is_interface_block, n_interface_point_in_block, interface_point_coords,&
       mapping, interpolate)

    ! GridSource: the Grid for the Source component
    ! GridTarget: the Grid for the Target component
    ! Router: the router to be set
    ! is_interface_block:
    ! Logical function which allows to skip the block if there is no !
    ! interface points in it. Optional, if not present then all the  !
    ! blocks are checked for the presence of the interface points    !

    ! n_interface_point_in_block:
    ! The function may be used for sparse 1D grids
    ! in which the real information is available only in the first
    ! n grid points, for the given plocks, the other points being
    ! unused. For this purpose, the function should return this
    ! number of used points. It is not applied, if the rurned value
    ! is larger than the total number of grid points, in the given
    ! plock.

    ! interface_point_coords:
    ! The subroutine which defines if the grid point is inside the
    ! interface layer. Optional, if not present, then all the grid
    ! points (at the source grid) are considered as the interface
    ! layer points

    ! mapping: transformation of the location coordinates between components

    ! interpolate: interpolation subroutine for the Target's grid

    type(LocalGridType),        intent(in) :: GridSource
    type(GridType),             intent(in) :: GridTarget
    type(RouterType),        intent(inout) :: Router
    interface
       logical function is_interface_block(iBlockLocal)
         implicit none
         integer,intent(in)::iBlockLocal
       end function is_interface_block
       integer function n_interface_point_in_block(iBlockLocal)

         ! To be applicable for the original purposes, the function
         ! should return 0 (or less) for local blocks in which there is
         ! no interface points. It may be also used for sparse 1D grids
         ! in which the real information is available only in the first
         ! n grid points, for the given plocks, the other points being
         ! unused. For this purpose, the function should return this
         ! number of used points. It is not applied, if the rurned value
         ! is larger than the total number of grid points, in the given
         ! plock.
         implicit none
         integer,intent(in) :: iBlockLocal
       end function n_interface_point_in_block
       subroutine interface_point_coords( &
            nDim, Xyz_D, nIndex, iIndex_I,&
            IsInterfacePoint)
         implicit none
         integer,intent(in)    :: nDim, nIndex
         logical,intent(out)   :: IsInterfacePoint
         real,   intent(inout) :: Xyz_D(nDim)
         integer,intent(inout) :: iIndex_I(nIndex)
       end subroutine interface_point_coords
       subroutine mapping(nDimIn, CoordIn_D, nDimOut, CoordOut_D, &
            IsInterfacePoint)
         ! this subroutine mapss coordinates between components
         integer, intent(in)  :: nDimIn
         real,    intent(in)  :: CoordIn_D(nDimIn)
         integer, intent(in)  :: nDimOut
         real,    intent(out) :: CoordOut_D(nDimOut)
         logical, intent(out) :: IsInterfacePoint
       end subroutine mapping
       subroutine interpolate(&
            nDim, Coord_D, Grid, &
            nIndex, iIndex_II, nImage, Weight_I)
         ! interpolation on Source's grid;
         ! provides PE and indices to access images of
         ! data location and interpolation weights
         use CON_grid_descriptor
         implicit none
         ! number of indices per data entry
         integer, intent(in)     :: nDim
         ! data location on Source
         real,    intent(inout)  :: Coord_D(nDim)
         ! grid descriptor
         type(GridType):: Grid
         ! indices of images, their number and interpolation weights
         integer, intent(in)     :: nIndex
         integer, intent(out)    :: iIndex_II(0:nIndex,2**nDim)
         integer, intent(out)    :: nImage
         real,    intent(out)    :: Weight_I(2**nDim)
       end subroutine interpolate
    end interface
    optional :: is_interface_block
    optional :: n_interface_point_in_block
    optional :: interface_point_coords
    optional :: mapping
    optional :: interpolate
    !--------------------------------------------------------------------------
    ! Return if the processor does not belong to the communicator
    if(.not.Router%IsProc)RETURN

    call set_router_dim(Router, GridTarget%nDim)
    if(is_proc(Router%iCompSource))&
         call set_semi_router_from_source(GridSource, GridTarget, Router,  &
         is_interface_block, n_interface_point_in_block,                   &
         interface_point_coords, mapping, interpolate)
    call synchronize_router_from_source(Router)
    if(is_proc(Router%iCompTarget))call update_semi_router_at_target(      &
         Router, GridTarget, interpolate)
  end subroutine construct_router_from_source
  !============================================================================
  subroutine set_semi_router_from_source(GridSource, GridTarget, Router,   &
       is_interface_block, n_interface_point_in_block,                     &
       interface_point_coords, mapping, interpolate)

    ! GridSource: the Grid for the Source component
    ! GridTarget: the Grid for the Target component
    ! Router: the router to be set
    ! is_interface_block:
    ! Logical function which allows to skip the block if there is no !
    ! interface points in it. Optional, if not present then all the  !
    ! blocks are checked for the presence of the interface points    !

    ! n_interface_point_in_block:
    ! The function may be used for sparse 1D grids
    ! in which the real information is available only in the first
    ! n grid points, for the given plocks, the other points being
    ! unused. For this purpose, the function should return this
    ! number of used points. It is not applied, if the rurned value
    ! is larger than the total number of grid points, in the given
    ! plock.

    ! interface_point_coords:
    ! The subroutine which defines if the grid point is inside the
    ! interface layer. Optional, if not present, then all the grid
    ! points (at the source grid) are considered as the interface
    ! layer points

    ! mapping: transformation of the location coordinates between components

    ! interpolate: interpolation subroutine for the Target's grid

    type(LocalGridType),   intent(in) :: GridSource
    type(GridType),        intent(in) :: GridTarget
    type(RouterType),   intent(inout) :: Router
    interface
       logical function is_interface_block(iBlockLocal)
         implicit none
         integer,intent(in)::iBlockLocal
       end function is_interface_block
       integer function n_interface_point_in_block(iBlockLocal)
         ! To be applicable for the original purposes, the function
         ! should return 0 (or less) for local blocks in which there is
         ! no interface points. It may be also used for sparse 1D grids
         ! in which the real information is available only in the first
         ! n grid points, for the given plocks, the other points being
         ! unused. For this purpose, the function should return this
         ! number of used points. It is not applied, if the rurned value
         ! is larger than the total number of grid points, in the given
         ! plock.
         implicit none
         integer,intent(in) :: iBlockLocal
       end function n_interface_point_in_block
       subroutine interface_point_coords(nDim, Coord_D, nIndex, iIndex_I, &
            IsInterfacePoint)
         implicit none
         integer,intent(in)    :: nDim, nIndex
         logical,intent(out)   :: IsInterfacePoint
         real,   intent(inout) :: Coord_D(nDim)
         integer,intent(inout) :: iIndex_I(nIndex)
       end subroutine interface_point_coords
       subroutine mapping(nDimIn, CoordIn_D, nDimOut, CoordOut_D, &
            IsInterfacePoint)
         ! this subroutine mapss coordinates between components
         integer, intent(in)  :: nDimIn
         real,    intent(in)  :: CoordIn_D(nDimIn)
         integer, intent(in)  :: nDimOut
         real,    intent(out) :: CoordOut_D(nDimOut)
         logical, intent(out) :: IsInterfacePoint
       end subroutine mapping
       subroutine interpolate(nDim, Coord_D, Grid, nIndex, iIndex_II,&
            nImage, Weight_I)
         ! interpolation on Source's grid;
         ! provides PE and indices to access images of
         ! data location and interpolation weights
         use CON_grid_descriptor
         implicit none
         ! number of indices per data entry
         integer, intent(in)     :: nDim
         ! data location on Source
         real,    intent(inout)  :: Coord_D(nDim)
         ! grid descriptor
         type(GridType):: Grid
         ! indices of images, their number and interpolation weights
         integer, intent(in)     :: nIndex
         integer, intent(out)    :: iIndex_II(0:nIndex,2**nDim)
         integer, intent(out)    :: nImage
         real,    intent(out)    :: Weight_I(2**nDim)
       end subroutine interpolate
    end interface
    optional :: is_interface_block
    optional :: n_interface_point_in_block
    optional :: interface_point_coords
    optional :: mapping
    optional :: interpolate

    ! The presence of optional parameters
    logical :: UseMappingFunction, DoCheckBlock, DoCheckBlockSize, &
         DoCheckPoint ,DoInterpolate
    ! MPI-related variables
    integer :: iProc
    logical :: DoCountOnly
    integer :: iProcTo, iBlockFrom, iLocalGridPoint
    integer :: iLocalPointLast, iPointInBlock, nInterfacePoint
    logical :: IsInterfacePoint
    integer :: nGetUbound_P(0:Router%nProc-1)
    integer :: nSendUbound

    real    :: CoordTarget_D(GridTarget%nDim)
    real    :: CoordSource_D(GridSource%nDim)
    integer :: iCell_D(GridSource%nDim)
    integer :: iIndexSend_I(Router%nIndexSource)
    ! components ids
    integer :: iCompTarget, iProc0Target, iProcLastTarget, &
         iProcStrideTarget
    ! dimensionality of components
    integer :: nDimSource, nDimTarget
    ! number of indices (e.g. cell indices) on components
    integer :: nIndexSource, nIndexTarget
    ! loop variables
    integer :: nBuffer, iBlockUsed
    ! offset in buffer
    integer :: nSendCumSum

    ! Misc
    integer :: iBlockMisc,  nAux
    character(len=*), parameter:: NameSub = 'set_semi_router_from_source'
    !--------------------------------------------------------------------------
    ! Return if the processor does not belong to the communicator
    if(.not.Router%IsProc) RETURN
    ! Stage 1:
    ! on Source PE determine the sources of data and
    ! how much data is to be sent from them
    if(.not.is_proc(Router%iCompSource))RETURN
    ! For given PE the index in the communicator is:
    iProc = Router % iProc
    ! identify components
    iCompTarget = Router%iCompTarget
    if(Router%IsLocal)then
       iProc0Target = 0;  iProcStrideTarget = 1
       iProcLastTarget = Router%nProc - 1
    else
       iProc0Target      = i_proc0(      iCompTarget)
       iProcLastTarget   = i_proc_last(  iCompTarget)
       iProcStrideTarget = i_proc_stride(iCompTarget)
    end if
    ! determine which optional actions should be taken
    UseMappingFunction = present(mapping)
    DoInterpolate      = present(interpolate)
    DoCheckBlock       = present(is_interface_block)
    DoCheckBlockSize   = present(n_interface_point_in_block)
    DoCheckPoint       = present(interface_point_coords)

    ! introduced for a better readability
    nDimTarget   = GridTarget%nDim
    nDimSource   = GridSource%nDim
    nIndexTarget = Router%nIndexTarget
    nIndexSource = Router%nIndexSource
    if(.not.UseMappingFunction .and. nDimTarget /= nDimSource )call CON_stop(&
         NameSub//': Mapping is needed for Target%nDim/=Source%nDim')

    ! some data will be sent to Source, determine amount:
    ! cell and block indexes are sent
    nAux = Router%nMappedPointIndex
    ! Check dimensions
    DoCountOnly=.true. ! To enter the loop
    do while(DoCountOnly)
       call check_router_allocation(Router)
       ! Store Upper bounds to control if the alllocated
       ! index arrays have sufficient size
       do iProcTo = iProc0Target, iProcLastTarget, iProcStrideTarget
          nGetUbound_P(iProcTo) = ubound(Router%iGet_P(iProcTo)%iCB_II,2)
       end do
       ! and also the buffer size
       nSendUbound = ubound(Router%BufferSource_II, 2)

       ! which processor holds a current image
       call check_size(1, [Router%nBufferSource], iBuffer_I = iProc_I)
       ! correct order of images in the send buffer
       call check_size(1, [Router%nBufferSource], iBuffer_I = iOrder_I)
       ! reset basic coupling information
       Router%nGet_P  = 0
       Router%nSend_P = 0
       ! reset index of a current data entry in the buffer
       nBuffer        = 0
       DoCountOnly    = .false.

       ! If the check shows that the allocated array is not
       ! sufficient, then DoCountOnly will be set to true. The loop
       ! then will be repeated for the second time

       ! Loop over local blocks
       BLOCKS: do iBlockUsed = 1, GridSource%nBlock
          ! Convert to the local block number
          iBlockFrom = GridSource%iIndex_IB(BLK_,iBlockUsed)
          ! Skip the block if desired: if there is known to be no interface!
          ! point in it                                                    !
          if( DoCheckBlock)then
             if(.not.is_interface_block(iBlockFrom))CYCLE BLOCKS
          end if
          ! Prepare a GlobalCellNumber Loop, for a given (octree) block
          ! 1. Global number of the last grid point in the previous
          ! global block, with the global block number = nBlockAll-1
          iLocalPointLast = GridSource%iIndex_IB(&
               GridPointFirst_,iBlockUsed) - 1
          nInterfacePoint = GridSource%nPointPerBlock
          if( DoCheckBlockSize)then
             nInterfacePoint = min(nInterfacePoint,&
                  n_interface_point_in_block(iBlockFrom))
          end if

          ! Now the subloop runs over the iPointBlock within
          ! the current local block
          POINTS:do iPointInBlock  = 1, nInterfacePoint

             ! Recover the global grid point number
             iLocalGridPoint =   iLocalPointLast + iPointInBlock

             ! Get cell index for this grid point in source.
             call global_i_grid_point_to_icb(GridSource,   &
                  iLocalGridPoint,&
                  iBlockMisc,     & !<=The output is iBlockUsed, which is
                  iCell_D)          ! the loop variable, cannot be affected

             ! Shape an index array for this grid point in source
             iIndexSend_I(nIndexSource) = iBlockFrom
             iIndexSend_I(1:nDimSource) = iCell_D

             ! Generalized coordinates of the Source
             ! grid point
             CoordSource_D = coord_grid_d(GridSource, iBlockUsed, iCell_D)

             ! Using the grid descriptor for the Source we may want
             ! 1. To select or deselect this point from the interface
             ! 2. To convert (with no change in the dimensions) the
             ! grid point coordinates (say, from generalized to
             ! Cartesian), that is why CoordSource_D has intent inout
             if( DoCheckPoint)then
                call interface_point_coords(&
                     nDimSource            ,&
                     CoordSource_D         ,&
                     nIndexSource          ,&
                     iIndexSend_I          ,&
                     IsInterfacePoint)
                if(.not.IsInterfacePoint)CYCLE POINTS

                ! Otherwise, no point is deselected,
                ! the generalized coordinates keep unchanged
             end if
             if(UseMappingFunction)then

                ! This is a mapping of the
                ! source point coordinates to the generalized
                ! coordinates of target
                call mapping(&
                     nDimSource,&
                     CoordSource_D,&
                     nDimTarget,&
                     CoordTarget_D, &
                     IsInterfacePoint)
                if(.not.IsInterfacePoint)CYCLE POINTS
             else

                ! Otherwise, the coordinates in the source
                ! are the same as those in the oputput from
                !
                CoordTarget_D=CoordSource_D
             end if

             ! Now, we put the point indexes to the router
             ! and its coordinates (and, probably,
             ! its global point number) to the buffer
             call put_point_to_semirouter
          end do POINTS   ! iGlobalPoints
       end do BLOCKS      ! iGlobalBlock
    end do                ! DoCountOnly
    ! fix the order of Buffer_I so contiguous chunks of data can be sent
    ! to the appropriate processors of Source,
    ! currently iOrder_I contains indices WITHIN these chunks
    nSendCumSum = Router%nSend_P(iProc)
    do iProcTo = iProc0Target, iProcLastTarget, iProcStrideTarget
       if(iProc==iProcTo)CYCLE
       where(iProc_I(1:nBuffer) == iProcTo)&
            iOrder_I( 1:nBuffer) = iOrder_I( 1:nBuffer) + nSendCumSum
       nSendCumSum = nSendCumSum + Router%nSend_P(iProcTo)
    end do

    ! the correct order is found, apply it
    Router%BufferSource_II(:,iOrder_I(1:nBuffer)) = &
         Router%BufferSource_II(:,1:nBuffer)
  contains
    !==========================================================================
    subroutine put_point_to_semirouter
      integer :: iIndexPut_II(0:Router%nIndexTarget,&
           2**GridTarget%nDim)
      real    :: Weight_I(2**GridTarget%nDim)
      real    :: CoordPass_D(GridTarget%nDim)
      integer :: iImage, nImage!, iProcDoNotAdd
      integer :: iProcLookUp_I(2**GridTarget%nDim)
      integer :: nProcToPut, iProcToPut
      integer, parameter:: iPointGlobal_ = 0,  &
           iPointInBlock_ = 2, iBlockAll_ = 1
      integer :: iAux_I(1:2)
      ! error message containers
      character(len=200):: &
           StringErrorFormat, StringErrorMessage
      !------------------------------------------------------------------------
      CoordPass_D = CoordTarget_D
      if( DoInterpolate)then
         call interpolate(  &
              nDimTarget,   &
              CoordPass_D,    &
              GridTarget,     &
              nIndexTarget, &
              iIndexPut_II, &
              nImage,       &
              Weight_I)
      else
         call nearest_grid_points(&
              nDim           = nDimTarget, &
              Coord_D          = CoordPass_D, &
              Grid             = GridTarget, &
              nIndex         = nIndexTarget, &
              iIndex_II      = iIndexPut_II, &
              nImage         = nImage, &
              Weight_I       = Weight_I)
      end if
      if(nImage < 1)then
         write(StringErrorFormat,'(a,i3,a)') '(a,',nDimTarget,'es15.7)'
         write(StringErrorMessage,StringErrorFormat)&
              NameSub//': Interpolation failed at location ', &
              CoordTarget_D
         call CON_stop(StringErrorMessage)
      end if

      ! go over the list of images and process result of interpolation
      IMAGES1:do iImage = 1, nImage
         iProcTo = iIndexPut_II(0, iImage)
         if(iImage==1)then
            iProcLookUp_I(1)=iProcTo
            nProcToPut=1
         else
            if(.not.any(iProcLookUp_I(1:nProcToPut)==iProcTo))then
               nProcToPut               = nProcToPut + 1
               iProcLookUp_I(nProcToPut) = iProcTo
            else
               CYCLE IMAGES1
            end if
         end if
         ! index of current data entry in the buffer
         Router%nGet_P(iProcTo)  = Router%nGet_P(iProcTo ) + 1
         Router%nSend_P(iProcTo) = Router%nSend_P(iProcTo) + 1

         DoCountOnly = DoCountOnly.or.&
              Router%nGet_P(iProcTo) > nGetUbound_P(iProcTo).or.&
              sum(Router%nSend_P(:)) > nSendUbound
      end do IMAGES1
      if(.not.DoCountOnly)then
         PROCTO:do iProcToPut = 1, nProcToPut
            iProcTo = iProcLookUp_I(iProcToPut)
            Router%iGet_P(iProcTo)%&
                 iCB_II(1:Router%nIndexSource, Router%nGet_P(iProcTo))&
                 = iIndexSend_I(1:Router%nIndexSource)
            Router%iGet_P(iProcTo)%iCB_II(0,Router%nGet_P(iProcTo)) = 1
            Router%Get_P(iProcTo)%Weight_I(Router%nGet_P(iProcTo))  = 1.0
            ! indices of the location where data has to be put
            ! store processor id for later use
            nBuffer = nBuffer + 1
            iProc_I(nBuffer) = iProcTo
            ! index of the image in the buffer FOR CURRENT
            ! source processor
            iOrder_I(nBuffer) = Router%nSend_P(iProcTo)

            ! fill the buffer to be sent
            ! coordinates on Source
            Router%BufferSource_II(&
                 1:Router%nDim, nBuffer)=&
                 CoordTarget_D
            if(nAux > 0)then
               iAux_I(iPointInBlock_) = iPointInBlock
               iAux_I(iBlockAll_    ) = GridSource%iIndex_IB(&
                    GlobalBlock_,iBlockUsed)
               Router%BufferSource_II(Router%iAuxStart:Router%iAuxEnd,&
                    nBuffer) = real(iAux_I(1:nAux))
            end if
         end do PROCTO
      end if
    end subroutine put_point_to_semirouter
    !==========================================================================
  end subroutine set_semi_router_from_source
  !============================================================================
  subroutine synchronize_router_from_source(Router)
    type(RouterType), intent(inout):: Router
    integer:: nRecvCumSum, nRecvCumSumMy, nSendCumSum

    ! MPI-related variables
    integer :: iStatusS_II(MPI_STATUS_SIZE, Router%nProc)
    integer :: iStatusR_II(MPI_STATUS_SIZE, Router%nProc)
    integer :: iRequestS_I(Router%nProc), iRequestR_I(Router%nProc)
    integer :: nRequestR, nRequestS, iError, iTag=0
    integer :: iProc, nProc
    integer :: iProcFrom, iProcTo
    integer :: iCompTarget,  iCompSource
    integer :: iProc0Source, iProcLastSource, iProcStrideSource
    integer :: iProc0Target, iProcLastTarget, iProcStrideTarget
    character(len=*), parameter:: NameSub = 'synchronize_router_from_source'
    !--------------------------------------------------------------------------
    ! Return if the processor does not belong to the communicator
    if(.not.Router%IsProc) RETURN
    ! For given PE the index in the communicator is:
    iProc = Router % iProc
    ! total number of processors and on components
    nProc       = Router%nProc
    iCompTarget = Router%iCompTarget
    iCompSource = Router%iCompSource
    if(Router%IsLocal)then
       iProc0Source = 0;  iProcStrideSource = 1
       iProcLastSource = nProc-1
       iProc0Target = 0;  iProcStrideTarget = 1
       iProcLastTarget = nProc-1
    else
       iProc0Target      = i_proc0(iCompTarget)
       iProcLastTarget   = i_proc_last(iCompTarget)
       iProcStrideTarget = i_proc_stride(iCompTarget)
       iProc0Source      = i_proc0(iCompSource)
       iProcLastSource   = i_proc_last(iCompSource)
       iProcStrideSource = i_proc_stride(iCompSource)
    end if
    ! Stage 2:
    ! send the router info to Target:
    ! first, send the amount of data to be received,
    ! then the info itself (stored in buffer on Source)

    ! post recvs
    nRequestR = 0
    if(is_proc(iCompTarget))then
       do iProcFrom = iProc0Source, iProcLastSource, iProcStrideSource
          ! Do not expect the message from self
          if(iProcFrom==iProc)CYCLE
          nRequestR = nRequestR + 1
          iTag = iProcFrom
          call MPI_Irecv(Router % nRecv_P(iProcFrom), 1, MPI_INTEGER,       &
               iProcFrom, iTag, Router%iComm, iRequestR_I(nRequestR), iError)
       end do
    end if
    ! Make sure all recv's are posted before using an rsend
    if(DoRSend)call MPI_BARRIER(Router%iCommUnion,  iError)
    ! post sends
    nRequestS = 0
    if(is_proc(iCompSource))then
       do iProcTo = iProc0Target, iProcLastTarget, iProcStrideTarget
          if(iProc==iProcTo)then
             ! Copy, if needed
             Router % nRecv_P(iProc) = Router % nSend_P(iProc)
             CYCLE
          end if
          nRequestS = nRequestS + 1
          iTag = iProc
          if(DoRSend)then
             call MPI_Rsend(Router % nSend_P(iProcTo), 1, MPI_INTEGER,&
                  iProcTo, iTag, Router%iComm, iError)
          else
             call MPI_Isend(Router % nSend_P(iProcTo), 1, MPI_INTEGER,&
                  iProcTo, iTag, Router%iComm, iRequestS_I(nRequestS), iError)
          end if
       end do
    end if
    call MPI_waitall(nRequestR, iRequestR_I, iStatusR_II, iError)
    if(.not.DoRSend)&
         call MPI_waitall(nRequestS, iRequestS_I, iStatusS_II, iError)
    ! send the actual router info
    ! post recvs
    nRequestR = 0
    if(is_proc(iCompTarget))then
       ! Allocate target buffer
       call check_router_allocation(Router)
       nRecvCumSum =  Router%nRecv_P(iProc)
       do iProcFrom = iProc0Source, iProcLastSource, iProcStrideSource
          if(Router%nRecv_P(iProcFrom) == 0)CYCLE
          ! Do not expect the message from self
          if(iProc==iProcFrom)CYCLE
          nRequestR = nRequestR + 1
          iTag = iProcFrom
          call MPI_Irecv(Router%BufferTarget_II(1:Router%nVar,         &
               1+nRecvCumSum:Router%nRecv_P(iProcFrom) + nRecvCumSum), &
               Router%nRecv_P(iProcFrom)*Router%nVar, MPI_REAL,        &
               iProcFrom, iTag, Router%iComm, iRequestR_I(nRequestR), iError)
          nRecvCumSum = nRecvCumSum + Router%nRecv_P(iProcFrom)
       end do
    end if  ! is_proc(iCompTarget)
    ! Make sure all recv's are posted before using an rsend
    if(DoRSend)call MPI_BARRIER(Router%iCommUnion,  iError)
    ! post sends
    nRequestS = 0
    if(is_proc(iCompSource))then
       ! Copy, if needed
       if(Router%nRecv_P(iProc) > 0)&
            Router%BufferTarget_II(:,1:Router%nRecv_P(iProc)) = &
            Router%BufferSource_II(:,1:Router%nSend_P(iProc))
       nSendCumSum = Router%nSend_P(iProc)
       do iProcTo = iProc0Target, iProcLastTarget, iProcStrideTarget
          if(Router % nSend_P(iProcTo) == 0) CYCLE
          if(iProc==iProcTo)CYCLE
          nRequestS = nRequestS + 1
          iTag = iProc
          if(DoRSend)then
             call MPI_Rsend(Router%BufferSource_II(1:Router%nVar,         &
                  1 + nSendCumSum:Router%nSend_P(iProcTo) + nSendCumSum), &
                  Router%nSend_P(iProcTo)*Router%nVar, MPI_REAL,          &
                  iProcTo, iTag, Router%iComm, iError)
          else
             call MPI_Isend(Router%BufferSource_II(1:Router%nVar,         &
                  1 + nSendCumSum:Router%nSend_P(iProcTo) + nSendCumSum), &
                  Router%nSend_P(iProcTo)*Router%nVar, MPI_REAL,          &
                  iProcTo, iTag, Router%iComm, iRequestS_I(nRequestS), iError)
          end if
          nSendCumSum = nSendCumSum + Router%nSend_P(iProcTo)
       end do
    end if
    ! Finalize transfer
    call MPI_waitall(nRequestR, iRequestR_I, iStatusR_II, iError)
    if(.not.DoRSend)&
         call MPI_waitall(nRequestS, iRequestS_I, iStatusS_II, iError)
  end subroutine synchronize_router_from_source
  !============================================================================
  subroutine update_semi_router_at_target(Router, GridTarget, interpolate)
    type(RouterType),   intent(inout):: Router
    type(GridType),     intent(in)   :: GridTarget
    interface
       subroutine interpolate(&
            nDim, Coord_D, Grid, &
            nIndex, iIndex_II, nImage, Weight_I)
         ! interpolation on Target's grid;
         ! provides PE and indices to access images of
         ! data location and interpolation weights
         use CON_grid_descriptor
         implicit none
         ! number of indices per data entry
         integer, intent(in):: nDim
         ! data location on Source
         real,    intent(inout):: Coord_D(nDim)
         ! grid descriptor
         type(GridType):: Grid
         ! indices of images, their number and interpolation weights
         integer, intent(in) :: nIndex
         integer, intent(out):: iIndex_II(0:nIndex,2**nDim)
         integer, intent(out):: nImage
         real,    intent(out):: Weight_I(2**nDim)
       end subroutine interpolate
    end interface
    optional:: interpolate
    integer :: nRecvCumSum
    integer :: iBuffer, iProcFrom
    integer :: iProc, nProc
    integer :: nDimTarget, nIndexTarget, nImageMax
    logical :: DoInterpolate
    integer :: iCompSource
    integer :: iProc0Source, iProcLastSource, iProcStrideSource
    ! Return if the processor does not belong to the communicator
    character(len=*), parameter:: NameSub = 'update_semi_router_at_target'
    !--------------------------------------------------------------------------
    if(.not.Router%IsProc) RETURN
    ! Return if the processor does not belong to the target
    if(.not.is_proc(Router%iCompTarget))RETURN
    ! identify components
    iCompSource = Router%iCompSource
    ! For given PE the index in the communicator is:
    iProc       = Router%iProc
    ! total number of processors and on components
    nProc       = Router%nProc
    if(Router%IsLocal)then
       iProc0Source = 0;  iProcStrideSource = 1
       iProcLastSource = nProc-1
    else
       iProc0Source      = i_proc0(iCompSource)
       iProcLastSource   = i_proc_last(iCompSource)
       iProcStrideSource = i_proc_stride(iCompSource)
    end if

    DoInterpolate = present(interpolate)

    nDimTarget       = GridTarget%nDim
    nIndexTarget     = Router%nIndexTarget

    nImageMax = 2**nDimTarget
    ! prepare containers for router information of Source side
    do iProcFrom = iProc0Source, iProcLastSource, iProcStrideSource
       Router%nPut_P(iProcFrom) = Router%nRecv_P(iProcFrom)*nImageMax
    end do
    call check_router_allocation(Router)
    Router%nPut_P = 0
    ! fill these containers

     ! First process the buffered info from this PE
    iProcFrom = iProc
    do iBuffer = 1, Router%nRecv_P(iProc)
       call put_point_to_semirouter
    end do
    nRecvCumSum = Router%nRecv_P(iProc)
    do iProcFrom = iProc0Source, iProcLastSource, iProcStrideSource
       if(iProcFrom == iProc)CYCLE
       do iBuffer = nRecvCumSum + 1, nRecvCumSum + Router%nRecv_P(iProcFrom)
          call put_point_to_semirouter
       end do    ! iBuffer
       ! increment the offset
       nRecvCumSum = nRecvCumSum + Router%nRecv_P(iProcFrom)
    end do       ! iProcFrom
  contains
    !==========================================================================
    subroutine put_point_to_semirouter
      ! interpolation-related variables
      integer :: iIndexPut_II(0:GridTarget%nDim+1,2**GridTarget%nDim)
      real    :: CoordTarget_D(GridTarget%nDim)
      real    :: Weight_I(2**GridTarget%nDim)
      real    :: CoordPass_D(GridTarget%nDim)
      integer :: iImage, nImage, nImagePart, iProcTo, iToPut
      !------------------------------------------------------------------------
      CoordTarget_D = Router%BufferTarget_II(1:Router%nDim, iBuffer)
      CoordPass_D   = CoordTarget_D
      if(DoInterpolate)then
         call interpolate(                   &
              nDim           = nDimTarget,   &
              Coord_D        = CoordPass_D,  &
              Grid           = GridTarget,   &
              nIndex         = nIndexTarget, &
              iIndex_II      = iIndexPut_II, &
              nImage         = nImage,       &
              Weight_I       = Weight_I)
      else
         call nearest_grid_points(nDimTarget, CoordPass_D, GridTarget, &
              nIndexTarget, iIndexPut_II, nImage, Weight_I)
      end if
      nImagePart = count(iIndexPut_II(0,1:nImage)==iProc)
      Router%nPut_P(iProcFrom) = Router%nPut_P(iProcFrom) + nImagePart
      if(nImagePart==0)&
           call CON_stop(NameSub//'No image on the receiving PE')
      ! indices
      do iImage = 1, nImage
         iProcTo = iIndexPut_II(0,iImage)
         if(iProc==iProcTo)then
            iToPut = Router%nPut_P(iProcFrom) + 1 - nImagePart
            Router%iPut_P(iProcFrom)%iCB_II(:,iToPut)&
                 = iIndexPut_II(:,iImage)
            Router%iPut_P(iProcFrom)%iCB_II(0,iToPut) = nImagePart
            Router%Put_P(iProcFrom)%Weight_I(iToPut)  = Weight_I(iImage)
            Router%DoAdd_P(iProcFrom)%DoAdd_I(iToPut) = .true.
            nImagePart = nImagePart - 1
         end if
      end do ! iImage
    end subroutine put_point_to_semirouter
    !==========================================================================
  end subroutine update_semi_router_at_target
  !============================================================================
  subroutine check_size(nRank, nSize_I, &
       Buffer_I, Buffer_II, PBuffer_I, PBuffer_II,&
       iBuffer_I, iBuffer_II, DoBuffer_I)
    ! check if size of a given array is sufficient;
    ! can call the subroutine for only one array at a time
    integer,                     intent(in)   :: nRank
    integer,                     intent(in)   :: nSize_I(nRank)
    real,   allocatable,optional,intent(inout)::  Buffer_I(:), Buffer_II(:,:)
    real,   pointer,    optional,intent(inout):: PBuffer_I(:),PBuffer_II(:,:)
    integer,allocatable,optional,intent(inout):: iBuffer_I(:),iBuffer_II(:,:)
    logical,allocatable,optional,intent(inout)::DoBuffer_I(:)
    integer:: iError
    logical:: IsPresent_I(7)
    character(len=*), parameter:: NameSub = 'check_size'
    !--------------------------------------------------------------------------
    IsPresent_I = [&
         present(  Buffer_I), present( Buffer_II), &
         present( PBuffer_I), present(PBuffer_II), &
         present( iBuffer_I), present(iBuffer_II), &
         present(DoBuffer_I)]

    if(count(IsPresent_I) /= 1)&
         call CON_stop(NameSub // ': incorrect call')

    select case(nRank)
    case(1)
       if(present(Buffer_I))then
          if(allocated(Buffer_I))then
             if(any(ubound(Buffer_I) < nSize_I))then
                deallocate(Buffer_I)
             else
                RETURN
             end if
          end if
          allocate(Buffer_I(nSize_I(1)),stat=iError)
          call check_allocate(iError,NameSub)
          RETURN
       elseif(present(PBuffer_I))then
          if(associated(PBuffer_I))then
             if(any(ubound(PBuffer_I) < nSize_I))then
                deallocate(PBuffer_I)
             else
                RETURN
             end if
          end if
          allocate(PBuffer_I(nSize_I(1)),stat=iError)
          call check_allocate(iError,NameSub)
          RETURN
       elseif(present(iBuffer_I))then
          if(allocated(iBuffer_I))then
             if(any(ubound(iBuffer_I) < nSize_I))then
                deallocate(iBuffer_I)
             else
                RETURN
             end if
          end if
          allocate(iBuffer_I(nSize_I(1)),stat=iError)
          call check_allocate(iError,NameSub)
          RETURN
       elseif(present(DoBuffer_I))then
          if(allocated(DoBuffer_I))then
             if(any(ubound(DoBuffer_I) < nSize_I))then
                deallocate(DoBuffer_I)
             else
                RETURN
             end if
          end if
          allocate(DoBuffer_I(nSize_I(1)),stat=iError)
          call check_allocate(iError,NameSub)
          RETURN
       end if
    case(2)
       if(present(Buffer_II))then
          if(allocated(Buffer_II))then
             if(ubound(Buffer_II,1)/=nSize_I(1).or.&
                  ubound(Buffer_II,2)<nSize_I(2))then
                deallocate(Buffer_II)
             else
                RETURN
             end if
          end if
          allocate(Buffer_II(nSize_I(1), nSize_I(2)),stat=iError)
          call check_allocate(iError,NameSub)
          RETURN
       elseif(present(iBuffer_II))then
          if(allocated(iBuffer_II))then
             if(ubound(iBuffer_II,1)/=nSize_I(1).or.&
                  ubound(iBuffer_II,2)<nSize_I(2))then
                deallocate(iBuffer_II)
             else
                RETURN
             end if
          end if
          allocate(iBuffer_II(nSize_I(1), nSize_I(2)),stat=iError)
          call check_allocate(iError,NameSub)
          RETURN
       elseif(present(PBuffer_II))then
          if(associated(PBuffer_II))then
             if(ubound(PBuffer_II,1)/=nSize_I(1).or.&
                  ubound(PBuffer_II,2)<nSize_I(2))then
                deallocate(PBuffer_II)
             else
                RETURN
             end if
          end if
          allocate(PBuffer_II(nSize_I(1), nSize_I(2)),stat=iError)
          call check_allocate(iError,NameSub)
          RETURN
       end if
    case default
       call CON_stop(NameSub // ': incorrect call')
    end select
  end subroutine check_size
  !============================================================================

end module CON_router
!==============================================================================

