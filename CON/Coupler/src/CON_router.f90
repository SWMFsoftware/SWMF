!^CFG COPYRIGHT UM                                              !
!BOP
!MODULE: CON_router - set the connection between the grids of different models
!INTERFACE:
Module CON_router
  !USES:
  use CON_grid_descriptor
  use ModUtilities, ONLY: check_allocate
  !DESCRIPTION:
!This file presents the class of routers between the grids, each!  
!of them can be either the uniformly spaced or Octree or Quadric!
!adaptive block grid                   .                        !
!
!The methods include: allocation, initialization, cleaner and   !
!two different constructors.                                    !
!
!REVISION HISTORY:
! Sokolov I.V.                                                  !
! 7.20.03-7.21.03                                               !
! igorsok@umich.edu                                             !
! phone(734)647-4705                                            !
!EOP
  implicit none

  logical,parameter:: UseUnionComm=.true.
!BOP
!DESCRIPTION:
!==========================DERIVED TYPES========================!
!\begin{verbatim}
  type DoAddPtrType
     logical,dimension(:),pointer::DoAdd_I
  end type DoAddPtrType
!\end{verbatim}
!---------------------------------------------------------------!
!See CON\_grid\_descriptor about iCB index. In the array iCB\_I !
!the second index enumerates the grid points belonging to some  !
!list, while the first one numerates the position of (0) PE, at !
!which the point is localized, (1:nDim) grid point indexes in   !
!the block and (nDim+1), if exists, stores the local block      !
!number.
!\begin{verbatim}
  type IndexPtrType
      integer,dimension(:,:),pointer::iCB_II
   end type IndexPtrType

  type WeightPtrType
     real,dimension(:),pointer::Weight_I
  end type WeightPtrType
!\end{verbatim}
!========================DERIVED TYPE===========================!
  type RouterType
!The router can be set between the grids of different dimensions!
     character(LEN=3)::Name

!For the router between LOCAL grids of a component we use the   !
!communicator of the model for sending-receiving the data,      !
!otherwise the global communicator                              !
     logical::IsLocal,IsProc
     integer::iProc,nProc,iComm

!If the union group is constructed, then for use with broadcast !
!we need the union communicator and the root PE ranks in this   !
!communicator
     integer::iCommUnion,iProc0Source,iProc0Target
!\begin{verbatim}
!As the default we use iCB indexes to construct the router,     !
!hence the grid point is characterized by the                   !
!GridDescriptor%nDim grid point indexes plus one more index for !
!the block number. Also we allow to use exactly                 ! 
!GridDescriptor%nDim indexes, without the block number which    !
!only seems to be on sence for the component which is localized !
!at one PE only, or which has exactly one block per PE          !
 
     integer::nIndexesSource,nIndexesTarget  

!The total amounts of the buffer segmentsto be sent-received    !
!to/from the PE. The total amounts of the grid points from which!
!the data should be got or to which the data should be put,some !
!data points may be counted more than one time                  !

     integer, dimension(:), pointer :: &
          nGet_P, nPut_P, nRecv_P, nSend_P
!iCB indexes and the weight coefficients for the points of the  !
!target and source grids, which are connected through the router!
     type(IndexPtrType), dimension(:), pointer :: iGet_P
     type(IndexPtrType), dimension(:), pointer :: iPut_P
     type(DoAddPtrType), dimension(:), pointer :: DoAdd_P
     type(WeightPtrType), dimension(:),pointer :: Get_P
     type(WeightPtrType), dimension(:),pointer :: Put_P
!\end{verbatim}
  end type RouterType
!EOP
!BOP
!PUBLIC MEMBER FUNCTIONS:
  private::allocate_get_arrays
  private::allocate_put_arrays
  private::check_router_allocation
!EOP
contains
!BOP
!IROUTINE: clean_router - destructor for the type
!INTERFACE:
  subroutine clean_router(Router)
    !INPUT ARGUMENTS:
    type(RouterType),intent(inout)::Router
!EOP
    integer::iPE,iError
    do iPE=0,Router%nProc-1
       if(associated(Router%iGet_P(iPE)%iCB_II))&
            deallocate(Router%iGet_P(iPE)%iCB_II)
       if(associated(Router%Get_P(iPE)%Weight_I))&
            deallocate(Router%Get_P(iPE)%Weight_I)
       if(associated(Router%iPut_P(iPE)%iCB_II))&
            deallocate(Router%iPut_P(iPE)%iCB_II)
       if(associated(Router%Put_P(iPE)%Weight_I))&
            deallocate(Router%Put_P(iPE)%Weight_I)
       if(associated(Router%DoAdd_P(iPE)%DoAdd_I))&
            deallocate(Router%DoAdd_P(iPE)%DoAdd_I)
    end do
    if(associated(Router%nGet_P))&
         deallocate(Router%nGet_P)
    if(associated(Router%nPut_P))&
         deallocate(Router%nPut_P)
    if(associated(Router%nSend_P))&
         deallocate(Router%nSend_P)
    if(associated(Router%nRecv_P))&
         deallocate(Router%nRecv_P)
    if(associated(Router%iGet_P))&
         deallocate(Router%iGet_P)
    if(associated(Router%Get_P))&
         deallocate(Router%Get_P)
    if(associated(Router%iPut_P))&
         deallocate(Router%iPut_P)
    if(associated(Router%Put_P))&
         deallocate(Router%Put_P)
    if(associated(Router%DoAdd_P))&
         deallocate(Router%DoAdd_P)
    if(UseUnionComm.and..not.Router%IsLocal.and.Router%IsProc)&
         call MPI_COMM_FREE(Router%iCommUnion,iError)
  end subroutine clean_router
!===============================================================!
!BOP
!IROUTINE: init_router - initialize the type
!INTERFACE:
  subroutine init_router(&
       GridDescriptorSource,&
       GridDescriptorTarget,&
       Router,&
       nIndexesSource,&
       nIndexesTarget)
    !INPUT ARGUMENTS:
    type(GridDescriptorType),intent(in)::&
         GridDescriptorSource,&
         GridDescriptorTarget
    type (RouterType),intent(out)::Router
    integer,intent(in),optional::nIndexesSource
    integer,intent(in),optional::nIndexesTarget
    !EOP
    integer::iPE,iError,LocalCompID_
    integer::nProc
    integer::iProc0Source,iProc0Target,iProcUnion
    integer::iGroupUnion,iGroupSource,iGroupTarget
    integer,dimension(1)::iRankGlobal_I,iRankUnion_I

!---------------------------------------------------------------!
!Check if the grids are both local or both global               !

    if(is_local_grid(GridDescriptorSource%DD%Ptr).and.&
         is_local_grid(GridDescriptorTarget%DD%Ptr))then
       Router%IsLocal=.true.
       LocalCompID_=compid_grid(GridDescriptorSource%DD%Ptr)
       if( LocalCompID_/=&
            compid_grid(GridDescriptorTarget%DD%Ptr))&
            call CON_stop(&
            'Do not couple Local grids of different components!')
       Router%iProc=i_proc(LocalCompID_)
       Router%nProc=n_proc(LocalCompID_)
       Router%iComm=i_comm(LocalCompID_)
       Router%iCommUnion=Router%iComm
       Router%iProc0Source=0
       Router%iProc0Target=0
       Router%IsProc=is_proc(LocalCompID_)
    elseif((.not.is_local_grid(GridDescriptorSource%DD%Ptr))&
         .and.(.not.is_local_grid(GridDescriptorTarget%DD%Ptr)))&
         then
       Router%IsLocal=.false.
       Router%iProc=i_proc()
       Router%nProc=n_proc()
       Router%iComm=i_comm()
       Router%IsProc=is_proc()
       iProc0Source=i_proc0(&
               compid_grid(GridDescriptorSource%DD%Ptr))
       iProc0Target=i_proc0(&
               compid_grid(GridDescriptorTarget%DD%Ptr))
       if(UseUnionComm)then
          iGroupSource=i_group(&
               compid_grid(GridDescriptorSource%DD%Ptr))
          iGroupTarget=i_group(&
               compid_grid(GridDescriptorTarget%DD%Ptr))
          if(iProc0Target>iProc0Source)then
             call MPI_GROUP_UNION(&
                  iGroupSource,&
                  iGroupTarget,&
                  iGroupUnion,&
                  iError)
             Router%iProc0Source=0
             iRankGlobal_I=iProc0Target
             call MPI_GROUP_TRANSLATE_RANKS(&
                  i_group(),&
                  1,&
                  iRankGlobal_I(1),&
                  iGroupUnion,&
                  iRankUnion_I(1),&
                  iError)
             Router%iProc0Target=iRankUnion_I(1)
          else
             call MPI_GROUP_UNION(&
                  iGroupTarget,&
                  iGroupSource,&
                  iGroupUnion,&
                  iError)
             Router%iProc0Target=0
             iRankGlobal_I=iProc0Source
             call MPI_GROUP_TRANSLATE_RANKS(&
                  i_group(),&
                  1,&
                  iRankGlobal_I(1),&
                  iGroupUnion,&
                  iRankUnion_I(1),&
                  iError)
             Router%iProc0Source=iRankUnion_I(1)
          end if
          call MPI_COMM_CREATE(&
               i_comm(),&
               iGroupUnion,&
               Router%iCommUnion,&
               iError)
          call MPI_group_rank(iGroupUnion,iProcUnion,iError)
          Router%IsProc=iProcUnion/=MPI_UNDEFINED
          if(iProcUnion/=Router%iProc0Target.and.&
               i_proc()==iProc0Target)call CON_stop(&
               'Wrongly defined Router%iProc0Target=',&
               Router%iProc0Target)
          if(iProcUnion/=Router%iProc0Source.and.&
               i_proc()==iProc0Source)call CON_stop(&
               'Wrongly defined Router%iProc0Source=',&
               Router%iProc0Source)
          call MPI_GROUP_FREE(iGroupUnion,iError)
       else
          Router%iCommUnion=Router%iComm
          Router%iProc0Source=iProc0Source
          Router%iProc0Target=iProc0Target
       end if
       write(*,*)
    else
       call CON_stop(&
            'Do not couple a Local grid with a global one')
    end if
    Router%nIndexesSource=ndim_grid(&
       GridDescriptorSource%DD%Ptr)+1
    if(present(nIndexesSource))then
       if(nIndexesSource>=Router%nIndexesSource-1)then
          Router%nIndexesSource=nIndexesSource
       else
          call CON_stop(&
               'nIndexesSource should be at least ',&
               Router%nIndexesSource-1)
       end if
    end if
    Router%nIndexesTarget=ndim_grid(&
       GridDescriptorTarget%DD%Ptr)+1

    if(present(nIndexesTarget))then
       if(nIndexesTarget>=Router%nIndexesTarget-1)then
          Router%nIndexesTarget=nIndexesTarget
       else
          call CON_stop(&
               'nIndexesTarget should be at least ',&
               Router%nIndexesTarget-1)
       end if
    end if
    nProc=Router%nProc
    allocate(Router%nGet_P(0:nProc-1),stat=iError)
    call check_allocate(iError,'nGet_P')
    allocate(Router%nPut_P(0:nProc-1),stat=iError)
    call check_allocate(iError,'nPut_P')
    allocate(Router%nSend_P(0:nProc-1),stat=iError)
    call check_allocate(iError,'nSend_P')
    allocate(Router%nRecv_P(0:nProc-1),stat=iError)
    call check_allocate(iError,'nRecv_P')
    allocate(Router%iGet_P(0:nProc-1),stat=iError)
    call check_allocate(iError,'iGet_P') 
    allocate(Router%Get_P(0:nProc-1),stat=iError)
    call check_allocate(iError,'Get_P') 
    do iPE=0,nProc-1
       nullify(Router%iGet_P(iPE)%iCB_II)
       nullify(Router%Get_P(iPE)%Weight_I)
       call allocate_get_arrays(Router,iPE,1)
    end do
    allocate(Router%iPut_P(0:nProc-1),stat=iError)
    call check_allocate(iError,'iPut_P') 
    allocate(Router%Put_P(0:nProc-1),stat=iError)
    call check_allocate(iError,'Put_P')
    allocate(Router%DoAdd_P(0:nProc-1),stat=iError)
    call check_allocate(iError,'DoAdd_P')

    do iPE=0,nProc-1
       nullify(Router%iPut_P(iPE)%iCB_II)
       nullify(Router%Put_P(iPE)%Weight_I)
       nullify(Router%DoAdd_P(iPE)%DoAdd_I)
       call allocate_put_arrays(Router,iPE,1)
    end do

    Router%nGet_P=0
    Router%nPut_P=0
    Router%nSend_P=0
    Router%nRecv_P=0
    do iPE=0,nProc-1
       Router%iGet_P(iPE)%iCB_II(:,:)=0
       Router%iPut_P(iPE)%iCB_II(:,:)=0
       Router%Put_P(iPE)%Weight_I(:)=cZero
       Router%Get_P(iPE)%Weight_I(:)=cZero
       Router%DoAdd_P(iPE)%DoAdd_I(:)=.false.
    end do
  end subroutine init_router
!============================PRIVATE============================!
  subroutine  allocate_get_arrays(Router,iPE,nLength)
    type(RouterType),intent(inout)::Router
    integer,intent(in)::iPE,nLength
    integer::iError
    if(associated(Router%iGet_P(iPE)%iCB_II))&
         deallocate(Router%iGet_P(iPE)%iCB_II)
    if(associated(Router%Get_P(iPE)%Weight_I))&
         deallocate(Router%Get_P(iPE)%Weight_I)
    allocate(Router%iGet_P(iPE)%iCB_II(&
         0:Router%nIndexesSource,nLength),stat=iError)
    call check_allocate(iError,'iGet_P%iCB_II')
    allocate(Router%Get_P(iPE)%Weight_I(nLength),stat=iError)
    call check_allocate(iError,'Get_P%Weight_I')
  end subroutine allocate_get_arrays
!============================PRIVATE============================!
  subroutine allocate_put_arrays(Router,iPE,nLength)
    type(RouterType),intent(inout)::Router
    integer,intent(in)::iPE,nLength
    integer::iError
    if(associated(Router%iPut_P(iPE)%iCB_II))&
         deallocate(Router%iPut_P(iPE)%iCB_II)
    if(associated(Router%DoAdd_P(iPE)%DoAdd_I))&
         deallocate(Router%DoAdd_P(iPE)%DoAdd_I)
    allocate(Router%iPut_P(iPE)%iCB_II(&
         0:Router%nIndexesTarget,nLength),stat=iError)
    call check_allocate(iError,'iPut_P%iCB_II') 
    allocate(Router%Put_P(iPE)%Weight_I(nLength),stat=iError)
    call check_allocate(iError,'Put_P%Weight_I')   
    allocate(Router%DoAdd_P(iPE)%DoAdd_I(nLength),stat=iError)
    call check_allocate(iError,'DoAdd_P%DoAdd_I')
  end subroutine allocate_put_arrays
!============================PRIVATE============================!
!
  subroutine check_router_allocation(Router,iProc,nProc)
    integer,intent(in)::iProc,nProc
    type(RouterType),intent(inout)::Router
    integer::iError,iPE,lLength
    do iPE=0,nProc-1
       if(ubound(&
            Router%iPut_P(iPE)%iCB_II,2)<&
            Router%nPut_P(iPE))then
          call allocate_put_arrays&
               (Router,iPE,Router%nPut_P(iPE))
       end if
       if(ubound(&
            Router%iGet_P(iPE)%iCB_II,2)<&
            Router%nGet_P(iPE))then
          call allocate_get_arrays&
               (Router,iPE,Router%nGet_P(iPE))
       end if
    end do
  end subroutine check_router_allocation
!===============================================================!
!Done up to this place 7.21.03                                  !

!===============================================================!
!===============================================================!
!BOP
!IROUTINE: set_router - work for a case of mapping FROM TARGET
!INTERFACE:
  subroutine set_router(&
!The Descriptor for Source grid points                          !
       GridDescriptorSource,& 
!The Descriptor for Target grid points                          !
       GridDescriptorTarget,&
!The Router to be set                                           !
       Router,&
!Logical function which allows to skip the block if there is no !
!interface points in it. Optional, if not present then all the  !
!blocks are checked for the presence of the interface points    !
       is_interface_block,&
!The subroutine which defines if the grid point is inside the   !
!interface layer. Optional, if not present, then all the grid   !
!points (at the target grid) are considered as the interface    !
!layer points                                                   !     
       interface_point_coords, &
!Mapping transformation which, in the treated case, maps the    !
!target grid point to an image point into the source domain     !    
       mapping,&
!First or second order interpolation procedures are available   !
!(see CON_grid_descriptor) to find the grid point at the source !
!grid and the interpoltion weights for the image point:         !
!nearest_grid_points and bilinear_interpolation.                !
!Optional, if not present then the nearest_grid_points is used  !
       interpolate) 
    !INPUT ARGUMENTS:
    interface
       logical function is_interface_block(lGlobalNode)
         implicit none
         integer,intent(in)::lGlobalNode 
       end function is_interface_block
       subroutine interface_point_coords(&
                  GridDescriptor,&
                  lGlobalTreeNode,&
                  Xyz_D,&
                  nIndexes,&
                  Index_I,&
                  IsInterfacePoint)
         use CON_grid_descriptor
         implicit none
         type(GridDescriptorType),intent(in)::GridDescriptor
         integer,intent(in)::lGlobalTreeNode,nIndexes
         logical,intent(out)::IsInterfacePoint
         real,dimension(GridDescriptor%nDim),intent(inout)::Xyz_D
         integer, dimension(nIndexes),intent(inout)::Index_I
       end subroutine interface_point_coords
       subroutine mapping(&
            nDimFrom,XyzFrom_D,nDimTo,XyzTo_D,IsInterfacePoint)
         implicit none                                                  
         integer,intent(in)::nDimFrom,nDimTo       
         real,dimension(nDimFrom),intent(in)::XyzFrom_D
         real,dimension(nDimTo),intent(out)::XyzTo_D
         logical,intent(out)::IsInterfacePoint
       end subroutine mapping
       subroutine interpolate(Xyz_D,&
                              GridDescriptor,&
                              nIndexes,&
                              Index_II,&
                              nImages,Weight_I)
         use CON_grid_descriptor
         implicit none
         type(GridDescriptorType)::GridDescriptor     
         real,dimension(GridDescriptor%nDim),&
              intent(inout)::Xyz_D 
         integer,intent(in)::nIndexes
         integer,dimension(&
              0:nIndexes,2**GridDescriptor%nDim)::Index_II
         integer,intent(out)::nImages
         real,dimension(2**GridDescriptor%nDim),&
              intent(out)::Weight_I
       end subroutine interpolate
    end interface

    optional::is_interface_block,interface_point_coords
    optional::mapping,interpolate

    type(GridDescriptorType),intent(in):: GridDescriptorSource
    type(GridDescriptorType),intent(in):: GridDescriptorTarget
    type(RouterType),intent(inout)::Router
    !EOP
    integer::iProc,nProc
    integer::lGlobalNode,lFound
    integer::iGlobalGridPoint,nGridPointsPerBlock
    logical::IsInterfacePoint
    integer::iImages,nImages,nImagesPart,iToGet
    integer::iProcTo,iProcFrom,iProcDoNotAdd,iPE
    integer,dimension(0:Router%nProc-1)::&
         nGetUbound_P,nPutUbound_P

    real,dimension(GridDescriptorTarget%nDim)::XyzTarget_D
    real,dimension(GridDescriptorSource%nDim)::XyzSource_D

    integer, dimension(Router%nIndexesTarget)::IndexRecv_I
    integer,dimension(0:Router%nIndexesSource,&
                      2**GridDescriptorSource%nDim)::&
                      IndexGet_II
    integer,dimension(2**GridDescriptorSource%nDim)::&
                      iProcLookUp_I
    integer::nProcToGet,iProcToGet
    logical::DoCountOnly,DoCountRed
    real,dimension(2**GridDescriptorSource%nDim)::Weight_I 
    integer::iError
    iProc=Router%iProc
    if(iProc<0)return
    if(.not.present(mapping).and.GridDescriptorTarget%nDim/=&
                                 GridDescriptorSource%nDim)&
            call CON_stop(&
            'Mapping is needed for Target%nDim/=Source%nDim')
    nProc=Router%nProc

    !Check dimensions

    DoCountOnly=.true. !To enter the loop
    do while(DoCountOnly)
       call check_router_allocation(Router,iProc,nProc)


!Store Upper bounds to control if the alllocated index array    !
!have sufficient size
       do iPE=0,nProc-1
          nGetUbound_P(iPE)=ubound(Router%iGet_P(iPE)%iCB_II,2)
          nPutUbound_P(iPE)=ubound(Router%iPut_P(iPE)%iCB_II,2)
       end do

       DoCountOnly=.false.
!If the check shows that the allocated array is not sufficient, ! 
!then DoCountOnly will be set to true. The loop then will be    !
!repeated for the second time                                   !



!Initialize the counters                                        !
       Router%nGet_P=0
       Router%nPut_P=0
       Router%nSend_P=0
       Router%nRecv_P=0

       nGridPointsPerBlock=n_grid_points_per_block(&
            GridDescriptorTarget)

!Block loop                                                     !
       do lGlobalNode=1,ntree_nodes(GridDescriptorTarget%DD%Ptr)

!Skip non-end octree nodes, if any                              !
          if(.not.is_used_block(GridDescriptorTarget%DD%Ptr,&
               lGlobalNode))CYCLE

!Skip the block if desired: if there is known to be no interface!
!point in it                                                    !
          if(present(is_interface_block))then
             if(.not.is_interface_block(lGlobalNode))CYCLE
          end if
!GlobalCellNumber Loop, for a given (octree) block              !
          do iGlobalGridPoint=&
                         1+nGridPointsPerBlock*(lGlobalNode-1),&
                         nGridPointsPerBlock*lGlobalNode

             call pe_and_blk(&
                  GridDescriptorTarget%DD%Ptr,lGlobalNode,&
                  iProcTo,IndexRecv_I(Router%nIndexesTarget))

             call global_i_grid_point_to_icb(&
                  GridDescriptorTarget,&
                  iGlobalGridPoint,&
                  lFound,& !Should be the same as lGlobalNode
                  IndexRecv_I(1:GridDescriptorTarget%nDim))
             XyzTarget_D=xyz_grid_d(&
                  GridDescriptorTarget,&
                  lGlobalNode,&
                  IndexRecv_I(1:GridDescriptorTarget%nDim))                    
             
             if(present(interface_point_coords))then
                call interface_point_coords(&
                     GridDescriptorTarget,&
                     lGlobalNode,&
                     XyzTarget_D,&
                     Router%nIndexesTarget,&
                     IndexRecv_I,&
                     IsInterfacePoint)
                if(.not.IsInterfacePoint)CYCLE 
             end if

             if(present(mapping))then
                call mapping(&
                     GridDescriptorTarget%nDim,&
                     XyzTarget_D,&
                     GridDescriptorSource%nDim,&
                     XyzSource_D,&
                     IsInterfacePoint)
                if(.not.IsInterfacePoint)CYCLE
             else
                XyzSource_D=XyzTarget_D
             end if

             if(present(interpolate))then
                call interpolate(&
                     XyzSource_D,&
                     GridDescriptorSource,&
                     Router%nIndexesSource,&
                     IndexGet_II,&
                     nImages,&
                     Weight_I)
             else
                call nearest_grid_points(&
                     XyzSource_D,&
                     GridDescriptorSource,&
                     Router%nIndexesSource,&
                     IndexGet_II,&
                     nImages,&
                     Weight_I)
             end if
!--------------------------------------------------------------!
!Lookup
             nImagesPart=0     !At all CPUs

                         
             do iImages=1,nImages
                iProcFrom=IndexGet_II(0,iImages)
!At the source PEs the number of terms in the partial sums are !
!found                                                         !
                if(iProc==iProcFrom)then
                   nImagesPart=nImagesPart+1
                   Router%nGet_P(iProcTo)=&
                        Router%nGet_P(iProcTo)+1
                   DoCountOnly=DoCountOnly.or.&
                        Router%nGet_P(iProcTo)>&
                        nGetUbound_P(iProcTo)
                end if
                
!At the target processor the PE list is defined which will send!
!partial sums                                                 !
                if(iProc==iProcTo)then
                   if(iImages==1)then
                      iProcLookUp_I(1)=iProcFrom
                      nProcToGet=1
                      Router%nPut_P(iProcFrom)=&
                           Router%nPut_P(iProcFrom)+1
                      Router%nRecv_P(iProcFrom)=&
                           Router%nRecv_P(iProcFrom)+1
                      DoCountOnly=DoCountOnly.or.&
                           Router%nPut_P(iProcFrom)>&
                           nPutUbound_P(iProcFrom)
                   else
                      if(.not.any(iProcLookUp_I(&
                           1:nProcToGet)==iProcFrom))then
                         nProcToGet=nProcToGet+1
                         iProcLookUp_I(nProcToGet)=iProcFrom
                         Router%nPut_P(iProcFrom)=&
                              Router%nPut_P(iProcFrom)+1
                         Router%nRecv_P(iProcFrom)=&
                              Router%nRecv_P(iProcFrom)+1
                         DoCountOnly=DoCountOnly.or.&
                              Router%nRecv_P(iProcFrom)>&
                              nPutUbound_P(iProcFrom)
                      end if
                   end if
                end if
             end do
             
             if(nImagesPart>0)Router%nSend_P(iProcTo)=&
                  Router%nSend_P(iProcTo)+1
             
             if(.not.DoCountOnly)then
                do iImages=1,nImages
                   iProcFrom=IndexGet_II(0,iImages)
                   if(iProc==iProcFrom)then
                      iToGet=Router%nGet_P(iProcTo)+1-nImagesPart
                      Router%iGet_P(iProcTo)%iCB_II(:,iToGet)&
                           =IndexGet_II(:,iImages)
                      Router%iGet_P(iProcTo)%iCB_II(0,iToGet)&
                           =nImagesPart
                      Router%Get_P(iProcTo)%Weight_I(iToGet)&
                           =Weight_I(iImages)
                      nImagesPart=nImagesPart-1
                   end if
                end do
                if(iProc==iProcTo)then
                   do iProcToGet=1,nProcToGet
                      iProcFrom=iProcLookUp_I(iProcToGet)
                      Router%iPut_P(iProcFrom)%&
                           iCB_II(1:Router%nIndexesTarget,&
                           Router%nPut_P(iProcFrom))&
                           =IndexRecv_I(1:Router%nIndexesTarget)
                      Router%iPut_P(iProcFrom)%&
                           iCB_II(0,Router%nPut_P(iProcFrom))&
                           =1
                      Router%Put_P(iProcFrom)%&
                           Weight_I(Router%nPut_P(iProcFrom))&
                           =cOne
                      Router%DoAdd_P(iProcFrom)%&
                           DoAdd_I(Router%nRecv_P(iProcFrom))=&
                           .true.
                   end do
                   if(any(iProcLookUp_I(&
                        1:nProcToGet)==iProcTo))then
                      iProcDoNotAdd=iProcTo
                   else
                      iProcDoNotAdd=minval(&
                           iProcLookUp_I(1:nProcToGet))
                   end if
                   Router%DoAdd_P(iProcDoNotAdd)%&
                        DoAdd_I(Router%nRecv_P(iProcDoNotAdd))=&
                        .false.
                end if
             end if
          end do !Global cell
       end do    !Target block


    end do       !Check if DoCountOnly

  end subroutine set_router
!===============================================================!
!===============================================================!
!This subroutine works for a case of mapping FROM SOURCE
!BOP
!IROUTINE: construct_router_from_source - work for a case of mapping FROM SOURCE
!INTERFACE: 
  subroutine construct_router_from_source(&
!The Descriptor for Source grid points   
       GridDescriptorSource,& 
!The Descriptor for Target grid points    
       GridDescriptorTarget,& 
!The Router to be set 
       Router,&
!Logical function which allows to skip the block if there is no !
!interface points in it. Optional, if not present then all the  !
!blocks are checked for the presence of the interface points    !
       is_interface_block,&
!The subroutine which defines if the grid point is inside the   !
!interface layer. Optional, if not present, then all the grid   !
!points (at the source grid) are considered as the interface    !
!layer points 
       interface_point_coords, &
!Mapping transformation which, in the treated case, maps the    !
!source grid point to an image point into the target domain     !    
       mapping,&
!First or second order interpolation procedures are available   !
!(see CON_grid_descriptor) to find the grid point at the target !
!grid and the interpoltion weights for the image point:         !
!nearest_grid_points and bilinear_interpolation.                !
!Optional, if not present then the nearest_grid_points is used  !
       interpolate) !First or second order interpolation
    !INPUT ARGUMENTS:
    interface
       logical function is_interface_block(lGlobalNode)
         implicit none
         integer,intent(in)::lGlobalNode 
       end function is_interface_block
       subroutine interface_point_coords(&
                  GridDescriptor,&
                  lGlobalTreeNode,&
                  Xyz_D,&
                  nIndexes,&
                  Index_I,&
                  IsInterfacePoint)
         use CON_grid_descriptor
         implicit none
         type(GridDescriptorType),intent(in)::GridDescriptor
         integer,intent(in)::lGlobalTreeNode,nIndexes
         logical,intent(out)::IsInterfacePoint
         real,dimension(GridDescriptor%nDim),intent(inout)::Xyz_D
         integer, dimension(nIndexes),intent(inout)::Index_I
       end subroutine interface_point_coords
       subroutine mapping(&
            nDimFrom,XyzFrom_D,nDimTo,XyzTo_D,IsInterfacePoint)
         implicit none                                                  
         integer,intent(in)::nDimFrom,nDimTo       
         real,dimension(nDimFrom),intent(in)::XyzFrom_D
         real,dimension(nDimTo),intent(out)::XyzTo_D
         logical,intent(out)::IsInterfacePoint
       end subroutine mapping
       subroutine interpolate(Xyz_D,&
                              GridDescriptor,&
                              nIndexes,&
                              Index_II,&
                              nImages,Weight_I)
         use CON_grid_descriptor
         implicit none
         type(GridDescriptorType)::GridDescriptor     
         real,dimension(GridDescriptor%nDim),&
              intent(inout)::Xyz_D 
         integer,intent(in)::nIndexes
         integer,dimension(&
              0:nIndexes,2**GridDescriptor%nDim)::Index_II
         integer,intent(out)::nImages
         real,dimension(2**GridDescriptor%nDim),&
              intent(out)::Weight_I
       end subroutine interpolate
    end interface

    optional::is_interface_block,interface_point_coords
    optional::mapping,interpolate

    type(GridDescriptorType),intent(in):: GridDescriptorSource
    type(GridDescriptorType),intent(in):: GridDescriptorTarget
    type(RouterType),intent(inout)::Router
    !EOP
    integer::iProc,nProc
    integer::lGlobalNode,lFound
    integer::iGlobalGridPoint,nGridPointsPerBlock
    logical::IsInterfacePoint
    integer::iImages,nImages,nImagesPart,iToPut
    integer::iProcTo,iProcFrom,iPE
    integer,dimension(0:Router%nProc-1)::&
         nGetUbound_P,nPutUbound_P

    real,dimension(GridDescriptorTarget%nDim)::XyzTarget_D
    real,dimension(GridDescriptorSource%nDim)::XyzSource_D

    integer, dimension(Router%nIndexesSource)::IndexGet_I
    integer,dimension(0:Router%nIndexesTarget,&
                      2**GridDescriptorTarget%nDim)::&
                      IndexPut_II
    integer,dimension(2**GridDescriptorTarget%nDim)::&
                      iProcLookUp_I
    integer::nProcToPut,iProcToPut
    logical::DoCountOnly,DoCountRed
    real,dimension(2**GridDescriptorTarget%nDim)::Weight_I 
    integer::iError

    iProc=Router%iProc
    if(iProc<0)return
    if(.not.present(mapping).and.GridDescriptorTarget%nDim/=&
                                 GridDescriptorSource%nDim)&
            call CON_stop(&
            'Mapping is needed for Target%nDim/=Source%nDim')
    nProc=Router%nProc

    !Check dimensions

    DoCountOnly=.true. !To enter the loop
    do while(DoCountOnly)
       call check_router_allocation(Router,iProc,nProc)


!Store Upper bounds to control if the alllocated index array    !
!have a sufficient length
       do iPE=0,nProc-1
          nGetUbound_P(iPE)=ubound(Router%iGet_P(iPE)%iCB_II,2)
          nPutUbound_P(iPE)=ubound(Router%iPut_P(iPE)%iCB_II,2)
       end do

       DoCountOnly=.false.
!If the check shows that the allocated array is not sufficient, ! 
!DoCountOnly will be set to true. The loop then will be repeated!
! for the second time



!Initialize the counters
       Router%nGet_P=0
       Router%nPut_P=0
       Router%nSend_P=0
       Router%nRecv_P=0

       nGridPointsPerBlock=n_grid_points_per_block(&
            GridDescriptorSource)

!Block loop
       do lGlobalNode=1,ntree_nodes(GridDescriptorSource%DD%Ptr)

!Skip non-end octree nodes, if any                              !
          if(.not.is_used_block(GridDescriptorSource%DD%Ptr,&
               lGlobalNode))CYCLE

!Skip the block if desired: if there is known to be no interface!
! point in it
          if(present(is_interface_block))then
             if(.not.is_interface_block(lGlobalNode))CYCLE
          end if

!GlobalCellNumber Loop, for a given (octree) block              !
          do iGlobalGridPoint=&
                         1+nGridPointsPerBlock*(lGlobalNode-1),&
                         nGridPointsPerBlock*lGlobalNode

             call pe_and_blk(&
                  GridDescriptorSource%DD%Ptr,lGlobalNode,&
                  iProcFrom,IndexGet_I(Router%nIndexesSource))

             call global_i_grid_point_to_icb(&
                  GridDescriptorSource,&
                  iGlobalGridPoint,&
                  lFound,& !Should be the same as lGlobalNode
                  IndexGet_I(1:GridDescriptorSource%nDim)) 
           
             XyzSource_D=xyz_grid_d(&
                  GridDescriptorSource,&
                  lGlobalNode,&
                  IndexGet_I(1:GridDescriptorSource%nDim))

             if(present(interface_point_coords))then
                call interface_point_coords(&
                     GridDescriptorSource,&
                     lGlobalNode,&
                     XyzSource_D,&
                     Router%nIndexesSource,&
                     IndexGet_I,&
                     IsInterfacePoint)
                if(.not.IsInterfacePoint)CYCLE 
             end if
 
             if(present(mapping))then
                call mapping(GridDescriptorSource%nDim,&
                             XyzSource_D,&
                             GridDescriptorTarget%nDim,&
                             XyzTarget_D,&
                             IsInterfacePoint)
             else
                XyzTarget_D=XyzSource_D
             end if


             if(present(interpolate))then
                call interpolate(&
                     XyzTarget_D,&
                     GridDescriptorTarget,&
                     Router%nIndexesTarget,&
                     IndexPut_II,&
                     nImages,&
                     Weight_I)
                else
                   call nearest_grid_points(&
                     XyzTarget_D,&
                     GridDescriptorTarget,&
                     Router%nIndexesTarget,&
                     IndexPut_II,&
                     nImages,&
                     Weight_I)
                end if

!See interface
!--------------------------------------------------------------!
!Lookup
             nImagesPart=0     !At all CPUs

             do iImages=1,nImages
                iProcTo=IndexPut_II(0,iImages)
!At the target PEs the number of terms in the partial sums are ! 
!found
                if(iProc==iProcTo)then
                   nImagesPart=nImagesPart+1
                   Router%nPut_P(iProcFrom)=&
                        Router%nPut_P(iProcFrom)+1
                   DoCountOnly=DoCountOnly.or.&
                        Router%nPut_P(iProcFrom)>&
                        nPutUbound_P(iProcFrom)
                end if

!At the source processor the PE list is defined which will get !
!the  partial sums
                if(iProc==iProcFrom)then
                   if(iImages==1)then
                      iProcLookUp_I(1)=iProcTo
                      nProcToPut=1
                      Router%nGet_P(iProcTo)=&
                           Router%nGet_P(iProcTo)+1
                      Router%nSend_P(iProcTo)=&
                           Router%nSend_P(iProcTo)+1
                      DoCountOnly=DoCountOnly.or.&
                           Router%nGet_P(iProcTo)>&
                           nGetUbound_P(iProcTo)
                   else
                      if(.not.any(iProcLookUp_I(&
                           1:nProcToPut)==iProcTo))then
                         nProcToPut=nProcToPut+1
                         iProcLookUp_I(nProcToPut)=iProcTo
                         Router%nGet_P(iProcTo)=&
                              Router%nGet_P(iProcTo)+1
                         Router%nSend_P(iProcTo)=&
                              Router%nSend_P(iProcTo)+1
                         DoCountOnly=DoCountOnly.or.&
                              Router%nGet_P(iProcTo)>&
                              nGetUbound_P(iProcTo)
                      end if
                   end if
                end if
             end do

             if(nImagesPart>0)Router%nRecv_P(iProcFrom)=&
                  Router%nRecv_P(iProcFrom)+1

             if(.not.DoCountOnly)then
                do iImages=1,nImages
                   iProcTo=IndexPut_II(0,iImages)
                   if(iProc==iProcTo)then
                      iToPut= Router%nPut_P(iProcFrom)+1&
                           -nImagesPart
                      Router%iPut_P(iProcFrom)%iCB_II(:,iToPut)&
                           =IndexPut_II(:,iImages)
                      Router%iPut_P(iProcFrom)%iCB_II(0,iToPut)&
                           =nImagesPart
                      Router%Put_P(iProcFrom)%Weight_I(iToPut)&
                           =Weight_I(iImages)
                      Router%DoAdd_P(iProcFrom)%&
                           DoAdd_I(iToPut)=.true.
                      nImagesPart=nImagesPart-1
                   end if
                end do
                if(iProc==iProcFrom)then
                   do iProcToPut=1,nProcToPut
                      iProcTo=iProcLookUp_I(iProcToPut)
                      Router%iGet_P(iProcTo)%&
                           iCB_II(1:Router%nIndexesSource,&
                           Router%nGet_P(iProcTo))&
                           =IndexGet_I(1:Router%nIndexesSource)
                      Router%iGet_P(iProcTo)%&
                           iCB_II(0,Router%nGet_P(iProcTo))&
                           =1
                      Router%Get_P(iProcTo)%&
                           Weight_I(Router%nGet_P(iProcTo))&
                           =cOne
                   end do
! The efficient way to properly find doadd's is not found yet
! For a while the target state vector which is filled in                 
! with this subroutine, should be nullified before applying the 
!global message pass  
                end if
             end if
          end do !Global cell
       end do    !Target block
    end do       !Check if DoCountOnly
  end subroutine construct_router_from_source
!====================END========================================!
end Module CON_router
!end Module CON_router
!\end{verbatim}                     !^CFG UNCOMMENT IF PRINTPS  !
!\end{document}                     !^CFG UNCOMMENT IF PRINTPS  !
!=============================LINE 912==========================!


