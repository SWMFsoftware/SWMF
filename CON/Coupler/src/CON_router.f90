!^CFG COPYRIGHT UM                                              !
!BOP
!MODULE: CON_router - set the connection between the grids of different models
!INTERFACE:
Module CON_router
  !USES:
  use CON_grid_descriptor
  use CON_global_vector
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
!\begin{verbatim}
     logical::IsLocal,IsProc
     integer::iProc,nProc,iComm
!\end{verbatim}
!If the union group is constructed, then for use with broadcast !
!we need the union communicator and the root PE ranks in this   !
!communicator
!\begin{verbatim}
     integer::iCommUnion,iProc0Source,iProc0Target
     integer,dimension(:),pointer::iTranslated_P
!\end{verbatim}
!As the default we use iCB indexes to construct the router,     !
!hence the grid point is characterized by the                   !
!GridDescriptor%nDim grid point indexes plus one more index for !
!the block number. Also we allow to use exactly                 ! 
!GridDescriptor%nDim indexes, without the block number which    !
!only seems to be of sence for the component which is localized !
!at one PE only, or which has exactly one block per PE          !
!\begin{verbatim}
      integer::nIndexesSource,nIndexesTarget  
!\end{verbatim}
!The total amounts of the buffer segments to be sent-received   !
!to/from the PE. The total amounts of the grid points from which!
!the data should be got or to which the data should be put,some !
!data points may be counted more than one time                  !
!\begin{verbatim}
     integer, dimension(:), pointer :: &
          nGet_P, nPut_P, nRecv_P, nSend_P
!\end{verbatim}
!iCB indexes and the weight coefficients for the points of the  !
!target and source grids, which are connected through the router!
!\begin{verbatim}
     type(IndexPtrType), dimension(:), pointer :: iGet_P
     type(IndexPtrType), dimension(:), pointer :: iPut_P
     type(DoAddPtrType), dimension(:), pointer :: DoAdd_P
     type(WeightPtrType), dimension(:),pointer :: Get_P
     type(WeightPtrType), dimension(:),pointer :: Put_P
!\end{verbatim}
  end type RouterType
!EOP
  integer,allocatable,dimension(:),save::iAux_P
!BOP
!PUBLIC MEMBER FUNCTIONS:
  private::allocate_get_arrays
  private::allocate_put_arrays
  private::check_router_allocation
  private::iAux_P
!EOP
contains
!BOP
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
          if(.not.allocated(iAux_P))then
             allocate(iAux_P(0:n_proc()-1),stat=iError)
             call check_allocate(iError,'iAux_P')
             do iPE=0,n_proc()-1
                iAux_P(iPE)=iPE
             end do
          end if
          nProc=Router%nProc
          allocate(Router%iTranslated_P(0:nProc-1),stat=iError)
          call check_allocate(iError,'iTranslated_P')
          
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
             call MPI_GROUP_TRANSLATE_RANKS(&
                  i_group(),&
                  n_proc(),&
                  iAux_P(0),&
                  iGroupUnion,&
                  Router%iTranslated_P,&
                  iError)
             Router%iProc0Target=&
                  Router%iTranslated_P(iProc0Target)
         else
             call MPI_GROUP_UNION(&
                  iGroupTarget,&
                  iGroupSource,&
                  iGroupUnion,&
                  iError)
             Router%iProc0Target=0
             call MPI_GROUP_TRANSLATE_RANKS(&
                  i_group(),&
                  n_proc(),&
                  iAux_P(0),&
                  iGroupUnion,&
                  Router%iTranslated_P,&
                  iError)
             Router%iProc0Source=&
                  Router%iTranslated_P(iProc0Source)
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
               'Wrongly defined Router%iProc0Target')
          if(iProcUnion/=Router%iProc0Source.and.&
               i_proc()==iProc0Source)call CON_stop(&
               'Wrongly defined Router%iProc0Source')
          call MPI_GROUP_FREE(iGroupUnion,iError)
       else
          Router%iCommUnion=Router%iComm
          Router%iProc0Source=iProc0Source
          Router%iProc0Target=iProc0Target
       end if
    else
       call CON_stop(&
            'Do not couple a Local grid with a global one')
    end if


    Router%nIndexesSource=ndim_grid(&
       GridDescriptorSource%DD%Ptr)+1
    if(present(nIndexesSource))then
       if(nIndexesSource>=Router%nIndexesSource-1&
            .or.nIndexesSource==1)then
          Router%nIndexesSource=nIndexesSource
       else
          write(*,*)'IndexMin=',Router%nIndexesSource-1
          call CON_stop('nIndexesSource should be at least IndexMin')
       end if
    end if

    Router%nIndexesTarget=ndim_grid(&
       GridDescriptorTarget%DD%Ptr)+1

    if(present(nIndexesTarget))then
       if(nIndexesTarget>=Router%nIndexesTarget-1&
            .or.nIndexesTarget==1)then
          Router%nIndexesTarget=nIndexesTarget
       else
          write(*,*)'IndexMin=',Router%nIndexesTarget-1
          call CON_stop('nIndexesTarget should be at least IndexMin')
       end if
    end if
    nProc=Router%nProc
    
    !Allocation:
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
  subroutine bcast_global_vector_in_router(&
       NameVector,&
       GD,&
       Router,&
       NameMask)
    character(LEN=*),intent(in)::NameVector,NameMask
    optional::NameMask
    type(GridDescriptorType),intent(in)::GD
    type(RouterType),intent(in)::Router

    if(Router%IsLocal)then
       call bcast_global_vector(&
            NameVector,&
            GD,&
            NameMask)
    else
       call bcast_global_vector(&
            NameVector,&
            GD,&
            Router%iTranslated_P,&
            Router%iCommUnion,&
            NameMask)
    end if
  end subroutine bcast_global_vector_in_router

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
!target grid point to an image point into the source domain,    !
!in case this mapping is implemented through a routine          !
       mapping,&
!Mapping transformation which, in the treated case, maps the    !
!target grid point to an image point into the source domain,    !
!in case this mapping is implemented through a global vector    !    
       NameMappingVector,&
!If mapping throught a global vector is used, some points can
!be avoided in setting the router using the named mask array
       NameMask,&
!First or second order interpolation procedures are available   !
!(see CON_grid_descriptor) to find the grid point at the source !
!grid and the interpoltion weights for the image point:         !
!nearest_grid_points and bilinear_interpolation.                !
!Optional, if not present then the nearest_grid_points is used  !
       interpolate)
    use ModIOUnit,ONLY:io_unit_new
    !INPUT ARGUMENTS:

    interface
       logical function is_interface_block(lGlobalNode)
         implicit none
         integer,intent(in)::lGlobalNode 
       end function is_interface_block

       subroutine interface_point_coords(&
                  GridDescriptor,&
                  lGlobalTreeNode,&
                  nDim,&
                  Xyz_D,&
                  nIndexes,&
                  Index_I,&
                  IsInterfacePoint)
         use CON_grid_descriptor
         implicit none
         type(GridDescriptorType),intent(in)::GridDescriptor
         integer,intent(in)::lGlobalTreeNode,nIndexes
         logical,intent(out)::IsInterfacePoint
         integer,intent(in)::nDim
         real,intent(inout)::Xyz_D(nDim)
         integer,intent(inout)::Index_I(nIndexes)
       end subroutine interface_point_coords

       subroutine mapping(&
            nDimFrom,XyzFrom_D,nDimTo,XyzTo_D,IsInterfacePoint)
         implicit none                                                  
         integer,intent(in)::nDimFrom,nDimTo       
         real,dimension(nDimFrom),intent(in)::XyzFrom_D
         real,dimension(nDimTo),intent(out)::XyzTo_D
         logical,intent(out)::IsInterfacePoint
       end subroutine mapping

       subroutine interpolate(&
            nDim,&
            Xyz_D,&
            GridDescriptor,&
            nIndexes,&
            Index_II,&
            nImages,&
            Weight_I)
         use CON_grid_descriptor
         implicit none
         integer,intent(in)::nDim
         real,intent(inout)::Xyz_D(nDim)
         type(GridDescriptorType)::GridDescriptor     
         integer,intent(in)::nIndexes
         integer,dimension(0:nIndexes,2**nDim)::Index_II
         integer,intent(out)::nImages
         real,dimension(2**nDim),intent(out)::Weight_I
       end subroutine interpolate
    end interface

    optional::is_interface_block,interface_point_coords
    optional::mapping,interpolate

    character(LEN=*),intent(in),optional::NameMappingVector
    character(LEN=*),intent(in),optional::NameMask

    type(GridDescriptorType),intent(in):: GridDescriptorSource
    type(GridDescriptorType),intent(in):: GridDescriptorTarget
    type(RouterType),intent(inout)::Router
    !EOP
    integer::iProc,nProc
    integer::lGlobalNode,iBlockAll
    integer::iGlobalGridPoint,nGridPointsPerBlock
    logical::IsInterfacePoint
    integer::iImages,nImages,nImagesPart,iToGet
    integer::iProcTo,iBlockTo,iProcFrom,iProcDoNotAdd,iPE
    integer,dimension(0:Router%nProc-1)::&
         nGetUbound_P,nPutUbound_P

    real,dimension(GridDescriptorTarget%nDim)::XyzTarget_D
    real,dimension(GridDescriptorSource%nDim)::&
                                   XyzSource_D,XyzStored_D
    integer,dimension(GridDescriptorTarget%nDim)::iCell_D
    integer, dimension(Router%nIndexesTarget)::IndexRecv_I
    integer,dimension(0:Router%nIndexesSource,&
                      2**GridDescriptorSource%nDim)::&
                      IndexGet_II
    integer,dimension(2**GridDescriptorSource%nDim)::&
                      iProcLookUp_I
    integer::nProcToGet,iProcToGet
    logical::DoCountOnly,DoCountRed
    real,dimension(2**GridDescriptorSource%nDim)::Weight_I
 
    real,dimension(:,:),pointer::XyzMapping_DI
   
    logical,dimension(:),pointer::Used_I
    logical::UseMask,UseMappingVector,UseMappingFunction

    logical::DoCheckBlock,DoCheckPoint,DoInterpolate
    integer::iError,iFile
    logical::DoTest,DoTestMe
    character(LEN=*),parameter::NameSub='Router'
    character(LEN=100):: NameFile
    !-------------------------
    
    !For given PE the number in the communicator is:
    iProc=Router%iProc
    !

    !Return if the processor does not belong to the communicator
    if(iProc<0)return

    DoTest=.false.; DoTestMe=.false.


    !Check a presence of mapping array.
    !Associate pointer if required.

    UseMappingVector=present(NameMappingVector)

    if(UseMappingVector)&
         call associate_with_global_vector(&
         XyzMapping_DI,NameMappingVector)

    UseMask=present(NameMask)

    if(UseMask)then
       call associate_with_global_mask(&
            Used_I,NameMask)
    end if
    if(UseMask)then
       NameFile = 'router_'//NameMask
       call CON_set_do_test(trim(NameFile),DoTest,DoTestMe)
      
    else
       NameFile ='router_'//NameMappingVector
       if(UseMappingVector)call CON_set_do_test(&
            trim(NameFile),DoTest,DoTestMe)
    end if
    DoTestMe=DoTest.and.iProc==Router%iProc0Target
    if(DoTestMe)write(*,*)'Router starts'
    UseMappingFunction=present(mapping)
  
    if(.not.(UseMappingFunction.or.UseMappingVector).and.&
         GridDescriptorTarget%nDim/=GridDescriptorSource%nDim)&
            call CON_stop(&
            'Mapping is needed for Target%nDim/=Source%nDim')
    nProc=Router%nProc

    DoCheckBlock=present(is_interface_block)
    DoCheckPoint=present(interface_point_coords)
    DoInterpolate=present(interpolate)

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
       if(DoTestMe)then
          iFile=io_unit_new()
          open(iFile,file=trim(NameFile),status='replace')
          write(iFile,*)'iPointGlobal Xyz_D'
          write(iFile,*)'iProcFrom   iCB indexes  Weitht  Sum(Weight)'//&
               'iImages '
       end if


!Initialize the counters                                        !
       Router%nGet_P=0
       Router%nPut_P=0
       Router%nSend_P=0
       Router%nRecv_P=0

       nGridPointsPerBlock=n_grid_points_per_block(&
            GridDescriptorTarget)

!Block loop                                                     !
       do iBlockAll=1,n_block_total(GridDescriptorTarget%DD%Ptr)

          lGlobalNode=i_global_node_a(&
               GridDescriptorTarget%DD%Ptr,iBlockAll)

          call pe_and_blk(&
               GridDescriptorTarget%DD%Ptr,lGlobalNode,&
               iProcTo,iBlockTo)
!Skip the block if desired: if there is known to be no interface!
!point in it                                                    !
          if( DoCheckBlock)then
             if(.not.is_interface_block(lGlobalNode))CYCLE
          end if
!GlobalCellNumber Loop, for a given (octree) block              !
          do iGlobalGridPoint=&
                         1+nGridPointsPerBlock*(iBlockAll-1),&
                         nGridPointsPerBlock*iBlockAll
             
             if(UseMask)then
                if(DoTestMe)&
                     write(iFile,*)'iGlobalPoint=',iGlobalGridPoint,&
                     ' Used_I=', Used_I(iGlobalGridPoint)
                if(.not.Used_I(iGlobalGridPoint))&
                     CYCLE
             end if
             
             IndexRecv_I(1)=iGlobalGridPoint
             if(Router%nIndexesTarget==1.and.&
                  UseMappingVector)then
                   
                XyzSource_D=XyzMapping_DI(&
                     1:GridDescriptorSource%nDim,&
                     iGlobalGridPoint)
             else 
                call global_i_grid_point_to_icb(&
                     GridDescriptorTarget,&
                     iGlobalGridPoint,&
                     lGlobalNode,& 
                     iCell_D)
                
                if(Router%nIndexesTarget/=1)then
                   IndexRecv_I(Router%nIndexesTarget)=iBlockTo
                   IndexRecv_I(1:GridDescriptorTarget%nDim)=&
                        iCell_D
                end if
                if(UseMappingVector)then
                   XyzSource_D=XyzMapping_DI(&
                        1:GridDescriptorSource%nDim,&
                        iGlobalGridPoint)
                else
                   XyzTarget_D=xyz_grid_d(&
                        GridDescriptorTarget,&
                        lGlobalNode,&
                        iCell_D)
                   if( DoCheckPoint)then
                      call interface_point_coords(&
                           GridDescriptorTarget,&
                           lGlobalNode,&
                           GridDescriptorTarget%nDim,&
                           XyzTarget_D,&
                           Router%nIndexesTarget,&
                           IndexRecv_I,&
                           IsInterfacePoint)
                      if(.not.IsInterfacePoint)CYCLE 
                   end if
                   if(UseMappingFunction)then
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
                end if
             end if
             ! call timing_start('set_router_interp')
             if(DoTestMe)then
                XyzStored_D=XyzSource_D
                write(iFile,*)iGlobalGridPoint,XyzSource_D
             end if
             if( DoInterpolate)then
                call interpolate(&
                     GridDescriptorSource%nDim,&
                     XyzSource_D,&
                     GridDescriptorSource,&
                     Router%nIndexesSource,&
                     IndexGet_II,&
                     nImages,&
                     Weight_I)
             else
                call nearest_grid_points(&
                     GridDescriptorSource%nDim,&
                     XyzSource_D,&
                     GridDescriptorSource,&
                     Router%nIndexesSource,&
                     IndexGet_II,&
                     nImages,&
                     Weight_I)
             end if
             if(nImages<1)then
                write(*,*)'nImages=', nImages
                call CON_stop('interpolation failed')
             end if
             if(DoTestMe)then
                do iImages=1,nImages
                   if(iImages==1)then
                      write(iFile,*)IndexGet_II(:,iImages),Weight_I(iImages),&
                           sum(Weight_I(1:nImages))
                   else
                      write(iFile,*)IndexGet_II(:,iImages),Weight_I(iImages),&
                           iImages
                   end if
                end do
                if(Router%nIndexesSource==&
                     GridDescriptorSource%nDim+1)then
                   XyzSource_D=cZero
                   do iImages=1,nImages
                      XyzSource_D=&
                           XyzSource_D+&
                           xyz_grid_d(GridDescriptorSource,&
                           i_global_node_bp(&
                           GridDescriptorSource%DD%Ptr,&
                           IndexGet_II(Router%nIndexesSource,iImages),&
                           IndexGet_II(0,iImages)),&
                           IndexGet_II(1:GridDescriptorSource%nDim,&
                           iImages))*Weight_I(iImages)
                      
                   end do
                   write(iFile,*)'Interpolated coordinate values=',&
                        XyzSource_D,' Error=',&
                        sqrt(sum((XyzSource_D-XyzStored_D)**2))
                end if
                write(iFile,*)
             end if
             ! call timing_stop('set_router_interp')
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

       if(DoTestMe)close(iFile)
    end do       !Check if DoCountOnly
    if(UseMappingVector)nullify(XyzMapping_DI)
    if(UseMask)nullify(Used_I)
  end subroutine set_router
!===============================================================!
!===============================================================!
  subroutine construct_router_from_source(&
       !The Descriptor for Source grid points:   
       GridDescriptorSource,& !<<<<<<<<<<<<<<
       !The Descriptor for Target grid points:   
       GridDescriptorTarget,& !<<<<<<<<<<<<<<
       !The Router to be set 
       Router,&               !<<<<<<<<<<<<<<
       !Logical function which allows to skip the block if there is no !
       !interface points in it. Optional, if not present then all the  !
       !blocks are checked for the presence of the interface points    !
       is_interface_block,&   !<<<<<<<<<<<<<<
       !The subroutine which defines if the grid point is inside the   !
       !interface layer. Optional, if not present, then all the grid   !
       !points (at the source grid) are considered as the interface    !
       !layer points 
       interface_point_coords, &  !<<<<<<<<<<
       !Mapping transformation which, in the treated case, maps the    !
       !source grid point to an image point into the target domain,    !
       !in case this mapping is implemented through a routine          !
       mapping,&              !<<<<<<<<<<<<<<
       !Mapping transformation which, in the treated case, maps the    !
       !source grid point to an image point into the target domain,    !
       !in case this mapping is implemented through a global vector    !    
       NameMappingVector,&    !<<<<<<<<<<<<<<
       !If mapping throught a global vector is used, some points can
       !be avoided in setting the router using the named mask array
       NameMask,&             !<<<<<<<<<<<<<<
       !First or second order interpolation procedures are available   !
       !(see CON_grid_descriptor) to find the grid point at the target !
       !grid and the interpoltion weights for the image point:         !
       !nearest_grid_points and bilinear_interpolation.                !
       !Optional, if not present then the nearest_grid_points is used  !
       interpolate)           !<<<<<<<<<<<<<<
    use ModIOUnit,ONLY:io_unit_new

    !INPUT ARGUMENTS:

    interface
       logical function is_interface_block(lGlobalNode)
         implicit none
         integer,intent(in)::lGlobalNode 
       end function is_interface_block
       
       subroutine interface_point_coords(&
                  GridDescriptor,&
                  lGlobalTreeNode,&
                  nDim,&
                  Xyz_D,&
                  nIndexes,&
                  Index_I,&
                  IsInterfacePoint)
         use CON_grid_descriptor
         implicit none
         type(GridDescriptorType),intent(in)::GridDescriptor
         integer,intent(in)::lGlobalTreeNode,nIndexes
         logical,intent(out)::IsInterfacePoint
         integer,intent(in)::nDim
         real,intent(inout)::Xyz_D(nDim)
         integer,intent(inout)::Index_I(nIndexes)
       end subroutine interface_point_coords

       subroutine mapping(&
            nDimFrom,XyzFrom_D,nDimTo,XyzTo_D,IsInterfacePoint)
         implicit none                                                  
         integer,intent(in)::nDimFrom,nDimTo       
         real,dimension(nDimFrom),intent(in)::XyzFrom_D
         real,dimension(nDimTo),intent(out)::XyzTo_D
         logical,intent(out)::IsInterfacePoint
       end subroutine mapping

       subroutine interpolate(&
            nDim,&
            Xyz_D,&
            GridDescriptor,&
            nIndexes,&
            Index_II,&
            nImages,&
            Weight_I)
         use CON_grid_descriptor
         implicit none
         integer,intent(in)::nDim
         real,intent(inout)::Xyz_D(nDim)
         type(GridDescriptorType)::GridDescriptor     
         integer,intent(in)::nIndexes
         integer,dimension(0:nIndexes,2**nDim)::Index_II
         integer,intent(out)::nImages
         real,dimension(2**nDim),intent(out)::Weight_I
       end subroutine interpolate
    end interface

    optional::is_interface_block,interface_point_coords
    optional::mapping,interpolate

    character(LEN=*),intent(in),optional::NameMappingVector
    character(LEN=*),intent(in),optional::NameMask

    type(GridDescriptorType),intent(in):: GridDescriptorSource
    type(GridDescriptorType),intent(in):: GridDescriptorTarget
    type(RouterType),intent(inout)::Router
    !EOP
    integer::iProc,nProc
    integer::lGlobalNode,iBlockAll
    integer::iGlobalGridPoint,nGridPointsPerBlock
    logical::IsInterfacePoint
    integer::iImages,nImages,nImagesPart,iToPut
    integer::iProcTo,iBlockFrom,iProcFrom,iPE
    integer,dimension(0:Router%nProc-1)::&
         nGetUbound_P,nPutUbound_P

    real,dimension(GridDescriptorTarget%nDim)::&
                                  XyzTarget_D, XyzStored_D
    real,dimension(GridDescriptorSource%nDim)::XyzSource_D

    integer,dimension(GridDescriptorSource%nDim)::iCell_D
    integer, dimension(Router%nIndexesSource)::IndexGet_I

    integer,dimension(0:Router%nIndexesTarget,&
                      2**GridDescriptorTarget%nDim)::&
                      IndexPut_II

    integer,dimension(2**GridDescriptorTarget%nDim)::&
                      iProcLookUp_I

    integer::nProcToPut,iProcToPut
    logical::DoCountOnly,DoCountRed

    real,dimension(2**GridDescriptorTarget%nDim)::Weight_I 

    real,dimension(:,:),pointer::XyzMapping_DI

    logical,dimension(:),pointer::Used_I
    logical::UseMask,UseMappingVector,UseMappingFunction
    logical::DoInterpolate,DoCheckBlock,DoCheckPoint
    integer::iError,iFile
    logical::DoTest,DoTestMe
    
    character(LEN=100):: NameFile
    !------------------------
    
    !For given PE the number in the communicator is:
    iProc=Router%iProc
    !

    !Return if the processor does not belong to the communicator
    if(iProc<0)return

    DoTest=.false.; DoTestMe=.false.


    !Check a presence of mapping array.
    !Associate pointer if required.
    
    UseMappingVector=present(NameMappingVector)

    if(UseMappingVector)&
         call associate_with_global_vector(&
         XyzMapping_DI,NameMappingVector)

    UseMask=present(NameMask)

    if(UseMask)then
       call associate_with_global_mask(&
            Used_I,NameMask)
    end if
    if(UseMask)then
       NameFile = 'router_'//NameMask
       call CON_set_do_test(trim(NameFile),DoTest,DoTestMe)
       NameFile = 'from_source_'//NameFile
    else
       NameFile ='router_'//NameMappingVector
       if(UseMappingVector)call CON_set_do_test(&
            trim(NameFile),DoTest,DoTestMe)
       NameFile = 'from_source_'//NameFile
    end if
    
    DoTestMe=DoTest.and.iProc==Router%iProc0Source
    if(DoTestMe)write(*,*)'Router from source starts'

    UseMappingFunction=present(mapping)
    
    if(.not.(UseMappingFunction.or.UseMappingVector).and.&
         GridDescriptorTarget%nDim/=GridDescriptorSource%nDim)&
            call CON_stop(&
            'Mapping is needed for Target%nDim/=Source%nDim')
    nProc=Router%nProc

    DoCheckBlock=present(is_interface_block)
    DoCheckPoint=present(interface_point_coords)
    DoInterpolate=present(interpolate)

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
       !DoCountOnly will be set to true. The loop then will be repeated!
       ! for the second time

       if(DoTestMe)then
          iFile=io_unit_new()
          open(iFile,file=NameFile,status='replace')
          write(iFile,*)'iPointGlobal Xyz_D'
          write(iFile,*)'iProcTo   iCB indexes  Weitht  Sum(Weight)'//&
               'iImages '
       end if


       !Initialize the counters
       Router%nGet_P=0
       Router%nPut_P=0
       Router%nSend_P=0
       Router%nRecv_P=0

       nGridPointsPerBlock=n_grid_points_per_block(&
            GridDescriptorSource)

       !Block loop
       do iBlockAll=1,n_block_total(GridDescriptorSource%DD%Ptr)

          !Skip non-end octree nodes, if any                              !
          lGlobalNode=i_global_node_a(&
               GridDescriptorSource%DD%Ptr,iBlockAll)

          call pe_and_blk(&
               GridDescriptorSource%DD%Ptr,lGlobalNode,&
               iProcFrom,iBlockFrom)
               

          !Skip the block if desired: if there is known to be no interface!
          ! point in it
          if( DoCheckBlock)then
             if(.not.is_interface_block(lGlobalNode))CYCLE
          end if

          !GlobalCellNumber Loop, for a given (octree) block              !
          do iGlobalGridPoint=&
                         1+nGridPointsPerBlock*(iBlockAll-1),&
                         nGridPointsPerBlock*iBlockAll
             if(UseMask)then
                if(DoTestMe)&
                     write(iFile,*)'iGlobalPoint=',iGlobalGridPoint,&
                     ' Used_I=', Used_I(iGlobalGridPoint)
                if(.not.Used_I(iGlobalGridPoint))&
                     CYCLE
             end if
             
             !Treat separately a case when
             !a single-index global vector 
             !is a source
             
             IndexGet_I(1) = iGlobalGridPoint
             !This value will be rewritten otherwise
             
             if(Router%nIndexesSource==1.and.&
                  UseMappingVector)then
                XyzTarget_D=XyzMapping_DI(&
                     1:GridDescriptorTarget%nDim,&
                     iGlobalGridPoint)
             else
                call global_i_grid_point_to_icb(&
                     GridDescriptorSource,&
                     iGlobalGridPoint,&
                     lGlobalNode, &
                     iCell_D) 
                if(Router%nIndexesSource/=1)then
                   IndexGet_I(Router%nIndexesSource)=iBlockFrom
                   IndexGet_I(1:GridDescriptorSource%nDim) =&
                        iCell_D
                end if
                
                
                if(UseMappingVector)then
                   XyzTarget_D=XyzMapping_DI(&
                        1:GridDescriptorTarget%nDim,&
                        iGlobalGridPoint) 
                else 
                   XyzSource_D=xyz_grid_d(&
                        GridDescriptorSource,&
                        lGlobalNode,&
                        IndexGet_I(1:GridDescriptorSource%nDim))
                   
                   if( DoCheckPoint)then
                      call interface_point_coords(&
                           GridDescriptorSource,&
                           lGlobalNode,&
                           GridDescriptorSource%nDim,&
                           XyzSource_D,&
                           Router%nIndexesSource,&
                           IndexGet_I,&
                           IsInterfacePoint)
                      if(.not.IsInterfacePoint)CYCLE 
                   end if
                   if(UseMappingFunction)then
                      call mapping(GridDescriptorSource%nDim,&
                           XyzSource_D,&
                           GridDescriptorTarget%nDim,&
                           XyzTarget_D,&
                           IsInterfacePoint)
                      if(.not.IsInterfacePoint)CYCLE 
                   else
                      XyzTarget_D=XyzSource_D
                   end if   !Mapping function
                end if  !MappingVector
             end if  !Global vector as a source
             if(DoTestMe)then
                XyzStored_D=XyzTarget_D
                write(iFile,*)iGlobalGridPoint,XyzTarget_D
             end if
             if( DoInterpolate)then
                call interpolate(&
                     GridDescriptorTarget%nDim,&
                     XyzTarget_D,&
                     GridDescriptorTarget,&
                     Router%nIndexesTarget,&
                     IndexPut_II,&
                     nImages,&
                     Weight_I)
             else
                call nearest_grid_points(&
                     GridDescriptorTarget%nDim,&   
                     XyzTarget_D,&
                     GridDescriptorTarget,&
                     Router%nIndexesTarget,&
                     IndexPut_II,&
                     nImages,&
                     Weight_I)
             end if
             if(nImages<1)then
                write(*,*)'nImages=', nImages
                call CON_stop('interpolation failed in router from source')
             end if
             if(DoTestMe)then
                do iImages=1,nImages
                   if(iImages==1)then
                      write(iFile,*)IndexPut_II(:,iImages),Weight_I(iImages),&
                           sum(Weight_I(1:nImages))
                   else
                      write(iFile,*)IndexPut_II(:,iImages),Weight_I(iImages),&
                           iImages
                   end if
                end do
                if(Router%nIndexesTarget==&
                     GridDescriptorSource%nDim+1)then
                   XyzTarget_D=cZero
                   do iImages=1,nImages
                      XyzTarget_D=&
                           XyzTarget_D+&
                           xyz_grid_d(GridDescriptorTarget,&
                           i_global_node_bp(&
                           GridDescriptorTarget%DD%Ptr,&
                           IndexPut_II(Router%nIndexesTarget,iImages),&
                           IndexPut_II(0,iImages)),&
                           IndexPut_II(1:GridDescriptorTarget%nDim,&
                           iImages))*Weight_I(iImages)
                      
                   end do
                   write(iFile,*)'Interpolated coordinate values=',&
                        XyzTarget_D,' Error=',&
                        sqrt(sum((XyzTarget_D-XyzStored_D)**2))
                end if
                write(iFile,*)
             end if
             !See interface
             !--------------------------------------!
             !Lookup

             nImagesPart=0     !At all CPUs

             do iImages=1,nImages
                iProcTo=IndexPut_II(0,iImages)
                !At the target PEs the number of terms 
                !in the partial sums are found

                if(iProc==iProcTo)then
                   nImagesPart=nImagesPart+1
                   Router%nPut_P(iProcFrom)=&
                        Router%nPut_P(iProcFrom)+1
                   DoCountOnly=DoCountOnly.or.&
                        Router%nPut_P(iProcFrom)>&
                        nPutUbound_P(iProcFrom)
                end if

                !At the source processor the PE list is defined !
                !which will get the  partial sums
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
       if(DoTestMe)close(iFile)
    end do       !Check if DoCountOnly
    if(UseMappingVector)nullify(XyzMapping_DI)
    if(UseMask)nullify(Used_I)
  end subroutine construct_router_from_source
!====================END========================================!
end Module CON_router
!end Module CON_router
!\end{verbatim}                     !^CFG UNCOMMENT IF PRINTPS  !
!\end{document}                     !^CFG UNCOMMENT IF PRINTPS  !
!=============================LINE 912==========================!


