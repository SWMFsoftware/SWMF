Module CON_test_global_message_pass
  use CON_global_message_pass
  implicit none
  private !Except
  public::test_global_message_pass
  integer,parameter::nGhostCells=2
  real, dimension(:,:,:,:), allocatable :: V0,V1
  real, dimension(:,:,:), allocatable :: V20,V21
  interface test_global_message_pass
     module procedure test_global_message_pass_dd
     module procedure test_global_message_pass_id
  end interface
contains
  subroutine test_global_message_pass_id(GridID_)
    integer,intent(in)::GridID_
    select case(ndim_grid(GridID_))
    case(2)
       call test_global_message_pass_2id(GridID_)
    case(3)
       call test_global_message_pass_3id(GridID_)
    case default
       call CON_stop(&
            'There is no means for testing the grid of dimensions of',&
            ndim_grid(GridID_))
    end select
  end subroutine test_global_message_pass_id
  subroutine test_global_message_pass_dd(DomainDecomposition)
    type(DomainDecompositionType),intent(inout),target::DomainDecomposition
    select case(ndim_grid(DomainDecomposition))
    case(2)
       call test_global_message_pass_2dd(DomainDecomposition)
    case(3)
       call test_global_message_pass_3dd(DomainDecomposition)
    case default
       call CON_stop(&
            'There is no means for testing the grid of dimensions of',&
            ndim_grid(DomainDecomposition))
    end select
  end subroutine test_global_message_pass_dd
!===============================================================!
! TESTS TESTS TESTS  TESTS TESTS TESTS  TESTS TESTS TESTS  TESTS!
!===============================================================!
! For three dimensional grid
  subroutine test_global_message_pass_3dd(Grid)
    type(DomainDecompositionType),intent(inout),target::Grid

    !Local variables
    integer :: i,j,k,iBlockStart, iBLK,nBLK
    integer :: iError,lOctree, iImages
    integer, parameter :: nDim=3
    real, dimension(nDim),parameter :: F_D=(/1.0,0.71,0.13/)
    real::Xyz_D(nDim),Weight_I(8),VMisc
    integer::nImages,lSearch
    integer,dimension(nDim)::iCells_D,nCells_D
    integer,dimension(0:4,8)::Indexes_II
    type(RouterType)::Router
    type(GridDescriptorType)::Cell2,Cell,Temp 
    !------------------------------------------
    if(.not.is_proc())return
    call set_standard_grid_descriptor(Grid,GridDescriptor=Cell)
    call set_standard_grid_descriptor(Grid,nGhostGridPoints=nGhostCells,&
         GridDescriptor=Cell2)
    call init_router&
         (Cell,Cell2,Router)
    nCells_D=ncells_decomposition_d(Grid)
    Xyz_D=cHalf*(xyz_max_d(Grid)+xyz_min_d(Grid))
    if(is_proc0())write(*,*)"XyzTest",Xyz_D
    call set_standard_grid_descriptor(Grid,GridDescriptor=Temp)
    call nearest_grid_points(Xyz_D,Temp,4,Indexes_II,nImages,Weight_I)
      if(is_proc0())then
       do iImages=1,nImages
          write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
       end do
       write(*,*)'Nearest Cell, nImages=',nImages,'Control sum:',&
            sum(Weight_I(1:nImages))
    end if
    Xyz_D=cHalf*(xyz_max_d(Grid)+xyz_min_d(Grid))
    call interpolation_fix_reschange(nDim,Xyz_D,Temp,4,Indexes_II,nImages,Weight_I)
      if(is_proc0())then
       do iImages=1,nImages
          write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
       end do
       write(*,*)'Bilinear, Cells, nImages=',nImages,'Control sum:',&
            sum(Weight_I(1:nImages))
    end if
    Xyz_D=cHalf*(xyz_max_d(Grid)+xyz_min_d(Grid))
    if(is_proc0())write(*,*)"XyzTest",Xyz_D
    call clean_grid_descriptor(Temp)
    call set_standard_grid_descriptor(Grid,&
         Standard_=Nodes_,&
         GridDescriptor=Temp)
    call nearest_grid_points(Xyz_D,Temp,4,Indexes_II,nImages,Weight_I)
      if(is_proc0())then
       do iImages=1,nImages
          write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
       end do
       write(*,*)'NearestNode, nImages=',nImages,'Control sum:',&
            sum(Weight_I(1:nImages))
    end if
    Xyz_D=cHalf*(xyz_max_d(Grid)+xyz_min_d(Grid))
    call interpolation_fix_reschange(nDim,Xyz_D,Temp,4,Indexes_II,nImages,Weight_I)
      if(is_proc0())then
       do iImages=1,nImages
          write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
       end do
       write(*,*)'Bilinear,nodes, nImages=',nImages,'Control sum:',&
            sum(Weight_I(1:nImages))
    end if
    call clean_grid_descriptor(Temp)
    call set_router(&
         GridDescriptorSource=Cell,&
         GridDescriptorTarget=Cell2,&
         Router=Router,&
         interface_point_coords=do_all_ghostcells) 
    nBLK=max_block_pe(Grid,i_proc())
    iBlockStart=min_block_pe(Grid,i_proc())
 
    if(is_proc0())then
       write(*,*)'For Grid of the component ',compid_grid(Grid)
       write(*,*)'Testing message_pass, PE=',i_proc(),'  Starting tests ...'
    end if
    allocate( V0(1-nGhostCells:&
         nCells_D(1)+nGhostCells,&
         1-nGhostCells:&
         nCells_D(2)+nGhostCells,&
         1-nGhostCells:&
         nCells_D(3)+nGhostCells,&
         iBlockStart:nBLK), stat=iError )
    call check_allocate(iError,'V0')
    allocate( V1(1-nGhostCells:&
         nCells_D(1)+nGhostCells,&
         1-nGhostCells:&
         nCells_D(2)+nGhostCells,&
         1-nGhostCells:&
         nCells_D(3)+nGhostCells,&
         iBlockStart:nBLK), stat=iError )
    call check_allocate(iError,'V0')

    V0=-999999.
    V1=V0
    do lOctree=1,ntree_nodes(Grid)
       if(.not.used_node(Grid,lOctree))CYCLE
       if(pe_decomposition(Grid,lOctree)/=i_proc())CYCLE
       iBLK=blk_decomposition(Grid,lOctree)
!V0 is assigned througout all cells
       do k=1-nGhostCells,nCells_D(3)+nGhostCells
          iCells_D(3)=k
          do j=1-nGhostCells,nCells_D(2)+nGhostCells
             iCells_D(2)=j
             do i=1-nGhostCells,nCells_D(1)+nGhostCells
                iCells_D(1)=i
                V0(i,j,k,iBLK)=&
                     dot_product(F_D,xyz_cell_d(Grid,lOctree,iCells_D))
             end do
          end do
       end do
!In sending messages in the ghostcells which cover coarser block
!message pass writes equal values to all the ghostcells which are
!covered by common coarser physical cell. To mimic this effect
!the data in the array V0 are coarsened too in such the ghostcells
       do k=1-nGhostCells,nCells_D(3)+nGhostCells,2
          iCells_D(3)=k
          do j=1-nGhostCells,nCells_D(2)+nGhostCells,2 
             iCells_D(2)=j
             do i=1-nGhostCells,nCells_D(1)+nGhostCells,2
                iCells_D(1)=i
                if(l_level_neighbor(Grid,lOctree,iCells_D)==1)then
                   V0(i:i+1,j:j+1,k:k+1,iBLK)=&
                        cEighth*sum(V0(i:i+1,j:j+1,k:k+1,iBLK))           
                end if
             end do
          end do
       end do
!V1 is assigned in true cells only
       do k=1,nCells_D(3)
          iCells_D(3)=k
          do j=1,nCells_D(2)
             iCells_D(2)=j
             do i=1,nCells_D(1)
                iCells_D(1)=i
                V1(i,j,k,iBLK)=&
                     dot_product(F_D,xyz_cell_d(Grid,lOctree,iCells_D))
             end do
          end do
       end do
       if (any(is_left_boundary_d(Grid,lOctree))&
            .or.any(is_right_boundary_d(Grid,lOctree)))&
            V1(:,:,:,iBLK)=V0(:,:,:,iBLK)       
    end do

!Ghost cell values for V1 are filled in using global message pass
    call global_message_pass(Router,&
                             1,&      !
                             getV1,&  !Wrapper for Source
                             putV1, & !Wrapper for Target
                             copyV1)  !Wrapper for Source and Target
!Compare V0 and V1
    if(max(maxval(V1-V0),maxval(V0-V1))>.01)then
       write(*,*)'Testing message_pass, PE=',i_proc(),&
            ' max difference=', &
            max(maxval(V1-V0),maxval(V0-V1)), &
            ' printing out values and exiting.'
       write(*,*)' '
       do lOctree=1,ntree_nodes(Grid)
          if(.not.used_node(Grid,lOctree))CYCLE
          if(pe_decomposition(Grid,lOctree)/=i_proc())CYCLE
          iBLK=blk_decomposition(Grid,lOctree)
          write(*,*)xyz_block_d(Grid,lOctree)
          do k=1-nGhostCells,nCells_D(3)+nGhostCells
             do j=1-nGhostCells,nCells_D(2)+nGhostCells 
                do i=1-nGhostCells,nCells_D(1)+nGhostCells
                   if(abs(V0(i,j,k,iBLK)-V1(i,j,k,iBLK))>.01)&
                        write(*,*)iBLK,i,j,k,V0(i,j,k,iBLK),V1(i,j,k,iBLK)
                end do
             end do
          end do
       end do
       write(*,*)' '
       call CON_stop('test_global_message_pass_failed')
    end if

    deallocate(V0)
    deallocate(V1)

    call clean_router(Router)
    call clean_grid_descriptor(Cell)
    call clean_grid_descriptor(Cell2)
    if(is_proc0())then
       write(*,*)'Testing message_pass, PE=',i_proc(), &
            '  All tests passed.'
       write(*,*)' '
    end if

  end subroutine test_global_message_pass_3dd
!===============================================================
! For three dimensional grid
  subroutine test_global_message_pass_3id(GridID_)
    integer,intent(in)::GridID_

    !Local variables
    integer :: i,j,k,iBlockStart, iBLK,nBLK
    integer :: iError,lOctree, iImages
    integer, parameter :: nDim=3
    real, dimension(nDim),parameter :: F_D=(/1.0,0.71,0.13/)
    real::Xyz_D(nDim),Weight_I(8),VMisc
    integer::nImages,lSearch,iPoint
    integer,dimension(nDim)::iCells_D,nCells_D
    integer,dimension(0:4,8)::Indexes_II
    type(RouterType)::Router
    type(GridDescriptorType)::Cell2,Cell,Temp 
    !------------------------------------------
    if(.not.is_proc())return
    call set_standard_grid_descriptor(GridID_,GridDescriptor=Cell)
    call set_standard_grid_descriptor(GridID_,nGhostGridPoints=nGhostCells,&
         GridDescriptor=Cell2)
    call init_router&
         (Cell,Cell2,Router)
    nCells_D=ncells_decomposition_d(GridID_)
    if(is_proc0())write(*,*)
    if(is_proc0())write(*,*)'Testing cell-centered grid'
    if(is_proc0())write(*,*)
    call set_standard_grid_descriptor(GridID_,GridDescriptor=Temp)
    do iPoint=1,35
       if(is_proc0())write(*,*)'Testing nearest_grid_points routine:'
       Xyz_D=((iPoint-2)*xyz_max_d(GridID_)+&
            (34-iPoint)*xyz_min_d(GridID_))/32
       if(is_proc0())write(*,*)"XyzTest",Xyz_D

       call nearest_grid_points(Xyz_D,Temp,4,Indexes_II,nImages,Weight_I)
       if(is_proc0())then
          do iImages=1,nImages
             write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
          end do
          write(*,*)'nImages=',nImages,'Control sum:',&
               sum(Weight_I(1:nImages))
       end if
    end do
    do iPoint=1,35
       if(is_proc0())write(*,*)'Testing bilinear_interpolation routine:'
       Xyz_D=((iPoint-2)*xyz_max_d(GridID_)+&
            (34-iPoint)*xyz_min_d(GridID_))/32
       if(is_proc0())write(*,*)"XyzTest",Xyz_D
       call bilinear_interpolation(nDim,Xyz_D,Temp,4,Indexes_II,nImages,Weight_I)
       if(is_proc0())then
          do iImages=1,nImages
             write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
          end do
          write(*,*)'nImages=',nImages,'Control sum:',&
               sum(Weight_I(1:nImages))
       end if
    end do

    do iPoint=1,35
       if(is_proc0())write(*,*)'Testing interpolation_fix_reschange routine:'
       Xyz_D=((iPoint-2)*xyz_max_d(GridID_)+&
            (34-iPoint)*xyz_min_d(GridID_))/32
       if(is_proc0())write(*,*)"XyzTest",Xyz_D
       call interpolation_fix_reschange(nDim,Xyz_D,Temp,4,Indexes_II,nImages,Weight_I)
       if(is_proc0())then
          do iImages=1,nImages
             write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
          end do
          write(*,*)'nImages=',nImages,'Control sum:',&
               sum(Weight_I(1:nImages))
       end if
    end do


    call clean_grid_descriptor(Temp)
    call set_standard_grid_descriptor(GridID_,&
         Standard_=Nodes_,&
         GridDescriptor=Temp)
    if(is_proc0())write(*,*)
    if(is_proc0())write(*,*)'Testing node grid'
    if(is_proc0())write(*,*)
    do iPoint=1,35
       if(is_proc0())write(*,*)'Testing nearest_grid_points routine:'
       Xyz_D=((iPoint-2)*xyz_max_d(GridID_)+&
            (34-iPoint)*xyz_min_d(GridID_))/32
       if(is_proc0())write(*,*)"XyzTest",Xyz_D
       call nearest_grid_points(Xyz_D,Temp,4,Indexes_II,nImages,Weight_I)
       if(is_proc0())then
          do iImages=1,nImages
             write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
          end do
          write(*,*)'nImages=',nImages,'Control sum:',&
               sum(Weight_I(1:nImages))
       end if
    end do
    do iPoint=1,35
       if(is_proc0())write(*,*)'Testing bilinear_interpolation routine:'
       Xyz_D=((iPoint-2)*xyz_max_d(GridID_)+&
            (34-iPoint)*xyz_min_d(GridID_))/32
       if(is_proc0())write(*,*)"XyzTest",Xyz_D
       call bilinear_interpolation(nDim,Xyz_D,Temp,4,Indexes_II,nImages,Weight_I)
       if(is_proc0())then
          do iImages=1,nImages
             write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
          end do
          write(*,*)'nImages=',nImages,'Control sum:',&
               sum(Weight_I(1:nImages))
       end if
    end do
    do iPoint=1,35
       if(is_proc0())write(*,*)'Testing interpolation_fix_reschange routine:'
       Xyz_D=((iPoint-2)*xyz_max_d(GridID_)+&
            (34-iPoint)*xyz_min_d(GridID_))/32 
       if(is_proc0())write(*,*)"XyzTest",Xyz_D
       call interpolation_fix_reschange(nDim,Xyz_D,Temp,4,Indexes_II,nImages,Weight_I)
       if(is_proc0())then
          do iImages=1,nImages
             write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
          end do
          write(*,*)'nImages=',nImages,'Control sum:',&
               sum(Weight_I(1:nImages))
       end if
    end do
    call clean_grid_descriptor(Temp)
    call set_router(&
         GridDescriptorSource=Cell,&
         GridDescriptorTarget=Cell2,&
         Router=Router,&
         interface_point_coords=do_all_ghostcells) 

    nBLK=max_block_pe(GridID_,i_proc())
    iBlockStart=min_block_pe(GridID_,i_proc())
    if(is_proc0())then
       write(*,*)'For Grid of the component ',compid_grid(GridID_)
       write(*,*)'Testing message_pass, PE=',i_proc(),'  Starting tests ...'
    end if
    allocate( V0(1-nGhostCells:&
         nCells_D(1)+nGhostCells,&
         1-nGhostCells:&
         nCells_D(2)+nGhostCells,&
         1-nGhostCells:&
         nCells_D(3)+nGhostCells,&
         iBlockStart:nBLK), stat=iError )
    call check_allocate(iError,'V0')
    allocate( V1(1-nGhostCells:&
         nCells_D(1)+nGhostCells,&
         1-nGhostCells:&
         nCells_D(2)+nGhostCells,&
         1-nGhostCells:&
         nCells_D(3)+nGhostCells,&
         iBlockStart:nBLK), stat=iError )
    call check_allocate(iError,'V0')

    V0=-999999.
    V1=V0
    do lOctree=1,ntree_nodes(GridID_)
       if(.not.used_node(GridID_,lOctree))CYCLE
       if(pe_decomposition(GridID_,lOctree)/=i_proc())CYCLE
       iBLK=blk_decomposition(GridID_,lOctree)
!V0 is assigned througout all cells
       do k=1-nGhostCells,nCells_D(3)+nGhostCells
          iCells_D(3)=k
          do j=1-nGhostCells,nCells_D(2)+nGhostCells
             iCells_D(2)=j
             do i=1-nGhostCells,nCells_D(1)+nGhostCells
                iCells_D(1)=i
                V0(i,j,k,iBLK)=&
                     dot_product(F_D,xyz_cell_d(GridID_,lOctree,iCells_D))
             end do
          end do
       end do
!In sending messages in the ghostcells which cover coarser block
!message pass writes equal values to all the ghostcells which are
!covered by common coarser physical cell. To mimic this effect
!the data in the array V0 are coarsened too in such the ghostcells
       do k=1-nGhostCells,nCells_D(3)+nGhostCells,2
          iCells_D(3)=k
          do j=1-nGhostCells,nCells_D(2)+nGhostCells,2 
             iCells_D(2)=j
             do i=1-nGhostCells,nCells_D(1)+nGhostCells,2
                iCells_D(1)=i
                if(l_level_neighbor(GridID_,lOctree,iCells_D)==1)then
                   V0(i:i+1,j:j+1,k:k+1,iBLK)=&
                        cEighth*sum(V0(i:i+1,j:j+1,k:k+1,iBLK))           
                end if
             end do
          end do
       end do
!V1 is assigned in true cells only
       do k=1,nCells_D(3)
          iCells_D(3)=k
          do j=1,nCells_D(2)
             iCells_D(2)=j
             do i=1,nCells_D(1)
                iCells_D(1)=i
                V1(i,j,k,iBLK)=&
                     dot_product(F_D,xyz_cell_d(GridID_,lOctree,iCells_D))
             end do
          end do
       end do
       if (any(is_left_boundary_d(GridID_,lOctree))&
            .or.any(is_right_boundary_d(GridID_,lOctree)))&
            V1(:,:,:,iBLK)=V0(:,:,:,iBLK)       
    end do

!Ghost cell values for V1 are filled in using global message pass
    call global_message_pass(Router,&
                             1,&      !
                             getV1,&  !Wrapper for Source
                             putV1, & !Wrapper for Target
                             copyV1)  !Wrapper for Source and Target
!Compare V0 and V1
    if(max(maxval(V1-V0),maxval(V0-V1))>.01)then
       write(*,*)'Testing message_pass, PE=',i_proc(),&
            ' max difference=', &
            max(maxval(V1-V0),maxval(V0-V1)), &
            ' printing out values and exiting.'
       write(*,*)' '
       do lOctree=1,ntree_nodes(GridID_)
          if(.not.used_node(GridID_,lOctree))CYCLE
          if(pe_decomposition(GridID_,lOctree)/=i_proc())CYCLE
          iBLK=blk_decomposition(GridID_,lOctree)
          write(*,*)xyz_block_d(GridID_,lOctree)
          do k=1-nGhostCells,nCells_D(3)+nGhostCells
             do j=1-nGhostCells,nCells_D(2)+nGhostCells 
                do i=1-nGhostCells,nCells_D(1)+nGhostCells
                   if(abs(V0(i,j,k,iBLK)-V1(i,j,k,iBLK))>.01)&
                        write(*,*)iBLK,i,j,k,V0(i,j,k,iBLK),V1(i,j,k,iBLK)
                end do
             end do
          end do
       end do
       write(*,*)' '
       call CON_stop('test_global_message_pass_failed')
    end if

    deallocate(V0)
    deallocate(V1)

    call clean_router(Router)
    call clean_grid_descriptor(Cell)
    call clean_grid_descriptor(Cell2)
    if(is_proc0())then
       write(*,*)'Testing message_pass, PE=',i_proc(), &
            '  All tests passed.'
       write(*,*)' '
    end if

  end subroutine test_global_message_pass_3id
!===============================================================
!Subroutines used for testing global message pass
  subroutine getV1(nPartial,&
                  iGetStart,&
                  Get,&
                  Weight,&
                  Buff_I,nVar)
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get    
    type(WeightPtrType),intent(in)::Weight
    real,dimension(nVar),intent(out)::Buff_I
    integer::iGet
    Buff_I(1)=V1(Get%iCB_II(1,iGetStart),&
                 Get%iCB_II(2,iGetStart),&
                 Get%iCB_II(3,iGetStart),&
                 Get%iCB_II(4,iGetStart))*&
                 Weight%Weight_I(iGetStart)
    do iGet=iGetStart+1,iGetStart+nPartial-1
       Buff_I(1)= Buff_I(1)+&
                  V1(Get%iCB_II(1,iGet),&
                     Get%iCB_II(2,iGet),&
                     Get%iCB_II(3,iGet),&
                     Get%iCB_II(4,iGet))*&
                     Weight%Weight_I(iGet)
    end do
  end subroutine getV1
!==============================================================
  subroutine putV1(nPartial,&
                   iPut,&
                   Put,&
                   Weight,&
                   DoAdd,Buff_I,nVar)
    integer,intent(in)::iPut,nPartial,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::Weight
    real,dimension(nVar),intent(in)::Buff_I
    logical,intent(in)::DoAdd
    if(DoAdd)then
       V1(Put%iCB_II(1,iPut),&
          Put%iCB_II(2,iPut),&
          Put%iCB_II(3,iPut),& 
          Put%iCB_II(4,iPut))=&
          V1(Put%iCB_II(1,iPut),&
             Put%iCB_II(2,iPut),&
             Put%iCB_II(3,iPut),& 
             Put%iCB_II(4,iPut))+Buff_I(1)
    else
       V1(Put%iCB_II(1,iPut),&
          Put%iCB_II(2,iPut),&
          Put%iCB_II(3,iPut),& 
          Put%iCB_II(4,iPut))=Buff_I(1)
    end if
  end subroutine putV1
!==============================================================
  subroutine copyV1(nPartialGet,&
                    iGetStart,  &
                    Get,&
                    Weight,&
                    nPartialPut,&
                    iPutStart,&
                    Put,Weight1,DoAdd)
    integer,intent(in)::nPartialGet,iGetStart
    integer,intent(in)::nPartialPut,iPutStart
    type(IndexPtrType),intent(in)::Get,Put
    type(WeightPtrType),intent(in)::Weight,Weight1
    logical,intent(in)::DoAdd
    integer::iGet
       V1(Put%iCB_II(1,iPutStart),&
          Put%iCB_II(2,iPutStart),&
          Put%iCB_II(3,iPutStart),& 
          Put%iCB_II(4,iPutStart))=V1(Get%iCB_II(1,iGetStart),&
                                    Get%iCB_II(2,iGetStart),&
                                    Get%iCB_II(3,iGetStart),&
                                    Get%iCB_II(4,iGetStart))*&
                                    Weight%Weight_I(iGetStart)
       do iGet=iGetStart+1,iGetStart+nPartialGet-1
          V1(Put%iCB_II(1,iPutStart),&
             Put%iCB_II(2,iPutStart),&
             Put%iCB_II(3,iPutStart),& 
             Put%iCB_II(4,iPutStart))=&
             V1(Put%iCB_II(1,iPutStart),&
                Put%iCB_II(2,iPutStart),&
                Put%iCB_II(3,iPutStart),& 
                Put%iCB_II(4,iPutStart))+&
                V1(Get%iCB_II(1,iGet),&
                   Get%iCB_II(2,iGet),&
                   Get%iCB_II(3,iGet),&
                   Get%iCB_II(4,iGet))*&
                   Weight%Weight_I(iGet)
       end do
  end subroutine copyV1
!=============================================================
! For two dimensional grid
  subroutine test_global_message_pass_2dd(Grid)
    use ModNumConst
    type(DomainDecompositionType),intent(inout),target::Grid

    !Local variables
    integer :: i,j, iBlockStart, iBLK,nBLK
    integer::iError,lOctree, iImages
    integer,parameter :: nDim=2
    real, dimension(nDim),parameter :: F_D=(/1.0,0.13/)
    real::Xyz_D(nDim),Weight_I(4),VMisc
    integer::nImages,lSearch
    integer,dimension(nDim)::iCells_D,nCells_D
    integer,dimension(0:3,4)::Indexes_II
    type(RouterType)::Router
    type(GridDescriptorType)::Cell2,Cell,Temp 
    !------------------------------------------
    if(.not.is_proc())return
    call set_standard_grid_descriptor(Grid,GridDescriptor=Cell)
    call set_standard_grid_descriptor(Grid,nGhostGridPoints=nGhostCells,&
         GridDescriptor=Cell2)
    call init_router&
         (Cell,Cell2,Router)
    nCells_D=ncells_decomposition_d(Grid)
    Xyz_D=cHalf*(xyz_max_d(Grid)+xyz_min_d(Grid))
    if(is_proc0())write(*,*)"XyzTest:=",Xyz_D
    call set_standard_grid_descriptor(Grid,GridDescriptor=Temp)
    call nearest_grid_points(Xyz_D,Temp,3,Indexes_II,nImages,Weight_I)
      if(is_proc0())then
       do iImages=1,nImages
          write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
       end do
       write(*,*)'Nearest Cell, nImages=',nImages,'Control sum:',&
            sum(Weight_I(1:nImages))
    end if
    Xyz_D=cHalf*(xyz_max_d(Grid)+xyz_min_d(Grid))
    call interpolation_fix_reschange(nDim,Xyz_D,Temp,3,Indexes_II,nImages,Weight_I)
      if(is_proc0())then
       do iImages=1,nImages
          write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
       end do
       write(*,*)'Bilinear, Cells, nImages=',nImages,'Control sum:',&
            sum(Weight_I(1:nImages))
    end if
    Xyz_D=cHalf*(xyz_max_d(Grid)+xyz_min_d(Grid))
    if(is_proc0())write(*,*)"XyzTest",Xyz_D
    call clean_grid_descriptor(Temp)
    call set_standard_grid_descriptor(Grid,&
         Standard_=Nodes_,&
         GridDescriptor=Temp)
    call nearest_grid_points(Xyz_D,Temp,3,Indexes_II,nImages,Weight_I)
      if(is_proc0())then
       do iImages=1,nImages
          write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
       end do
       write(*,*)'NearestNode, nImages=',nImages,'Control sum:',&
            sum(Weight_I(1:nImages))
    end if
    Xyz_D=cHalf*(xyz_max_d(Grid)+xyz_min_d(Grid))
    call interpolation_fix_reschange(nDim,Xyz_D,Temp,3,Indexes_II,nImages,Weight_I)
      if(is_proc0())then
       do iImages=1,nImages
          write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
       end do
       write(*,*)'Bilinear,nodes, nImages=',nImages,'Control sum:',&
            sum(Weight_I(1:nImages))
    end if
    call clean_grid_descriptor(Temp)
    call set_router(&
         GridDescriptorSource=Cell,&
         GridDescriptorTarget=Cell2,&
         Router=Router,&
         interface_point_coords=do_all_ghostcells) 

    nBLK=max_block_pe(Grid,i_proc())
    iBlockStart=min_block_pe(Grid,i_proc())

    write(*,*)'For Grid of the component ',compid_grid(Grid)
    write(*,*)'Testing message_pass, PE=',i_proc(),'  Starting tests ...'
    allocate( V20(1-nGhostCells:&
         nCells_D(1)+nGhostCells,&
         1-nGhostCells:&
         nCells_D(2)+nGhostCells,&
         iBlockStart:nBLK), stat=iError )
    call check_allocate(iError,'V20')
    allocate( V21(1-nGhostCells:&
         nCells_D(1)+nGhostCells,&
         1-nGhostCells:&
         nCells_D(2)+nGhostCells,&
         iBlockStart:nBLK), stat=iError )
    call check_allocate(iError,'V21')

    V20=-999999.
    V21=V20
    do lOctree=1,ntree_nodes(Grid)
       if(.not.used_node(Grid,lOctree))CYCLE
       if(pe_decomposition(Grid,lOctree)/=i_proc())CYCLE
       iBLK=blk_decomposition(Grid,lOctree)
!V20 is assigned througout all cells
       do j=1-nGhostCells,nCells_D(2)+nGhostCells 
          iCells_D(2)=j
          do i=1-nGhostCells,nCells_D(1)+nGhostCells
             iCells_D(1)=i
             V20(i,j,iBLK)= &
                  dot_product(F_D,xyz_cell_d(Grid,lOctree,iCells_D))
          end do
       end do
!In sending messages in the ghostcells which cover coarser block
!message pass writes equal values to all the ghostcells which are
!covered by common coarser physical cell. To mimic this effect
!the data in the array V0 are coarsened too in such the ghostcells
       do j=1-nGhostCells,nCells_D(2)+nGhostCells,2 
          iCells_D(2)=j
          do i=1-nGhostCells,nCells_D(1)+nGhostCells,2
             iCells_D(1)=i
             if(l_level_neighbor(Grid,lOctree,iCells_D)==1)&
                  V20(i:i+1,j:j+1,iBLK)=&
                  cQuarter*sum(V20(i:i+1,j:j+1,iBLK))          
          end do
       end do
!V21 is assigned in true cells only
       do j=1,nCells_D(2)
          iCells_D(2)=j
          do i=1,nCells_D(1)
             iCells_D(1)=i
             V21(i,j,iBLK)= &
                  dot_product(F_D,xyz_cell_d(Grid,lOctree,iCells_D))
          end do
       end do
       if (any(is_left_boundary_d(Grid,lOctree))&
            .or.any(is_right_boundary_d(Grid,lOctree)))&
            V21(:,:,iBLK)=V20(:,:,iBLK)       
    end do

!Ghost cell values for V21 are filled in using global message pass
    call global_message_pass(Router,&
                             1,&      !
                             get2V1,&  !Wrapper for Source
                             put2V1, & !Wrapper for Target
                             copy2V1)  !Wrapper for Source and Target
!Compare V20 and V21
    if(max(maxval(V21-V20),maxval(V20-V21))>.01)then
       write(*,*)'Testing message_pass, PE=',i_proc(),&
            ' max difference=', &
            max(maxval(V21-V20),maxval(V20-V21)), &
            ' printing out values and exiting.'
       write(*,*)' '
       do lOctree=1,ntree_nodes(Grid)
          if(.not.used_node(Grid,lOctree))CYCLE
          if(pe_decomposition(Grid,lOctree)/=i_proc())CYCLE
          iBLK=blk_decomposition(Grid,lOctree)
          write(*,*)xyz_block_d(Grid,lOctree)
          do j=1-nGhostCells,nCells_D(2)+nGhostCells 
             do i=1-nGhostCells,nCells_D(1)+nGhostCells
                if(abs(V20(i,j,iBLK)-V21(i,j,iBLK))>.01)&
                     write(*,*)iBLK,i,j,V20(i,j,iBLK),V21(i,j,iBLK)
             end do
          end do
       end do
       write(*,*)' '
       call CON_stop('test_global_message_pass_failed')
    end if

    call global_barrier

    deallocate(V20)
    deallocate(V21)

    call clean_router(Router)
    call clean_grid_descriptor(Cell2)
    call clean_grid_descriptor(Cell)

    write(*,*)'Testing message_pass, PE=',i_proc(), &
         '  All tests passed.'
    write(*,*)' '

  end subroutine test_global_message_pass_2dd
  !========================================================
  subroutine test_global_message_pass_2id(GridID_)
    use ModNumConst
    integer,intent(in)::GridID_

    !Local variables
    integer :: i,j, iBlockStart, iBLK,nBLK
    integer::iError,lOctree, iImages
    integer,parameter :: nDim=2
    real, dimension(nDim),parameter :: F_D=(/1.0,0.13/)
    real::Xyz_D(nDim),Weight_I(4),VMisc
    integer::nImages,lSearch
    integer,dimension(nDim)::iCells_D,nCells_D
    integer,dimension(0:3,4)::Indexes_II
    type(RouterType)::Router
    type(GridDescriptorType)::Cell2,Cell,Temp 
    !------------------------------------------
    if(.not.is_proc())return
    call set_standard_grid_descriptor(GridID_,GridDescriptor=Cell)
    call set_standard_grid_descriptor(GridID_,nGhostGridPoints=nGhostCells,&
         GridDescriptor=Cell2)
    call init_router&
         (Cell,Cell2,Router)
    nCells_D=ncells_decomposition_d(GridID_)
    Xyz_D=cHalf*(xyz_max_d(GridID_)+xyz_min_d(GridID_))
    if(is_proc0())write(*,*)"XyzTest",Xyz_D
    call set_standard_grid_descriptor(GridID_,GridDescriptor=Temp)
    call nearest_grid_points(Xyz_D,Temp,3,Indexes_II,nImages,Weight_I)
      if(is_proc0())then
       do iImages=1,nImages
          write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
       end do
       write(*,*)'Nearest Cell, nImages=',nImages,'Control sum:',&
            sum(Weight_I(1:nImages))
    end if
    Xyz_D=cHalf*(xyz_max_d(GridID_)+xyz_min_d(GridID_))
    call interpolation_fix_reschange(nDim,Xyz_D,Temp,3,Indexes_II,nImages,Weight_I)
      if(is_proc0())then
       do iImages=1,nImages
          write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
       end do
       write(*,*)'Bilinear, Cells, nImages=',nImages,'Control sum:',&
            sum(Weight_I(1:nImages))
    end if
    Xyz_D=cHalf*(xyz_max_d(GridID_)+xyz_min_d(GridID_))
    if(is_proc0())write(*,*)"XyzTest",Xyz_D
    call clean_grid_descriptor(Temp)
    call set_standard_grid_descriptor(GridID_,&
         Standard_=Nodes_,&
         GridDescriptor=Temp)
    call nearest_grid_points(Xyz_D,Temp,3,Indexes_II,nImages,Weight_I)
      if(is_proc0())then
       do iImages=1,nImages
          write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
       end do
       write(*,*)'NearestNode, nImages=',nImages,'Control sum:',&
            sum(Weight_I(1:nImages))
    end if
    Xyz_D=cHalf*(xyz_max_d(GridID_)+xyz_min_d(GridID_))
    call interpolation_fix_reschange(nDim,Xyz_D,Temp,3,Indexes_II,nImages,Weight_I)
      if(is_proc0())then
       do iImages=1,nImages
          write(*,*) Indexes_II(:,iImages),Weight_I(iImages)
       end do
       write(*,*)'Bilinear,nodes, nImages=',nImages,'Control sum:',&
            sum(Weight_I(1:nImages))
    end if
    call clean_grid_descriptor(Temp)
    call set_router(&
         GridDescriptorSource=Cell,&
         GridDescriptorTarget=Cell2,&
         Router=Router,&
         interface_point_coords=do_all_ghostcells) 

    nBLK=max_block_pe(GridID_,i_proc())
    iBlockStart=min_block_pe(GridID_,i_proc())

    write(*,*)'For Grid of the component ',compid_grid(GridID_)
    write(*,*)'Testing message_pass, PE=',i_proc(),'  Starting tests ...'
    allocate( V20(1-nGhostCells:&
         nCells_D(1)+nGhostCells,&
         1-nGhostCells:&
         nCells_D(2)+nGhostCells,&
         iBlockStart:nBLK), stat=iError )
    call check_allocate(iError,'V20')
    allocate( V21(1-nGhostCells:&
         nCells_D(1)+nGhostCells,&
         1-nGhostCells:&
         nCells_D(2)+nGhostCells,&
         iBlockStart:nBLK), stat=iError )
    call check_allocate(iError,'V21')

    V20=-999999.
    V21=V20
    do lOctree=1,ntree_nodes(GridID_)
       if(.not.used_node(GridID_,lOctree))CYCLE
       if(pe_decomposition(GridID_,lOctree)/=i_proc())CYCLE
       iBLK=blk_decomposition(GridID_,lOctree)
!V20 is assigned througout all cells
       do j=1-nGhostCells,nCells_D(2)+nGhostCells 
          iCells_D(2)=j
          do i=1-nGhostCells,nCells_D(1)+nGhostCells
             iCells_D(1)=i
             V20(i,j,iBLK)= &
                  dot_product(F_D,xyz_cell_d(GridID_,lOctree,iCells_D))
          end do
       end do
!In sending messages in the ghostcells which cover coarser block
!message pass writes equal values to all the ghostcells which are
!covered by common coarser physical cell. To mimic this effect
!the data in the array V0 are coarsened too in such the ghostcells
       do j=1-nGhostCells,nCells_D(2)+nGhostCells,2 
          iCells_D(2)=j
          do i=1-nGhostCells,nCells_D(1)+nGhostCells,2
             iCells_D(1)=i
             if(l_level_neighbor(GridID_,lOctree,iCells_D)==1)&
                  V20(i:i+1,j:j+1,iBLK)=&
                  cQuarter*sum(V20(i:i+1,j:j+1,iBLK))          
          end do
       end do
!V21 is assigned in true cells only
       do j=1,nCells_D(2)
          iCells_D(2)=j
          do i=1,nCells_D(1)
             iCells_D(1)=i
             V21(i,j,iBLK)= &
                  dot_product(F_D,xyz_cell_d(GridID_,lOctree,iCells_D))
          end do
       end do
       if (any(is_left_boundary_d(GridID_,lOctree))&
            .or.any(is_right_boundary_d(GridID_,lOctree)))&
            V21(:,:,iBLK)=V20(:,:,iBLK)       
    end do

!Ghost cell values for V21 are filled in using global message pass
    call global_message_pass(Router,&
                             1,&      !
                             get2V1,&  !Wrapper for Source
                             put2V1, & !Wrapper for Target
                             copy2V1)  !Wrapper for Source and Target
!Compare V20 and V21
    if(max(maxval(V21-V20),maxval(V20-V21))>.01)then
       write(*,*)'Testing message_pass, PE=',i_proc(),&
            ' max difference=', &
            max(maxval(V21-V20),maxval(V20-V21)), &
            ' printing out values and exiting.'
       write(*,*)' '
       do lOctree=1,ntree_nodes(GridID_)
          if(.not.used_node(GridID_,lOctree))CYCLE
          if(pe_decomposition(GridID_,lOctree)/=i_proc())CYCLE
          iBLK=blk_decomposition(GridID_,lOctree)
          write(*,*)xyz_block_d(GridID_,lOctree)
          do j=1-nGhostCells,nCells_D(2)+nGhostCells 
             do i=1-nGhostCells,nCells_D(1)+nGhostCells
                if(abs(V20(i,j,iBLK)-V21(i,j,iBLK))>.01)&
                     write(*,*)iBLK,i,j,V20(i,j,iBLK),V21(i,j,iBLK)
             end do
          end do
       end do
       write(*,*)' '
       call CON_stop('test_global_message_pass_failed')
    end if

    call global_barrier

    deallocate(V20)
    deallocate(V21)

    call clean_router(Router)
    call clean_grid_descriptor(Cell2)
    call clean_grid_descriptor(Cell)

    write(*,*)'Testing message_pass, PE=',i_proc(), &
         '  All tests passed.'
    write(*,*)' '

  end subroutine test_global_message_pass_2id
!===============================================================
!===============================================================
!Subroutines used for testing global message pass
  subroutine get2V1(nPartial,&
                  iGetStart,&
                  Get,&
                  Weight,&
                  Buff_I,nVar)
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get    
    type(WeightPtrType),intent(in)::Weight
    real,dimension(nVar),intent(out)::Buff_I
    integer::iGet
    Buff_I(1)=V21(Get%iCB_II(1,iGetStart),&
                 Get%iCB_II(2,iGetStart),&
                 Get%iCB_II(3,iGetStart))*&
                 Weight%Weight_I(iGetStart)
    do iGet=iGetStart+1,iGetStart+nPartial-1
       Buff_I(1)= Buff_I(1)+&
                  V21(Get%iCB_II(1,iGet),&
                     Get%iCB_II(2,iGet),&
                     Get%iCB_II(3,iGet))*&
                     Weight%Weight_I(iGet)
    end do
  end subroutine get2V1
!==============================================================
  subroutine put2V1(nPartial,&
                   iPut,&
                   Put,&
                   Weight,&
                   DoAdd,Buff_I,nVar)
    integer,intent(in)::iPut,nPartial,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::Weight
    real,dimension(nVar),intent(in)::Buff_I
    logical,intent(in)::DoAdd
    if(DoAdd)then
       V21(Put%iCB_II(1,iPut),&
          Put%iCB_II(2,iPut),&
          Put%iCB_II(3,iPut))=&
          V21(Put%iCB_II(1,iPut),&
             Put%iCB_II(2,iPut),&
             Put%iCB_II(3,iPut))+Buff_I(1)
    else
       V21(Put%iCB_II(1,iPut),&
          Put%iCB_II(2,iPut),&
          Put%iCB_II(3,iPut))=Buff_I(1)
    end if
  end subroutine put2V1
!==============================================================
  subroutine copy2V1(nPartialGet,&
                    iGetStart,  &
                    Get,&
                    Weight,&
                    nPartialPut,&
                    iPutStart,&
                    Put,Weight1,DoAdd)
    integer,intent(in)::nPartialGet,iGetStart
    integer,intent(in)::nPartialPut,iPutStart
    type(IndexPtrType),intent(in)::Get,Put
    type(WeightPtrType),intent(in)::Weight,Weight1
    logical,intent(in)::DoAdd
    integer::iGet
       V21(Put%iCB_II(1,iPutStart),&
          Put%iCB_II(2,iPutStart),&
          Put%iCB_II(3,iPutStart))=V21(Get%iCB_II(1,iGetStart),&
                                    Get%iCB_II(2,iGetStart),&
                                    Get%iCB_II(3,iGetStart))*&
                                    Weight%Weight_I(iGetStart)
       do iGet=iGetStart+1,iGetStart+nPartialGet-1
          V21(Put%iCB_II(1,iPutStart),&
             Put%iCB_II(2,iPutStart),&
             Put%iCB_II(3,iPutStart))=&
             V21(Put%iCB_II(1,iPutStart),&
                Put%iCB_II(2,iPutStart),&
                Put%iCB_II(3,iPutStart))+&
                V21(Get%iCB_II(1,iGet),&
                   Get%iCB_II(2,iGet),&
                   Get%iCB_II(3,iGet))*&
                   Weight%Weight_I(iGet)
       end do
  end subroutine copy2V1
!==================================================================
  subroutine do_all_ghostcells(&
       GridDescriptor,&
       lGlobalTreeNode,&
       nDim,&
       Xyz_D,&
       nIndexes,&
       Index_I,&
       IsInterfacePoint)
    type(GridDescriptorType),intent(in):: GridDescriptor
    integer,intent(in)::lGlobalTreeNode,nIndexes
    logical,intent(out)::IsInterfacePoint
    integer,intent(in)::nDim
    real,dimension(nDim),intent(inout)::Xyz_D
    integer, dimension(nIndexes),intent(inout)::Index_I
    logical,dimension(nDim)::&
         IsLeftFace_D,IsRightFace_D
    !----------------------------------------------------------------------

    IsLeftFace_D=Index_I(1:nDim)<1
    IsRightFace_D=Index_I(1:nDim)>&
         ncells_decomposition_d(GridDescriptor%DD%Ptr)
    IsInterfacePoint=&
         any(IsLeftFace_D.or.IsRightFace_D).and.&
         (.not.any(IsLeftFace_D.and.is_left_boundary_d(&
         GridDescriptor%DD%Ptr,lGlobalTreeNode))).and.&  
         (.not.any(IsRightFace_D.and.is_right_boundary_d(&
         GridDescriptor%DD%Ptr,lGlobalTreeNode)))
  end subroutine do_all_ghostcells
end Module CON_test_global_message_pass
