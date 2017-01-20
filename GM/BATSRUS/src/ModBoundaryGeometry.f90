!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!============================================================================

module ModBoundaryGeometry

  implicit none
  SAVE

  ! iBoundary_GB contains the index of the boundary that the cell belongs to.
  integer, allocatable :: iBoundary_GB(:,:,:,:)

  ! Cells inside domain have index domain_ that is smaller than smallest 
  ! boundary index (body2_ = -2)
  integer, parameter :: domain_ = -10

contains

  !============================================================================
  subroutine init_mod_boundary_cells

    use ModSize, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, MaxBlock
    !--------------------------------------------------------------------------

    if(.not. allocated(iBoundary_GB)) then
       allocate(iBoundary_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
       iBoundary_GB = domain_
    end if

  end subroutine init_mod_boundary_cells

  !============================================================================

  subroutine fix_block_geometry(iBlock, DoSolveSolidIn)

    use ModMain, ONLY: body1_, body2_, ExtraBc_, Solid_, &
         UseBody2, UseExtraBoundary, UseSolidState, &
         ProcTest, BlkTest, iTest, jTest, kTest
    use ModGeometry, ONLY: &
         MinFaceBoundary, MaxFaceBoundary, IsBoundaryBlock_IB, &
         IsBoundaryCell_GI, R2_BLK, Rmin2_BLK, Body_BLK, &
         far_field_BCs_BLK, true_blk, true_cell
    use ModPhysics, ONLY : xBody2,yBody2,zBody2
    use BATL_lib, ONLY: &
         iProc, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nI, nJ, nK, nG, &
         Xyz_DGB, CellSize_DB, &
         CoordMin_DB, CoordMax_DB, CoordMin_D, CoordMax_D, IsPeriodic_D

    implicit none

    integer, intent(in) :: iBlock
    logical, intent(in), optional :: DoSolveSolidIn

    integer :: i, j, k, iBoundary
    logical :: DoSolveSolid
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='fix_block_geometry'
    !--------------------------------------------------------------------------

    if(iBlock==BlkTest .and. iProc==ProcTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    DoSolveSolid = .true.
    if(present(DoSolveSolidIn)) DoSolveSolid = DoSolveSolidIn

    if (UseBody2) then
       ! calculate the radius as measured from the second body
       ! Note that the second body can move
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          R2_BLK(i,j,k,iBlock) = sqrt( &
               sum( (Xyz_DGB(:,i,j,k,iBlock)- (/xBody2, yBody2, zBody2/))**2 ))
       end do; end do; end do
       Rmin2_BLK(iBlock) = minval(R2_BLK(:,:,:,iBlock))
    else
       Rmin2_BLK(iBlock) = 0.0
    end if

    far_field_BCs_BLK(iBlock) = any( .not. IsPeriodic_D .and.   &
         (    CoordMin_DB(:,iBlock) <= CoordMin_D + 1e-12       &
         .or. CoordMax_DB(:,iBlock) >= CoordMax_D - 1e-12) )

    if(DoTestMe)then
       write(*,*)NameSub,': far_field_bcs_blk=',far_field_bcs_BLK(iBlock)
       write(*,*)NameSub,': CoordMin_DB=', CoordMin_DB(:,iBlock)
       write(*,*)NameSub,': CellSize_D = ',CellSize_DB(:,iBlock)
       write(*,*)NameSub,': CoordMin_D =', CoordMin_D
       write(*,*)NameSub,': CoordMax_D =', CoordMax_D
       write(*,*)NameSub,': MinFaceBoundary, MaxFaceBoundary=',  &
            MinFaceBoundary, MaxFaceBoundary
    end if

    !\
    ! TRUE_CELL: if not inside a body or outside the outer face boundary
    !/
    true_cell(:,:,:,iBlock)=.true.
    IsBoundaryBlock_IB(:,iBlock)=.false.
    do iBoundary = MinFaceBoundary, MaxFaceBoundary
       IsBoundaryBlock_IB(iBoundary,iBlock) = .true.
    end do
    IsBoundaryBlock_IB(ExtraBc_,iBlock) = UseExtraBoundary
    if(DoSolveSolid)then
       IsBoundaryBlock_IB(Solid_,iBlock) = .false.
    else
       IsBoundaryBlock_IB(Solid_,iBlock) = UseSolidState
    end if

    ! set IsBoundaryCell_GI to be used for true_cell
    call set_boundary_cells(iBlock)

    do iBoundary = MinFaceBoundary, max(MaxFaceBoundary,Body1_)
       IsBoundaryBlock_IB(iBoundary,iBlock) = &
            any(IsBoundaryCell_GI(:,:,:,iBoundary))
       true_cell(:,:,:,iBlock) = &
            true_cell(:,:,:,iBlock).and..not.IsBoundaryCell_GI(:,:,:,iBoundary)
    end do

    if(DoTestMe)write(*,*) NameSub,&
         ' after call set_boundary_cells: true_cell(iTest-nG:iTest+nG)=',&
         true_cell(iTest-nG:iTest+nG,jTest,kTest,iBlock)

    IsBoundaryCell_GI(:,:,:,ExtraBc_) = &
         UseExtraBoundary .and. IsBoundaryCell_GI(:,:,:,ExtraBc_)

    if(DoSolveSolid)then
       IsBoundaryCell_GI(:,:,:,Solid_) = .false.
    else
       IsBoundaryCell_GI(:,:,:,Solid_) = &
            UseSolidState .and. IsBoundaryCell_GI(:,:,:,Solid_)
    end if

    ! Copying  the IsBoundaryCell_GI into the format for iBoundary_GB
    iBoundary_GB(:,:,:,iBlock) = domain_
    do iBoundary = Solid_, 6
       where(IsBoundaryCell_GI(:,:,:,iBoundary))
          iBoundary_GB(:,:,:,iBlock) = iBoundary
       end where
    end do

    ! Alow other places to set true_cell
    true_cell(1:nI,1:nJ,1:nK,iBlock) = true_cell(1:nI,1:nJ,1:nK,iBlock) &
         .and. iBoundary_GB(1:nI,1:nJ,1:nK,iBlock) == domain_

    ! body_BLK: if any cell INCLUDING ghost cells is inside body(ies)
    body_BLK(iBlock) = .not. all(true_cell(:,:,:,iBlock))

    ! TRUE_BLK: if all cells EXCLUDING ghost cells are outside body(ies)
    true_BLK(iBlock) = all(true_cell(1:nI,1:nJ,1:nK,iBlock))

    if(DoTestMe)then
       write(*,*) NameSub,&
            ' finished with iBoundary_GB(iTest-nG:iTest+nG)=',&
            iBoundary_GB(iTest-nG:iTest+nG,jTest,kTest,iBlock)

       write(*,*) NameSub,&
            ' finished with true_cell(iTest-nG:iTest+nG)=',&
            true_cell(iTest-nG:iTest+nG,jTest,kTest,iBlock)
    end if

  end subroutine fix_block_geometry

  !============================================================================

  subroutine fix_geometry(DoSolveSolidIn)

    use BATL_lib, ONLY: nBlock

    logical, intent(in) :: DoSolveSolidIn

    integer :: iBlock
    !--------------------------------------------------------------------------

    do iBlock = 1, nBlock
       call fix_block_geometry(iBlock, DoSolveSolidIn)
    end do

  end subroutine fix_geometry

  !============================================================================

  subroutine set_boundary_cells(iBlock)

    use ModMain, ONLY: iTest, jTest, kTest, BlkTest, ProcTest, &
         Body1, Body1_, UseBody2, Body2_, ExtraBc_, Solid_
    use ModPhysics,  ONLY: Rbody
    use ModGeometry, ONLY: R_BLK, IsBoundaryBlock_IB, IsBoundaryCell_GI
    use ModPhysics,  ONLY: Rbody2
    use ModGeometry, ONLY: R2_BLK
    use ModGeometry, ONLY: x1, x2, y1, y2, z1, z2
    use BATL_lib,    ONLY: iProc, Xyz_DGB, x_, y_, z_, nG
    use ModUserInterface ! user_set_boundary_cells

    implicit none
    integer, intent(in)::iBlock

    logical:: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='set_boundary_cells'
    !--------------------------------------------------------------------------
    if(iBlock==BlkTest .and. iProc==ProcTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    if(DoTestMe)write(*,*) NameSub,' IsBoundaryBlock_IB=', &
         IsBoundaryBlock_IB(:,iBlock)

    IsBoundaryCell_GI=.false.  

    if(UseBody2 .and. IsBoundaryBlock_IB(Body2_,iBlock)) &
         IsBoundaryCell_GI(:,:,:,Body2_) = R2_BLK(:,:,:,iBlock) < rBody2  

    if(IsBoundaryBlock_IB(Body1_,iBlock)) &
         IsBoundaryCell_GI(:,:,:,Body1_) = &
         body1    .and. R_BLK(:,:,:,iBlock) < Rbody

    if(IsBoundaryBlock_IB(ExtraBc_,iBlock).or.IsBoundaryBlock_IB(Solid_,iBlock))&
         call user_set_boundary_cells(iBlock)

    if(IsBoundaryBlock_IB(1,iBlock)) &
         IsBoundaryCell_GI(:,:,:,1)= Xyz_DGB(x_,:,:,:,iBlock) < x1  
    if(IsBoundaryBlock_IB(2,iBlock)) &
         IsBoundaryCell_GI(:,:,:,2)= Xyz_DGB(x_,:,:,:,iBlock) > x2 
    if(IsBoundaryBlock_IB(3,iBlock)) &
         IsBoundaryCell_GI(:,:,:,3)= Xyz_DGB(y_,:,:,:,iBlock) < y1
    if(IsBoundaryBlock_IB(4,iBlock)) &
         IsBoundaryCell_GI(:,:,:,4)= Xyz_DGB(y_,:,:,:,iBlock) > y2 
    if(IsBoundaryBlock_IB(5,iBlock)) &
         IsBoundaryCell_GI(:,:,:,5)= Xyz_DGB(z_,:,:,:,iBlock) < z1
    if(IsBoundaryBlock_IB(6,iBlock)) &
         IsBoundaryCell_GI(:,:,:,6)= Xyz_DGB(z_,:,:,:,iBlock) > z2  

  end subroutine set_boundary_cells

end module ModBoundaryGeometry

!=============================================================================

subroutine fix_boundary_ghost_cells

  ! Recalculate true_cell information in ghost cells if grid changed. 

  use ModBoundaryGeometry, ONLY: iBoundary_GB, domain_
  use ModMain, ONLY : nBlock, Unused_B, iNewGrid, iNewDecomposition, &
       body2_, BlkTest, iTest, jTest, kTest, iteration_number, Solid_
  use ModGeometry, ONLY: true_cell, body_BLK, IsBoundaryBlock_IB
  !use ModProcMH, ONLY: iProc
  use BATL_lib, ONLY: message_pass_cell

  implicit none

  integer:: iBlock, iBoundary
  integer:: iGridHere = -1, iDecompositionHere = -1

  logical:: DoTest, DoTestMe
  character(len=*), parameter:: NameSub = 'fix_boundary_ghost_cells'
  !----------------------------------------------------------------------------
  if(iGridHere==iNewGrid .and. iDecompositionHere==iNewDecomposition) RETURN

  iGridHere = iNewGrid; iDecompositionHere = iNewDecomposition

  call set_oktest(NameSub, DoTest, DoTestMe)

  if(DoTestMe) write(*,*)NameSub,' starting with true_cell(i-2:i+2)=', &
       true_cell(iTest-2:iTest+2,jTest,kTest,BlkTest)

  ! DoResChangeOnly=true works as long as the ghost cells are correctly set
  ! away from resolution changes. This usually holds, but not if the 
  ! boundary cells are set based on the state variables read from a restart
  ! file that has no ghost cell information saved. This can only happen
  ! at the very beginning of a run when iteration_number == 0.

  call message_pass_cell(iBoundary_GB, &
       nProlongOrderIn=1, nCoarseLayerIn=2, &
       DoSendCornerIn=.true., DoRestrictFaceIn=.true., &
       DoResChangeOnlyIn=iteration_number>0, NameOperatorIn='max')

  if(DoTestMe) write(*,*) NameSub,': iBoundary_GB(i-2:i+2)=', &
       iBoundary_GB(iTest-2:iTest+2,jTest,kTest,BlkTest)

  do iBlock = 1, nBlock
     if(Unused_B(iBlock)) CYCLE

     true_cell(:,:,:,iBlock) = true_cell(:,:,:,iBlock)  &
          .and. (iBoundary_GB(:,:,:,iBlock) == domain_)

     body_BLK(iBlock) = .not. all(true_cell(:,:,:,iBlock))   

     do iBoundary = Solid_, 6
        IsBoundaryBlock_IB(iBoundary,iBlock)= &
             any(iBoundary_GB(:,:,:,iBlock) == iBoundary)
     end do
  end do

  if(DoTestMe) write(*,*) NameSub,' finished with true_cell(i-2:i+2)=', &
       true_cell(iTest-2:iTest+2,jTest,kTest,BlkTest)

end subroutine fix_boundary_ghost_cells
