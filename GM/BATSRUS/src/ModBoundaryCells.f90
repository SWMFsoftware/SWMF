!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!============================================================================

module ModBoundaryCells

  implicit none
  SAVE

  ! iBoundary_GB contains the index of the boundary that the cell belongs to.
  integer, allocatable :: iBoundary_GB(:,:,:,:)

  ! Cells inside domain have index domain_ that is smaller than smallest 
  ! boundary index (body2_ = -2)
  integer, parameter :: domain_ = -10

  ! The DomainOp sets the operator applied for boundary indexes
  ! in ghost cells that are restricted. So the largest index is kept.
  character(len=*), parameter :: DomainOp='max'

contains

  !===========================================================================
  subroutine init_mod_boundary_cells

    use ModSize, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, MaxBlock
    !-------------------------------------------------------------------------

    if(.not. allocated(iBoundary_GB)) then
       allocate(iBoundary_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
       iBoundary_GB = domain_
    end if

  end subroutine init_mod_boundary_cells

end module ModBoundaryCells

!=============================================================================

subroutine fix_boundary_ghost_cells

  ! Recalculate true_cell information in ghost cells if grid changed. 

  use ModBoundaryCells, ONLY: iBoundary_GB, DomainOp, domain_
  use ModMain, ONLY : nBlock, Unused_B, iNewGrid, iNewDecomposition, &
       body2_, BlkTest, iTest, jTest, kTest
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

  call message_pass_cell(iBoundary_GB, &
       nProlongOrderIn=1, nCoarseLayerIn=2, &
       DoSendCornerIn=.true., DoRestrictFaceIn=.true., &
       NameOperatorIn=DomainOp)

  if(DoTestMe) write(*,*) NameSub,': iBoundary_GB(i-2:i+2)=', &
       iBoundary_GB(iTest-2:iTest+2,jTest,kTest,BlkTest)

  do iBlock = 1, nBlock
     if(Unused_B(iBlock)) CYCLE

     true_cell(:,:,:,iBlock) = true_cell(:,:,:,iBlock)  &
          .and. (iBoundary_GB(:,:,:,iBlock) == domain_)

     body_BLK(iBlock) = .not. all(true_cell(:,:,:,iBlock))   

     do iBoundary = body2_, 6
        IsBoundaryBlock_IB(iBoundary,iBlock)= &
             any(iBoundary_GB(:,:,:,iBlock) == iBoundary)
     end do
  end do

  if(DoTestMe) write(*,*) NameSub,' finished with true_cell(i-2:i+2)=', &
       true_cell(iTest-2:iTest+2,jTest,kTest,BlkTest)

end subroutine fix_boundary_ghost_cells
