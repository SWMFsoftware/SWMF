module ModBoundaryCells
  implicit none
  SAVE
  logical::SaveBoundaryCells=.false.
  integer::MinBoundarySaved=-777,MaxBoundarySaved=-777,nBoundarySaved=0
  logical,allocatable,dimension(:,:,:,:,:)::IsBoundaryCell_IGB
contains
  subroutine allocate_boundary_cells
    use ModSize
    use ModGeometry,ONLY:MinBoundary,MaxBoundary,ExtraBC_
    use ModMain,ONLY:UseExtraBoundary
    use ModProcMH
    if(.not.SaveBoundaryCells)then
       if(iProc==0)write(*,*)&
            'Do not allocate boundary cell array, if SaveBoundaryCells=',&
            SaveBoundaryCells
       return
    end if
    if(allocated(IsBoundaryCell_IGB))return
   
    if(UseExtraBoundary)then
       MinBoundarySaved=ExtraBC_
    else
       MinBoundarySaved=ExtraBC_+1
    end if
    MaxBoundarySaved=MaxBoundary
    if(MaxBoundarySaved<MinBoundarySaved)then
       if(iProc==0)write(*,*)'Only Extra and Outer boundary cells can be saved'
       SaveBoundaryCells=.false.
    end if
    nBoundarySaved=MaxBoundarySaved-MinBoundarySaved+1
    allocate(IsBoundaryCell_IGB(MinBoundarySaved:MaxBoundarySaved,&
         1-gcn:nI+gcn, 1-gcn:nJ+gcn,  1-gcn:nK+gcn, nBLK))
  end subroutine allocate_boundary_cells
end module ModBoundaryCells
!=============================================================================
subroutine fix_boundary_ghost_cells(UseMonotoneRestrict)
  use ModBoundaryCells
  use ModMain, ONLY : unusedBLK 
  use ModMain, ONLY : iNewGrid, iNewDecomposition
  use ModGeometry
  use ModProcMH
  implicit none
  !Subroutine arguements
  logical :: UseMonotoneRestrict
  integer::iBlock,iBoundary
  integer,save::iGridHere=-1,iDecompositionHere=-1
!-----------------------------------------------------------------------------
  if(iGridHere==iNewGrid.and.iDecompositionHere==iNewDecomposition)return
!  if(iProc==0)write(*,*)'Start fix boundary cells'
  iGridHere=iNewGrid; iDecompositionHere=iNewDecomposition
  call message_pass_boundary_cells(UseMonotoneRestrict)
  do iBlock=1,nBLK
     if(unusedBLK(iBlock))CYCLE
     do iBoundary=MinBoundarySaved,MaxBoundarySaved
        IsBoundaryBlock_IB(iBoundary,iBlock)=&
             any(IsBoundaryCell_IGB(iBoundary,:,:,:,iBlock))
        true_cell(:,:,:,iBlock) = true_cell(:,:,:,iBlock) &
             .and. .not.IsBoundaryCell_IGB(iBoundary,:,:,:,iBlock)
     end do
     !\
     ! body_BLK: if any cell INCLUDING ghost cells is inside body(ies)
     !/
     body_BLK(iBlock) = .not. all(true_cell(:,:,:,iBlock))
  end do
end subroutine fix_boundary_ghost_cells
