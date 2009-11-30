module ModSetOuterBC

  ! Notation: 1g - first ghost cell,    2g - second ghost cell
  !           1p - first physical cell, 2p - second physical cell

  integer :: imin1g,imax1g,imin2g,imax2g,imin1p,imax1p,imin2p,imax2p
  integer :: jmin1g,jmax1g,jmin2g,jmax2g,jmin1p,jmax1p,jmin2p,jmax2p
  integer :: kmin1g,kmax1g,kmin2g,kmax2g,kmin1p,kmax1p,kmin2p,kmax2p
  integer :: iBLK

end module ModSetOuterBC

!=============================================================================

module ModFaceBc

  ! Variables shared by subroutines set_BCs.f90 and user_face_bcs in ModUser
  
  use ModVarIndexes, ONLY: nVar

  implicit none

  SAVE

  ! True if only boundaries at resolution changes are updated
  logical:: DoResChangeOnly

  ! The type and index of the boundary
  character(len=20) :: TypeBc

  ! Negative iBoundary indicates which body we are computing for.
  ! Positive iBoundary numerates the sides for the outer BCs.
  integer:: iBoundary

  ! Index of the face
  integer:: iFace, jFace, kFace, iBlockBc

  ! The side of the cell defined with respect to the cell inside the domain
  integer :: iSide

  ! The values on the physical side and the ghost cell side of the boudary
  real:: VarsTrueFace_V(nVar), VarsGhostFace_V(nVar)

  ! The coordinates of the face center and the B0 field at that point
  real :: FaceCoords_D(3), B0Face_D(3)

  ! The time at which the (time dependent) boundary condition is calculated
  real :: TimeBc

end module ModFaceBc

!============================================================================

module ModBoundaryCells

  implicit none
  SAVE

  logical::SaveBoundaryCells=.false.
  logical::ResetBody2Cells=.false.
  integer::MinBoundarySaved=-777,MaxBoundarySaved=-777,nBoundarySaved=0
  logical,allocatable,dimension(:,:,:,:,:)::IsBoundaryCell_IGB

contains

  !===========================================================================
  subroutine allocate_boundary_cells

    use ModSize
    use ModGeometry, ONLY: MaxBoundary, ExtraBC_
    use ModMain,     ONLY: UseExtraBoundary
    use ModProcMH
    !-------------------------------------------------------------------------
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
    if(MaxBoundarySaved < MinBoundarySaved)then
       if(iProc==0)write(*,*)'Only Extra and Outer boundary cells can be saved'
       SaveBoundaryCells=.false.
    end if
    nBoundarySaved = MaxBoundarySaved - MinBoundarySaved + 1
    allocate(IsBoundaryCell_IGB(MinBoundarySaved:MaxBoundarySaved,&
         1-gcn:nI+gcn, 1-gcn:nJ+gcn,  1-gcn:nK+gcn, nBLK))

  end subroutine allocate_boundary_cells

end module ModBoundaryCells

!=============================================================================

subroutine fix_boundary_ghost_cells(UseMonotoneRestrict)

  use ModBoundaryCells, ONLY: MinBoundarySaved, MaxBoundarySaved, &
       IsBoundaryCell_IGB
  use ModMain, ONLY : nBlock, UnusedBlk, iNewGrid, iNewDecomposition
  use ModGeometry, ONLY: true_cell, body_BLK, IsBoundaryBlock_IB
  !use ModProcMH, ONLY: iProc

  implicit none

  logical, intent(in):: UseMonotoneRestrict

  integer:: iBlock, iBoundary
  integer:: iGridHere=-1, iDecompositionHere=-1
!-----------------------------------------------------------------------------
  if(iGridHere==iNewGrid .and. iDecompositionHere==iNewDecomposition) RETURN

!  if(iProc==0)write(*,*)'Start fix boundary cells'

  iGridHere=iNewGrid; iDecompositionHere=iNewDecomposition
  call message_pass_boundary_cells(UseMonotoneRestrict)

  do iBlock = 1, nBlock

     if(unusedBLK(iBlock)) CYCLE

     do iBoundary = MinBoundarySaved, MaxBoundarySaved
        IsBoundaryBlock_IB(iBoundary,iBlock)=&
             any(IsBoundaryCell_IGB(iBoundary,:,:,:,iBlock))
        true_cell(:,:,:,iBlock) = true_cell(:,:,:,iBlock) &
             .and. .not.IsBoundaryCell_IGB(iBoundary,:,:,:,iBlock)
     end do

     ! body_BLK: if any cell INCLUDING ghost cells is outside the comp. domain
     body_BLK(iBlock) = .not. all(true_cell(:,:,:,iBlock))

  end do

end subroutine fix_boundary_ghost_cells
