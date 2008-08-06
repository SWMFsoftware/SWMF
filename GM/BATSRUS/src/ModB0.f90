module ModB0

  use ModSize, ONLY: nDim, nI, nJ, nK, MaxBlock, &
       East_, West_, South_, North_, Bot_, Top_

  use ModParallel, ONLY: &
       NeiLEast, NeiLWest, NeiLSouth, NeiLNorth, NeiLBot, NeiLTop

  ! B0 is the analytic part of the full magnetic field: FullB = B0 + B1
  ! B0 is divergence free, and usually curl free (unless UseCurlB0 is true). 
  ! B0 may be time dependent.
  ! The cell centered B0 is stored for all blocks and all (ghost) cells.
  ! The face centered B0 is calculated as the average of the cell centered 
  ! values except at resolution changes, where a more complicated approximation
  ! is used.

  implicit none

  private ! except

  ! Variables for all blocks
  real, public, allocatable, dimension(:,:,:,:,:):: &
       B0_DGB, B0ResChange_DXSB, B0ResChange_DYSB, B0ResChange_DZSB

  ! Variables for a particular block
  real, public, allocatable, dimension(:,:,:,:) :: B0_DX, B0_DY, B0_DZ

  public:: init_mod_b0, set_b0_face

contains

  !============================================================================

  subroutine init_mod_b0

    if(allocated(B0_DGB)) RETURN

    allocate( &
         B0_DGB(nDim, -1:nI+2, -1:nJ+2, -1:nK+2, MaxBlock), &
         B0ResChange_DXSB(nDim, nJ, nK, East_ :West_ , MaxBlock), &
         B0ResChange_DYSB(nDim, nI, nK, South_:North_, MaxBlock), &
         B0ResChange_DZSB(nDim, nI, nJ, Bot_  :Top_  , MaxBlock), &
         B0_DX(nDim, nI+1, nJ, nK),                               &
         B0_DY(nDim, nI, nJ+1, nK),                               &
         B0_DZ(nDim, nI, nJ, nK+1) )

  end subroutine init_mod_b0

  !============================================================================

  subroutine set_b0_face(iBlock)

    ! Set the face centered magnetic field for block iBlock
    integer, intent(in) :: iBlock
    !-------------------------------------------------------------------------

    ! Average 
    B0_DX = 0.5*(B0_DGB(:,0:nI  ,1:nJ,1:nK,iBlock) &
         +       B0_DGB(:,1:nI+1,1:nJ,1:nK,iBlock))
    B0_DY = 0.5*(B0_DGB(:,1:nI,0:nJ  ,1:nK,iBlock) &
         +       B0_DGB(:,1:nI,1:nJ+1,1:nK,iBlock))
    B0_DZ = 0.5*(B0_DGB(:,1:nI,1:nJ,0:nK,  iBlock) &
         +       B0_DGB(:,1:nI,1:nJ,1:nK+1,iBlock))

    ! Resolution changes
    if(NeiLEast(iBlock) == -1) &
         B0_DX(:,   1,1:nJ,1:nK) = B0ResChange_DXSB(:,:,:,East_ ,iBlock)
    if(NeiLWest(iBlock) == -1) &
         B0_DX(:,nI+1,1:nJ,1:nK) = B0ResChange_DXSB(:,:,:,West_ ,iBlock)
    if(NeiLSouth(iBlock) == -1) &
         B0_DY(:,1:nI,   1,1:nK) = B0ResChange_DYSB(:,:,:,South_,iBlock)
    if(NeiLNorth(iBlock) == -1) &
         B0_DY(:,1:nI,nJ+1,1:nK) = B0ResChange_DYSB(:,:,:,North_,iBlock)
    if(NeiLBot(iBlock) == -1) &
         B0_DZ(:,1:nI,1:nJ,   1) = B0ResChange_DZSB(:,:,:,Bot_,  iBlock)
    if(NeiLTop(iBlock) == -1) &
         B0_DZ(:,1:nI,1:nJ,nK+1) = B0ResChange_DZSB(:,:,:,Top_,  iBlock)

  end subroutine set_b0_face

end module ModB0
