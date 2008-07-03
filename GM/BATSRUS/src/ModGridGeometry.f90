module ModGridGeometry

  ! Set coordinate information for the block adaptive grid

  use ModSize
  use ModOctreeNew

  implicit none

  private ! except

  real, public :: &
       DomainMin_D(MaxDim), & ! minimum (gen.) coordinates of the domain
       DomainMax_D(MaxDim), & ! maximum (gen.) coordinates of the domain

  real, public, dimension(MaxDim, MaxBlock) :: &
       CornerMin_DB(MaxDim, MaxBlock), & ! starting (gen.) coordinates 
       CellSize_DB(MaxDim, MaxBlock), & ! cell size in (gen.) coords
       InvVolume_DB(MaxDim, MaxBlock) ! inverse volume of cells in real space

  public:: set_block_geometry

contains

  subroutine set_block_geometry(iBlock)

    ! Set block coordinates for local block index iBlock
    ! lookup global block index first

    integer, intent(in):: iBlock   ! local block index

    iBlockAll = iBlockAll_BP(iBlock, iProc)

  end subroutine set_block_geometry

end module ModGridGeometry
