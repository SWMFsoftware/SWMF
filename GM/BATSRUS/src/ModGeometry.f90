!^CFG COPYRIGHT UM
module ModGeometry
  use ModSize
  use ModMain,       ONLY: UseBody2, body2_, ExtraBc_
  use ModIO,         ONLY: iUnitOut, write_prefix
  use ModProcMH,     ONLY: iProc
  use ModCovariant                

  implicit none
  SAVE

  !\
  ! Geometry parameters.
  !/
  real :: x1, x2, y1, y2, z1, z2
  real :: DomainVolume = -1.0
  real :: dxyz(3)
  real :: xyzStart(3)
  real :: xyzStart_BLK(3,MaxBlock)
  real :: XyzMin_D(3)
  real :: XyzMax_D(3)

  ! Mirror symmetry in the 3 directions
  integer:: nMirror_D(3) = 1

  !\
  ! Other block solution and geometry parameters.
  !/
  real :: minDXvalue, maxDXvalue

  real, dimension(MaxBlock) :: dx_BLK, dy_BLK, dz_BLK, Rmin_BLK
  real, dimension(MaxBlock) :: Rmin2_BLK
  real, dimension(MaxBlock) :: fAx_BLK, fAy_BLK, fAz_BLK, cV_BLK

  real, allocatable :: vInv_CB(:,:,:,:)
 

  ! Variables describing cells inside boundaries
 
  logical,dimension(MaxBlock) :: BodyFlg_B 
  logical,dimension(MaxBlock) :: DoFixExtraBoundary_B                          

  !true when at least one cell in the block (including ghost cells) is not true
  logical :: body_BLK(MaxBlock)

  ! true when all cells in block (not including ghost cells) are true_cells 
  logical :: true_BLK(MaxBlock)

  ! true cells are cells that are not inside a body
  logical, allocatable :: true_cell(:,:,:,:)
  logical,dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,body2_:Top_):: &
       IsBoundaryCell_GI
  logical,dimension(body2_:Top_,MaxBlock):: IsBoundaryBlock_IB 
  integer :: MinBoundary=Top_, MaxBoundary=body2_                    
  logical :: far_field_BCs_BLK(MaxBlock)

  ! Block cell coordinates
  real, allocatable :: x_BLK(:,:,:,:)
  real, allocatable :: y_BLK(:,:,:,:)
  real, allocatable :: z_BLK(:,:,:,:)
  real, allocatable :: R_BLK(:,:,:,:)
  real, allocatable :: R2_BLK(:,:,:,:)

contains
  !============================================================================
  subroutine init_mod_geometry

    if(allocated(vInv_CB)) return
    allocate(vInv_CB(nI,nJ,nK,MaxBlock))
    allocate(true_cell(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(x_BLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(y_BLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(z_BLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(R_BLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    if(UseBody2) allocate(R2_BLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    if(iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'init_mod_geometry allocated arrays'
    end if

  end subroutine init_mod_geometry
  !============================================================================
  subroutine clean_mod_geometry

    if(.not.allocated(vInv_CB)) RETURN

    deallocate(vInv_CB)
    if(allocated(true_cell)) deallocate(true_cell)
    if(allocated(x_BLK))     deallocate(x_BLK)
    if(allocated(y_BLK))     deallocate(y_BLK)
    if(allocated(z_BLK))     deallocate(z_BLK)
    if(allocated(R_BLK))     deallocate(R_BLK)
    if(allocated(R2_BLK))    deallocate(R2_BLK)

    if(iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'clean_mod_geometry deallocated arrays'
    end if

  end subroutine clean_mod_geometry

end module ModGeometry
