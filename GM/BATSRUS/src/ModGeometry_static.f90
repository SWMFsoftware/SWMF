!^CFG COPYRIGHT UM
Module ModGeometry
  use ModSize
  use ModMain,ONLY:body2_,ExtraBc_
  use ModIO,         ONLY: iUnitOut, write_prefix
  use ModProcMH,     ONLY: iProc
  use ModCovariant                

  implicit none
  SAVE

  ! Logical parameter indicating static vs. dynamic allocation
  logical, parameter :: IsDynamicGeometry = .false.

  !\
  ! Geometry parameters.
  !/
  real  ::    x1, x2, y1, y2, z1, z2
  real :: dxyz(3), xyzStart(3), xyzStart_BLK(3,nBLK),XyzMin_D(3),XyzMax_D(3)

  ! Mirror symmetry in the 3 directions
  integer:: nMirror_D(3) = 1

  !\
  ! Other block solution and geometry parameters.
  !/
  real :: minDXvalue, maxDXvalue

  real, dimension(nBLK) :: dx_BLK, dy_BLK, dz_BLK, Rmin_BLK
  real, dimension(nBLK) :: Rmin2_BLK                                       !^CFG IF SECONDBODY
  real, dimension(nBLK) :: fAx_BLK, fAy_BLK, fAz_BLK, cV_BLK

  real,dimension(nI,nJ,nK,nBLK):: vInv_CB         
 

  ! Variables describing cells inside boundaries

 
  logical,dimension(nBLK) :: BodyFlg_B 
  logical,dimension(nBLK) :: DoFixExtraBoundary_B                          

  !true when at least one cell in the block (including ghost cells) is not true
  logical :: body_BLK(nBLK)

  ! true when all cells in block (not including ghost cells) are true_cells 
  logical :: true_BLK(nBLK)

  ! true cells are cells that are not inside a body
  logical,dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK) ::true_cell
  logical,dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,body2_:Top_) ::&          
       IsBoundaryCell_GI
  logical,dimension(body2_:Top_,nBLK):: IsBoundaryBlock_IB 
  integer :: MinBoundary=Top_, MaxBoundary=body2_                    
  logical :: far_field_BCs_BLK(nBLK)                                                

  ! Block cell coordinates
  real,  dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn,nBLK) :: &
       x_BLK,y_BLK,z_BLK,R_BLK
  real,  dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn,nBLK) :: &   !^CFG IF SECONDBODY
       R2_BLK                                                            !^CFG IF SECONDBODY
contains
  !============================================================================
  subroutine init_mod_geometry

    if(IsDynamicGeometry .and. iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'init_mod_geometry allocated arrays'
    end if

  end subroutine init_mod_geometry
  !============================================================================
  subroutine clean_mod_geometry

    if(IsDynamicGeometry .and. iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'clean_mod_geometry deallocated arrays'
    end if

  end subroutine clean_mod_geometry

end module ModGeometry
