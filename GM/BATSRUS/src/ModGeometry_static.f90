!^CFG COPYRIGHT UM
Module ModGeometry
  use ModSize
  use ModMain,ONLY:body2_,ExtraBc_
  use ModIO,         ONLY: iUnitOut, write_prefix
  use ModProcMH,     ONLY: iProc

  implicit none
  SAVE

  ! Logical parameter indicating static vs. dynamic allocation
  logical, parameter :: IsDynamicGeometry = .false.

  !\
  ! Geometry parameters.
  !/
  real  ::    x1, x2, y1, y2, z1, z2
  real :: dxyz(3), xyzStart(3), xyzStart_BLK(3,nBLK),XyzMin_D(3),XyzMax_D(3)

  !\
  ! Other block solution and geometry parameters.
  !/
  real :: minDXvalue, maxDXvalue

  real, dimension(nBLK) :: dx_BLK, dy_BLK, dz_BLK, Rmin_BLK
  real, dimension(nBLK) :: Rmin2_BLK                                       !^CFG IF SECONDBODY
  real, dimension(nBLK) :: fAx_BLK, fAy_BLK, fAz_BLK, cV_BLK               !^CFG IF CARTESIAN

  real,dimension(nI,nJ,nK,nBLK):: vInv_CB         
  !^CFG IF NOT CARTESIAN BEGIN

  logical::UseCovariant=.false.


  !For a vertex-based logically cartesian (spherical, cylindircal) grid 
  !(UseVertexBasedGrid=.true.) the node coordinates are defined
  !in terms of an arbitrary pointwide transformation of nodes of an 
  !original cartesian (spherical,cylindrical) block adaptive grid.
  !Advantage: the possiblity to use the arbitrary transformation.
  !Disadvantages: the cell center coordinates can not be definied unambigously
  !and the difference of the state variables across the face does not evaluate
  !the gradient in the direction, normal to this face (stricly speaking).
  !Cell-centered grids are used if UseVertexBasedGrid=.false. (default value)
  !Advantage: for some particular geometries (spherical, cylindrical) the 
  !control volumes are the Voronoy cells (any face is perpendicular to the line
  !connecting the centers of the neighboring cells). 
  !Disadvantages: even in these particular cases it is not easy to properly define 
  !the face area vectors at the resolution change. More general cell-centered grid 
  !either is not logically cartesian, or does not consist of the Voronoy cells only.
  !
  logical :: UseVertexBasedGrid=.false.
  character (len=20) ::TypeGeometry='cartesian'                            !^CFG END CARTESIAN

 

  ! Variables describing cells inside boundaries

 
  logical,dimension(nBLK) :: BodyFlg_B 
  logical,dimension(nBLK) :: DoFixExtraBoundary_B                          !^CFG IF FACEOUTERBC

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
  real,allocatable,dimension(:,:,:,:,:):: &            !^CFG IF NOT CARTESIAN
        FaceAreaI_DFB,FaceAreaJ_DFB,FaceAreaK_DFB      !^CFG IF NOT CARTESIAN
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
  !^CFG IF NOT CARTESIAN BEGIN
  subroutine allocate_face_area_vectors
    allocate(FaceAreaI_DFB(nDim,2-gcn:nI+gcn,0:nJ+1,0:nK+1,nBLK))
    allocate(FaceAreaJ_DFB(nDim,0:nI+1,2-gcn:nJ+gcn,0:nK+1,nBLK))
    allocate(FaceAreaK_DFB(nDim,0:nI+1,0:nJ+1,2-gcn:nK+gcn,nBLK))
  end subroutine allocate_face_area_vectors
    !^CFG END CARTESIAN    

end module ModGeometry
