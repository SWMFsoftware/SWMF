!^CFG COPYRIGHT UM
!^CFG FILE NOT CARTESIAN
module ModCovariant
  use ModSize
  implicit none
  save
  ! introduced in such a manner that the face 
  ! area vector is equal to    
  ! ((x,y,z)_BLK(i,j,k,iBLK)-(x,y,z)_BLK(i-1,j,k,iBLK))&
  !                                   *FaceAreaI(i,j,k,iBLK) 
  !real :: FaceAreaI_FB(0:nI+2,0:nJ+1,0:nK+1,nBLK)                             
  !real :: FaceAreaJ_FB(0:nI+1,0:nJ+2,0:nK+1,nBLK)                       
  !real :: FaceAreaK_FB(0:nI+1,0:nJ+1,0:nK+2,nBLK)
  real,dimension(nBLK) :: &
       FaceArea2MinI_B, FaceArea2MinJ_B,FaceArea2MinK_B
  
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
  character (len=20) ::TypeGeometry='cartesian'                            
  real,allocatable,dimension(:,:,:,:,:):: &            
        FaceAreaI_DFB,FaceAreaJ_DFB,FaceAreaK_DFB      
contains
  subroutine allocate_face_area_vectors
    if(allocated(FaceAreaI_DFB))return
    allocate(FaceAreaI_DFB(nDim,2-gcn:nI+gcn,0:nJ+1,0:nK+1,nBLK))
    allocate(FaceAreaJ_DFB(nDim,0:nI+1,2-gcn:nJ+gcn,0:nK+1,nBLK))
    allocate(FaceAreaK_DFB(nDim,0:nI+1,0:nJ+1,2-gcn:nK+gcn,nBLK))
  end subroutine allocate_face_area_vectors
  subroutine allocate_old_levels
  end subroutine allocate_old_levels
  subroutine calc_face_area_i(i,j,k,iBLK,FaceAreaVector_D)
    integer,intent(in)::i,j,k,iBLK
    real,dimension(3),intent(out)::FaceAreaVector_D
    FaceAreaVector_D=FaceAreaI_DFB(:,i,j,k,iBLK)
  end subroutine calc_face_area_i
  subroutine calc_face_area_j(i,j,k,iBLK,FaceAreaVector_D)
    integer,intent(in)::i,j,k,iBLK
    real,dimension(3),intent(out)::FaceAreaVector_D
    FaceAreaVector_D=FaceAreaJ_DFB(:,i,j,k,iBLK)
  end subroutine calc_face_area_j
  subroutine calc_face_area_k(i,j,k,iBLK,FaceAreaVector_D)
    integer,intent(in)::i,j,k,iBLK
    real,dimension(3),intent(out)::FaceAreaVector_D
    FaceAreaVector_D=FaceAreaK_DFB(:,i,j,k,iBLK)
  end subroutine calc_face_area_k
end module ModCovariant

!=============End ModCovariant.f90=========================
