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
  real :: FaceAreaI_FB(0:nI+2,0:nJ+1,0:nK+1,nBLK)                             
  real :: FaceAreaJ_FB(0:nI+1,0:nJ+2,0:nK+1,nBLK)                       
  real :: FaceAreaK_FB(0:nI+1,0:nJ+1,0:nK+2,nBLK)
  real,dimension(nBLK) :: &
       FaceArea2MinI_B, FaceArea2MinJ_B,FaceArea2MinK_B
contains
  subroutine calc_faceareaI(i,j,k,iBLK,FaceAreaVector_D)
    use ModGeometry,ONLY:x_BLK,y_BLK,z_BLK
    integer,intent(in)::i,j,k,iBLK
    real,dimension(3),intent(out)::FaceAreaVector_D
    FaceAreaVector_D(1)=(x_BLK(i,j,k,iBLK)-x_BLK(i-1,j,k,iBLK))&
         *FaceAreaI_FB(i,j,k,iBLK)
    FaceAreaVector_D(2)=(y_BLK(i,j,k,iBLK)-y_BLK(i-1,j,k,iBLK))&
         *FaceAreaI_FB(i,j,k,iBLK)
    FaceAreaVector_D(3)=(z_BLK(i,j,k,iBLK)-z_BLK(i-1,j,k,iBLK))&
         *FaceAreaI_FB(i,j,k,iBLK)
  end subroutine calc_faceareaI
  subroutine calc_faceareaJ(i,j,k,iBLK,FaceAreaVector_D)
    use ModGeometry,ONLY:x_BLK,y_BLK,z_BLK
    integer,intent(in)::i,j,k,iBLK
    real,dimension(3),intent(out)::FaceAreaVector_D
    FaceAreaVector_D(1)=(x_BLK(i,j,k,iBLK)-x_BLK(i,j-1,k,iBLK))&
         *FaceAreaJ_FB(i,j,k,iBLK)
    FaceAreaVector_D(2)=(y_BLK(i,j,k,iBLK)-y_BLK(i,j-1,k,iBLK))&
         *FaceAreaJ_FB(i,j,k,iBLK)
    FaceAreaVector_D(3)=(z_BLK(i,j,k,iBLK)-z_BLK(i,j-1,k,iBLK))&
         *FaceAreaJ_FB(i,j,k,iBLK)
  end subroutine calc_faceareaJ
  subroutine calc_faceareaK(i,j,k,iBLK,FaceAreaVector_D)
    use ModGeometry,ONLY:x_BLK,y_BLK,z_BLK
    integer,intent(in)::i,j,k,iBLK
    real,dimension(3),intent(out)::FaceAreaVector_D
    FaceAreaVector_D(1)=(x_BLK(i,j,k,iBLK)-x_BLK(i,j,k-1,iBLK))&
         *FaceAreaK_FB(i,j,k,iBLK)
    FaceAreaVector_D(2)=(y_BLK(i,j,k,iBLK)-y_BLK(i,j,k-1,iBLK))&
         *FaceAreaK_FB(i,j,k,iBLK)
    FaceAreaVector_D(3)=(z_BLK(i,j,k,iBLK)-z_BLK(i,j,k-1,iBLK))&
         *FaceAreaK_FB(i,j,k,iBLK)
  end subroutine calc_faceareaK
end module ModCovariant

!=============End ModCovariant.f90=========================
