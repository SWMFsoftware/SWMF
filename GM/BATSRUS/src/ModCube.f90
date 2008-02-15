!^CFG COPYRIGHT UM
!
!=================================
module ModCube
  use ModSize
  use ModParallel,ONLY:NOBLK
  implicit none
  integer,parameter,dimension(3)::nCell2_D=(/nI/2,nJ/2,nK/2/)
  private::nCell2_D,i_child
 !/////////////////////////////////////////////////////////////////////////////
  !
  !    North                      8--------------7
  !    Phi,j  Top               / |            / |
  !     ^     z,k,Theta       /   | face iZ=1/   |
  !     |   /               /     |        /     |
  !     | /     West      5-------+------6       |
  !     |-----> R,i       |   e -1|      |       |
  !                       |  c =  1------+-------2
  !                       | a iX/        |     /
  !                       |f  /          |   /
  !                       | /            | /
  !                       4--------------3 --(front, right, and bottom 3D box corner)
  !
  !  Point 7 is:  +x, +y, +z, West, North, Top
  !
  !  Point 4 is:  -x, -y, -z, East, South, Bottom
  !=================================================================================!

  !The shift of the child corner with respect to the parent corner
  integer, parameter, dimension(8,3):: iShiftChild_ID = reshape(&
       !1    2     3     4    5     6     7     8        child index
       (/0,  nI/2, nI/2, 0,   0,    nI/2, nI/2, 0,     & !i shift
       0,    0,    0,    0,   nJ/2, nJ/2, nJ/2, nJ/2,  & !j shift
       nK/2, nK/2, 0,    0,   0,    0,    nK/2, nK/2/),& !k shift
       (/8,3/))
contains
  !==============================================================!
  integer function i_child(iX,iY,iZ)
    integer,intent(in)::iX,iY,iZ       !should be 0 or 1
    integer,dimension(8),parameter::Bin2Child_I=(/4,1,5,8,3,2,6,7/)
    !Transforms the binary number ix|iy|iz to the child number,
    !such that iShiftChild_ID(iChild,:)/nCell2_D=(iX,iY,iZ)
    i_child=Bin2Child_I(4*iX+2*iY+iZ+1)
  end function i_child
  !==============================================================!
  subroutine refined_children_list(iX,iY,iZ,iChildren_I)
    integer,intent(in)::iX,iY,iZ 
    !Allowed values: -1,0,1. 
    !(iX,iY,iZ) is the direction from the center of the coarser block
    !towards face (edge, corners).
    
    integer,dimension(4),intent(out)::iChildren_I
    !The child number list for neighboring blocks
    !-----------------------------------------------------
    integer::i,j,k,lSubF
    lSubf=0;iChildren_I=NOBLK
    do i=-min(0,iX),& !for iX=-1 all the nighbors are in the western part
                     1-max(0,iX)! For iX=+1, in the eastern part
       do j=-min(0,iY),&       ! For iY=-1, in the nothern part
                     1-max(0,iY)! For iY=+1, in the southern part
          do k=-min(0,iZ),&    ! For iZ=-1, in the top part
                     1-max(0,iZ)! For iZ=-1, in the bottom part
             !Part here means the part of the common parent octree block for 
             !all the finer neighbors in the direction (iX,iY,iZ)
             lSubF=lSubF+1; iChildren_I(lSubF)=i_child(i,j,k)
          end do
       end do
    end do
  end subroutine refined_children_list
  !============================================================================!
  !==============================================================
  !===============================TESTING PROCEDURES============!
  !=test verifies that subfaces arrays constructed above is in 
  !=accordance with tree_neighbors_across_pole=!
  subroutine test_refined_neighbor(iX,iY,iZ,iChild_I)
    use ModParallel,ONLY:NOBLK
    integer,intent(in)::iX,iY,iZ,iChild_I(4)
    integer::iChildrenTest_I(4)
    !--------------------------------------------------------------------------!
    call refined_children_list(iX,iY,iZ,iChildrenTest_I)
    if(any(iChildrenTest_I/=iChild_I))then
       write(*,*)'Wrong refined children list in tree_neighbor_across_pole'
       write(*,*)'The positions of children in the list should be:',&
            iChildrenTest_I
       write(*,*)'Actual children list',iChild_I
       call stop_mpi('Failed')
    end if
  end subroutine test_refined_neighbor
  !============================================================================!
end module ModCube
