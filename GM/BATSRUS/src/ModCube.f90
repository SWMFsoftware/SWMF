!^CFG COPYRIGHT UM
!
!=================================
module ModCube
  use ModSize
  use ModParallel,ONLY:NOBLK
  implicit none
  Private
  integer,parameter,dimension(3)::nCell2_D=(/nI/2,nJ/2,nK/2/)
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
  
  interface child1_at_face
     module procedure child1_at_face_vec
     module procedure child1_at_face_3s
  end interface

   interface is_not_at_face
     module procedure is_not_at_face_vec
     module procedure is_not_at_face_3s
  end interface
                                 !  vec(iDirC2F_D,iChild,IsAtFace,
  interface get_position_at_face !...3s(iX,iY,iZ,iChild,IsAtFace,...
     module procedure get_position_shift_3s  ! ....,iShift_D)
     module procedure get_position_shift_vec ! ....,iShift_D)
     module procedure get_position_subf_3s   ! ....,iSubF)
     module procedure get_position_subf_vec  ! ....,iSubF)
  end interface
  public::child1_at_face         !iDirC2F_D=>Child1
  public::is_not_at_face         !iDirC2F_D,Child=>is_not_at_face
  public::set_indices            !iDirC2F_D=>cell indices for send/receive
  !-----------------------------------------------------------------------!
  public::get_position_at_face   !iDirC2F,Child=>iSubface or index shift
  public::refined_children_list  !iDirC2F_D=>list of children at face
  public::iShiftChild_ID         
  public::test_refined_neighbor
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
  !=============is_not_at_face, child1_at_face===================!
  !In the following routines "face" means face, edge or corner    
  !iDirC2F_D is the direction from the center of the coarser block
  !towards the face, child number is for finer blocks, having the
  !common face with the coarser one. This definition of iDirC2F_D 
  !does not work for the case when the face is the the pole. For
  !this case a  more general definition for iDirC2F_D should be 
  !applied,specifically iDirC2F_D is the direction from the face 
  !to the center of the parent octree block In a 3-scalar version, 
  !iDirC2F=(/iX,iY,iZ/)
  !!=========================child1 at face =====================!
  !!At faces (except for corner, some cell indexes may change 
  !along the face. Child1 is the number of child for the block which
  !includes the starting poit for such variable indices. 
  !==============================================================!
  integer function child1_at_face_vec(iDirC2F_D)
    integer,dimension(3),intent(in)::iDirC2F_D
    child1_at_face_vec=i_child( -min(0,iDirC2F_D(1)),&
        -min(0,iDirC2F_D(2)),-min(0,iDirC2F_D(3)))
  end function child1_at_face_vec
  !==============================================================!
  integer function child1_at_face_3s(iX,iY,iZ)
    integer,intent(in)::iX,iY,iZ
    child1_at_face_3s=i_child( -min(0,iX),-min(0,iY),-min(0,iZ))
  end function child1_at_face_3s
  !==============================================================!
  logical function is_not_at_face_vec(iDirC2F_D,iChild)
    integer,dimension(3),intent(in)::iDirC2F_D
    integer,intent(in)::iChild
    is_not_at_face_vec=any(iDirC2F_D/=0.and.&
         iShiftChild_ID(iChild,:)- &
         iShiftChild_ID(child1_at_face_vec(iDirC2F_D),:)/=0)
  end function is_not_at_face_vec
  !==============================================================!
  logical function is_not_at_face_3s(iX,iY,iZ,iChild)
    integer,intent(in)::iX,iY,iZ,iChild
    is_not_at_face_3s=any((/iX,iY,iZ/)/=0.and.&
         iShiftChild_ID(iChild,:)- &
         iShiftChild_ID(child1_at_face_3s(iX,iY,iZ),:)/=0)
  end function is_not_at_face_3s
  !==============================================================!
  !==============Combined routines:  ============================!
  subroutine get_position_shift_3s(&
       iX,iY,iZ,iChild,IsAtFace,iShift_D)
    integer,intent(in)::iX,iY,iZ,iChild
    logical,intent(out)::IsAtFace
    integer,intent(out),dimension(3)::iShift_D
    !-----------------------------------------!
 
    iShift_D=0
    IsAtFace=.not.is_not_at_face_3s(iX,iY,iZ,iChild)
    if(IsAtFace)iShift_D=iShiftChild_ID(iChild,:) &
         -iShiftChild_ID(Child1_at_face_3s(iX,iY,iZ),:)
  end subroutine get_position_shift_3s
  !==============================================================!
  subroutine get_position_shift_vec(&
       iDirC2F_D,iChild,IsAtFace,iShift_D)
    integer,intent(in)::iDirC2F_D(3),iChild
    logical,intent(out)::IsAtFace
    integer,intent(out),dimension(3)::iShift_D
    !-----------------------------------------!
 
    iShift_D=0
    IsAtFace=.not.is_not_at_face_vec(iDirC2F_D,iChild)
    if(IsAtFace)iShift_D=iShiftChild_ID(iChild,:) &
         -iShiftChild_ID(Child1_at_face_vec(iDirC2F_D),:)
  end subroutine get_position_shift_vec
  !==============================================================
  !==============================================================!
  !=====================SET_INDICES==============================!
  !This is the only version of set_indices which does not use a  !
  !magic subface ordering indes. For sending messages across the !
  !pole the indexes along the direction of iDirPole should be    !
  !flipped
  subroutine set_indices(&
       iDirS2R_D,& !Direction from S(end) to R(ecv) blocks
       nLayerS,nLayerR,& !nLayers for physical- and ghost cells 
       iMinS_D,iMaxS_D,& !Indexes to be sent, for sending block
       iMinR_D,iMaxR_D,& !Indexes to be sent, for receiving block
       iLevelR,        & !Level of receiving block
       iChild,         & !Child number for finer block, if needed
       iChild1,        & !Child1 defined as described above
       nExtraCell)       !
    integer,intent(in)::iDirS2R_D(3),nLayerS,nLayerR
    integer,dimension(3),intent(out)::iMinS_D,iMaxS_D
    integer,dimension(3),intent(out)::iMinR_D,iMaxR_D
    integer,intent(in)::iLevelR
    integer,intent(in),optional::iChild,iChild1,nExtraCell
    !------------------------------------------------------------!
    integer::iShift_D(3),iDirC2F_D(3)
    logical::IsFaseIndex_D(3)
    where(iDirS2R_D==+1)   !send UP
       iMaxS_D=nCells   ; iMinS_D=iMaxS_D+1-nLayerS
       iMaxR_D=0         ; iMinR_D=iMaxR_D+1-nLayerR
    end where
    where(iDirS2R_D==-1)   !send DOWN
       iMinS_D=1         ; iMaxS_D=iMinS_D-1+nLayerS
       iMinR_D=nCells+1 ; iMaxR_D=iMinR_D-1+nLayerR
    end where
    select case(iLevelR) ! Indices along the face
    case(0)              ! Recv block is at the same level
       where(iDirS2R_D==0)
          iMinS_D=1      ; iMaxS_D=nCells
          iMinR_D=1      ; iMaxR_D=nCells
       end where
    case(-1)             ! Recv block is finer
       if(.not.present(iChild1))then
          iShift_D=iShiftChild_ID(iChild,:)-&
               iShiftChild_ID(child1_at_face_vec(iDirS2R_D),:)
       else
          iShift_D=iShiftChild_ID(iChild,:)-&
               iShiftChild_ID(iChild1,:)
       end if
       if(.not.present(nExtraCell))then
          where(iDirS2R_D==0)
             iMinS_D=1 + iShift_D ; iMaxS_D=nCell2_D + iShift_D
             iMinR_D=1            ; iMaxR_D=nCells
          end where
       else
          where(iDirS2R_D==0.and.iShift_D==0)
             iMinS_D=1            ; iMaxS_D=nCell2_D + nExtraCell
             iMinR_D=1            ; iMaxR_D=nCells +2*nExtraCell
          end where
          where(iDirS2R_D==0.and.iShift_D/=0)
             iMinS_D=1 - nExtraCell + nCell2_D   ; iMaxS_D=nCells 
             iMinR_D=1 - 2 * nExtraCell          ; iMaxR_D=nCells 
          end where
       end if
    case(+1)    !Receiveng block is coarser
       if(.not.present(iChild1))then
          iDirC2F_D=-iDirS2R_D
          iShift_D=iShiftChild_ID(iChild,:)-&
               iShiftChild_ID(child1_at_face_vec(iDirC2F_D),:)
       else
          iShift_D=iShiftChild_ID(iChild,:)-&
               iShiftChild_ID(iChild1,:)
       end if
       where(iDirS2R_D==0)
          iMinS_D=1             ; iMaxS_D=nCells 
          iMinR_D=1 + iShift_D  ; iMaxR_D=nCell2_D + iShift_D
       end where
    end select
  end subroutine set_indices
  !===================================================================!
  !The further routines are based on the order of children in the list!
  !created by routines for finding neighbors                          !
  subroutine refined_children_list(iX,iY,iZ,iChildren_I)
    integer,intent(in)::iX,iY,iZ 
    !Allowed values: -1,0,1. 
    !(iX,iY,iZ) is the direction from the center of the COARSER block
    !towards face (edge, corners), i.e. towards FINER neighbors.
    
    integer,dimension(4),intent(out)::iChildren_I
    !The child number list for neighboring blocks
    !-----------------------------------------------------
    integer::i,j,k,lSubF,iChild2
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
    if(.not.(iZ**2==1.and.iX**2+iY**2==0))return
    !consider separatelely  Top_ or Bot_ face
    iChild2=iChildren_I(2)
    iChildren_I(2)=iChildren_I(3)
    iChildren_I(3)=iChild2
  end subroutine refined_children_list
  !==============================================================!
  subroutine get_position_subf_3s(iX,iY,iZ,iChild,IsAtFace,iSubF)
    integer,intent(in)::iX,iY,iZ,iChild
    logical,intent(out)::IsAtFace
    integer,intent(out)::iSubF
    !-----------------------------------------!
    integer,dimension(4)::iChild_I
    integer,dimension(1)::iLoc
    iSubf=0
    call refined_children_list(iX,iY,iZ,iChild_I)
    IsAtFace=any(iChild_I==iChild)
    if(.not.IsAtFace)return
    iLoc=minloc(iChild_I,MASK=iChild_I==iChild)
    iSubF=iLoc(1)
  end subroutine get_position_subf_3s
  !=============================================================!
  subroutine get_position_subf_vec(iDirC2F_D,iChild,IsAtFace,iSubF)
    integer,intent(in)::iDirC2F_D(3),iChild
    logical,intent(out)::IsAtFace
    integer,intent(out)::iSubF
    !-----------------------------------------!
    call get_position_subf_3s(iDirC2F_D(1),iDirC2F_D(2),&
     iDirC2F_D(3),iChild,IsAtFace,iSubF)
  end subroutine get_position_subf_vec
  !==============================================================
  !===============================TESTING PROCEDURES============!
  !=test verifies that subfaces arrays constructed above is in 
  !=accordance with tree_neighbors_across_pole=!
  subroutine test_refined_neighbor(iX,iY,iZ,iChild_I)
    use ModParallel,ONLY:NOBLK
    integer,intent(in)::iX,iY,iZ,iChild_I(4)
    integer::iChildrenTest_I(4)
    !------------------------------------------------------------------------!
    call refined_children_list(iX,iY,iZ,iChildrenTest_I)
    if(any(iChildrenTest_I/=iChild_I))then
       write(*,*)'Wrong refined children list in tree_neighbor_list'
       write(*,*)'The positions of children in the list should be:',&
            iChildrenTest_I
       write(*,*)'Actual children list',iChild_I
       call stop_mpi('Failed')
    end if
  end subroutine test_refined_neighbor
end module ModCube
  
