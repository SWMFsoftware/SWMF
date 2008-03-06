!^CFG COPYRIGHT UM
!
!=================================
module ModCube
  use ModSize
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
  integer,dimension(8),parameter::iBin2Child_I=(/4,1,5,8,3,2,6,7/)
  public::is_not_at_face         !iDirC2F_D,Child=>is_not_at_face
  public::set_indices            !iDirC2F_D=>cell indices for send/receive
  !-----------------------------------------------------------------------!
  public::iShiftChild_ID
  public::get_children_list
contains
  !==============================================================!
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
  logical function is_not_at_face(iDirC2F_D,iChild)
    integer,dimension(3),intent(in)::iDirC2F_D
    integer,intent(in)::iChild
    integer::iBin,iChild1
    iBin=1-&
         (4*min(0,iDirC2F_D(1))+&
          2*min(0,iDirC2F_D(2))+&
            min(0,iDirC2F_D(3)))
    iChild1=iBin2Child_I(iBin)
    is_not_at_face=any(iDirC2F_D/=0.and.&
         iShiftChild_ID(iChild,:)- &
         iShiftChild_ID(iChild1,:)/=0)
  end function is_not_at_face
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
    integer::iShift_D(3),iChild1Here,iBin
    logical::IsFaseIndex_D(3)
    where(iDirS2R_D==+1)   !send UP
       iMaxS_D=nCells    
       iMinS_D=iMaxS_D+1-nLayerS
       iMaxR_D=0         
       iMinR_D=iMaxR_D+1-nLayerR
    end where
    where(iDirS2R_D==-1)   !send DOWN
       iMinS_D=1         
       iMaxS_D=iMinS_D-1+nLayerS
       iMinR_D=nCells+1  
       iMaxR_D=iMinR_D-1+nLayerR
    end where
    select case(iLevelR) ! Indices along the face
    case(0)              ! Recv block is at the same level
       where(iDirS2R_D==0)
          iMinS_D=1      
          iMaxS_D=nCells
          iMinR_D=1      
          iMaxR_D=nCells
       end where
    case(-1)             ! Recv block is finer
       if(.not.present(iChild1))then
          iBin=1-(&
               4*min(0,iDirS2R_D(1))+&
               2*min(0,iDirS2R_D(2))+&
               min(0,iDirS2R_D(3)))
          iChild1Here=iBin2Child_I(iBin)
          iShift_D=iShiftChild_ID(iChild,:)-&
               iShiftChild_ID(iChild1Here,:)
       else
          iShift_D=iShiftChild_ID(iChild,:)-&
               iShiftChild_ID(iChild1,:)
       end if
       if(.not.present(nExtraCell))then
          where(iDirS2R_D==0)
             iMinS_D=1 + iShift_D 
             iMaxS_D=nCell2_D + iShift_D
             iMinR_D=1            
             iMaxR_D=nCells
          end where
       else
          where(iDirS2R_D==0.and.iShift_D==0)
             iMinS_D=1            
             iMaxS_D=nCell2_D + nExtraCell
             iMinR_D=1            
             iMaxR_D=nCells +2*nExtraCell
          end where
          where(iDirS2R_D==0.and.iShift_D/=0)
             iMinS_D=1 - nExtraCell + nCell2_D   
             iMaxS_D=nCells 
             iMinR_D=1 - 2 * nExtraCell          
             iMaxR_D=nCells 
          end where
       end if
    case(+1)    !Receiveng block is coarser
       if(.not.present(iChild1))then
          iBin=1-(&
               4*min(0,-iDirS2R_D(1))+&
               2*min(0,-iDirS2R_D(2))+&
                 min(0,-iDirS2R_D(3)))
          iChild1Here=iBin2Child_I(iBin)
          iShift_D=iShiftChild_ID(iChild,:)-&
               iShiftChild_ID(iChild1Here,:)
       else
          iShift_D=iShiftChild_ID(iChild,:)-&
               iShiftChild_ID(iChild1,:)
       end if
       where(iDirS2R_D==0)
          iMinS_D=1             
          iMaxS_D=nCells 
          iMinR_D=1 + iShift_D  
          iMaxR_D=nCell2_D + iShift_D
       end where
    end select
  end subroutine set_indices
!===================================================================!
  !The further routines are based on the order of children in the list!
  !created by routines for finding neighbors                          !
  subroutine get_children_list(iX,iY,iZ,iChildren_I)
    use ModParallel,ONLY:NOBLK
    integer,intent(in)::iX,iY,iZ 
    !Allowed values: -1,0,1. 
    !(iX,iY,iZ) is the direction from the center of the COARSER block
    !towards face (edge, corners), i.e. towards FINER neighbors.
    
    integer,dimension(4),intent(out)::iChildren_I
    !The child number list for neighboring blocks
    !-----------------------------------------------------
    integer::i,j,k,lSubF,iChild2,iBin
    lSubf=0;iChildren_I=NOBLK
    do i=-min(0,iX),& !for iX=-1 all the nighbors are in the western part
                     1-max(0,iX)! For iX=+1, in the eastern part
       do j=-min(0,iY),&       ! For iY=-1, in the nothern part
                     1-max(0,iY)! For iY=+1, in the southern part
          do k=-min(0,iZ),&    ! For iZ=-1, in the top part
                     1-max(0,iZ)! For iZ=+1, in the bottom part
             !Part here means the part of the common parent octree block for 
             !all the finer neighbors in the direction (iX,iY,iZ)
             lSubF= lSubF + 1
             iBin = 1 + 4*i + 2*j + k
             iChildren_I(lSubF)=iBin2Child_I(iBin)
          end do
       end do
    end do
    if(.not.(iZ**2==1.and.iX**2+iY**2==0))return
    !consider separatelely  Top_ or Bot_ face
    iChild2=iChildren_I(2)
    iChildren_I(2)=iChildren_I(3)
    iChildren_I(3)=iChild2
  end subroutine get_children_list
end module ModCube
  
