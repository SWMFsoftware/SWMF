!^CFG COPYRIGHT UM
!
!==========================================================================!

module ModPolarNeighbor
  use ModCube
  use ModNumConst
  implicit none
  save
 
  interface check_pole_inside_b
     module procedure check_pole_inside_blk
     module procedure check_pole_inside_blkpe
     module procedure check_pole_inside_octreeblk
  end interface
contains
  !===============================finding polar neighbors==================!
  !=====tree_neighbor_fixed: inoptimal, provided here to demonstrate the way 
  !to apply tree_neighbor_across_pole==========! 
  subroutine tree_neighbor_fixed(iProc,iBlock,iX,iY,iZ,&
       iPE_I, iBLK_I, iChild_I, iLevel)
    integer,intent(in)::iProc,iBlock,iX,iY,iZ
    integer,intent(out),dimension(4)::iPE_I, iBLK_I, iChild_I
    integer,intent(out)::iLevel
    !--------------------------------------------------------!
    integer::iDir_D(3),iDirPole,iLoopPole
    logical::IsPole

    call check_pole_inside_b(iBlock,iProc,IsPole,iDirPole,iLoopPole)
 
    iDir_D=(/iX,iY,iZ/)
    
    if(IsPole.and.iDir_D(iDirPole)==iLoopPole)then
    
       call tree_neighbor_across_pole(iProc,iBlock,iX,iY,iZ,iDirPole,&
            iPE_I, iBLK_I, iChild_I, iLevel)
    else
       
       call treeNeighbor(iProc,iBlock,iX,iY,iZ, &
       iPE_I, iBLK_I, iChild_I, iLevel)
    end if
  end subroutine tree_neighbor_fixed
  !=============================================================================================!
  !=====tree_neighbor_across_pole, requires that the direction of iX,iY,iZ is towards the pole=!
  subroutine tree_neighbor_across_pole(iProc,iBlock,iX,iY,iZ,iDirPole,&
       iPE_I, iBLK_I, iChild_I, iLevel)
    use ModParallel
    integer,intent(in)::iProc,iBlock,iX,iY,iZ,iDirPole
    integer,intent(out),dimension(4)::iPE_I, iBLK_I, iChild_I
    integer,intent(out)::iLevel
    logical::IsPole
    integer::iDirPoleCheck,iLoopPole,iDir_D(3)
    integer::i,iPut
    
    iDir_D=(/iX,iY,iZ/)
    !Check if the direction to the neighbor is across the pole
    call check_pole_inside_blkpe(iBlock,iProc,IsPole,iDirPoleCheck,iLoopPole)
    if(.not.(IsPole.and.iDirPoleCheck==iDirPole&
         .and.iDir_D(iDirPole)==iLoopPole))call stop_mpi(&
         'tree_neighbor_across_pole:no pole from a given BLK to a given dir')
    iDir_D(iDirPole)=0
  
    call treeNeighbor(iProc,iBlock, iDir_D(1), iDir_D(2), iDir_D(3), &
         iPE_I, iBLK_I, iChild_I,iLevel)
  
    select case(iLevel)
    case(NOBLK)
       return !The neighbor is out of tree
    case(0,1)
       call find_axial_neighbor(iPE_I(1),iBLK_I(1),iPE_I(1),iBLK_I(1))
       return
    case(-1)
       iDir_D=(/iX,iY,iZ/)
       !In spherical geometry (iDirPole==3) finer block may only occur at
       !iX=+-1, in axial one (iDirPole==1) at iZ==+-1
       if(.not.iDir_D(4-iDirPole)**2==1)call stop_mpi(&
            'Failed:  tree_neighbor_across_pole')
       !Sort out subfaces of the blocks which do not touch the pole
 
       do i=1,4
          if(iBLK_I(i)==NOBLK)exit
          call check_pole_inside_blkpe(&
               iBLK_I(i),iPE_I(i),IsPole,iDirPoleCheck,iLoopPole)
          if(.not.(IsPole.and.iDirPoleCheck==iDirPole&
               .and.iDir_D(iDirPole)==iLoopPole))then
             iPE_I(i)=NOBLK; iBLK_I(i)=NOBLK; iChild_I(i)=NOBLK
          else
             call find_axial_neighbor(iPE_I(i),iBLK_I(i),iPE_I(i),iBLK_I(i))
             !Move the earlier created NOBLK lines to top position,if needed
             iPut=i
             do while(iPut>1)
                if(iBLK_I(iPut-1)/=NOBLK)exit
                iPut=iPut-1
             end do
             if(iPut/=i)then
                iPE_I(iPut)=iPE_I(i)
                iBLK_I(iPut)=iBLK_I(i)
                iChild_I(iPut)=iChild_I(i)
                iPE_I(i)=NOBLK; iBLK_I(i)=NOBLK; iChild_I(i)=NOBLK
             end if
          end if
       end do
    end select
  end subroutine tree_neighbor_across_pole
  !=============================================================================================!
  !Simplified version, gives only iPEOut and iBLKOut for the block which has the common section !
  !of the pole with iPEIn,iBLKIn (face neighbor, the face being degenerated and having zero area!
  !=============================================================================================! 
  subroutine find_axial_neighbor(iPEIn,iBLKIn,iPEOut,iBLKout)
    use ModParallel, ONLY : proc_dims
    use ModOctree
    implicit none
    integer,intent(in)::iPEIn,iBLKIn
    integer,intent(out)::iPEOut,iBLKOut
    integer,dimension(99)::iChild_I
    integer::i,j,k,iRoot,jRoot,kRoot,iLevel,iLevelIn
    type (adaptive_block_ptr) :: InBlkPtr,OutBlkPtr

    nullify (InBlkPtr % ptr)
    InBlkPtr%ptr=> global_block_ptrs(iBlkIn,iPEIn+1)%ptr
    iLevelIn=InBlkPtr%ptr%LEV

    do iLevel=iLevelIn,1,-1
       iChild_I(iLevel)=InBlkPtr%ptr%child_number
       InBlkPtr%ptr=>InBlkPtr%ptr%parent%ptr
    end do
    !Found root cell, find neighbor root cell

    !Found root cell, find neighbor root cell
    iRoot = InBlkPtr%ptr%iRoot 
    jRoot = InBlkPtr%ptr%jRoot
    kRoot = InBlkPtr%ptr%kRoot
    !Check
    if (.not.associated(InBlkPtr%ptr,&
       octree_roots(iRoot,jRoot,kRoot)%ptr))&
       call stop_mpi('Failure in the algorithm for finding octree root')
    
    jRoot=1+mod(jRoot-1+proc_dims(2)/2,proc_dims(2))

    nullify (OutBlkPtr%ptr)
    
    OutBlkPtr%ptr=>octree_roots(iRoot,jRoot,kRoot)%ptr
    do iLevel=1,iLevelIn
       
       OutBlkPtr%ptr => OutBlkPtr%ptr%child(iChild_I(iLevel))%ptr
       if (.not.associated(OutBlkPtr%ptr))&
            call stop_mpi('Wrong Axial Neighbor is found') ! ensure block is allocated
    end do
    iPEOut  = outBlkPtr%ptr%PE
    iBlkOut = outBlkPtr%ptr%BLK
  end subroutine find_axial_neighbor
  !================================CHECK if block is polar===================================!
  logical function is_pole_inside_b(iBlock)
    use ModGeometry,ONLY:IsPole_B=>DoFixExtraBoundary_B,is_axial_geometry
    implicit none
    integer,intent(in)::iBlock
    if(.not.is_axial_geometry())then
       is_pole_inside_b=.false.
       return
    end if
    is_pole_inside_b=IsPole_B(iBlock)
  end function is_pole_inside_b
  !==========================================================================================!
  subroutine check_pole_inside_blk(iBlock,IsPole,iDir,iLoop)
    use ModOctree
    use ModGeometry,ONLY:dz_BLK,XyzStart_BLK,TypeGeometry
    use ModMain,ONLY:Theta_
    implicit none
    integer,intent(in) ::iBlock
    logical,intent(out)::IsPole   !True if the pole is inside
    integer,intent(out)::iDir     !iDir=1 pole is in the direction 'i'
    integer,intent(out)::iLoop    !iLoop=+1(-1) in positive (negative) 
    !direction, from a given block
    !--------------------------------------------------------------!
    real::dTheta

    IsPole=is_pole_inside_b(iBlock)
    if(IsPole)then
       if(index(TypeGeometry,'spherical')>0)then
          !calculate a minimal non-zero face area squared along theta direction
          dTheta=dz_BLK(iBlock); iDir=3
          if(XyzStart_BLK(Theta_,iBlock)-dTheta < -cHalfPi)then
             iLoop=-1
          else
             iLoop=+1
          end if
       else !Axial geometry
          iDir=1;iLoop=-1
       end if
    else
       iDir=2;iLoop=0  !no pole
    end if
  end subroutine check_pole_inside_blk
  !==========================================================================================!
  subroutine check_pole_inside_blkpe(iBlock,iProc,IsPole,iDir,iLoop)
    use ModOctree
    implicit none
    integer,intent(in)::iBlock,iProc
    logical,intent(out)::IsPole   !True if the pole is inside
    integer,intent(out)::iDir     !iDir=1 pole is in the direction 'i'
    integer,intent(out)::iLoop    !iLoop=+1(-1) in positive (negative) 
    type(adaptive_block_ptr)::BlkAux
    !direction, from a given block
    if(.not.associated(global_block_ptrs(iBlock, iProc+1) % ptr))then
       write(*,*)'Not associated global block pointer for iBlock=',iBlock,&
            ' at PE=',iProc
       call stop_mpi('Stop in check_pole_inside_blk')
    end if
    BlkAux%ptr=>global_block_ptrs(iBlock,iProc+1)%ptr
    call check_pole_inside_octreeblk(BlkAux,&
         IsPole,iDir,iLoop)
  end subroutine check_pole_inside_blkpe
  !=========================================================================!
  subroutine check_pole_inside_octreeblk(OctreeBlk,IsPole,iDir,iLoop)
    use ModCovariant,ONLY:TypeGeometry,is_axial_geometry
    use ModParallel, ONLY : proc_dims
    use ModOctree
    implicit none
    type(adaptive_block_ptr),intent(in)::OctreeBlk
    type(adaptive_block_ptr)::BlkAux
    logical,intent(out)::IsPole   !True if the pole is inside
    integer,intent(out)::iDir     !iDir=1 pole is in the direction 'i'
    integer,intent(out)::iLoop    !iLoop=+1(-1) in positive (negative) 
    !direction, from a given block

    if(.not.OctreeBlk%ptr % used)&
         call stop_mpi('Attempt to check pole inside unused octree block')
    IsPole=.false.;iDir=2;iLoop=0
    if(.not.is_axial_geometry())return
    if(.not.OctreeBlk%ptr%IsExtraBoundaryOrPole)return
    if(index(TypeGeometry,'spherical')>0)then
       if(proc_dims(3)==1.and.OctreeBlk%ptr%LEV==0)call stop_mpi(&
            'Can not check pole: the block covers both hemispheres')

       !for        iz=-1      iy=0       ix=0

       call find_tree_neighbor(OctreeBlk,BlkAux,0,0,-1,IsPole)
       if(IsPole)then
          iDir=3; iLoop=-1; return
       end if

       !for        iz=+1      iy=0       ix=0  

       call find_tree_neighbor(OctreeBlk,BlkAux,0,0,+1,IsPole)   
       if(IsPole)then
          iDir=3; iLoop=+1; return
       end if
    elseif(index(TypeGeometry,'cylindrical')>0.or.&
         index(TypeGeometry,'axial')>0)then

       !for        iz=0       iy=0       ix=-1

       call find_tree_neighbor(OctreeBlk,BlkAux,-1,0,0,IsPole)   
       if(IsPole)then
          iDir=1; iLoop=-1 ;return
       end if
    else
       call stop_mpi('Failed check_pole_inside_octreeblk')
    end if
    call stop_mpi('Failed check_pole_inside_octreeblk')
  end subroutine check_pole_inside_octreeblk
end module ModPolarNeighbor
