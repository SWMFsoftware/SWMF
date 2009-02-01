!^CFG COPYRIGHT UM
module ModCovariant
  use ModSize
  use ModNumConst
  use ModMain,ONLY:unusedBLK
  use ModCoordTransform
  use ModParallel,ONLY:BLKneighborLEV,NOBLK
  use ModUtilities,ONLY: check_allocate
  implicit none
  save
 
  real,dimension(nBLK) :: &
       FaceArea2MinI_B, FaceArea2MinJ_B, FaceArea2MinK_B
  
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
  logical :: UseVertexBasedGrid=.true.
  character (len=20) ::TypeGeometry='cartesian'                            
  real,allocatable,dimension(:,:,:,:,:):: &            
        FaceAreaI_DFB,FaceAreaJ_DFB,FaceAreaK_DFB
  integer,allocatable,dimension(:,:,:,:)::OldLevel_IIIB
  logical,dimension(-1:1,-1:1,-1:1)::IsNotCorner_III
  
  !Parameters of 
  real::rTorusLarge=6.0,rTorusSmall=0.50
  integer,parameter::nToroidalBoundaryPoints=1024
  !
  real,dimension(0:nToroidalBoundaryPoints)::TorusSurface_I
  logical::IsInitializedTorusGeometry=.false.
contains
  subroutine allocate_face_area_vectors
    integer::iError
    if(allocated(FaceAreaI_DFB))return
    allocate(FaceAreaI_DFB(nDim,1:nI+1,1:nJ,1:nK,nBLK),stat=iError)
    call check_allocate(iError,'Face Areas are not allocated')
    allocate(FaceAreaJ_DFB(nDim,1:nI,1:nJ+1,1:nK,nBLK),stat=iError)
    call check_allocate(iError,'Face Areas are not allocated')
    allocate(FaceAreaK_DFB(nDim,1:nI,1:nJ,1:nK+1,nBLK),stat=iError)
    call check_allocate(iError,'Face Areas are not allocated')
    FaceAreaI_DFB=cOne; FaceAreaJ_DFB=cOne; FaceAreaK_DFB=cOne
  end subroutine allocate_face_area_vectors
!---------------------------------------------------------------------------------
  subroutine allocate_old_levels
    if(allocated(OldLevel_IIIB))return
    allocate(OldLevel_IIIB(-1:1,-1:1,-1:1,nBLK))
    OldLevel_IIIB=NOBLK
    IsNotCorner_III=.true.
    IsNotCorner_III(-1,-1,-1)=.false.
    IsNotCorner_III(+1,-1,-1)=.false.
    IsNotCorner_III(-1,+1,-1)=.false.
    IsNotCorner_III(+1,+1,-1)=.false.
    IsNotCorner_III(-1,-1,+1)=.false.
    IsNotCorner_III(+1,-1,+1)=.false.
    IsNotCorner_III(-1,+1,+1)=.false.
    IsNotCorner_III(+1,+1,+1)=.false.
  end subroutine allocate_old_levels
!---------------------------------------------------------------------------------
  logical function do_fix_geometry_at_reschange(iBlock)
    integer,intent(in)::iBlock
    if(unusedBLK(iBlock).or.(.not.UseCovariant).or.(.not.UseVertexBasedGrid))then
       do_fix_geometry_at_reschange=.false.
       return
    end if
    do_fix_geometry_at_reschange=any(IsNotCorner_III(:,:,:).and.&
         OldLevel_IIIB(:,:,:,iBlock)/=BLKneighborLEV(:,:,:,iBlock).and.&
         (OldLevel_IIIB(:,:,:,iBlock)==-1.or.BLKneighborLEV(:,:,:,iBlock)==-1))
  end function do_fix_geometry_at_reschange

!---------------------------------------------------------------------------------
  subroutine get_face_area_i(&
       XyzNodes_DIII,&                      !(in) Cartesian coordinates of nodes
       iStart,iMax,jStart,jMax,kStart,kMax,&!(in) FACE indexes. 
       FaceAreaI_DF)                        !(out)Face area vectors
    !-----------------------------------------------------------------------------
    integer,intent(in)::iStart,iMax,jStart,jMax,kStart,kMax
    real,intent(in),dimension(nDim,iStart:iMax,jStart:jMax+1,kStart:kMax+1)::&
         XyzNodes_DIII
    real,intent(out),dimension(nDim,iStart:iMax,jStart:jMax,kStart:kMax)::&
         FaceAreaI_DF
    !-----------------------------------------------------------------------------
    integer::i,j,k
    do k=kStart,kMax; do j=jStart,jMax; do i=iStart,iMax
       FaceAreaI_DF(:,i,j,k)=cHalf*&
            cross_product(XyzNodes_DIII(:,i,j+1,k+1)-XyzNodes_DIII(:,i,j  ,k),&
                          XyzNodes_DIII(:,i,j  ,k+1)-XyzNodes_DIII(:,i,j+1,k))
    end do; end do; end do
  end subroutine get_face_area_i
!---------------------------------------------------------------------------------
  subroutine get_face_area_j(&
       XyzNodes_DIII,&                      !(in) Cartesian coordinates of nodes
       iStart,iMax,jStart,jMax,kStart,kMax,&!(in) FACE indexes. 
       FaceAreaJ_DF)                        !(out)Face area vectors
    !-----------------------------------------------------------------------------
    integer,intent(in)::iStart,iMax,jStart,jMax,kStart,kMax
    real,intent(in),dimension(nDim,iStart:iMax+1,jStart:jMax,kStart:kMax+1)::&
         XyzNodes_DIII
    real,intent(out),dimension(nDim,iStart:iMax,jStart:jMax,kStart:kMax)::&
         FaceAreaJ_DF
    !-----------------------------------------------------------------------------
    integer::i,j,k
    do k=kStart,kMax; do j=jStart,jMax; do i=iStart,iMax
       FaceAreaJ_DF(:,i,j,k)=cHalf*&
            cross_product(XyzNodes_DIII(:,i+1,j,k+1)-XyzNodes_DIII(:,i,j  ,k),&
                          XyzNodes_DIII(:,i+1,j  ,k)-XyzNodes_DIII(:,i,j,k+1))
    end do; end do; end do
  end subroutine get_face_area_j
!---------------------------------------------------------------------------------
  subroutine get_face_area_k(&
       XyzNodes_DIII,&                      !(in) Cartesian coordinates of nodes
       iStart,iMax,jStart,jMax,kStart,kMax,&!(in) FACE indexes. 
       FaceAreaK_DF)                        !(out)Face area vectors
    !-----------------------------------------------------------------------------
    integer,intent(in)::iStart,iMax,jStart,jMax,kStart,kMax
    real,intent(in),dimension(nDim,iStart:iMax+1,jStart:jMax+1,kStart:kMax)::&
         XyzNodes_DIII
    real,intent(out),dimension(nDim,iStart:iMax,jStart:jMax,kStart:kMax)::&
         FaceAreaK_DF
    !-----------------------------------------------------------------------------
    integer::i,j,k
    do k=kStart,kMax; do j=jStart,jMax; do i=iStart,iMax
       FaceAreaK_DF(:,i,j,k)=cHalf*&
            cross_product(XyzNodes_DIII(:,i+1,j+1,k)-XyzNodes_DIII(:,i,j  ,k),&
                          XyzNodes_DIII(:,i,j+1  ,k)-XyzNodes_DIII(:,i+1,j,k))
    end do; end do; end do
  end subroutine get_face_area_k
!-------------------------------------------------------------!
  logical function is_axial_geometry()
    is_axial_geometry=index(TypeGeometry,'spherical')  >0.or.&
                      index(TypeGeometry,'cylindrical')>0.or.&
                      index(TypeGeometry,'axial')>0 
  end function is_axial_geometry
end module ModCovariant


