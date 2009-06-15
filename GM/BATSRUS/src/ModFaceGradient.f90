!^CFG COPYRIGHT UM
!==============================================================================
module ModFaceGradient

  use ModSize, ONLY: nDim, nI, nJ, nK

  implicit none
  save

  private ! except

  ! Public methods
  public :: set_block_field
  public :: calc_face_gradient
  public :: calc_face_curl

  ! Jacobian matrix for covariant grid: Dcovariant/Dcartesian
  real :: DcoordDxyz_DDFD(nDim,nDim,1:nI+1,1:nJ+1,1:nK+1,nDim)

contains

  !============================================================================

  subroutine set_block_field(iBlock, nVar, Field1_VG, Field_VG)

    ! correct the ghostcells of the given scalar/vector field on iBlock

    use ModParallel, ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot, BlkNeighborLev, NOBLK

    integer, intent(in) :: iBlock, nVar
    real, intent(inout) :: Field1_VG(nVar,0:nI+1,0:nJ+1,0:nK+1)
    real, intent(inout) :: Field_VG(nVar,-1:nI+2,-1:nJ+2,-1:nK+2)

    real, parameter :: c0 = 0.5, p0 = 1./6., F1 = 1./3.

    integer :: i1, j1, k1, i2, j2, k2, iC, jC, kC, iSide, jSide, kSide
    integer :: iL, iR, jL, jR, kL, kR
    integer :: ip, jp, kp

    logical :: IsEqualLevel_G(0:nI+1,0:nJ+1,0:nK+1)
    !--------------------------------------------------------------------------

    Field1_VG = Field_VG(:,0:nI+1,0:nJ+1,0:nK+1)

    do kSide = -1,1; do jSide = -1,1; do iSide = -1,1
       if(iSide==0)then
          iL = 1; iR = nI
       elseif(iSide==1)then
          iL = nI+1; iR = iL
       else
          iL = 0; iR = iL
       end if
       if(jSide==0)then
          jL = 1; jR = nJ
       elseif(jSide==1)then
          jL = nJ+1; jR = jL
       else
          jL = 0; jR = jL
       end if
       if(kSide==0)then
          kL = 1; kR = nK
       elseif(kSide==1)then
          kL = nK+1; kR = kL
       else
          kL = 0; kR = kL
       end if
       if( BlkNeighborLev(iSide, jSide, kSide, iBlock) == 0 )then
          IsEqualLevel_G(iL:iR,jL:jR,kL:kR) = .true.
       else
          IsEqualLevel_G(iL:iR,jL:jR,kL:kR) = .false.
       end if
    end do; end do; end do

    ! Do six faces
    if(NeiLeast(iBlock) == 1)then
       do k1=1, nK, 2; do j1=1, nJ, 2; do k2 = k1,k1+1; do j2 = j1,j1+1
          jp = 3*j2 - 2*j1 -1
          kp = 3*k2 - 2*k1 -1
          if(IsEqualLevel_G(0,jp,kp))then
             Field_VG(:,0,j2,k2) = c0*Field1_VG(:,0,j2,k2) &
                  + 0.25*Field1_VG(:,0,jp,kp) + 0.25*Field_VG(:,1,j2,k2)
          else
             Field_VG(:,0,j2,k2) = c0*Field1_VG(:,0,j2,k2) &
                  + p0*Field1_VG(:,0,jp,kp) + F1*Field_VG(:,1,j2,k2)
          end if
       end do; end do; end do; end do
    end if

    if(NeiLwest(iBlock) == 1)then
       do k1=1, nK, 2; do j1=1, nJ, 2; do k2 = k1,k1+1; do j2 = j1,j1+1
          jp = 3*j2 - 2*j1 -1
          kp = 3*k2 - 2*k1 -1
          if(IsEqualLevel_G(nI+1,jp,kp))then
             Field_VG(:,nI+1,j2,k2) = c0*Field1_VG(:,nI+1,j2,k2) &
                  + 0.25*Field1_VG(:,nI+1,jp,kp) + 0.25*Field_VG(:,nI,j2,k2)
          else
             Field_VG(:,nI+1,j2,k2) = c0*Field1_VG(:,nI+1,j2,k2) &
                  + p0*Field1_VG(:,nI+1,jp,kp) + F1*Field_VG(:,nI,j2,k2)
          end if
       end do; end do; end do; end do
    end if

    if(NeiLsouth(iBlock) == 1)then
       do k1=1, nK, 2; do i1=1, nI, 2; do k2 = k1,k1+1; do i2 = i1,i1+1
          ip = 3*i2 - 2*i1 -1
          kp = 3*k2 - 2*k1 -1
          if(IsEqualLevel_G(ip,0,kp))then
             Field_VG(:,i2,0,k2) = c0*Field1_VG(:,i2,0,k2) &
                  + 0.25*Field1_VG(:,ip,0,kp) + 0.25*Field_VG(:,i2,1,k2)
          else
             Field_VG(:,i2,0,k2) = c0*Field1_VG(:,i2,0,k2) &
                  + p0*Field1_VG(:,ip,0,kp) + F1*Field_VG(:,i2,1,k2)
          end if
       end do; end do; end do; end do
    end if

    if(NeiLnorth(iBlock) == 1)then
       do k1=1, nK, 2; do i1=1, nI, 2; do k2 = k1,k1+1; do i2 = i1,i1+1
          ip = 3*i2 - 2*i1 -1
          kp = 3*k2 - 2*k1 -1
          if(IsEqualLevel_G(ip,nJ+1,kp))then
             Field_VG(:,i2,nJ+1,k2) = c0*Field1_VG(:,i2,nJ+1,k2) &
                  + 0.25*Field1_VG(:,ip,nJ+1,kp) + 0.25*Field_VG(:,i2,nJ,k2)
          else
             Field_VG(:,i2,nJ+1,k2) = c0*Field1_VG(:,i2,nJ+1,k2) &
                  + p0*Field1_VG(:,ip,nJ+1,kp) + F1*Field_VG(:,i2,nJ,k2)
          end if
       end do; end do; end do; end do
    end if

    if(NeiLbot(iBlock) == 1)then
       do j1=1, nJ, 2; do i1=1, nI, 2; do j2 = j1,j1+1; do i2 = i1,i1+1
          ip = 3*i2 - 2*i1 -1
          jp = 3*j2 - 2*j1 -1
          if(IsEqualLevel_G(ip,jp,0))then
             Field_VG(:,i2,j2,0) = c0*Field1_VG(:,i2,j2,0) &
                  + 0.25*Field1_VG(:,ip,jp,0) + 0.25*Field_VG(:,i2,j2,1)
          else
             Field_VG(:,i2,j2,0) = c0*Field1_VG(:,i2,j2,0) &
                  + p0*Field1_VG(:,ip,jp,0) + F1*Field_VG(:,i2,j2,1)
          end if
       end do; end do; end do; end do
    end if

    if(NeiLtop(iBlock) == 1)then
       do j1=1, nJ, 2; do i1=1, nI, 2; do j2 = j1,j1+1; do i2 = i1,i1+1
          ip = 3*i2 - 2*i1 -1
          jp = 3*j2 - 2*j1 -1
          if(IsEqualLevel_G(ip,jp,nK+1))then
             Field_VG(:,i2,j2,nK+1) = c0*Field1_VG(:,i2,j2,nK+1) &
                  + 0.25*Field1_VG(:,ip,jp,nK+1) + 0.25*Field_VG(:,i2,j2,nK)
          else
             Field_VG(:,i2,j2,nK+1) = c0*Field1_VG(:,i2,j2,nK+1) &
                  + p0*Field1_VG(:,ip,jp,nK+1) + F1*Field_VG(:,i2,j2,nK)
          end if
       end do; end do; end do; end do
    end if

    ! Do 12 edges
    ! 4 X edges
    do kSide = -1,1,2; do jSide = -1,1,2
       if(  BlkNeighborLev(0, jSide, kSide, iBlock) /= 1 .and. .not. ( &
            BlkNeighborLev(0, jSide, kSide, iBlock) == NOBLK .and. ( &
            BlkNeighborLev(0, jSide, 0, iBlock) == 1 .or. &
            BlkNeighborLev(0, 0, kSide, iBlock) == 1))) CYCLE

       j1=1; if(jSide==1) j1=nJ; jC = j1+jSide
       k1=1; if(kSide==1) k1=nK; kC = k1+kSide
       do i1 = 1,nI,2; do i2 = i1, i1+1
          ip = 3*i2 - 2*i1 -1
          if(IsEqualLevel_G(ip,jC,kC))then
             Field_VG(:,i2,jC,kC) = c0*Field1_VG(:,i2,jC,kC) &
                  + 0.25*Field1_VG(:,ip,jC,kC) + 0.25*Field_VG(:,i2,j1,k1)
          else
             Field_VG(:,i2,jC,kC) = c0*Field1_VG(:,i2,jC,kC) &
                  + p0*Field1_VG(:,ip,jC,kC) + F1*Field_VG(:,i2,j1,k1)
          end if
       end do; end do
    end do; end do
    ! 4 Y edges
    do kSide = -1,1,2; do iSide = -1,1,2
       if(  BlkNeighborLev(iSide, 0, kSide, iBlock) /= 1 .and. .not. ( &
            BlkNeighborLev(iSide, 0, kSide, iBlock) == NOBLK .and. ( &
            BlkNeighborLev(iSide, 0, 0, iBlock) == 1 .or. &
            BlkNeighborLev(0, 0, kSide, iBlock) == 1))) CYCLE

       i1=1; if(iSide==1) i1=nI; iC = i1+iSide
       k1=1; if(kSide==1) k1=nK; kC = k1+kSide
       do j1 = 1, nJ, 2; do j2 = j1, j1+1
          jp = 3*j2 - 2*j1 -1
          if(IsEqualLevel_G(iC,jp,kC))then
             Field_VG(:,iC,j2,kC) = c0*Field1_VG(:,iC,j2,kC) &
                  + 0.25*Field1_VG(:,iC,jp,kC) + 0.25*Field_VG(:,i1,j2,k1)
          else
             Field_VG(:,iC,j2,kC) = c0*Field1_VG(:,iC,j2,kC) &
                  + p0*Field1_VG(:,iC,jp,kC) + F1*Field_VG(:,i1,j2,k1)
          end if
       end do; end do
    end do; end do
    ! 4 Z edges
    do jSide = -1,1,2; do iSide = -1,1,2
       if(  BlkNeighborLev(iSide, jSide, 0, iBlock) /= 1.and. .not. ( &
            BlkNeighborLev(iSide, jSide, 0, iBlock) == NOBLK .and. ( &
            BlkNeighborLev(iSide, 0, 0, iBlock) == 1 .or. &
            BlkNeighborLev(0, jSide, 0, iBlock) == 1))) CYCLE

       i1=1; if(iSide==1) i1=nI; iC = i1+iSide
       j1=1; if(jSide==1) j1=nJ; jC = j1+jSide
       do k1 = 1, nK, 2 ; do k2 = k1, k1 + 1
          kp = 3*k2 - 2*k1 -1
          if(IsEqualLevel_G(iC,jC,kp))then
             Field_VG(:,iC,jC,k2) = c0*Field1_VG(:,iC,jC,k2) &
                  + 0.25*Field1_VG(:,iC,jC,kp) + 0.25*Field_VG(:,i1,j1,k2)
          else
             Field_VG(:,iC,jC,k2) = c0*Field1_VG(:,iC,jC,k2) &
                  + p0*Field1_VG(:,iC,jC,kp) + F1*Field_VG(:,i1,j1,k2)
          end if
       end do; end do         
    end do; end do

  end subroutine set_block_field

  !============================================================================

  subroutine set_block_jacobian_face(iBlock)

    use ModMain, ONLY: x_, y_, z_
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, Dx_Blk, Dy_Blk, Dz_Blk
    use ModCoordTransform, ONLY: inverse_matrix

    integer, intent(in):: iBlock

    ! Dxyz/Dcoord matrix for one cell
    real:: DxyzDcoord_DD(nDim,nDim)

    ! Transverse gradients
    real:: TransGrad_DDG(nDim,nDim,-1:nI+2,-1:nJ+2,-1:nK+2)

    ! Cell center coordinates for this block
    real:: Xyz_DG(nDim,-1:nI+2,-1:nJ+2,-1:nK+2)

    ! Inverse of cell size
    real :: InvDx, InvDy, InvDz

    ! Indexes
    integer:: i, j, k

    !coeff of Ui+2 and Ui+1 to get normal derivative
    real, parameter:: fp2 = -1./24.0, fp1 = 9.0/8.0 
    !coeff of Ui+2 and Ui+1 for transverse derivatives
    real, parameter:: dp2 = -1./12.0, dp1 = 2.0/3.0 
    !coeff to average transverse derivatives
    real, parameter:: ap2 = -1./16.0, ap1 = 9.0/16. 

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub='set_block_jacobian_face'
    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    ! Calculate the dCovariant/dCartesian matrix

    InvDx = 1.0/Dx_Blk(iBlock)
    InvDy = 1.0/Dy_Blk(iBlock)
    InvDz = 1.0/Dz_Blk(iBlock)

    Xyz_DG(x_,:,:,:) = x_BLK(:,:,:,iBlock)
    Xyz_DG(y_,:,:,:) = y_BLK(:,:,:,iBlock)
    Xyz_DG(z_,:,:,:) = z_BLK(:,:,:,iBlock)

    do k=-1,nK+2; do j=-1,nJ+2; do i=1,nI
       TransGrad_DDG(:,1,i,j,k)=  &
            ( dp1* (Xyz_DG(:,i+1,j,k) - Xyz_DG(:,i-1,j,k)) &
            + dp2* (Xyz_DG(:,i+2,j,k) - Xyz_DG(:,i-2,j,k)))
    end do; end do; end do

    do k=-1,nK+2; do j=1,nJ; do i=-1,nI+2
       TransGrad_DDG(:,2,i,j,k)=  &
            ( dp1* (Xyz_DG(:,i,j+1,k) - Xyz_DG(:,i,j-1,k)) &
            + dp2* (Xyz_DG(:,i,j+2,k) - Xyz_DG(:,i,j-2,k)))
    end do; end do; end do

    do k=1,nK; do j=-1,nJ+2; do i=-1,nI+2
       TransGrad_DDG(:,3,i,j,k)=  &
            ( dp1* (Xyz_DG(:,i,j,k+1) - Xyz_DG(:,i,j,k-1)) &
            + dp2* (Xyz_DG(:,i,j,k+2) - Xyz_DG(:,i,j,k-2)))
    end do; end do; end do

    ! coord1 face
    do k=1,nK; do j=1,nJ; do i=1,nI+1
       ! DxyzDcoord along coord1 face
       DxyzDcoord_DD(:,1) = InvDx* &
            (  fp1*(Xyz_DG(:,i  ,j,k) - Xyz_DG(:,i-1,j,k)) &
            +  fp2*(Xyz_DG(:,i+1,j,k) - Xyz_DG(:,i-2,j,k)))
       DxyzDcoord_DD(:,2) = InvDy* &
            ( ap1*( TransGrad_DDG(:,2,i  ,j,k) + TransGrad_DDG(:,2,i-1,j,k)) &
            + ap2*( TransGrad_DDG(:,2,i+1,j,k) + TransGrad_DDG(:,2,i-2,j,k)))
       DxyzDcoord_DD(:,3) = InvDz* &
            ( ap1*( TransGrad_DDG(:,3,i  ,j,k) + TransGrad_DDG(:,3,i-1,j,k)) &
            + ap2*( TransGrad_DDG(:,3,i+1,j,k) + TransGrad_DDG(:,3,i-2,j,k)))
       DcoordDxyz_DDFD(:,:,i,j,k,1) = &
            inverse_matrix(DxyzDcoord_DD, DoIgnoreSingular=.true.)
    end do; end do; end do

    ! coord2 face
    do k=1,nK; do j=1,nJ+1; do i=1,nI
       ! DxyzDcoord along coord2 face
       DxyzDcoord_DD(:,1) = InvDx* &
            ( ap1*( TransGrad_DDG(:,1,i,j  ,k) + TransGrad_DDG(:,1,i,j-1,k)) &
            + ap2*( TransGrad_DDG(:,1,i,j+1,k) + TransGrad_DDG(:,1,i,j-2,k)))
       DxyzDcoord_DD(:,2) = InvDy* &
            (  fp1*(Xyz_DG(:,i,j  ,k) - Xyz_DG(:,i,j-1,k)) &
            +  fp2*(Xyz_DG(:,i,j+1,k) - Xyz_DG(:,i,j-2,k)))
       DxyzDcoord_DD(:,3) = InvDz* &
            ( ap1*( TransGrad_DDG(:,3,i,j  ,k) + TransGrad_DDG(:,3,i,j-1,k)) &
            + ap2*( TransGrad_DDG(:,3,i,j+1,k) + TransGrad_DDG(:,3,i,j-2,k)))

       DcoordDxyz_DDFD(:,:,i,j,k,2) = &
            inverse_matrix(DxyzDcoord_DD, DoIgnoreSingular=.true.)
    end do; end do; end do

    ! coord3 face
    do k=1,nK+1; do j=1,nJ; do i=1,nI
       ! DxyzDcoord along coord3 face
       DxyzDcoord_DD(:,1) = InvDx* &
            ( ap1*( TransGrad_DDG(:,1,i,j,k  ) + TransGrad_DDG(:,1,i,j,k-1)) &
            + ap2*( TransGrad_DDG(:,1,i,j,k+1) + TransGrad_DDG(:,1,i,j,k-2)))
       DxyzDcoord_DD(:,2) = InvDy* &
            ( ap1*( TransGrad_DDG(:,2,i,j,k  ) + TransGrad_DDG(:,2,i,j,k-1)) &
            + ap2*( TransGrad_DDG(:,2,i,j,k+1) + TransGrad_DDG(:,2,i,j,k-2)))
       DxyzDcoord_DD(:,3) = InvDz* &
            (  fp1*(Xyz_DG(:,i,j,k  ) - Xyz_DG(:,i,j,k-1)) &
            +  fp2*(Xyz_DG(:,i,j,k+1) - Xyz_DG(:,i,j,k-2)))

       DcoordDxyz_DDFD(:,:,i,j,k,3) = &
            inverse_matrix(DxyzDcoord_DD, DoIgnoreSingular=.true.)
    end do; end do; end do

  end subroutine set_block_jacobian_face

  !============================================================================

  subroutine calc_face_gradient(iDir, i, j, k, iBlock, Scalar_G, IsNewBlock, &
       FaceGrad_D) 

    ! calculate the cell face gradient of Scalar_G

    use ModAdvance,    ONLY: State_VGB
    use ModCovariant,  ONLY: UseCovariant
    use ModGeometry,   ONLY: Dx_Blk, Dy_Blk, Dz_Blk
    use ModMain,       ONLY: x_, y_, z_
    use ModParallel,   ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot, BlkNeighborLev

    integer, intent(in) :: iDir, i, j, k, iBlock
    real, intent(inout) :: Scalar_G(-1:nI+2,-1:nJ+2,-1:nK+2)
    logical, intent(inout) :: IsNewBlock
    real, intent(out) :: FaceGrad_D(3)

    integer :: iL, iR, jL, jR, kL, kR
    real :: Ax, Bx, Cx, Ay, By, Cy, Az, Bz, Cz
    Real :: InvDx, InvDy, InvDz
    real :: Scalar1_G(0:nI+1,0:nJ+1,0:nK+1)
    !--------------------------------------------------------------------------
    InvDx = 1.0/Dx_Blk(iBlock)
    InvDy = 1.0/Dy_Blk(iBlock)
    InvDz = 1.0/Dz_Blk(iBlock)

    if(IsNewBlock)then
       call set_block_field(iBlock, 1, Scalar1_G, Scalar_G)
       if(UseCovariant) call set_block_jacobian_face(iBlock)

       IsNewBlock = .false.
    end if

    ! Central difference with averaging in orthogonal direction
    iR = i+1; iL = i-1; 
    jR = j+1; jL = j-1; 
    kR = k+1; kL = k-1; 

    Ax = -0.25*InvDx; Bx = 0.0; Cx = +0.25*InvDx
    Ay = -0.25*InvDy; By = 0.0; Cy = +0.25*InvDy
    Az = -0.25*InvDz; Bz = 0.0; Cz = +0.25*InvDz

    if(i==1)then
       if(NeiLeast(iBlock)==-1 &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev(-1,-1, 0,iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev(-1, 1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev(-1, 0,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev(-1, 0, 1,iBlock)==-1)) &
            )then
          iL = i+1; iR = i+2; Ax=InvDx; Bx=-0.75*InvDx; Cx=-0.25*InvDx
       end if
    elseif(i==nI)then
       if(NeiLwest(iBlock)==-1 &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev( 1,-1, 0,iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev( 1, 1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev( 1, 0,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev( 1, 0, 1,iBlock)==-1)) &
            )then
          iL = i-1; iR = i-2; Ax=-InvDx; Bx=0.75*InvDx; Cx=0.25*InvDx
       end if
    end if

    if(j==1)then
       if(NeiLsouth(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1,-1, 0,iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1,-1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev( 0,-1,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev( 0,-1, 1,iBlock)==-1)) &
            )then
          jL = j+1; jR = j+2; Ay=InvDy; By=-0.75*InvDy; Cy=-0.25*InvDy
       end if
    elseif(j==nJ)then
       if(NeiLnorth(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1, 1, 0,iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1, 1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev( 0, 1,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev( 0, 1, 1,iBlock)==-1))&
            )then
          jL = j-1; jR = j-2; Ay=-InvDy; By=0.75*InvDy; Cy=0.25*InvDy
       end if
    end if

    if(k==1)then
       if(NeiLbot(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1, 0,-1,iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1, 0,-1,iBlock)==-1)) &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev( 0,-1,-1,iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev( 0, 1,-1,iBlock)==-1)) &
            )then
          kL = k+1; kR = k+2; Az=InvDz; Bz=-0.75*InvDz; Cz=-0.25*InvDz
       end if
    elseif(k==nK)then
       if(NeiLtop(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1, 0, 1,iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1, 0, 1,iBlock)==-1)) &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev( 0,-1, 1,iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev( 0, 1, 1,iBlock)==-1)) &
            )then
          kL = k-1; kR = k-2; Az=-InvDz; Bz=0.75*InvDz; Cz=0.25*InvDz
       end if
    end if

    ! Use central difference to get gradient at face
    select case(iDir)
    case(x_)
       FaceGrad_D(x_) = InvDx*(Scalar_G(i,j,k) - Scalar_G(i-1,j,k))
       FaceGrad_D(y_) = &
            + Ay*(Scalar_G(i-1,jL,k) + Scalar_G(i,jL,k)) &
            + By*(Scalar_G(i-1,j ,k) + Scalar_G(i,j ,k)) &
            + Cy*(Scalar_G(i-1,jR,k) + Scalar_G(i,jR,k))
       FaceGrad_D(z_) = &
            + Az*(Scalar_G(i-1,j,kL) + Scalar_G(i,j,kL)) &
            + Bz*(Scalar_G(i-1,j,k ) + Scalar_G(i,j,k )) &
            + Cz*(Scalar_G(i-1,j,kR) + Scalar_G(i,j,kR))
    case(y_)
       FaceGrad_D(x_) = &
            + Ax*(Scalar_G(iL,j-1,k) + Scalar_G(iL,j,k)) &
            + Bx*(Scalar_G(i ,j-1,k) + Scalar_G(i ,j,k)) &
            + Cx*(Scalar_G(iR,j-1,k) + Scalar_G(iR,j,k))
       FaceGrad_D(y_) = InvDy*(Scalar_G(i,j,k) - Scalar_G(i,j-1,k))
       FaceGrad_D(z_) = &
            + Az*(Scalar_G(i,j-1,kL) + Scalar_G(i,j,kL)) &
            + Bz*(Scalar_G(i,j-1,k ) + Scalar_G(i,j,k )) &
            + Cz*(Scalar_G(i,j-1,kR) + Scalar_G(i,j,kR))
    case(z_)
       FaceGrad_D(x_) = &
            + Ax*(Scalar_G(iL,j,k-1) + Scalar_G(iL,j,k)) &
            + Bx*(Scalar_G(i ,j,k-1) + Scalar_G(i ,j,k)) &
            + Cx*(Scalar_G(iR,j,k-1) + Scalar_G(iR,j,k))
       FaceGrad_D(y_) = &
            + Ay*(Scalar_G(i,jL,k-1) + Scalar_G(i,jL,k))  &
            + By*(Scalar_G(i,j ,k-1) + Scalar_G(i,j ,k))  &
            + Cy*(Scalar_G(i,jR,k-1) + Scalar_G(i,jR,k))
       FaceGrad_D(z_) = InvDz*(Scalar_G(i,j,k) - Scalar_G(i,j,k-1))
    case default
       write(*,*)'Error in calc_face_gradient: iDir=',iDir
       call stop_mpi('DEBUG')
    end select

    ! multiply by the coordinate transformation matrix to obtain the
    ! cartesian gradient from the partial derivatives dScalar/dCovariant
    if(UseCovariant) &
         FaceGrad_D = matmul(FaceGrad_D,DcoordDxyz_DDFD(:,:,i,j,k,iDir))

  end subroutine calc_face_gradient

  !============================================================================

  subroutine calc_face_curl(iDir, i, j, k, iBlock, Vector_DG, IsNewBlock, &
       FaceCurl_D)

    use ModMain,      ONLY: x_, y_, z_
    use ModGeometry,  ONLY: Dx_BLK, Dy_BLK, Dz_BLK
    use ModCovariant, ONLY: UseCovariant                
    use ModParallel,  ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot, BlkNeighborLev

    integer, intent(in) :: iDir, i, j, k, iBlock
    real, intent(inout) :: Vector_DG(3,-1:nI+2,-1:nJ+2,-1:nK+2)
    logical, intent(inout) :: IsNewBlock
    real, intent(out)  :: FaceCurl_D(3)

    integer :: iL, iR, jL, jR, kL, kR
    real :: Ax, Bx, Cx, Ay, By, Cy, Az, Bz, Cz
    Real :: InvDx, InvDy, InvDz
    real :: Vector1_DG(3,0:nI+1,0:nJ+1,0:nK+1)
    !--------------------------------------------------------------------------

    InvDx = 1.0/dx_Blk(iBlock)
    InvDy = 1.0/dy_Blk(iBlock)
    InvDz = 1.0/dz_Blk(iBlock)

    if(IsNewBlock)then
       call set_block_field(iBlock, 3, Vector1_DG, Vector_DG)
       if(UseCovariant) call set_block_jacobian_face(iBlock)

       IsNewBlock = .false.
    end if

    ! Central difference with averaging in orthogonal direction
    iR = i+1; iL = i-1; 
    jR = j+1; jL = j-1; 
    kR = k+1; kL = k-1; 

    Ax = -0.25*InvDx; Bx = 0.0; Cx = +0.25*InvDx
    Ay = -0.25*InvDy; By = 0.0; Cy = +0.25*InvDy
    Az = -0.25*InvDz; Bz = 0.0; Cz = +0.25*InvDz

    if(i==1)then
       if(NeiLeast(iBlock)==-1 &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev(-1,-1, 0,iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev(-1, 1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev(-1, 0,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev(-1, 0, 1,iBlock)==-1)) &
            )then
          iL = i+1; iR = i+2; Ax=InvDx; Bx=-0.75*InvDx; Cx=-0.25*InvDx
       end if
    elseif(i==nI)then
       if(NeiLwest(iBlock)==-1 &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev( 1,-1, 0,iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev( 1, 1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev( 1, 0,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev( 1, 0, 1,iBlock)==-1)) &
            )then
          iL = i-1; iR = i-2; Ax=-InvDx; Bx=0.75*InvDx; Cx=0.25*InvDx
       end if
    end if

    if(j==1)then
       if(NeiLsouth(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1,-1, 0,iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1,-1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev( 0,-1,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev( 0,-1, 1,iBlock)==-1)) &
            )then
          jL = j+1; jR = j+2; Ay=InvDy; By=-0.75*InvDy; Cy=-0.25*InvDy
       end if
    elseif(j==nJ)then
       if(NeiLnorth(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1, 1, 0,iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1, 1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev( 0, 1,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev( 0, 1, 1,iBlock)==-1)) &
            )then
          jL = j-1; jR = j-2; Ay=-InvDy; By=0.75*InvDy; Cy=0.25*InvDy
       end if
    end if

    if(k==1)then
       if(NeiLbot(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1, 0,-1,iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1, 0,-1,iBlock)==-1)) &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev( 0,-1,-1,iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev( 0, 1,-1,iBlock)==-1)) &
            )then
          kL = k+1; kR = k+2; Az=InvDz; Bz=-0.75*InvDz; Cz=-0.25*InvDz
       end if
    elseif(k==nK)then
       if(NeiLtop(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1, 0, 1,iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1, 0, 1,iBlock)==-1)) &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev( 0,-1, 1,iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev( 0, 1, 1,iBlock)==-1)) &
            )then
          kL = k-1; kR = k-2; Az=-InvDz; Bz=0.75*InvDz; Cz=0.25*InvDz
       end if
    end if

    if(UseCovariant)then               
       call calc_covariant_curl
    else
       call calc_cartesian_curl
    end if                             

  contains

    !==========================================================================

    subroutine calc_cartesian_curl

      !------------------------------------------------------------------------
      select case(iDir)
      case(x_)
         FaceCurl_D(y_) = &
              -InvDx*(Vector_DG(z_,i,j,k) - Vector_DG(z_,i-1,j,k)) &
              + Az*(Vector_DG(x_,i-1,j,kL) + Vector_DG(x_,i,j,kL)) &
              + Bz*(Vector_DG(x_,i-1,j,k ) + Vector_DG(x_,i,j,k )) &
              + Cz*(Vector_DG(x_,i-1,j,kR) + Vector_DG(x_,i,j,kR))

         FaceCurl_D(z_) = &
              +InvDx*(Vector_DG(y_,i,j,k) - Vector_DG(y_,i-1,j,k)) &
              - Ay*(Vector_DG(x_,i-1,jL,k) + Vector_DG(x_,i,jL,k)) &
              - By*(Vector_DG(x_,i-1,j ,k) + Vector_DG(x_,i,j ,k)) &
              - Cy*(Vector_DG(x_,i-1,jR,k) + Vector_DG(x_,i,jR,k)) 

         FaceCurl_D(x_) = &
              + Ay*(Vector_DG(z_,i-1,jL,k) + Vector_DG(z_,i,jL,k )) &
              + By*(Vector_DG(z_,i-1,j ,k) + Vector_DG(z_,i,j ,k )) &
              + Cy*(Vector_DG(z_,i-1,jR,k) + Vector_DG(z_,i,jR,k )) &
              - Az*(Vector_DG(y_,i-1,j,kL) + Vector_DG(y_,i,j ,kL)) &
              - Bz*(Vector_DG(y_,i-1,j,k ) + Vector_DG(y_,i,j ,k )) &
              - Cz*(Vector_DG(y_,i-1,j,kR) + Vector_DG(y_,i,j ,kR))

      case(y_)
         FaceCurl_D(x_) = &
              +InvDy*(Vector_DG(z_,i,j,k) - Vector_DG(z_,i,j-1,k)) &
              - Az*(Vector_DG(y_,i,j-1,kL) + Vector_DG(y_,i,j,kL)) &
              - Bz*(Vector_DG(y_,i,j-1,k ) + Vector_DG(y_,i,j,k )) &
              - Cz*(Vector_DG(y_,i,j-1,kR) + Vector_DG(y_,i,j,kR))

         FaceCurl_D(z_) = &
              -InvDy*(Vector_DG(x_,i,j,k) - Vector_DG(x_,i,j-1,k)) &
              + Ax*(Vector_DG(y_,iL,j-1,k) + Vector_DG(y_,iL,j,k)) &
              + Bx*(Vector_DG(y_,i ,j-1,k) + Vector_DG(y_,i ,j,k)) &
              + Cx*(Vector_DG(y_,iR,j-1,k) + Vector_DG(y_,iR,j,k))

         FaceCurl_D(y_) = &
              + Az*(Vector_DG(x_,i,j-1,kL) + Vector_DG(x_,i,j,kL)) &
              + Bz*(Vector_DG(x_,i,j-1,k ) + Vector_DG(x_,i,j,k )) &
              + Cz*(Vector_DG(x_,i,j-1,kR) + Vector_DG(x_,i,j,kR)) &
              - Ax*(Vector_DG(z_,iL,j-1,k) + Vector_DG(z_,iL,j,k)) &
              - Bx*(Vector_DG(z_,i ,j-1,k) + Vector_DG(z_,i ,j,k)) &
              - Cx*(Vector_DG(z_,iR,j-1,k) + Vector_DG(z_,iR,j,k))

      case(z_)
         FaceCurl_D(x_) = &
              -InvDz*(Vector_DG(y_,i,j,k) - Vector_DG(y_,i,j,k-1)) & 
              + Ay*(Vector_DG(z_,i,jL,k-1) + Vector_DG(z_,i,jL,k)) &
              + By*(Vector_DG(z_,i,j ,k-1) + Vector_DG(z_,i,j ,k)) &
              + Cy*(Vector_DG(z_,i,jR,k-1) + Vector_DG(z_,i,jR,k))

         FaceCurl_D(y_) = &
              +InvDz*(Vector_DG(x_,i,j,k) - Vector_DG(x_,i,j,k-1)) &
              - Ax*(Vector_DG(z_,iL,j,k-1) + Vector_DG(z_,iL,j,k)) &
              - Bx*(Vector_DG(z_,i ,j,k-1) + Vector_DG(z_,i ,j,k)) &
              - Cx*(Vector_DG(z_,iR,j,k-1) + Vector_DG(z_,iR,j,k))

         FaceCurl_D(z_) = &
              + Ax*(Vector_DG(y_,iL,j,k-1) + Vector_DG(y_,iL,j,k)) &
              + Bx*(Vector_DG(y_,i ,j,k-1) + Vector_DG(y_,i ,j,k)) &
              + Cx*(Vector_DG(y_,iR,j,k-1) + Vector_DG(y_,iR,j,k)) &
              - Ay*(Vector_DG(x_,i,jL,k-1) + Vector_DG(x_,i,jL,k)) &
              - By*(Vector_DG(x_,i,j ,k-1) + Vector_DG(x_,i,j ,k)) &
              - Cy*(Vector_DG(x_,i,jR,k-1) + Vector_DG(x_,i,jR,k))

      case default
         write(*,*)'Error in calc_cartesian_curl: iDir=',iDir
         call stop_mpi('DEBUG')
      end select

    end subroutine calc_cartesian_curl

    !==========================================================================

    subroutine calc_covariant_curl

      real :: DvectorDcoord_DD(nDim,nDim)
      !------------------------------------------------------------------------

      ! Calculate the partial derivatives dVector/dCovariant
      select case(iDir)
      case(x_)
         DvectorDcoord_DD(:,1) = &
              InvDx*(Vector_DG(:,i,j,k) - Vector_DG(:,i-1,j,k))
         DvectorDcoord_DD(:,2) = &
              + Ay*(Vector_DG(:,i-1,jL,k) + Vector_DG(:,i,jL,k)) &
              + By*(Vector_DG(:,i-1,j ,k) + Vector_DG(:,i,j ,k)) &
              + Cy*(Vector_DG(:,i-1,jR,k) + Vector_DG(:,i,jR,k))
         DvectorDcoord_DD(:,3) = &
              + Az*(Vector_DG(:,i-1,j,kL) + Vector_DG(:,i,j,kL)) &
              + Bz*(Vector_DG(:,i-1,j,k ) + Vector_DG(:,i,j,k )) &
              + Cz*(Vector_DG(:,i-1,j,kR) + Vector_DG(:,i,j,kR))
         
      case(y_)
         DvectorDcoord_DD(:,1) = &
              + Ax*(Vector_DG(:,iL,j-1,k) + Vector_DG(:,iL,j,k)) &
              + Bx*(Vector_DG(:,i ,j-1,k) + Vector_DG(:,i ,j,k)) &
              + Cx*(Vector_DG(:,iR,j-1,k) + Vector_DG(:,iR,j,k))
         DvectorDcoord_DD(:,2) = &
              InvDy*(Vector_DG(:,i,j,k) - Vector_DG(:,i,j-1,k))
         DvectorDcoord_DD(:,3) = &
              + Az*(Vector_DG(:,i,j-1,kL) + Vector_DG(:,i,j,kL)) &
              + Bz*(Vector_DG(:,i,j-1,k ) + Vector_DG(:,i,j,k )) &
              + Cz*(Vector_DG(:,i,j-1,kR) + Vector_DG(:,i,j,kR))

      case(z_)
         DvectorDcoord_DD(:,1) = &
              + Ax*(Vector_DG(:,iL,j,k-1) + Vector_DG(:,iL,j,k)) &
              + Bx*(Vector_DG(:,i ,j,k-1) + Vector_DG(:,i ,j,k)) &
              + Cx*(Vector_DG(:,iR,j,k-1) + Vector_DG(:,iR,j,k))
         DvectorDcoord_DD(:,2) = &
              + Ay*(Vector_DG(:,i,jL,k-1) + Vector_DG(:,i,jL,k)) &
              + By*(Vector_DG(:,i,j ,k-1) + Vector_DG(:,i,j ,k)) &
              + Cy*(Vector_DG(:,i,jR,k-1) + Vector_DG(:,i,jR,k))
         DvectorDcoord_DD(:,3) = &
              InvDz*(Vector_DG(:,i,j,k) - Vector_DG(:,i,j,k-1))
      end select

      ! Curl_x = Dvector_z/Dy - Dvector_y/Dz
      FaceCurl_D(x_) = &
           + sum(DvectorDcoord_DD(z_,:)*DcoordDxyz_DDFD(:,y_,i,j,k,iDir)) &
           - sum(DvectorDcoord_DD(y_,:)*DcoordDxyz_DDFD(:,z_,i,j,k,iDir))

      ! Curl_y = Dvector_x/Dz - Dvector_z/Dx
      FaceCurl_D(y_) = &
           + sum(DvectorDcoord_DD(x_,:)*DcoordDxyz_DDFD(:,z_,i,j,k,iDir)) &
           - sum(DvectorDcoord_DD(z_,:)*DcoordDxyz_DDFD(:,x_,i,j,k,iDir))

      ! Curl_z = Dvector_y/Dx - Dvector_x/Dy
      FaceCurl_D(z_) = &
           + sum(DvectorDcoord_DD(y_,:)*DcoordDxyz_DDFD(:,x_,i,j,k,iDir)) &
           - sum(DvectorDcoord_DD(x_,:)*DcoordDxyz_DDFD(:,y_,i,j,k,iDir))

    end subroutine calc_covariant_curl

  end subroutine calc_face_curl

end module ModFaceGradient
