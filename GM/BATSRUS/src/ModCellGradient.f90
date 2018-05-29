module ModCellGradient

  use BATL_lib, ONLY: &
       test_start, test_stop

  ! Calculate cell centered gradient, divergence and curl

  use BATL_lib, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nI, nJ, nK, &
       CellSize_DB,  CellFace_DB, FaceNormal_DDFB, CellVolume_GB, &
       IsCartesian, IsCartesianGrid, Unused_B, nBlock, MaxBlock,&
       nDim, jDim_, kDim_, x_, y_, z_, Dim1_, Dim2_, Dim3_
  use ModGeometry, ONLY: body_blk, true_cell

  implicit none

  SAVE
  private ! except

  public:: calc_divergence ! calculate divergence
  public:: calc_gradient   ! calculate gradient
  public:: calc_gradient_ghost  ! set gradient in ghost cells too
  public:: calc_cell_curl_ghost ! calculate curl with 1 layer of ghost cells
  
  real, public, allocatable :: GradVar_DGB(:,:,:,:,:)

  ! Local variables -------------

  integer, allocatable:: iTrue_G(:,:,:)

  interface calc_gradient
     module procedure calc_gradient1, calc_gradient3
  end interface

contains
  !============================================================================

  subroutine calc_divergence(iBlock, Var_DG, nG, Div_G, UseBodyCellIn)

    ! Calculate divergence of Var_DG and return it in Div_G
    ! Only the physical cells are calculated but Div_G can have
    ! nG ghost cells.
    ! Body (false) cells are ignored unless UseBodyCellIn is set to true.

    use ModGeometry, ONLY: body_blk, true_cell

    ! Calculate divergence Div_G of Var_DG on a Cartesian grid

    integer, intent(in):: iBlock
    real, intent(in) :: Var_DG(nDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    integer, intent(in):: nG ! number of ghost cells in Div_G
    real, intent(inout):: &  ! preserve ghost cell values!
         Div_G(1-nG:nI+nG,1-nG*jDim_:nJ+nG*jDim_,1-nG*kDim_:nK+nG*kDim_)
    logical, intent(in), optional:: UseBodyCellIn

    logical:: UseBodyCell

    real:: InvDxHalf, InvDyHalf, InvDzHalf, InvDx, InvDy, InvDz
    real:: VarR_D(nDim), VarL_D(nDIm)
    integer :: i, j, k, iL, iR, jL, jR, kL, kR

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_divergence'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    UseBodyCell = .false.
    if(present(UseBodyCellIn)) UseBodyCell = UseBodyCellIn

    if(UseBodyCell .or. .not. body_blk(iBlock)) then
       if(IsCartesian)then

          ! Simple central differencing

          InvDxHalf = 0.5/CellSize_DB(1,iBlock)
          InvDyHalf = 0.5/CellSize_DB(2,iBlock)
          InvDzHalf = 0.5/CellSize_DB(3,iBlock)

          do k=1,nK; do j=1,nJ; do i=1,nI
             Div_G(i,j,k) = &
                  InvDxHalf*(Var_DG(Dim1_,i+1,j,k) - Var_DG(Dim1_,i-1,j,k))
             if(nJ > 1) Div_G(i,j,k) = Div_G(i,j,k) + &
                  InvDyHalf*(Var_DG(Dim2_,i,j+1,k) - Var_DG(Dim2_,i,j-1,k))
             if(nK > 1) Div_G(i,j,k) = Div_G(i,j,k) + &
                  InvDzHalf*(Var_DG(Dim3_,i,j,k+1) - Var_DG(Dim3_,i,j,k-1))
          end do; end do; end do

       else
          ! div(Var) = Sum(AreaNormal.VarFace_D)/Volume for the 2*nDim faces
          ! VarFace_D = (VarCenter_D + VarNeighbor_D)/2
          ! The contributions from the cell center cancel out
          do k=1,nK; do j=1,nJ; do i=1,nI
             Div_G(i,j,k) = &
                  +sum(FaceNormal_DDFB(:,1,i+1,j,k,iBlock)*Var_DG(:,i+1,j,k)) &
                  -sum(FaceNormal_DDFB(:,1,i  ,j,k,iBlock)*Var_DG(:,i-1,j,k))
             if(nJ > 1) Div_G(i,j,k) = Div_G(i,j,k) &
                  +sum(FaceNormal_DDFB(:,2,i,j+1,k,iBlock)*Var_DG(:,i,j+1,k)) &
                  -sum(FaceNormal_DDFB(:,2,i,j  ,k,iBlock)*Var_DG(:,i,j-1,k))
             if(nK > 1) Div_G(i,j,k) = Div_G(i,j,k) &
                  +sum(FaceNormal_DDFB(:,3,i,j,k+1,iBlock)*Var_DG(:,i,j,k+1)) &
                  -sum(FaceNormal_DDFB(:,3,i,j,k  ,iBlock)*Var_DG(:,i,j,k-1))
             Div_G(i,j,k) = Div_G(i,j,k)*0.5/CellVolume_GB(i,j,k,iBlock)
          end do; end do; end do

       end if
    else
       ! One sided differences if the neighbor is inside a body

       ! Set iTrue_G to 1 in used cells and 0 in non-used cells
       if(.not.allocated(iTrue_G)) allocate(iTrue_G(0:nI+1,0:nJ+1,0:nK+1))
       where(true_cell(0:nI+1,0:nJ+1,0:nK+1,iBlock))
          iTrue_G = 1
       elsewhere
          iTrue_G = 0
       end where

       if(IsCartesian)then
          ! Where .not.true_cell, all the gradients are zero
          ! In true_cell the input to gradient from the face neighbor
          ! is ignored, if the face neighbor is .not.true_cell, the input
          ! from the opposite cell is doubled in this case

          InvDx = 1/CellSize_DB(1,iBlock)
          InvDy = 1/CellSize_DB(2,iBlock)
          InvDz = 1/CellSize_DB(3,iBlock)

          do k=1,nK; do j=1,nJ; do i=1,nI
             Div_G(i,j,k) = 0.0
             if(.not.true_cell(i,j,k,iBlock)) CYCLE

             ! Shift to left and right if the neighbor is a true cell
             iL = i - iTrue_G(i-1,j,k); iR = i + iTrue_G(i+1,j,k)

             ! If at least one neighbor is true, calculate derivative
             if(iL /= iR) Div_G(i,j,k) = &
                  InvDx*(Var_DG(Dim1_,iR,j,k) - Var_DG(Dim1_,iL,j,k))/(iR - iL)

             if(nJ == 1) CYCLE
             jL = j - iTrue_G(i,j-1,k); jR = j + iTrue_G(i,j+1,k)
             if(jL /= jR)Div_G(i,j,k) = Div_G(i,j,k) + &
                  InvDy*(Var_DG(Dim2_,i,jR,k) - Var_DG(Dim2_,i,jL,k))/(jR - jL)

             if(nK == 1) CYCLE
             kL = k - iTrue_G(i,j,k-1); kR = k + iTrue_G(i,j,k+1)
             if(kL /= kR)Div_G(i,j,k) = Div_G(i,j,k) + &
                  InvDz*(Var_DG(Dim3_,i,j,kR) - Var_DG(Dim3_,i,j,kR))/(kR - kL)

          end do; end do; end do
       else

          ! For false cells set div(Var) = 0. For true cells
          ! div(Var) = Sum(AreaNormal.VarFace_D)/Volume for the 2*nDim faces
          ! where VarFace_D = (VarNeighbor_D + VarCenter_D)/2
          ! Since sum(AreaNormal) = 0, the VarCenter_D contributions vanish.
          ! For false neighbor on one side extrapolate from the other side.
          ! For false neighbors on both sides set VarNeigbbor_D = VarCenter_D.

          do k=1,nK; do j=1,nJ; do i=1,nI
             Div_G(i,j,k) = 0.0
             if(.not.true_cell(i,j,k,iBlock)) CYCLE

             ! Shift to left and right if the neighbor is a true cell
             iL = i - iTrue_G(i-1,j,k); iR = i + iTrue_G(i+1,j,k)
             ! Set neighbor values
             VarL_D = Var_DG(:,iL,j,k); VarR_D = Var_DG(:,iR,j,k)
             ! Extrapolate from other side for false neighbors
             if(iL == i .and. iR > i) VarL_D = 2*VarL_D - VarR_D
             if(iR == i .and. iL < i) VarR_D = 2*VarR_D - VarL_D
             Div_G(i,j,k) = sum(FaceNormal_DDFB(:,1,i+1,j,k,iBlock)*VarR_D) &
                  -         sum(FaceNormal_DDFB(:,1,i  ,j,k,iBlock)*VarL_D)

             if(nJ > 1)then
                jL = j - iTrue_G(i,j-1,k); jR = j + iTrue_G(i,j+1,k)
                VarL_D = Var_DG(:,i,jL,k); VarR_D = Var_DG(:,i,jR,k)
                if(jL == j .and. jR > j) VarL_D = 2*VarL_D - VarR_D
                if(jR == j .and. jL < j) VarR_D = 2*VarR_D - VarL_D
                Div_G(i,j,k) = Div_G(i,j,k) &
                     + sum(FaceNormal_DDFB(:,2,i,j+1,k,iBlock)*VarR_D) &
                     - sum(FaceNormal_DDFB(:,2,i,j  ,k,iBlock)*VarL_D)
             end if
             if(nK > 1)then
                kL = k - iTrue_G(i,j,k-1); kR = k + iTrue_G(i,j,k+1)
                VarL_D = Var_DG(:,i,j,kL); VarR_D = Var_DG(:,i,j,kR)
                if(kL == k .and. kR > k) VarL_D = 2*VarL_D - VarR_D
                if(kR == k .and. kL < k) VarR_D = 2*VarR_D - VarL_D
                Div_G(i,j,k) = Div_G(i,j,k) &
                     + sum(FaceNormal_DDFB(:,3,i,j,k+1,iBlock)*VarR_D) &
                     - sum(FaceNormal_DDFB(:,3,i,j,k  ,iBlock)*VarL_D)
             end if
             Div_G(i,j,k) = Div_G(i,j,k)*0.5/CellVolume_GB(i,j,k,iBlock)
          end do; end do; end do
       end if
    end if

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine calc_divergence
  !============================================================================

  subroutine calc_gradient1(iBlock, Var_G, nG, Grad_DG, UseBodyCellIn)

    ! Calculate gradient of Var_G and return it in Grad_DG
    ! Only the physical cells are calculated but Grad_DG can have
    ! nG ghost cells.
    ! Body (false) cells are ignored unless UseBodyCellIn is set to true.

    integer, intent(in) :: iBlock
    real,    intent(in) :: Var_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    integer, intent(in) :: nG  ! number of ghost cells in Grad_DG
    real,    intent(inout):: & ! preserve ghost cell values
         Grad_DG(nDim,1-nG:nI+nG,1-nG*jDim_:nJ+nG*jDim_,1-nG*kDim_:nK+nG*kDim_)
    logical, intent(in), optional:: UseBodyCellIn

    logical:: UseBodyCell

    real:: InvDxHalf, InvDyHalf, InvDzHalf, InvDx, InvDy, InvDz, VarL, VarR
    integer :: i, j, k, iL, iR, jL, jR, kL, kR

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_gradient1'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    UseBodyCell = .false.
    if(present(UseBodyCellIn)) UseBodyCell = UseBodyCellIn

    if(UseBodyCell .or. .not. body_blk(iBlock)) then
       if(IsCartesian)then

          ! Simple central differencing

          InvDxHalf = 0.5/CellSize_DB(1,iBlock)
          InvDyHalf = 0.5/CellSize_DB(2,iBlock)
          InvDzHalf = 0.5/CellSize_DB(3,iBlock)

          do k=1,nK; do j=1,nJ; do i=1,nI
             Grad_DG(Dim1_,i,j,k) = &
                  InvDxHalf*(Var_G(i+1,j,k) - Var_G(i-1,j,k))
             if(nJ > 1) Grad_DG(Dim2_,i,j,k) = &
                  InvDyHalf*(Var_G(i,j+1,k) - Var_G(i,j-1,k))
             if(nK > 1) Grad_DG(Dim3_,i,j,k) = &
                  InvDzHalf*(Var_G(i,j,k+1) - Var_G(i,j,k-1))
          end do; end do; end do
       else
          ! grad(Var) = Sum(AreaNormal*VarFace)/Volume for the 2*nDim faces
          ! VarFace = (VarCenter + VarNeighbor)/2
          ! The contributions from the cell center cancel out
          do k=1,nK; do j=1,nJ; do i=1,nI
             Grad_DG(:,i,j,k) = &
                  FaceNormal_DDFB(:,1,i+1,j,k,iBlock)*Var_G(i+1,j,k) - &
                  FaceNormal_DDFB(:,1,i  ,j,k,iBlock)*Var_G(i-1,j,k)
             if(nJ > 1) Grad_DG(:,i,j,k) = Grad_DG(:,i,j,k) + &
                  FaceNormal_DDFB(:,2,i,j+1,k,iBlock)*Var_G(i,j+1,k) - &
                  FaceNormal_DDFB(:,2,i,j  ,k,iBlock)*Var_G(i,j-1,k)
             if(nK > 1) Grad_DG(:,i,j,k) = Grad_DG(:,i,j,k) + &
                  FaceNormal_DDFB(:,3,i,j,k+1,iBlock)*Var_G(i,j,k+1) - &
                  FaceNormal_DDFB(:,3,i,j,k  ,iBlock)*Var_G(i,j,k-1)
             Grad_DG(:,i,j,k) = Grad_DG(:,i,j,k) &
                  *0.5/CellVolume_GB(i,j,k,iBlock)
          end do; end do; end do
       end if
    else
       ! One sided differences if the neighbor is inside a body

       ! Set iTrue_G to 1 in used cells and 0 in non-used cells
       if(.not.allocated(iTrue_G)) allocate(iTrue_G(0:nI+1,0:nJ+1,0:nK+1))
       where(true_cell(0:nI+1,0:nJ+1,0:nK+1,iBlock))
          iTrue_G = 1
       elsewhere
          iTrue_G = 0
       end where

       if(IsCartesian)then
          InvDx = 1/CellSize_DB(1,iBlock)
          InvDy = 1/CellSize_DB(2,iBlock)
          InvDz = 1/CellSize_DB(3,iBlock)

          do k=1,nK; do j=1,nJ; do i=1,nI
             Grad_DG(:,i,j,k) = 0.0
             if(.not.true_cell(i,j,k,iBlock)) CYCLE

             ! Shift to left and right if the neighbor is a true cell
             ! If at least one neighbor is true, calculate derivative
             iL = i - iTrue_G(i-1,j,k); iR = i + iTrue_G(i+1,j,k)
             if(iL /= iR) Grad_DG(Dim1_,i,j,k) = &
                  InvDx*(Var_G(iR,j,k) - Var_G(iL,j,k))/(iR - iL)

             if(nJ == 1) CYCLE
             jL = j - iTrue_G(i,j-1,k); jR = j + iTrue_G(i,j+1,k)
             if(jL /= jR) Grad_DG(Dim2_,i,j,k) = &
                  InvDy*(Var_G(i,jR,k) - Var_G(i,jL,k))/(jR - jL)

             if(nK == 1) CYCLE
             kL = k - iTrue_G(i,j,k-1); kR = k + iTrue_G(i,j,k+1)
             if(kL /= kR) Grad_DG(Dim3_,i,j,k) = &
                  InvDz*(Var_G(i,j,kR) - Var_G(i,j,kR))/(kR - kL)

          end do; end do; end do
       else

          ! For false cells grad(Var) = 0. For true cells
          ! grad(Var) = Sum(AreaNormal*VarFace)/Volume for the 2*nDim faces
          ! where VarFace = (VarNeighbor + VarCenter)/2.
          ! Since sum(AreaNormal) = 0, the VarCenter contributions vanish.
          ! For false neighbor on one side extrapolate from the other side.
          ! For false neighbors on both sides set VarNeigbbor = VarCenter.

          do k=1,nK; do j=1,nJ; do i=1,nI
             Grad_DG(:,i,j,k) = 0.0
             if(.not.true_cell(i,j,k,iBlock)) CYCLE

             ! Shift to left and right if the neighbor is a true cell
             iL = i - iTrue_G(i-1,j,k); iR = i + iTrue_G(i+1,j,k)
             ! Set neighbor values
             VarR = Var_G(iR,j,k); VarL = Var_G(iL,j,k)
             ! Extrapolate from other side for false neighbors
             if(iL == i .and. iR > i) VarL = 2*VarL - VarR
             if(iR == i .and. iL < i) VarR = 2*VarR - VarL
             Grad_DG(:,i,j,k) = FaceNormal_DDFB(:,1,i+1,j,k,iBlock)*VarR &
                  -             FaceNormal_DDFB(:,1,i  ,j,k,iBlock)*VarL

             if(nJ > 1)then
                jL = j - iTrue_G(i,j-1,k); jR = j + iTrue_G(i,j+1,k)
                VarR = Var_G(i,jR,k); VarL = Var_G(i,jL,k)
                if(jL == j .and. jR > j) VarL = 2*VarL - VarR
                if(jR == j .and. jL < j) VarR = 2*VarR - VarL
                Grad_DG(:,i,j,k) = Grad_DG(:,i,j,k) &
                     + FaceNormal_DDFB(:,2,i,j+1,k,iBlock)*VarR &
                     - FaceNormal_DDFB(:,2,i,j  ,k,iBlock)*VarL
             end if
             if(nK > 1)then
                kL = k - iTrue_G(i,j,k-1); kR = k + iTrue_G(i,j,k+1)
                VarR = Var_G(i,j,kR); VarL = Var_G(i,j,kL)
                if(kL == k .and. kR > k) VarL = 2*VarL - VarR
                if(kR == k .and. kL < k) VarR = 2*VarR - VarL
                Grad_DG(:,i,j,k) = Grad_DG(:,i,j,k) &
                     + FaceNormal_DDFB(:,3,i,j,k+1,iBlock)*VarR &
                     - FaceNormal_DDFB(:,3,i,j,k  ,iBlock)*VarL
             end if
             Grad_DG(:,i,j,k) = Grad_DG(:,i,j,k) &
                  *0.5/CellVolume_GB(i,j,k,iBlock)
          end do; end do; end do
       end if
    end if

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine calc_gradient1
  !============================================================================

  subroutine calc_gradient3(iBlock, Var_G, GradX_C, GradY_C, GradZ_C)

    ! Calculate the
    ! This is an interface to cartesian or gencoord_gradient.

    integer, intent(in):: iBlock
    real,    intent(in):: Var_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(out):: GradX_C(nI,nJ,nK), GradY_C(nI,nJ,nK), GradZ_C(nI,nJ,nK)

    real:: OneTrue_G(0:nI+1,0:nJ+1,0:nK+1)
    real:: VInvHalf
    real :: FaceArea_DS(3,6), Difference_S(6)
    integer :: i, j, k
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_gradient3'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    if(IsCartesianGrid)then
       if(.not.body_blk(iBlock)) then
          do k=1,nK; do j=1,nJ; do i=1,nI
             VInvHalf = 0.5/CellVolume_GB(i,j,k,iBlock)

             GradX_C(i,j,k) = CellFace_DB(x_,iBlock)*&
                  (Var_G(i+1,j,k) - Var_G(i-1,j,k))*VInvHalf
             if(nJ == 1)then
                GradY_C(i,j,k) = 0.0
             else
                GradY_C(i,j,k) = CellFace_DB(y_,iBlock)*&
                     (Var_G(i,j+1,k) - Var_G(i,j-1,k))*VInvHalf
             end if
             if(nK == 1)then
                GradZ_C(i,j,k) = 0.0
             else
                GradZ_C(i,j,k) = CellFace_DB(z_,iBlock)*&
                     (Var_G(i,j,k+1) - Var_G(i,j,k-1))*VInvHalf
             end if
          end do; end do; end do
       else
          where(true_cell(0:nI+1,0:nJ+1,0:nK+1,iBlock))
             OneTrue_G = 1.0
          elsewhere
             OneTrue_G = 0.0
          end where
          !
          !\
          ! Where .not.true_cell, all the gradients are zero
          ! In true_cell the input to gradient from the face neighbor
          ! is ignored, if the face neighbor is .not.true_cell, the input
          ! from the opposite cell is doubled in this case
          !/
          !
          do k=1,nK; do j=1,nJ; do i=1,nI
             VInvHalf = 0.5/CellVolume_GB(i,j,k,iBlock)

             GradX_C(i,j,k) = CellSize_DB(x_,iBlock)*&
                  OneTrue_G(i,j,k)*(&
                  (Var_G(i+1,j,k) - Var_G(i,j,k))*&
                  OneTrue_G(i+1,j,k)*&
                  (2.0 - OneTrue_G(i-1,j,k)) + &
                  (Var_G(i,j,k)-Var_G(i-1,j,k))*&
                  OneTrue_G(i-1,j,k)*&
                  (2.0 - OneTrue_G(i+1,j,k)) )*VInvHalf

             if(nJ==1)then
                GradY_C(i,j,k) = 0.0
             else
                GradY_C(i,j,k) = CellSize_DB(y_,iBlock)*&
                     OneTrue_G(i,j,k)*(&
                     (Var_G(i,j+1,k) - Var_G(i,j,k))*&
                     OneTrue_G(i,j+1,k)*&
                     (2.0 - OneTrue_G(i,j-1,k))+&
                     (Var_G(i,j,k)-Var_G(i,j-1,k))*&
                     OneTrue_G(i,j-1,k)*&
                     (2.0 - OneTrue_G(i,j+1,k)) )*VInvHalf
             end if
             if(nK == 1)then
                GradZ_C(i,j,k) = 0.0
             else
                GradZ_C(i,j,k) = CellSize_DB(z_,iBlock)*&
                     OneTrue_G(i,j,k)*(&
                     (Var_G(i,j,k+1)-Var_G(i,j,k))*&
                     OneTrue_G(i,j,k+1)*&
                     (2.0 - OneTrue_G(i,j,k-1))+&
                     (Var_G(i,j,k)-Var_G(i,j,k-1))*&
                     OneTrue_G(i,j,k-1)*&
                     (2.0 - OneTrue_G(i,j,k+1)) )*VInvHalf
             end if
          end do; end do; end do
       end if
    else
       if(.not.body_BLK(iBlock)) then
          do k=1,nK; do j=1,nJ; do i=1,nI
             VInvHalf = 0.5/CellVolume_GB(i,j,k,iBlock)

             FaceArea_DS(:,1:2) = FaceNormal_DDFB(:,1,i:i+1,j,k,iBlock)
             FaceArea_DS(:,3:4) = FaceNormal_DDFB(:,2,i,j:j+1,k,iBlock)
             FaceArea_DS(:,5:6) = FaceNormal_DDFB(:,3,i,j,k:k+1,iBlock)

             Difference_S(1) = -(Var_G(i-1,j,k) + Var_G(i,j,k))
             Difference_S(2) = +(Var_G(i+1,j,k) + Var_G(i,j,k))
             Difference_S(3) = -(Var_G(i,j-1,k) + Var_G(i,j,k))
             Difference_S(4) = +(Var_G(i,j+1,k) + Var_G(i,j,k))
             Difference_S(5) = -(Var_G(i,j,k-1) + Var_G(i,j,k))
             Difference_S(6) = +(Var_G(i,j,k+1) + Var_G(i,j,k))

             GradX_C(i,j,k) = sum(FaceArea_DS(x_,:)*Difference_S)*VInvHalf
             GradY_C(i,j,k) = sum(FaceArea_DS(y_,:)*Difference_S)*VInvHalf
             GradZ_C(i,j,k) = sum(FaceArea_DS(z_,:)*Difference_S)*VInvHalf
          end do; end do; end do
       else
          where(true_cell(0:nI+1,0:nJ+1,0:nK+1,iBlock))
             OneTrue_G = 1.0
          elsewhere
             OneTrue_G = 0.0
          end where
          do k=1,nK;  do j=1,nJ; do i=1,nI
             if(.not.true_cell(i,j,k,iBlock))then
                GradX_C(i,j,k) = 0.0
                GradY_C(i,j,k) = 0.0
                GradZ_C(i,j,k) = 0.0
                CYCLE
             end if

             FaceArea_DS(:,1:2) = FaceNormal_DDFB(:,1,i:i+1,j,k,iBlock)
             FaceArea_DS(:,3:4) = FaceNormal_DDFB(:,2,i,j:j+1,k,iBlock)
             FaceArea_DS(:,5:6) = FaceNormal_DDFB(:,3,i,j,k:k+1,iBlock)

             Difference_S(1) =  OneTrue_G(i-1,j,k)*&
                  (Var_G(i,j,k)-Var_G(i-1,j,k))+&
                  (1.0 - OneTrue_G(i-1,j,k))*&
                  OneTrue_G(i+1,j,k)*&
                  (Var_G(i+1,j,k)-Var_G(i,j,k))

             Difference_S(2) =  OneTrue_G(i+1,j,k)*&
                  (Var_G(i+1,j,k)-Var_G(i,j,k))+&
                  (1.0 - OneTrue_G(i+1,j,k))*&
                  OneTrue_G(i-1,j,k)*&
                  (Var_G(i,j,k)-Var_G(i-1,j,k))

             Difference_S(3)=  OneTrue_G(i,j-1,k)*&
                  (Var_G(i,j,k)-Var_G(i,j-1,k))+&
                  (1.0 - OneTrue_G(i,j-1,k))*&
                  OneTrue_G(i,j+1,k)*&
                  (Var_G(i,j+1,k)-Var_G(i,j,k))

             Difference_S(4)=  OneTrue_G(i,j+1,k)*&
                  (Var_G(i,j+1,k)-Var_G(i,j,k))+&
                  (1.0 - OneTrue_G(i,j+1,k))*&
                  OneTrue_G(i,j-1,k)*&
                  (Var_G(i,j,k)-Var_G(i,j-1,k))

             Difference_S(5)  =  OneTrue_G(i,j,k-1)*&
                  (Var_G(i,j,k)-Var_G(i,j,k-1))+&
                  (1.0 - OneTrue_G(i,j,k-1))*&
                  OneTrue_G(i,j,k+1)*&
                  (Var_G(i,j,k+1)-Var_G(i,j,k))

             Difference_S(6)  =  OneTrue_G(i,j,k+1)*&
                  (Var_G(i,j,k+1)-Var_G(i,j,k))+&
                  (1.0 - OneTrue_G(i,j,k+1))*&
                  OneTrue_G(i,j,k-1)*&
                  (Var_G(i,j,k)-Var_G(i,j,k-1))

             VInvHalf = 0.5/CellVolume_GB(i,j,k,iBlock)

             GradX_C(i,j,k) = sum(FaceArea_DS(x_,:)*Difference_S)*VInvHalf
             GradY_C(i,j,k) = sum(FaceArea_DS(y_,:)*Difference_S)*VInvHalf
             GradZ_C(i,j,k) = sum(FaceArea_DS(z_,:)*Difference_S)*VInvHalf
          end do; end do; end do
       end if
    end if
    call test_stop(NameSub, DoTest, iBlock)
  end subroutine calc_gradient3
  !============================================================================
  subroutine calc_gradient_ghost(Var_CB, GradVarInOut_DGB)

    use ModFaceGradient, ONLY: set_block_field2
    use BATL_lib, ONLY: j0_, nJp1_, k0_, nKp1_, message_pass_cell

    ! Unless GradVarInOut_DGB is provided, the result is stored
    ! in the public array GradVar_DGB, which can be used for interpolation
    ! outside of this module
    ! CAUTION: you may omit GradVarInOut_DGB ONLY(! ) if you can guarantee that
    ! calc_gradient_ghost isn't called again before you use GradVar_DGB

    real, intent(in)    :: Var_CB(1:nI,1:nJ,1:nK,1:MaxBlock)

    real, optional, intent(inout) :: &
         GradVarInOut_DGB(nDim,0:nI+1,j0_:nJp1_,k0_:nKp1_,MaxBlock)

    real:: Var_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)
    real:: VarCopy_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    integer:: iBlock

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_gradient_ghost'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    Var_GB = 0.0; VarCopy_G = 0.0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock))CYCLE
       Var_GB(1:nI,1:nJ,1:nK,iBlock) = Var_CB(:,:,:,iBlock)
    end do

    ! set scalar variables in the ghost cells
    call message_pass_cell(Var_GB, nProlongOrderIn = 1)
    if(.not.present(GradVarInOut_DGB))then
       if(.not. allocated(GradVar_DGB)) &
            allocate(GradVar_DGB(nDim,0:nI+1,j0_:nJp1_,k0_:nKp1_,MaxBlock))
    end if

    do iBlock = 1, nBlock
       if(Unused_B(iBlock))CYCLE
       ! correct ghost cells
       call set_block_field2(iBlock, nVar=1, &
            Field1_VG = VarCopy_G, Field_VG=Var_GB(:,:,:,iBlock))

       ! compute gradient in cell centers
       if(present(GradVarInOut_DGB))then
          call calc_gradient1(iBlock, Var_GB(:,:,:,iBlock), 1, &
               GradVarInOut_DGB(:,:,:,:,iBlock))
       else
          call calc_gradient1(iBlock, Var_GB(:,:,:,iBlock), 1, &
               GradVar_DGB(:,:,:,:,iBlock))
       end if
    end do

    ! fill ghost cells for gradient via message pass
    if(present(GradVarInOut_DGB))then
       call message_pass_cell(&
            nVar            = nDim, &
            nG              = 1,    &
            State_VGB       = GradVarInOut_DGB, &
            nProlongOrderIn = 1)
    else
       call message_pass_cell(&
            nVar            = nDim, &
            nG              = 1,    &
            State_VGB       = GradVar_DGB, &
            nProlongOrderIn = 1)
    end if

    call test_stop(NameSub, DoTest)
  end subroutine calc_gradient_ghost
  !============================================================================

!  subroutine calc_cell_curl_ghost(iBlock, Var_DG, nG, curl_DG, UseBodyCellIn)
!
!    ! Calculate curl of Var_DG and return it in Curl_DG
!    ! Physical cells and 1 layer of ghost cells are calculated (nG>=2)
!    ! Body (false) cells are ignored unless UseBodyCellIn is set to true.
!    ! Calculate curl of Var_DG on a Cartesian grid. 
!    ! Need to include general coordinates calculation in the future.
!
!    use ModGeometry, ONLY: body_blk, true_cell
!    !hyzhou test
!    use ModAdvance,  ONLY: State_VGB, Bx_, By_, Bz_
!
!
!    integer, intent(in):: iBlock
!    real, intent(in) :: Var_DG(3,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
!    integer, intent(in):: nG ! number of ghost cells in curl_DG
!    real, intent(inout):: &  
!         curl_DG(3,1-nG:nI+nG,1-nG*jDim_:nJ+nG*jDim_,1-nG*kDim_:nK+nG*kDim_)
!    logical, intent(in), optional:: UseBodyCellIn
!
!    logical:: UseBodyCell
!
!    real:: InvDxHalf, InvDyHalf, InvDzHalf
!    integer :: i, j, k
!
!    logical:: DoTest
!    character(len=*), parameter:: NameSub = 'calc_cell_curl'
!    !--------------------------------------------------------------------------
!    call test_start(NameSub, DoTest, iBlock)
!    UseBodyCell = .false.
!    if(present(UseBodyCellIn)) UseBodyCell = UseBodyCellIn
!
!    if(UseBodyCell .or. .not. body_blk(iBlock)) then
!       if(IsCartesian)then
!
!          ! For test comparison with get_current
!          !Var_DG = State_VGB(Bx_:Bz_,:,:,:,iBlock)
!
!          ! Simple central differencing
!          
!          InvDxHalf = 0.5/CellSize_DB(1,iBlock)
!          InvDyHalf = 0.5/CellSize_DB(2,iBlock)
!          InvDzHalf = 0.5/CellSize_DB(3,iBlock)
!
!          ! Calculate curl for 1 layer of ghost cells
!          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
!             curl_DG(x_,i,j,k) = &
!                  InvDyHalf*( Var_DG(z_,i,j+1,k) - Var_DG(z_,i,j-1,k) ) - &
!                  InvDzHalf*( Var_DG(y_,i,j,k+1) - Var_DG(y_,i,j,k-1) )
!             
!             curl_DG(y_,i,j,k) = &
!                  InvDzHalf*( Var_DG(x_,i,j,k+1) - Var_DG(x_,i,j,k-1) ) - &
!                  InvDxHalf*( Var_DG(z_,i+1,j,k) - Var_DG(z_,i-1,j,k) )
!             
!             curl_DG(z_,i,j,k) = &
!                  InvDxHalf*( Var_DG(y_,i+1,j,k) - Var_DG(y_,i-1,j,k) ) - &
!                  InvDyHalf*( Var_DG(x_,i,j+1,k) - Var_DG(x_,i,j-1,k) )
!          end do; end do; end do  
!       else
!          ! Need to be implemented for general coordinates
!       end if
!    end if
!
!    call test_stop(NameSub, DoTest, iBlock)    
!  end subroutine calc_cell_curl_ghost
  !============================================================================

  ! idea copy from get_current
  subroutine calc_cell_curl_ghost(i,j,k,iBlock,Vector_DG,Curl_D)

    use BATL_lib, ONLY: IsCartesianGrid, IsRzGeometry, Xyz_DGB, CellSize_DB
    use ModParallel, ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot
    use ModSize,     ONLY: nI, nJ, nK, x_, y_, z_
    use ModGeometry, ONLY: True_Cell, true_BLK
    
    integer, intent(in) :: i, j, k, iBlock
    real,    intent(in) :: Vector_DG(3,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real,    intent(out):: Curl_D(3)
    
    integer:: iL, iR, jL, jR, kL, kR
    real   :: Ax, Bx, Cx, Ay, By, Cy, Az, Bz, Cz
    real   :: InvDx2, InvDy2, InvDz2
    !--------------------------------------------------------------------------
!    if(.not.True_Cell(i,j,k,iBlock))then
!       Curl_D = 0.0
!       RETURN
!    endif

    InvDx2 = 0.5/CellSize_DB(x_,iBlock)
    InvDy2 = 0.5/CellSize_DB(y_,iBlock)
    InvDz2 = 0.5/CellSize_DB(z_,iBlock)

    ! Central difference
    iR = i+1; iL = i-1;
    Ax = -InvDx2; Bx = 0.0; Cx = +InvDx2

    jR = j+1; jL = j-1;
    Ay = -InvDy2; By = 0.0; Cy = +InvDy2

    kR = k+1; kL = k-1
    Az = -InvDz2; Bz = 0.0; Cz = +InvDz2


    if(IsCartesianGrid)then
       call calc_cartesian_curl
    else
       call calc_gencoord_curl
    end if

  contains
    !==========================================================================
    subroutine calc_cartesian_curl

      !------------------------------------------------------------------------
      Curl_D(x_) = &
           + Ay*Vector_DG(z_,i,jL,k ) &
           + By*Vector_DG(z_,i,j ,k ) &
           + Cy*Vector_DG(z_,i,jR,k ) &
           - Az*Vector_DG(y_,i,j ,kL) &
           - Bz*Vector_DG(y_,i,j ,k ) &
           - Cz*Vector_DG(y_,i,j ,kR)

      Curl_D(y_) = &
           + Az*Vector_DG(x_,i ,j,kL) &
           + Bz*Vector_DG(x_,i ,j,k ) &
           + Cz*Vector_DG(x_,i ,j,kR) &
           - Ax*Vector_DG(z_,iL,j,k ) &
           - Bx*Vector_DG(z_,i ,j,k ) &
           - Cx*Vector_DG(z_,iR,j,k )

      Curl_D(z_) = &
           + Ax*Vector_DG(y_,iL,j ,k) &
           + Bx*Vector_DG(y_,i ,j ,k) &
           + Cx*Vector_DG(y_,iR,j ,k) &
           - Ay*Vector_DG(x_,i ,jL,k) &
           - By*Vector_DG(x_,i ,j ,k) &
           - Cy*Vector_DG(x_,i ,jR,k)

      ! Correct curl for rz-geometry: Curl(z) = Curl(z) + Bphi/radius
      if(IsRzGeometry) Curl_D(x_) = Curl_D(x_) &
           + Vector_DG(z_,i,j,k)/Xyz_DGB(y_,i,j,k,iBlock)

    end subroutine calc_cartesian_curl
    !==========================================================================

    subroutine calc_gencoord_curl

      use ModCoordTransform, ONLY: inverse_matrix

      real :: DxyzDcoord_DD(3,3), DcoordDxyz_DD(3,3), DbDcoord_DD(3,3)

      ! Get the dCartesian/dGencoord matrix with central difference
      !------------------------------------------------------------------------
      DxyzDcoord_DD(:,1) = InvDx2 &
           *(Xyz_DGB(:,i+1,j,k,iBlock) - Xyz_DGB(:,i-1,j,k,iBlock))

      DxyzDcoord_DD(:,2) = InvDy2 &
           *(Xyz_DGB(:,i,j+1,k,iBlock) - Xyz_DGB(:,i,j-1,k,iBlock))

      if(nK > 1)then
         DxyzDcoord_DD(:,3) = InvDz2 &
              *(Xyz_DGB(:,i,j,k+1,iBlock) - Xyz_DGB(:,i,j,k-1,iBlock))
      else
         DxyzDcoord_DD(:,3) = (/ 0., 0., 1./)
      end if

      DcoordDxyz_DD = inverse_matrix(DxyzDcoord_DD, DoIgnoreSingular=.true.)

      ! Calculate the partial derivatives dB/dGencoord
      DbDcoord_DD(:,1) = &
           + Ax*Vector_DG(x_:z_,iL,j,k) &
           + Bx*Vector_DG(x_:z_,i ,j,k) &
           + Cx*Vector_DG(x_:z_,iR,j,k)

      DbDcoord_DD(:,2) = &
           + Ay*Vector_DG(x_:z_,i,jL,k) &
           + By*Vector_DG(x_:z_,i,j ,k) &
           + Cy*Vector_DG(x_:z_,i,jR,k)

      if(nK > 1)then
         DbDcoord_DD(:,3) = &
              + Az*Vector_DG(x_:z_,i,j,kL) &
              + Bz*Vector_DG(x_:z_,i,j,k ) &
              + Cz*Vector_DG(x_:z_,i,j,kR)
      else
         DbDcoord_DD(:,3) = 0.0
      end if

      ! curl_x = Dbz/Dy - Dby/Dz = Dbz/Dcoord.Dcoord/Dy - DBy/Dcoord.Dccord/dz
      Curl_D(x_) = &
           + sum(DbDcoord_DD(z_,:)*DcoordDxyz_DD(:,y_)) &
           - sum(DbDcoord_DD(y_,:)*DcoordDxyz_DD(:,z_))

      ! curl_y = Dbx/Dz - Dbz/Dx
      Curl_D(y_) = &
           + sum(DbDcoord_DD(x_,:)*DcoordDxyz_DD(:,z_)) &
           - sum(DbDcoord_DD(z_,:)*DcoordDxyz_DD(:,x_))

      ! curl_z = Dby/Dx - Dbx/Dy
      Curl_D(z_) = &
           + sum(DbDcoord_DD(y_,:)*DcoordDxyz_DD(:,x_)) &
           - sum(DbDcoord_DD(x_,:)*DcoordDxyz_DD(:,y_))

    end subroutine calc_gencoord_curl
    !==========================================================================
  end subroutine calc_cell_curl_ghost
  !============================================================================
  
end module ModCellGradient
!==============================================================================
