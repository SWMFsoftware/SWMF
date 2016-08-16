module ModCellGradient

  ! Calculate cell centered gradient and divergence

  use BATL_lib, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nI, nJ, nK, &
       CellSize_DB,  CellFace_DB, FaceNormal_DDFB, CellVolume_GB, &
       IsCartesian, IsCartesianGrid, &
       nDim, jDim_, kDim_, x_, y_, z_, Dim1_, Dim2_, Dim3_
  use ModGeometry, ONLY: body_blk, true_cell

  implicit none

  SAVE
  private ! except

  public:: calc_divergence
  public:: calc_gradient

  integer, allocatable:: iTrue_G(:,:,:)


  interface calc_gradient
     module procedure calc_gradient1, calc_gradient3
  end interface

contains

  !==========================================================================
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

    character(len=*), parameter:: NameSub = 'calc_divergence'
    !------------------------------------------------------------------------
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

  end subroutine calc_divergence

  !==========================================================================
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

    character(len=*), parameter:: NameSub = 'calc_gradient1'
    !-------------------------------------------------------------------------
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

  end subroutine calc_gradient1

  !==========================================================================
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
    !------------------------------------------------------------------------
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
  end subroutine calc_gradient3

end module ModCellGradient
