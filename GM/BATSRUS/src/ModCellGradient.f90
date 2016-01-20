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

  real, allocatable:: iTrue_G(:,:,:)

contains

  !==========================================================================
  subroutine calc_divergence(iBlock, Var_DG, nG, Div_G)

    ! This is an interface to cartesian or gen. coord divergence

    integer, intent(in):: iBlock
    real, intent(in) :: Var_DG(nDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    integer, intent(in):: nG ! number of ghost cells in Div_G
    real, intent(out):: &
         Div_G(1-nG:nI+nG,1-nG*jDim_:nJ+nG*jDim_,1-nG*kDim_:nK+nG*kDim_)

    character(len=*), parameter:: NameSub = 'calc_divergence'
    !------------------------------------------------------------------------
    if(IsCartesian)then
       call calc_divergence_cart(iBlock, Var_DG, nG, Div_G)
    else
       call stop_mpi(NameSub//': non-cartesian to be implemented')
    end if

  end subroutine calc_divergence

  !==========================================================================
  subroutine calc_gradient(iBlock, Var_G, GradX_C, GradY_C, GradZ_C)

    ! This is an interface to cartesian or gencoord_gradient.

    integer, intent(in):: iBlock
    real,    intent(in):: Var_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(out):: GradX_C(nI,nJ,nK), GradY_C(nI,nJ,nK), GradZ_C(nI,nJ,nK)
    !------------------------------------------------------------------------
    if(IsCartesianGrid)then
       call cartesian_gradient(iBlock, Var_G, GradX_C, GradY_C, GradZ_C)
    else
       call gencoord_gradient(iBlock, Var_G, GradX_C, GradY_C, GradZ_C)  
    end if

  end subroutine calc_gradient

  !==========================================================================
  subroutine cartesian_gradient(iBlock, Var_G, GradX_C,GradY_C,GradZ_C)

    use ModGeometry, ONLY: body_blk, true_cell

    integer,intent(in) :: iBlock

    real, intent(in) :: Var_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK)
    real, intent(out):: GradX_C(nI,nJ,nK), GradY_C(nI,nJ,nK), GradZ_C(nI,nJ,nK)

    real:: OneTrue_G(0:nI+1,0:nJ+1,0:nK+1)
    real:: VInvHalf
    integer :: i, j, k
    !------------------------------------------------------------------------

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

  end subroutine cartesian_gradient
  !==========================================================================
  subroutine gencoord_gradient(iBlock, Var_G, GradX_C, GradY_C, GradZ_C)

    integer,intent(in) :: iBlock

    real, intent(in) :: &
         Var_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(out):: &
         GradX_C(nI,nJ,nK), GradY_C(nI,nJ,nK), GradZ_C(nI,nJ,nK)

    real:: OneTrue_G(0:nI+1,0:nJ+1,0:nK+1)

    integer :: i, j, k

    real :: FaceArea_DS(3,6), Difference_S(6)

    real::VInvHalf
    !------------------------------------------------------------------------
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

  end subroutine gencoord_gradient
  !==========================================================================

  subroutine calc_divergence_cart(iBlock, Var_DG, nG, Div_G)

    ! Calculate divergence Div_G of Var_DG on a Cartesian grid

    use ModGeometry, ONLY: body_blk, true_cell

    integer,intent(in) :: iBlock
    real, intent(in) :: Var_DG(nDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    integer, intent(in):: nG ! number of ghost cells in Div_G
    real, intent(out):: &
         Div_G(1-nG:nI+nG,1-nG*jDim_:nJ+nG*jDim_,1-nG*kDim_:nK+nG*kDim_)

    real:: InvDxHalf, InvDyHalf, InvDzHalf, InvDx, InvDy, InvDz
    integer :: i, j, k, iL, iR, jL, jR, kL, kR
    !------------------------------------------------------------------------
    if(.not.body_blk(iBlock)) then

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
       InvDx = 1/CellSize_DB(1,iBlock)
       InvDy = 1/CellSize_DB(2,iBlock)
       InvDz = 1/CellSize_DB(3,iBlock)

       ! Set iTrue_G to 1 in used cells and 0 in non-used cells
       if(.not.allocated(iTrue_G)) allocate(iTrue_G(0:nI+1,0:nJ+1,0:nK+1))
       where(true_cell(0:nI+1,0:nJ+1,0:nK+1,iBlock)) 
          iTrue_G = 1
       elsewhere
          iTrue_G = 0
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
          Div_G(i,j,k) = 0.0
          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          ! Shift to left and right if the neighbor is a true cell
          iL = i - iTrue_G(i-1,j,k); iR = i + iTrue_G(i+1,j,k)

          ! If at least one neighbor is true, calculate derivative
          if(iL == iR) Div_G(i,j,k) = &
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
    end if

  end subroutine calc_divergence_cart

end module ModCellGradient
