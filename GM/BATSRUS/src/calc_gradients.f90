!^CFG COPYRIGHT UM
!==============================================================================
subroutine grad1D(iObsolete, iBlock, Var_G, &
      DifferenceX_G, DifferenceY_G, DifferenceZ_G, TypeObsolete, jObsolete)

  ! This is an interface to central_difference or covariant_gradient 
  ! programs.
  ! For using in amr_criteria only

  use ModMain, ONLY : nI,nJ,nK,gcn
  use ModGeometry,ONLY: UseCovariant              
  implicit none

  integer,           intent(in) :: iBlock
  integer, intent(in) :: iObsolete,jObsolete               !Obsolete
  real, intent(in), dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn) :: Var_G
  real, intent(out), dimension( 0:nI+1, 0:nJ+1, 0:nK+1) :: &
        DifferenceX_G, DifferenceY_G, DifferenceZ_G
  character (LEN=*), intent(in) :: TypeObsolete            !Obsolete

  !--------------------------------------------------------------------------
  if(UseCovariant)then                               
     call covariant_gradient(iBlock, Var_G,&  
          DifferenceX_G, DifferenceY_G, DifferenceZ_G)  
  else                                                
     call central_differences(iBlock, Var_G,&         
          DifferenceX_G, DifferenceY_G, DifferenceZ_G)
     continue
  end if                                             
end subroutine grad1D

!==============================================================================
subroutine central_differences(iBlock, Var_G,&     
     DifferenceX_G, DifferenceY_G, DifferenceZ_G)
  use ModSize
  use ModGeometry,ONLY:body_blk, true_cell, &
       fAX_BLK, fAY_BLK, fAZ_BLK, &
       vInv_CB
  use ModNumConst
  implicit none

  integer,intent(in) :: iBlock

  real, dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn),&
       intent(in) :: Var_G
  real, dimension(0:nI+1, 0:nJ+1, 0:nK+1),&
       intent(out) ::  DifferenceX_G, DifferenceY_G, DifferenceZ_G

  real, dimension(0:nI+1, 0:nJ+1, 0:nK+1) :: OneTrue_G

  integer :: i, j, k
  !---------------------------------------------------

  !To fill in the ghostcells
  DifferenceX_G = cZero
  DifferenceY_G = cZero
  DifferenceZ_G = cZero

  if(.not.body_blk(iBlock)) then
     do k=1,nK; do j=1,nJ; do i=1,nI
        DifferenceX_G(i,j,k) = cHalf*fAX_BLK(iBlock)*&
             (Var_G(i+1, j, k)-Var_G(i-1,j,k))*vInv_CB(i,j,k,iBlock)
        DifferenceY_G(i,j,k) = cHalf*fAY_BLK(iBlock)*&
             (Var_G(i, j+1, k)-Var_G(i,j-1,k))*vInv_CB(i,j,k,iBlock)
        DifferenceZ_G(i,j,k) = cHalf*fAZ_BLK(iBlock)*&
             (Var_G(i, j, k+1)-Var_G(i, j, k-1))*vInv_CB(i,j,k,iBlock)
     end do; end do; end do
  else
     where(true_cell(0:nI+1, 0:nJ+1, 0:nK+1,iBlock)) 
        OneTrue_G=cOne
     elsewhere
        OneTrue_G=cZero
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
        DifferenceX_G(i,j,k) = cHalf*fAX_BLK(iBlock)*&
             OneTrue_G(i,j,k)*(&
             (Var_G(i+1,j,k)-Var_G(i,j,k))*&
             OneTrue_G(i+1,j,k)*&
             (2.0 - OneTrue_G(i-1,j,k)) + &
             (Var_G(i,j,k)-Var_G(i-1,j,k))*&
             OneTrue_G(i-1,j,k)*&
             (2.0 - OneTrue_G(i+1,j,k)) )*vInv_CB(i,j,k,iBlock)

        DifferenceY_G(i,j,k) = cHalf*fAY_BLK(iBlock)*&
             OneTrue_G(i,j,k)*(&
             (Var_G(i,j+1,k)-Var_G(i,j,k))*&
             OneTrue_G(i,j+1,k)*&
             (2.0 - OneTrue_G(i,j-1,k))+&
             (Var_G(i,j,k)-Var_G(i,j-1,k))*&
             OneTrue_G(i,j-1,k)*&
             (2.0 - OneTrue_G(i,j+1,k)) )*vInv_CB(i,j,k,iBlock)

        DifferenceZ_G(i,j,k) = cHalf*fAZ_BLK(iBlock)*&
             OneTrue_G(i,j,k)*(&
             (Var_G(i,j,k+1)-Var_G(i,j,k))*&
             OneTrue_G(i,j,k+1)*&
             (2.0 - OneTrue_G(i,j,k-1))+&
             (Var_G(i,j,k)-Var_G(i,j,k-1))*&
             OneTrue_G(i,j,k-1)*&
             (2.0 - OneTrue_G(i,j,k+1)) )*vInv_CB(i,j,k,iBlock)
     end do; end do; end do
  end if
end subroutine central_differences

