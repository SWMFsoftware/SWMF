!^CFG COPYRIGHT UM
!==============================================================================
subroutine fix_block_geometry(iBlock)

  use ModMain, ONLY: body1_, body2_, ExtraBc_,&
       UseExtraBoundary, ProcTest, BlkTest   
  use ModMain, ONLY: UseBody2                    !^CFG IF SECONDBODY
  use ModGeometry ! , ONLY: true_cell, dx_blk, dy_blk, dz_blk
  use ModPhysics, ONLY : xBody2,yBody2,zBody2 !^CFG IF SECONDBODY
  use ModParallel, ONLY : periodic3D
  use ModBoundaryCells
  use BATL_lib, ONLY: Xyz_DGB
  implicit none

  integer, intent(in) :: iBlock

  integer :: i,j,k, iBoundary

  logical:: DoTest, DoTestMe
  character(len=*), parameter :: NameSub='fix_block_geometry'
  !---------------------------------------------------------------------------

  if(iBlock==BlkTest .and. iProc==ProcTest)then
     call set_oktest(NameSub, DoTest, DoTestMe)
  else
     DoTest = .false.; DoTestMe = .false.
  end if

  if (UseBody2) then                        !^CFG IF SECONDBODY BEGIN
     ! calculate the radius as measured from the second body
     ! Note that the second body can move
     R2_BLK(:,:,:,iBlock) = sqrt( &
          (Xyz_DGB(x_,:,:,:,iBlock)-xBody2)**2 + &
          (Xyz_DGB(y_,:,:,:,iBlock)-yBody2)**2 + &
          (Xyz_DGB(z_,:,:,:,iBlock)-zBody2)**2)
     Rmin2_BLK(iBlock) = minval(R2_BLK(:,:,:,iBlock))
  else
     Rmin2_BLK(iBlock) = 0.0
  end if                                    !^CFG END SECONDBODY

  far_field_BCs_BLK(iBlock) = &
       (((xyzStart_BLK(1,iBlock)-dx_BLK(iBlock))<XyzMin_D(1).or.&
       (xyzStart_BLK(1,iBlock)+nI*dx_BLK(iBlock))>XyzMax_D(1)) &
       .and. .not.periodic3D(1)) .or. &
       (((xyzStart_BLK(2,iBlock)-dy_BLK(iBlock))<XyzMin_D(2).or.&
       (xyzStart_BLK(2,iBlock)+nJ*dy_BLK(iBlock))>XyzMax_D(2)) &
       .and. .not.periodic3D(2)).or. &
       (((xyzStart_BLK(3,iBlock)-dz_BLK(iBlock))<XyzMin_D(3).or.&
       (xyzStart_BLK(3,iBlock)+nK*dz_BLK(iBlock))>XyzMax_D(3)) &
       .and. .not.periodic3D(3))  


  if(DoTestMe)then
     write(*,*)NameSub,': far_field_bcs_blk=',far_field_bcs_BLK(iBlock)
     write(*,*)NameSub,': xyzStart_BLK=',xyzStart_BLK(:,BlkTest)
     write(*,*)NameSub,': dx,dy,dz= ',dx_BLK(BlkTest),dy_BLK(iBlock),dz_BLK(iBlock)
     write(*,*)NameSub,': xyzmin=',XyzMin_D(:)
     write(*,*)NameSub,': xyzmax=',XyzMax_D(:)
  end if

  !\
  ! TRUE_CELL: if not inside a body or outside the outer face boundary
  !/
  true_cell(:,:,:,iBlock)=.true.
  IsBoundaryBlock_IB(:,iBlock)=.false.
  do iBoundary = MinBoundary, MaxBoundary
     IsBoundaryBlock_IB(iBoundary,iBlock) = .true.
  end do
  IsBoundaryBlock_IB(ExtraBc_,iBlock) = UseExtraBoundary

  ! set true_cell array (seting IsBoundaryCell_GI)
  call set_boundary_cells(iBlock)

  do iBoundary = MinBoundary, min(MaxBoundary,Body1_)
     IsBoundaryBlock_IB(iBoundary,iBlock)=any(IsBoundaryCell_GI(:,:,:,iBoundary))
     true_cell(:,:,:,iBlock) = &
          true_cell(:,:,:,iBlock) .and. .not.IsBoundaryCell_GI(:,:,:,iBoundary)
  end do

  IsBoundaryCell_GI(:,:,:,ExtraBc_) = &
       UseExtraBoundary .and. IsBoundaryCell_GI(:,:,:,ExtraBc_)

  ! Copying  the IsBoundaryCell_GI into the format for iBoundary_GB
  iBoundary_GB(:,:,:,iBlock) = domain_
  do iBoundary = body2_, 6
     where(IsBoundaryCell_GI(:,:,:,iBoundary))
        iBoundary_GB(:,:,:,iBlock) = iBoundary
     end where
  end do

  ! Alow other places to set true_cell
  true_cell(1:nI,1:nJ,1:nK,iBlock) = true_cell(1:nI,1:nJ,1:nK,iBlock) &
       .and. iBoundary_GB(1:nI,1:nJ,1:nK,iBlock) == domain_

  ! body_BLK: if any cell INCLUDING ghost cells is inside body(ies)
  body_BLK(iBlock) = .not. all(true_cell(:,:,:,iBlock))

  ! TRUE_BLK: if all cells EXCLUDING ghost cells are outside body(ies)
  true_BLK(iBlock) = all(true_cell(1:nI,1:nJ,1:nK,iBlock))

end subroutine fix_block_geometry

!=============================================================================

subroutine set_boundary_cells(iBlock)

  use ModProcMH
  use ModMain
  use ModPhysics,  ONLY: Rbody
  use ModGeometry, ONLY: R_BLK, IsBoundaryBlock_IB, IsBoundaryCell_GI
  use ModPhysics,  ONLY: Rbody2                            !^CFG IF SECONDBODY
  use ModGeometry, ONLY: R2_BLK                            !^CFG IF SECONDBODY
  use ModGeometry,ONLY:x1,x2,y1,y2,z1,z2
  use BATL_lib, ONLY: Xyz_DGB
  use ModUser, ONLY: user_set_boundary_cells

  implicit none
  integer, intent(in)::iBlock
  !----------------------------------------------------------------------------

  IsBoundaryCell_GI=.false.  
  !^CFG IF SECONDBODY BEGIN
  if(UseBody2 .and. IsBoundaryBlock_IB(Body2_,iBlock)) &
       IsBoundaryCell_GI(:,:,:,Body2_) = R2_BLK(:,:,:,iBlock) < rBody2  
  !^CFG END SECONDBODY

  if(IsBoundaryBlock_IB(Body1_,iBlock)) &
       IsBoundaryCell_GI(:,:,:,Body1_) = &
       body1    .and. R_BLK(:,:,:,iBlock) < Rbody

  if(IsBoundaryBlock_IB(ExtraBc_,iBlock))&             
       call user_set_boundary_cells(iBlock)

  if(IsBoundaryBlock_IB(1,iBlock)) &
       IsBoundaryCell_GI(:,:,:,1)= Xyz_DGB(x_,:,:,:,iBlock) < x1  
  if(IsBoundaryBlock_IB(2,iBlock)) &
       IsBoundaryCell_GI(:,:,:,2)= Xyz_DGB(x_,:,:,:,iBlock) > x2 
  if(IsBoundaryBlock_IB(3,iBlock)) &
       IsBoundaryCell_GI(:,:,:,3)= Xyz_DGB(y_,:,:,:,iBlock) < y1
  if(IsBoundaryBlock_IB(4,iBlock)) &
       IsBoundaryCell_GI(:,:,:,4)= Xyz_DGB(y_,:,:,:,iBlock) > y2 
  if(IsBoundaryBlock_IB(5,iBlock)) &
       IsBoundaryCell_GI(:,:,:,5)= Xyz_DGB(z_,:,:,:,iBlock) < z1
  if(IsBoundaryBlock_IB(6,iBlock)) &
       IsBoundaryCell_GI(:,:,:,6)= Xyz_DGB(z_,:,:,:,iBlock) > z2  


end subroutine set_boundary_cells
