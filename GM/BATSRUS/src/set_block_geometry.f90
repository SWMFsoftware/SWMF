!^CFG COPYRIGHT UM
!==============================================================================
subroutine fix_block_geometry(iBLK)

  use ModMain, ONLY:body1_,body2_,ExtraBc_,&
       UseExtraBoundary,ProcTest,BlkTest   
  use ModMain, ONLY: UseBody2                    !^CFG IF SECONDBODY
  use ModNodes
  use ModGeometry
  use ModNumConst
  use ModPhysics, ONLY : xBody2,yBody2,zBody2 !^CFG IF SECONDBODY
  use ModParallel, ONLY : periodic3D
  use ModBoundaryCells
  implicit none

  integer, intent(in) :: iBLK

  integer :: i,j,k, iBoundary

  logical:: DoTest, DoTestMe
  character(len=*), parameter :: NameSub='fix_block_geometry'
  !---------------------------------------------------------------------------

  if(iBlk==BlkTest .and. iProc==ProcTest)then
     call set_oktest(NameSub, DoTest, DoTestMe)
  else
     DoTest = .false.; DoTestMe = .false.
  end if

  if (UseBody2) then                        !^CFG IF SECONDBODY BEGIN
     ! calculate the radius as measured from the second body
     ! Note that the second body can move
     R2_BLK(:,:,:,iBLK) = sqrt( &
          (x_BLK(:,:,:,iBLK)-xBody2)**2 + &
          (y_BLK(:,:,:,iBLK)-yBody2)**2 + &
          (z_BLK(:,:,:,iBLK)-zBody2)**2)
     Rmin2_BLK(iBLK) = minval(R2_BLK(:,:,:,iBLK))
  else
     Rmin2_BLK(iBLK) = 0.0
  end if                                    !^CFG END SECONDBODY

  far_field_BCs_BLK(iBLK) = &
       (((xyzStart_BLK(1,iBLK)-dx_BLK(iBLK))<XyzMin_D(1).or.&
       (xyzStart_BLK(1,iBLK)+nI*dx_BLK(iBLK))>XyzMax_D(1)) &
       .and. .not.periodic3D(1)) .or. &
       (((xyzStart_BLK(2,iBLK)-dy_BLK(iBLK))<XyzMin_D(2).or.&
       (xyzStart_BLK(2,iBLK)+nJ*dy_BLK(iBLK))>XyzMax_D(2)) &
       .and. .not.periodic3D(2)).or. &
       (((xyzStart_BLK(3,iBLK)-dz_BLK(iBLK))<XyzMin_D(3).or.&
       (xyzStart_BLK(3,iBLK)+nK*dz_BLK(iBLK))>XyzMax_D(3)) &
       .and. .not.periodic3D(3))  


  if(DoTestMe)then
     write(*,*)NameSub,': far_field_bcs_blk=',far_field_bcs_BLK(iBlk)
     write(*,*)NameSub,': xyzStart_BLK=',xyzStart_BLK(:,BlkTest)
     write(*,*)NameSub,': dx,dy,dz= ',dx_BLK(BlkTest),dy_BLK(iBLK),dz_BLK(iBLK)
     write(*,*)NameSub,': xyzmin=',XyzMin_D(:)
     write(*,*)NameSub,': xyzmax=',XyzMax_D(:)
  end if

  !\
  ! TRUE_CELL: if not inside a body or outside the outer face boundary
  !/
  true_cell(:,:,:,iBLK)=.true.
  IsBoundaryBlock_IB(:,iBLK)=.false.
  do iBoundary = MinBoundary, MaxBoundary
     IsBoundaryBlock_IB(iBoundary,iBLK) = .true.
  end do
  IsBoundaryBlock_IB(ExtraBc_,iBLK) = UseExtraBoundary

  ! set true_cell array (seting IsBoundaryCell_GI)
  call set_boundary_cells(iBLK)

  do iBoundary = MinBoundary, min(MaxBoundary,Body1_)
     IsBoundaryBlock_IB(iBoundary,iBLK)=any(IsBoundaryCell_GI(:,:,:,iBoundary))
     true_cell(:,:,:,iBLK) = &
          true_cell(:,:,:,iBLK) .and. .not.IsBoundaryCell_GI(:,:,:,iBoundary)
  end do

  IsBoundaryCell_GI(:,:,:,ExtraBc_) = &
       UseExtraBoundary .and. IsBoundaryCell_GI(:,:,:,ExtraBc_)

  ! Copying  the IsBoundaryCell_GI into the format for iBoundary_GB
  iBoundary_GB(:,:,:,iBlk) = domain_
  do iBoundary = body2_, 6
     where(IsBoundaryCell_GI(:,:,:,iBoundary))
        iBoundary_GB(:,:,:,iBlk) = iBoundary
     end where
  end do

  ! Alow other places to set true_cell
  true_cell(1:nI,1:nJ,1:nK,iBLK) = true_cell(1:nI,1:nJ,1:nK,iBLK) &
       .and. iBoundary_GB(1:nI,1:nJ,1:nK,iBlk) == domain_

  ! body_BLK: if any cell INCLUDING ghost cells is inside body(ies)
  body_BLK(iBLK) = .not. all(true_cell(:,:,:,iBLK))

  ! TRUE_BLK: if all cells EXCLUDING ghost cells are outside body(ies)
  true_BLK(iBLK) = all(true_cell(1:nI,1:nJ,1:nK,iBLK))

end subroutine fix_block_geometry

!=============================================================================

subroutine set_boundary_cells(iBLK)

  use ModProcMH
  use ModMain
  use ModPhysics,  ONLY: Rbody
  use ModGeometry, ONLY: R_BLK, IsBoundaryBlock_IB, IsBoundaryCell_GI
  use ModPhysics,  ONLY: Rbody2                            !^CFG IF SECONDBODY
  use ModGeometry, ONLY: R2_BLK                            !^CFG IF SECONDBODY
  use ModGeometry,ONLY:x1,x2,y1,y2,z1,z2,x_BLK,y_BLK,z_BLK 
  use ModUser, ONLY: user_set_boundary_cells

  implicit none
  integer,intent(in)::iBLK
  !----------------------------------------------------------------------------

  IsBoundaryCell_GI=.false.  
  !^CFG IF SECONDBODY BEGIN
  if(UseBody2 .and. IsBoundaryBlock_IB(Body2_,iBLK)) &
       IsBoundaryCell_GI(:,:,:,Body2_) = R2_BLK(:,:,:,iBLK) < rBody2  
  !^CFG END SECONDBODY

  if(IsBoundaryBlock_IB(Body1_,iBLK)) &
       IsBoundaryCell_GI(:,:,:,Body1_) = &
       body1    .and. R_BLK(:,:,:,iBLK) < Rbody


  if(IsBoundaryBlock_IB(ExtraBc_,iBLK))&             
       call user_set_boundary_cells(iBLK)



  if(IsBoundaryBlock_IB(1,iBLK)) &
       IsBoundaryCell_GI(:,:,:,1)=x_BLK(:,:,:,iBLK)<x1  
  if(IsBoundaryBlock_IB(2,iBLK)) &
       IsBoundaryCell_GI(:,:,:,2)=x_BLK(:,:,:,iBLK)>x2 
  if(IsBoundaryBlock_IB(3,iBLK)) &
       IsBoundaryCell_GI(:,:,:,3)= y_BLK(:,:,:,iBLK)<y1
  if(IsBoundaryBlock_IB(4,iBLK)) &
       IsBoundaryCell_GI(:,:,:,4)= y_BLK(:,:,:,iBLK)>y2 
  if(IsBoundaryBlock_IB(5,iBLK)) &
       IsBoundaryCell_GI(:,:,:,5)= z_BLK(:,:,:,iBLK)<z1
  if(IsBoundaryBlock_IB(6,iBLK)) &
       IsBoundaryCell_GI(:,:,:,6)= z_BLK(:,:,:,iBLK)>z2  


end subroutine set_boundary_cells
