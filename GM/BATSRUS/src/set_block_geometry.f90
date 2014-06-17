!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!This code is a copyright protected software (c) 2002- University of Michigan
!==============================================================================
subroutine fix_block_geometry(iBlock)

  use ModMain, ONLY: body1_, body2_, ExtraBc_,&
       UseExtraBoundary, ProcTest, BlkTest   
  use ModMain, ONLY: UseBody2
  use ModGeometry, ONLY: &
       MinFaceBoundary, MaxFaceBoundary, IsBoundaryBlock_IB, &
       IsBoundaryCell_GI, R2_BLK, Rmin2_BLK, Body_BLK, &
       far_field_BCs_BLK, true_blk, true_cell
  use ModPhysics, ONLY : xBody2,yBody2,zBody2
  use ModBoundaryCells
  use BATL_lib, ONLY: iProc, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nI, nJ, nK, &
       Xyz_DGB, CellSize_DB, &
       CoordMin_DB, CoordMax_DB, CoordMin_D, CoordMax_D, IsPeriodic_D

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

  if (UseBody2) then
     ! calculate the radius as measured from the second body
     ! Note that the second body can move
     do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
        R2_BLK(i,j,k,iBlock) = sqrt( &
          sum( (Xyz_DGB(:,i,j,k,iBlock)- (/xBody2, yBody2, zBody2/))**2 ))
     end do; end do; end do
     Rmin2_BLK(iBlock) = minval(R2_BLK(:,:,:,iBlock))
  else
     Rmin2_BLK(iBlock) = 0.0
  end if

  far_field_BCs_BLK(iBlock) = any( .not. IsPeriodic_D .and.   &
       (    CoordMin_DB(:,iBlock) <= CoordMin_D + 1e-12       &
       .or. CoordMax_DB(:,iBlock) >= CoordMax_D - 1e-12) )

  if(DoTestMe)then
     write(*,*)NameSub,': far_field_bcs_blk=',far_field_bcs_BLK(iBlock)
     write(*,*)NameSub,': CoordMin_DB=', CoordMin_DB(:,iBlock)
     write(*,*)NameSub,': CellSize_D = ',CellSize_DB(:,iBlock)
     write(*,*)NameSub,': CoordMin_D =', CoordMin_D
     write(*,*)NameSub,': CoordMax_D =', CoordMax_D
  end if

  !\
  ! TRUE_CELL: if not inside a body or outside the outer face boundary
  !/
  true_cell(:,:,:,iBlock)=.true.
  IsBoundaryBlock_IB(:,iBlock)=.false.
  do iBoundary = MinFaceBoundary, MaxFaceBoundary
     IsBoundaryBlock_IB(iBoundary,iBlock) = .true.
  end do
  IsBoundaryBlock_IB(ExtraBc_,iBlock) = UseExtraBoundary

  ! set true_cell array (seting IsBoundaryCell_GI)
  call set_boundary_cells(iBlock)

  do iBoundary = MinFaceBoundary, min(MaxFaceBoundary,Body1_)
     IsBoundaryBlock_IB(iBoundary,iBlock) = &
          any(IsBoundaryCell_GI(:,:,:,iBoundary))
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
  use ModPhysics,  ONLY: Rbody2
  use ModGeometry, ONLY: R2_BLK
  use ModGeometry,ONLY:x1,x2,y1,y2,z1,z2
  use BATL_lib, ONLY: Xyz_DGB
  use ModUserInterface ! user_set_boundary_cells

  implicit none
  integer, intent(in)::iBlock
  !----------------------------------------------------------------------------

  IsBoundaryCell_GI=.false.  

  if(UseBody2 .and. IsBoundaryBlock_IB(Body2_,iBlock)) &
       IsBoundaryCell_GI(:,:,:,Body2_) = R2_BLK(:,:,:,iBlock) < rBody2  


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
