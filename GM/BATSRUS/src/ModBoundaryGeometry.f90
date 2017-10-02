!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!============================================================================

module ModBoundaryGeometry

  implicit none
  SAVE

  ! iBoundary_GB contains the index of the boundary that the cell belongs to.
  integer, allocatable :: iBoundary_GB(:,:,:,:)

  ! Cells inside domain have index domain_ that is smaller than smallest 
  ! boundary index (SolidBc_ = -3)
  integer, parameter :: domain_ = -10

  ! Boundary parameter for SolidState
  character(len=15) :: TypeSolidGeometry = 'none'
  real :: rSolid = -1
  real :: SolidLimitDt = 0

contains

  !============================================================================
  subroutine init_mod_boundary_cells

    use ModSize, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, MaxBlock
    !--------------------------------------------------------------------------

    if(.not. allocated(iBoundary_GB)) then
       allocate(iBoundary_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
       iBoundary_GB = domain_
    end if

  end subroutine init_mod_boundary_cells

  !============================================================================

  subroutine read_boundary_geometry_param(NameCommand)

    use ModMain,      ONLY: UseSolidState, TypeCellBc_I, TypeFaceBc_I
    use ModMain,      ONLY: Coord1MinBc_, xMinBc_, solidBc_
    use ModReadParam, ONLY: read_var
    use BATL_size,    ONLY: nDim
    
    integer :: i

    character(len=*), intent(in):: NameCommand
    character(len=*), parameter:: NameSub = 'read_boundary_geometry_param'
    !--------------------------------------------------------------------------

    select case(NameCommand)
    case("#SOLIDSTATE")
       call read_var('UseSolidState', UseSolidState)
       if(UseSolidState)then
          call read_var('TypeBcSolid',TypeFaceBc_I(solidBc_))
          call read_var('TypeSolidGeometry',TypeSolidGeometry)

          select case(TypeSolidGeometry)
          case("sphere")
             call read_var('rSolid',rSolid)
          case("user")
             !?
          case default
             call stop_mpi(NameSub//': Command='//NameCommand//&
                  ', Geometry type='//TypeSolidGeometry//&
                  ' has not been implemented!')
          end select

          call read_var('SolidLimitDt',SolidLimitDt)
       end if
    case("#OUTERBOUNDARY")
       do i = Coord1MinBc_, 2*nDim
          call read_var('TypeCellBc', TypeCellBc_I(i))
       end do
    case("#BOXBOUNDARY")
       do i = xMinBc_, xMinBc_-1+2*nDim
          call read_var('TypeFaceBc', TypeFaceBc_I(i))
       end do

    case default
       call stop_mpi(NameSub//' unknown command='//NameCommand)
    end select

  end subroutine read_boundary_geometry_param

  !============================================================================

  subroutine fix_block_geometry(iBlock, DoSolveSolidIn)

    use ModMain, ONLY: Body1, body1_, body2_, &
         UseBody2, UseExtraBoundary, UseSolidState, &
         ProcTest, BlkTest, iTest, jTest, kTest, TypeFaceBc_I, &
         xMinBc_, xMaxBc_, yMinBc_, yMaxBc_, zMinBc_, zMaxBc_, solidBc_
    use ModGeometry, ONLY: &
         R_BLK, R2_BLK, Rmin2_BLK, Body_BLK, &
         far_field_BCs_BLK, true_blk, true_cell,&
         x1, x2, y1, y2, z1, z2
    use ModPhysics, ONLY : xBody2,yBody2,zBody2, rbody, rBody2
    use BATL_lib, ONLY: &
         iProc, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nI, nJ, nK, nG, &
         Xyz_DGB, CellSize_DB, x_, y_, z_, &
         CoordMin_DB, CoordMax_DB, CoordMin_D, CoordMax_D, IsPeriodic_D
    use ModAdvance, ONLY: time_BLK

    implicit none

    integer, intent(in) :: iBlock
    logical, intent(in), optional :: DoSolveSolidIn

    integer :: i, j, k
    logical :: DoSolveSolid
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='fix_block_geometry'
    !--------------------------------------------------------------------------

    if(iBlock==BlkTest .and. iProc==ProcTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    DoSolveSolid = .true.
    if(present(DoSolveSolidIn)) DoSolveSolid = DoSolveSolidIn

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


    ! TRUE_CELL: if not inside a body or outside the outer face boundary
    ! so true cells are those where iBoundary_GB==domain_
    true_cell(:,:,:,iBlock)=.true.

    ! Reset for every level of refinement                               
    iBoundary_GB(:,:,:,iBlock) = domain_

    ! Extra boundary
    if(UseExtraBoundary) call user_set_boundary_cells(iBlock)

    ! Set boundary type and timestep inside body for solidBC
    if(UseSolidState) then
       if(.not. DoSolveSolid) then
          select case(TypeSolidGeometry)
          case('sphere')
             where( R_BLK(:,:,:,iBlock) < rSolid)&
                  iBoundary_GB(:,:,:,iBlock) = solidBc_
             where( R_BLK(1:nI,1:nJ,1:nK,iBlock) < rSolid)&
                  time_BLK(:,:,:,iBlock) = 0.0        
          end select
       else ! DoSolveSolid=.TRUE.
          select case(TypeSolidGeometry)
          case('sphere')
             where( R_BLK(1:nI,1:nJ,1:nK,iBlock) < rSolid)&
                  time_BLK(:,:,:,iBlock) = SolidLimitDt
          end select
       end if
    end if
    
    ! Set iBoundary_GB for body1 and body2 (always face boundary)            
    if(UseBody2) then
       where( R2_BLK(:,:,:,iBlock) < rbody2) &
          iBoundary_GB(:,:,:,iBlock) = body2_
    end if
    if(body1) then
       where( R_BLK(:,:,:,iBlock) < rbody ) &
         iBoundary_GB(:,:,:,iBlock) = body1_
    end if

    ! No face BC is applied for TypeFaceBc_I(1:6) if not set in PARAM.in   
    if(.not. all(TypeFaceBc_I(xMinBc_:zMaxBc_)=='none') )then
       where( Xyz_DGB(x_,:,:,:,iBlock) < x1 ) &
            iBoundary_GB(:,:,:,iBlock) = xMinBc_
       where( Xyz_DGB(x_,:,:,:,iBlock) > x2 ) &
            iBoundary_GB(:,:,:,iBlock) = xMaxBc_
       where( Xyz_DGB(y_,:,:,:,iBlock) < y1 ) &
            iBoundary_GB(:,:,:,iBlock) = yMinBc_
       where( Xyz_DGB(y_,:,:,:,iBlock) > y2 ) &
            iBoundary_GB(:,:,:,iBlock) = yMaxBc_
       where( Xyz_DGB(z_,:,:,:,iBlock) < z1 ) &
            iBoundary_GB(:,:,:,iBlock) = zMinBc_
       where( Xyz_DGB(z_,:,:,:,iBlock) > z2 ) &
            iBoundary_GB(:,:,:,iBlock) = zMaxBc_
    end if

    true_cell(:,:,:,iBlock) = (iBoundary_GB(:,:,:,iBlock)==domain_)

    ! body_BLK: if any cell INCLUDING ghost cells is inside body(ies)
    body_BLK(iBlock) = .not. all(true_cell(:,:,:,iBlock))

    ! TRUE_BLK: if all cells EXCLUDING ghost cells are outside body(ies)
    true_BLK(iBlock) = all(true_cell(1:nI,1:nJ,1:nK,iBlock))

    if(DoTestMe)then
       write(*,*) NameSub,&
            ' finished with iBoundary_GB(iTest-nG:iTest+nG)=',&
            iBoundary_GB(iTest-nG:iTest+nG,jTest,kTest,iBlock)

       write(*,*) NameSub,&
            ' finished with true_cell(iTest-nG:iTest+nG)=',&
            true_cell(iTest-nG:iTest+nG,jTest,kTest,iBlock)
    end if

  end subroutine fix_block_geometry

  !============================================================================

  subroutine fix_geometry(DoSolveSolidIn)

    use BATL_lib, ONLY: nBlock

    logical, intent(in) :: DoSolveSolidIn

    integer :: iBlock
    !--------------------------------------------------------------------------

    do iBlock = 1, nBlock
       call fix_block_geometry(iBlock, DoSolveSolidIn)
    end do

  end subroutine fix_geometry

end module ModBoundaryGeometry

!=============================================================================

subroutine fix_boundary_ghost_cells

  ! Recalculate true_cell information in ghost cells if grid changed. 

  use ModBoundaryGeometry, ONLY: iBoundary_GB, domain_
  use ModMain, ONLY : nBlock, Unused_B, iNewGrid, iNewDecomposition, &
       BlkTest, iTest, jTest, kTest, iteration_number, nOrderProlong
  use ModGeometry, ONLY: true_cell, body_BLK
  !use ModProcMH, ONLY: iProc
  use BATL_lib, ONLY: message_pass_cell

  implicit none

  integer:: iBlock
  integer:: iGridHere = -1, iDecompositionHere = -1

  logical:: DoTest, DoTestMe
  character(len=*), parameter:: NameSub = 'fix_boundary_ghost_cells'
  !----------------------------------------------------------------------------
  if(iGridHere==iNewGrid .and. iDecompositionHere==iNewDecomposition) RETURN

  iGridHere = iNewGrid; iDecompositionHere = iNewDecomposition

  call set_oktest(NameSub, DoTest, DoTestMe)

  if(DoTestMe) write(*,*)NameSub,' starting with true_cell(i-2:i+2)=', &
       true_cell(iTest-2:iTest+2,jTest,kTest,BlkTest)

  ! DoResChangeOnly=true works as long as the ghost cells are correctly set
  ! away from resolution changes. This usually holds, but not if the 
  ! boundary cells are set based on the state variables read from a restart
  ! file that has no ghost cell information saved. This can only happen
  ! at the very beginning of a run when iteration_number == 0.

  if(nOrderProlong > 1)then
     call message_pass_cell(iBoundary_GB, &
          DoResChangeOnlyIn=iteration_number>0, NameOperatorIn='max')
  else
     call message_pass_cell(iBoundary_GB, &
          nProlongOrderIn=1, nCoarseLayerIn=2, &
          DoSendCornerIn=.true., DoRestrictFaceIn=.true., &
          DoResChangeOnlyIn=iteration_number>0, NameOperatorIn='max')
  end if

  if(DoTestMe) write(*,*) NameSub,': iBoundary_GB(i-2:i+2)=', &
       iBoundary_GB(iTest-2:iTest+2,jTest,kTest,BlkTest)

  do iBlock = 1, nBlock
     if(Unused_B(iBlock)) CYCLE

     true_cell(:,:,:,iBlock) = true_cell(:,:,:,iBlock)  &
          .and. (iBoundary_GB(:,:,:,iBlock) == domain_)

     body_BLK(iBlock) = .not. all(true_cell(:,:,:,iBlock))   
  end do

  if(DoTestMe) write(*,*) NameSub,' finished with true_cell(i-2:i+2)=', &
       true_cell(iTest-2:iTest+2,jTest,kTest,BlkTest)

end subroutine fix_boundary_ghost_cells
