module BATL_grid

  use BATL_size
  use BATL_tree
  use BATL_geometry, ONLY: IsCartesian, TypeGeometry

  implicit none

  private ! except

  public :: init_grid
  public :: clean_grid
  public :: create_grid
  public :: create_grid_block
  public :: test_grid

  real, public:: &
       CoordMin_D(MaxDim),      & ! Min gen. coords of domain
       CoordMax_D(MaxDim)         ! Max gen. coordinates of domain

  real, public, allocatable::   &
       CoordMin_DB(:,:),        & ! Min gen. coordinates of a block
       CoordMax_DB(:,:),        & ! Max gen. coordinates of a block
       CellSize_DB(:,:),        & ! Cell size in gen. coordinates
       CellFace_DB(:,:),        & ! Cell faces for Cartesian grids
       CellFace_DFB(:,:,:,:,:), & ! Cell faces for general grids
       CellVolume_B(:),         & ! Cell volume for Cartesian grids
       CellVolume_GB(:,:,:,:),  & ! Cell volume for general grids
       Xyz_DGB(:,:,:,:,:)         ! Cartesian cell centers coords

  ! Local variables

  logical :: DoInitializeGrid = .true.
  

contains
  !============================================================================
  subroutine init_grid(CoordMinIn_D, CoordMaxIn_D)

    real, intent(in):: CoordMinIn_D(nDim), CoordMaxIn_D(nDim)
    !-------------------------------------------------------------------------
    if(.not. DoInitializeGrid) RETURN

    DoInitializeGrid = .false.

    ! Make sure that the thickness is unity in the ignored dimensions
    CoordMin_D = -0.5
    CoordMax_D = +0.5
    CoordMin_D(1:nDim) = CoordMinIn_D
    CoordMax_D(1:nDim) = CoordMaxIn_D

    allocate(CoordMin_DB(MaxDim,MaxBlock))
    allocate(CoordMax_DB(MaxDim,MaxBlock))
    allocate(CellSize_DB(MaxDim,MaxBlock))

    allocate(CellFace_DB(MaxDim,MaxBlock))
    if(.not.IsCartesian) &
         allocate(CellFace_DFB(MaxDim,1:nI+1,1:nJ+1,1:nK+1,MaxBlock))

    allocate(CellVolume_B(MaxBlock))
    allocate(CellVolume_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(Xyz_DGB(MaxDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))

  end subroutine init_grid
  !===========================================================================
  subroutine clean_grid

    if(DoInitializeGrid) RETURN

    DoInitializeGrid = .true.

    deallocate(CoordMin_DB, CoordMax_DB, CellSize_DB, CellFace_DB, &
         CellVolume_B, Xyz_DGB)
    if(allocated(CellFace_DFB)) deallocate(CellFace_DFB)
    if(allocated(CellVolume_GB))deallocate(CellVolume_GB)

    CoordMin_D =  0.0
    CoordMax_D = -1.0

  end subroutine clean_grid

  !===========================================================================

  subroutine create_grid_block(iBlock, iNodeIn)

    ! Create geometrical information for block iBlock on the local PE

    integer, intent(in):: iBlock

    ! In case iNode_B is not set, iNodeIn can provide the node info
    integer, optional, intent(in):: iNodeIn

    character(len=*), parameter:: NameSub = 'create_grid_block'

    real :: PositionMin_D(MaxDim), PositionMax_D(MaxDim)
    integer :: iNode, i, j, k
    !----------------------------------------------------------------------
    if(present(iNodeIn))then
       iNode = iNodeIn
    else
       iNode = iNode_B(iBlock)
    end if
    call get_tree_position(iNode, PositionMin_D, PositionMax_D)

    CoordMin_DB(:,iBlock)= CoordMin_D + (CoordMax_D - CoordMin_D)*PositionMin_D
    CoordMax_DB(:,iBlock)= CoordMin_D + (CoordMax_D - CoordMin_D)*PositionMax_D

    CellSize_DB(:,iBlock) = (CoordMax_DB(:,iBlock) - CoordMin_DB(:,iBlock)) &
         / nIJK_D

    if(IsCartesian .or. TypeGeometry == 'rz')then
       ! For RZ geometry this is true in generalized coordinate sense only
       CellVolume_B(iBlock) = product(CellSize_DB(:,iBlock))
       CellFace_DB(:,iBlock) = CellVolume_B(iBlock) / CellSize_DB(:,iBlock)

       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          Xyz_DGB(:,i,j,k,iBlock) = CoordMin_DB(:,iBlock) + &
               ( (/i, j, k/) - 0.5 ) * CellSize_DB(:,iBlock)
       end do; end do; end do

       if(TypeGeometry == 'rz')then
          do j = MinJ, MaxJ
             CellVolume_GB(:,j,:,iBlock) = &
                  CellVolume_B(iBlock)*abs(Xyz_DGB(2,1,j,1,iBlock))
          end do
          do j = 1, nJ
             CellFace_DFB(1,:,j,1:nK,iBlock) = &
                  CellFace_DB(1,iBlock)*abs(Xyz_DGB(2,1,j,1,iBlock))
          end do
          do j = 1, nJ+1
             ! Could use node coordinate here !!!
             CellFace_DFB(2,1:nI,j,1:nK,iBlock) = CellFace_DB(2,iBlock) &
                  *0.5*sum(abs(Xyz_DGB(2,1,j-1:j,1,iBlock)))
          end do
       else
          ! Also useful for Cartesian to keep code simple
          CellVolume_GB(:,:,:,iBlock) = CellVolume_B(iBlock)
       end if
    else
       call CON_stop(NameSub//': '//TypeGeometry// &
            ' geometry is not yet implemented')
    end if

  end subroutine create_grid_block

  !===========================================================================
  subroutine create_grid

    integer:: iBlock
    !------------------------------------------------------------------------
    do iBlock = 1, nBlock
       if(Unused_B(iBlock))CYCLE
       call create_grid_block(iBlock)
    end do

  end subroutine create_grid
  !===========================================================================

  subroutine show_grid_block(iBlock)

    use BATL_mpi, ONLY: iProc

    integer, intent(in):: iBlock

    ! Show grid information for block iBlock

    character(len=*), parameter:: NameSub = 'show_grid_block'
    !------------------------------------------------------------------------
    if(Unused_B(iBlock))then
       write(*,*) NameSub//' WARNING unused block ',iBlock,' on proc',iProc
       RETURN
    end if
    write(*,*)'show_grid_block for iProc, iBlock=',iProc, iBlock
    write(*,*)'CoordMin  =', CoordMin_DB(:,iBlock)
    write(*,*)'CoordMax  =', CoordMax_DB(:,iBlock)
    write(*,*)'CellSize  =', CellSize_DB(:,iBlock)
    if(IsCartesian)then
       write(*,*)'CellFace  =', CellFace_DB(:,iBlock)
       write(*,*)'CellVolume=', CellVolume_B(iBlock)
    else
       write(*,*)'CellFace(1, 1, 1)  =', CellFace_DFB(1:nDim,1,1,1,iBlock)
       write(*,*)'CellVolume(1, 1, 1)=', CellVolume_GB(1,1,1,iBlock)
    end if
    write(*,*)'Xyz( 1, 1, 1)=', Xyz_DGB(:, 1, 1, 1,iBlock)
    write(*,*)'Xyz(nI, 1, 1)=', Xyz_DGB(:,nI, 1, 1,iBlock)
    write(*,*)'Xyz( 1,nJ, 1)=', Xyz_DGB(:, 1,nJ, 1,iBlock)
    write(*,*)'Xyz( 1, 1,nK)=', Xyz_DGB(:, 1, 1,nK,iBlock)
    write(*,*)'Xyz(nI,nJ,nK)=', Xyz_DGB(:,nI,nJ,nK,iBlock)

  end subroutine show_grid_block

  !===========================================================================

  subroutine show_grid

    use BATL_mpi, ONLY: iProc, nProc, barrier_mpi

    ! Show all blocks sequentially on all processors, ie. show_grid 
    ! must be called from all processors of the MPI communicator iComm!

    integer:: iBlock, iPe
    !------------------------------------------------------------------------

    call barrier_mpi
    do iPe = 0, nProc - 1
       if(iPe == iProc) then
          do iBlock = 1, nBlock
             if(Unused_B(iBlock)) CYCLE
             call show_grid_block(iBlock)
          end do
       end if
       call barrier_mpi
    end do

  end subroutine show_grid

  !===========================================================================

  subroutine test_grid

    use BATL_mpi, ONLY: iProc
    use BATL_geometry, ONLY: init_geometry
    use ModNumConst, ONLY: i_DD

    integer :: iBlock, nBlockAll, Int_D(MaxDim)

    integer, parameter:: MaxBlockTest            = 50
    integer, parameter:: nRootTest_D(MaxDim)     = (/3,2,1/)
    logical, parameter:: IsPeriodicTest_D(MaxDim)= (/.true., .true., .false./)
    real:: DomainMin_D(MaxDim) = (/ 3.0, 2.0, 1.0 /)
    real:: DomainMax_D(MaxDim) = (/ 9.0, 6.0, 4.0 /)

    real, parameter:: Tolerance = 1e-6

    integer:: i, j, k, Di, Dj, iDim
    real:: Radius
    real, allocatable:: CellVolumeCart_B(:), CellFaceCart_DB(:,:)


    logical:: DoTestMe
    character(len=*), parameter :: NameSub = 'test_grid'
    !-----------------------------------------------------------------------
    DoTestMe = iProc == 0

    if(DoTestMe) write(*,*)'Testing init_grid'
    if(DoTestMe) write(*,*)'nDimAmr, nIJK_D=', nDimAmr, nIJK_D

    ! Set Cartesian grid geometry before initializing tree and grid
    call init_geometry( IsPeriodicIn_D = IsPeriodicTest_D(1:nDim) )
    call init_tree(MaxBlockTest)
    call init_grid( DomainMin_D(1:nDim), DomainMax_D(1:nDim) )
    call set_tree_root( nRootTest_D(1:nDim))

    call refine_tree_node(3)
    call distribute_tree(.true.)
    if(DoTestMe) call show_tree('After distribute_tree')

    if(DoTestMe) write(*,*)'Testing create_grid'
    call create_grid

    call show_grid

    if(nDim == 2)then
       if(DoTestMe) write(*,*)'Testing create_grid in RZ geometry'

       ! Store Cartesian values for checking
       allocate(CellVolumeCart_B(MaxBlock), CellFaceCart_DB(MaxDim,MaxBlock))
       CellFaceCart_DB = CellFace_DB
       CellVolumeCart_B= CellVolume_B

       ! Clean Cartesian grid
       call clean_grid

       ! Initialize RZ grid
       call init_geometry(TypeGeometryIn = 'rz')
       call init_grid( DomainMin_D(1:nDim), DomainMax_D(1:nDim) )
       call create_grid
       call show_grid

       ! Check relative to Cartesian
       do iBlock = 1, nBlock
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             if(abs( CellVolume_GB(i,j,k,iBlock) &
                  - abs(Xyz_DGB(2,i,j,k,iBlock))*CellVolumeCart_B(iBlock)) &
                  < Tolerance) CYCLE
             write(*,*)NameSub,' ERROR: incorrect cell volume=', &
                  CellVolume_GB(i,j,k,iBlock),' should be', &
                  abs(Xyz_DGB(2,i,j,k,iBlock))*CellVolumeCart_B(iBlock), &
                  ' at i,j,k,iBlock,iProc=', i, j, k, iBlock, iProc
          end do; end do; end do
          do iDim = 1, nDim
             Di = i_DD(1,iDim); Dj = i_DD(2,iDim)
             do k = 1, nK; do j = 1, nJ+Dj; do i = 1, nI+Di
                Radius = 0.5*sum(abs(Xyz_DGB(2,i-Di:i,j-Dj:j,k,iBlock)))
                if(abs(CellFace_DFB(iDim,i,j,k,iBlock) - &
                     Radius*CellFaceCart_DB(iDim,iBlock)) &
                      < Tolerance) CYCLE
                write(*,*)NameSub,' ERROR: incorrect face area=', &
                     CellFace_DFB(iDim,i,j,k,iBlock),' should be', &
                     Radius*CellFaceCart_DB(iDim,iBlock), &
                     ' at iDim,i,j,k,iBlock,iProc=', &
                     iDim, i, j, k, iBlock, iProc

             end do; end do; end do
          end do
       end do
    end if

    if(DoTestMe) write(*,*)'Testing clean_grid'
    call clean_grid
    call clean_tree
    
  end subroutine test_grid

end module BATL_grid
