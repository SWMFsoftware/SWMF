module BATL_lib

  ! Collection of all public methods and data that an application can access

  use BATL_size
  use BATL_mpi
  use BATL_tree
  use BATL_geometry
  use BATL_grid
  use BATL_pass_cell
  use BATL_pass_face

  implicit none

  private ! except

  ! Public methods of this module
  public:: init_batl
  public:: clean_batl
  public:: init_grid_batl
  public:: regrid_batl

  ! Inherited from BATL_size
  public:: MaxDim, nDim, nDimAmr
  public:: MaxBlock, nBlock
  public:: nI, nJ, nK, nIJK, nIJK_D
  public:: MinI, MaxI, MinJ, MaxJ, MinK, MaxK

  ! Inherited from BATL_mpi
  public:: init_mpi, clean_mpi, barrier_mpi
  public:: iComm, nProc, iProc          

  ! Inherited from BATL_tree
  public:: nNodeUsed
  public:: Unused_B
  public:: iNode_B
  public:: DiLevelNei_IIIB, iNodeNei_IIIB
  public:: iStatusNew_A, Refine_, Coarsen_
  public:: iTree_IA, MinLevel_, MaxLevel_, Unset_

  ! Inherited from BATL_geometry
  public:: TypeGeometry, IsCartesian, IsRzGeometry, IsSpherical, &
       IsCylindrical, IsPeriodic_D, x_, y_, z_, r_

  ! Inherited from BATL_grid
  public:: CoordMin_D, CoordMax_D, CoordMin_DB, CoordMax_DB, CellSize_DB
  public:: Xyz_DGB
  public:: CellFace_DB, CellFace_DFB
  public:: CellVolume_B, CellVolume_GB

  ! Inherited from BATL_pass_cell
  public:: message_pass_cell

  ! Inherited from BATL_pass_face
  public:: message_pass_face
  public:: store_face_flux
  public:: apply_flux_correction

contains
  !============================================================================
  subroutine init_batl(&
       CoordMinIn_D, CoordMaxIn_D, MaxBlockIn, &
       TypeGeometryIn, IsPeriodicIn_D, nRootIn_D)

    integer, intent(in):: MaxBlockIn         ! max number of blocks/processor
    real,    intent(in):: CoordMinIn_D(nDim) ! min (gen) coordinates of domain
    real,    intent(in):: CoordMaxIn_D(nDim) ! max (gen) coordinates of domain

    ! Number of blocks per dimension on the tree root level
    integer,          optional, intent(in):: nRootIn_D(nDim)
    ! Grid geometry type (cartesian, spherical, etc.)
    character(len=*), optional, intent(in):: TypeGeometryIn
    ! Periodicity of grid boundaries per dimention
    logical,          optional, intent(in):: IsPeriodicIn_D(nDim)
    !-------------------------------------------------------------------------
    call init_tree(MaxBlockIn)
    call init_geometry(TypeGeometryIn, IsPeriodicIn_D)
    call init_grid(CoordMinIn_D, CoordMaxIn_D)
    call set_tree_root(nRootIn_D)
    call distribute_tree(.true.)
    call create_grid

  end subroutine init_batl
  !============================================================================
  subroutine clean_batl
    call clean_grid
    call clean_tree
  end subroutine clean_batl

  !============================================================================

  subroutine init_grid_batl

    use BATL_tree, ONLY: adapt_tree, distribute_tree
    use BATL_grid, ONLY: create_grid
    !------------------------------------------------------------------------

    call adapt_tree
    call distribute_tree(DoMove=.true.)
    call create_grid

  end subroutine init_grid_batl

  !============================================================================

  subroutine regrid_batl(nVar, State_VGB, DoBalanceEachLevelIn)

    use BATL_tree, ONLY: adapt_tree, distribute_tree, move_tree, &
         iTree_IA, Level_
    use BATL_amr,  ONLY: do_amr
    use BATL_grid, ONLY: create_grid

    integer, intent(in)   :: nVar
    real,    intent(inout):: &
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    logical, intent(in), optional:: DoBalanceEachLevelIn

    logical:: DoBalanceEachLevel
    !------------------------------------------------------------------------
    DoBalanceEachLevel = .false.
    if(present(DoBalanceEachLevelIn)) DoBalanceEachLevel = DoBalanceEachLevelIn

    ! Coarsen and refine the tree nodes
    call adapt_tree
    ! Load balance the tree
    if(DoBalanceEachLevel)then
       call distribute_tree(DoMove=.false.,iTypeNode_A=iTree_IA(Level_,:)+1)
    else
       call distribute_tree(DoMove=.false.)
    end if
    ! Coarsen, refine and load balance the flow variables
    call do_amr(nVar,State_VGB)

    ! Finalize the tree information
    call move_tree

  end subroutine regrid_batl

end module BATL_lib
