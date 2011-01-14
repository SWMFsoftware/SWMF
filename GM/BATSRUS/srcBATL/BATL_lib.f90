module BATL_lib

  ! Collection of all public methods and data that an application can access

  use BATL_size
  use BATL_mpi
  use BATL_tree
  use BATL_geometry
  use BATL_grid
  use BATL_amr
  use BATL_amr_criteria
  use BATL_pass_cell
  use BATL_pass_face

  implicit none

  private ! except

  ! Public methods and variables of this module
  public:: init_batl
  public:: clean_batl
  public:: init_grid_batl
  public:: regrid_batl

  logical, public:: IsBatlInitialized = .false.

  ! Inherited from BATL_size
  public:: MaxDim, nDim, nDimAmr, iDimAmr_D
  public:: MaxBlock, nBlock
  public:: nI, nJ, nK, nIJK, nIJK_D
  public:: MinI, MaxI, MinJ, MaxJ, MinK, MaxK

  ! Inherited from BATL_mpi
  public:: init_mpi, clean_mpi, barrier_mpi
  public:: iComm, nProc, iProc          

  ! Inherited from BATL_tree
  public:: nNodeUsed
  public:: Unused_B, Unused_BP
  public:: iNode_B, iMortonNode_A
  public:: DiLevelNei_IIIB, iNodeNei_IIIB
  public:: iStatusNew_A, Refine_, Coarsen_
  public:: iTree_IA, MinLevel_, MaxLevel_, Unset_
  public:: write_tree_file, read_tree_file

  ! Inherited from BATL_geometry
  public:: TypeGeometry, IsCartesian, IsRzGeometry, IsSpherical, &
       IsCylindrical, IsPeriodic_D, x_, y_, z_, r_

  ! Inherited from BATL_grid
  public:: CoordMin_D, CoordMax_D, CoordMin_DB, CoordMax_DB, CellSize_DB
  public:: Xyz_DGB
  public:: CellFace_DB, CellFace_DFB
  public:: CellVolume_B, CellVolume_GB
  public:: find_grid_block, interpolate_grid

  ! Inherited from BATL_amr
  public:: iAmrChange_B
  public:: AmrRemoved_, AmrUnchanged_, AmrMoved_, AmrRefined_, AmrCoarsened_

  ! Inherited from BATL_amr_criteria
  public:: set_amr_criteria, clean_amr_criteria, read_amr_criteria_param
  public:: AmrCrit_IB, nAmrCrit
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

    ! Initialize the block-adaptive tree and the domain. 
    !
    ! CoordMin_D and CoordMax_D are the domain boundaries. For non-Cartesian
    ! grids this is meant in generalized coordinates.
    ! 
    ! MaxBlockIn gives the maximum number of blocks per processor
    ! during the whole run.
    ! At the base level the domain is decomposed into the root blocks.
    !
    ! The optional nRootIn_D argument provides the number of root blocks 
    ! in each dimension. The default is a single root block at the base level.
    !
    ! The optional TypeGeometry describes the geometry of the coordinate
    ! system. Default is "Cartesian". Currently the only other option is
    ! RZ geometry. 
    ! 
    ! The optional IsPeriodicIn_D argument tells if a certain direction
    ! is periodic or not. 
    !
    ! At the completion of this subroutine, the domain and the block-tree 
    ! are initialized, the root blocks are distributed and their coordinates
    ! cell volumes, face areas, etc. are all set.
    !
    !-------------------------------------------------------------------------
    if(IsBatlInitialized) RETURN

    call init_tree(MaxBlockIn)
    call init_geometry(TypeGeometryIn, IsPeriodicIn_D)
    call init_grid(CoordMinIn_D, CoordMaxIn_D)
    call set_tree_root(nRootIn_D)
    call distribute_tree(.true.)
    call create_grid

    IsBatlInitialized = .true.

  end subroutine init_batl
  !============================================================================
  subroutine clean_batl

    call clean_amr_criteria
    ! Free up memory by cleaning the grid and tree data
    call clean_grid
    call clean_tree

    IsBatlInitialized = .false.

  end subroutine clean_batl

  !============================================================================

  subroutine init_grid_batl(DoRefine_B)

    logical, optional, intent(in):: DoRefine_B(MaxBlock)

    ! Initialize the grid. The grid does not contain any data at this point.
    ! The optional DoRefine_B argument allows initial refinement of the grid.
    ! After completion, the grid nodes are distributed over the processors,
    ! using Morton ordering, and the coordinates and related variables
    ! (cell volume, face areas, etc.) are all set.

    integer:: iBlock
    !------------------------------------------------------------------------

    if(present(DoRefine_B))then
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          if(DoRefine_B(iBlock)) iStatusNew_A(iNode_B(iBlock)) = Refine_
       end do
       call adapt_tree
    end if
    call distribute_tree(DoMove=.true.)
    call create_grid

  end subroutine init_grid_batl

  !============================================================================

  subroutine regrid_batl(nVar, State_VGB, Dt_B, DoRefine_B, DoCoarsen_B, &
       DoBalanceEachLevelIn, DoTestIn)

    integer, intent(in)   :: nVar                         ! number of variables
    real,    intent(inout):: &                            ! state variables
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    real, intent(inout), optional:: Dt_B(MaxBlock)        ! time step limit

    logical, intent(in), optional:: DoRefine_B(MaxBlock)  ! request to refine
    logical, intent(in), optional:: DoCoarsen_B(MaxBlock) ! request to coarsen
    logical, intent(in), optional:: DoBalanceEachLevelIn  ! balance per level?
    logical, intent(in), optional:: DoTestIn              ! print test info

    ! Refine, coarsen and load balance the blocks containing the nVar 
    ! state variables in State_VGB. Use second order accurate conservative
    ! restriction and prolongation operators. Load balance by Morton ordering.
    ! If DoBalanceEachLevelIn is true, balance each AMR level independently.
    !
    ! Refinement and coarsening is primarily based on the iStatusNew_A array 
    ! (available via the BATL_lib module) that is indexed by nodes 
    ! (so a processor can request refinement for a non-local block). 
    ! It should be set to the values Refine_ and Coarsen_. The AMR algorithm
    ! checks whether the requests can be done while keeping resolution
    ! changes at most a factor of 2 (proper nesting). Refinement requests
    ! take precedence over coarsening requests when multiple processors
    ! set the iStatusNew_A for the same node. 
    !
    ! The optional arguments DoRefine_B and DoCoarsen_B are provided 
    ! for convenience. They can be used to request refinement and 
    ! coarsening for local blocks only.
    !
    ! The optional Dt_B argument contains the time step limit per block. 
    ! This information is moved together with the block during load balance,
    ! The time step limit is divided by 2 for prolonged blocks, and 
    ! multiplied by 2 (and minimum is taken over the children) for coarsened
    ! blocks.
    !
    ! Note that this estimate of the time step limit for prolonged/restricted
    ! blocks is only an approximate solution, and for partial AMR
    ! it is not correct at all. One can always recalculate the time step
    ! limit for the coarsened and prolonged blocks.
    !
    ! The iAmrChange_B array provides information about the changes due to AMR
    ! in the local blocks. Its value is set to AmrRemoved_ for blocks that
    ! became unused, to AmrUnchanged_ for no change, to AmrMoved_ for blocks
    ! that moved from one processor to another, AmrRefined_ for blocks that
    ! were just refined, and AmrCoarsened_ for blocks that were just coarsened.

    logical:: DoBalanceEachLevel
    integer:: iBlock

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'regrid_batl'
    !------------------------------------------------------------------------
    DoTest = .false.
    if(present(DoTestIn)) DoTest = DoTestIn

    if(DoTest)write(*,*) NameSub,' starting with nVar=', nVar

    DoBalanceEachLevel = .false.
    if(present(DoBalanceEachLevelIn)) DoBalanceEachLevel = DoBalanceEachLevelIn

    if(present(DoCoarsen_B))then
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          if(DoCoarsen_B(iBlock)) iStatusNew_A(iNode_B(iBlock)) = Coarsen_
       end do
    end if

    if(present(DoRefine_B))then
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          if(DoRefine_B(iBlock)) iStatusNew_A(iNode_B(iBlock)) = Refine_
       end do
    end if

    ! Coarsen and refine the tree nodes
    if(DoTest)write(*,*) NameSub,' call adapt_tree'
    call adapt_tree

    ! Load balance the tree
    if(DoTest)write(*,*) NameSub, &
         ' call distribute_tree with DoBalanceEachLevel=', DoBalanceEachLevel
         
    if(DoBalanceEachLevel)then
       call distribute_tree(DoMove=.false., iTypeNode_A=iTree_IA(Level_,:)+1)
    else
       call distribute_tree(DoMove=.false.)
    end if

    ! Coarsen, refine and load balance the flow variables, and set Dt_B.
    if(DoTest)write(*,*) NameSub,' call do_amr'
    call do_amr(nVar, State_VGB, Dt_B, DoTestIn=DoTestIn)

    ! Finalize the tree information
    if(DoTest)write(*,*) NameSub,' call move_tree'
    call move_tree

    if(DoTest)write(*,*) NameSub,' finished'

  end subroutine regrid_batl

end module BATL_lib
