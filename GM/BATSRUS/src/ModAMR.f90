!^CFG COPYRIGHT UM
module ModAMR

  use BATL_lib, ONLY:  unusedBlock_BP => Unused_BP

  use ModCube
  implicit none
  save

  ! Array converting CHILD number and FACE index to SUBFACE index
  integer, dimension(8,6) :: child2subface, child2subedge
  data child2subface  / &
       2,0,0,1,3,0,0,4, &     !East face
       0,2,1,0,0,3,4,0, &     !West face
       2,4,3,1,0,0,0,0, &     !South face
       0,0,0,0,1,3,4,2, &     !North face
       0,0,2,1,3,4,0,0, &     !Bot face
       1,2,0,0,0,0,4,3/       !Top face

  ! Array converting CHILD number and FACE index to SUBEDGE index
  data child2subedge  / &
       1,0,0,1,2,0,0,2, &     !East face  Y edges
       0,1,1,0,0,2,2,0, &     !West face  Y edges
       2,2,1,1,0,0,0,0, &     !South face Z edgee
       0,0,0,0,1,1,2,2, &     !North face Z edges
       0,0,2,1,1,2,0,0, &     !Bot face   X edges
       1,2,0,0,0,0,2,1/       !Top face   X edges
  
  !\
  ! Local and global refinement/coarsening and neighbor parameters.
  !/
  integer, dimension(:,:), allocatable :: availableBLKs

  integer :: local_cube(8), local_cubeBLK(8)

  real, dimension(:,:,:), allocatable :: refine_criteria_IBP

  logical, dimension(:,:), allocatable :: refine_list
  logical, dimension(:,:), allocatable :: coarsen_list

  integer, dimension(:,:), allocatable :: ListToCoarsen
  integer, dimension(:,:), allocatable :: ListToRefine
  integer, dimension(:,:), allocatable :: SortB
  integer, dimension(:,:), allocatable :: SortP
  real, dimension(:,:), allocatable :: SortC

  !\
  ! Refinement criteria parameters
  !/
  integer            :: nAmrCriteria, MaxTotalBlocks
  real               :: PercentCoarsen, PercentRefine
  character (len=20),dimension(:), allocatable:: RefineCrit, TypeTransient_I
  real, dimension(:), allocatable :: CoarsenLimit_I, RefineLimit_I

  ! Refine for criterion n only if it is above RefineCritMin_I(n)
  real, dimension(:), allocatable :: RefineCritMin_I
  real,dimension(:,:), allocatable :: AmrCriteria_IB
  ! Coarsen only if the rescaled (0.0 to 1.0) criterion is below CoarsenCritMax
  real               :: CoarsenCritMax     = 2.0

  integer, dimension(:), allocatable :: SortIndex_I
  
  !\
  ! Refinement parameters.
  !/
  integer :: initial_refine_levels, nRefineLevelIC, nRefineLevel
  integer :: DnAmr = -1
  integer :: min_block_level, max_block_level
  real    :: min_cell_dx, max_cell_dx
  real    :: DtAmr = -1.0
  logical :: automatic_refinement, fix_body_level
  logical :: DoAmr = .false.
  logical :: DoProfileAmr = .false.
  ! Needed by amr_criteria
  integer :: nCritGeo = 0, nCritPhys = 0
  !\
  ! Variables controlling grid resolution in different areas
  !/

  integer, parameter :: MaxArea = 100, lNameArea = 20
  integer            :: nArea   = 0

  type AreaType
     character(len=lNameArea) :: Name
     real                     :: Resolution
     real, dimension(3)       :: Center_D,  Size_D
     real                     :: Radius1
     logical                  :: DoRotate
     real, dimension(3,3)     :: Rotate_DD
  end type AreaType

  type(AreaType) :: Area_I(MaxArea)

contains
  !============================================================================
  subroutine init_mod_amr(nCrit)

    use ModMain, ONLY: MaxBlock
    use ModProcMH, ONLY: nProc

    integer, intent(in) :: nCrit
    !-----------------------------------------------------------------------

    ! clean for each time we have new refinment criteia
    call clean_mod_amr

    allocate(RefineCrit(nCrit),&
         TypeTransient_I(nCrit),&
         CoarsenLimit_I(NCrit+1),&
         RefineLimit_I(NCrit+1),&
         RefineCritMin_I(nCrit),&
         AmrCriteria_IB(nCrit+1,MaxBlock),&
         refine_criteria_IBP(nCrit+1,MaxBlock,nProc),&
         SortIndex_I(nCrit+1))
    CoarsenLimit_I      = -1.0
    RefineLimit_I       = -1.0
    RefineCritMin_I     = -1.0
    refine_criteria_IBP = 0.00

  end subroutine init_mod_amr
  !============================================================================
  subroutine clean_mod_amr

    if(allocated(RefineCrit)) then
       deallocate(RefineCrit,&
            TypeTransient_I,&
            CoarsenLimit_I,&
            RefineLimit_I,&
            RefineCritMin_I,&
            AmrCriteria_IB,&
            refine_criteria_IBP,&
            SortIndex_I)
    end if
  end subroutine clean_mod_amr

end module ModAMR
