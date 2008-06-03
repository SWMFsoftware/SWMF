!^CFG COPYRIGHT UM
module ModAMR
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

  real, dimension(:,:,:), allocatable :: refine_criteria_list

  logical, dimension(:,:), allocatable :: refine_list
  logical, dimension(:,:), allocatable :: coarsen_list

  integer, dimension(:,:), allocatable :: ListToCoarsen
  integer, dimension(:,:), allocatable :: ListToRefine
  integer, dimension(:,:), allocatable :: SortB
  integer, dimension(:,:), allocatable :: SortP
  real, dimension(:,:), allocatable :: SortC

  logical, dimension(:,:),     allocatable :: unusedBlock_BP

  !\
  ! Refinement criteria parameters
  !/
  integer            :: nRefineCrit, MaxTotalBlocks
  real               :: PercentCoarsen, PercentRefine
  character (len=20) :: RefineCrit(3), TypeTransient_I(3)

  ! Refine for criterion n only if it is above RefineCritMin_I(n)
  real               :: RefineCritMin_I(3) = (/-1.0,-1.0,-1.0/)
  ! Coarsen only if the rescaled (0.0 to 1.0) criterion is below CoarsenCritMax
  real               :: CoarsenCritMax     = 2.0

  !\
  ! Refinement parameters.
  !/
  integer :: initial_refine_levels, nRefineLevelIC, nRefineLevel
  integer :: dn_refine
  integer :: min_block_level, max_block_level
  real    :: min_cell_dx, max_cell_dx
  logical :: automatic_refinement, fix_body_level


  !\
  ! Variables controlling grid resolution in different areas
  !/

  integer, parameter :: MaxArea = 100, lNameArea = 20
  integer            :: nArea   = 0

  type AreaType
     character(len=lNameArea) :: Name
     real                     :: Resolution
     real, dimension(3)       :: Center_D,  Size_D
     real                     :: Radius1, Radius2
     logical                  :: DoRotate
     real, dimension(3,3)     :: Rotate_DD
  end type AreaType

  type(AreaType) :: Area_I(MaxArea)

end module ModAMR
