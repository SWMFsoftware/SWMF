!^CFG COPYRIGHT UM
module ModAMR

  use BATL_lib, ONLY:  unusedBlock_BP => Unused_BP

  implicit none
  save

  ! Local and global refinement/coarsening and neighbor parameters.

  ! Refinement criteria parameters
  integer          :: nAmrCriteria, MaxTotalBlocks
  real             :: PercentCoarsen, PercentRefine
  real, allocatable:: CoarsenLimit_I(:), RefineLimit_I(:)
  character(len=20), allocatable:: RefineCrit(:), TypeTransient_I(:)

  real, allocatable:: AmrCriteria_IB(:,:)

  ! Refinement parameters.
  integer :: initial_refine_levels, nRefineLevelIC, nRefineLevel
  integer :: DnAmr = -1
  integer :: min_block_level, max_block_level
  real    :: min_cell_dx, max_cell_dx
  real    :: DtAmr = -1.0
  logical :: automatic_refinement
  logical :: DoAmr = .false.
  logical :: DoProfileAmr = .false.

  ! Needed by amr_criteria
  integer :: nCritGeo = 0, nCritPhys = 0

contains
  !============================================================================
  subroutine init_mod_amr(nCrit)

    use ModMain, ONLY: MaxBlock

    integer, intent(in) :: nCrit
    !-----------------------------------------------------------------------

    ! clean for each time we have new refinment criteia
    call clean_mod_amr

    allocate(RefineCrit(nCrit),              &
         TypeTransient_I(nCrit),             &
         CoarsenLimit_I(nCrit+1),            &
         RefineLimit_I(nCrit+1),             &
         AmrCriteria_IB(nCrit+1,MaxBlock))

    CoarsenLimit_I      = -1.0
    RefineLimit_I       = -1.0

  end subroutine init_mod_amr
  !============================================================================
  subroutine clean_mod_amr

    if(allocated(RefineCrit)) deallocate( &
         RefineCrit,&
         TypeTransient_I,&
         CoarsenLimit_I,&
         RefineLimit_I,&
         AmrCriteria_IB)

  end subroutine clean_mod_amr

end module ModAMR
