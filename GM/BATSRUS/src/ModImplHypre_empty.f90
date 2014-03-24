!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModImplHypre

  ! Interface to the EMPTY version of HYPRE.
  ! This module is needed when HYPRE is not configured/available.

  implicit none

  private
  public:: hypre_read_param
  public:: hypre_initialize
  public:: hypre_finalize
  public:: hypre_set_matrix_block
  public:: hypre_set_matrix
  public:: hypre_preconditioner

  logical, public, parameter:: IsHypreAvailable = .false. ! signals empty

contains

  !==========================================================================
  subroutine hypre_read_param

    call stop_mpi('hypre_read_param: empty! Use Config.pl -hypre')

  end subroutine hypre_read_param

  !==========================================================================
  subroutine hypre_initialize

    call stop_mpi('hypre_initialize: empty! Use Config.pl -hypre')

  end subroutine hypre_initialize

  !==========================================================================
  subroutine hypre_finalize

    call stop_mpi('hypre_finalize: empty! Use Config.pl -hypre')

  end subroutine hypre_finalize

  !===========================================================================

  subroutine hypre_set_matrix_block(iImplBlock, Jacobian_CI)

    integer, intent(in):: iImplBlock
    real,    intent(inout):: Jacobian_CI(1,1,1,1)
    !------------------------------------------------------------------------
    call stop_mpi('hypre_set_matrix_block: empty! Use Config.pl -hypre')

  end subroutine hypre_set_matrix_block

  !============================================================================

  subroutine hypre_set_matrix(DoInitHypreAmg)

    logical, intent(in):: DoInitHypreAmg
    !------------------------------------------------------------------------
    call stop_mpi('hypre_set_matrix: empty! Use Config.pl -hypre')

  end subroutine hypre_set_matrix

  !============================================================================

  subroutine hypre_preconditioner(n, y_I)

    integer, intent(in):: n
    real, intent(inout):: y_I(n)
    !-------------------------------------------------------------------------

    call stop_mpi('hypre_preconditioner: empty! Use Config.pl -hypre')

  end subroutine hypre_preconditioner

end module ModImplHypre
