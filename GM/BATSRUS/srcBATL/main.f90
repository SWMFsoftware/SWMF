!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program BATL_main

  use BATL_amr,        ONLY: test_amr
  use BATL_mpi,        ONLY: init_mpi
  use BATL_tree,       ONLY: test_tree
  use BATL_geometry,   ONLY: test_geometry
  use BATL_grid,       ONLY: test_grid
  use BATL_pass_cell,  ONLY: test_pass_cell
  use BATL_pass_node,  ONLY: test_pass_node
  use BATL_pass_face,  ONLY: test_pass_face
  use BATL_amr_criteria, ONLY: test_amr_criteria
  use ModMpi

  implicit none

  integer:: iError
  !--------------------------------------------------------------------------
  call MPI_init(iError)
  call init_mpi(MPI_COMM_WORLD)

  call test_tree
  call test_geometry
  call test_grid
  call test_pass_cell
  call test_pass_face
  call test_pass_node
  call test_amr
  call test_amr_criteria

  call MPI_finalize(iError)
  
end program BATL_main
!=============================================================================

include 'external_routines.f90'
