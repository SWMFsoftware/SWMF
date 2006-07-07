!^CFG COPYRIGHT UM
!========================================================================
module ModUser
  ! This is the default user module which contains empty methods defined
  ! in ModUserEmpty.f90
  !
  ! Please see the documentation, and the files ModUserEmpty.f90 and 
  ! srcUser/ModUserExamples.f90 for information about what the different user
  ! subroutines do and how to implement them for your specific problem.

  use ModUserEmpty, ONLY:               &
       user_read_inputs,                &
       user_init_session,               &
       user_set_ics,                    &
       user_initial_perturbation,       &
       user_set_boundary_cells,         &
       user_face_bcs,                   &
       user_set_outerbcs,               &
       user_specify_initial_refinement, &
       user_amr_criteria,               &
       user_write_progress,             &
       user_get_log_var,                &
       user_set_plot_var,               &
       user_calc_sources,               &
       user_heat_source,                &
       user_get_b0,                     &
       user_update_states

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = 'DEFAULT EMPTY ROUTINES'

end module ModUser
