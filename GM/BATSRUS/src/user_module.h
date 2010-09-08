  implicit none  

  SAVE ! save all module variables

  private ! except  
  public :: VersionUserModule  
  public :: NameUserModule  
  public :: user_set_boundary_cells
  public :: user_face_bcs
  public :: user_set_outerbcs
  public :: user_initial_perturbation
  public :: user_set_ics
  public :: user_init_session
  public :: user_specify_refinement
  public :: user_amr_criteria
  public :: user_read_inputs
  public :: user_write_progress
  public :: user_get_log_var
  public :: user_set_plot_var
  public :: user_calc_sources
  public :: user_init_point_implicit
  public :: user_get_b0
  public :: user_update_states
  public :: user_io_units
  public :: user_normalization
  public :: user_set_resistivity
  public :: user_material_properties
  public :: user_calc_sources_adjoint
  public :: user_update_states_adjoint
