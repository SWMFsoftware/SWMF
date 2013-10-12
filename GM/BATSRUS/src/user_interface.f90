!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!This code is a copyright protected software (c) 2002- University of Michigan
! We provide external subroutine interfaces to ModUser routines 
! to avoid circular dependencies

!=====================================================================
subroutine user_set_boundary_cells(iBlock)

  use ModUser, ONLY: user_sub => user_set_boundary_cells
  implicit none

  integer,intent(in)::iBlock
  !-------------------------------------------------------------------
  call user_sub(iBlock)

end subroutine user_set_boundary_cells

!=====================================================================
subroutine user_set_face_boundary(VarsGhostFace_V)

  use ModAdvance, ONLY: nVar
  use ModUser, ONLY: user_sub => user_set_face_boundary
  implicit none

  real, intent(out):: VarsGhostFace_V(nVar)
  !-------------------------------------------------------------------
  call user_sub(VarsGhostFace_V)

end subroutine user_set_face_boundary

!=====================================================================
subroutine user_set_cell_boundary(iBlock, iSide, TypeBc, IsFound)

  use ModUser, ONLY: user_sub => user_set_cell_boundary
  implicit none

  integer,          intent(in)  :: iBlock, iSide
  character(len=*), intent(in)  :: TypeBc
  logical,          intent(out) :: IsFound
  !-------------------------------------------------------------------
  call user_sub(iBlock,iSide, TypeBc, IsFound)

end subroutine user_set_cell_boundary

!=====================================================================
subroutine user_initial_perturbation

  use ModUser, ONLY: user_sub => user_initial_perturbation
  implicit none
  !-------------------------------------------------------------------
  call user_sub

end subroutine user_initial_perturbation

!=====================================================================
subroutine user_set_ics(iBlock)

  use ModUser, ONLY: user_sub => user_set_ics
  implicit none

  integer, intent(in) :: iBlock
  !-------------------------------------------------------------------
  call user_sub(iBlock)

end subroutine user_set_ics

!=====================================================================
subroutine user_init_session

  use ModUser, ONLY: user_sub => user_init_session
  implicit none
  !-------------------------------------------------------------------
  call user_sub

end subroutine user_init_session

!=====================================================================
subroutine user_action(NameAction)

  use ModUser, ONLY: user_sub => user_action
  implicit none

  character(len=*), intent(in):: NameAction
  !-------------------------------------------------------------------
  call user_sub(NameAction)

end subroutine user_action
!=====================================================================
subroutine user_specify_refinement(iBlock, iArea, DoRefine)

  use ModUser, ONLY: user_sub => user_specify_refinement
  implicit none

  integer, intent(in) :: iBlock, iArea
  logical,intent(out) :: DoRefine

  !-------------------------------------------------------------------
  call user_sub(iBlock, iArea, DoRefine)

end subroutine user_specify_refinement

!=====================================================================
subroutine user_amr_criteria(iBlock, UserCriteria, TypeCriteria, IsFound)

  use ModUser, ONLY: user_sub => user_amr_criteria
  implicit none

  integer, intent(in)          :: iBlock
  real, intent(out)            :: UserCriteria
  character(len=*), intent(in) :: TypeCriteria
  logical, intent(inout)       :: IsFound
  !-------------------------------------------------------------------
  call user_sub(iBlock, UserCriteria, TypeCriteria, IsFound)

end subroutine user_amr_criteria

!=====================================================================
subroutine user_read_inputs

  use ModUser, ONLY: user_sub => user_read_inputs
  implicit none
  !-------------------------------------------------------------------
  call user_sub

end subroutine user_read_inputs

!=====================================================================
subroutine user_get_log_var(VarValue, TypeVar, Radius)

  use ModUser, ONLY: user_sub => user_get_log_var
  implicit none

  real, intent(out)           :: VarValue
  character(len=*), intent(in):: TypeVar
  real, intent(in)            :: Radius
  !-------------------------------------------------------------------
  call user_sub(VarValue, TypeVar, Radius)
  
end subroutine user_get_log_var

!====================================================================

subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
     PlotVar_G, PlotVarBody, UsePlotVarBody, &
     NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

  use ModSize, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
  use ModUser, ONLY: user_sub => user_set_plot_var
  implicit none

  integer,          intent(in)   :: iBlock
  character(len=*), intent(in)   :: NameVar
  logical,          intent(in)   :: IsDimensional
  real,             intent(out)  :: &
       PlotVar_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
  real,             intent(out)  :: PlotVarBody
  logical,          intent(out)  :: UsePlotVarBody
  character(len=*), intent(inout):: NameTecVar
  character(len=*), intent(inout):: NameTecUnit
  character(len=*), intent(inout):: NameIdlUnit
  logical,          intent(out)  :: IsFound

  !-------------------------------------------------------------------
  call user_sub(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

end subroutine user_set_plot_var

!====================================================================

subroutine user_calc_sources(iBlock)

  use ModUser, ONLY: user_sub => user_calc_sources
  implicit none

  integer, intent(in) :: iBlock
  !-------------------------------------------------------------------
  call user_sub(iBlock)

end subroutine user_calc_sources

!=====================================================================

subroutine user_init_point_implicit

  use ModUser, ONLY: user_sub => user_init_point_implicit
  implicit none
  !-------------------------------------------------------------------
  call user_sub

end subroutine user_init_point_implicit

!=====================================================================

subroutine user_get_b0(x, y, z, B0_D)

  use ModUser, ONLY: user_sub => user_get_b0
  implicit none

  real, intent(in) :: x, y, z
  real, intent(out):: B0_D(3)
  !-------------------------------------------------------------------
  call user_sub(x, y, z, B0_D)

end subroutine user_get_b0

!=====================================================================
subroutine user_update_states(iStage, iBlock)

  use ModUser, ONLY: user_sub => user_update_states
  implicit none
 
  integer,intent(in):: iStage, iBlock
  !-------------------------------------------------------------------
  call user_sub(iStage, iBlock)

end subroutine user_update_states

!=====================================================================
subroutine user_normalization

  use ModUser, ONLY: user_sub => user_normalization
  implicit none
  !-------------------------------------------------------------------
  call user_sub

end subroutine user_normalization

!=====================================================================
subroutine user_io_units

  use ModUser, ONLY: user_sub => user_io_units
  implicit none
  !-------------------------------------------------------------------
  call user_sub
  
end subroutine user_io_units

!=====================================================================
subroutine user_set_resistivity(iBlock, Eta_G)

  use ModSize, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
  use ModUser, ONLY: user_sub => user_set_resistivity
  implicit none

  integer, intent(in) :: iBlock
  real,    intent(out):: Eta_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK) 
  !-------------------------------------------------------------------
  call user_sub(iBlock, Eta_G)

end subroutine user_set_resistivity
!=====================================================================
subroutine user_material_te(State_V, TeOut)

  use ModVarIndexes, ONLY: nVar
  use ModUser, ONLY: user_material_properties

  real, intent(in) :: State_V(nVar)
  real, intent(out) :: TeOut
  !-------------------------------------------------------------------
  call user_material_properties(State_V, TeOut=TeOut)

end subroutine user_material_te
