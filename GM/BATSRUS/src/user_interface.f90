!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

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
subroutine user_specify_region(iArea, iBlock, nValue, NameLocation, &
     IsInside, IsInside_I, Value_I)

  use ModUser, ONLY: user_sub => user_specify_region
  implicit none

  integer,   intent(in):: iArea        ! area index in BATL_region
  integer,   intent(in):: iBlock       ! block index
  integer,   intent(in):: nValue       ! number of output values
  character, intent(in):: NameLocation ! c, g, x, y, z, or n

  logical, optional, intent(out) :: IsInside
  logical, optional, intent(out) :: IsInside_I(nValue)
  real,    optional, intent(out) :: Value_I(nValue)

  character(len=*), parameter :: NameSub = 'user_specify_region'
  !-------------------------------------------------------------------
  call user_sub(iArea, iBlock, nValue, NameLocation, &
     IsInside, IsInside_I, Value_I)

end subroutine user_specify_region

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
  real, optional, intent(in)  :: Radius
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

  real, intent(in)   :: x, y, z
  real, intent(inout):: B0_D(3)
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
!============================================================================
subroutine user_material_properties(State_V, i, j, k, iBlock, iDir, &
     EinternalIn, TeIn, NatomicOut, AverageIonChargeOut, &
     EinternalOut, TeOut, PressureOut, &
     CvOut, GammaOut, HeatCondOut, IonHeatCondOut, TeTiRelaxOut, &
     OpacityPlanckOut_W, OpacityEmissionOut_W, OpacityRosselandOut_W, &
     PlanckOut_W)

  use ModVarIndexes, ONLY: nWave, nVar
  use ModUser, ONLY: user_sub => user_material_properties
  implicit none

  ! The State_V vector is in normalized units, all other physical
  ! quantities are in SI.
  !
  ! If the electron energy is used, then EinternalIn, EinternalOut,
  ! PressureOut, CvOut refer to the electron internal energies,
  ! electron pressure, and electron specific heat, respectively.
  ! Otherwise they refer to the total (electron + ion) internal energies,
  ! total (electron + ion) pressure, and the total specific heat.

  real, intent(in) :: State_V(nVar)
  integer, optional, intent(in):: i, j, k, iBlock, iDir  ! cell/face index
  real, optional, intent(in)  :: EinternalIn             ! [J/m^3]
  real, optional, intent(in)  :: TeIn                    ! [K]
  real, optional, intent(out) :: NatomicOut              ! [1/m^3]
  real, optional, intent(out) :: AverageIonChargeOut     ! dimensionless
  real, optional, intent(out) :: EinternalOut            ! [J/m^3]
  real, optional, intent(out) :: TeOut                   ! [K]
  real, optional, intent(out) :: PressureOut             ! [Pa]
  real, optional, intent(out) :: CvOut                   ! [J/(K*m^3)]
  real, optional, intent(out) :: GammaOut                ! dimensionless
  real, optional, intent(out) :: HeatCondOut             ! [J/(m*K*s)]
  real, optional, intent(out) :: IonHeatCondOut          ! [J/(m*K*s)]
  real, optional, intent(out) :: TeTiRelaxOut            ! [1/s]
  real, optional, intent(out) :: OpacityPlanckOut_W(nWave)    ! [1/m]
  real, optional, intent(out) :: OpacityEmissionOut_W(nWave)  ! [1/m]
  real, optional, intent(out) :: OpacityRosselandOut_W(nWave) ! [1/m]

  ! Multi-group specific interface. The variables are respectively:
  !  Group Planckian spectral energy density
  real, optional, intent(out) :: PlanckOut_W(nWave)      ! [J/m^3]
  !---------------------------------------------------------------------------
  call user_sub(State_V, i, j, k, iBlock, iDir, &
       EinternalIn, TeIn, NatomicOut, AverageIonChargeOut, &
       EinternalOut, TeOut, PressureOut, &
       CvOut, GammaOut, HeatCondOut, IonHeatCondOut, TeTiRelaxOut, &
       OpacityPlanckOut_W, OpacityEmissionOut_W, OpacityRosselandOut_W, &
       PlanckOut_W)

end subroutine user_material_properties
