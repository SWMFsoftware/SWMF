!^CFG COPYRIGHT UM
!==============================================================================
module ModUserEmpty

  ! This module contains empty user routines.  They should be "used" 
  ! (included) in the srcUser/ModUser*.f90 files for routines that the user 
  ! does not wish to implement.

  implicit none

  private :: stop_user

contains

  !=====================================================================
  subroutine user_set_boundary_cells(iBLK)

    integer,intent(in)::iBLK

    character (len=*), parameter :: NameSub = 'user_set_boundary_cells'
    !-------------------------------------------------------------------
    call stop_user(NameSub)
  end subroutine user_set_boundary_cells

  !=====================================================================
  subroutine user_face_bcs(VarsGhostFace_V)

    use ModAdvance, ONLY: nVar

    real, intent(out):: VarsGhostFace_V(nVar)

    character (len=*), parameter :: NameSub = 'user_face_bcs' 
    !-------------------------------------------------------------------
    call stop_user(NameSub)
  end subroutine user_face_bcs

  !=====================================================================
  subroutine user_set_outerbcs(iBlock,iSide, TypeBc, IsFound)

    integer,          intent(in)  :: iBlock, iSide
    character(len=20),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    character (len=*), parameter :: NameSub = 'user_set_outerbcs'
    !-------------------------------------------------------------------
    IsFound = .false.
    call stop_user(NameSub)
  end subroutine user_set_outerbcs

  !=====================================================================
  subroutine user_initial_perturbation
    use ModMain,ONLY: nBlockMax
    character (len=*), parameter :: NameSub = 'user_initial_perturbation'
    integer::iBlock
    !-------------------------------------------------------------------
    !The routine is called once and should be applied for all blocks, the 
    !do-loop should be present. Another distinction from user_set_ics is that 
    !user_initial_perturbation can be applied after restart, while
    !user_set_ICs cannot.

    do iBlock = 1, nBlockMax
    
       call stop_user(NameSub)
    end do
  end subroutine user_initial_perturbation

  !=====================================================================
  subroutine user_set_ics
    use ModMain,ONLY: globalBLK

    character (len=*), parameter :: NameSub = 'user_set_ics'
    integer::iBlock
    !-------------------------------------------------------------------
    !The routine is called for each block, the number of block should
    !be passed to the routine using globalBLK

    iBlock = globalBLK  

    call stop_user(NameSub)
  end subroutine user_set_ics

  !=====================================================================
  subroutine user_init_session

    character (len=*), parameter :: NameSub = 'user_init_session'
    !-------------------------------------------------------------------
    call stop_user(NameSub)
  end subroutine user_init_session

  !=====================================================================
  subroutine user_specify_refinement(iBlock, iArea, DoRefine)

    integer, intent(in) :: iBlock, iArea
    logical,intent(out) :: DoRefine

    character (len=*), parameter :: NameSub = 'user_specify_refinement'
    !-------------------------------------------------------------------
    call stop_user(NameSub)
  end subroutine user_specify_refinement

  !=====================================================================
  subroutine user_amr_criteria(iBLK, UserCriteria, TypeCriteria, IsFound)

    integer, intent(in)          :: iBLK
    real, intent(out)            :: userCriteria
    character (len=*),intent(in) :: TypeCriteria
    logical, intent(inout)       :: IsFound

    character (len=*), parameter :: NameSub = 'user_amr_criteria'
    !-------------------------------------------------------------------
    call stop_user(NameSub)
  end subroutine user_amr_criteria

  !=====================================================================
  subroutine user_read_inputs

    character (len=*), parameter :: NameSub = 'user_read_inputs'
    !-------------------------------------------------------------------
    call stop_user(NameSub)
  end subroutine user_read_inputs

  !=====================================================================
  subroutine user_write_progress

    character (len=*), parameter :: NameSub = 'user_write_progress'
    !-------------------------------------------------------------------
    call stop_user(NameSub)
  end subroutine user_write_progress

  !=====================================================================
  subroutine user_get_log_var(VarValue, TypeVar, Radius)

    real, intent(out)            :: VarValue
    character (len=*), intent(in):: TypeVar
    real, intent(in), optional :: Radius

    character (len=*), parameter :: NameSub = 'user_get_log_var'
    !-------------------------------------------------------------------
    call stop_user(NameSub)
  end subroutine user_get_log_var

  !====================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModSize, ONLY: nI, nJ, nK

    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(-1:nI+2, -1:nJ+2, -1:nK+2)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    !-------------------------------------------------------------------
    call stop_user(NameSub)

  end subroutine user_set_plot_var

  !====================================================================

  subroutine user_calc_sources

    character (len=*), parameter :: NameSub = 'user_calc_sources'
    !-------------------------------------------------------------------
    call stop_user(NameSub)
  end subroutine user_calc_sources

  !=====================================================================

  subroutine user_init_point_implicit

    character (len=*), parameter :: NameSub = 'user_init_point_implicit'
    !-------------------------------------------------------------------
    call stop_user(NameSub)

  end subroutine user_init_point_implicit

  !=====================================================================

  subroutine user_get_b0(xx,yy,zz,B0)

    real, intent(in):: xx,yy,zz
    real, intent(out), dimension(3):: B0

    character (len=*), parameter :: NameSub = 'user_get_b0'
    !-------------------------------------------------------------------
    call stop_user(NameSub)
  end subroutine user_get_b0

  !=====================================================================
  subroutine user_update_states(iStage,iBlock)

    integer,intent(in)::iStage,iBlock

    character (len=*), parameter :: NameSub = 'user_update_states'
    !-------------------------------------------------------------------
    call stop_user(NameSub)
  end subroutine user_update_states

  !=====================================================================
  subroutine user_normalization
    use ModPhysics 
    
    character (len=*), parameter :: NameSub = 'user_normalization'
    !-------------------------------------------------------------------
    call stop_user(NameSub)
  end subroutine user_normalization

  !=====================================================================
  subroutine user_io_units
    use ModPhysics 
    
    character (len=*), parameter :: NameSub = 'user_io_units'
    !-------------------------------------------------------------------
    call stop_user(NameSub)
  end subroutine user_io_units

  !=====================================================================
  subroutine user_set_resistivity(iBlock, Eta_G)
    ! This subrountine set the eta for every block
    use ModSize

    integer, intent(in) :: iBlock
    real,    intent(out):: Eta_G(-1:nI+2,-1:nJ+2,-1:nK+2) 
    character (len=*), parameter :: NameSub = 'user_set_resistivity'

    !-------------------------------------------------------------------
    call stop_mpi(NameSub)

  end subroutine user_set_resistivity

  !===========================================================================

  subroutine user_material_properties(State_V, EinternalSiIn, &
       TeSiIn, EinternalSiOut, TeSiOut, PressureSiOut, CvSiOut, GammaOut, &
       AbsorptionOpacitySiOut, RosselandMeanOpacitySiOut, &
       HeatConductionCoefSiOut)

    ! The State_V vector is in normalized units

    use ModVarIndexes, ONLY: nVar

    real, intent(in) :: State_V(nVar)
    real, optional, intent(in)  :: EinternalSiIn             ! [J/m^3]
    real, optional, intent(in)  :: TeSiIn                    ! [K]
    real, optional, intent(out) :: EinternalSiOut            ! [J/m^3]
    real, optional, intent(out) :: TeSiOut                   ! [K]
    real, optional, intent(out) :: PressureSiOut             ! [Pa]
    real, optional, intent(out) :: CvSiOut                   ! [J/(K*m^3)]
    real, optional, intent(out) :: GammaOut                  ! dimensionless
    real, optional, intent(out) :: AbsorptionOpacitySiOut    ! [1/m]
    real, optional, intent(out) :: RosselandMeanOpacitySiOut ! [1/m]
    real, optional, intent(out) :: HeatConductionCoefSiOut   ! [Jm^2/(Ks)]

    character (len=*), parameter :: NameSub = 'user_material_properties'
    !------------------------------------------------------------------------
    call stop_mpi(NameSub)

  end subroutine user_material_properties

  !=====================================================================
  subroutine stop_user(NameSub)
    ! Note that this routine is not a user routine but just a routine
    ! which warns the user if they try to use an unimplemented user routine.

    character (len=*), intent(in) :: NameSub
    !-------------------------------------------------------------------
    call stop_mpi('You are trying to call the empty user routine '//   &
         NameSub//'. Please implement the routine in src/ModUser.f90')
  end subroutine stop_user


end module ModUserEmpty
!==============================================================================
