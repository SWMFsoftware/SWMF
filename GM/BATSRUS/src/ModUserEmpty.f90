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

    character (len=*), parameter :: Name='user_set_boundary_cells'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_set_boundary_cells

  !=====================================================================
  subroutine user_face_bcs(iFace, jFace, kFace, iBlock, iSide, iBoundary,&
       iter, time_now, FaceCoords_D, &
       VarsTrueFace_V, VarsGhostFace_V, &
       B0Face_D,  UseIonosphereHere, UseCorotationHere)
    use ModSize, ONLY: nDim
    use ModAdvance, ONLY: nFaceValueVars

    integer,intent(in)::iFace,jFace,kFace,iBlock,iSide,iBoundary,iter
    real, intent(in):: time_now
    real, dimension(nDim), intent(in) :: FaceCoords_D, B0Face_D
    real, dimension(nFaceValueVars), intent(in) :: VarsTrueFace_V
    real, dimension(nFaceValueVars), intent(out):: VarsGhostFace_V
    logical, intent(in):: UseIonosphereHere, UseCorotationHere

    character (len=*), parameter :: Name='user_face_bcs' 
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_face_bcs

  !=====================================================================
  subroutine user_set_outerbcs(iBlock,iSide, TypeBc,found)

    integer,intent(in)::iBlock, iSide
    logical,intent(out) :: found
    character (len=20),intent(in) :: TypeBc

    character (len=*), parameter :: Name='user_set_outerbcs'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_set_outerbcs

  !=====================================================================
  subroutine user_initial_perturbation

    character (len=*), parameter :: Name='user_initial_perturbation'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_initial_perturbation

  !=====================================================================
  subroutine user_set_ics

    character (len=*), parameter :: Name='user_set_ics'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_set_ics

  !=====================================================================
  subroutine user_init_session

    character (len=*), parameter :: Name='user_init_session'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_init_session

  !=====================================================================
  subroutine user_specify_initial_refinement(iBLK,refineBlock,lev,DxBlock, &
       xCenter,yCenter,zCenter,rCenter,                        &
       minx,miny,minz,minR,maxx,maxy,maxz,maxR,found)

    logical,intent(out) :: refineBlock, found
    integer, intent(in) :: lev
    real, intent(in)    :: DxBlock
    real, intent(in)    :: xCenter,yCenter,zCenter,rCenter
    real, intent(in)    :: minx,miny,minz,minR
    real, intent(in)    :: maxx,maxy,maxz,maxR
    integer, intent(in) :: iBLK

    character (len=*), parameter :: Name='user_specify_initial_refinement'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_specify_initial_refinement

  !=====================================================================
  subroutine user_amr_criteria(iBLK, UserCriteria, TypeCriteria, IsFound)

    integer, intent(in)          :: iBLK
    real, intent(out)            :: userCriteria
    character (len=*),intent(in) :: TypeCriteria
    logical, intent(inout)       :: IsFound

    character (len=*), parameter :: Name='user_amr_criteria'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_amr_criteria

  !=====================================================================
  subroutine user_read_inputs

    character (len=*), parameter :: Name='user_read_inputs'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_read_inputs

  !=====================================================================
  subroutine user_write_progress

    character (len=*), parameter :: Name='user_write_progress'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_write_progress

  !=====================================================================
  subroutine user_get_log_var(VarValue, TypeVar)

    real, intent(out)            :: VarValue
    character (len=*), intent(in):: TypeVar

    character (len=*), parameter :: Name='user_get_log_var'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_get_log_var

  !====================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModSize, ONLY: nI, nJ, nK

    integer,          intent(in) :: iBlock
    character(len=*), intent(in) :: NameVar
    logical,          intent(in) :: IsDimensional
    real,             intent(out):: PlotVar_G(-1:nI+2, -1:nJ+2, -1:nK+2)
    real,             intent(out):: PlotVarBody
    logical,          intent(out):: UsePlotVarBody
    character(len=*), intent(out):: NameTecVar
    character(len=*), intent(out):: NameTecUnit
    character(len=*), intent(out):: NameIdlUnit
    logical,          intent(out):: IsFound

    character (len=*), parameter :: Name='user_set_plot_var'
    !-------------------------------------------------------------------
    call stop_user(Name)

  end subroutine user_set_plot_var

  !====================================================================

  subroutine user_calc_sources

    character (len=*), parameter :: Name='user_calc_sources'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_calc_sources

  !=====================================================================
  subroutine user_get_b0(xx,yy,zz,B0)

    real, intent(in):: xx,yy,zz
    real, intent(out), dimension(3):: B0

    character (len=*), parameter :: Name='user_get_b0'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_get_b0

  !=====================================================================
  subroutine user_update_states(iStage,iBlock)

    integer,intent(in)::iStage,iBlock

    character (len=*), parameter :: Name='user_update_states'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_update_states

  !=====================================================================
  subroutine user_io_units
    use ModPhysics 
    
    character (len=*), parameter :: Name='user_io_units'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_io_units

  !=====================================================================
  subroutine stop_user(Name)
    ! Note that this routine is not a user routine but just a routine
    ! which warns the user if they try to use an unimplemented user routine.

    character (len=*), intent(in) :: Name
    !-------------------------------------------------------------------
    call stop_mpi('You are trying to call the empty user routine '//   &
         Name//'. Please implement the routine in src/ModUser.f90')
  end subroutine stop_user

end module ModUserEmpty
!==============================================================================
