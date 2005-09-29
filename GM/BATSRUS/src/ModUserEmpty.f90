!^CFG COPYRIGHT UM
!^CFG FILE USERFILES
!========================================================================
!========================================================================
!
! User routines let the user make customizations to BATSRUS without 
! having to make changes in the "kernel" of the code.
!
! This file contains empty user routines.  They should be "used" 
! (included) in the srcUser/ModUser*.f90 files for routines that the user 
! does not wish to implement.
!
! Please see the documentation, and the files ModUserDefault.f90 and 
! srcUser/ModUserExamples.f90 for information about what the different user
! subroutines do and how to implement them for your specific problem.
!
!========================================================================
!========================================================================


!========================================================================
module ModUserEmpty
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
  subroutine user_face_bcs(iFace,jFace,kFace,iBlock,iSide,iBoundary,&
       iter,time_now, FaceCoords_D,&
       VarsTrueFace_V,VarsGhostFace_V,&
       B0Face_D,  UseIonosphereHere, UseCorotationHere)
    use ModSize, ONLY: nDim
    use ModAdvance, ONLY: nFaceValueVars

    integer,intent(in)::iFace,jFace,kFace,iBlock,iSide,iBoundary,iter
    real,intent(in)::time_now
    real,dimension(nDim),intent(in)::FaceCoords_D,B0Face_D
    real,dimension(nFaceValueVars),intent(in)::VarsTrueFace_V
    real,dimension(nFaceValueVars),intent(out)::VarsGhostFace_V
    logical,intent(in)::UseIonosphereHere,UseCorotationHere

    character (len=*), parameter :: Name='user_face_bcs' 
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_face_bcs

  !=====================================================================
  subroutine user_set_outerBCs(iBlock,iSide, TypeBc,found)

    integer,intent(in)::iBlock, iSide
    logical,intent(out) :: found
    character (len=20),intent(in) :: TypeBc

    character (len=*), parameter :: Name='user_set_outerBCs'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_set_outerBCs

  !=====================================================================
  subroutine user_initial_perturbation

    character (len=*), parameter :: Name='user_initial_perturbation'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_initial_perturbation

  !=====================================================================
  subroutine user_set_ICs

    character (len=*), parameter :: Name='user_set_ICs'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_set_ICs

  !=====================================================================
  subroutine user_specify_initial_refinement(iBLK,refineBlock,lev,dxBlock,   &
       xCenter,yCenter,zCenter,rCenter,                        &
       minx,miny,minz,minR,maxx,maxy,maxz,maxR,found)

    logical,intent(out) :: refineBlock, found
    integer, intent(in) :: lev
    real, intent(in) :: dxBlock
    real, intent(in) :: xCenter,yCenter,zCenter,rCenter
    real, intent(in) :: minx,miny,minz,minR
    real, intent(in) :: maxx,maxy,maxz,maxR
    integer,intent(in) :: iBLK

    character (len=*), parameter :: Name='user_specify_initial_refinement'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_specify_initial_refinement

  !=====================================================================
  subroutine user_amr_criteria(iBLK, userCriteria, TypeCriteria, found)

    logical ,intent(inout):: found
    integer, intent(in) :: iBLK
    real, intent(out) :: userCriteria
    character (len=20),intent(in) :: TypeCriteria

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
  subroutine user_get_log_var(VarValue,TypeVar)

    real, intent(out):: VarValue
    character (LEN=10), intent(in):: TypeVar

    character (len=*), parameter :: Name='user_get_log_var'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_get_log_var

  !=====================================================================
  subroutine user_calc_sources

    character (len=*), parameter :: Name='user_calc_sources'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_calc_sources

  !=====================================================================
  subroutine user_heat_source

    character (len=*), parameter :: Name='user_heat_source'
    !-------------------------------------------------------------------
    call stop_user(Name)
  end subroutine user_heat_source

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
  subroutine stop_user(Name)
    ! Note that this routine is not a user routine but just a routine
    ! which warns the user if they try to use an unimplemented user routine.

    character (len=*), intent(in) :: Name
    !-------------------------------------------------------------------
    call stop_mpi('You are trying to call a user supplied version of '//   &
         'the routine '//name//'but you have implmented only the '//       &
         'empty version. Please implement the routine in your user file'// &
         'located in srcUser/ModUser*.f90')
  end subroutine stop_user

end module ModUserEmpty
