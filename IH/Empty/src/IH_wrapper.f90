!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module IH_wrapper

  ! Wrapper for an "empty" Inner Heliosphere (IH) component

  implicit none

  private ! except

  ! CON wrapper
  public:: IH_set_param
  public:: IH_init_session
  public:: IH_run
  public:: IH_save_restart
  public:: IH_finalize

  ! Global buffer coupling
  public:: IH_get_for_global_buffer

  ! Coupling toolkit
  public:: IH_synchronize_refinement
  public:: IH_get_for_mh
  public:: IH_get_for_mh_with_xyz
  public:: IH_put_from_mh

  ! Coupling with SC
  public:: IH_set_buffer_grid
  public:: IH_set_buffer_grid_get_info
  public:: IH_save_global_buffer
  public:: IH_match_ibc

  ! Coupling with SP
  public:: IH_get_for_sp
  public:: IH_get_a_line_point

  ! Coupling with GM
  public:: IH_get_for_gm

contains
  !==========================================================================
  subroutine IH_set_param(CompInfo, TypeAction)

    use CON_comp_info

    character (len=*), parameter :: NameSub='IH_set_param'

    ! Arguments
    type(CompInfoType), intent(inout):: CompInfo   ! Information for this comp.
    character (len=*), intent(in)    :: TypeAction ! What to do
    !-------------------------------------------------------------------------
    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,&
            Use        =.false., &
            NameVersion='Empty', &
            Version    =0.0)

    case default
       call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')
    end select

  end subroutine IH_set_param

  !============================================================================

  subroutine IH_init_session(iSession, TimeSimulation)

    !INPUT PARAMETERS:
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='IH_init_session'

    call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')

  end subroutine IH_init_session

  !============================================================================

  subroutine IH_finalize(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='IH_finalize'

    call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')

  end subroutine IH_finalize

  !============================================================================

  subroutine IH_save_restart(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='IH_save_restart'

    call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')

  end subroutine IH_save_restart

  !============================================================================

  subroutine IH_run(TimeSimulation,TimeSimulationLimit)

    !INPUT/OUTPUT ARGUMENTS:
    real, intent(inout):: TimeSimulation   ! current time of component

    !INPUT ARGUMENTS:
    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    character(len=*), parameter :: NameSub='IH_run'

    call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')

  end subroutine IH_run

  !===============================================================

  subroutine IH_synchronize_refinement(iProc0,iCommUnion)

    integer, intent(in) ::iProc0,iCommUnion
    character(len=*), parameter :: NameSub='IH_synchronize_refinement'

    call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')

  end subroutine IH_synchronize_refinement

  !===============================================================

  subroutine IH_get_for_gm(&
       nPartial,iGetStart,Get,W,State_V,nVar,TimeCoupling)

    use CON_router, ONLY: IndexPtrType, WeightPtrType

    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::W
    real,dimension(nVar),intent(out)::State_V
    real,intent(in)::TimeCoupling

    character(len=*), parameter :: NameSub='IH_get_for_gm'

    call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')

  end subroutine IH_get_for_gm
  !===================================================================!
  subroutine IH_get_for_mh(&
       nPartial,iGetStart,Get,W,State_V,nVar)

    use CON_router, ONLY: IndexPtrType, WeightPtrType

    integer,intent(in)              ::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)   ::Get
    type(WeightPtrType),intent(in)  ::W
    real,dimension(nVar),intent(out)::State_V

    character(len=*), parameter :: NameSub='IH_get_for_mh'

    call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')
  end subroutine IH_get_for_mh
  !===================================================================!
  subroutine IH_get_for_mh_with_xyz(&
       nPartial,iGetStart,Get,W,State_V,nVar)

    use CON_router, ONLY: IndexPtrType, WeightPtrType

    integer,intent(in)               :: nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)    :: Get
    type(WeightPtrType),intent(in)   :: W
    real,dimension(nVar),intent(out) ::State_V

    character(len=*), parameter :: NameSub='IH_get_for_mh_with_xyz'

    call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')
  end subroutine IH_get_for_mh_with_xyz
  !===================================================================!
  subroutine IH_set_buffer_grid(DD)

    use CON_router, ONLY: DomainDecompositionType
    type(DomainDecompositionType), intent(out)::DD

    character(len=*), parameter :: NameSub='IH_set_buffer_grid'

    call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')
  end subroutine IH_set_buffer_grid
  !===================================================================!
  subroutine IH_get_for_sp(&
       nPartial,iGetStart,Get,W,State_V,nVar)
    use CON_router, ONLY: IndexPtrType, WeightPtrType

    !INPUT ARGUMENTS:
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::W
    real,dimension(nVar),intent(out)::State_V

    character(len=*), parameter :: NameSub='IH_get_for_sp'

    call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')
  end subroutine IH_get_for_sp
  !===================================================================!
  subroutine IH_get_a_line_point(&
       nPartial,iGetStart,Get,W,State_V,nVar)

    use CON_router, ONLY: IndexPtrType, WeightPtrType

    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::W
    real,dimension(nVar),intent(out)::State_V

    character(len=*), parameter :: NameSub='IH_get_a_line_point'

    call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')
  end subroutine IH_get_a_line_point
  !===================================================================!
  subroutine IH_put_from_mh(nPartial,&
       iPutStart,&
       Put,& 
       Weight,&
       DoAdd,&
       StateSI_V,&
       nVar)
    use CON_router,    ONLY: IndexPtrType, WeightPtrType

    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::Weight
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::StateSI_V

    character (len=*), parameter :: NameSub='IH_put_from_mh.f90'

    call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')
  end subroutine IH_put_from_mh
  !===================================================================!
  subroutine IH_match_ibc

    character(len=*), parameter :: NameSub='IH_match_ibc'

    call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')

  end subroutine IH_match_ibc
  !===============================================================
  subroutine IH_set_buffer_grid_get_info(CompID_, &
       nR, nPhi, nTheta, BufferMinMax_DI)

    integer, intent(in)     :: CompID_
    integer, intent(out)    :: nR, nPhi, nTheta
    real, intent(out)       :: BufferMinMax_DI(3,2)

    character(len=*), parameter :: NameSub = 'IH_set_buffer_grid_get_info'
    ! ---------------------------------------------------------------

    call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')

  end subroutine IH_set_buffer_grid_get_info
  !===============================================================
  subroutine IH_save_global_buffer(nVar, nR, nPhi, nTheta, BufferIn_VG)

    integer,intent(in) :: nVar, nR, nPhi, nTheta
    real,intent(in)    :: BufferIn_VG(nVar, nR, nPhi, nTheta)

    character(len=*), parameter :: NameSub = 'IH_save_global_buffer'
    !-------------------------------------------------------------

    call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')

  end subroutine IH_save_global_buffer
  !===============================================================
  subroutine IH_get_for_global_buffer(&
       nR, nPhi,nTheta, BufferMinMax_DI, &
       TimeCoupling, iCompSource, iCompTarget, Buffer_VG)

    ! Buffer size and limits
    integer,intent(in) :: nR, nPhi, nTheta
    real, intent(in)   :: TimeCoupling
    real, intent(in)   :: BufferMinMax_DI(3,2)
    integer,intent(in) :: iCompSource, iCompTarget

    ! State variables to be fiiled in all buffer grid points
    real, intent(out):: Buffer_VG(:,:,:,:)

    character (len=*), parameter :: NameSub='IH_get_for_buffer_grid'

    call CON_stop(NameSub//': IH_ERROR: empty version cannot be used!')

  end subroutine IH_get_for_global_buffer

end module IH_wrapper
