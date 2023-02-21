!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module UA_wrapper

  ! Wrapper for Upper Atmosphere (UA) component

  use ModUtilities, ONLY: CON_stop

  implicit none

  private ! except

  public:: UA_set_param
  public:: UA_init_session
  public:: UA_run
  public:: UA_save_restart
  public:: UA_finalize

  ! IE Coupler:
  public :: UA_get_info_for_ie
  public :: UA_get_for_ie
  public :: UA_put_from_ie

  ! GM Coupler
  public :: UA_find_points
  public :: UA_get_for_gm
  public :: UA_get_grid_info
  
contains
  !============================================================================
  subroutine UA_set_param(CompInfo, TypeAction)

    use CON_comp_info

    character (len=*), parameter :: NameSub='UA_set_param'

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
       call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')
    end select

  end subroutine UA_set_param
  !============================================================================
  subroutine UA_init_session(iSession, TimeSimulation)

    !INPUT PARAMETERS:
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='UA_init_session'

    call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')

  end subroutine UA_init_session
  !============================================================================
  subroutine UA_finalize(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='UA_finalize'

    call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')

  end subroutine UA_finalize
  !============================================================================
  subroutine UA_save_restart(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='UA_save_restart'

    call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')

  end subroutine UA_save_restart
  !============================================================================
  subroutine UA_run(TimeSimulation,TimeSimulationLimit)

    !INPUT/OUTPUT ARGUMENTS:
    real, intent(inout):: TimeSimulation   ! current time of component

    !INPUT ARGUMENTS:
    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    character(len=*), parameter :: NameSub='UA_run'

    call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')

  end subroutine UA_run

  !============================================================================
  subroutine UA_get_info_for_ie(nVar, NameVar_V, nMagLat, nMagLon)

    !OUTPUT ARGUMENTS:
    integer, intent(out) :: nVar
    integer, intent(out), optional :: nMagLat, nMagLon
    character(len=*), intent(out), optional :: NameVar_V(:)

    character(len=*), parameter :: NameSub='UA_get_info_for_ie'

    call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')

  end subroutine UA_get_info_for_ie

  !============================================================================
  subroutine UA_put_from_ie(Buffer_IIV, iSizeIn, jSizeIn, nVarIn, &
       NameVarIn_V, iBlock)

    !INPUT/OUTPUT ARGUMENTS:
    integer, intent(in)           :: iSizeIn, jSizeIn, nVarIn, iBlock
    real, intent(in)              :: Buffer_IIV(iSizeIn,jSizeIn,nVarIn)
    character (len=*),intent(in)  :: NameVarIn_V(nVarIn)

    character (len=*), parameter :: NameSub='UA_put_from_ie'

    call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')

  end subroutine UA_put_from_ie

  !============================================================================
  subroutine UA_get_for_ie(BufferOut_IIBV, nMltIn, nLatIn, nVarIn, NameVarIn_V)

    ! INPUT ARGUMENTS:
    integer,          intent(in) :: nMltIn, nLatIn, nVarIn
    character(len=3), intent(in) :: NameVarIn_V(nVarIn)

    ! OUTPUT ARGUMENTS:
    real, intent(out) :: BufferOut_IIBV(nMltIn, nLatIn, 2, nVarIn)

    character (len=*), parameter :: NameSub='UA_get_for_ie'

    call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')

  end subroutine UA_get_for_ie

  !============================================================================
  subroutine UA_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

    integer, intent(in) :: nDimIn
    integer, intent(in) :: nPoint
    real,    intent(in) :: Xyz_DI(nDimIn,nPoint)
    integer, intent(out):: iProc_I(nPoint)

    character(len=*), parameter:: NameSub = 'UA_find_points'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')
    
  end subroutine UA_find_points

  !============================================================================
  subroutine UA_get_for_gm(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
       Data_VI)

    logical,          intent(in) :: IsNew
    character(len=*), intent(in) :: NameVar
    integer,          intent(in) :: nVarIn
    integer,          intent(in) :: nDimIn
    integer,          intent(in) :: nPoint
    real,             intent(in) :: Xyz_DI(nDimIn,nPoint)
    real, intent(out):: Data_VI(nVarIn,nPoint)

    character(len=*), parameter :: NameSub='UA_get_for_gm'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')

  end subroutine UA_get_for_gm

  !============================================================================
  subroutine UA_get_grid_info(nDimOut, iGridOut, iDecompOut)

    integer, intent(out):: nDimOut
    integer, intent(out):: iGridOut
    integer, intent(out):: iDecompOut

    character(len=*), parameter :: NameSub = 'UA_get_grid_info'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')

  end subroutine UA_get_grid_info
  
end module UA_wrapper
!==============================================================================

! The following subroutines are empty versions of those in UA/GITM2/src/
! The call to these routines is commented out in CON_couple_ie_ua.f90
!
!subroutine UA_fill_electrodynamics(UAr2_fac, UAr2_ped, UAr2_hal, &
!     UAr2_lats, UAr2_mlts)
!
!  character(len=*), parameter :: NameSub='UA_fill_electrodynamics'
!
!  call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')
!
!end subroutine UA_fill_electrodynamics
!
!!===========================================================================
!
!subroutine UA_calc_electrodynamics(UAi_nMLTs, UAi_nLats)
!
!  character(len=*), parameter :: NameSub='UA_calc_electrodynamics'
!
!  call CON_stop(NameSub//': UA_ERROR: empty version cannot be used!')
!
!end subroutine UA_calc_electrodynamics


