!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module IE_wrapper

  use ModUtilities, ONLY: CON_stop

  implicit none

contains
  ! Wrapper for the "empty" Ionosphere Electrodynamics (IE) component
  !==========================================================================
  subroutine IE_set_param(CompInfo, TypeAction)

    use CON_comp_info

    character (len=*), parameter :: NameSub='IE_set_param'

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
       call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')
    end select

  end subroutine IE_set_param
  !============================================================================
  subroutine IE_init_session(iSession, TimeSimulation)

    !INPUT PARAMETERS:
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='IE_init_session'

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_init_session
  !============================================================================
  subroutine IE_finalize(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='IE_finalize'

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_finalize
  !============================================================================
  subroutine IE_save_restart(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='IE_save_restart'

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_save_restart
  !============================================================================
  subroutine IE_run(TimeSimulation,TimeSimulationLimit)

    !INPUT/OUTPUT ARGUMENTS:
    real, intent(inout) :: TimeSimulation   ! current time of component

    !INPUT ARGUMENTS:
    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    character(len=*), parameter :: NameSub='IE_run'

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_run
  !============================================================================
  subroutine IE_get_for_gm(Buffer_IIV, iSize, jSize, nVar, NameVar_I, &
       tSimulation)

    integer,          intent(in) :: iSize, jSize, nVar
    real,             intent(out):: Buffer_IIV(iSize,jSize,nVar)
    character(len=*), intent(in) :: NameVar_I(nVar)
    real,             intent(in) :: tSimulation

    character (len=*),parameter :: NameSub = 'IE_get_for_gm'

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_get_for_gm
  !============================================================================
  subroutine IE_put_from_gm(Buffer_IIV, iSize, jSize, nVar)

    integer,          intent(in) :: iSize, jSize, nVar
    real,             intent(in) :: Buffer_IIV(iSize,jSize,nVar)

    character (len=*),parameter :: NameSub='IE_put_from_gm'

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_put_from_gm
  !============================================================================
  subroutine IE_get_for_pw(Buffer_IIV, iSize, jSize, nVar, Name_V, NameHem,&
       tSimulation)

    character (len=*),parameter :: NameSub='IE_get_for_pw'

    integer, intent(in)           :: iSize, jSize, nVar
    real, intent(out)             :: Buffer_IIV(iSize,jSize,nVar)
    character (len=*),intent(in)  :: NameHem
    character (len=*),intent(in)  :: Name_V(nVar)
    real,             intent(in)  :: tSimulation

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_get_for_pw
  !============================================================================
  subroutine IE_get_for_rb(Buffer_IIV, iSize, jSize, nVar, Name_V, NameHem,&
       tSimulation)

    character (len=*),parameter :: NameSub='IE_get_for_rb'

    integer, intent(in)           :: iSize, jSize, nVar
    real, intent(out)             :: Buffer_IIV(iSize,jSize,nVar)
    character (len=*),intent(in)  :: NameHem
    character (len=*),intent(in)  :: Name_V(nVar)
    real,             intent(in)  :: tSimulation

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_get_for_rb
  !============================================================================
  subroutine IE_get_for_ps(Buffer_II, iSize, jSize, tSimulation)

    integer, intent(in) :: iSize, jSize
    real, intent(out)   :: Buffer_II(iSize,jSize)
    real, intent(in)    :: tSimulation

    character (len=*),parameter :: NameSub='IE_get_for_ps'

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_get_for_ps
  !============================================================================
  subroutine IE_get_info_for_im(use_ua, nEngInput, nVarImIe)

    integer,          intent(out) :: nVarImIe
    integer,          intent(in) :: nEngInput
    logical,          intent(in) :: use_ua

    character (len=*), parameter :: NameSub='IE_get_info_for_im'

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_get_info_for_im
  !============================================================================
  subroutine IE_get_for_im(nPoint,iPointStart,Index,Weight,Buff_V,nVar)

    use CON_router,   ONLY: IndexPtrType, WeightPtrType

    integer,intent(in)            :: nPoint, iPointStart, nVar
    real,intent(out)              :: Buff_V(nVar)
    type(IndexPtrType),intent(in) :: Index
    type(WeightPtrType),intent(in):: Weight

    character (len=*),parameter :: NameSub='IE_get_for_im'

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_get_for_im
  !============================================================================
  subroutine IE_put_from_UA(Buffer_IIBV, nMLTs, nLats, nVarIn, NameVarUaIn_V)

    integer,          intent(in) :: nMlts, nLats, nVarIn
    character(len=3), intent(in) :: NameVarUaIn_V(nVarIn)
    real,             intent(in) :: Buffer_IIBV(nMlts, nLats, 2, nVarIn)

    character (len=*),parameter :: NameSub='IE_put_from_UA'

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_put_from_UA
  !============================================================================
  subroutine IE_get_info_for_ua(nVar, NameVar_V)

    integer,          intent(out)           :: nVar
    character(len=*), intent(out), optional :: NameVar_V(:)

    character(len=*), parameter :: NameSub='IE_get_info_for_ua'

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_get_info_for_ua

  !============================================================================

  subroutine IE_get_for_ua(Buffer_IIV,iSize,jSize,nVarIn,NameVar_V, &
       iBlock,tSimulation)

    integer,          intent(in)  :: iSize,jSize, nVarIn, iBlock
    real,             intent(out) :: Buffer_IIV(iSize,jSize,nVarIn)
    character (len=*),intent(in)  :: NameVar_V(nVarIn)
    real,             intent(in)  :: tSimulation

    character (len=*),parameter :: NameSub='IE_get_for_ua'

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_get_for_ua
  !============================================================================
  subroutine IE_setnMlts(iComponent, nMLTsIn, iError)

    integer, intent(in)  :: iComponent, nMLTsIn
    integer, intent(out) :: iError

    character (len=*), parameter :: NameSub='IE_setnMlts'

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_setnMlts
  !============================================================================
  subroutine IE_setnLats(iComponent, nLatsIn, iError)

    integer, intent(in)  :: iComponent, nLatsIn
    integer, intent(out) :: iError

    character (len=*), parameter :: NameSub='IE_setnLats'

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_setnLats
  !============================================================================
  subroutine IE_setgrid(iComponent, MLTsIn, LatsIn, iError)

    integer, intent(in) :: iComponent
    real, intent(in) :: MLTsIn,LatsIn
    integer, intent(out) :: iError

    character (len=*), parameter :: NameSub='IE_setgrid'

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_setgrid
  !============================================================================
  subroutine IE_put_from_im(nPoint,iPointStart,Index,Weight,DoAdd,Buff_V,nVar)

    use CON_router,   ONLY: IndexPtrType, WeightPtrType

    character(len=*), parameter   :: NameSub='IE_put_from_im'
    integer,intent(in)            :: nPoint, iPointStart, nVar
    real, intent(in)              :: Buff_V(nVar)
    type(IndexPtrType),intent(in) :: Index
    type(WeightPtrType),intent(in):: Weight
    logical,intent(in)            :: DoAdd

    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_put_from_im
  !============================================================================
  subroutine IE_put_from_im_complete

    write(*,*)"What is IE_put_from_im_complete really supposed to do?"

  end subroutine IE_put_from_im_complete
  !============================================================================
end module IE_wrapper
!============================================================================
subroutine SPS_put_into_ie(Buffer_II, iSize, jSize, NameVar, iBlock)

  use ModUtilities, ONLY: CON_stop
  implicit none

  integer, intent(in)           :: iSize,jSize
  real, intent(in)              :: Buffer_II(iSize,jSize)
  character (len=*),intent(in)  :: NameVar
  integer,intent(in)            :: iBlock

  character (len=*), parameter :: NameSub='SPS_put_into_ie'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine SPS_put_into_ie
!============================================================================
subroutine initialize_ie_ua_buffers(iOutputError)

  use ModUtilities, ONLY: CON_stop
  implicit none

  integer :: iOutputError

  character (len=*),parameter :: NameSub='initialize_ie_ua_buffers'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine initialize_ie_ua_buffers
!============================================================================

