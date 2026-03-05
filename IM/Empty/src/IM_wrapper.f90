!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module IM_wrapper

  ! Wrapper for empty Internal Magnetosphere (IM) component

  use ModUtilities, ONLY: CON_stop

  implicit none

  private ! except

  public:: IM_set_param
  public:: IM_init_session
  public:: IM_run
  public:: IM_save_restart
  public:: IM_finalize

  ! Coupling with IE
  public:: IM_get_info_for_ie
  public:: IM_get_for_ie
  public:: IM_put_from_ie_mpi
  public:: IM_put_from_ie
  public:: IM_put_from_ie_complete

  ! Coupling with GM
  public:: IM_get_for_gm
  public:: IM_put_from_gm
  public:: IM_put_from_gm_line
  public:: IM_put_from_gm_crcm
  public:: IM_put_sat_from_gm

contains
  !==========================================================================
  subroutine IM_set_param(CompInfo, TypeAction)

    use CON_comp_info

    character (len=*), parameter :: NameSub='IM_set_param'

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
       call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')
    end select

  end subroutine IM_set_param
  !============================================================================
  subroutine IM_init_session(iSession, TimeSimulation)

    !INPUT PARAMETERS:
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='IM_init_session'

    call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

  end subroutine IM_init_session
  !============================================================================
  subroutine IM_finalize(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='IM_finalize'

    call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

  end subroutine IM_finalize
  !============================================================================
  subroutine IM_save_restart(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='IM_save_restart'

    call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

  end subroutine IM_save_restart
  !============================================================================
  subroutine IM_run(TimeSimulation,TimeSimulationLimit)

    !INPUT/OUTPUT ARGUMENTS:
    real, intent(inout):: TimeSimulation   ! current time of component

    !INPUT ARGUMENTS:
    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    character(len=*), parameter :: NameSub='IM_run'

    call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

  end subroutine IM_run
  !============================================================================
  subroutine IM_put_from_ie_mpi(nTheta, nPhi, Potential_II)

    integer, intent(in):: nTheta, nPhi
    real,    intent(in):: Potential_II(nTheta, nPhi)

    character(len=*), parameter   :: NameSub='IM_put_from_ie_mpi'

    call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

  end subroutine IM_put_from_ie_mpi
  !============================================================================
  subroutine IM_put_from_ie(nPoint,iPointStart,Index,Weight,DoAdd,Buff_V,nVar)

    use CON_router,   ONLY: IndexPtrType, WeightPtrType

    integer,intent(in)            :: nPoint, iPointStart, nVar
    real, intent(in)              :: Buff_V(nVar)
    type(IndexPtrType),intent(in) :: Index
    type(WeightPtrType),intent(in):: Weight
    logical,intent(in)            :: DoAdd

    character(len=*), parameter   :: NameSub = 'IM_put_from_ie'

    call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

  end subroutine IM_put_from_ie
  !============================================================================
  subroutine IM_put_from_ie_complete

    character(len=*), parameter   :: NameSub='IM_put_from_ie_complete'

    call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

  end subroutine IM_put_from_ie_complete
  !============================================================================
  subroutine IM_put_from_gm_line(nRadiusIn, nLonIn, Map_DSII, &
       nVarLineIn, nPointLineIn, BufferLine_VI, NameVar)

    integer, intent(in) :: nRadiusIn, nLonIn
    real,    intent(in) :: Map_DSII(3,2,nRadiusIn,nLonIn)
    integer, intent(in) :: nVarLineIn, nPointLineIn
    real,    intent(in) :: BufferLine_VI(nVarLineIn,nPointLineIn)
    character(len=*), intent(in) :: NameVar

    character (len=*),parameter :: NameSub='IM_put_from_gm_line'

    call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

  end subroutine IM_put_from_gm_line
  !============================================================================
  subroutine IM_put_from_gm(Buffer_IIV,BufferKp,iSizeIn,jSizeIn,nVarIn,NameVar)

    integer, intent(in) :: iSizeIn,jSizeIn,nVarIn
    real,    intent(in) :: BufferKp
    real, dimension(iSizeIn,jSizeIn,nVarIn), intent(in) :: Buffer_IIV
    character (len=*),intent(in)       :: NameVar

    character (len=*),parameter :: NameSub='IM_put_from_gm'

    call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

  end subroutine IM_put_from_gm
  !============================================================================
  subroutine IM_get_for_gm(Buffer_IIV, iSizeIn, jSizeIn, nVar, NameVar)

    character (len=*),parameter :: NameSub='IM_get_for_gm'

    integer, intent(in) :: iSizeIn, jSizeIn, nVar
    real,    intent(out):: Buffer_IIV(iSizeIn,jSizeIn,nVar)
    character (len=*),intent(in):: NameVar

    call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

  end subroutine IM_get_for_gm
  !===========================================================================
  subroutine IM_put_from_gm_crcm( &
       Integral_IIV, BufferKp, BufferAe, iSizeIn, jSizeIn,&
       nIntegralIn, BufferLine_VI, nVarLine, nPointLine, NameVar, &
       SolarWind_V, tSimulation)

    integer, intent(in) :: iSizeIn, jSizeIn, nIntegralIn
    real,    intent(in) :: Integral_IIV(iSizeIn,jSizeIn,nIntegralIn)
    real,    intent(in) :: BufferKp, BufferAe
    integer, intent(in) :: nVarLine, nPointLine
    real,    intent(in) :: BufferLine_VI(nVarLine,nPointLine)
    real,    intent(in) :: SolarWind_V(8)
    real,    intent(in) :: tSimulation
    character (len=*), intent(in) :: NameVar

    character (len=*), parameter :: NameSub='IM_put_from_gm_crcm'

    call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

  end subroutine IM_put_from_gm_crcm
  !============================================================================
  subroutine IM_put_sat_from_gm(nSats, Buffer_I, Buffer_III)

    character (len=*), parameter :: NameSub='IM_put_sat_from_gm'

    integer, intent(in)            :: nSats
    real, intent(in)               :: Buffer_III(3,2,nSats)
    character(len=100), intent(in) :: Buffer_I(nSats)

    call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

  end subroutine IM_put_sat_from_gm
  !============================================================================
  subroutine IM_get_info_for_ie(nEngIM)

    character(len=*), parameter :: NameSub='IM_get_info_for_ie'

    integer, intent(out) :: nEngIM

    call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

  end subroutine IM_get_info_for_ie
  !============================================================================
  subroutine IM_get_for_ie(nPoint,iPointStart,Index,Weight,Buff_V,nVar)

    use CON_router,   ONLY: IndexPtrType, WeightPtrType

    character(len=*), parameter :: NameSub='IM_get_for_ie'

    integer,intent(in)            :: nPoint, iPointStart, nVar
    real,intent(out)              :: Buff_V(nVar)
    type(IndexPtrType),intent(in) :: Index
    type(WeightPtrType),intent(in):: Weight

    call CON_stop(NameSub//': IM_ERROR: empty version cannot be used!')

  end subroutine IM_get_for_ie
  !============================================================================
end module IM_wrapper
!==============================================================================
