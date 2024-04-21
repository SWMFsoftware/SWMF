!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModTime
  use ModTimeConvert, ONLY: time_int_to_julian, time_real_to_julian
  use ModKind,    ONLY: Real8_
  use ModUtilities, ONLY: CON_stop

  implicit none

  SAVE
  ! Iteration and time in SP
  integer      :: iIter   = 0
  real         :: SPTime  = 0.0
  ! Time of the last data output
  real         :: DataInputTime = 0.0
  logical      :: IsSteadyState = .false.
  ! SPTime is the time in the SP model, DataInputTime is the time
  ! of last data input. These two varables are set/modified in the
  ! following ways
  ! INITIALLY

  ! 1. SPTime and DataInputTime are set to zero in the current file
  !    which does not matter since in ANY case these values will be
  !    overwritten.
  ! 2. They both may be set with the command #TIMESIMULATION in the
  !    PARAM.in command. Which also does not matter since in ANY
  !    case these values will be overwritten.
  ! 3. If the SP model works within the SWMF in the first
  !    init_session command in the SP_wrapper both SPTime and
  !    DataInputTime are assigned to the SWMF simulation time
  ! 4. Both within the SWMF and in a stand-alone mode if DoReadMhData
  !    is set to .true. during the reading THE FIRST file the
  !    DataInputTime is read from the first MhData file AND SPTime is
  !    assigned to be equal to this read value too.
  ! IN THE RUN TIME
  ! 1. During each coupling new value is assigned to DataInputTime
  !    equal to the SWMF simulation time at the moment of coupling
  !    The SPTime keeps unchanged and is LESS THAN DataInputTile.
  ! 2. While reading MhData file a new value is assigned to
  !    DataInputTime which is written in the said file. The SPTime
  !    keeps unchanged and is LESS THAN DataInputTile.
  ! 3. In the course of SP_run the difference in values of density
  !    related to "old" state at SPTime and new one at DataInputTime
  !    is used to calculate the Lagrangian time derivatives.
  ! 4. With these Lagrangian derivatives the coefficients in the
  !    kinetic equation are used to solve the kinetic equation
  !    and to advance the momentum distribution function in time
  !    from SPTime to DataInputTime using advance routine. At the
  !    end of advance routine the SPTime is synchronized with
  !    DataInputTime.

  ! StartTime converted to real
  real(Real8_) :: StartTime
  real(Real8_) :: StartTimeJulian
  ! In the stand-alone mode StartTime is read from PARAM.in file
  ! in the format iYear, iMonth, iDay, iHour, iMinute, iSecond
  ! Function time_int_to_real is used then to assign the value
  ! of Real8_StartTime. Alternatively, by applying the inverse
  ! function, time_real_to_int to the total of Real8_StartTime +
  ! SPTime, we can obtain the date and time, which may be used
  ! for example, in the file names to label the exact time for
  ! the data sets.

  ! This is the same default value as in the SWMF
  integer:: iStartTime_I(7) = [2000,3,21,10,45,0,0]

contains
  !============================================================================
  subroutine read_param(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand ! From PARAM.in
    logical :: DoTimeAccurate = .true.
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
       case('#NSTEP')
          call read_var('nStep',iIter)
       case('#TIMESIMULATION')
          call read_var('SPTime',SPTime)
          ! The last time input occurred the same time
          DataInputTime = SPTime
       case("#STARTTIME", "#SETREALTIME")
          call read_var('iYear'  ,iStartTime_I(1))
          call read_var('iMonth' ,iStartTime_I(2))
          call read_var('iDay'   ,iStartTime_I(3))
          call read_var('iHour'  ,iStartTime_I(4))
          call read_var('iMinute',iStartTime_I(5))
          call read_var('iSecond',iStartTime_I(6))
          iStartTime_I(7) = 0
       case("#TIMEACCURATE")
          call read_var('DoTimeAccurate', DoTimeAccurate)
          IsSteadyState = .not. DoTimeAccurate
    case default
       call CON_stop(NameSub//' Unknown command '//NameCommand)
    end select

  end subroutine read_param
  !============================================================================
  subroutine init
    use ModTimeConvert, ONLY: time_int_to_real
    !--------------------------------------------------------------------------
    call time_int_to_real(iStartTime_I, StartTime)
    ! also save the start time in Julian days;
    call time_int_to_julian(iStartTime_I, StartTimeJulian)

  end subroutine init
  !============================================================================
end module PT_ModTime
!==============================================================================
