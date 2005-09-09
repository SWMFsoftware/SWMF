!^CMP COPYRIGHT UM
!BOP
!MODULE: SWMF interface - a set of external subroutines to drive the SWMF
!DESCRIPTION:
! The subroutines in this file provide an interface to the SWMF library, 
! when it does not run in stand alone mode.
! These are all external subroutines, so that the external application 
! does not have to compile any SWMF modules 
! (which avoids a lot of compilation problems).
! All subroutines return an error code, which is 0 on success.

!REVISION HISTORY:
! 09/01/05 G.Toth - initial version
!EOP

!BOP ==========================================================================
!ROUTINE: SWMF_initialize - initialize the SWMF (for the first session)
!INTERFACE:
subroutine SWMF_initialize(iComm, iTimeStart_I, TimeSim, TimeStop, &
     IsLastSession, iError)

  !USES:
  use CON_main,      ONLY: initialize
  use CON_variables, ONLY: iErrorSwmf
  use CON_time,      ONLY: TimeStart, tSimulation, tSimulationMax, &
       DoTimeAccurate, MaxIteration
  use ModTimeConvert,ONLY: time_int_to_real, time_real_to_int
  use ModKind,       ONLY: Real8_

  implicit none

  !INPUT ARGUMENTS:
  integer, intent(in) :: iComm           ! The MPI communicator for the SWMF
  integer, intent(in) :: iTimeStart_I(7) ! Start time (year ... millisec)
  real(Real8_), intent(in) :: TimeSim    ! Simulation time (0.0 unless restart)
  real(Real8_), intent(in) :: TimeStop   ! Final simulation time

  !OUTPUT ARGUMENTS:
  integer, intent(out):: iError         ! SWMF error code (0 on success)
  logical, intent(out):: IsLastSession  ! True if there is only one session

  !DESCRIPTION:
  ! Obtains the MPI communicator for the whole SWMF.
  ! Initializes the SWMF for the first session with a start date and time,
  ! current and final simulation time. The output arguments indicate
  ! if there are more than one sessions according to the PARAM.in file.
  ! The current and final simulation times are passed as 8 byte reals,
  ! because this is available on all platforms.
  ! \newpage
  !EOP
  !BOC ------------------------------------------------------------------------
  ! Initialize SWMF
  call initialize(iComm)
  iError = iErrorSwmf
  if(iError /= 0) RETURN

  ! Set the real*8 field of CON_time::StartTime
  call time_int_to_real(iTimeStart_I, TimeStart % Time)
  ! Set the integer and string parts of CON_time::StartTime
  call time_real_to_int(TimeStart)

  ! Set the current simulation time
  tSimulation = TimeSim

  ! Set the final simulation time
  tSimulationMax = TimeStop

  ! Set time accurate mode
  DoTimeAccurate = .true.
  MaxIteration   = -1

  ! Initialize first session
  call SWMF_initialize_session(IsLastSession, iError)
  !EOC
end subroutine SWMF_initialize

!BOP ==========================================================================
!ROUTINE: SWMF_initialize_session - read parameters and initialize session
!INTERFACE:
subroutine SWMF_initialize_session(IsLastSession, iError)
  !USES:
  use CON_io,        ONLY: read_inputs
  use CON_session,   ONLY: init_session
  use CON_variables, ONLY: iErrorSwmf
  implicit none
  !OUTPUT ARGUMENTS:
  logical, intent(out):: IsLastSession ! True if this is the last session
  integer, intent(out):: iError        ! Error code, 0 on success
  !DESCRIPTION:
  ! Read input parameters and initialize a session.
  ! The IsLastSession argument indicates if this is the last session to run
  ! according to the PARAM.in file.
  !EOP
  !BOC ------------------------------------------------------------------------
  call read_inputs(IsLastSession)
  if(iErrorSwmf == 0) call init_session
  iError = iErrorSwmf
  !EOC
end subroutine SWMF_initialize_session

!BOP ==========================================================================
!ROUTINE: SWMF_run - run the SWMF
!INTERFACE:
subroutine SWMF_run(DoStop, iError)
  !USES:
  use CON_session,   ONLY: do_session
  use CON_variables, ONLY: iErrorSwmf
  implicit none
  !OUTPUT ARGUMENTS:
  logical, intent(out):: DoStop ! True if the SWMF requested a stop
  integer, intent(out):: iError ! Error code, 0 on success
  !DESCRIPTION:
  ! Run the SWMF until a stop condition is reached.
  ! The DoStop argument indicates if a final stop has been requested 
  ! by the SWMF.
  !EOP
  !BOC ------------------------------------------------------------------------
  call do_session(DoStop)
  iError = iErrorSwmf
  !EOC
end subroutine SWMF_run

!BOP ==========================================================================
!ROUTINE: SWMF_finalize - finalize the SWMF
!INTERFACE:
subroutine SWMF_finalize(iError)
  !USES:
  use CON_main,      ONLY: finalize
  use CON_variables, ONLY: iErrorSwmf
  implicit none
  !OUTPUT ARGUMENTS:
  integer, intent(out):: iError ! Error code, 0 on success
  !DESCRIPTION:
  ! Finalize the SWMF after the last session is done.
  !EOP
  !BOC ------------------------------------------------------------------------
  call finalize
  iError = iErrorSwmf
  !EOC
end subroutine SWMF_finalize
