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
subroutine SWMF_initialize(iComm, iTimeStart_I, TimeSim, IsLastSession, iError)

  !USES:
  use CON_main,      ONLY: initialize
  use CON_session,   ONLY: init_session
  use CON_variables, ONLY: iErrorSwmf
  use CON_time,      ONLY: TimeStart, tSimulation
  use ModTimeConvert,ONLY: time_int_to_real, time_real_to_int
  use ModKind,       ONLY: Real4_

  implicit none

  !INPUT ARGUMENTS:
  integer, intent(in) :: iComm           ! The MPI communicator for the SWMF
  integer, intent(in) :: iTimeStart_I(7) ! Start time (year ... millisec)
  real(Real4_), intent(in) :: TimeSim    ! Simulation time (0.0 unless restart)

  !OUTPUT ARGUMENTS:
  integer, intent(out):: iError         ! SWMF error code (0 on success)
  logical, intent(out):: IsLastSession  ! True if there is only one session

  !DESCRIPTION:
  ! Obtains the MPI communicator for the whole SWMF.
  ! Initializes the SWMF for the first session.
  ! Indicates if there are more than one sessions according to the PARAM.in.
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

  ! Set the simulation time CON_time::tSimulation
  tSimulation = TimeSim

  ! Initialize for first session
  call init_session(IsLastSession)
  iError = iErrorSwmf
  !EOC
end subroutine SWMF_initialize

!BOP ==========================================================================
!ROUTINE: SWMF_initialize_session - initialize SWMF for 2nd, 3rd etc. sessions
!INTERFACE:
subroutine SWMF_initialize_session(IsLastSession, iError)
  !USES:
  use CON_session,   ONLY: init_session
  use CON_variables, ONLY: iErrorSwmf
  implicit none
  !OUTPUT ARGUMENTS:
  logical, intent(out):: IsLastSession ! True if this is the last session
  integer, intent(out):: iError        ! Error code, 0 on success
  !DESCRIPTION:
  ! Initialize for 2nd, 3rd etc. sessions in a multi-session run.
  ! Indicates if this is the last session to run.
  !EOP
  !BOC ------------------------------------------------------------------------
  call init_session(IsLastSession)
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
  ! Run the SWMF for a certain time. To be implemented !!!
  ! Indicate if a stop has been requested.
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
