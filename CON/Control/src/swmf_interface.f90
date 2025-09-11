!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! The subroutines in this file provide an interface to the SWMF library,
! when it does not run in stand alone mode.
! These are all external subroutines, so that the external application
! does not have to compile any SWMF modules
! (which avoids a lot of compilation problems).
! All subroutines return an error code, which is 0 on success.

! revision history:
! 09/01/05 G.Toth - initial version
! 02/06/06          modified SWMF_run to pass coupling time
!                   added SWMF_couple to transfer 2D grid data
subroutine SWMF_initialize(iComm, iTimeStart_I, TimeSim, TimeStop, &
     IsLastSession, iError)

  use CON_main, ONLY: initialize
  use CON_variables, ONLY: iErrorSwmf
  use CON_time, ONLY: TimeStart, tSimulation, tSimulationMax, &
       DoTimeAccurate, MaxIteration
  use ModTimeConvert, ONLY: time_int_to_real, time_real_to_int
  use ModKind, ONLY: Real8_

  implicit none

  integer, intent(in) :: iComm           ! The MPI communicator for the SWMF
  integer, intent(in) :: iTimeStart_I(7) ! Start time (year ... millisec)
  real(Real8_), intent(in) :: TimeSim    ! Simulation time (0.0 unless restart)
  real(Real8_), intent(in) :: TimeStop   ! Final simulation time

  integer, intent(out):: iError         ! SWMF error code (0 on success)
  logical, intent(out):: IsLastSession  ! True if there is only one session

  ! initialize the SWMF (for the first session)
  ! Obtains the MPI communicator for the whole SWMF.
  ! Initializes the SWMF for the first session with a start date and time,
  ! current and final simulation time. The output arguments indicate
  ! if there are more than one sessions according to the PARAM.in file.
  ! The current and final simulation times are passed as 8 byte reals,
  ! because this is available on all platforms.
  ! \newpage
  ! Initialize SWMF
  !----------------------------------------------------------------------------
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

end subroutine SWMF_initialize
!==============================================================================
subroutine SWMF_initialize_session(IsLastSession, iError)

  ! read parameters and initialize session

  use CON_io, ONLY: read_inputs
  use CON_session, ONLY: init_session
  use CON_variables, ONLY: iErrorSwmf
  implicit none
  logical, intent(out):: IsLastSession ! True if this is the last session
  integer, intent(out):: iError        ! Error code, 0 on success
  ! Read input parameters and initialize a session.
  ! The IsLastSession argument indicates if this is the last session to run
  ! according to the PARAM.in file.
  !----------------------------------------------------------------------------
  call read_inputs(IsLastSession)
  if(iErrorSwmf == 0) call init_session
  iError = iErrorSwmf

end subroutine SWMF_initialize_session
!==============================================================================
subroutine SWMF_run(NameComp, tCouple, tSimulationOut, DoStop, iError)

  use CON_session, ONLY: do_session
  use CON_variables, ONLY: iErrorSwmf
  use CON_time, ONLY: tSimulation
  use CON_comp_param, ONLY: MaxComp, lNameComp, i_comp_name
  use ModKind, ONLY: Real8_

  implicit none

  character(len=lNameComp), intent(in):: NameComp ! Component to couple with
  real(Real8_),             intent(in):: tCouple  ! Next coupling time
  real(Real8_), intent(out):: tSimulationOut ! Current SWMF simulation time
  logical,      intent(out):: DoStop         ! True if SWMF requested a stop
  integer,      intent(out):: iError         ! Error code, 0 on success

  ! Run the SWMF until the coupling time or a stop condition is reached.
  ! If NameComp is the name of one of the SWMF components then the coupling
  ! time tCouple affects only the processors that run component NameComp.
  ! If NameComp is set to ** then all the SWMF processors synchronize
  ! and return at the coupling time. The tSimulationOut contains the
  ! simulation time at return. The DoStop argument indicates if a
  ! final stop has been requested by the SWMF.
  ! local variables

  real    :: tCouple_C(MaxComp)
  !----------------------------------------------------------------------------
  tCouple_C = -1.0
  if(tCouple >= 0.0)then
     if(NameComp == '**')then
        ! couple with all the components
        tCouple_C = tCouple
     else
        ! couple with a single component
        tCouple_C(i_comp_name(NameComp)) = tCouple
     end if
  end if

  call do_session(DoStop, tCouple_C)

  tSimulationOut = tSimulation
  iError = iErrorSwmf

end subroutine SWMF_run
!==============================================================================
subroutine SWMF_finalize(iError)

  use CON_main, ONLY: finalize
  use CON_variables, ONLY: iErrorSwmf

  implicit none

  integer, intent(out):: iError ! Error code, 0 on success

  ! Finalize the SWMF after the last session is done.
  !----------------------------------------------------------------------------
  call finalize
  iError = iErrorSwmf

end subroutine SWMF_finalize
!==============================================================================
subroutine SWMF_couple(NameFrom, NameTo, NameCoord, &
     nVar, nX, nY, xMin, xMax, yMin, yMax, Data_VII, iError)

  use ModKind, ONLY: Real8_
  use CON_comp_param, ONLY: lNameComp, i_comp_name
  use CON_world, ONLY: i_comm
  use ModMpiOrig

  use GM_wrapper, ONLY: GM_put_from_ih_buffer

  implicit none
  character(len=*), intent(in) :: NameFrom   ! Provider component
  character(len=*), intent(in) :: NameTo     ! Receiver component
  character(len=3), intent(in) :: NameCoord  ! Coordinate system description
  integer,          intent(in) :: nVar       ! Number of variables
  integer,          intent(in) :: nX, nY     ! Grid size
  real(Real8_),     intent(in) :: xMin, xMax ! X range
  real(Real8_),     intent(in) :: yMin, yMax ! Y range

  real(Real8_),  intent(inout) :: Data_VII(nVar, nX, nY) ! grid data pointer

  integer, intent(out):: iError ! Error code, 0 on success

  ! SWMF coupling with an external code
  ! The coupling interface of the SWMF when coupled to an external code.
  ! This subroutine can be used to couple with various components.
  ! The couplings are identified by the NameFrom and NameTo strings.
  ! One of the strings must be the name of an SWMF component, the
  ! other string identifies the external component and it should not
  ! coincide with any of the SWMF component names.
  ! The data is transferred on a uniform 2D grid. The grid dimension,
  ! coordinate system and coordinate ranges are determined by the external
  ! code. The data array is passed to/from the appropriate SWMF component.
  ! The coordinate ranges and the data must be in SI units.
  ! The coordinate system description is a 3 character string (e.g. GSM)
  ! which is used to do the coordinate transformation between the SWMF
  ! component and the buffer grid.
  integer :: iComm
  character(len=*), parameter:: NameSub = 'SWMF_couple'
  !----------------------------------------------------------------------------
  iError = 1

  ! Do the appropriate coupling depending on the NameFrom-NameTo pair
  select case(NameFrom)
  case('ESMF_IPE')
     select case(NameTo)
     case('IE')
        if(nVar /= 2)then
           write(*,*)NameSub//' ERROR: '// &
                'coupling to IE requires 2 variables, not nVar=',nVar
           RETURN
        end if

        ! Put data into IE. Convert to the default real used by SWMF.
        ! call IE_put_from_ih_buffer(NameCoord, nX, nY, &
        !     real(xMin), real(xMax), real(yMin), real(yMax), real(Data_VII))
     case default
        write(*,*)NameSub//' ERROR: '// &
             'coupling to '//NameTo//' is not implemented'
        RETURN
     end select
  case default
     write(*,*)NameSub//' ERROR: '// &
          'coupling from '//NameFrom//' is not implemented'
     RETURN
  end select

  iError = 0

end subroutine SWMF_couple
!==============================================================================

