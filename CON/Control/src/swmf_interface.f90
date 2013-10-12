! !  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
! !  For more information, see http://csem.engin.umich.edu/tools/swmf
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
! 02/06/06          modified SWMF_run to pass coupling time
!                   added SWMF_couple to transfer 2D grid data
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
subroutine SWMF_run(NameComp, tCouple, tSimulationOut, DoStop, iError)
  !USES:
  use CON_session,    ONLY: do_session
  use CON_variables,  ONLY: iErrorSwmf
  use CON_time,       ONLY: tSimulation
  use CON_comp_param, ONLY: MaxComp, lNameComp, i_comp_name
  use ModKind,        ONLY: Real8_
  implicit none
  !INPUT ARGUMENTS:
  character(len=lNameComp), intent(in):: NameComp ! Component to couple with
  real(Real8_),             intent(in):: tCouple  ! Next coupling time
  !OUTPUT ARGUMENTS:
  real(Real8_), intent(out):: tSimulationOut ! Current SWMF simulation time
  logical,      intent(out):: DoStop         ! True if SWMF requested a stop
  integer,      intent(out):: iError         ! Error code, 0 on success
  !DESCRIPTION:
  ! Run the SWMF until the coupling time or a stop condition is reached.
  ! If NameComp is the name of one of the SWMF components then the coupling 
  ! time tCouple affects only the processors that run component NameComp.
  ! If NameComp is set to ** then all the SWMF processors synchronize
  ! and return at the coupling time. The tSimulationOut contains the 
  ! simulation time at return. The DoStop argument indicates if a 
  ! final stop has been requested by the SWMF.
  !EOP
  !LOCAL VARIABLES:
  real    :: tCouple_C(MaxComp)
  !BOC ------------------------------------------------------------------------
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

!BOP ==========================================================================
!ROUTINE: SWMF_couple - SWMF coupling with an external code
!INTERFACE:
subroutine SWMF_couple(NameFrom, NameTo, NameCoord, &
     nVar, nX, nY, xMin, xMax, yMin, yMax, Data_VII, iError)

  !USES:
  use ModKind, ONLY: Real8_
  use CON_comp_param, ONLY: lNameComp, i_comp_name
  use CON_world,      ONLY: i_comm
  use ModMpi

  implicit none
  !INPUT ARGUMENTS:
  character(len=*), intent(in) :: NameFrom   ! Provider component
  character(len=*), intent(in) :: NameTo     ! Receiver component
  character(len=3), intent(in) :: NameCoord  ! Coordinate system description
  integer,          intent(in) :: nVar       ! Number of variables
  integer,          intent(in) :: nX, nY     ! Grid size
  real(Real8_),     intent(in) :: xMin, xMax ! X range
  real(Real8_),     intent(in) :: yMin, yMax ! Y range

  !INPUT/OUTPUT ARGUMENTS:
  real(Real8_),  intent(inout) :: Data_VII(nVar, nX, nY) ! grid data pointer

  !OUTPUT ARGUMENTS:
  integer, intent(out):: iError ! Error code, 0 on success
  !DESCRIPTION:
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
  !EOP
  integer :: iComm
  character(len=*), parameter :: NameSub='SWMF_receive'
  !BOC ------------------------------------------------------------------------
  iError = 1

  ! Do the appropriate coupling depending on the NameFrom-NameTo pair
  select case(NameFrom)
  case('ESMF_IH')
     select case(NameTo)
     case('GM')
        if(nVar /= 8)then
           write(*,*)NameSub//' ERROR: '// &
                'coupling to GM requires 8 variables, not nVar=',nVar
           return
        end if

        ! Broadcast information from the root processor 
        ! to the rest of the component
        iComm = i_comm(NameTo)
        call MPI_bcast(Data_VII, nVar*nX*nY, MPI_DOUBLE_PRECISION, 0, &
             iComm, iError)
        if(iError /= 0)then
           write(*,*)NameSub//' ERROR: '// &
                'MPI_bcast failed for component '//NameTo
           return
        end if

        ! Put data into GM. Convert to the default real used by SWMF.
        call GM_put_from_ih_buffer(NameCoord, nX, nY, &
             real(xMin), real(xMax), real(yMin), real(yMax), real(Data_VII))
     case default
        write(*,*)NameSub//' ERROR: '// &
             'coupling to '//NameTo//' is not implemented'
        return
     end select
  case default
     write(*,*)NameSub//' ERROR: '// &
          'coupling from '//NameFrom//' is not implemented'
     return
  end select

  iError = 0
  !EOC

end subroutine SWMF_couple

