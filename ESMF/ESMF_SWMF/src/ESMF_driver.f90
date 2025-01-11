program ESMF_driver

  !  Application Driver for the coupled ESMF-SWMF system.
  !  Creates the top ESMF_SWMF Gridded Component and calls the
  !  Initialize, Run, and Finalize routines for it.
  !
  !  The top Gridded Component creates and manages the ESMF and SWMF
  !  subcomponents internally. The SWMF is treated as a single component
  !  which is coupled to (some of) the ESMF component(s) periodically.

  ! ESMF module, defines all ESMF data types and procedures
  use ESMF

  ! Top ESMF-SWMF Gridded Component registration routines
  use ESMF_grid_comp, ONLY: ESMF_set_services

  ! Various variables
  use ESMFSWMF_variables, ONLY: NameParamFile, iProc, nProc, &
       Year_, Month_, Day_, Hour_, Minute_, Second_, MilliSec_, &
       iStartTime_I, iFinishTime_I, TimeSimulation, &
       read_esmf_swmf_input, write_log, write_error

  implicit none

  ! Components
  type(ESMF_GridComp) :: EsmfSwmfComp
  ! States, Virtual Machines, Layouts and processor index
  type(ESMF_VM)       :: DefaultVM

  ! A clock, starting and stop times and timestep
  type(ESMF_Clock)        :: Clock
  type(ESMF_Time)         :: StartTime
  type(ESMF_Time)         :: StopTime
  type(ESMF_TimeInterval) :: TimeStep

  ! Variables for the clock

  ! Simulation time (non-zero for restart)
  type(ESMF_Time)         :: CurrentTime
  type(ESMF_TimeInterval) :: SimTime
  integer(ESMF_KIND_I4)   :: iSecond, iMillisec

  ! Store default values before reading

  ! Return code for error checks
  integer :: iError

  ! Initialize the ESMF Framework
  !----------------------------------------------------------------------------
  call ESMF_Initialize(configFileName=trim(NameParamFile), &
       defaultCalkind=ESMF_CALKIND_GREGORIAN, rc=iError)
  if (iError /= ESMF_SUCCESS) stop 'ESMF_Initialize FAILED'

  call write_log("ESMF-SWMF Driver start")

  ! Get the default VM which contains all PEs this job was started on.
  call ESMF_VMGetGlobal(defaultVM, rc=iError)
  if(iError /= ESMF_SUCCESS) call my_error('ESMF_VMGetGlobal failed')

  call ESMF_VMGet(defaultVM, petcount=nProc, localpet=iProc, rc=iError)
  if(iError /= ESMF_SUCCESS) call my_error('ESMF_VMGet failed')

  ! Create section

  ! Create the top Gridded component, passing in the default layout.
  EsmfSwmfComp = ESMF_GridCompCreate(name="ESMF-SWMF Component", rc=iError)
  if(iError /= ESMF_SUCCESS) call my_error('ESMF_GridCompCreate')
  call write_log("Component Create finished")

  ! Register section
  call ESMF_GridCompSetServices(EsmfSwmfComp, &
       userRoutine=ESMF_set_services, rc=iError)
  if(iError /= ESMF_SUCCESS) call my_error('ESMF_GridCompSetServices')

  ! Read input paramterers
  call read_esmf_swmf_input(iError)
  if(iError /= ESMF_SUCCESS) call my_error('call read_esmf_swmf_input failed')

  ! Create and initialize a clock
  ! Based on values from the Config file, create a Clock.

  call ESMF_TimeSet(StartTime,    &
       yy=iStartTime_I(Year_),    &
       mm=iStartTime_I(Month_),   &
       dd=iStartTime_I(Day_),     &
       h =iStartTime_I(Hour_),    &
       m =iStartTime_I(Minute_),  &
       s =iStartTime_I(Second_),  &
       ms=iStartTime_I(Millisec_),&
       rc=iError)

  if(iError /= ESMF_SUCCESS) call my_error('ESMF_TimeSet StartTime')

  call ESMF_TimeSet(StopTime,     &
       yy=iFinishTime_I(Year_),   &
       mm=iFinishTime_I(Month_),  &
       dd=iFinishTime_I(Day_),    &
       h =iFinishTime_I(Hour_),   &
       m =iFinishTime_I(Minute_), &
       s =iFinishTime_I(Second_), &
       ms=iFinishTime_I(Millisec_))

  if(iError /= ESMF_SUCCESS) call my_error('ESMF_TimeSet StopTime')

  TimeStep= StopTime-StartTime
  Clock = ESMF_ClockCreate(TimeStep, StartTime, stopTime=StopTime, &
       runDuration=StopTime-StartTime, name="application Clock", rc=iError)
  if(iError /= ESMF_SUCCESS) call my_error('ESMF_ClockCreate')

  if(TimeSimulation /= 0.0)then
     iSecond   = int(TimeSimulation)
     iMillisec = nint(1000*(TimeSimulation-iSecond))
     call ESMF_TimeIntervalSet(SimTime, s=iSecond, ms=iMillisec, rc=iError)
     if(iError /= ESMF_SUCCESS) call my_error('SMF_TimeIntervalSet SimTime')

     CurrentTime = StartTime + SimTime
     call ESMF_ClockSet(Clock, currtime = CurrentTime, rc=iError)
     if(iError /= ESMF_SUCCESS) call my_error('ESMF_ClockSet Clock')
  end if

  !  Init, Run, and Finalize section
  call ESMF_GridCompInitialize(EsmfSwmfComp, clock=Clock, rc=iError)
  if (iError /= ESMF_SUCCESS) call my_error('ESMF_GridCompInitialize')

  call ESMF_GridCompRun(EsmfSwmfComp, clock=Clock, rc=iError)
  if (iError /= ESMF_SUCCESS) call my_error('ESMF_GridCompRun')

  call ESMF_GridCompFinalize(EsmfSwmfComp, clock=Clock, rc=iError)
  if (iError /= ESMF_SUCCESS) call my_error('ESMF_GridCompFinalize')

  ! Clean up
  call ESMF_ClockDestroy(clock, rc=iError)
  if (iError /= ESMF_SUCCESS) call my_error('SMF_ClockDestroy')

  call ESMF_GridCompDestroy(EsmfSwmfComp, rc=iError)
  if (iError /= ESMF_SUCCESS) call my_error('ESMF_GridCompDestroy')

  call ESMF_Finalize

contains
  !============================================================================
  subroutine my_error(String)

    character(len=*), intent(in) :: String

    !--------------------------------------------------------------------------
    call write_error('ESMF_driver '//String)

  end subroutine my_error
  !============================================================================
end program ESMF_driver
!==============================================================================
