program ESMF_SWMF_Driver

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
  use ESMF_SWMF_GridCompMod, only : SetServices => ESMF_SWMF_SetServices

  ! Various variables
  use ESMF_SWMF_Mod

  implicit none

  ! Local variables
  character (len=*), parameter :: NameParamFile = "ESMF_SWMF.input"

  ! Components
  type(ESMF_GridComp) :: EsmfSwmfComp
  ! States, Virtual Machines, Layouts and processor index
  type(ESMF_VM)       :: DefaultVM
  integer             :: iProc, nProc

  ! The grid used to pass MHD state at the SWFM/GM inflow boundary
  type(ESMF_Grid) :: grid

  ! A clock, starting and stop times and timestep
  type(ESMF_Clock)        :: Clock
  type(ESMF_Time)         :: StartTime
  type(ESMF_Time)         :: StopTime
  type(ESMF_TimeInterval) :: TimeStep

  ! Variables for the clock
  integer :: iTime                              ! Index for date-time arrays

  ! Simulation time (non-zero for restart)
  type(ESMF_Time)         :: CurrentTime
  type(ESMF_TimeInterval) :: SimTime
  integer(ESMF_KIND_I4)   :: iSecond, iMillisec

  ! Store default values before reading
  integer :: iDefaultTmp                        ! Temporary default integer
  real(ESMF_KIND_R8)    :: DefaultTmp           ! Temporary default real

  ! Return code for error checks
  integer :: rc
  !----------------------------------------------------------------------------

  ! Initialize the ESMF Framework
  call ESMF_Initialize(defaultCalkind=ESMF_CALKIND_GREGORIAN, rc=rc)
  if (rc /= ESMF_SUCCESS) stop 'ESMF_Initialize FAILED'

  call ESMF_LogWrite("ESMF-SWMF Driver start", ESMF_LOGMSG_INFO)

  ! Get the default VM which contains all PEs this job was started on.
  call ESMF_VMGetGlobal(defaultVM, rc=rc)
  if(rc /= ESMF_SUCCESS) call my_error('ESMF_VMGetGlobal failed')

  ! Get the processor number
  call ESMF_VMGet(defaultVM, petcount=nProc, localpet=iProc, rc=rc)
  if(rc /= ESMF_SUCCESS) call my_error('ESMF_VMGet failed')

  ! Read input paramterers
  call read_esmf_swmf_input(nProc, iProc, rc)
  if(rc /= ESMF_SUCCESS) call my_error('call read_esmf_swmf_input failed')

  ! Create section

  ! Create the top Gridded component, passing in the default layout.
  EsmfSwmfComp = ESMF_GridCompCreate(name="ESMF-SWMF Component", rc=rc)

  call ESMF_LogWrite("Component Create finished", ESMF_LOGMSG_INFO)

  ! Register section
  !!! DOES NOT COMPILE???
!!!  call ESMF_GridCompSetServices(EsmfSwmfComp, userRoutine=SetServices, rc=rc)
!  if (ESMF_LogFoundError(rcToCheck=rc, msg='Registration failed', &
!       line=__LINE__, file=__FILE__)) &
!       call ESMF_Finalize

  !  Create and initialize a clock, and a grid.

  ! Based on values from the Config file, create a Clock.  

  call ESMF_TimeIntervalSet(TimeStep, s=iCoupleFreq, rc=rc)

  if(rc /= ESMF_SUCCESS) then
     if(iProc==0)write(*,*)'ESMF_SWMF ERROR: iCoupleFreq=', iCoupleFreq
     call my_error('ESMF_TimeIntervalSet(s=iCoupleFreq) failed')
  end if

  call ESMF_TimeSet(startTime,    &
       yy=iStartTime_I(Year_),    &
       mm=iStartTime_I(Month_),   &
       dd=iStartTime_I(Day_),     &
       h =iStartTime_I(Hour_),    &
       m =iStartTime_I(Minute_),  &
       s =iStartTime_I(Second_),  &
       ms=iStartTime_I(Millisec_),&
       rc=rc)

  if(rc /= ESMF_SUCCESS) then
     if(iProc == 0)write(*,*) 'ESMF_SWMF ERROR: ',&
          'Setting start time failed:', iStartTime_I
     call ESMF_Finalize
  end if

  call ESMF_TimeSet(stopTime,     &
       yy=iFinishTime_I(Year_),   &
       mm=iFinishTime_I(Month_),  &
       dd=iFinishTime_I(Day_),    &
       h =iFinishTime_I(Hour_),   &
       m =iFinishTime_I(Minute_), &
       s =iFinishTime_I(Second_), &
       ms=iFinishTime_I(Millisec_))

  if(rc /= ESMF_SUCCESS) then
     if(iProc == 0)write(*,*) 'ESMF_SWMF ERROR: ',&
          'Setting finish time failed:',iFinishTime_I
     call ESMF_Finalize
  end if

  Clock = ESMF_ClockCreate(TimeStep, StartTime, stopTime=StopTime, &
       name="pplication Clock", rc=rc)

  if(rc /= ESMF_SUCCESS)then
     if(iProc == 0)write(*,*) 'ESMF_SWMF ERROR: ',&
          'Setting clock failed, start time=',iStartTime_I, &
          ' finish time=',iFinishTime_I,' coupling freq=',iCoupleFreq
     call my_error('ESMF_ClockCreate failed')
  end if

  if(TimeSimulation /= 0.0)then
     iSecond   = int(TimeSimulation)
     iMillisec = nint(1000*(TimeSimulation-iSecond))
     call ESMF_TimeIntervalSet(SimTime, s=iSecond, ms=iMillisec, rc=rc)
     if(rc /= ESMF_SUCCESS)then
        if(iProc == 0)write(*,*) 'ESMF_SWMF ERROR: ',&
             'Setting time interval failed for simulation time=',&
             TimeSimulation,' s=',iSecond,' ms=',iMillisec
        call my_error('ESMF_TimeIntervalSet failed')
     end if

     CurrentTime = StartTime + SimTime
     call ESMF_ClockSet(Clock, currtime = CurrentTime, rc=rc)

     if(rc /= ESMF_SUCCESS)then
        if(iProc == 0)write(*,*) 'ESMF_SWMF ERROR: ',&
             'Setting clock to Simulation time=',TimeSimulation,' failed'
        call my_error('ESMF_ClockSet failed')
     end if
  end if

  !  Init, Run, and Finalize section
  call ESMF_GridCompInitialize(EsmfSwmfComp, clock=Clock, rc=rc)
  if (rc /= ESMF_SUCCESS) call my_error('EsmfSwmfComp:init failed')

  call ESMF_GridCompRun(EsmfSwmfComp, clock=Clock, rc=rc)
  if (rc /= ESMF_SUCCESS) call my_error('EsmfSwmfComp:run failed')

  call ESMF_GridCompFinalize(EsmfSwmfComp, clock=Clock, rc=rc)
  if (rc /= ESMF_SUCCESS) call my_error('EsmfSwmfComp:finalize failed')

  ! Clean up

  call ESMF_ClockDestroy(clock, rc=rc)

  call ESMF_GridCompDestroy(EsmfSwmfComp, rc=rc)

  call ESMF_Finalize

contains
  !============================================================================
  subroutine my_error(String)

    character(len=*), intent(in) :: String

    write(*,*)'ERROR in ESMF_SWMF_Driver, iProc=', iProc
    write(*,*) String
    call ESMF_Finalize
    stop

  end subroutine my_error
  !============================================================================
end program ESMF_SWMF_Driver

