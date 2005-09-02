!------------------------------------------------------------------------------
! MODULE ESMF_SWMF_Driver.f90 - main program for coupled ESMF-SWMF application
!
! !DESCRIPTION:
!  Application Driver for the coupled ESMF-SWMF system.  
!  Creates the top ESMF_SWMF Gridded Component and calls the 
!  Initialize, Run, and Finalize routines for it.  
!
!  The top Gridded Component creates and manages the ESMF and SWMF 
!  subcomponents internally. The SWMF is treated as a single component
!  which is coupled to (some of) the ESMF component(s) periodically.
!
!BOP
!\begin{verbatim}

program ESMF_SWMF_Driver

  ! ESMF module, defines all ESMF data types and procedures
  use ESMF_Mod

  ! Top ESMF-SWMF Gridded Component registration routines
  use ESMF_SWMF_GridCompMod, only : SetServices => ESMF_SWMF_SetServices

  ! Named indexes for integer time arrays
  use ModTimeArray

  implicit none

  ! Local variables
  character (len=*), parameter :: NameParamFile = "ESMF_SWMF.input"

  ! Components
  type(ESMF_GridComp) :: compGridded

  ! States, Virtual Machines, and Layouts and processor index
  type(ESMF_VM)       :: defaultvm
  type(ESMF_DELayout) :: defaultlayout
  type(ESMF_State)    :: defaultstate
  integer             :: iProc

  ! Configuration information
  type(ESMF_Config) :: config

  ! A common grid
  type(ESMF_Grid) :: grid

  ! Variables related to the (unused) grid
  integer :: i_max=10, j_max=10
  real(ESMF_KIND_R8) :: x_min=0.0, x_max=1.0, y_min=0.0, y_max=1.0

  ! A clock, starting and stop times and timestep
  type(ESMF_Clock)        :: clock
  type(ESMF_Time)         :: startTime
  type(ESMF_Time)         :: stopTime
  type(ESMF_TimeInterval) :: timeStep

  ! Variables for the clock
  integer :: iTime                              ! Index for date-time arrays
  integer :: iStartTime_I(Year_:MilliSec_)  = & ! Start date-time
       (/2000, 3, 21, 10, 45, 0, 0/)            !   with defaults
  integer :: iFinishTime_I(Year_:MilliSec_) = & ! Finish date-time
       (/2000, 3, 21, 10, 45, 0, 0/)            !   with defaults
  integer :: iDefaultTmp                        ! Temporary for default value

  integer :: iCoupleFreq = 1                    ! Coupling frequency in seconds

  ! Labels used in the input file
  character (len=*), parameter :: StringStart  = 'Start '
  character (len=*), parameter :: StringFinish = 'Finish '
  character (len=9), parameter :: StringTime_I(Year_:MilliSec_) = (/ &
       'Year:    ',&
       'Month:   ',&
       'Day:     ',&
       'Hour:    ',&
       'Minute:  ',&
       'Second:  ',&
       'Millisec:' /)

  ! Simulation time (non-zero for restart)
  type(ESMF_Time)         :: CurrentTime
  type(ESMF_TimeInterval) :: SimTime
  real(ESMF_KIND_R4)      :: TimeSimulation = 0.0
  integer(ESMF_KIND_I4)   :: iSecond, iMillisec

  ! Return code for error checks
  integer :: rc

  !----------------------------------------------------------------------------
  !  Initialize the ESMF Framework
  !----------------------------------------------------------------------------

  call ESMF_Initialize(defaultCalendar=ESMF_CAL_GREGORIAN, rc=rc)
  if (rc /= ESMF_SUCCESS) stop 'ESMF_Initialize FAILED'


  call ESMF_LogWrite("ESMF-SWMF Driver start", ESMF_LOG_INFO)

  ! Get the default VM which contains all PEs this job was started on.
  call ESMF_VMGetGlobal(defaultvm, rc)

  ! Get the processor number
  call ESMF_VMGet(defaultvm, localPET=iProc)

  !
  ! Read in Configuration information from a default config file
  !

  config = ESMF_ConfigCreate(rc)

  call ESMF_ConfigLoadFile(config, NameParamFile, rc = rc)
  if(rc /= ESMF_SUCCESS) then
     if(iProc == 0)write(*,*) &
          'ESMF_ConfigLoadFile FAILED for file '//NameParamFile
     call ESMF_Finalize
  endif

  ! Get the start time, stop time, and running intervals
  ! for the main time loop.
  !
  do iTime = Year_, Millisec_
     iDefaultTmp = iStartTime_I(iTime)
     call ESMF_ConfigGetAttribute(config, iStartTime_I(iTime),&
          StringStart//trim(StringTime_I(iTime)), rc=rc)
     if(rc /= ESMF_SUCCESS) then
        if(iProc == 0) write(*,*) &
             'Did not read ',StringStart//trim(StringTime_I(iTime)), &
             ' setting default value= ', iDefaultTmp
        iStartTime_I(iTime) = iDefaultTmp
     end if
  end do
  do iTime = Year_, Millisec_
     call ESMF_ConfigGetAttribute(config, iFinishTime_I(iTime),&
          StringFinish//trim(StringTime_I(iTime)), rc=rc)

     if(rc /= ESMF_SUCCESS) then
        if(iProc == 0)write(*,*) &
             'Did not read ',StringFinish//trim(StringTime_I(iTime)), &
             ' setting default value= ', iDefaultTmp
        iStartTime_I(iTime) = iDefaultTmp
     end if
  end do

  call ESMF_ConfigGetAttribute(config, &
       TimeSimulation, 'Simulation Time:', rc=rc)
  if(rc /= ESMF_SUCCESS) then
     if(iProc == 0)write(*,*) &
          'Did not read Simulation Time: setting default value = 0.0'
     TimeSimulation = 0.0
  end if

  iDefaultTmp = iCoupleFreq
  call ESMF_ConfigGetAttribute(config, &
       iCoupleFreq, 'Coupling Frequency:', rc=rc)
  if(rc /= ESMF_SUCCESS) then
     if(iProc == 0)write(*,*) &
          'Did not read Coupling Frequency: setting default value= ', &
          iDefaultTmp
     iCoupleFreq = iDefaultTmp
  end if

  !----------------------------------------------------------------------------
  !    Create section
  !----------------------------------------------------------------------------

  ! Create the top Gridded component, passing in the default layout.
  compGridded = ESMF_GridCompCreate(defaultvm, "ESMF Gridded Component", rc=rc)

  call ESMF_LogWrite("Component Create finished", ESMF_LOG_INFO)


  !----------------------------------------------------------------------------
  !  Register section
  !----------------------------------------------------------------------------

  call ESMF_GridCompSetServices(compGridded, SetServices, rc)
  if (ESMF_LogMsgFoundError(rc, "Registration failed", rc)) goto 10

  !----------------------------------------------------------------------------
  !  Create and initialize a clock, and a grid.
  !----------------------------------------------------------------------------

  ! Based on values from the Config file, create a Clock.  

  call ESMF_TimeIntervalSet(timeStep, s=iCoupleFreq, rc=rc)

  if(rc /= ESMF_SUCCESS) then
     if(iProc == 0) write(*,*)'Setting coupling frequency failed:',iCoupleFreq
     call ESMF_Finalize
  end if

  call ESMF_TimeSet(startTime, &
       yy=iStartTime_I(Year_), &
       mm=iStartTime_I(Month_), &
       dd=iStartTime_I(Day_), &
       h =iStartTime_I(Hour_), &
       m =iStartTime_I(Minute_), &
       s =iStartTime_I(Second_), &
       ms=iStartTime_I(Millisec_),&
       rc=rc)

  if(rc /= ESMF_SUCCESS) then
     if(iProc == 0)write(*,*)'Setting start time failed:',iStartTime_I
     call ESMF_Finalize
  end if

  call ESMF_TimeSet(stopTime, &
       yy=iFinishTime_I(Year_), &
       mm=iFinishTime_I(Month_), &
       dd=iFinishTime_I(Day_), &
       h =iFinishTime_I(Hour_), &
       m =iFinishTime_I(Minute_), &
       s =iFinishTime_I(Second_), &
       ms=iFinishTime_I(Millisec_))

  if(rc /= ESMF_SUCCESS) then
     if(iProc == 0)write(*,*)'Setting finish time failed:',iFinishTime_I
     call ESMF_Finalize
  end if

  clock = ESMF_ClockCreate("Application Clock", timeStep, startTime, &
       stopTime, rc=rc)

  if(rc /= ESMF_SUCCESS)then
     if(iProc == 0)write(*,*)&
          'Setting clock failed, start time=',iStartTime_I, &
          ' finish time=',iFinishTime_I,' coupling freq=',iCoupleFreq
     call ESMF_Finalize
  end if

  if(TimeSimulation /= 0.0)then
     iSecond   = int(TimeSimulation)
     iMillisec = nint(1000*(TimeSimulation-iSecond))
     call ESMF_TimeIntervalSet(SimTime, s=iSecond, ms=iMillisec, rc=rc)
     if(rc /= ESMF_SUCCESS)then
        if(iProc == 0)write(*,*)&
             'Setting time interval failed for simulation time=',&
             TimeSimulation,' s=',iSecond,' ms=',iMillisec
        call ESMF_Finalize
     end if

     CurrentTime = StartTime + SimTime
     call ESMF_ClockSet(clock, currtime = CurrentTime, rc=rc)

     if(rc /= ESMF_SUCCESS)then
        if(iProc == 0)write(*,*)&
             'Setting clock to Simulation time=',TimeSimulation,' failed'
        call ESMF_Finalize
     end if
  end if

  ! Setup a grid. This is not used for any actual data transfer yet.
  ! Get a default layout based on the VM.
  defaultlayout = ESMF_DELayoutCreate(defaultvm, rc=rc)

  grid = ESMF_GridCreateHorzXYUni(counts=(/i_max, j_max/), &
       minGlobalCoordPerDim=(/x_min, y_min/), &
       maxGlobalCoordPerDim=(/x_max, y_max/), &
       name="ESMF-SWMF common grid", rc=rc)
  call ESMF_GridDistribute(grid, delayout=defaultlayout, rc=rc)

  ! Attach the Grid to the Component
  call ESMF_GridCompSet(compGridded, grid=grid, rc=rc)


  !----------------------------------------------------------------------------
  !  Create and initialize a State to use for both import and export.
  !----------------------------------------------------------------------------

  defaultstate = ESMF_StateCreate("Default Gridded State", rc=rc)

  !----------------------------------------------------------------------------
  !  Init, Run, and Finalize section
  !----------------------------------------------------------------------------

  call ESMF_GridCompInitialize(compGridded, defaultstate, defaultstate, &
       clock, rc=rc)
  if (ESMF_LogMsgFoundError(rc, "Initialize failed", rc)) goto 10

  call ESMF_GridCompRun(compGridded, defaultstate, defaultstate, &
       clock, rc=rc)
  if (ESMF_LogMsgFoundError(rc, "Run failed", rc)) goto 10

  call ESMF_GridCompFinalize(compGridded, defaultstate, defaultstate, &
       clock, rc=rc)
  if (ESMF_LogMsgFoundError(rc, "Finalize failed", rc)) goto 10


  !----------------------------------------------------------------------------
  !     Destroy section
  !----------------------------------------------------------------------------

  ! Clean up

  call ESMF_ClockDestroy(clock, rc)

  call ESMF_StateDestroy(defaultstate, rc)

  call ESMF_GridCompDestroy(compGridded, rc)

  call ESMF_DELayoutDestroy(defaultLayout, rc)

  !----------------------------------------------------------------------------

10 continue

  call ESMF_Finalize(rc=rc)

end program ESMF_SWMF_Driver

!\end{verbatim}    
!EOP
