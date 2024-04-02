module SWMF_grid_comp

  ! This is the SWMF Gridded Component, which acts as an interface
  ! to the SWMF. It does not actually have a grid, but it does
  ! the initialization, running and finalization of the SWMF.

  ! ESMF Framework module
  use ESMF

  use ESMFSWMF_variables, ONLY: &
       DoRunSwmf, DoBlockAllSwmf, &
       NameSwmfComp, iProc0SwmfComp, iProcLastSwmfComp, &
       NameFieldEsmf_V, nVarEsmf, &
       Year_, Month_, Day_, Hour_, Minute_, Second_, MilliSec_, &
       write_log, write_error

  ! Size of ionosphere grid in SWMF/IE model
  use IE_ModSize, ONLY: nLat => IONO_nTheta, nLon => IONO_nPsi
  ! Conversion to radians
  use ModNumConst, ONLY: cDegToRad

  implicit none

  private

  public:: SetServices

  ! Local variables

  ! Grid related to the SWMF Ridley Ionosphere Model
  ! This is a 2D spherical grid representing the height integrated ionosphere.
  ! RIM is in SM coordinates (aligned with Sun-Earth direction):
  ! +Z points to north magnetic dipole and the Sun is in the +X-Z halfplane.

  ! Coordinate arrays
  real(ESMF_KIND_R8), pointer, save:: Lon_I(:), Lat_I(:)

contains
  !============================================================================
  subroutine SetServices(gComp, rc)

    type(ESMF_GridComp) :: gComp
    integer, intent(out):: rc

    call ESMF_GridCompSetEntryPoint(gComp, ESMF_METHOD_INITIALIZE, &
         userRoutine=my_init, rc=rc)
    call ESMF_GridCompSetEntryPoint(gComp, ESMF_METHOD_RUN, &
         userRoutine=my_run, rc=rc)
    call ESMF_GridCompSetEntryPoint(gComp, ESMF_METHOD_FINALIZE, &
         userRoutine=my_final, rc=rc)

  end subroutine SetServices
  !============================================================================
  subroutine my_init(gComp, ImportState, ExportState, ExternalClock, rc)

    type(ESMF_GridComp) :: gComp
    type(ESMF_State) :: ImportState
    type(ESMF_State) :: ExportState
    type(ESMF_Clock) :: ExternalClock
    integer, intent(out):: rc

    logical          :: IsLastSession ! true if SWMF has a single session
    type(ESMF_VM)    :: Vm
    type(ESMF_Grid)  :: Grid
    integer          :: iComm, iProc
    type(ESMF_Time)  :: StartTime
    integer          :: iStartTime_I(Year_:Millisec_)
    type(ESMF_TimeInterval) :: SimTime, RunDuration
    integer(ESMF_KIND_I4)   :: iSecond, iMilliSec
    real(ESMF_KIND_R8)      :: TimeSim, TimeStop

    integer:: i
    !--------------------------------------------------------------------------
    call write_log("SWMF_grid_comp:init routine called")
    rc = ESMF_FAILURE

    ! Obtain the VM for the SWMF gridded component
    call ESMF_GridCompGet(gComp, vm=Vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_GridCompGet')

    ! Obtain the MPI communicator for the VM
    call ESMF_VMGet(Vm, mpiCommunicator=iComm,  rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_VMGet')
    
    ! Obtain the start time from the clock 
    call ESMF_ClockGet(externalclock, startTime=StartTime, &
         currSimTime=SimTime, runDuration=RunDuration, rc=rc)
    if(rc /= ESMF_SUCCESS) call	my_error('ESMF_ClockGet')

    call ESMF_TimeGet(StartTime,   &
         yy=iStartTime_I(Year_),   &
         mm=iStartTime_I(Month_),  &
         dd=iStartTime_I(Day_),    &
         h =iStartTime_I(Hour_),   &
         m =iStartTime_I(Minute_), &
         s =iStartTime_I(Second_), &
         ms=iStartTime_I(Millisec_), &
         rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_TimeGet')

    ! Obtain the simulation time from the clock
    call ESMF_TimeIntervalGet(SimTime, s=iSecond, ms=iMillisec, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_TimeIntervalGet Sim')
    TimeSim = iSecond + iMillisec/1000.0

    ! Obtain the final simulation time from the clock
    call ESMF_TimeIntervalGet(RunDuration, s=iSecond, ms=iMillisec, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_TimeIntervalGet Run')
    TimeStop = iSecond + iMillisec/1000.0

    ! Initialze the SWMF with this MPI communicator and start time
    call write_log("SWMF_initialize routine called")
    call SWMF_initialize(iComm, iStartTime_I, &
         TimeSim, TimeStop, IsLastSession, rc)
    call write_log("SWMF_initialize routine returned")
    if(rc /= 0)call my_error('SWMF_initialize')

    rc = ESMF_SUCCESS
    call write_log("SWMF_grid_comp:init routine returned")

  end subroutine my_init
  !============================================================================
  subroutine my_run(gComp, ImportState, ExportState, Clock, rc)

    type(ESMF_GridComp):: gComp
    type(ESMF_State):: ImportState
    type(ESMF_State):: ExportState
    type(ESMF_Clock):: Clock
    integer, intent(out):: rc

    ! Access to the data
    real(ESMF_KIND_R8), pointer     :: Ptr(:,:)
    real(ESMF_KIND_R8), allocatable :: Data_VII(:,:,:)
    integer                         :: iVar

    ! Access to time
    type(ESMF_TimeInterval) :: SimTime, TimeStep
    integer(ESMF_KIND_I4)   :: iSec, iMilliSec

    ! Parameters for the SWMF_run interface
    logical            :: DoStop            ! true if SWMF requests a stop
    real(ESMF_KIND_R8) :: tCurrent, tCouple ! Current and next coupling times
    real(ESMF_KIND_R8) :: tSimSwmf          ! SWMF Simulation time

    !--------------------------------------------------------------------------
    call write_log("SWMF_grid_comp:run routine called")
    rc = ESMF_FAILURE

    ! Get the current time from the clock
    call ESMF_ClockGet(Clock, CurrSimTime=SimTime, TimeStep=TimeStep, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_ClockGet')
    call ESMF_TimeIntervalGet(SimTime, s=iSec, ms=iMilliSec, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_TimeIntervalGet current')
    tCurrent = iSec + 0.001*iMilliSec

    ! Calculate simulation time for next coupling
    SimTime = SimTime + TimeStep
    call ESMF_TimeIntervalGet(SimTime, s=iSec, ms=iMilliSec, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_TimeIntervalGet couple')
    tCouple = iSec + 0.001*iMilliSec
    
    call write_log("SWMF_run routine called!")
    write(*,*)'SWMF_run starts  with tCouple =',tCouple
    if(.not.DoRunSwmf)then
       ! Pretend that SWMF reached the coupling time
       tSimSwmf = tCouple
       rc = 0
    elseif(DoBlockAllSwmf)then
       call SWMF_run('**', tCouple, tSimSwmf, DoStop, rc)
    else
       call SWMF_run(NameSwmfComp, tCouple, tSimSwmf, DoStop, rc)
    end if
    write(*,*)'SWMF_run returns with tSimSwmf=', tSimSwmf
    call write_log("SWMF_run routine returned!")
    if(rc /= 0)call my_error('SWMF_run')
    
    call write_log("SWMF_grid_comp:run routine returned")

    rc = ESMF_SUCCESS

  end subroutine my_run
  !============================================================================
  subroutine my_final(gComp, ImportState, ExportState, ExternalClock, rc)

    type(ESMF_GridComp) :: gcomp
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    type(ESMF_Clock) :: externalclock
    integer, intent(out) :: rc

    type(ESMF_VM)    :: vm
    integer          :: iProc
    !--------------------------------------------------------------------------
    call write_log("SWMF_finalize routine called")

    call SWMF_finalize(rc)
    call write_log("SWMF_finalize routine returned")
    if(rc /= 0) call my_error("SWMF_finalize FAILED")

  end subroutine my_final
  !============================================================================
  subroutine my_error(String)

    ! Write out error message and stop
    
    character(len=*), intent(in) :: String
    !--------------------------------------------------------------------------
    call write_error('SWMF_grid_comp '//String)
    
  end subroutine my_error
  !============================================================================
end module SWMF_grid_comp
