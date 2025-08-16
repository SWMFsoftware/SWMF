module ESMFSWMF_variables

  ! Various entities needed for the ESMF-SWMF coupling

  use ESMF
  use NUOPC
  implicit none

  private
  public:: read_esmf_swmf_input, add_fields, write_log, write_error

  ! main configuration file
  character (len=*), public, parameter :: NameParamFile = "ESMF_SWMF.input"

  ! named integer indexes for integer time arrays
  integer, public, parameter :: &
       Year_=1, Month_=2, Day_=3, Hour_=4, Minute_=5, Second_=6, MilliSec_=7

  ! number of ESMF variables and their names to be sent to SWMF
  integer, public, parameter :: nVarEsmf = 2
  character(len=4), public, parameter :: NameFieldEsmf_V(nVarEsmf) = &
       [ 'Hall', 'Ped ' ]

  ! number of SWMF variables and their names to be sent to ESMF
  integer, public, parameter :: nVarSwmf = 4
  character(len=4), public, parameter :: NameFieldSwmf_V(nVarSwmf) = &
       [ 'jFac', 'Epot', 'Aver', 'Diff' ]

  ! Time related variables
  integer, public:: iStartTime_I(Year_:MilliSec_)  = & ! Start date-time
       [2000, 3, 21, 10, 45, 0, 0]                     !   with defaults
  integer, public:: iFinishTime_I(Year_:MilliSec_) = & ! Finish date-time
       [2000, 3, 21, 10, 45, 0, 0]                     !   with defaults
  real(ESMF_KIND_R8), public:: TimeSimulation = 0.0

  ! Root processor and number of processors on global VM
  integer, public :: iProc=0, nProc=1

  ! SWMF component to couple with
  character(len=2), public :: NameSwmfComp = 'IE'

  ! The ESMF communicates with SwmfComp within the SWMF.
  ! The processors used by SwmfComp are obtained from the PARAM.in file.
  ! These indexes are relative to the SWMF MPI communicator
  integer, public:: iProc0SwmfComp=0, iProcLastSwmfComp=0

  ! Testing
  logical, public :: DoTest = .false.

  ! Field values and coordinate coefficients for testing
  real, public, parameter:: FieldTest_V(nVarEsmf) = [3.0, 5.0]
  real, public, parameter:: CoordCoefTest = 0.1

  ! Change of Hall field during ESMF run
  real, public, parameter:: dHallPerDtTest = 0.4

contains
  !============================================================================
  subroutine read_esmf_swmf_input(iError)

    integer, intent(out) :: iError    ! error code

    ! Labels used in the input file for the start and finish times
    character (len=*), parameter :: StringStart  = 'Start '
    character (len=*), parameter :: StringFinish = 'Finish '
    character (len=9), parameter :: StringTime_I(Year_:MilliSec_) = [ &
         'Year:    ',&
         'Month:   ',&
         'Day:     ',&
         'Hour:    ',&
         'Minute:  ',&
         'Second:  ',&
         'Millisec:' ]

    integer :: iTime

    type(ESMF_Config) :: Config
    character (len=*), parameter :: NameParamFile = "ESMF_SWMF.input"

    ! Store default values before reading
    integer               :: iDefaultTmp          ! Temporary default integer
    character             :: StringTmp            ! Temporary string
    !--------------------------------------------------------------------------
    call write_log("ESMF_SWMF_Mod:read_esmf_swmf_input called")
    iError = ESMF_FAILURE

    ! Read in Configuration information from NameParamFile
    Config = ESMF_ConfigCreate(rc=iError)
    if(iError /= ESMF_SUCCESS) RETURN

    call ESMF_ConfigLoadFile(Config, NameParamFile, rc=iError)
    if(iError /= ESMF_SUCCESS) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF ERROR: ',&
            'ESMF_ConfigLoadFile FAILED for file '//trim(NameParamFile)
       RETURN
    endif

    ! Get start date/time
    do iTime = Year_, Millisec_
       iDefaultTmp = iStartTime_I(iTime)
       call ESMF_ConfigGetAttribute(Config, iStartTime_I(iTime),&
            label=StringStart//trim(StringTime_I(iTime)), rc=iError)
       if(iError /= ESMF_SUCCESS) then
          if(iProc == 0) write(*,*) 'ESMF_SWMF did not read ', &
               StringStart//trim(StringTime_I(iTime)), &
               ' setting default value= ', iDefaultTmp
          iStartTime_I(iTime) = iDefaultTmp
       end if
    end do

    ! Get stop date/time
    do iTime = Year_, Millisec_
       call ESMF_ConfigGetAttribute(Config, iFinishTime_I(iTime),&
            label=StringFinish//trim(StringTime_I(iTime)), rc=iError)
       if(iError /= ESMF_SUCCESS) then
          if(iProc == 0)write(*,*) 'ESMF_SWMF did not read ',&
               StringFinish//trim(StringTime_I(iTime)), &
               ' setting default value= ', iDefaultTmp
          iStartTime_I(iTime) = iDefaultTmp
       end if
    end do

    ! Get current simulation time (non-zero for restart)
    call ESMF_ConfigGetAttribute(Config, &
         TimeSimulation, label='Simulation Time:', rc=iError)
    if(iError /= ESMF_SUCCESS) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF did not read ',&
            'Simulation Time: setting default value = 0.0'
       TimeSimulation = 0.0
    end if
    if(TimeSimulation < 0.0)then
       if(iProc == 0) write(*,*) 'ESMF_SWMF ERROR: ', &
            'TimeSimulation =',TimeSimulation,' should not be negative!'
       call ESMF_Finalize
    end if

    ! Read the SWMF component name that is coupled
    call ESMF_ConfigGetAttribute(Config, NameSwmfComp, &
         label='SWMF Component:', rc=iError)
    if(iError /= ESMF_SUCCESS) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF: ', &
            'Setting default for SWMF Component: IE'
       NameSwmfComp = 'IE'
    end if

    ! Read testing option
    call ESMF_ConfigGetAttribute(Config, DoTest, &
         label='DoTest:', rc=iError)
    if(iError /= ESMF_SUCCESS) then
       if(iProc == 0) write(*,*) 'ESMF_SWMF did not read ', &
            'DoTest: setting default value = true'
       DoTest = .true.
    else
       if(iProc == 0) write(*,*) 'ESMF_SWMF ', &
            'DoTest is set to ', DoTest
    end if

    call ESMF_ConfigDestroy(Config, rc=iError)
    if(iError /= ESMF_SUCCESS) RETURN

    call read_swmf_layout(iError)
    if(iError /= ESMF_SUCCESS) RETURN

    iError = ESMF_SUCCESS
    call write_log("ESMF_SWMF_Mod:read_esmf_swmf_input returned")

  end subroutine read_esmf_swmf_input
  !============================================================================
  subroutine read_swmf_layout(iError)

    ! Get the processors of the SWMF component to be coupled with
    ! from the PARAM.in file

    integer, intent(out) :: iError ! error code

    integer :: iUnit
    character(len=100) :: String
    logical :: DoRead
    !--------------------------------------------------------------------------
    call CON_io_unit_new_ext(iUnit)
    open(unit=iUnit, file='PARAM.in', iostat=iError)
    if(iError /= 0)then
       if(iProc==0)write(*,*)'ESMF_SWMF ERROR: '// &
            'could not open PARAM.in file'
       iError = ESMF_FAILURE; RETURN
    end if

    ! Read the PARAM.in file
    DoRead = .false.
    READLAYOUT: do
       read(iUnit,'(a)',iostat=iError) String

       if(iError /= 0)then
          if(iProc==0)write(*,*)'ESMF_SWMF ERROR: '// &
               'could not read PARAM.in file'
          iError = ESMF_FAILURE; RETURN
       end if
       if(String(1:13) == '#COMPONENTMAP') DoRead = .true.
       if(.not.DoRead) CYCLE
       if(String(1:2) == NameSwmfComp)then
          read(String, *, iostat=iError)NameSwmfComp, &
               iProc0SwmfComp, iProcLastSwmfComp
          if(iError/=0)then
             if(iProc==0)write(*,*)'ESMF_SWMF ERROR: '// &
                  'could not read iProcRoot for '//NameSwmfComp// &
                  ' from line=',trim(String),' in PARAM.in'
             iError = ESMF_FAILURE; RETURN
          end if
          EXIT READLAYOUT
       end if
       if(String == '')then
          if(iProc==0)write(*,*)'ESMF_SWMF ERROR: '// &
               'could not find component '//NameSwmfComp//' in #COMPONENTMAP'
          iError = ESMF_FAILURE; RETURN
       end if
    end do READLAYOUT
    close(iUnit)

    if(iProc==0) write(*,*)'ESMF_SWMF: '//NameSwmfComp//' Root, Last=', &
         iProc0SwmfComp, iProcLastSwmfComp

  end subroutine read_swmf_layout
  !============================================================================
  subroutine add_fields(Grid, State, IsFromEsmf, iError)

    type(ESMF_Grid), intent(inout) :: Grid
    type(ESMF_State),    intent(inout) :: State
    logical,             intent(in)    :: IsFromEsmf
    integer,             intent(out)   :: iError

    type(ESMF_Field)     :: Field
    type(ESMF_ArraySpec) :: ArraySpec
    integer              :: iVar
    character(len=4)     :: NameField
    ! real(ESMF_KIND_R8), pointer :: Ptr_C(:,:)
    type(ESMF_VM)      :: Vm

    ! Name of the component
    character(len=100) :: Name="UNKNOWN"
    !--------------------------------------------------------------------------
    call write_log("ESMF_SWMF_Mod:add_fields called")
    iError = ESMF_FAILURE

    call ESMF_VMGetCurrent(Vm, rc=iError)
    if (iError /= ESMF_SUCCESS) call my_error('ESMF_VMGetCurrent')
    call ESMF_VMGet(Vm, localPet=iProc, rc=iError)
    if (iError /= ESMF_SUCCESS) call my_error('ESMF_VMGet')

    call ESMF_GridGet(Grid, name=Name)
    if (iError /= ESMF_SUCCESS) call my_error('ESMF_GridGet')

    call ESMF_ArraySpecSet(ArraySpec, rank=2, typekind=ESMF_TYPEKIND_R8)
    if (iError /= ESMF_SUCCESS) call my_error('ESMF_ArraySpecSet')

    if(IsFromEsmf)then
       do iVar = 1, nVarEsmf
          NameField = NameFieldEsmf_V(iVar)
          if (NUOPC_IsConnected(State, fieldName=trim(NameField))) then
             write(*,*) iProc,' Adding ESMF field=', NameField,' to ',trim(Name)
             Field = ESMF_FieldCreate(Grid, arrayspec=ArraySpec, &
                  staggerloc=ESMF_STAGGERLOC_CORNER, name=NameField, rc=iError)
             if(iError /= ESMF_SUCCESS) call my_error('ESMF_FieldCreate ' &
                  //NameField//' for '//trim(Name))
             call NUOPC_Realize(State, field=Field, rc=iError)
             if(iError /= ESMF_SUCCESS) call my_error( &
                  'NUOPC_Realize '//NameField//' to '//trim(Name))
          else
             call ESMF_StateRemove(State, [ trim(NameField) ], rc=iError)
             if(iError /= ESMF_SUCCESS) call my_error( &
                  'ESMF_StateRemove '//NameField)
          end if
       end do
    else
       do iVar = 1, nVarSwmf
          NameField = NameFieldSwmf_V(iVar)
          if (NUOPC_IsConnected(State, fieldName=trim(NameField))) then
             write(*,*)'Adding SWMF field=', NameField,' to ',trim(Name)
             Field = ESMF_FieldCreate(Grid, typekind=ESMF_TYPEKIND_R8, &
                  staggerloc=ESMF_STAGGERLOC_CORNER, name=NameField, rc=iError)
             if(iError /= ESMF_SUCCESS) call my_error('ESMF_FieldCreate ' &
                  //trim(NameField)//' for '//trim(Name))
             call NUOPC_Realize(State, field=Field, rc=iError)
             if(iError /= ESMF_SUCCESS) call my_error( &
                  'NUOPC_Realize '//NameField//' to '//trim(Name))
          else
             call ESMF_StateRemove(State, [ trim(NameField) ], rc=iError)
             if(iError /= ESMF_SUCCESS) call my_error( &
                  'ESMF_StateRemove '//NameField)
          end if
       end do
    endif
    iError = ESMF_SUCCESS
    call write_log("ESMF_SWMF_Mod:add_fields returned")

  end subroutine add_fields
  !============================================================================
  subroutine my_error(String)

    ! Write out error message and stop

    character(len=*), intent(in) :: String
    !--------------------------------------------------------------------------
    call write_error('ESMFSWMF_variables '//String)

  end subroutine my_error
  !============================================================================
  subroutine write_log(String)

    ! write out the log message in String and flush the file

    character(len=*), intent(in):: String
    !--------------------------------------------------------------------------
    call ESMF_LogWrite(String, ESMF_LOGMSG_INFO)
    call ESMF_LogFlush()

  end subroutine write_log
  !============================================================================
  subroutine write_error(String)

    ! Write out processor index and the error message String and stop

    character(len=*), intent(in):: String
    !--------------------------------------------------------------------------
    write(*,'(a,i4,a,a)') 'ERROR (iProc=',iProc,'):', String
    call ESMF_finalize
    stop

  end subroutine write_error
  !============================================================================
end module ESMFSWMF_variables
!==============================================================================
