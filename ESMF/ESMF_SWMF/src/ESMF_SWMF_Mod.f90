module ESMF_SWMF_Mod

  ! Various entities needed for the ESMF-SWMF coupling

  use ESMF
  implicit none

  private
  public:: read_esmf_swmf_input, add_fields, log_write

  type(ESMF_Sync_Flag), public:: SyncFlag = ESMF_SYNC_BLOCKING

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
  integer, public :: iCoupleFreq = 1   ! Coupling frequency in seconds

  ! SWMF runs on processor ranks iProcRootSwmf to iProcLastSwmf
  ! ESMF runs on processor ranks iProcRootEsmf to iProcLastEsmf
  integer, public :: iProcRootSwmf=1, iProcLastSwmf=-1, nProcSwmf
  integer, public :: iProcRootEsmf=0, iProcLastEsmf=0

  ! SWMF component to couple with
  character(len=2), public :: NameSwmfComp = 'IE'

  ! The ESMF communicates with SwmfComp within the SWMF.
  ! The processors used by SwmfComp are obtained from the PARAM.in file.
  ! These indexes are relative to the SWMF MPI communicator
  integer, public:: iProc0SwmfComp=0, iProcLastSwmfComp=0, nProcSwmfComp=1

  ! If DoRunSwmf is true, run the SWMF for real, otherwise just pretend
  logical, public:: DoRunSwmf = .true.

  ! When SWMF and ESMF are coupled, one can block the whole SWMF
  ! or only the component the ESMF is communicating with. The latter
  ! is more efficient but it can result in a dead-lock if the ESMF and
  ! SWMF overlap.
  logical, public:: DoBlockAllSwmf=.false.

  ! Testing
  logical, public, parameter:: DoTest = .true.

  ! Field values and coordinate coefficients for testing
  real, public, parameter:: FieldTest_V(nVarEsmf) = [3.0, 5.0]
  real, public, parameter:: CoordCoefTest = 0.1

  ! Change of Hall field during ESMF run
  real, public, parameter:: dHallPerDtTest = 0.4
  
contains
  !============================================================================
  subroutine read_esmf_swmf_input(nProc, iProc, rc)

    integer, intent(in)  :: nProc ! total number of processors
    integer, intent(in)  :: iProc ! processor rank for this CPU
    integer, intent(out) :: rc    ! error code

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
    real(ESMF_KIND_R8)    :: DefaultTmp           ! Temporary default real
    character             :: StringTmp            ! Temporary string
    !--------------------------------------------------------------------------
    call log_write("ESMF_SWMF_Mod:read_esmf_swmf_input called")
    rc = ESMF_FAILURE

    ! Read in Configuration information from NameParamFile
    Config = ESMF_ConfigCreate(rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    call ESMF_ConfigLoadFile(Config, NameParamFile, rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF ERROR: ',&
            'ESMF_ConfigLoadFile FAILED for file '//trim(NameParamFile)
       RETURN
    endif

    ! Get start date/time
    do iTime = Year_, Millisec_
       iDefaultTmp = iStartTime_I(iTime)
       call ESMF_ConfigGetAttribute(Config, iStartTime_I(iTime),&
            label=StringStart//trim(StringTime_I(iTime)), rc=rc)
       if(rc /= ESMF_SUCCESS) then
          if(iProc == 0) write(*,*) 'ESMF_SWMF did not read ', &
               StringStart//trim(StringTime_I(iTime)), &
               ' setting default value= ', iDefaultTmp
          iStartTime_I(iTime) = iDefaultTmp
       end if
    end do

    ! Get stop date/time
    do iTime = Year_, Millisec_
       call ESMF_ConfigGetAttribute(Config, iFinishTime_I(iTime),&
            label=StringFinish//trim(StringTime_I(iTime)), rc=rc)
       if(rc /= ESMF_SUCCESS) then
          if(iProc == 0)write(*,*) 'ESMF_SWMF did not read ',&
               StringFinish//trim(StringTime_I(iTime)), &
               ' setting default value= ', iDefaultTmp
          iStartTime_I(iTime) = iDefaultTmp
       end if
    end do

    ! Get current simulation time (non-zero for restart)
    call ESMF_ConfigGetAttribute(Config, &
         TimeSimulation, label='Simulation Time:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF did not read ',&
            'Simulation Time: setting default value = 0.0'
       TimeSimulation = 0.0
    end if
    if(TimeSimulation < 0.0)then
       if(iProc == 0) write(*,*) 'ESMF_SWMF ERROR: ', &
            'TimeSimulation =',TimeSimulation,' should not be negative!'
       call ESMF_Finalize
    end if

    ! Get coupling frequency
    iDefaultTmp = iCoupleFreq
    call ESMF_ConfigGetAttribute(Config, &
         iCoupleFreq, label='Coupling Frequency:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF did not read ',&
            'Coupling Frequency: setting default value= ', &
            iDefaultTmp
       iCoupleFreq = iDefaultTmp
    end if
    if(iCoupleFreq < 1)then
       if(iProc == 0) write(*,*) 'ESMF_SWMF ERROR: ', &
            'iCoupleFreq =', iCoupleFreq,' should be positive!'
       call ESMF_Finalize
    end if

    ! Read in layout information
    config = ESMF_ConfigCreate(rc=rc)
    call ESMF_ConfigLoadFile(Config, NameParamFile, rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF ERROR: ', &
            'ESMF_ConfigLoadFile FAILED for file '//NameParamFile
       RETURN
    endif

    ! Read root PE for the SWMF (default is 1)
    call ESMF_ConfigGetAttribute(Config, iProcRootSwmf, &
         label='SWMF Root PE:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       iProcRootSwmf = min(nProc - 1, 1)
    elseif(iProcRootSwmf < 0) then
       iProcRootSwmf = max(0, iProcRootSwmf + nProc)
    else
       iProcRootSwmf = min(iProcRootSwmf, nProc - 1)
    end if

    ! Read last PE for the SWMF (default is nProc - 1)
    call ESMF_ConfigGetAttribute(Config, iProcLastSwmf, &
         label='SWMF Last PE:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       iProcLastSwmf = nProc-1
    elseif(iProcLastSwmf < 0) then
       iProcLastSwmf = max(iProcRootSwmf, iProcLastSwmf + nProc)
    else
       iProcLastSwmf = max(iProcRootSwmf, min(iProcLastSwmf, nProc - 1))
    end if

    ! Number of PEs used by the SWMF
    nProcSwmf = iProcLastSwmf - iProcRootSwmf + 1

    ! Read root PE for the ESMF (default is 0)
    call ESMF_ConfigGetAttribute(Config, iProcRootEsmf, &
         label='ESMF Root PE:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       iProcRootEsmf = 0
    elseif(iProcRootEsmf < 0) then
       iProcRootEsmf = max(0, iProcRootEsmf + nProc)
    else
       iProcRootEsmf = min(nProc - 1, iProcRootEsmf)
    end if

    ! Read last PE for the ESMF (default is iProcRootEsmf)
    call ESMF_ConfigGetAttribute(Config, iProcLastEsmf, &
         label='ESMF Last PE:', rc=rc)
    if(rc /= ESMF_SUCCESS)then
       iProcLastEsmf = iProcRootEsmf
    elseif(iProcLastEsmf < 0)then
       iProcLastEsmf = max(iProcRootEsmf, iProcLastEsmf + nProc)
    else
       iProcLastEsmf = max(iProcRootEsmf, min(nProc - 1, iProcLastEsmf))
    end if

    write(*,*)'iProcRootEsmf, iProcLastEsmf=', iProcRootEsmf, iProcLastEsmf
    write(*,*)'iProcRootSwmf, iProcLastSwmf=', iProcRootSwmf, iProcLastSwmf

    call ESMF_ConfigGetAttribute(Config, StringTmp, &
         label='Run the SWMF [y/n]:', rc=rc)
    if(rc == ESMF_SUCCESS) then
       DoRunSwmf = StringTmp == 'y' .or. StringTmp == 'Y' .or. &
            StringTmp == 't' .or. StringTmp == 'T'
    endif
    write(*,*)'DoRunSwmf=', DoRunSwmf
    
    call ESMF_ConfigGetAttribute(Config, StringTmp, &
         label='Block all SWMF [y/n]:', rc=rc)
    
    if(rc == ESMF_SUCCESS) then
       DoBlockAllSwmf = StringTmp == 'y' .or. StringTmp == 'Y' .or. &
            StringTmp == 't' .or. StringTmp == 'T'
       write(*,*)'Block all SWMF [y/n]=',StringTmp,'! DoBlock=',DoBlockAllSwmf
    else
       ! If there is an overlap between the SWMF and ESMF processors
       ! it is a good idea to block the whole SWMF during coupling
       if(iProcRootSwmf <= iProcRootEsmf)then
          DoBlockAllSwmf = iProcLastSwmf >= iProcRootEsmf
       else
          DoBlockAllSwmf = iProcLastEsmf >= iProcRootSwmf
       end if
       if(iProc == 0)then
          if(DoBlockAllSwmf)then
             write(*,*) 'ESMF_SWMF: ', &
                  'ESMF and SWMF layouts overlap, ',&
                  'setting DoBlockAllSwmf=T'
          else
             write(*,*) 'ESMF_SWMF: ', &
                  'ESMF and SWMF layouts are disjoint, ',&
                  'setting DoBlockAllSwmf=F'
          end if
       end if
    end if

    ! Read the SWMF component name that is coupled
    call ESMF_ConfigGetAttribute(Config, NameSwmfComp, &
         label='SWMF Component:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF: ', &
            'Setting default for SWMF Component: IE'
       NameSwmfComp = 'IE'
    end if

    call ESMF_ConfigDestroy(Config, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    call read_swmf_layout(iProc, rc)
    if(rc /= ESMF_SUCCESS) RETURN

    rc = ESMF_SUCCESS
    call log_write("ESMF_SWMF_Mod:read_esmf_swmf_input returned")

  end subroutine read_esmf_swmf_input
  !============================================================================
  subroutine read_swmf_layout(iProc, rc)

    ! Get the root processor for the SWMF component to be coupled with
    ! from the PARAM.in file

    integer, intent(in)  :: iProc  ! rank of processor
    integer, intent(out) :: rc     ! error code

    integer :: iUnit
    character(len=100) :: String
    logical :: DoRead
    !--------------------------------------------------------------------------
    call CON_io_unit_new_ext(iUnit)
    open(unit=iUnit, file='PARAM.in', iostat=rc)
    if(rc /= 0)then
       if(iProc==0)write(*,*)'ESMF_SWMF ERROR: '// &
            'could not open PARAM.in file'
       rc = ESMF_FAILURE; RETURN
    end if

    ! Read the PARAM.in file
    DoRead = .false.
    READLAYOUT: do
       read(iUnit,'(a)',iostat=rc) String
    
       if(rc /= 0)then
          if(iProc==0)write(*,*)'ESMF_SWMF ERROR: '// &
               'could not read PARAM.in file'
          rc = ESMF_FAILURE; RETURN
       end if
       if(String(1:13) == '#COMPONENTMAP') DoRead = .true.
       if(.not.DoRead) CYCLE
       if(String(1:2) == NameSwmfComp)then
          read(String, *, iostat=rc)NameSwmfComp, &
               iProc0SwmfComp, iProcLastSwmfComp
          if(rc/=0)then
             if(iProc==0)write(*,*)'ESMF_SWMF ERROR: '// &
                  'could not read iProcRoot for '//NameSwmfComp// &
                  ' from line=',trim(String),' in PARAM.in'
             rc = ESMF_FAILURE; RETURN
          end if
          if(iProc0SwmfComp < 0) then
             ! Negative value is relative to the end
             iProc0SwmfComp = max(0, iProc0SwmfComp + nProcSwmf)
          else
             ! Positive value is limited by nProcSwmf
             iProc0SwmfComp = min(0, nProcSwmf - 1)
          end if
          if(iProcLastSwmfComp < 0) then
             iProcLastSwmfComp = iProcLastSwmfComp + nProcSwmf
          else
             iProcLastSwmfComp = min(iProcLastSwmfComp, nProcSwmf - 1)
          end if
          ! Number of processors used by SwmfComp
          nProcSwmfComp = iProcLastSwmfComp - iProc0SwmfComp + 1
          EXIT READLAYOUT
       end if
       if(String == '')then
          if(iProc==0)write(*,*)'ESMF_SWMF ERROR: '// &
               'could not find component '//NameSwmfComp//' in #COMPONENTMAP'
          rc = ESMF_FAILURE; RETURN
       end if
    end do READLAYOUT
    close(iUnit)

    if(iProc == 0)write(*,*)'ESMF_SWMF: '//NameSwmfComp// &
         ' Root, Last, nProc=', &
         iProc0SwmfComp, iProcLastSwmfComp, nProcSwmfComp

  end subroutine read_swmf_layout
  !============================================================================
  subroutine add_fields(Grid, State, IsFromEsmf, rc)

    type(ESMF_Grid), intent(inout) :: Grid
    type(ESMF_State),    intent(inout) :: State
    logical,             intent(in)    :: IsFromEsmf
    integer,             intent(out)   :: rc

    type(ESMF_Field)     :: Field
    type(ESMF_ArraySpec) :: ArraySpec
    integer              :: iVar, iPet
    character(len=4)     :: NameField
    ! real(ESMF_KIND_R8), pointer :: Ptr_C(:,:)
    type(ESMF_VM)      :: Vm
    
    ! Name of the component
    character(len=100) :: Name="UNKNOWN"
    !--------------------------------------------------------------------------
    call log_write("ESMF_SWMF_Mod:add_fields called")
    rc = ESMF_FAILURE

    call ESMF_VMGetCurrent(Vm, rc=rc)
    if (rc /= ESMF_SUCCESS) call my_error('ESMF_VMGetCurrent')
    call ESMF_VMGet(Vm, localPet=iPet, rc=rc)
    if (rc /= ESMF_SUCCESS) call my_error('ESMF_VMGet')

    call ESMF_GridGet(Grid, name=Name)
    if (rc /= ESMF_SUCCESS) call my_error('ESMF_GridGet')
    
    call ESMF_ArraySpecSet(ArraySpec, rank=2, typekind=ESMF_TYPEKIND_R8)
    if (rc /= ESMF_SUCCESS) call my_error('ESMF_ArraySpecSet')
    
    if(IsFromEsmf)then
       do iVar = 1, nVarEsmf
          NameField = NameFieldEsmf_V(iVar)
          write(*,*) iPet,' Adding ESMF field=', NameField,' to ',trim(Name)
          Field = ESMF_FieldCreate(Grid, arrayspec=ArraySpec, &
               staggerloc=ESMF_STAGGERLOC_CORNER, name=NameField, rc=rc)
          if(rc /= ESMF_SUCCESS) call my_error('ESMF_FieldCreate ' &
               //NameField//' for '//trim(Name))
          ! nullify(Ptr_C)
          ! call ESMF_FieldGet(Field, farrayPtr=Ptr_C, rc=rc)
          ! if(rc /= ESMF_SUCCESS) call my_error("ESMF_FieldGet for "//NameField// &
          !      ' for '//trim(Name))
          ! if(size(Ptr_C) == 0)then
          !    write(*,*)'!!! size(Ptr_C)=', size(Ptr_C)
          !    call my_error("ESMF_FieldGet "//NameField//' for  '//trim(Name))
          ! end if
          ! Ptr_C = -7.77
          call ESMF_StateAdd(State, [Field], rc=rc)
          if(rc /= ESMF_SUCCESS) call my_error('ESMF_StateAdd '//NameField// &
                ' to '//trim(Name))
       end do
    else
       do iVar = 1, nVarSwmf
          NameField = NameFieldSwmf_V(iVar)
          write(*,*)'Adding SWMF field=', NameField,' to ',trim(Name)
          Field = ESMF_FieldCreate(Grid, typekind=ESMF_TYPEKIND_R8, &
               staggerloc=ESMF_STAGGERLOC_CORNER, name=NameField, rc=rc)
          if(rc /= ESMF_SUCCESS) call my_error('ESMF_FieldCreate ' &
               //trim(NameField)//' for '//trim(Name))
          call ESMF_StateAdd(State, [Field], rc=rc)
          if(rc /= ESMF_SUCCESS) call my_error('ESMF_StateAdd ' &
               //trim(NameField)//' to '//trim(Name))
       end do
    endif
    rc = ESMF_SUCCESS
    call log_write("ESMF_SWMF_Mod:add_fields returned")

  end subroutine add_fields
  !============================================================================
  subroutine my_error(String)

    character(len=*), intent(in) :: String

    write(*,*)'ERROR in ESMF_SWMF_Mod:', String
    call ESMF_Finalize
    stop

  end subroutine my_error
  !============================================================================
  subroutine log_write(String)

    character(len=*), intent(in):: String
    !--------------------------------------------------------------------------
    call ESMF_LogWrite(String, ESMF_LOGMSG_INFO)
    call ESMF_LogFlush()

  end subroutine log_write
    
end module ESMF_SWMF_Mod
!==============================================================================
