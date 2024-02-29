module ESMF_SWMF_Mod

  ! Various entities needed for the ESMF-SWMF coupling

  use ESMF
  use IE_ModSize, ONLY: nLatSwmf => IONO_nTheta, nLonSwmf => IONO_nPsi
  implicit none

  private
  public:: read_esmf_swmf_input, add_fields

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

  ! Variables related to the layout information
  ! SWMF runs on processor ranks iProcRootSwmf to iProcLastSwmf,
  ! ESMF runs on processor ranks iProcRootEsmf to iProcLastEsmf
  integer, public :: iProcRootSwmf, iProcLastSwmf, nProcSwmf
  integer, public :: iProcRootEsmf, iProcLastEsmf

  ! SWMF component to couple with
  character(len=2), public :: NameSwmfComp = 'IE'

  ! The ESMF communicates with SwmfComp within the SWMF.
  ! The processors used by SwmfComp are obtained from the PARAM.in file.
  integer, public:: iProc0SwmfComp=0, iProcLastSwmfComp=0, nProcSwmfComp=1

  ! If DoRunSwmf is true, run the SWMF for real, otherwise just pretend
  logical, public:: DoRunSwmf = .true.

  ! When SWMF and ESMF are coupled, one can block the whole SWMF
  ! or only the component the ESMF is communicating with. The latter
  ! is more efficient but it can result in a dead-lock if the ESMF and
  ! SWMF overlap.
  logical, public:: DoBlockAllSwmf=.false.

  ! Variables related to the grid used between the ESMF and SWMF components.
  ! This is a 2D spherical grid representing the height integrated ionosphere.
  ! In RIM it is in SM coordinates:
  ! +Z points to north magnetic dipole and the Sun is in the +X-Z halfplane.
  public:: nLonSwmf, nLatSwmf                 ! SWMF grid size from IE_ModSize
  integer, public:: nLonEsmf=48, nLatEsmf=120 ! Default ESMF grid size
  real(ESMF_KIND_R8), parameter, public:: LonMin = 0.0      ! Minimum longitude
  real(ESMF_KIND_R8), parameter, public:: LonMax = 360.0    ! Maximum longitude
  real(ESMF_KIND_R8), parameter, public:: LatMinSwmf = -90.0 ! Min SWMF lat
  real(ESMF_KIND_R8), parameter, public:: LatMaxSwmf = +90.0 ! Max SWMF lat
  real(ESMF_KIND_R8), public:: LatMinEsmf = -90.0 ! Min latitude of ESMF grid
  real(ESMF_KIND_R8), public:: LatMaxEsmf = +90.0 ! Max latitude of ESMF grid

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
    call ESMF_LogWrite("ESMF_SWMF_Mod:read_esmf_swmf_input called", &
         ESMF_LOGMSG_INFO)
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

    ! Read root PE for the SWMF
    call ESMF_ConfigGetAttribute(Config, iProcRootSwmf, &
         label='SWMF Root PE:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF: ', &
            'Settind default for SWMF Root PE: 0'
       iProcRootSwmf = 0
    end if
    if(iProcRootSwmf < 0) then
       if(iProc == 0)write(*,*) 'WARNING in ESMF_SWMF: ', &
            'SWMF Root PE rank negative! Setting it to 0'
       iProcRootSwmf = 0
    end if
    if(iProcRootSwmf >= nProc) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF: ', &
            'SWMF Root PE rank too large, setting nProc-1=',nProc-1
       iProcRootSwmf = nProc-1
    end if

    ! Read last PE for the SWMF
    call ESMF_ConfigGetAttribute(Config, iProcLastSwmf, &
         label='SWMF Last PE:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF: ', &
            'Setting default for SWMF Last PE: nProc-1=',nProc-1
       iProcLastSwmf = nProc-1
    end if
    if(iProcLastSwmf < iProcRootSwmf) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF: ', &
            'SWMF Last PE rank too small, setting it to iProcRootSwmf=',&
            iProcRootSwmf
       iProcLastSwmf = iProcRootSwmf
    end if
    if(iProcLastSwmf >= nProc) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF: ', &
            'SWMF Last PE rank too large, setting it to nProc-1=', nProc-1
       iProcLastSwmf = nProc-1
    end if
    ! Number of PEs used by the SWMF
    nProcSwmf = iProcLastSwmf - iProcRootSwmf + 1

    ! Read root PE for the ESMF
    call ESMF_ConfigGetAttribute(Config, iProcRootEsmf, &
         label='ESMF Root PE:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF: ', &
            'Setting default for ESMF Root PE: 0'
       iProcRootEsmf = 0
    end if
    if(iProcRootEsmf < 0) then
       if(iProc == 0)write(*,*) 'WARNING in ESMF_SWMF: ', &
            'Negative rank for ESMF Root PE! Setting it to 0'
       iProcRootEsmf = 0
    end if
    if(iProcRootEsmf >= nProc) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF: ', &
            'ESMF Root PE rank too large, setting it to nProc-1=', nProc-1
       iProcRootEsmf = nProc-1
    end if

    ! Read last PE for the ESMF
    call ESMF_ConfigGetAttribute(Config, iProcLastEsmf, &
         label='ESMF Last PE:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF: ', &
            'Setting default for ESMF Last PE: nProc-1= ',nProc-1
       iProcLastEsmf = nProc-1
    end if
    if(iProcLastEsmf < iProcRootEsmf) then
       if(iProc == 0)write(*,*) 'WARNING in ESMF_SWMF: ', &
            'ESMF Last PE rank too small! Setting it to iProcRootEsmf=', &
            iProcRootEsmf
       iProcLastEsmf = iProcRootEsmf
    end if
    if(iProcLastEsmf >= nProc) then
       if(iProc == 0)write(*,*) 'ESMF_SWMF: ', &
            'ESMF Last PE rank too large, setting it to nProc-1=', nProc-1
       iProcLastEsmf = nProc-1
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

    ! Get grid size
    iDefaultTmp = nLonEsmf
    call ESMF_ConfigGetAttribute(Config, nLonEsmf, label='nLon:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0) write(*,*) 'ESMF_SWMF did not read nLon, ', &
            'setting default value= ', iDefaultTmp
       nLonEsmf = iDefaultTmp
    end if
    if(nLonEsmf < 1)then
       if(iProc == 0) write(*,*) 'ESMF_SWMF ERROR: ', &
            'nLon =',nLonEsmf,' should be positive!'
       rc = ESMF_FAILURE; RETURN
    end if

    iDefaultTmp = nLatEsmf
    call ESMF_ConfigGetAttribute(Config, nLatEsmf, label='nLat:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0) write(*,*) 'ESMF_SWMF did not read nLat, ', &
            'setting default value= ', iDefaultTmp
       nLatEsmf = iDefaultTmp
    end if
    if(nLatEsmf < 1)then
       if(iProc == 0) write(*,*) 'ESMF_SWMF ERROR: ', &
            'nLat =',nLatEsmf,' should be positive!'
       rc = ESMF_FAILURE; RETURN
    end if

    DefaultTmp = LatMinEsmf
    call ESMF_ConfigGetAttribute(Config, LatMinEsmf, label='LatMin:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0) write(*,*) 'ESMF_SWMF did not read LatMin, '// &
            'setting default value= ', DefaultTmp
       LatMinEsmf = DefaultTmp
    end if

    DefaultTmp = LatMaxEsmf
    call ESMF_ConfigGetAttribute(Config, LatMaxEsmf, label='LatMax:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0) write(*,*) 'ESMF_SWMF did not read LatMax, '// &
            'setting default value= ', DefaultTmp
       LatMaxEsmf = DefaultTmp
    end if
    if(LatMaxEsmf <= LatMinEsmf)then
       if(iProc == 0) write(*,*) 'ESMF_SWMF ERROR: LatMaxEsmf=', &
            LatMaxEsmf, ' should be larger than LatMinEsmf=', LatMinEsmf
       rc = ESMF_FAILURE; RETURN
    end if

    call ESMF_ConfigDestroy(Config, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    call read_swmf_layout(iProc, rc)
    if(rc /= ESMF_SUCCESS) RETURN

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ESMF_SWMF_Mod:read_esmf_swmf_input returned", &
         ESMF_LOGMSG_INFO)

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
  subroutine add_fields(GridComp, State, IsFromEsmf, rc)

    type(ESMF_GridComp), intent(inout) :: GridComp
    type(ESMF_State),    intent(inout) :: State
    logical,             intent(in)    :: IsFromEsmf
    integer,             intent(out)   :: rc

    type(ESMF_Grid)      :: Grid
    type(ESMF_Field)     :: Field
    type(ESMF_ArraySpec) :: ArraySpec
    integer              :: iVar
    character(len=4)     :: NameField
    ! real(ESMF_KIND_R8), pointer :: Ptr_C(:,:)

    ! Name of the component
    character(len=100) :: Name
    !--------------------------------------------------------------------------
    call ESMF_LogWrite("ESMF_SWMF_Mod:add_fields called", ESMF_LOGMSG_INFO)
    rc = ESMF_FAILURE

    ! Get default layout for the VM of GridComp
    call ESMF_GridCompGet(GridComp, grid=Grid, name=Name, rc=rc)
    if(rc /= ESMF_SUCCESS) call my_error('ESMF_GridCompGet failed')

    call ESMF_ArraySpecSet(ArraySpec, rank=2, typekind=ESMF_TYPEKIND_R8)
    if (rc /= ESMF_SUCCESS) call my_error('ESMF_ArraySpecSet')
    
    if(IsFromEsmf)then
       do iVar = 1, nVarEsmf
          NameField = NameFieldEsmf_V(iVar)
          write(*,*)'Adding ESMF field=', NameField,' to ',trim(Name)
          Field = ESMF_FieldCreate(Grid, arrayspec=ArraySpec, &
               staggerloc=ESMF_STAGGERLOC_CENTER, name=NameField, rc=rc)
          if(rc /= ESMF_SUCCESS) call my_error('ESMF_FieldCreate ' &
               //trim(NameField)//' for '//trim(Name))
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
               staggerloc=ESMF_STAGGERLOC_CENTER, name=NameField, rc=rc)
          if(rc /= ESMF_SUCCESS) call my_error('ESMF_FieldCreate ' &
               //trim(NameField)//' for '//trim(Name))
          call ESMF_StateAdd(State, [Field], rc=rc)
          if(rc /= ESMF_SUCCESS) call my_error('ESMF_StateAdd ' &
               //trim(NameField)//' to '//trim(Name))
       end do
    endif
    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ESMF_SWMF_Mod:add_fields returned", &
         ESMF_LOGMSG_INFO)

  end subroutine add_fields
  !============================================================================
  subroutine my_error(String)

    character(len=*), intent(in) :: String

    write(*,*)'ERROR in ESMF_SWMF_Mod:', String
    call ESMF_Finalize
    stop

  end subroutine my_error
  !============================================================================

end module ESMF_SWMF_Mod
!==============================================================================
