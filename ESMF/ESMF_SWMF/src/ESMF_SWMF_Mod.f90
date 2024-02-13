module ESMF_SWMF_Mod

  ! Various entities needed for the ESMF-SWMF coupling

  use ESMF
  implicit none

  ! named integer indexes for integer time arrays
  integer, parameter :: &
       Year_=1, Month_=2, Day_=3, Hour_=4, Minute_=5, Second_=6, MilliSec_=7

  ! number of MHD variables and their names
  integer, parameter :: nVar = 8
  character(len=3), dimension(nVar), parameter :: NameField_V = &
       [ 'Rho', 'Ux ', 'Uy ', 'Uz ', 'Bx ', 'By ', 'Bz ', 'P  ' ]

  ! Time related variables
  integer :: iStartTime_I(Year_:MilliSec_)  = & ! Start date-time
       [2000, 3, 21, 10, 45, 0, 0]            !   with defaults
  integer :: iFinishTime_I(Year_:MilliSec_) = & ! Finish date-time
       [2000, 3, 21, 10, 45, 0, 0]            !   with defaults
  real(ESMF_KIND_R8)      :: TimeSimulation = 0.0
  integer :: iCoupleFreq = 1                    ! Coupling frequency in seconds

  ! Variables related to the layout information
  ! SWMF runs on processor ranks iProcRootSwmf to iProcLastSwmf,
  ! ESMF runs on processor ranks iProcRootEsmf to iProcLastEsmf
  integer :: iProcRootSwmf, iProcLastSwmf
  integer :: iProcRootEsmf, iProcLastEsmf

  ! SWMF component to couple with
  character(len=2) :: NameSwmfComp = 'GM'

  ! The ESMF communicates with iProcCoupleSwmf within the SWMF layout
  ! This variable is determined from NameSwmfComp and the LAYOUT.in file.
  integer :: iProcCoupleSwmf=0

  ! When SWMF and ESMF are coupled, the one can block the whole SWMF
  ! or only the component the ESMF is communicating with. The latter
  ! is more efficient but it can result in a dead-lock if the ESMF and
  ! SWMF overlap.
  logical :: DoBlockAllSwmf=.false.

  ! Variables related to the grid used between the ESMF and SWMF components.
  ! Distance units are Earth radii, and the 2D uniform grid is defined
  ! in the Geo-Solar-Magnetic (GSM) coordinate system. It is orthogonal
  ! to the GSM-X axis.
  real(ESMF_KIND_R8), parameter :: rEarth = 6378000.0 ! Earth radius in SI unit
  integer            :: iMax=9, jMax=9                ! Default grid size
  real(ESMF_KIND_R8) :: yMin=-128.0                   ! Minimum Y in rEarth
  real(ESMF_KIND_R8) :: yMax=+128.0                   ! Maximum Y in rEarth
  real(ESMF_KIND_R8) :: zMin=-128.0                   ! Minimum Z in rEarth
  real(ESMF_KIND_R8) :: zMax=+128.0                   ! Maximum Z in rEarth

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
            'ESMF_ConfigLoadFile FAILED for file '//NameParamFile
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

    call ESMF_ConfigGetAttribute(Config, StringTmp, &
         label='Block all SWMF [y/n]:', rc=rc)
    if(rc == ESMF_SUCCESS) then
       DoBlockAllSwmf = StringTmp == 'y' .or. StringTmp == 'Y' .or. &
            StringTmp == 't' .or. StringTmp == 'T'
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
            'Setting default for SWMF Component: GM'
       NameSwmfComp = 'GM'
    end if

    ! Get grid size
    iDefaultTmp = iMax
    call ESMF_ConfigGetAttribute(Config, iMax, label='iMax:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0) write(*,*) 'ESMF_SWMF did not read iMax, ', &
            'setting default value= ', iDefaultTmp
       iMax = iDefaultTmp
    end if
    if(iMax < 1)then
       if(iProc == 0) write(*,*) 'ESMF_SWMF ERROR: ', &
            'iMax =',iMax,' should be positive!'
       rc = ESMF_FAILURE; RETURN
    end if

    iDefaultTmp = jMax
    call ESMF_ConfigGetAttribute(Config, jMax, label='jMax:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0) write(*,*) 'ESMF_SWMF did not read jMax, ', &
            'setting default value= ', iDefaultTmp
       jMax = iDefaultTmp
    end if
    if(jMax < 1)then
       if(iProc == 0) write(*,*) 'ESMF_SWMF ERROR: ', &
            'jMax =',jMax,' should be positive!'
       rc = ESMF_FAILURE; RETURN
    end if

    ! Get grid coordinate range
    DefaultTmp = yMin
    call ESMF_ConfigGetAttribute(Config, yMin, label='yMin:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0) write(*,*) 'ESMF_SWMF did not read yMin, '// &
            'setting default value= ', DefaultTmp
       yMin = DefaultTmp
    end if
    if(yMin >= 0.0)then
       if(iProc == 0) write(*,*) 'ESMF_SWMF ERROR: ', &
            'yMin =',yMin,' should be negative!'
       rc = ESMF_FAILURE; RETURN
    end if

    DefaultTmp = yMax
    call ESMF_ConfigGetAttribute(Config, yMax, label='yMax:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0) write(*,*) 'ESMF_SWMF did not read yMax, '// &
            'setting default value= ', DefaultTmp
       yMax = DefaultTmp
    end if
    if(yMax <= 0.0)then
       if(iProc == 0) write(*,*) 'ESMF_SWMF ERROR: ', &
            'yMax =',yMax,' should be positive!'
       rc = ESMF_FAILURE; RETURN
    end if

    DefaultTmp = zMin
    call ESMF_ConfigGetAttribute(Config, zMin, label='zMin:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0) write(*,*) 'ESMF_SWMF did not read zMin, '// &
            'setting default value= ', DefaultTmp
       zMin = DefaultTmp
    end if
    if(zMin >= 0.0)then
       if(iProc == 0) write(*,*) 'ESMF_SWMF ERROR: ', &
            'zMin =',zMin,' should be negative!'
       rc = ESMF_FAILURE; RETURN
    end if

    DefaultTmp = zMax
    call ESMF_ConfigGetAttribute(Config, zMax, label='zMax:', rc=rc)
    if(rc /= ESMF_SUCCESS) then
       if(iProc == 0) write(*,*) 'ESMF_SWMF did not read zMax, '// &
            'setting default value= ', DefaultTmp
       zMax = DefaultTmp
    end if
    if(zMax <= 0.0)then
       if(iProc == 0) write(*,*) 'ESMF_SWMF ERROR: ', &
            'zMax =',zMax,' should be positive!'
       rc = ESMF_FAILURE; RETURN
    end if

    ! Convert to SI units
    yMin = yMin*rEarth
    yMax = yMax*rEarth
    zMin = zMin*rEarth
    zMax = zMax*rEarth

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

    integer, intent(in)  :: iProc  ! rank of processor
    integer, intent(out) :: rc     ! error code

    integer :: iUnit
    character(len=100) :: String
    logical :: DoRead
    ! Get the root processor for the SWMF component to be coupled with
    ! from the LAYOUT.in file
    !--------------------------------------------------------------------------
    call CON_io_unit_new_ext(iUnit)
    open(unit=iUnit, file='LAYOUT.in', iostat=rc)
    if(rc /= 0)then
       if(iProc==0)write(*,*)'ESMF_SWMF ERROR: '// &
            'could not open LAYOUT.in file'
       rc = ESMF_FAILURE; RETURN
    end if

    ! Read the LAYOUT.in file
    DoRead = .false.
    READLAYOUT: do
       read(iUnit,'(a)',iostat=rc)String
       write(*,*)'!!! rc, String=', rc, String
    
       if(rc/=0)then
          if(iProc==0)write(*,*)'ESMF_SWMF ERROR: '// &
               'could not read LAYOUT.in file'
          rc = ESMF_FAILURE; RETURN
       end if
       if(String(1:13) == '#COMPONENTMAP') DoRead = .true.
       if(.not.DoRead) CYCLE
       if(String(1:2) == NameSwmfComp)then
          read(String, *, iostat=rc)NameSwmfComp, iProcCoupleSwmf
          write(*,*)'!!! rc, iProcCoupleSwmf=', rc, iProcCoupleSwmf
          if(rc/=0)then
             if(iProc==0)write(*,*)'ESMF_SWMF ERROR: '// &
                  'could not read iProcRoot for '//NameSwmfComp// &
                  ' from line=',trim(String),' in LAYOUT.in'
             rc = ESMF_FAILURE; RETURN
          end if
          ! Limit iProcCoupleSwmf by the number of SWMF processors
          iProcCoupleSwmf = min(iProcCoupleSwmf, iProcLastSwmf - iProcRootSwmf)
          EXIT READLAYOUT
       end if
       if(String == '#END')then
          if(iProc==0)write(*,*)'ESMF_SWMF ERROR: '// &
               'could not find component '//NameSwmfComp//' in LAYOUT.in file'
          rc = ESMF_FAILURE; RETURN
       end if
    end do READLAYOUT
    close(iUnit)

    if(iProc==0)write(*,*)'ESMF_SWMF: '// &
         NameSwmfComp//'   Root PE rank in SWMF=',iProcCoupleSwmf

  end subroutine read_swmf_layout
  !============================================================================
  subroutine add_mhd_fields(GridComp, State, InitialValue, rc)

    type(ESMF_GridComp), intent(inout) :: GridComp
    type(ESMF_State),    intent(inout) :: State
    real,                intent(in)    :: InitialValue
    integer,             intent(out)   :: rc

    type(ESMF_VM)        :: Vm
    type(ESMF_DELayout)  :: Layout
    type(ESMF_Grid)      :: Grid
    type(ESMF_ArraySpec) :: ArraySpec

    type(ESMF_Field)     :: MhdField
    integer              :: iVar
    ! Access to the MHD data
    real(ESMF_KIND_R8), pointer :: MyData(:,:)

    ! Name of the component
    character(len=100) :: Name

    ! Number of PE-s in SWMF
    integer              :: nProc

    ! Number of grid cells per DE=PE in SWMF
    integer, allocatable :: nCell_P(:)
    !--------------------------------------------------------------------------
    call ESMF_LogWrite("ESMF_SWMF_Mod:add_mhd_fields called", ESMF_LOGMSG_INFO)
    rc = ESMF_FAILURE

    ! Get default layout for the VM of GridComp
    call ESMF_GridCompGet(GridComp, vm=vm, grid=Grid, name=Name, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    Layout = ESMF_DELayoutCreate(vm, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    ! call ESMF_DELayoutPrint(Layout, rc=rc)
    ! if(rc /= ESMF_SUCCESS) return

    ! Distribute grid over the local vm
    ! This should be done in GridCreate
    !if(index(Name, 'SWMF') > 0)then
    !   ! For the SWMF all the data is on iProcCoupleSwmf
    !   call ESMF_VMGet(VM, petcount=nProc, rc=rc)
    !   if(rc /= ESMF_SUCCESS) RETURN
    !
    !   ! Create cell count array with all iMax cells on PE iProcCoupleSwmf
    !   allocate(nCell_P(nProc))
    !   nCell_P = 0
    !   nCell_P(iProcCoupleSwmf+1) = iMax
    !
    !   ! write(*,*)'!!! nProc, iProcCoupleSwmf, nCell_P=', &
    !   !     nProc, iProcCoupleSwmf, nCell_P
    !
    !   ! The SMWF grid is on a single processor
    !   call ESMF_GridDistribute(Grid, delayout=Layout, &
    !        countsPerDEDim1=nCell_P, &
    !        countsPerDEDim2=[ jMax ], rc=rc)
    !
    !   deallocate(nCell_P)
    !else
    !   ! For the ESMF the data is distributed
    !   call ESMF_GridDistribute(Grid, delayout=Layout, rc=rc)
    !end if
    !if(rc /= ESMF_SUCCESS) RETURN

    call ESMF_GridCompSet(GridComp, grid=Grid, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    ! Add MHD fields using the grid layout from DELayout
    call ESMF_ArraySpecSet(ArraySpec, rank=2, &
         typekind=ESMF_TYPEKIND_R8, rc=rc)
    if(rc /= ESMF_SUCCESS) RETURN

    do iVar=1, nVar
       MhdField = ESMF_FieldCreate(Grid, ArraySpec, &
            staggerloc=ESMF_STAGGERLOC_CENTER, &
            name=NameField_V(iVar), rc=rc)
       if(rc /= ESMF_SUCCESS) RETURN

       ! call ESMF_FieldGetDataPointer(MhdField, MyData)
       if(rc /= ESMF_SUCCESS) RETURN

       if(size(MyData) > 0) MyData = InitialValue

       !call ESMF_StateAddField(State, MhdField, rc=rc)
       !if(rc /= ESMF_SUCCESS) RETURN

       ! write(*,*)'iVar,shape(MyData),value=',iVar,shape(MyData),InitialValue
    end do

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ESMF_SWMF_Mod:add_mhd_fields returned", &
         ESMF_LOGMSG_INFO)

  end subroutine add_mhd_fields
  !============================================================================
end module ESMF_SWMF_Mod
!==============================================================================
