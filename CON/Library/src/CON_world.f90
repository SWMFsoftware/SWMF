!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module  CON_world

  ! CON_world reads, stores and provides information about all the
  ! registered components. Components are registered by being listed
  ! in the #COMPONENTMAP/#LAYOUT command in the PARAM.in file,
  ! which is read by the world_setup subroutine.
  ! The information is stored in the private CompInfo_C array
  ! of derived type CompInfoType which is defined in CON_comp_info.
  ! The array is indexed from 1 to the number of registered components,
  ! which can be obtained with the n_comp public function.
  !
  ! The variable holding the index for the registered components should be
  ! named lComp where l stands for "listed".
  ! The fixed component ID, on the other hand, is referred to as iComp.
  ! Conversion from lComp to iComp can be done with the i_comp public
  ! function. A typical loop accessing all registered and used components
  ! is like this
  !
  ! do lComp = 1,n_comp()
  !     iComp = i_comp(lComp)
  !     if(.not.use_comp(iComp)) CYCLE
  !     if(is_proc0(iComp)) write(*,*) &
  !          'Number of processors for component ',lComp,&
  !          ' is ',n_proc(iComp)
  ! end do
  !
  ! The above example shows three access functions, use_comp,
  ! is_proc0 and n_proc which are called with the fixed component
  ! ID argument.
  ! All the similar functions, such as i_comm, i_group,
  ! i_proc, i_proc0 and {is_proc can be called with the
  ! fixed integer component ID, the two-character component name,
  ! or no argument. The no argument call returns the value relevant for
  ! the whole framework. Examples:
  !
  ! write(*,*) n_proc(IE_)   ! number of processors used by the IE component
  ! write(*,*) n_proc('GM')  ! number of processors used by the GM component
  ! write(*,*) n_proc(CON_)  ! number of processors used by active components
  ! write(*,*) n_proc('CON') ! number of processors used by active components
  ! write(*,*) n_proc()      ! number of processors used by SWMF.
  !
  ! While functions provide a convenient access to most of the information
  ! about the components, the general get_comp_info and
  ! put_comp_info subroutines provide access to any subset of
  ! the available information with a single call.

  use ModMpi
  use ModIoUnit, ONLY : UNITTMP_, io_unit_clean
  use ModUtilities, ONLY: CON_stop
  use CON_comp_param
  use CON_comp_info, &
       comp_info_init  => init, &
       comp_info_get   => get,  &
       comp_info_put   => put,  &
       comp_info_clean => clean

  implicit none

  private ! except

  public :: MaxComp         ! max number of components from CON_comp_param
  public :: world_init      ! constructor initializes registry for components
  public :: world_setup     ! set registry values (reads from PARAM.in)
  public :: world_clean     ! destructor cleans up
  public :: world_used      ! setup group for PEs used by active components
  public :: n_comp          ! return number of components
  public :: use_comp        ! return true if component is used
  public :: i_comp          ! return component ID for list index
  public :: get_comp_info   ! get info about a component
  public :: put_comp_info   ! put info about a component not set by setup
  public :: is_proc         ! return true if PE is used by component
  public :: is_proc0        ! retrun true if PE is the root for component
  public :: i_proc          ! return processor index in the component group
  public :: n_proc          ! return number of processors used by component
  public :: i_proc0         ! return world processor index of component's root
  public :: is_thread       ! return true if PE is used by a component thread
  public :: i_thread        ! return thread index for a component
  public :: n_thread        ! return number of threads used by a component
  public :: i_comm          ! return communicator used by component
  public :: i_group         ! return group used by component
  public :: i_proc_stride   ! return the PE stride
  public :: i_proc_last     ! return the last PE rank used by comp

  public :: check_i_comp    ! stop with error if component ID is out of range
  !                           this method is inherited from CON_comp_param
  integer, public:: nComp=0 ! Number of registered components
  integer, public:: &
       lComp_I(MaxComp) = -1! Convert named index to list index

  integer, parameter, public:: &
       CON_ = 0           ! Index of MPI group of PEs used by active components

  ! local variables

  integer :: iProcWorld   = -1 ! global PE rank
  integer :: nProcWorld   =  0 ! global number of PEs
  integer :: iCommWorld   = -1 ! global communicator
  integer :: iGroupWorld  = -1 ! global group
  integer :: iProc0World  =  0 ! root processor

  integer :: iProcUsed  = -1 ! PE rank within the used PE group
  integer :: nProcUsed  = -1 ! number of active PEs
  integer :: iCommUsed  = -1 ! communicator for active PEs
  integer :: iGroupUsed = -1 ! group for active PEs
  integer :: iProc0Used = -1 ! global index of root of active PEs

  integer :: MaxThread    =  1 ! maximum number of OpenMP threads

  integer, save :: iComp_C(MaxComp) ! Convert list index to named index

  type(CompInfoType), save, pointer :: CompInfo_C(:) ! Component info array

  interface world_setup; module procedure &
       setup_from_file
  end interface world_setup
  interface use_comp; module procedure &
       use_comp_name, use_comp_id
  end interface use_comp
  interface l_comp; module procedure &
       l_comp_name, l_comp_id
  end interface l_comp
  interface i_comp; module procedure &
       i_comp_name, i_comp_id
  end interface i_comp
  interface get_comp_info; module procedure &
       get_comp_info_name, get_comp_info_id
  end interface get_comp_info
  interface put_comp_info; module procedure &
       put_comp_info_name, put_comp_info_id
  end interface put_comp_info
  interface is_proc; module procedure &
       is_proc_name, is_proc_id, is_proc_world
  end interface is_proc
  interface is_proc0; module procedure &
       is_proc0_name, is_proc0_id, is_proc0_world
  end interface is_proc0
  interface i_proc; module procedure &
       i_proc_name, i_proc_id, i_proc_world
  end interface i_proc
  interface i_proc0; module procedure &
       i_proc0_name, i_proc0_id, i_proc0_world
  end interface i_proc0
  interface n_proc; module procedure &
       n_proc_name, n_proc_id, n_proc_world
  end interface n_proc
  interface i_proc_stride
     module procedure i_proc_stride_world
     module procedure i_proc_stride_id
     module procedure i_proc_stride_name
  end interface i_proc_stride
  interface i_proc_last
     module procedure i_proc_last_world
     module procedure i_proc_last_name
     module procedure i_proc_last_id
  end interface i_proc_last
  interface is_thread; module procedure &
       is_thread_name, is_thread_id
  end interface is_thread
  interface i_thread; module procedure &
       i_thread_name, i_thread_id
  end interface i_thread
  interface n_thread; module procedure &
       n_thread_name, n_thread_id
  end interface n_thread
  interface i_comm; module procedure &
       i_comm_name, i_comm_id, i_comm_world
  end interface i_comm
  interface i_group; module procedure &
       i_group_name, i_group_id, i_group_world
  end interface i_group

  ! revision history:
  !
  !  June 2003 - O. Volberg <volov@umich.edu> - initial version
  !  this code was inspired by Multi Program-Comp Handshaking (MPH) Utility
  !  written by Yun He and Chris Ding, NERSC/LBL, January 2001, and by
  !  MCTWorld class in the Model Coupling Tooolkit (MCT) written by
  !  J.W. Larson and R. Jacob at the Argonne National Laboratory,
  !  release 1.0.6, 2002.
  !
  !  July 12 2003 - G. Toth <gtoth@umich.edu> - simplify and add useful methods
  !  2019 - H. Zhou and G. Toth added thread related methods

  character(len=*),parameter :: NameMod='CON_world'

contains
  !============================================================================

  subroutine world_init(iComm)

    ! Initalize and start up MPI for the application.
    ! If the optional argument is present, use it as the communicator
    ! for the SWMF, and do not call MPI_init.
    ! If the optional argument is not present, use MPI_COMM_WORLD
    ! and call MPI_init.
    use omp_lib

    integer, intent(in), optional :: iComm

    integer :: iError

    character(len=*), parameter:: NameSub = 'world_init'
    !--------------------------------------------------------------------------

    if(present(iComm))then
       iCommWorld = iComm
    else
       iCommWorld = MPI_COMM_WORLD
    end if
    call MPI_COMM_RANK (iCommWorld, iProcWorld,  iError)
    call MPI_COMM_SIZE (iCommWorld, nProcWorld,  iError)
    call MPI_COMM_GROUP(iCommWorld, iGroupWorld, iError)
    nComp = 0
    nullify(CompInfo_C)

    !$ MaxThread = omp_get_max_threads()

  end subroutine world_init
  !============================================================================
  subroutine world_clean
    integer :: iError

    character(len=*), parameter:: NameSub = 'world_clean'
    !--------------------------------------------------------------------------
    call io_unit_clean

  end subroutine world_clean
  !============================================================================
  logical function is_proc_world()
    !--------------------------------------------------------------------------
    is_proc_world = iProcWorld>=0
  end function is_proc_world
  !============================================================================
  logical function is_proc0_world()
    !--------------------------------------------------------------------------
    is_proc0_world = iProcWorld==0
  end function is_proc0_world
  !============================================================================
  integer function i_proc_world()
    !--------------------------------------------------------------------------
    i_proc_world = iProcWorld
  end function i_proc_world
  !============================================================================
  integer function i_proc0_world()
    !--------------------------------------------------------------------------
    i_proc0_world = iProc0World
  end function i_proc0_world
  !============================================================================
  integer function n_proc_world()
    !--------------------------------------------------------------------------
    n_proc_world = nProcWorld
  end function n_proc_world
  !============================================================================
  integer function i_comm_world()
    !--------------------------------------------------------------------------
    i_comm_world = iCommWorld
  end function i_comm_world
  !============================================================================
  integer function i_group_world()
    !--------------------------------------------------------------------------
    i_group_world = iGroupWorld
  end function i_group_world
  !============================================================================
  integer function i_proc_stride_world()
    !--------------------------------------------------------------------------
    i_proc_stride_world=1
  end function i_proc_stride_world
  !============================================================================
  integer function i_proc_last_world()
    !--------------------------------------------------------------------------
    i_proc_last_world=nProcWorld-1
  end function i_proc_last_world
  !============================================================================
  subroutine setup_from_file

    ! Reads PARAM.in and registers the listed components. The MPI
    ! parameters are initialized using the init subroutine from
    ! CON_comp_info. The format is the following:
    !
    ! Name      First   Last    Stride   nThread
    ! ============================================
    ! #COMPONENTMAP
    ! GM       0        9999    8        8   ! GM runs with 8 threads
    ! PW       0        9999   -1       -1   ! PW runs with MaxThread threads
    ! IE       -3      -2       1            ! IE runs on nProc-3:nProc-2
    ! IM       -1      -1       1            ! IM runs on nProc-1
    !
    ! The layout information starts with the #COMPONENTMAP command
    ! (or the alternative name #LAYOUT) and ends with an empty line.
    ! The first column contains the
    ! two-character component name, the 2nd, 3rd, and 4th columns
    ! are integers containing the ranks of the first and last processors
    ! used by the component, and the stride between the processor ranks.
    ! The optional 4th parameter defines the number of OpenMP threads per MPI
    ! process for models that can use OpenMP. The Stride and nThread are
    ! typically equal.
    ! If the first (or last) rank is negative, it is evaluated as nProc-VALUE.
    ! If the last rank exceeds the number of processors nProc used by SWMF,
    ! it is reduced to the rank of the last processor (nProc-1) and a warning
    ! message is printed to STDOUT.
    ! If the stride or the number of threads are negative, they are replaced
    ! by the maximum number of threads MaxThread divided by their absolute
    ! values. If the number of threads nThread exceeds MaxThread, it gets
    ! reduced to MaxThread.

    ! One line of input
    character (len=100) :: String

    ! String starting and ending comp layout
    character (len=*), parameter  :: StringStart  = "#LAYOUT"
    character (len=*), parameter  :: StringStart2 = "#COMPONENTMAP"
    character (len=*), parameter  :: StringEnd    = "#END"
    character (len=*), parameter  :: StringEnd2   = " "

    ! Current line number
    integer :: nLine

    ! Input fields in the line
    character(len=16) :: Name
    integer :: iProcZero, iProcLast, iProcStride, nThread

    ! Temporary storage for processor ranges
    integer :: iProcRange_IC(ProcZero_:ProcStride_,MaxComp)

    ! Temporary storage for number of threads
    integer :: nThread_C(MaxComp)

    ! Temporary storage for component names
    character (len=lNameComp) :: Name_C(MaxComp)

    logical :: IsExisting
    integer :: lComp, iComp, iError

    ! Processor 0 looks for StringStart in NameMapFile
    character(len=*), parameter:: NameSub = 'setup_from_file'
    !--------------------------------------------------------------------------
    if(iProcWorld==0)then
       inquire(FILE=NameMapFile, EXIST=IsExisting)
       if(.not.IsExisting) NameMapFile = "PARAM.in"

       open(UNITTMP_,file=NameMapFile,iostat=iError,status="old",action="read")
       if(iError/=0) then
          write(*,'(3a)') NameSub,' SWMF_ERROR: cannot open file ',NameMapFile
          call CON_stop('Please check the file permissions of ' &
               //NameMapFile)
       endif

       nLine = 0
       do
          nLine = nLine + 1
          read(UNITTMP_, '(a)',IOSTAT=iError) String
          if(iError/=0)then
             close (UNITTMP_)
             write(*,'(a)') NameSub//' SWMF_ERROR: could not find '// &
                  StringStart//' or '//StringStart2//' in file '//NameMapFile
             call CON_stop('Please edit file '//NameMapFile)
          end if
          if(  String(1:len(StringStart))  == StringStart .or. &
               String(1:len(StringStart2)) == StringStart2) EXIT
       end do
    end if

    RECLOOP: do
       ! Processor zero reads and broadcasts lines
       if(iProcWorld==0) then
          read(UNITTMP_,'(a)',IOSTAT=iError) String
          if (iError/=0) then
             write(*,'(a,i2,a)')  NameSub// &
                  ' SWMF_ERROR: cannot read line after',&
                  nLine,' lines from the file '//NameMapFile
             call CON_stop('Please edit '//NameMapFile)
          end if
          nLine = nLine + 1
       end if

       call MPI_bcast(String,len(String),MPI_CHARACTER,0,iCommWorld,iError)

       if(String == StringEnd .or. String == StringEnd2) then
          if(nComp == 0) call CON_stop(NameSub // &
               'SWMF_ERROR: no components specified in the file='//NameMapFile)
          EXIT RECLOOP
       endif

       ! Process line in parallel
       read(String,*,iostat=iError) Name, iProcZero, iProcLast, iProcStride, &
            nThread

       if(iError/=0)then
          nThread = 1

          read(String,*,iostat=iError) Name, iProcZero, iProcLast, iProcStride
          if(iError/=0 .and. iProcWorld==0) then
             write(*,'(a,i2,a)')  NameSub// &
                  ' SWMF_ERROR: cannot read component name and PE range '// &
                  'from "'&
                  //trim(String)//'" at line ',nLine, &
                  ' from the file '//NameMapFile
             call CON_stop('Please edit '//NameMapFile)
          end if

       end if

       Name = adjustl(Name)
       if(is_valid_comp_name(Name)) then
          iComp=i_comp_name(Name)
       else
          if(iProcWorld==0) then

             write(*,'(a,i2,a)')  NameSub// &
                  ' SWMF_ERROR: invalid component name '//trim(Name)// &
                  ' at line ',nLine,' in the file '//NameMapFile
             call CON_stop('Please edit '//NameMapFile)
          endif
       endif
       nComp = nComp + 1

       ! Check if compoenent has been registered already
       if(lComp_I(iComp)>0 .and. iProcWorld==0)then
          write(*,'(a,i2,a,i2,a,i2,a)') &
               NameSub//' SWMF_ERROR component '//trim(Name)//&
               ' has been already registered as component ',&
               lComp_I(iComp),' and now again as component ',&
               nComp,' at line ',nLine,' in file '//NameMapFile
          call CON_stop(NameSub//' SWMF_ERROR Please edit '//NameMapFile)
          CYCLE RECLOOP
       end if

       ! Correct and store parameters read from this line
       lComp_I(iComp)                   = nComp  ! list index for named index
       iComp_C(nComp)                   = iComp  ! named index for list index
       Name_C(nComp)                    = Name

       ! Negative iProcZero is interpreted as counting back from the end
       if(iProcZero < 0) iProcZero = max(0, nProcWorld + iProcZero)
       ! Negative iProcLast is interpreted as counting back from the end
       if(iProcLast < 0) iProcLast = max(0, nProcWorld + iProcLast)

       iProcRange_IC(ProcZero_,nComp) = min(iProcZero, iProcLast, nProcWorld-1)
       iProcRange_IC(ProcLast_,nComp) = min(iProcLast, nProcWorld-1)
       ! Negative iProcStride is interpreted as MaxThread divided by its abs
       ! value
       if(iProcStride < 0) iProcStride = MaxThread/abs(iProcStride)
       iProcRange_IC(ProcStride_,nComp) = max(iProcStride,1)

       ! Negative nThread is interpreted as MaxThread divided by its abs value
       if(nThread < 0) nThread = MaxThread/abs(nThread)
       ! Store thread info
       nThread_C(nComp) = max(min(MaxThread,nThread),1)

       if(iProcWorld==0)then
          ! Report adjustments
          if(iProcRange_IC(ProcZero_, nComp) /= iProcZero) &
               write(*,'(a,i5,a,i5)')NameSub//&
               ' SWMF_WARNING for '//trim(Name)// &
               ' iProcZero=',iProcZero,&
               ' changed to ',iProcRange_IC(ProcZero_,  nComp)
          if(iProcRange_IC(ProcLast_, nComp) /= iProcLast) &
               write(*,'(a,i5,a,i5)')NameSub//&
               ' SWMF_WARNING for '//trim(Name)// &
               ' iProcLast=',iProcLast,&
               ' changed to ',iProcRange_IC(ProcLast_,  nComp)
          if(iProcRange_IC(ProcStride_,nComp) /= iProcStride) &
               write(*,'(a,i5,a,i5)')NameSub//&
               ' SWMF_WARNING for '//trim(Name)// &
               ' iProcStride=',iProcStride,&
               ' changed to ',iProcRange_IC(ProcStride_,nComp)
       end if
    end do RECLOOP
    if(iProcWorld == 0) close (UNITTMP_)

    allocate(CompInfo_C(CON_:nComp))
    CompInfo_C(CON_) % Name = 'CON'  ! Group of used PEs

    ! Check component names
    do lComp = 1, nComp
       call comp_info_init(CompInfo_C(lComp),&
            Name_C(lComp), iGroupWorld, iCommWorld, iProcRange_IC(:,lComp),&
            nThread_C(lComp), iError)
       if(iError > 0)then
          if(iProcWorld==0)write(*,'(a,3i5)')NameSub//&
               ' SWMF_ERROR: cannot produce layout for component '// &
               Name_C(lComp)//' using processor range ',&
               iProcRange_IC(:,lComp)
          call CON_stop(&
               'Please edit layout or increase number of PEs')
       end if
    end do

  end subroutine setup_from_file
  !============================================================================
  subroutine world_used(IsVerbose)

    logical, intent(in), optional:: IsVerbose

    ! Create communicator for all PEs used by the active componenats
    ! at the beginning of a session.
    ! OpenMP or future sessions may use the rest.
    ! Write out basic information if IsVerbose is true

    integer:: lComp, iError, iGroup

    character(len=*), parameter:: NameSub = 'world_used'
    !--------------------------------------------------------------------------
    if(nProcUsed /= -1)then
       call MPI_group_free(iGroupUsed, iError)
       if(iProcUsed /= MPI_UNDEFINED) call MPI_comm_free(iCommUsed, iError)
    end if

    iProc0Used = nProcWorld      ! initial value for the root of used PEs
    iGroupUsed = MPI_GROUP_EMPTY ! initial value for the group of used PEs
    do lComp = 1, nComp
       if(.not. CompInfo_C(lComp) % Use) CYCLE
       iGroup = iGroupUsed
       ! Add new group. Try to make the root the PE with smallest global rank
       if(iProc0Used < CompInfo_C(lComp)%iMpiParam_I(ProcZero_))then
          call MPI_group_union(iGroup, CompInfo_C(lComp)%iMpiParam_I(Group_), &
               iGroupUsed, iError)
       else
          iProc0Used = CompInfo_C(lComp)%iMpiParam_I(ProcZero_)
          call MPI_group_union(CompInfo_C(lComp)%iMpiParam_I(Group_), iGroup, &
               iGroupUsed, iError)
       end if
       if(iGroup /= MPI_GROUP_EMPTY) call MPI_group_free(iGroup, iError)
    end do
    if(iProc0Used == nProcWorld) &
         call CON_stop(NameSub//': no used components!')

    ! Calculate rank and size for active PE group
    call MPI_group_size(iGroupUsed, nProcUsed, iError)
    call MPI_group_rank(iGroupUsed, iProcUsed, iError)

    ! Create communicator of used PEs
    call MPI_comm_create(iCommWorld, iGroupUsed, iCommUsed, iError)

    ! Store MPI info for CON_ group
    CompInfo_C(CON_) % iMpiParam_I = &
         [iProc0Used, -1, -1, iProcUsed, nProcUsed, iCommUsed, iGroupUsed]

    if(.not.present(IsVerbose)) RETURN

    if(.not.IsVerbose) RETURN

    if(iProcWorld==0)write(*,*) NameSub,' finished with iProc0Used=', &
         iProc0Used,' nProcUsed=', nProcUsed

  end subroutine world_used
  !============================================================================
  integer function n_comp()
    character(len=*), parameter:: NameSub = 'n_comp'
    !--------------------------------------------------------------------------
    n_comp =  nComp
    if(nComp==0) &
         write(*,'(a)') NameSub//' SWMF_ERROR: no components defined '
  end function n_comp
  !============================================================================
  integer function l_comp_id(iComp, NameSub, DoCheckRegistered)

    ! Convert named index into list index and check everything

    integer, intent(in)          :: iComp
    character(len=*), intent(in) :: NameSub ! Name of the caller
    logical, intent(in), optional:: DoCheckRegistered
    integer :: lComp
    logical :: DoCheck

    !--------------------------------------------------------------------------
    if(iComp == CON_)then
       l_comp_id = CON_
       RETURN
    end if

    call check_i_comp(iComp,NameSub)

    lComp = lComp_I(iComp)

    DoCheck = .true.; if(present(DoCheckRegistered)) DoCheck=DoCheckRegistered
    if(DoCheck)then
       if(lComp < 1 .or. lComp>nComp) then
          write(*,'(a,i3,a,i3)')NameSub//' SWMF_ERROR iComp '// &
               NameComp_I(iComp)//' is not found among the ',&
               nComp,' registered components, lComp=',lComp
          call CON_stop('Error in the caller method')
       end if
    end if

    l_comp_id = lComp

  end function l_comp_id
  !============================================================================
  integer function l_comp_name(NameComp, NameSub, DoCheckRegistered)

    ! Convert the component name into list index and check everything

    character(len=*), intent(in) :: NameComp ! Name of the component
    character(len=*), intent(in) :: NameSub  ! Name of the caller
    logical, intent(in), optional:: DoCheckRegistered
    integer :: lComp
    logical :: DoCheck

    !--------------------------------------------------------------------------
    if(NameComp == 'CON')then
       l_comp_name = CON_
       RETURN
    end if

    if(.not.is_valid_comp_name(NameComp))then
       write(*,'(a)')NameSub//' SWMF_ERROR name '//NameComp// &
            ' is not a valid component name'
       call CON_stop('Error in the caller method')
    end if

    lComp = lComp_I(i_comp(NameComp))

    DoCheck = .true.; if(present(DoCheckRegistered)) DoCheck=DoCheckRegistered
    if(DoCheck)then
       if(lComp < 1 .or. lComp>nComp) then
          write(*,'(a,i3,a,i3)')NameSub//' SWMF_ERROR component '// &
               NameComp//' is not found among the ',&
               nComp,' registered components, lComp=',lComp
          call CON_stop('Error in the caller method')
       end if
    end if

    l_comp_name = lComp

  end function l_comp_name
  !============================================================================
  integer function i_comp_id(lComp)
    ! return named index for list index lComp
    ! Note: i_comp_name is implemented in CON_comp_param

    integer, intent(in) :: lComp
    integer :: iComp
    character(len=*), parameter:: NameSub = 'i_comp_id'
    !--------------------------------------------------------------------------
    if(lComp<1.or.lComp>nComp)then
       write(*,'(a,i3,a,i3)')NameSub//' SWMF_ERROR lComp=',lComp,&
            ' is out of range 1..nComp=',nComp
       call CON_stop('Error in caller method')
    end if
    iComp=iComp_C(lComp)
    if(iComp<1.or.iComp>MaxComp)then
       write(*,'(a,i3,a,i3,a,i3)')NameSub//' SWMF_ERROR lComp=',lComp,&
            ' gives iComp=',iComp,&
            ' which is out of range 1..MaxComp=',MaxComp
       call CON_stop(NameSub//'iComp_C array is not set correctly')
    end if
    i_comp_id = iComp

  end function i_comp_id
  !============================================================================
  subroutine get_comp_info_name(NameComp, &
       iProc, nProc, iComm, iGroup, iProcZero, iProcLast, iProcStride, &
       nThread, iUnitOut, Use, Name, NameVersion, Version, CompInfo)

    character(len=*), intent(in) :: NameComp
    integer, optional, intent(out) :: &
         iProc, nProc, iComm, iGroup, iProcZero,  iProcLast, iProcStride, &
         nThread, iUnitOut
    logical,                     optional, intent(out) :: use
    character(len=lNameComp),    optional, intent(out) :: Name
    character(len=lNameVersion), optional, intent(out) :: NameVersion
    real,                        optional, intent(out) :: Version
    type(CompInfoType),          optional, intent(out) :: CompInfo
    character(len=*), parameter:: NameSub = 'get_comp_info_name'
    !--------------------------------------------------------------------------
    if(present(CompInfo)) CompInfo = CompInfo_C(l_comp(NameComp,NameSub))

    call comp_info_get(CompInfo_C(l_comp(NameComp,NameSub)), &
         iProc, nProc, iComm, iGroup, iProcZero, iProcLast, iProcStride, &
         nThread, iUnitOut, Use, Name, NameVersion, Version)

  end subroutine get_comp_info_name
  !============================================================================
  subroutine get_comp_info_id(iComp, &
       iProc, nProc, iComm, iGroup, iProcZero, iProcLast, iProcStride, &
       nThread, iUnitOut, Use, Name, NameVersion, Version, CompInfo)

    integer, intent(in) :: iComp

    integer, optional, intent(out) :: &
         iProc, nProc, iComm, iGroup, iProcZero,  iProcLast, iProcStride, &
         nThread, iUnitOut
    logical,                     optional, intent(out) :: use
    character(len=lNameComp),    optional, intent(out) :: Name
    character(len=lNameVersion), optional, intent(out) :: NameVersion
    real,                        optional, intent(out) :: Version
    type(CompInfoType),          optional, intent(out) :: CompInfo

    ! Obtain any information about the component which is identified
    ! with its ID iComp, or name (using get_comp_info_name).

    character(len=*), parameter:: NameSub = 'get_comp_info_id'
    !--------------------------------------------------------------------------
    if(present(CompInfo)) CompInfo = CompInfo_C(l_comp(iComp,NameSub))

    call comp_info_get(CompInfo_C(l_comp(iComp,NameSub)), &
         iProc, nProc, iComm, iGroup, iProcZero, iProcLast, iProcStride, &
         nThread, iUnitOut, Use, Name, NameVersion, Version)

  end subroutine get_comp_info_id
  !============================================================================
  subroutine put_comp_info_name(NameComp, iUnitOut, Use, NameVersion, Version)

    character(len=*), intent(in) :: NameComp
    integer,          optional, intent(in) :: iUnitOut    ! for STDOUT
    logical,          optional, intent(in) :: use         ! is used
    character(len=*), optional, intent(in) :: NameVersion ! version name
    real,             optional, intent(in) :: Version     ! version number

    ! Set component information which is not set by setup.
    ! The first argument is either the two-character component name
    ! NameComp or the integer component ID (using put_comp_info_id).

    character(len=*), parameter:: NameSub = 'put_comp_info_name'
    !--------------------------------------------------------------------------
    call comp_info_put(CompInfo_C(l_comp(NameComp,NameSub)), &
         iUnitOut, Use, NameVersion, Version)

  end subroutine put_comp_info_name
  !============================================================================
  subroutine put_comp_info_id(iComp, iUnitOut, Use, NameVersion, Version)

    integer, intent(in) :: iComp
    integer,          optional, intent(in) :: iUnitOut
    logical,          optional, intent(in) :: use
    character(len=*), optional, intent(in) :: NameVersion
    real,             optional, intent(in) :: Version
    character(len=*), parameter:: NameSub = 'put_comp_info_id'
    !--------------------------------------------------------------------------
    call comp_info_put(CompInfo_C(l_comp(iComp,NameSub)), &
         iUnitOut, Use, NameVersion, Version)

  end subroutine put_comp_info_id
  !============================================================================
  logical function use_comp_name(NameComp)

    ! Return true if component is used
    character(len=*), intent(in) :: NameComp
    integer :: lComp
    character(len=*), parameter:: NameSub = 'use_comp_name'
    !--------------------------------------------------------------------------
    lComp = l_comp(NameComp, NameSub, .false.)
    if(lComp > 0)then
       use_comp_name = CompInfo_C(lComp) % Use
    else
       use_comp_name = .false.
    end if

  end function use_comp_name
  !============================================================================
  logical function use_comp_id(iComp)

    ! Return true if component is used
    integer, intent(in) :: iComp
    integer :: lComp
    character(len=*), parameter:: NameSub = 'use_comp_id'
    !--------------------------------------------------------------------------
    lComp = l_comp(iComp, NameSub, .false.)
    if(lComp > 0)then
       use_comp_id = CompInfo_C(lComp) % Use
    else
       use_comp_id = .false.
    end if

  end function use_comp_id
  !============================================================================
  logical function is_proc_name(Name)
    character(len=*), intent(in) :: Name
    integer :: lComp
    character(len=*), parameter:: NameSub = 'is_proc_name'
    !--------------------------------------------------------------------------
    lComp = l_comp(Name, NameSub, .false.)
    if(lComp >= CON_)then
       is_proc_name = CompInfo_C(lComp) % iMpiParam_I(Proc_) >= 0
    else
       is_proc_name=.false.
    end if
  end function is_proc_name
  !============================================================================
  logical function is_proc_id(iComp)
    integer, intent(in) :: iComp
    integer :: lComp
    character(len=*), parameter:: NameSub = 'is_proc_id'
    !--------------------------------------------------------------------------
    lComp = l_comp(iComp,NameSub,.false.)
    if(lComp >= CON_)then
       is_proc_id = CompInfo_C(lComp) % iMpiParam_I(Proc_) >= 0
    else
       is_proc_id = .false.
    end if

  end function is_proc_id
  !============================================================================
  logical function is_proc0_name(Name)
    character(len=*), intent(in) :: Name

    character(len=*), parameter:: NameSub = 'is_proc0_name'
    !--------------------------------------------------------------------------
    is_proc0_name = CompInfo_C(l_comp(Name,NameSub)) % iMpiParam_I(Proc_)==0
  end function is_proc0_name
  !============================================================================
  logical function is_proc0_id(iComp)
    integer, intent(in) :: iComp

    character(len=*), parameter:: NameSub = 'is_proc0_id'
    !--------------------------------------------------------------------------
    is_proc0_id = CompInfo_C(l_comp(iComp,NameSub)) % iMpiParam_I(Proc_)==0
  end function is_proc0_id
  !============================================================================
  integer function i_proc_name(Name)
    character(len=*), intent(in) :: Name

    character(len=*), parameter:: NameSub = 'i_proc_name'
    !--------------------------------------------------------------------------
    i_proc_name = CompInfo_C(l_comp(Name,NameSub)) % iMpiParam_I(Proc_)
  end function i_proc_name
  !============================================================================
  integer function i_proc_id(iComp)
    integer, intent(in) :: iComp

    character(len=*), parameter:: NameSub = 'i_proc_id'
    !--------------------------------------------------------------------------
    i_proc_id = CompInfo_C(l_comp(iComp,NameSub)) % iMpiParam_I(Proc_)
  end function i_proc_id
  !============================================================================
  integer function n_proc_name(Name)
    character(len=*), intent(in) :: Name

    character(len=*), parameter:: NameSub = 'n_proc_name'
    !--------------------------------------------------------------------------
    n_proc_name = CompInfo_C(l_comp(Name,NameSub)) % iMpiParam_I(nProc_)
  end function n_proc_name
  !============================================================================
  integer function n_proc_id(iComp)
    integer, intent(in) :: iComp

    character(len=*), parameter:: NameSub = 'n_proc_id'
    !--------------------------------------------------------------------------
    n_proc_id = CompInfo_C(l_comp(iComp,NameSub)) % iMpiParam_I(nProc_)
  end function n_proc_id
  !============================================================================
  integer function i_proc0_name(Name)
    character(len=*), intent(in) :: Name

    character(len=*), parameter:: NameSub = 'i_proc0_name'
    !--------------------------------------------------------------------------
    i_proc0_name = CompInfo_C(l_comp(Name,NameSub)) % iMpiParam_I(ProcZero_)
  end function i_proc0_name
  !============================================================================
  integer function i_proc0_id(iComp)
    integer, intent(in) :: iComp

    character(len=*), parameter:: NameSub = 'i_proc0_id'
    !--------------------------------------------------------------------------
    i_proc0_id = CompInfo_C(l_comp(iComp,NameSub)) % iMpiParam_I(ProcZero_)
  end function i_proc0_id
  !============================================================================
  logical function is_thread_name(Name)
    character(len=*), intent(in) :: Name
    integer :: lComp
    character(len=*), parameter:: NameSub = 'is_thread_name'
    !--------------------------------------------------------------------------
    lComp = l_comp(Name, NameSub, .false.)
    if(lComp > 0)then
       is_thread_name = CompInfo_C(lComp) % iThread >= 0
    else
       is_thread_name = .false.
    end if
  end function is_thread_name
  !============================================================================
  logical function is_thread_id(iComp)
    integer, intent(in) :: iComp
    integer :: lComp
    character(len=*), parameter:: NameSub = 'is_thread_id'
    !--------------------------------------------------------------------------
    lComp = l_comp(iComp, NameSub, .false.)
    if(lComp > 0)then
       is_thread_id = CompInfo_C(lComp) % iThread >= 0
    else
       is_thread_id = .false.
    end if
  end function is_thread_id
  !============================================================================
  integer function i_thread_name(Name)
    character(len=*), intent(in) :: Name

    character(len=*), parameter:: NameSub = 'i_thread_name'
    !--------------------------------------------------------------------------
    i_thread_name =  CompInfo_C(l_comp(Name,NameSub)) % iThread
  end function i_thread_name
  !============================================================================
  integer function i_thread_id(iComp)
    integer, intent(in) :: iComp

    character(len=*), parameter:: NameSub = 'i_thread_id'
    !--------------------------------------------------------------------------
    i_thread_id =  CompInfo_C(l_comp(iComp,NameSub)) % iThread
  end function i_thread_id
  !============================================================================
  integer function n_thread_name(Name)
    character(len=*), intent(in) :: Name

    character(len=*), parameter:: NameSub = 'n_thread_name'
    !--------------------------------------------------------------------------
    n_thread_name =  CompInfo_C(l_comp(Name,NameSub)) % nThread
  end function n_thread_name
  !============================================================================
  integer function n_thread_id(iComp)
    integer, intent(in) :: iComp

    character(len=*), parameter:: NameSub = 'n_thread_id'
    !--------------------------------------------------------------------------
    n_thread_id =  CompInfo_C(l_comp(iComp,NameSub)) % nThread
  end function n_thread_id
  !============================================================================
  integer function i_comm_name(Name)
    character(len=*), intent(in) :: Name

    character(len=*), parameter:: NameSub = 'i_comm_name'
    !--------------------------------------------------------------------------
    i_comm_name =  CompInfo_C(l_comp(Name,NameSub)) % iMpiParam_I(Comm_)
  end function i_comm_name
  !============================================================================
  integer function i_comm_id(iComp)
    integer, intent(in) :: iComp

    character(len=*), parameter:: NameSub = 'i_comm_id'
    !--------------------------------------------------------------------------
    i_comm_id =  CompInfo_C(l_comp(iComp,NameSub)) % iMpiParam_I(Comm_)
  end function i_comm_id
  !============================================================================
  integer function i_group_name(Name)
    character(len=*), intent(in) :: Name

    character(len=*), parameter:: NameSub = 'i_group_name'
    !--------------------------------------------------------------------------
    i_group_name =  CompInfo_C(l_comp(Name,NameSub)) % iMpiParam_I(Group_)
  end function i_group_name
  !============================================================================
  integer function i_group_id(iComp)
    integer, intent(in) :: iComp

    character(len=*), parameter:: NameSub = 'i_group_id'
    !--------------------------------------------------------------------------
    i_group_id =  CompInfo_C(l_comp(iComp,NameSub)) % iMpiParam_I(Group_)
  end function i_group_id
  !============================================================================
  integer function i_proc_stride_name(Name)
    character(len=*), intent(in) :: Name

    character(len=*), parameter:: NameSub = 'i_proc_stride_name'
    !--------------------------------------------------------------------------
    i_proc_stride_name =  &
         CompInfo_C(l_comp(Name,NameSub)) % iMpiParam_I(ProcStride_)
  end function i_proc_stride_name
  !============================================================================
  integer function i_proc_stride_id(iComp)
    integer, intent(in) :: iComp

    character(len=*), parameter:: NameSub = 'i_proc_stride_id'
    !--------------------------------------------------------------------------
    i_proc_stride_id =  &
         CompInfo_C(l_comp(iComp,NameSub)) % iMpiParam_I(ProcStride_)
  end function i_proc_stride_id
  !============================================================================
  integer function i_proc_last_name(Name)
    character(len=*),intent(in)::Name
    character(len=*), parameter:: NameSub = 'i_proc_last_name'
    !--------------------------------------------------------------------------
    i_proc_last_name = &
         CompInfo_C(l_comp(Name,NameSub)) % iMpiParam_I(ProcLast_)
  end function i_proc_last_name
  !============================================================================
  integer function i_proc_last_id(iComp)
    integer,intent(in)::iComp
    character(len=*), parameter:: NameSub = 'i_proc_last_id'
    !--------------------------------------------------------------------------
    i_proc_last_id=CompInfo_C(l_comp(iComp,NameSub)) % iMpiParam_I(ProcLast_)
  end function i_proc_last_id
  !============================================================================

end module CON_world
!==============================================================================
