!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!               Space Weather Modeling Framework (SWMF)                !
!    Center for Space Environment Modeling, The University of Michigan !
!-----------------------------------------------------------------------
!
!BOP -------------------------------------------------------------------------
!
!MODULE: CON_world - component registration and handling component info
! 
!DESCRIPTION: 
! CON\_world reads, stores and provides information about all the
! registered components. Components are registered by being listed
! in the LAYOUT.in file wich is read by the {\bf setup} subroutine.
! The information is stored in the private {\bf CompInfo\_C} array
! of derived type CompInfoType which is defined in CON\_comp\_info.
! The array is indexed from 1 to the number of registered components
! which can be obtained with the {\bf n\_comp} public function.
!
! The variable holding the index for the registered components should be 
! named {\bf lComp} where {\bf l} stands for {\it listed}.
! The fixed component ID, on the other hand, is referred to as {\bf iComp}.
! Conversion from {\bf lComp} to {\bf iComp} can be done with the
! {\bf i\_comp} public function. A typical loop accessing all 
! registered and used components is like this
! \begin{verbatim}
! do lComp = 1,n_comp()
!     iComp = i_comp(lComp)
!     if(.not.use_comp(iComp)) CYCLE
!     if(is_proc0(iComp)) write(*,*) &
!          'Number of processors for component ',lComp,&
!          ' is ',n_proc(iComp)
! end do
! \end{verbatim}
! The above example shows three access functions, {\bf use\_comp},
! {\bf is\_proc0} and {\bf n\_proc} which are called with the fixed component 
! ID argument.
! All the similar functions, such as {\bf i\_comm, i\_group,
! i\_proc, i\_proc0} and {is\_proc} can be called with the
! fixed integer component ID, the two-character component name,
! or no argument. The no argument call returns the value relevant for
! the whole framework. Examples:
! \begin{verbatim}
!    write(*,*) n_proc(IE_)  ! number of processors used by the IE component
!    write(*,*) n_proc('GM') ! number of processors used by the GM component
!    write(*,*) n_proc()     ! number of processors used by SWMF.
! \end{verbatim}
! While functions provide a convenient access to most of the information
! about the components, the general {\bf get\_comp\_info} and 
! {\bf put\_comp\_info} subroutines provide access to any subset of
! the available information with a single call.

!INTERFACE:	
!
module  CON_world
  !
  !USES:
  !
  use ModMpi
  use ModIoUnit, ONLY : UNITTMP_, io_unit_clean
  use CON_comp_param
  use CON_comp_info, &
       comp_info_init  => init, &
       comp_info_get   => get,  &
       comp_info_put   => put,  &
       comp_info_clean => clean

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:
  public :: MaxComp
  public :: world_init      ! constructor initializes registry for components
  public :: world_setup     ! set registry values (reads from LAYOUT.in)
  public :: world_clean     ! destructor cleans up
  public :: world_abort     ! abort execution
  public :: n_comp          ! return number of components
  public :: use_comp        ! return true if component is used
  public :: i_comp          ! return component ID for list index
  public :: is_proc0_world  ! return true if processor is root of world group
  public :: i_proc_world    ! return world processor index
  public :: n_proc_world    ! return number of processors in world
  public :: i_comm_world    ! return world communicator
  public :: i_group_world   ! return world group
  public :: get_comp_info   ! get info about a component
  public :: put_comp_info   ! put info about a component not set by setup
  public :: is_proc         ! return true if PE is used by component
  public :: is_proc0        ! retrun true if PE is the root for component
  public :: i_proc          ! return processor index in the component group
  public :: n_proc          ! return number of processors used by component
  public :: i_proc0         ! return world processor index of component's root
  public :: i_comm          ! return communicator used by component
  public :: i_group         ! return group used by component
  public :: i_proc_stride   ! return the PE stride
  public :: i_proc_last     ! return the last PE rank used by comp

  public :: check_i_comp    ! stop with error if component ID is out of range
  !                           this method is inherited from CON_comp_param

  !LOCAL VARIABLES:
  integer :: iProcWorld=-1  ! global PE rank
  integer :: nProcWorld=0   ! global number of PE-s
  integer :: iCommWorld=-1  ! global communicator
  integer :: iGroupWorld=-1 ! global group
  integer :: iProc0World=0

  integer :: nComp=0        ! Number of registered components

  integer, save :: lComp_I(MaxComp) ! Convert named index to list index
  integer, save :: iComp_C(MaxComp) ! Convert list index to named index

  type(CompInfoType), save, pointer :: CompInfo_C(:) ! Component info array

  interface world_setup; module procedure &
       setup_from_file  
  end interface
  interface use_comp; module procedure &
       use_comp_name, use_comp_id
  end interface
  interface l_comp; module procedure &
       l_comp_name, l_comp_id
  end interface
  interface i_comp; module procedure &
       i_comp_name, i_comp_id
  end interface
  interface get_comp_info; module procedure &
       get_comp_info_name, get_comp_info_id
  end interface
  interface put_comp_info; module procedure &
       put_comp_info_name, put_comp_info_id
  end interface
  interface is_proc; module procedure &
       is_proc_name, is_proc_id, is_proc_world
  end interface
  interface is_proc0; module procedure &
       is_proc0_name, is_proc0_id, is_proc0_world
  end interface
  interface i_proc; module procedure &
       i_proc_name, i_proc_id, i_proc_world
  end interface
  interface i_proc0; module procedure &
       i_proc0_name, i_proc0_id, i_proc0_world
  end interface
  interface n_proc; module procedure &
       n_proc_name, n_proc_id, n_proc_world
  end interface
  interface i_comm; module procedure &
       i_comm_name, i_comm_id, i_comm_world
  end interface
  interface i_group; module procedure &
       i_group_name, i_group_id, i_group_world
  end interface
  interface i_proc_stride
     module procedure i_proc_stride_world
     module procedure i_proc_stride_id
     module procedure i_proc_stride_name
  end interface
  interface i_proc_last
     module procedure i_proc_last_world
     module procedure i_proc_last_name
     module procedure i_proc_last_id
  end interface

  !REVISION HISTORY: 
  !
  !  June 2003 - O. Volberg <volov@umich.edu> - initial version
  !  this code was inspired by Multi Program-Comp Handshaking (MPH) Utility
  !  written by Yun He and Chris Ding, NERSC/LBL, January 2001, and by 
  !  MCTWorld class in the Model Coupling Tooolkit (MCT) written by 
  !  J.W. Larson and R. Jacob at the Argonne National Laboratory, 
  !  release 1.0.6, 2002.
  !
  !  July 12 2003 - G. Toth <gtoth@umich.edu> - simplify and add useful methods
  !EOP ------------------------------------------------------------------------

  character(len=*),parameter :: NameMod='CON_world'

contains

  subroutine world_init(iComm)
    !\  
    ! Initalize and start up MPI for the application.
    ! If the optional argument is present, use it as the communicator
    ! for the SWMF, and do not call MPI_init.
    ! If the optional argument is not present, use MPI_COMM_WORLD
    ! and call MPI_init.
    !/
    integer, intent(in), optional :: iComm 

    integer :: iError

    character(len=*), parameter :: NameSub = NameMod//'::world_init'
    !-----------------------------------------------------------------

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

  end subroutine world_init
  !============================================================================
  subroutine world_clean
    character(len=*),parameter :: NameSub=NameMod//'::world_clean'
    integer :: iError
    !------------------------------------------------------------------------
    call io_unit_clean

  end subroutine world_clean
  !===============================================================!
  subroutine world_abort(String, iValue, Value)
    character (len=*), intent(in) :: String
    integer, optional, intent(in) :: iValue
    real,    optional, intent(in) :: Value
    integer :: iError,nError
    !-------------------------------------------------------------------------
    write(*,'(a)',ADVANCE='NO') String
    if(present(iValue)) write(*,*) iValue
    if(present(Value))  write(*,*)  Value
    write(*,*)
    write(*,*)'!!! Stopping execution !!! requested by processor ', i_proc()
    call io_unit_clean
    call MPI_abort(i_comm(), nError, iError)
    stop
  end subroutine world_abort
  !===========================================================================
  logical function is_proc_world()
    is_proc_world = iProcWorld>=0
  end function is_proc_world
  !===========================================================================
  logical function is_proc0_world()
    is_proc0_world = iProcWorld==0
  end function is_proc0_world
  !===========================================================================
  integer function i_proc_world()
    i_proc_world = iProcWorld
  end function i_proc_world
  !===========================================================================
  integer function i_proc0_world()
    i_proc0_world = iProc0World
  end function i_proc0_world
  !===========================================================================
  integer function n_proc_world()
    n_proc_world = nProcWorld
  end function n_proc_world
  !===========================================================================
  integer function i_comm_world()
    i_comm_world = iCommWorld
  end function i_comm_world
  !===========================================================================
  integer function i_group_world()
    i_group_world = iGroupWorld
  end function i_group_world
  !===========================================================================
  integer function i_proc_stride_world()
    i_proc_stride_world=1
  end function i_proc_stride_world
  !===========================================================================
  integer function i_proc_last_world()
    i_proc_last_world=nProcWorld-1
  end function i_proc_last_world
  !===========================================================================
  !BOP =======================================================================
  !IROUTINE: setup_from_file - read LAYOUT.in and register components
  !INTERFACE:
  subroutine setup_from_file

    !DESCRIPTION:
    ! Reads LAYOUT.in and registers the listed components. The MPI
    ! parameters are initialized using the {\bf init} subroutine from
    ! CON\_comp\_info. The LAYOUT.in file has the following format:
    ! \begin{verbatim}
    ! remarks
    ! 
    ! Name	First	Last	Stride
    ! ==============================
    ! #COMPONENTMAP
    ! GM       0        9999    1        ! GM runs on all PE-s
    ! IE       0        2       2        ! IE runs on PE 0 and 2
    ! IM       1        1       1        ! IM runs of PE 1
    ! #END
    !
    ! remarks
    !\end{verbatim}
    ! The layout information is between the {\bf \#COMPONENTMAP}
    ! and {\bf \#END} commands. The first column contains the 
    ! two-character componentn name, the 2nd, 3rd, and 4th columns
    ! are integers containing the ranks of the first and last processors
    ! used by the component, and the stride between the processor ranks.
    ! If the last rank exceeds the number of processors used by SWMF,
    ! it is reduced to the rank of the last processor and a warning
    ! message is printed to STDOUT.
    !EOP

    character(len=*), parameter :: NameSub = NameMod//'::setup_from_file' 

    ! One line of input
    character (len=100) :: String

    ! String starting and ending comp layout 
    character (len=*), parameter  :: StringStart = "#COMPONENTMAP"
    character (len=*), parameter  :: StringEnd   = "#END" 

    ! Current line number
    integer :: nLine 

    ! Input fields in the line
    character(len=16) :: Name
    integer :: iProcZero, iProcLast, iProcStride

    ! Temporary storage for processor ranges
    integer :: iProcRange_IC(ProcZero_:ProcStride_,MaxComp)

    ! Temporary storage for component names
    character (len=lNameComp) :: Name_C(MaxComp)

    logical :: IsExisting
    integer :: lComp, iComp, iError
    !-------------------------------------------------------------------------

    ! Processor 0 looks for StringStart in NameMapFile
    if(iProcWorld==0)then
       inquire(FILE=NameMapFile, EXIST=IsExisting)
       if(.not. IsExisting) then
          write(*,'(4a)')  NameSub,' SWMF_ERROR: the map file ',NameMapFile, &
               ' does not exist in the current directory'
          call world_abort('Please create file '//NameMapFile)
       endif

       open(UNITTMP_,file=NameMapFile,iostat=iError,status="old",action="read")
       if(iError/=0) then
          write(*,'(3a)') NameSub,' SWMF_ERROR: can not open file ',NameMapFile
          call world_abort('Please check the file permissions of ' &
               //NameMapFile)
       endif

       nLine = 0
       do
          nLine = nLine + 1
          read(UNITTMP_, '(a)',IOSTAT=iError) String
          if(iError/=0)then
             close (UNITTMP_)
             write(*,'(a)') NameSub//' SWMF_ERROR: could not find '// &
                  StringStart//' in file '//NameMapFile
             call world_abort('Please edit file '//NameMapFile)
          end if
          if (String == StringStart) EXIT
       end do
    end if

    RECLOOP: do 
       ! Processor zero reads and broadcasts lines
       if(iProcWorld==0) then
          read(UNITTMP_,'(a)',IOSTAT=iError) String
          if (iError/=0) then
             write(*,'(a,i2,a)')  NameSub// &
                  ' SWMF_ERROR: can not read line after',&
                  nLine,' lines from the file '//NameMapFile
             call world_abort('Please edit '//NameMapFile)
          end if
          nLine          = nLine + 1
       end if

       call MPI_bcast(String,len(String),MPI_CHARACTER,0,iCommWorld,iError)

       if(String == StringEnd) then
          if (nComp == 0) call world_abort(NameSub // &
               'SWMF_ERROR: no components specified in the file='//NameMapFile)
          exit RECLOOP
       endif

       ! Process line in parallel
       read(String,*,iostat=iError) Name, iProcZero, iProcLast, iProcStride
       if (iError/=0 .and. iProcWorld==0) then
          write(*,'(a,i2,a)')  NameSub// &
               ' SWMF_ERROR: can not read component name and PE range from "'&
               //trim(String)//'" at line ',nLine, &
               ' from the file '//NameMapFile
          call world_abort('Please edit '//NameMapFile)
       end if

       Name = adjustl(Name)
       if (is_valid_comp_name(Name)) then
          iComp=i_comp_name(Name) 
       else
          if (iProcWorld==0) then

             write(*,'(a,i2,a)')  NameSub// &
                  ' SWMF_ERROR: invalid component name '//trim(Name)// &
                  ' at line ',nLine,' in the file '//NameMapFile
             call world_abort('Please edit '//NameMapFile)
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
          call world_abort(NameSub//' SWMF_ERROR Please edit '//NameMapFile)
          CYCLE RECLOOP
       end if

       ! Correct and store parameters read from this line
       lComp_I(iComp)                   = nComp  ! list index for named index
       iComp_C(nComp)                   = iComp  ! named index for list index
       Name_C(nComp)                    = Name

       iProcRange_IC(ProcZero_,  nComp) = min(iProcZero,iProcLast,nProcWorld-1)
       iProcRange_IC(ProcLast_,  nComp) = min(iProcLast,nProcWorld-1)
       iProcRange_IC(ProcStride_,nComp) = max(iProcStride,1)

       if(iProcWorld==0)then
          ! Report adjustments
          if(iProcRange_IC(ProcZero_,  nComp) /= iProcZero) &
               write(*,'(a,i4,a,i4)')NameSub//&
               ' SWMF_WARNING for '//trim(Name)// &
               ' iProcZero=',iProcZero,&
               ' changed to ',iProcRange_IC(ProcZero_,  nComp)
          if(iProcRange_IC(ProcLast_,  nComp) /= iProcLast) &
               write(*,'(a,i4,a,i4)')NameSub//&
               ' SWMF_WARNING for '//trim(Name)// &
               ' iProcLast=',iProcLast,&
               ' changed to ',iProcRange_IC(ProcLast_,  nComp)
          if(iProcRange_IC(ProcStride_,nComp) /= iProcStride) &
               write(*,'(a,i4,a,i4)')NameSub//&
               ' SWMF_WARNING for '//trim(Name)// &
               ' iProcStride=',iProcStride,&
               ' changed to ',iProcRange_IC(ProcStride_,nComp)
       end if
    end do RECLOOP
    if(iProcWorld == 0) close (UNITTMP_)

    allocate(CompInfo_C(nComp))

    ! check component names
    do lComp=1,nComp
       call comp_info_init(CompInfo_C(lComp),&
            Name_C(lComp),iGroupWorld, &
            iCommWorld, iProcRange_IC(:,lComp),iError)
       if(iError>0)then
          if(iProcWorld==0)write(*,'(a,3i4)')NameSub//&
               ' SWMF_ERROR: cannot produce layout for component '// &
               Name_C(lComp)//' using processor range ',&
               iProcRange_IC(:,lComp)
          call world_abort(&
               'Please edit '//NameMapFile//' or increase number of PEs')
       end if
    end do

  end subroutine setup_from_file
  !============================================================================
  integer function n_comp()
    character(len=*), parameter :: NameSub = NameMod//'::n_comp'
    !-------------------------------------------------------------
    n_comp =  nComp
    if(nComp==0) &
         write(*,'(a)') NameSub//' SWMF_ERROR: no components defined '
  end function n_comp
  !============================================================================
  integer function l_comp_id(iComp,NameSub,DoCheckRegistered)

    ! Convert named index into list index and check everything

    integer, intent(in)          :: iComp
    character(len=*), intent(in) :: NameSub ! Name of the caller
    logical, intent(in), optional:: DoCheckRegistered
    integer :: lComp
    logical :: DoCheck
    !-------------------------------------------------------------------------

    call check_i_comp(iComp,NameSub)

    lComp = lComp_I(iComp)

    DoCheck = .true.; if(present(DoCheckRegistered)) DoCheck=DoCheckRegistered
    if(DoCheck)then
       if(lComp < 1 .or. lComp>nComp) then
          write(*,'(a,i3,a,i3)')NameSub//' SWMF_ERROR iComp '// &
               NameComp_I(iComp)//' is not found among the ',&
               nComp,' registered components, lComp=',lComp
          call world_abort('Error in the caller method')
       end if
    end if

    l_comp_id = lComp

  end function l_comp_id
  !============================================================================
  integer function l_comp_name(NameComp,NameSub,DoCheckRegistered)

    ! Convert the component name into list index and check everything

    character(len=*), intent(in) :: NameComp ! Name of the component
    character(len=*), intent(in) :: NameSub  ! Name of the caller
    logical, intent(in), optional:: DoCheckRegistered
    integer :: lComp
    logical :: DoCheck
    !-------------------------------------------------------------------------
    if(.not.is_valid_comp_name(NameComp))then
       write(*,'(a)')NameSub//' SWMF_ERROR name '//NameComp// &
            ' is not a valid component name'
       call world_abort('Error in the caller method')
    end if

    lComp = lComp_I(i_comp(NameComp))

    DoCheck = .true.; if(present(DoCheckRegistered)) DoCheck=DoCheckRegistered
    if(DoCheck)then
       if(lComp < 1 .or. lComp>nComp) then
          write(*,'(a,i3,a,i3)')NameSub//' SWMF_ERROR component '// &
               NameComp//' is not found among the ',&
               nComp,' registered components, lComp=',lComp
          call world_abort('Error in the caller method')
       end if
    end if

    l_comp_name = lComp

  end function l_comp_name
  !============================================================================
  integer function i_comp_id(lComp)
    ! return named index for list index lComp
    ! Note: i_comp_name is implemented in CON_comp_param

    character(len=*), parameter :: NameSub=NameMod//'::i_comp_id'

    integer, intent(in) :: lComp
    integer :: iComp
    !--------------------------------------------------------------------------
    if(lComp<1.or.lComp>nComp)then
       write(*,'(a,i3,a,i3)')NameSub//' SWMF_ERROR lComp=',lComp,&
            ' is out of range 1..nComp=',nComp
       call world_abort('Error in caller method')
    end if
    iComp=iComp_C(lComp)
    if(iComp<1.or.iComp>MaxComp)then
       write(*,'(a,i3,a,i3,a,i3)')NameSub//' SWMF_ERROR lComp=',lComp,&
            ' gives iComp=',iComp,&
            ' which is out of range 1..MaxComp=',MaxComp
       call world_abort(NameSub//'iComp_C array is not set correctly')
    end if
    i_comp_id = iComp

  end function i_comp_id
  !==============================================================
  subroutine get_comp_info_name(NameComp, &
       iProc, nProc, iComm, iGroup, iProcZero, iProcLast, iProcStride, &
       iUnitOut, Use, Name, NameVersion, Version, CompInfo)

    character(len=*), parameter :: NameSub = NameMod//'::get_comp_info_name'

    character(len=*), intent(in) :: NameComp
    integer, optional, intent(out) :: &
         iProc, nProc, iComm, iGroup, iProcZero,  iProcLast, iProcStride, &
         iUnitOut
    logical,                     optional, intent(out) :: Use
    character(len=lNameComp),    optional, intent(out) :: Name
    character(len=lNameVersion), optional, intent(out) :: NameVersion
    real,                        optional, intent(out) :: Version
    type(CompInfoType),          optional, intent(out) :: CompInfo
    !--------------------------------------------------------------------------
    if(present(CompInfo)) CompInfo = CompInfo_C(l_comp(NameComp,NameSub))

    call comp_info_get(CompInfo_C(l_comp(NameComp,NameSub)), &
         iProc, nProc, iComm, iGroup, iProcZero, iProcLast, iProcStride, &
         iUnitOut, Use, Name, NameVersion, Version)

  end subroutine get_comp_info_name
  !BOP =======================================================================
  !IROUTINE: get_comp_info - get component information
  !INTERFACE:
  subroutine get_comp_info_id(iComp, &
       iProc, nProc, iComm, iGroup, iProcZero, iProcLast, iProcStride, &
       iUnitOut, Use, Name, NameVersion, Version, CompInfo)

    !INPUT ARGUMENTS:
    integer, intent(in) :: iComp

    !OUTPUT ARGUMENTS:
    integer, optional, intent(out) :: &
         iProc, nProc, iComm, iGroup, iProcZero,  iProcLast, iProcStride, &
         iUnitOut
    logical,                     optional, intent(out) :: Use
    character(len=lNameComp),    optional, intent(out) :: Name
    character(len=lNameVersion), optional, intent(out) :: NameVersion
    real,                        optional, intent(out) :: Version
    type(CompInfoType),          optional, intent(out) :: CompInfo

    !DESCRIPTION:
    ! Obtain any information about the component which is identified
    ! with its ID {\bf iComp}, or name (using {\bf get\_comp\_info\_name}).
    !EOP

    character(len=*), parameter :: NameSub = NameMod//'::get_comp_info_id'
    !---------------------------------------------------------------------
    if(present(CompInfo)) CompInfo = CompInfo_C(l_comp(iComp,NameSub))

    call comp_info_get(CompInfo_C(l_comp(iComp,NameSub)), &
         iProc, nProc, iComm, iGroup, iProcZero, iProcLast, iProcStride, &
         iUnitOut, Use, Name, NameVersion, Version)

  end subroutine get_comp_info_id
  !BOP ========================================================================
  !IROUTINE: put_comp_info - put component information into registry
  !INTERFACE: 
  subroutine put_comp_info_name(NameComp, iUnitOut, Use, NameVersion, Version)

    !INPUT ARGUMENTS:
    character(len=*), intent(in) :: NameComp
    integer,          optional, intent(in) :: iUnitOut    ! for STDOUT
    logical,          optional, intent(in) :: Use         ! is used
    character(len=*), optional, intent(in) :: NameVersion ! version name
    real,             optional, intent(in) :: Version     ! version number

    !DESCRIPTION:
    ! Set component information which is not set by {\bf setup}.
    ! The first argument is either the two-character component name 
    ! {\bf NameComp} or the integer component ID (using 
    ! {\bf put\_comp\_info\_id}.
    !EOP
    character(len=*), parameter :: NameSub = NameMod//'::put_comp_info_name'

    !--------------------------------------------------------------------------
    call comp_info_put(CompInfo_C(l_comp(NameComp,NameSub)), &
         iUnitOut, Use, NameVersion, Version)

  end subroutine put_comp_info_name
  !============================================================================
  subroutine put_comp_info_id(iComp, iUnitOut, Use, NameVersion, Version)

    character(len=*), parameter :: NameSub = NameMod//'::put_comp_info_id'

    integer, intent(in) :: iComp
    integer,          optional, intent(in) :: iUnitOut
    logical,          optional, intent(in) :: Use
    character(len=*), optional, intent(in) :: NameVersion
    real,             optional, intent(in) :: Version
    !--------------------------------------------------------------------------
    call comp_info_put(CompInfo_C(l_comp(iComp,NameSub)), &
         iUnitOut, Use, NameVersion, Version)

  end subroutine put_comp_info_id
  !============================================================================
  logical function use_comp_name(NameComp)
    ! Return true if component is used
    character(len=*), intent(in) :: NameComp
    character(len=*), parameter :: NameSub=NameMod//'::use_comp_name'
    integer :: lComp
    !-------------------------------------------------------------------------
    lComp = l_comp(NameComp,NameSub,.false.)
    if(lComp>0)then
       use_comp_name = CompInfo_C(lComp) % Use
    else
       use_comp_name = .false.
    end if

  end function use_comp_name
  !============================================================================
  logical function use_comp_id(iComp)

    ! Return true if component is used
    integer, intent(in) :: iComp
    character(len=*), parameter :: NameSub=NameMod//'::use_comp_id'
    integer :: lComp
    !-------------------------------------------------------------------------
    lComp = l_comp(iComp,NameSub,.false.)
    if(lComp>0)then
       use_comp_id = CompInfo_C(lComp) % Use
    else
       use_comp_id = .false.
    end if

  end function use_comp_id
  !============================================================================
  logical function is_proc_name(Name)
    character(len=*), intent(in) :: Name
    character(len=*), parameter :: NameSub=NameMod//'::is_proc_name'
    integer :: lComp
    !-----------------------------------------------------------------------
    lComp = l_comp(Name,NameSub,.false.)
    if(lComp>0)then
       is_proc_name = CompInfo_C(lComp) % iMpiParam_I(Proc_) >= 0
    else
       is_proc_name=.false.
    end if
  end function is_proc_name
  !============================================================================
  logical function is_proc_id(iComp)
    integer, intent(in) :: iComp
    character(len=*), parameter :: NameSub=NameMod//'::is_proc_id'
    integer :: lComp
    !-----------------------------------------------------------------------
    lComp = l_comp(iComp,NameSub,.false.)
    if(lComp>0)then
       is_proc_id = CompInfo_C(lComp) % iMpiParam_I(Proc_) >= 0
    else
       is_proc_id=.false.
    end if
  end function is_proc_id
  !===========================================================================
  logical function is_proc0_name(Name)
    character(len=*), intent(in) :: Name
    character(len=*), parameter :: NameSub=NameMod//'::is_proc0_name'
    !-------------------------------------------------------------------------
    is_proc0_name =  CompInfo_C(l_comp(Name,NameSub)) % iMpiParam_I(Proc_)==0
  end function is_proc0_name
  !===========================================================================
  logical function is_proc0_id(iComp)
    integer, intent(in) :: iComp
    character(len=*), parameter :: NameSub=NameMod//'::is_proc0_id'
    !-------------------------------------------------------------------------
    is_proc0_id =  CompInfo_C(l_comp(iComp,NameSub)) % iMpiParam_I(Proc_)==0
  end function is_proc0_id
  !===========================================================================
  integer function i_proc_name(Name)
    character(len=*), intent(in) :: Name
    character(len=*), parameter :: NameSub=NameMod//'::i_proc_name'
    !-------------------------------------------------------------------------
    i_proc_name =  CompInfo_C(l_comp(Name,NameSub)) % iMpiParam_I(Proc_)
  end function i_proc_name
  !===========================================================================
  integer function i_proc_id(iComp)
    integer, intent(in) :: iComp
    character(len=*), parameter :: NameSub=NameMod//'::i_proc_id'
    !-------------------------------------------------------------------------
    i_proc_id =  CompInfo_C(l_comp(iComp,NameSub)) % iMpiParam_I(Proc_)
  end function i_proc_id
  !===========================================================================
  integer function n_proc_name(Name)
    character(len=*), intent(in) :: Name
    character(len=*), parameter :: NameSub=NameMod//'::n_proc_name'
    !-------------------------------------------------------------------------
    n_proc_name =  CompInfo_C(l_comp(Name,NameSub)) % iMpiParam_I(nProc_)
  end function n_proc_name
  !===========================================================================
  integer function n_proc_id(iComp)
    integer, intent(in) :: iComp
    character(len=*), parameter :: NameSub=NameMod//'::n_proc_id'
    !-------------------------------------------------------------------------
    n_proc_id =  CompInfo_C(l_comp(iComp,NameSub)) % iMpiParam_I(nProc_)
  end function n_proc_id
  !===========================================================================
  integer function i_proc0_name(Name)
    character(len=*), intent(in) :: Name
    character(len=*), parameter :: NameSub=NameMod//'::i_proc0_name'
    !-------------------------------------------------------------------------
    i_proc0_name =  CompInfo_C(l_comp(Name,NameSub)) % iMpiParam_I(ProcZero_)
  end function i_proc0_name
  !===========================================================================
  integer function i_proc0_id(iComp)
    integer, intent(in) :: iComp
    character(len=*), parameter :: NameSub=NameMod//'::i_proc0_id'
    !-------------------------------------------------------------------------
    i_proc0_id =  CompInfo_C(l_comp(iComp,NameSub)) % iMpiParam_I(ProcZero_)
  end function i_proc0_id
  !===========================================================================
  integer function i_comm_name(Name)
    character(len=*), intent(in) :: Name
    character(len=*), parameter :: NameSub=NameMod//'::i_comm_name'
    !-------------------------------------------------------------------------
    i_comm_name =  CompInfo_C(l_comp(Name,NameSub)) % iMpiParam_I(Comm_)
  end function i_comm_name
  !===========================================================================
  integer function i_comm_id(iComp)
    integer, intent(in) :: iComp
    character(len=*), parameter :: NameSub=NameMod//'::i_comm_id'
    !-------------------------------------------------------------------------
    i_comm_id =  CompInfo_C(l_comp(iComp,NameSub)) % iMpiParam_I(Comm_)
  end function i_comm_id
  !===========================================================================
  integer function i_group_name(Name)
    character(len=*), intent(in) :: Name
    character(len=*), parameter :: NameSub=NameMod//'::i_group_name'
    !-------------------------------------------------------------------------
    i_group_name =  CompInfo_C(l_comp(Name,NameSub)) % iMpiParam_I(Group_)
  end function i_group_name
  !===========================================================================
  integer function i_group_id(iComp)
    integer, intent(in) :: iComp
    character(len=*), parameter :: NameSub=NameMod//'::i_group_id'
    !-------------------------------------------------------------------------
    i_group_id =  CompInfo_C(l_comp(iComp,NameSub)) % iMpiParam_I(Group_)
  end function i_group_id
  !===========================================================================
  integer function i_proc_stride_name(Name)
    character(len=*), intent(in) :: Name
    character(len=*), parameter :: NameSub=NameMod//'::i_proc_stride_name'
    !-------------------------------------------------------------------------
    i_proc_stride_name =  &
         CompInfo_C(l_comp(Name,NameSub)) % iMpiParam_I(ProcStride_)
  end function i_proc_stride_name
  !===========================================================================
  integer function i_proc_stride_id(iComp)
    integer, intent(in) :: iComp
    character(len=*), parameter :: NameSub=NameMod//'::i_proc_stride_id'
    !-------------------------------------------------------------------------
    i_proc_stride_id =  &
         CompInfo_C(l_comp(iComp,NameSub)) % iMpiParam_I(ProcStride_)
  end function i_proc_stride_id
  !===========================================================================
  integer function i_proc_last_name(Name)
    character(len=*),intent(in)::Name
    character(len=*),parameter :: NameSub=NameMod//'::i_proc_last_name'
    i_proc_last_name = &
         CompInfo_C(l_comp(Name,NameSub)) % iMpiParam_I(ProcLast_)
  end function i_proc_last_name
  !===========================================================================
  integer function i_proc_last_id(iComp)
    integer,intent(in)::iComp
    character(len=*),parameter :: NameSub=NameMod//'::i_proc_last_id'
    i_proc_last_id=CompInfo_C(l_comp(iComp,NameSub)) % iMpiParam_I(ProcLast_)
  end function i_proc_last_id

end module CON_world
