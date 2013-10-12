!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!               Space Weather Modeling Framework (SWMF)                !
!    Center for Space Environment Modeling, The University of Michigan !
!-----------------------------------------------------------------------
!BOP -------------------------------------------------------------------------
!
!MODULE: CON_comp_info - Component information creation and handling
! 
!DESCRIPTION: 
! This module provides a public derived type {\bf CompInfoType}
! to store component information such as name, version name and number, 
! whether the component is used, the unit number for (redirected) STDOUT, 
! and 7 MPI parameters, such as the MPI communicator, MPI group, 
! number of processors used by the component, the global rank of the 
! first and last processors used by the component, the stride between
! the used processors, and the rank of the current PE 
! within the component group.
!
! The MPI information is created with the {\bf init} subroutine.
!
! Other data can be set with the {\bf put} subroutine.
!
! All data can be accessed with the {\bf get} subroutine.
!
! The data can be released with the {\bf clean} subroutine.

!INTERFACE:	
!
module CON_comp_info
  !
  !USES:
  !
  use CON_comp_param
  use ModMpi
  use ModIoUnit, ONLY : STDOUT_
  implicit none

  private ! except

  !
  !PUBLIC TYPES:
  !
  public :: CompInfoType ! The class data structure

  type CompInfoType
     character (len=lNameComp)        :: Name
     character (len=lNameVersion)     :: NameVersion
     real                             :: Version
     logical                          :: Use
     integer                          :: iUnitOut
     integer, dimension(nMpiParam)    :: iMpiParam_I
  end type CompInfoType

  !
  !PUBLIC MEMBER FUNCTIONS:
  public :: init   ! constructor creates layout information from PE range
  public :: clean  ! destructor 
  public :: get    ! get component information
  public :: put    ! put component information not set by init

  !
  !REVISION HISTORY: 
  !
  !  June 2003 - O. Volberg <volov@umich.edu> - initial version
  !  July 2003 - G. Toth <gtoth@umich.edu>    - major rewrite
  !  Aug  2003 - G. Toth <gtoth@umich.edu>    - improved description
  !
  !EOP ------------------------------------------------------------------------

  character(len=*),parameter :: NameMod='CON_comp_info'

contains

  subroutine init(Info, Name, iGroupWorld, iCommWorld, iProcRange_I, iError)

    type(CompInfoType), intent(inout) :: Info
    character (len=*), intent(in)     :: Name
    integer, intent(in)               :: iGroupWorld,iCommWorld
    ! Processor range contains first rank, last rank and stride
    ! indicating ranks in group or processors to be included in the new group 
    integer, intent(in)               :: iProcRange_I(ProcZero_:ProcStride_) 
    integer, intent(out)              :: iError

    ! local variables
    integer :: iProcWorld, nProcWorld

    character(len=*), parameter :: NameSub = NameMod//'::init'
    !------------------------------------------------------------------------

    iError = 0

    ! Check if the name of the component is valid
    if (is_valid_comp_name(Name)) then
       Info % Name = Name
    else
       iError = 1
       write(*,'(a)')NameSub//' SWMF_ERROR: invalid component name '//Name
       RETURN
    end if

    ! get the number of PE-s in the world group
    call MPI_group_size(iGroupWorld, nProcWorld, iError)
    if(iError /= 0) then
       write(*,'(a)')NameSub//' MPI_ERROR: in MPI_group_size'
       RETURN
    endif

    ! Check if the range is correct and fits into the world group
    if ( iProcRange_I(ProcZero_) < 0 .or. &
         iProcRange_I(ProcZero_) > iProcRange_I(ProcLast_) .or. &
         iProcRange_I(ProcLast_) > nProcWorld - 1 .or. &
         iProcRange_I(ProcStride_) < 0  ) then
       iError = 1
       write(*,'(a,3i4,a,i4,a)')NameSub// &
            ' SWMF_ERROR: incorrect iProcRange_I = ',iProcRange_I,&
            ' when running on ',nProcWorld,' processor(s)'
       RETURN
    endif

    !Create the group from the processor range
    call MPI_group_range_incl(iGroupWorld, 1, iProcRange_I, &
         Info%iMpiParam_I(Group_), iError)
    if(iError /= 0) then
       write(*,'(a)')NameSub//' MPI_ERROR: in MPI_group_range_incl'
       RETURN
    endif

    ! Store the processor range
    Info%iMpiParam_I(ProcZero_:ProcStride_) = iProcRange_I

    ! Get the rank of a PE withing this group
    call MPI_group_rank(Info%iMpiParam_I(Group_), &
         Info%iMpiParam_I(Proc_), iError)

    if(iError /= 0) then
       write(*,'(a)')NameSub//' MPI_ERROR: in  MPI_group_rank'
       RETURN
    endif

    ! Create the communicator for this group
    call MPI_comm_create(iCommWorld, Info%iMpiParam_I(Group_), &
         Info%iMpiParam_I(Comm_), iError)

    if(iError /= 0) then
       write(*,'(a)')NameSub//' MPI_ERROR: in MPI_comm_create'
       RETURN
    endif

    ! get  the rank in the world communicator 
    call MPI_comm_rank(iCommWorld, iProcWorld, iError)

    ! Let the root PE of the group calculate the size of the group
    if (iProcWorld == Info%iMpiParam_I(ProcZero_) ) then
       call MPI_comm_size(Info%iMpiParam_I(Comm_), &
            Info%iMpiParam_I(nProc_), iError)

       if(iError /= 0) then
          write(*,'(a)')NameSub//'MPI_ERROR:  in  MPI_comm_size '
          RETURN
       endif
    end if

    ! Broadcast the group size to all
    call MPI_bcast(Info%iMpiParam_I(nProc_), 1, &
         MPI_INTEGER, Info%iMpiParam_I(ProcZero_), iCommWorld, iError)

    if(iError /= 0) then
       write(*,'(a)')NameSub//' MPI_ERROR: in MPI_bcast '
       RETURN
    endif

    Info % iUnitOut    = STDOUT_
    Info % Use         = .true.
    Info % NameVersion = 'unset'
    Info % Version     = -1.0

  end subroutine init

  !============================================================
  subroutine clean(Info, iError)

    type(CompInfoType), intent(inout) :: Info
    integer, intent(out) :: iError ! return code
    character(len=*), parameter :: NameSub = NameMod//':: clean_' 
    !------------------------------------------------------------------------
    call MPI_group_free(Info%iMpiParam_I(Group_), iError)
    if (iError /=  MPI_SUCCESS) then
       write(*,'(a)')NameSub// &
            ' MPI_ERROR: in MPI_group_free - can not free the group '
       RETURN
    end if
    Info%Name       = ""
    Info%NameVersion = ""
    Info%Version     = 0.0
    Info%iMpiParam_I = MPI_UNDEFINED

  end subroutine clean

  !===========================================================================
  subroutine get(Info, iProc, nProc, iComm, iGroup, iFirst, iLast, iStride, &
       iUnitOut, Use, Name, NameVersion, Version)

    type(CompInfoType), intent(in) :: Info
    integer,                     optional, intent(out) :: &
         iProc, nProc, iComm, iGroup,iFirst,  iLast, iStride, iUnitOut
    logical,                     optional, intent(out) :: Use
    character(len=lNameComp),    optional, intent(out) :: Name
    character(len=lNameVersion), optional, intent(out) :: NameVersion
    real,                        optional, intent(out) :: Version
    !-------------------------------------------------------------------------
    if( present(iProc)   ) iProc   = Info%iMpiParam_I(Proc_) 
    if( present(nProc)   ) nProc   = Info%iMpiParam_I(nProc_)
    if( present(iComm)   ) iComm   = Info%iMpiParam_I(Comm_)
    if( present(iGroup)  ) iGroup  = Info%iMpiParam_I(Group_)
    if( present(iFirst)  ) iFirst  = Info%iMpiParam_I(ProcZero_)
    if( present(iLast)   ) iLast   = Info%iMpiParam_I(ProcLast_)
    if( present(iStride) ) iStride = Info%iMpiParam_I(ProcStride_)
    if( present(iUnitOut)    ) iUnitOut     = Info%iUnitOut
    if( present(Use)         ) Use          = Info%Use
    if( present(Name)        ) Name         = Info%Name
    if( present(NameVersion) ) NameVersion  = Info%NameVersion
    if( present(Version)     ) Version      = Info%Version
  end subroutine get
  !===========================================================================
  subroutine put(Info, iUnitOut, Use, NameVersion, Version)

    type(CompInfoType), intent(inout) :: Info
    integer,          optional, intent(in) :: iUnitOut
    logical,          optional, intent(in) :: Use
    character(len=*), optional, intent(in) :: NameVersion
    real,             optional, intent(in) :: Version
    !-------------------------------------------------------------------------
    if( present(iUnitOut)    ) Info%iUnitOut    = iUnitOut
    if( present(Use)         ) Info%Use         = Use
    if( present(NameVersion) ) Info%NameVersion = NameVersion
    if( present(Version)     ) Info%Version     = Version

  end subroutine put

end module CON_comp_info
