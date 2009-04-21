!^CFG COPYRIGHT UM
!
!BOP
!
!MODULE: ModUtilities - Simple Methods for CON and Components
!
!DESCRIPTION:
! Simple methods which are used by CON and can be used 
! by the science components too. 
!
! This module is almost self contained. Only ModIoUnit, ModMpi, ModKind 
! and the external subroutine CON\_stop are used.
!
! F77 and C++ codes will need some F90 interface to access these utilities.

!REVISION HISTORY:
!  17Mar04 - Gabor Toth <gtoth@umich.edu> - initial prototype/prolog/code
!            Based on the self-contained part of SWMF_methods.
!EOP

module ModUtilities

  logical :: DoFlush = .true.

contains
  !BOP ========================================================================
  !ROUTINE: check_dir - check if a directory exists
  !INTERFACE:
  subroutine check_dir(NameDir)

    !USES:
    use ModIoUnit, ONLY: UNITTMP_
    implicit none

    !INPUT ARGUMENTS:
    character(len=*), intent(in) :: NameDir

    !DESCRIPTION:
    ! Check if a directory exists by trying to open a file in it.
    ! Die with an error message if the directory does not exist.
    ! Directory names are cached, so multiple calls with the same
    ! name do not repeat the check.
    !
    ! {\bf This subroutine should be called by the root PE of the component 
    ! only!}
    ! Calling the subroutine from multiple PE-s may result in a fatal error,
    ! namely one PE may delete the file written by the other PE, so the
    ! other PE thinks that the directory does not exist.
    !EOP

    character(len=*), parameter :: NameSub='check_dir'
    integer, parameter :: MaxDir=100, lNameDir=100
    integer, save :: nDir=0

    character(len=lNameDir), save :: NameDir_I(MaxDir)
    integer :: iDir, iError
    !--------------------------------------------------------------------------
    ! Only directory names shorter than lNameDir can be stored
    if(len_trim(NameDir) <= lNameDir) then
       ! Check if this directory has been checked already
       do iDir=1,nDir
          if(NameDir_I(iDir)==NameDir) RETURN
       end do

       ! Increase counter for different directory names
       nDir=nDir+1

       ! Store new name if possible. 
       if(nDir <= MaxDir) NameDir_I(nDir)=NameDir

       ! If not, warn once, and keep checking...
       if(nDir == MaxDir+1) write(*,'(a)')NameSub // &
            ' SWMF_WARNING: too many different directories!'
    end if

    ! Try to open a file in this directory
    open(UNITTMP_, file=trim(NameDir)//'.test', status='unknown', &
         iostat = iError)

    if (iError /= 0) then
       write(*,'(a,i4)')NameSub//&
            ' SWMF_ERROR: could not open file in directory '//trim(NameDir)//&
            ', iError=',iError
       call CON_stop(NameSub//' ERROR: Cannot find/write into directory '&
            //trim(NameDir))
    else
       close(UNITTMP_, status = 'DELETE')
    endif

  end subroutine check_dir

  !BOP ========================================================================
  !ROUTINE: fix_dir_name - add a slash to the end of the directory name
  !INTERFACE:
  subroutine fix_dir_name(NameDir)

    implicit none

    !INPUT/OUTPUT ARGUMENTS:
    character(len=*), intent(inout) :: NameDir

    !DESCRIPTION:
    ! Append a '/' at the end of the directory name if it is not there
    ! and the directory name is not zero length (empty string).
    !
    ! {\bf This subroutine should be called by all PE-s of the component!}
    !EOP

    character(len=*), parameter :: NameSub='fix_dir_name'
    integer :: i
    !--------------------------------------------------------------------------
    i = len_trim(NameDir)
    if(i>0 .and. NameDir(i:i) /= '/')then
       if(i >= len(NameDir)) call CON_stop(NameSub// &
            "ERROR cannot append / for directory name "//NameDir)
       NameDir(i+1:i+1) = '/'
    end if

  end subroutine fix_dir_name

  !BOP ========================================================================
  !ROUTINE: flush_unit - flush output
  !INTERFACE:
  subroutine flush_unit(iUnit)

    !USES:
#ifdef compNAGF95
    use F90_UNIX_IO,only: flush 
#endif
    implicit none

    !INPUT ARGUMENTS:
    integer, intent(in) :: iUnit

    !DESCRIPTION:
    ! Do a flush in an operating system dependent manner if DoFlush is true.
    !EOP

    integer :: iError
    !-------------------------------------------------------------------------
    if(.not.DoFlush) RETURN

#ifdef sysIRIX64
    call flush(iUnit,iError) 
#endif

#ifdef sysAIX
    call flush_(iUnit)       
#endif

#ifdef compNAGF95
    call flush(iUnit,iError)
#endif
 
#ifdef compPGF90
    call flush(iUnit)
#endif

#ifdef compXLF90
    call flush(iUnit)
#endif

#ifdef compifort
    call flush(iUnit) 
#endif

#ifdef compmpif90
    call flush(iUnit) 
#endif

#ifdef sysOSF1
    call flush(iUnit) 
#endif

#ifdef syslf95
    call flush(iUnit) 
#endif

  end subroutine flush_unit

  !BOP ========================================================================
  !ROUTINE: split_string - split space separated list into string array
  !INTERFACE:
  subroutine split_string(String, MaxString, String_I, nString, StringSepIn)

    implicit none

    !INPUT ARGUMENTS:
    character (len=*), intent(in):: String
    integer, intent(in) :: MaxString
    character, optional, intent(in):: StringSepIn
    !OUTPUT ARGUMENTS:
    character (len=*), intent(out):: String_I(MaxString)
    integer, intent(out):: nString

    !DESCRIPTION:
    ! Cut the input string into words. The words are separated with 1 or more
    ! spaces. Leading and trailing spaces are ignored. For example
    !\begin{verbatim}
    ! ' IE  GM ' --> nString=2, String\_I=(/'IE','GM'/)
    !\end{verbatim}
    !EOP

    character(len=*), parameter :: NameSub = 'split_string'

    character:: StringSep

    character(len=len(String)+1) :: StringTmp

    integer :: i,l
    !--------------------------------------------------------------------------
    if(present(StringSepIn))then
       StringSep = StringSepIn
    else
       StringSep = ' '
    endif
    nString   = 0
    StringTmp = String
    l         = len_trim(StringTmp)
    StringTmp = trim(StringTmp) // StringSep
    do
       StringTmp = adjustl(StringTmp)       ! Remove leading spaces
       i = index(StringTmp,StringSep)       ! Find end of first part   
       if(i <= 1) RETURN                    ! Nothing before the separator
       nString = nString +1                 ! Count parts

       String_I(nString) = StringTmp(1:i-1) ! Put part into string array
       StringTmp=StringTmp(i+1:l+1)         ! Delete part+separator from string

       if(nString == MaxString) RETURN      ! Check for maximum number of parts
    end do
  end subroutine split_string

  !BOP ========================================================================
  !ROUTINE: upper_case - convert string to all upper case
  !INTERFACE:
  subroutine upper_case(String)

    implicit none

    !INPUT/OUTPUT ARGUMENTS:
    character (len=*), intent(inout) :: String

    !DESCRIPTION:
    ! Change characters to upper case in String
    !EOP

    integer, parameter :: iA=ichar('a'), iZ=ichar('z'), Di=ichar('A')-iA
    integer :: i, iC
    !--------------------------------------------------------------------------
    do i = 1, len_trim(String)
       iC=ichar(String(i:i))
       if(iC>=iA.and.iC<=iZ) String(i:i)=char(iC+Di)
    end do

  end subroutine upper_case

  !BOP ========================================================================
  !ROUTINE: lower_case - convert string to all lower case
  !INTERFACE:
  subroutine lower_case(String)

    implicit none

    !INPUT/OUTPUT ARGUMENTS:
    character (len=*), intent(inout) :: String

    !DESCRIPTION:
    ! Change characters to lower case in String
    !EOP

    integer, parameter :: iA=ichar('A'), iZ=ichar('Z'), Di=ichar('a')-iA
    integer :: i, iC
    !--------------------------------------------------------------------------
    do i = 1, len_trim(String)
       iC=ichar(String(i:i))
       if(iC>=iA.and.iC<=iZ) String(i:i)=char(iC+Di)
    end do

  end subroutine lower_case

  !BOP ========================================================================
  !ROUTINE: sleep - sleep a given number of seconds
  !INTERFACE:
  subroutine sleep(DtCpu)
    !USES
    use ModMpi, ONLY : MPI_wtime
    use ModKind
    implicit none
    !INPUT ARGUMENTS:
    real, intent(in) :: DtCpu  ! CPU time to sleep (in seconds)
    !LOCAL VARIABLES:
    real(Real8_) :: tCpu0
    !DESCRIPTION:
    ! This subroutine returns after the number of seconds
    ! given in its argument.
    !EOP
    !BOC
    tCpu0 = MPI_WTIME()
    do
       if(MPI_WTIME() > tCpu0 + DtCpu) RETURN
    end do
    !EOC
  end subroutine sleep

  !BOP ========================================================================
  !ROUTINE: check_allocate - check and stop for allocation errors
  !INTERFACE:
  subroutine check_allocate(iError,NameArray)

    !INPUT ARGUMENTS:
    integer,intent(in)::iError
    character(LEN=*),intent(in)::NameArray
    !EOP
    !BOC
    if (iError > 0) call CON_stop('check_allocate F90_ERROR '// &
         'Could not allocate array '//NameArray)
    !EOC
  end subroutine check_allocate

  !============================================================================
  subroutine test_mod_utility

  ! Test split_string, read a string containing separators then print substrings
  ! Do this multiple times with various settings

   
    integer, parameter :: MaxString = 20
    integer :: nString
    character(len=30) :: String_I(MaxString)  
  
    integer, parameter :: nTest = 4
    character :: StringSep_I(nTest) = (/ ' ', ',', ':', ';'/)
    character(len=*), parameter :: String = &
            'SC EE;IH:SP,GM IM ;RB ,IE UA '
    
    !character(len=*), parameter :: NameSub = 'test_split_string'
    integer :: iTest, iString

    write(*,*) '[string]', String
    do iTest = 1, nTest
       call split_string(String, MaxString, String_I, nString, StringSep_I(iTest))
       write(*,*) '[test',iTest, ']'
       write(*,*) 'with separator ', StringSep_I(iTest), &
       		   ' split to', nString, 'part(s) as'
       do iString = 1, nString
       	  write(*,*) String_I(iString)	 
       end do	  	  
    end do
  end subroutine test_mod_utility

end module ModUtilities
