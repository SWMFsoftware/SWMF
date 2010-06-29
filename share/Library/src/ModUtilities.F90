!-*- mode: f90 -*- 
!^CFG COPYRIGHT UM
!
!BOP
!
!MODULE: ModUtilities - Simple Methods for CON and Components
!INTERFACE:
module ModUtilities

  !DESCRIPTION:
  ! Simple methods which are used by CON and can be used 
  ! by the science components too. 
  !
  ! This module is almost self contained. Only ModIoUnit, ModMpi, ModKind 
  ! and the external subroutine CON\_stop are used.
  !
  ! F77 and C++ codes need an F90 interface to access these utilities.
  !EOP

  implicit none

  private ! except

  public:: check_dir
  public:: fix_dir_name
  public:: flush_unit
  public:: split_string
  public:: join_string
  public:: upper_case
  public:: lower_case
  public:: sleep
  public:: check_allocate
  public:: test_mod_utility

  logical, public :: DoFlush = .true.

contains
  !BOP ========================================================================
  !ROUTINE: check_dir - check if a directory exists
  !INTERFACE:
  subroutine check_dir(NameDir)

    !USES:
    use ModIoUnit, ONLY: UNITTMP_

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
    if(i == 0) RETURN
    if(NameDir(i:i) == '/') RETURN

    if(i >= len(NameDir)) call CON_stop(NameSub// &
         "ERROR cannot append / to directory name "//NameDir)

    NameDir(i+1:i+1) = '/'

  end subroutine fix_dir_name

  !BOP ========================================================================
  !ROUTINE: flush_unit - flush output
  !INTERFACE:
  subroutine flush_unit(iUnit)

    !USES:
#ifdef compNAGF95
    use F90_UNIX_IO,only: flush 
#endif

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
  !ROUTINE: split_string - split string into array of substrings
  !INTERFACE:
  subroutine split_string(String, MaxString, String_I, nString, &
       StringSepIn, UseArraySyntaxIn)

    !INPUT ARGUMENTS:
    character(len=*),    intent(in):: String    ! string to be split
    integer,             intent(in):: MaxString ! maximum array size

    !OPTIONAL ARGUMENTS
    character, optional, intent(in):: StringSepIn      ! separator string
    logical,   optional, intent(in):: UseArraySyntaxIn ! expand Var(10:20:2)

    !OUTPUT ARGUMENTS:
    character (len=*), intent(out):: String_I(MaxString) ! array of substrings
    integer,           intent(out):: nString             ! number of substrings

    !DESCRIPTION:
    ! Cut the input string into an array of substrings. The separator
    ! character is either StringSepIn or space (default). 
    ! Multiple consecutive separator characters are treated as one.
    ! Leading and trailing spaces are ignored. For example
    !\begin{verbatim}
    ! ' IE  GM ' --> nString=2, String\_I=(/'IE','GM'/)
    !\end{verbatim}
    ! When UseArraySyntax is present, then expand strings containing
    ! parens into an array of substrings ending with numbers, e.g.
    !\begin{verbatim}
    ! 'Var(4)'      --> nString=4,  String\_I=(/'Var1','Var2','Var3','Var4'/)
    ! 'Var(11)'     --> nString=11, String\_I=(/'Var01','Var02',...,'Var11'/)
    ! 'Var(3:5)'    --> nString=3,  String\_I=(/'Var3','Var4','Var5'/)
    ! 'Var(7:11:2)' --> nString=3,  String\_I=(/'Var07','Var09','Var11'/)
    !\end{verbatim}
    !EOP

    character:: StringSep
    logical:: UseArraySyntax

    character(len=len(String)+1) :: StringTmp

    integer :: i,l

    character(len=*), parameter :: NameSub = 'split_string'
    !--------------------------------------------------------------------------
    StringSep = ' '
    if(present(StringSepIn)) StringSep = StringSepIn

    UseArraySyntax = .false.
    if(present(UseArraySyntaxIn)) UseArraySyntax = UseArraySyntaxIn
    
    nString   = 0
    StringTmp = String
    l         = len_trim(StringTmp)
    StringTmp = trim(StringTmp) // StringSep
    do
       StringTmp = adjustl(StringTmp)       ! Remove leading spaces
       i = index(StringTmp, StringSep)      ! Find end of first part   
       if(i <= 1) RETURN                    ! Nothing before the separator
       nString = nString +1                 ! Count parts

       String_I(nString) = StringTmp(1:i-1) ! Put part into string array
       StringTmp=StringTmp(i+1:l+1)         ! Delete part+separator from string

       if(UseArraySyntax) call expand_array(String_I(nString))

       if(nString == MaxString) RETURN      ! Check for maximum number of parts
    end do

  contains
    !========================================================================
    subroutine expand_array(String1)

      ! Expand String1 if it contains array syntax, e.g.
      ! "Var(04)"     to   "Var01", "Var02", "Var03", "Var04"
      ! "Var(2:4)"    to   "Var2", "Var3", "Var4"
      ! "Var(8:12:2)" to   "Var08", "Var10", "Var12"

      character(len=*), intent(inout):: String1

      character(len=len(String1)) :: String2
      character(len=6):: StringFormat

      integer:: j, k, l, m, lNum, iFirst, iLast, Di, iNum, iError
      !---------------------------------------------------------------------
      ! Find the opening paren if any
      j = index(String1,'(')
      if(j < 1) RETURN
      k = index(String1,')')

      if(k < j) &
           call CON_stop(NameSub//' missing closing paren in String='//String)

      ! Check for colon
      l = index(String1,':')
      if(l > j)then
         ! read initial index value before the first colon        
         read(String1(j+1:l-1),*,IOSTAT=iError) iFirst
         if(iError /= 0 .or. iFirst < 1) call CON_stop(NameSub// &
              ' invalid initial index value in String='//String)
      else
         iFirst = 1
         l = j
      end if
      
      ! Check for a second colon
      m = index(String1,':',back=.true.)
      if(m > l)then
         ! read index stride value after the seecond colon        
         read(String1(m+1:k-1),*,IOSTAT=iError) Di
         if(iError /= 0 .or. Di < 1) call CON_stop(NameSub// &
              ' invalid index stride value in String='//String)
      else
         Di = 1
         m  = k
      end if

      ! read the last index value between the l and m positions
      read(String1(l+1:m-1),*,IOSTAT=iError) iLast
      if(iError /= 0 .or. iLast < iFirst) call CON_stop(NameSub// &
           ' invalid maximum index value in String='//String)

      ! Set length of numerical string and the format string
      lNum = m - l - 1
      write(StringFormat,'(a,i1,a,i1,a)') "(i",lNum,".",lNum,")"

      ! Set the beginning part of the string to the variable name
      String2 = ''
      String2(1:j-1) = String1(1:j-1)

      ! Expand variable names by repating name and adding numerical value
      nString = nString - 1
      do iNum = iFirst, iLast, Di
 
         write(String2(j:j+lNum),StringFormat) iNum
         nString = nString + 1
         String_I(nString) = String2

         if(nString == MaxString) RETURN
         
      end do
    end subroutine expand_array

  end subroutine split_string

  !BOP ========================================================================
  !ROUTINE: join_string - join string array into a single string
  !INTERFACE:

  subroutine join_string(nString, String_I, String, StringSepIn)

    !INPUT ARGUMENTS:
    integer,             intent(in):: nString
    character(len=*),    intent(in):: String_I(nString)
    character, optional, intent(in):: StringSepIn
    !OUTPUT ARGUMENTS:
    character (len=*), intent(out):: String

    !DESCRIPTION:
    ! Join the input string array into one string. The parts are joined
    ! with spaces or the optional StringSepIn.
    !EOP

    character(len=*), parameter :: NameSub = 'join_string'

    character:: StringSep

    integer :: i
    !--------------------------------------------------------------------------
    if(present(StringSepIn))then
       StringSep = StringSepIn
    else
       StringSep = ' '
    endif
    String = String_I(1)
    do i = 2, nString
      String = trim(String) // StringSep // String_I(i)
    end do

  end subroutine join_string

  !BOP ========================================================================
  !ROUTINE: upper_case - convert string to all upper case
  !INTERFACE:
  subroutine upper_case(String)

    !INPUT/OUTPUT ARGUMENTS:
    character (len=*), intent(inout) :: String

    !DESCRIPTION:
    ! Change characters to upper case in String
    !EOP

    integer, parameter :: iA=ichar('a'), iZ=ichar('z'), Di=ichar('A')-iA
    integer :: i, iC
    !--------------------------------------------------------------------------
    do i = 1, len_trim(String)
       iC = ichar(String(i:i))
       if(iC >= iA .and. iC <= iZ) String(i:i) = char(iC+Di)
    end do

  end subroutine upper_case

  !BOP ========================================================================
  !ROUTINE: lower_case - convert string to all lower case
  !INTERFACE:
  subroutine lower_case(String)

    !INPUT/OUTPUT ARGUMENTS:
    character (len=*), intent(inout) :: String

    !DESCRIPTION:
    ! Change characters to lower case in String
    !EOP

    integer, parameter :: iA=ichar('A'), iZ=ichar('Z'), Di=ichar('a')-iA
    integer :: i, iC
    !--------------------------------------------------------------------------
    do i = 1, len_trim(String)
       iC = ichar(String(i:i))
       if(iC >= iA .and. iC <= iZ) String(i:i) = char(iC+Di)
    end do

  end subroutine lower_case

  !BOP ========================================================================
  !ROUTINE: sleep - sleep a given number of seconds
  !INTERFACE:
  subroutine sleep(DtCpu)
    !USES
    use ModMpi, ONLY : MPI_wtime
    use ModKind

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

    ! Test split_string, read a string containing separators 
    ! then print substrings
    ! Do this multiple times with various settings
   
    character(len=500):: String
    integer, parameter :: MaxString = 20
    integer :: nString
    character(len=30) :: String_I(MaxString)  
    integer :: iString

    character(len=*), parameter :: NameSub = 'test_mod_utility'
    !-----------------------------------------------------------------------
    write(*,'(a)') 'testing check_dir'
    write(*,'(a)') 'check directory "."'
    call check_dir('.')
    write(*,'(a)') 'check_dir returned successfully'
    write(*,'(a)') 'check directory "xxx/"'
    call check_dir('xxx/')
    write(*,*)
    write(*,'(a)') 'testing fix_dir_name'
    String = ''
    call fix_dir_name(String)
    write(*,'(a)') 'fixed empty string=' // trim(String)
    String = 'GM/BATSRUS/data'
    write(*,'(a)') 'original    string=' // trim(String)
    call fix_dir_name(String)
    write(*,'(a)') 'fixed first string=' // trim(String)
    call fix_dir_name(String)
    write(*,'(a)') 'fixed again string=' // trim(String)

    write(*,*)
    write(*,'(a)') 'testing split_string'
    String = '  a(3)  bb(04:06) c,ddd ee,ff gg(8:12:2) '
    write(*,'(a)') 'String=' // trim(String)

    call split_string(String, MaxString, String_I, nString)
    write(*,'(a,i3,a)') 'with space separator split to', nString, ' parts:'
    do iString = 1, nString
       write(*,'(a)') trim(String_I(iString))
    end do

    call split_string(String, MaxString, String_I, nString, ',')
    write(*,'(a,i3,a)') 'with comma separator split to', nString, ' parts:'
    do iString = 1, nString
       write(*,'(a)') trim(String_I(iString))
    end do

    call split_string(String, MaxString, String_I, nString, &
         UseArraySyntaxIn=.true.)
    write(*,'(a,i3,a)') 'with UseArraySyntax split to', nString,' parts:'
    do iString = 1, nString
       write(*,'(a)') trim(String_I(iString))
    end do

    write(*,*)
    write(*,'(a)') 'testing join_string'
    call join_string(nString, String_I, String, ' ')
    write(*,'(a)') 'joined string='//trim(String)

    write(*,*)
    write(*,'(a)') 'testing upper_case and lower_case'
    String = 'abCD 123:'
    write(*,'(a)') 'mixed case string='//trim(String)
    call upper_case(String)
    write(*,'(a)') 'upper case string='//trim(String)
    call lower_case(String)
    write(*,'(a)') 'lower case string='//trim(String)

  end subroutine test_mod_utility

end module ModUtilities
