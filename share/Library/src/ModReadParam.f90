!BOP
!MODULE: ModReadParam - read, include, broadcast and distribute parameters
!INTERFACE:
module ModReadParam

  !DESCRIPTION:
  ! This is a library for reading parameters and distribute them between
  ! the components.
  !
  ! {\bf Subroutine read\_file(NameFile)} reads the text from file NameFile
  ! which may include further files. The files can be included with the
  ! \begin{verbatim}
  ! #INCLUDE
  ! some/filename
  ! \end{verbatim}
  ! command. Include files can be nested up to MaxNestedFile=10 levels.
  !
  ! The text either ends at the end of the file, or at the 
  ! \begin{verbatim}
  ! #END
  ! \end{verbatim}
  ! command. After reading, the text is broadcast to all processors.
  ! The text buffer contains at most MaxLine=1000 lines, which are at most
  ! lStringLine=100 character long. Normally only the control component
  ! calls {\bf read\_file}.
  !
  ! {\bf Subroutine read\_init} can select a part of the text buffer
  ! starting from line iLine+1 to the last line nLine.
  ! It also sets the session number, the name of the component which should
  ! read the selected text, and the output unit used by the compnent
  ! for echoing the parameters. Normally only the control component 
  ! calls {\bf read\_init}.
  !
  ! {\bf Function read\_line} reads the next line from the selected
  ! part of the text. If there are no more lines in the selected part,
  ! the logical function returns .false. The line itself and the line number
  ! are provided in optional output arguments.
  !
  ! \bigskip
  ! For the control component as well as for many physical components 
  ! the parameters are given in from of command lines followed by parameter 
  ! lines. The commands start with a \# character, which is usually 
  ! followed by upper case letters and numbers. Anything after a space
  ! or TAB (char(9)) is ignored. The number of parameter lines is determined 
  ! by the command.  Each parameter line contains either a logical, an
  ! integer, a real or a string variable. Comments placed after at least 
  ! 3 spaces or one TAB character are ignored. The lines which do not contain 
  ! a command or the corresponding parameters are ignored, and can be used 
  ! for remarks.  For example
  ! \begin{verbatim}
  ! #COMMAND1
  ! param1      Description1
  ! param2      Description2
  !
  ! some remark
  !
  ! #COMMAND2 some comment here
  ! param3      Description3
  ! \end{verbatim}
  ! {\bf Function read\_command} returns .true. if a command is found
  ! in the previously read line and it provides the name of the command
  ! in it output argument by stripping off anything behind a blank.
  !
  ! {\bf Subroutine read\_var} reads the parameter from the parameter line.
  ! The name of the parameter is provided as an input argument, the value
  ! is obtained from the output argument.
  ! 
  ! {\bf Subroutine read\_echo\_set} can be used to set the logical,
  ! which determines if the input parameters should be echoed back. 
  ! The default is no echo.
  !
  ! \bigskip
  ! The typical way to read parameters in form of commands and parameters
  ! in a component is the following
  ! \begin{verbatim}
  !   use ModReadParam
  !   implicit none
  !   character (len=lStringLine) :: NameCommand
  !   do
  !       if(.not.read_line() ) EXIT
  !       if(.not.read_command(NameCommand)) CYCLE
  !       select case(NameCommand)
  !       case("#IONODIR")
  !           call read_var("NameIonoDir",NameIonoDir)
  !       case("#IONOSPHERE")
  !             call read_var('iConductanceModel',iConductanceModel)
  !             call read_var('UseFullCurrent' ,UseFullCurrent)
  !             call read_var('F10.7 Flux',Flux107)
  !       ...
  !       end select
  !   end do
  ! \end{verbatim}
  ! If the component does not use the \#COMMAND parameters format,
  ! the lines can be read line by line like this
  ! \begin{verbatim}
  !   use ModReadParam
  !   implicit none
  !   character (len=lStringLine) :: StringLine
  !   do
  !       if(.not.read_line(StringLine) ) EXIT
  !       ! process StringLine
  !       read(StringLine,*) MyVariables
  !       ...
  !   end do
  ! \end{verbatim}
  ! If the component cannot process the parameters line by line,
  ! the following methods are available to read the text into a string array:
  !
  ! {\bf Function i\_line\_read()} returns the number of the line
  ! before the first line. 
  !
  ! {\bf Function n\_line\_read()} returns the line number of the last line.
  !
  ! {\bf Subroutine  read\_text} returns the selected part of the 
  ! text buffer in its output argument. 
  ! 
  ! \bigskip
  ! These methods can be used like this
  ! \begin{verbatim}
  !     use ModReadParam
  !     implicit none
  !     character(len=lStringLine), allocatable :: StringLine_I(:)
  !     allocate(StringLine_I(i_line_read()+1:n_line_read()))
  !     call read_text(StringLine_I)
  !     ! process text buffer
  !     ...
  ! \end{verbatim}

  !USES:
  use ModMpi
  use ModIoUnit, ONLY: io_unit_new, STDOUT_

  implicit none

  private ! except

  !PUBLIC DATA MEMBERS:
  integer, parameter, public :: lStringLine=100 ! Max length of input lines

  !PUBLIC MEMBER FUNCTIONS:
  public :: read_file         ! Read text string from parameter file and bcast
  public :: read_init         ! Select the appropriate section of the text
  public :: read_line         ! Read next line, return false at the end
  public :: read_command      ! Read command, return false if not a command
  public :: read_var          ! Read scalar variable of any type
  public :: read_echo_set     ! Set if parameters should be echoed
  public :: i_line_read       ! Return the current line number
  public :: n_line_read       ! Return the last line number in the selected
  public :: i_session_read    ! Return the session number
  public :: read_text         ! Provide the full text in the output argument

  !EOP

  character(len=*), parameter :: NameMod='ModReadParam'

  ! Text buffer to hold all the input lines
  integer, parameter :: MaxLine=1000
  character (len=lStringLine), save :: StringLine_I(MaxLine)

  character(len=lStringLine)        :: StringLine

  character(len=2)  :: NameComp    = '  '  ! Name of the component
  character(len=3)  :: StringPrefix= '   ' ! Prefix for echo
  integer           :: iIoUnit     = -1    ! Unit number for echo
  integer           :: iLine=0             ! Current line number
  integer           :: nLine=0             ! Last line number
  integer           :: iSession=0          ! Session number
  logical           :: DoEcho = .false.    ! Do we echo parameters?

  interface read_var
     module procedure read_var_c, read_var_l, read_var_r, read_var_i
  end interface

contains

  !BOP =======================================================================
  !IROUTINE: read_file - read parameter file
  !INTERFACE:
  subroutine read_file(NameFile,iComm,NameRestartFile)

    !INPUT ARGUMENTS:
    character (len=*), intent(in):: NameFile ! Name of the base param file
    integer,           intent(in):: iComm    ! MPI communicator for broadcast

    ! Name of the restart file to be read if a #RESTART command is found
    character (len=*), intent(in), optional :: NameRestartFile 

    !EOP
    integer, parameter :: MaxNestedFile = 10

    character(len=*), parameter :: NameSub=NameMod//'::read_file'

    character (len=lStringLine) :: NameCommand

    integer :: iUnit_I(MaxNestedFile)

    integer :: iFile, i, iError, iProc

    logical :: IsFound

    logical :: Done=.false., DoInclude
    !-----------------------------------------------------------------------
    if(Done)call CON_stop(NameSub//&
         ' ERROR: the parameter file should be read only once!')

    ! Get processor rank
    call MPI_comm_rank(iComm,iProc,iError)

    !\
    ! Read all input file(s) into memory and broadcast
    !/
    if(iProc==0)then
       nLine=0
       inquire(file=NameFile,EXIST=IsFound)
       if(.not.IsFound)call CON_stop(NameSub//' SWMF_ERROR: '//&
            trim(NameFile)//" cannot be found")
       iFile=1
       iUnit_I(iFile)=io_unit_new()
       open(iUnit_I(iFile),file=NameFile,status="old")
       do
          read(iUnit_I(iFile),'(a)',ERR=100,END=100) StringLine
          NameCommand=StringLine
          i=index(NameCommand,' '); 
          if(i>0)NameCommand(i:len(NameCommand))=' '
          i=index(NameCommand,char(9)); 
          if(i>0)NameCommand(i:len(NameCommand))=' '
          if(NameCommand=='#INCLUDE')then
             ! Include text from file following the #INCLUDE command
             read(iUnit_I(iFile),'(a)')StringLine
             DoInclude = .true.
          elseif(present(NameRestartFile) .and. NameCommand=='#RESTART')then
             ! If #RESTART command is followed by true then include
             ! the file named NameRestartFile
             read(iUnit_I(iFile),*,IOSTAT=iError)DoInclude
             if(iError>0)then
                write(*,*) NameSub,&
                     " ERROR: could not read logical after #RESTART command",&
                     " at line ",nLine+1
                call CON_stop("Correct "//trim(NameFile))
             end if
             if(DoInclude)then
                StringLine = NameRestartFile
             else
                StringLine = ' ' ! remove #RESTART command
             end if
          else
             DoInclude = .false.
          end if
          if(DoInclude)then
             iFile = iFile + 1
             if(iFile > MaxNestedFile)call CON_stop(NameSub// &
                  " SWMF_ERROR: more than MaxNestedFile nested files")
             inquire(file=StringLine,EXIST=IsFound)
             if(.not.IsFound)call CON_stop(NameSub// &
                  " SWMF_ERROR: include file cannot be found, name="//&
                  trim(StringLine))
             iUnit_I(iFile) = io_unit_new()
             open(iUnit_I(iFile),FILE=StringLine,STATUS="old")
             CYCLE
          else if(NameCommand/='#END')then
             ! Store line into buffer
             nLine=nLine+1
             if(nLine>maxline)call CON_stop(NameSub// &
                  " SWMF_ERROR: too many lines of input")
             StringLine_I(nLine)=StringLine
             CYCLE
          end if

100       continue
          close (iUnit_I(iFile))
          if(iFile > 1)then
             ! Continue reading the calling file
             iFile = iFile - 1
             CYCLE
          else
             ! The base file ended, stop reading
             EXIT
          end if
       end do
       if(nLine==0)call CON_stop(NameSub// &
            " SWMF_ERROR: no lines of input read")
    end if
    ! Broadcast the number of lines and the text itself to all processors
    call MPI_Bcast(nLine,1,MPI_INTEGER,0,iComm,iError)

    if(iError>0)call CON_stop(NameSub// &
         " MPI_ERROR: number of lines could not be broadcast")

    call MPI_Bcast(StringLine_I,lStringLine*nLine,MPI_CHARACTER,&
         0,iComm,iError)

    if(iError>0)call CON_stop(NameSub// &
         " MPI_ERROR: text could not be broadcast")

    if(iProc==0)write(*,'(a,i4,a)') NameSub// &
         ': read and broadcast nLine=',nLine,' lines of text'

    Done = .true.

  end subroutine read_file

  !===========================================================================
  subroutine read_init(NameCompIn, iSessionIn, iLineIn, nLineIn, iIoUnitIn)

    ! Initialize module variables

    character (len=2), intent(in) :: NameCompIn
    integer,           intent(in) :: iSessionIn, iLineIn
    integer, optional, intent(in) :: nLineIn, iIoUnitIn
    !------------------------------------------------------------------------
    NameComp     = NameCompIn
    iSession     = iSessionIn
    iLine        = iLineIn
    if(present(nLineIn))then
       nLine     = nLineIn
    else
       nLine     = size(StringLine_I)
    end if
    if(present(iIoUnitIn))then
       iIoUnit   = iIoUnitIn
    else
       iIoUnit   = STDOUT_
    end if
    if(iIoUnit==STDOUT_ .and. len_trim(NameComp)>0 )then
       StringPrefix = NameComp//': '
    else
       StringPrefix = ''
    end if
    
  end subroutine read_init

  !===========================================================================
  subroutine read_echo_set(DoEchoIn)

    logical, intent(in) :: DoEchoIn

    DoEcho = DoEchoIn

  end subroutine read_echo_set

  !BOP =======================================================================
  !IROUTINE: read_line - read the next line from the text buffer
  !INTERFACE:
  logical function read_line(StringLineOut, iLineOut)

    !OUTPUT ARGUMENTS:
    character (len=*), optional, intent(out) :: StringLineOut
    integer, optional, intent(out)           :: iLineOut

    !DESCRIPTION:
    ! Read the current line from StringLine\_I into StringLine,
    ! set the optional StringLineOut and iLineOut arguments.
    ! Return .true. if successful, otherwise (if there are
    ! no more lines in the selected part of the text buffer)
    ! return .false. and an empty string in StringLineOut if present.
    !EOP

    iLine=iLine+1
    if(present(iLineOut)) iLineOut = iLine
    if(iLine <= nLine)then
       StringLine = StringLine_I(iLine)
       if(present(StringLineOut)) StringLineOut = StringLine
       read_line  = .true.
    else
       if(present(StringLineOut)) StringLineOut = ''
       read_line = .false.
    endif

  end function read_line

  !BOP =======================================================================
  !IROUTINE: read_command - read the name of the command from the current line
  !INTERFACE:
  logical function read_command(NameCommand)

    !OUTPUT ARGUMENTS:
    character (len=*), intent(out) :: NameCommand

    !DESCRIPTION:
    ! If the current line contains a command name (starting with \#),
    ! return true, and put the name of the command into the 
    ! output argument. Otherwise return .false. and an empty string.
    !EOP

    integer :: i

    !------------------------------------------------------------------------
    if(StringLine(1:1)=="#")then

       if(DoEcho)then
          write(iIoUnit,'(a)')trim(StringPrefix)
          write(iIoUnit,'(a)')trim(StringPrefix)//trim(StringLine)
       end if

       ! Remove anything after a space or TAB
       i=index(StringLine,' ');     if(i>0)StringLine(i:len(StringLine))=' '
       i=index(StringLine,char(9)); if(i>0)StringLine(i:len(StringLine))=' '

       NameCommand = StringLine
       read_command = .true.
    else
       NameCommand  = ''
       read_command = .false.
    endif

  end function read_command

  !===========================================================================
  subroutine read_line_param(Type,Name,iError)

    ! read next line from text

    character (len=*), intent(in) :: Type, Name
    integer, optional, intent(out):: iError
    !------------------------------------------------------------------------

    if(present(iError))iError=0
    iLine=iLine+1
    if(iLine>nLine)then
       if(present(iError))then
          iError=1
       else
          call CON_stop(&
               'Unexpected end of text after line='//StringLine)
       end if
    end if
    StringLine=StringLine_I(iLine)

    if(len_trim(StringLine)==0)then
       if(present(iError))then
          iError=2
       else
          write(*,'(a,i8,a,a,a,i2)') 'Error at line iLine=',iLine, &
               ' in component ',NameComp,' in session',iSession
          call CON_stop(&
               'Empty line found while reading '//Type//' variable '//Name)
       end if
    endif

  end subroutine read_line_param

  !===========================================================================
  subroutine read_echo(IsSame,Name)

    logical, intent(in) :: IsSame
    character (len=*), intent(in)    :: Name
    integer :: nTrimLine
    !-------------------------------------------------------------------------

    if(index(StringLine,Name)<1)then
       nTrimLine=len_trim(StringLine)
       StringLine=StringLine(1:nTrimLine)//char(9)//char(9)//Name
    end if
    if(IsSame)then
       nTrimLine=len_trim(StringLine)
       StringLine=StringLine(1:nTrimLine)//' (default/unchanged)'
    end if
    nTrimLine=len_trim(StringLine)
    write(iIoUnit,'(a)') trim(StringPrefix)//StringLine(1:nTrimLine)
       
  end subroutine read_echo

  !===========================================================================
  subroutine read_error(Type,Name,iError)

    ! Print error message for reading error of variable named Name of type Type

    character (len=*), intent(in) :: Type, Name
    integer, optional, intent(out):: iError

    !------------------------------------------------------------------------

    if(present(iError))then
       select case(Type)
       case('integer')
          iError = 3
       case('logical')
          iError = 4
       case('real')
          iError = 5
       case('character')
          iError = 6
       case default
          iError = -1
       end select
    else
       write(*,'(a,i3)')'Error in component '//NameComp//' in session',iSession
       call CON_stop('Error reading '//Type//' variable '//Name// &
            ' from line='//StringLine)
    end if

  end subroutine read_error
  !===========================================================================
  subroutine read_var_c(Name,StringVar,iError)
    
    ! Read a string variable described by Name
    
    implicit none

    ! Arguments
    character (len=*), intent(in)   :: Name
    character (len=*), intent(inout):: StringVar
    integer, optional, intent(out)  :: iError

    ! Local variable
    integer :: i,j

    !-------------------------------------------------------------------------

    call read_line_param('character',Name,iError)

    if(DoEcho)call read_echo(StringLine(1:len(StringVar))==StringVar,Name)

    ! Get rid of trailing comments after a TAB character or 3 spaces
    i=index(StringLine,char(9))
    j=index(StringLine,'   ')
    if(i>1.and.j>1)then
       i = min(i,j)      ! Both TAB and 3 spaces found, take the closer one
    else
       i = max(i,j)      ! Take the one that was found (if any)
    end if

    if(i>0)StringLine(i:len(StringLine))=' '

    StringVar=StringLine(1:len(StringVar))

  end subroutine read_var_c

  !BOP ========================================================================
  !IROUTINE: read_var - read a variable following the command.
  !INTERFACE:
  subroutine read_var_i(Name,IntVar,iError)

    !INPUT ARGUMENTS:
    character (len=*), intent(in)   :: Name
    !OUTPUT ARGUMENTS:
    integer,           intent(inout):: IntVar
    integer, optional, intent(out)  :: iError

    !DESCRIPTION:
    ! Read a variable from the next line in the buffer.
    ! The variable name is given by the string Name, which is used in the
    ! echoing of the parameters as well as in error messages. 
    ! The value is returned in the non-optional output argument. 
    ! If the optional argument iError is present, it returns a non-zero
    ! value in case an error occurs. If iError is not present, all errors
    ! result in an error message and an abort of the run.
    ! There are four variants of this subroutine: for integer, real,
    ! character string and logical variable types.
    !EOP

    ! Local variable
    integer :: IntTmp

    !-------------------------------------------------------------------------
    call read_line_param('integer',Name,iError)

    read(StringLine,*,err=1) IntTmp
    if(DoEcho)call read_echo(IntTmp==IntVar,Name)
    IntVar=IntTmp

    return
1   continue
    call read_error('integer',Name,iError)

  end subroutine read_var_i

  !===========================================================================
  subroutine read_var_r(Name,RealVar,iError)

    ! Read a real variable described by Name
    
    implicit none

    ! Arguments
    character (len=*), intent(in)   :: Name
    real, intent(inout)             :: RealVar
    integer, optional, intent(out)  :: iError

    ! Local variable
    real :: RealTmp

    !-------------------------------------------------------------------------
    call read_line_param('real',Name,iError)

    read(StringLine,*,err=1) RealTmp
    if(DoEcho)call read_echo(RealTmp==RealVar,Name)
    RealVar=RealTmp

    return
1   continue
    call read_error('real',Name,iError)

  end subroutine read_var_r

  !===========================================================================
  subroutine read_var_l(Name,IsLogicVar,iError)

    ! Read a logical variable described by Name

    ! Arguments
    character (len=*), intent(in)   :: Name
    logical, intent(inout)          :: IsLogicVar
    integer, optional, intent(out)  :: iError

    ! Local variable
    logical :: IsLogicTmp

    !-------------------------------------------------------------------------
    call read_line_param('logical',Name,iError)

    read(StringLine,*,err=1) IsLogicTmp
    if(DoEcho)call read_echo(IsLogicTmp.eqv.IsLogicVar,Name)
    IsLogicVar=IsLogicTmp

    return
1   continue
    call read_error('logical',Name,iError)

  end subroutine read_var_l
  !===========================================================================
  integer function i_line_read()
    i_line_read = iLine
  end function i_line_read
  !===========================================================================
  integer function n_line_read()
    n_line_read = nLine
  end function n_line_read
  !===========================================================================
  integer function i_session_read()
    i_session_read = iSession
  end function i_session_read
  !BOP =======================================================================
  !IROUTINE: read_text - obtain selected text buffer 
  !INTERFACE:
  subroutine read_text(String_I)
    !OUTPUT ARGUMENTS:
    character(len=lStringLine), intent(out) :: String_I(iLine+1:nLine)
    !EOP
    !BOC
    String_I = StringLine_I(iLine+1:nLine)
    !EOC
  end subroutine read_text
  !===========================================================================
end module ModReadParam
