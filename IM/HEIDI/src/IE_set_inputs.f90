
subroutine IE_set_inputs(StringInputLines)

  use ModIE_Interface
  use ModFiles
  use ModExtras

  implicit none

  character (len=100), dimension(*), intent(in) :: StringInputLines
  character (len=100) :: StringLine
  logical :: IsDone
  integer :: iLine, iDebugProc

  iLine = 1

  IsDone = .false.

  do while (.not.IsDone)

     StringLine = StringInputLines(iLine)

     if (StringLine(1:1) == "#") then

        if (index(StringLine,"#BACKGROUND") > 0) then

           call read_in_string(IE_NameOfModelDir)
           call read_in_string(IE_NameOfEFieldModel)
           call read_in_string(IE_NameOfAuroralModel)
           call read_in_string(IE_NameOfSolarModel)

           if (index(IE_NameOfAuroralModel,'IHP') > 0) &
                IE_NameOfAuroralModel = 'ihp'
           if (index(IE_NameOfAuroralModel,'PEM') > 0) &
                IE_NameOfAuroralModel = 'pem'

           if (index(IE_NameOfEFieldModel,'AMIE') > 0) &
                IE_NameOfEFieldModel = 'amie'

           if (index(IE_NameOfEFieldModel,'weimer01') > 0) &
                IE_NameOfEFieldModel = 'weimer01'
           if (index(IE_NameOfEFieldModel,'Weimer01') > 0) &
                IE_NameOfEFieldModel = 'weimer01'
           if (index(IE_NameOfEFieldModel,'WEIMER01') > 0) &
                IE_NameOfEFieldModel = 'weimer01'

           if (index(IE_NameOfEFieldModel,'weimer') > 0 .and. &
                index(IE_NameOfEFieldModel,'01') == 0) &
                IE_NameOfEFieldModel = 'weimer96'
           if (index(IE_NameOfEFieldModel,'Weimer') > 0 .and. &
                index(IE_NameOfEFieldModel,'01') == 0) &
                IE_NameOfEFieldModel = 'weimer96'
           if (index(IE_NameOfEFieldModel,'WEIMER') > 0 .and. &
                index(IE_NameOfEFieldModel,'01') == 0) &
                IE_NameOfEFieldModel = 'weimer96'

           if (index(IE_NameOfEFieldModel,'weimer96') > 0) &
                IE_NameOfEFieldModel = 'weimer96'
           if (index(IE_NameOfEFieldModel,'Weimer96') > 0) &
                IE_NameOfEFieldModel = 'weimer96'
           if (index(IE_NameOfEFieldModel,'WEIMER96') > 0) &
                IE_NameOfEFieldModel = 'weimer96'

           if (index(IE_NameOfEFieldModel,'SAMIE') > 0) &
                IE_NameOfEFieldModel = 'samie'

        endif

        if (index(StringLine,"#AMIEFILES") > 0) then
           call read_in_string(AMIEFileNorth)
           call read_in_string(AMIEFileSouth)
           IE_NameOfEFieldModel = 'amie'
           IE_NameOfAuroralModel = 'amie'
        endif

        if (index(StringLine,"#DEBUG") > 0) then
           call read_in_int(iDebugLevel)
           call read_in_int(iDebugProc)
           if (iDebugProc >= 0 .and. iProc /= iDebugProc) then
              iDebugLevel = -1
           endif
        endif

        if (index(StringLine,"#END") > 0) then
           IsDone = .true.
        endif

        if (iLine >= MaxInputLines) then
           IsDone = .true.
        endif

     else

        iLine = iLine + 1

     endif

  enddo

contains

  subroutine read_in_int(variable)
    integer, intent(out) :: variable
    iline = iline + 1
    read(StringInputLines(iline),*) variable
  end subroutine read_in_int

  subroutine read_in_logical(variable)
    logical, intent(out) :: variable
    iline = iline + 1
    read(StringInputLines(iline),*) variable
  end subroutine read_in_logical

  subroutine read_in_string(variable)
    character (len=100), intent(out) :: variable
    iline = iline + 1
    variable = StringInputLines(iline)
  end subroutine read_in_string

  subroutine read_in_real(variable)
    real :: variable
    iline = iline + 1
    read(StringInputLines(iline),*) variable
  end subroutine read_in_real

  subroutine read_in_time(variable)
    real*8 :: variable
    integer, dimension(7) :: itime
    iline = iline + 1
    read(StringInputLines(iline),*) itime(1:6)
    itime(7) = 0
    call time_int_to_real(itime, variable)
  end subroutine read_in_time

end subroutine IE_set_inputs
