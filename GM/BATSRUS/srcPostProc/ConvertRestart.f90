!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program convert_restart_file_format

  implicit none

  integer, parameter :: Real4_=selected_real_kind(6,30)
  integer, parameter :: Real8_=selected_real_kind(12,100)

  integer, parameter :: iUnitTmp = 10, iUnitBlk = 11, iUnitOne = 12
  character(len=*), parameter :: StringRestartExt=".rst"
  character(len=*), parameter :: NameHeaderFile="restart.H"
  character(len=*), parameter :: NameDataFile="data.rst"

  integer :: nByteReal = -1                 ! Real precision in restart files
  integer :: nI=-1, nJ=-1, nK=-1, nBlock=-1 ! Grid size
  integer :: nVar=8                         ! Number of variables
  real    :: TimeSimulation = -1.0          ! Simulation time

  logical :: IsPrevStateSaved = .false.     ! True if previous step is saved
  logical :: IsBfaceSaved     = .false.     ! True if face centered B is saved

  ! One can use 'block' or 'direct' format for input and output restart files
  character (len=10)  :: TypeRestartInFile ='block'

  character(len=100) :: NameFile

  ! Variables to read/write arbitrary precision data files
  real (Real8_) :: Dt8, Time8, Dxyz8_D(3), Xyz8_D(3)
  real (Real4_) :: Dt4, Time4, Dxyz4_D(3), Xyz4_D(3)
  real (Real8_), allocatable :: State8_VC(:,:,:,:), StatePrev8_CV(:,:,:,:)
  real (Real4_), allocatable :: State4_VC(:,:,:,:), StatePrev4_CV(:,:,:,:)
  real (Real4_), allocatable :: B4_X(:,:,:), B4_Y(:,:,:), B4_Z(:,:,:)
  real (Real4_), allocatable :: B8_X(:,:,:), B8_Y(:,:,:), B8_Z(:,:,:)

  !----------------------------------------------------------------------------

  write(*,*)'ConvertRestart v1.0 (G. Toth, 2006) starting'

  call read_restart_header

  write(*,*)'Converting from TypeRestart =',trim(TypeRestartInFile)
  write(*,*)'Precision or reals          =',nByteReal
  write(*,*)'Number of state variables   =',nVar
  write(*,*)'Number of cells in a block  =',nI,'x',nJ,'x',nK
  write(*,*)'Number of blocks            =',nBlock
  write(*,*)'Simulation time             =',TimeSimulation
  if(IsBfaceSaved)write(*,*)'Face centered magnetic field was saved'
  if(IsPrevStateSaved)write(*,*)'Previous state was saved'
  write(*,*)'Conversion starting ...'

  select case(TypeRestartInFile)
  case('block')
     call write_one_restart_file
  case('one')
     call read_one_restart_file
  case default
     call my_stop('Unknown TypeRestartInFile='//TypeRestartInFile)
  end select

  write(*,*)'ConvertRestart finished'

contains

  subroutine read_restart_header

    character(len=100) :: NameCommand, NameEquation
    integer :: iError
    character(len=*), parameter:: NameSub = 'read_restart_header'
    !--------------------------------------------------------------------------

    open(iUnitTmp,file=NameHeaderFile,status='old',iostat=iError)
    if(iError /= 0) call my_stop('Could not open '//NameHeaderFile)
    do
       read(iUnitTmp,'(a)',iostat=iError) NameCommand
       if(iError /= 0) EXIT
       select case(NameCommand)
       case("#TIMESIMULATION")
          read(iUnitTmp,*) TimeSimulation
       case("#NEWRESTART")
          read(iUnitTmp,*,iostat=iError) IsBfaceSaved
          ! If constrained transport is configured out 
          ! the logical parameter is not present
          if(iError /= 0) IsBfaceSaved = .false.
       case("#PRECISION")
          read(iUnitTmp,*) nByteReal
       case("#RESTARTINFILE")
          read(iUnitTmp,*) TypeRestartInFile
       case("#EQUATION")
          read(iUnitTmp,*) NameEquation
          read(iUnitTmp,*) nVar
       case("#CHECKGRIDSIZE")
          read(iUnitTmp,*) nI
          read(iUnitTmp,*) nJ
          read(iUnitTmp,*) nK
          read(iUnitTmp,*) nBlock
       case("#NPREVIOUS")
          IsPrevStateSaved = .true.
       case default
          CYCLE
       end select
    end do
    close(iUnitTmp)

    if(TimeSimulation < 0.0)then
       write(*,*)'Invalid TimeSimulation=',TimeSimulation
       call my_stop('Missing #TIMESIMULATION command in '//NameHeaderFile)
    end if

    if(nI<0 .or. nJ<0 .or. nK < 0 .or. nBlock <0)then
       write(*,*)'Invalid nI,nJ,nK,nBlock=',nI,nJ,nK,nBlock
       call my_stop('Missing #CHECKGRIDSIZE command in '//NameHeaderFile)
    end if

    if(nByteReal==4)then
       Time4 = TimeSimulation
       allocate(State4_VC(nVar,nI,nJ,nK))
       if(IsPrevStateSaved) &
            allocate(StatePrev4_CV(nI,nJ,nK,nVar))
       if(IsBfaceSaved) &
            allocate(B4_X(nI+1,nJ,nK), B4_Y(nI,nJ+1,nK), B4_Z(nI,nJ,nK+1))
    elseif(nByteReal==8)then
       Time8 = TimeSimulation
       allocate(State8_VC(nVar,nI,nJ,nK))
       if(IsPrevStateSaved) &
            allocate(StatePrev8_CV(nI,nJ,nK,nVar))
       if(IsBfaceSaved) &
            allocate(B8_X(nI+1,nJ,nK), B8_Y(nI,nJ+1,nK), B8_Z(nI,nJ,nK+1))
    else
       write(*,*)'nByteReal =',nByteReal,' should be 4 or 8 !!!'
       call my_stop('Missing #PRECISION command in '//NameHeaderFile)
    end if

  end subroutine read_restart_header

  !============================================================================

  subroutine read_restart_file(iBlock)

    integer, intent(in) :: iBlock

    integer   :: iError, iVar
    character :: StringDigit

    character (len=*), parameter :: NameSub='read_restart_file'
    !--------------------------------------------------------------------
    write(StringDigit,'(i1)') max(5,1+int(alog10(real(iBlock))))

    write(NameFile,'(a,i'//StringDigit//'.'//StringDigit//',a)') &
         'blk',iBlock,StringRestartExt

    open(iUnitBlk, file=NameFile, status='old', form='UNFORMATTED',&
         iostat = iError)

    if(iError /= 0)call my_stop(NameSub//' could not open: '//trim(NameFile))

    if(nByteReal == 8)then
       read(iUnitBlk, iostat = iError) Dt8, Time8
       read(iUnitBlk, iostat = iError) Dxyz8_D, Xyz8_D
       read(iUnitBlk, iostat = iError) (State8_VC(iVar,:,:,:),iVar=1,nVar)
       if(IsBfaceSaved) &
            read(iUnitBlk, iostat = iError) B8_X, B8_Y, B8_Z
       if(IsPrevStateSaved) &
            read(iUnitBlk, iostat = iError) StatePrev8_CV
    else
       read(iUnitBlk, iostat = iError) Dt4, Time4
       read(iUnitBlk, iostat = iError) Dxyz4_D, Xyz4_D
       read(iUnitBlk, iostat = iError) (State4_VC(iVar,:,:,:),iVar=1,nVar)
       if(IsBfaceSaved) &
            read(iUnitBlk, iostat = iError) b4_X, b4_Y, b4_Z
       if(IsPrevStateSaved) &
            read(iUnitBlk, iostat = iError) StatePrev4_CV
    endif

    if(iError /= 0) call my_stop(NameSub//' could not read from '//NameFile)

    close(iUnitBlk)

  end subroutine read_restart_file

  !===========================================================================

  subroutine write_restart_file(iBlock)

    integer, intent(in) :: iBlock

    character (len=*), parameter :: NameSub='write_restart_file'
    integer:: iError, iVar
    character:: StringDigit
    !--------------------------------------------------------------------

    write(StringDigit,'(i1)') max(5,int(1+alog10(real(iBlock))))

    write(NameFile,'(a,i'//StringDigit//'.'//StringDigit//',a)') &
         'blk',iBlock,StringRestartExt

    open(iUnitBlk, file=NameFile, status="replace", form='UNFORMATTED', &
         iostat=iError)
    if(iError /= 0) call my_stop(NameSub//' failed opening '//NameFile)

    if(nByteReal == 8)then
       write(iUnitBlk, iostat = iError) Dt8, Time8
       write(iUnitBlk, iostat = iError) Dxyz8_D, Xyz8_D
       write(iUnitBlk, iostat = iError) (State8_VC(iVar,:,:,:),iVar=1,nVar)
       if(IsBfaceSaved) &
            write(iUnitBlk, iostat = iError) B8_X, B8_Y, B8_Z
       if(IsPrevStateSaved) &
            write(iUnitBlk, iostat = iError) StatePrev8_CV
    else
       write(iUnitBlk, iostat = iError) Dt4, Time4
       write(iUnitBlk, iostat = iError) Dxyz4_D, Xyz4_D
       write(iUnitBlk, iostat = iError) (State4_VC(iVar,:,:,:),iVar=1,nVar)
       if(IsBfaceSaved) &
            write(iUnitBlk, iostat = iError) B4_X, B4_Y, B4_Z
       if(IsPrevStateSaved) &
            write(iUnitBlk, iostat = iError) StatePrev4_CV
    endif

    if(iError /= 0)call my_stop(NameSub//' failed writing into '//NameFile)

    close(iUnitBlk)

  end subroutine write_restart_file

  !============================================================================

  subroutine open_one_restart_file(DoRead)

    logical, intent(in) :: DoRead

    integer :: lRecord, l, iError
    character(len=*), parameter :: NameSub='open_one_restart_file'
    !-------------------------------------------------------------------------

    ! Calculate the record length for the first block
    inquire (IOLENGTH = lRecord ) Dt8, Dxyz8_D, Xyz8_D, State8_VC
    if(IsBfaceSaved)then
       inquire (IOLENGTH = l) B8_X, B8_Y, B8_Z
       lRecord = lRecord + l
    end if
    if(IsPrevStateSaved)then
       inquire (IOLENGTH = l) StatePrev8_CV
       lRecord = lRecord + l
    end if

    if(nByteReal == 4) lRecord = lRecord/2

    if(DoRead)then
       open(iUnitOne, file=NameDataFile, &
            RECL = lRecord, ACCESS = 'direct', FORM = 'unformatted', &
            status = 'old', iostat=iError)
    else
       open(iUnitOne, file=NameDataFile, &
            RECL = lRecord, ACCESS = 'direct', FORM = 'unformatted', &
            status = 'replace', iostat=iError)
    end if
    if(iError /= 0)then
       write(*,*) NameSub,': ERROR for DoRead=',DoRead
       call my_stop(NameSub//': could not open file='//NameDataFile)
    end if

  end subroutine open_one_restart_file

  !============================================================================

  subroutine read_one_restart_file

    character (len=*), parameter :: NameSub='read_one_restart_file'
    integer :: iBlock
    !-------------------------------------------------------------------------

    call open_one_restart_file(DoRead = .true.)

    do iBlock = 1, nBlock

       if(nByteReal == 4)then
          if(IsBfaceSaved)then
             ! Read face centered magnetic field 
             read(iUnitOne, rec=iBlock)  &
                  Dt4, Dxyz4_D, Xyz4_D, State4_VC, &
                  B4_X, B4_Y, B4_Z
          elseif(IsPrevStateSaved)then
             ! Read previous time step for sake of BDF2 scheme
             read(iUnitOne, rec=iBlock)  &
                  Dt4, Dxyz4_D, Xyz4_D, State4_VC, &
                  StatePrev4_CV
          else
             read(iUnitOne, rec=iBlock)  &
                  Dt4, Dxyz4_D, Xyz4_D, State4_VC
          end if

       else
          if(IsBfaceSaved)then
             ! Read face centered magnetic field 
             read(iUnitOne, rec=iBlock)  &
                  Dt8, Dxyz8_D, Xyz8_D, State8_VC, &
                  B8_X, B8_Y, B8_Z
          elseif(IsPrevStateSaved)then
             ! Read previous time step for sake of BDF2 scheme
             read(iUnitOne, rec=iBlock)  &
                  Dt8, Dxyz8_D, Xyz8_D, State8_VC, &
                  StatePrev8_CV
          else
             read(iUnitOne, rec=iBlock)  &
                  Dt8, Dxyz8_D, Xyz8_D, State8_VC
          end if

       end if

       ! write it into a block file
       call write_restart_file(iBlock)

    end do

    close(iUnitOne)

  end subroutine read_one_restart_file

  !============================================================================

  subroutine write_one_restart_file

    character (len=*), parameter :: NameSub='write_one_restart_file'
    integer :: iBlock
    !--------------------------------------------------------------------

    call open_one_restart_file(DoRead = .false.)

    do iBlock = 1, nBlock

       ! read data from a block file
       call read_restart_file(iBlock)

       if(nByteReal == 4)then
          if(IsBfaceSaved)then
             ! Read face centered magnetic field 
             write(iUnitOne, rec=iBlock)  &
                  Dt4, Dxyz4_D, Xyz4_D, State4_VC, &
                  B4_X, B4_Y, B4_Z
          elseif(IsPrevStateSaved)then
             ! Read previous time step for sake of BDF2 scheme
             write(iUnitOne, rec=iBlock)  &
                  Dt4, Dxyz4_D, Xyz4_D, State4_VC, &
                  StatePrev4_CV
          else
             write(iUnitOne, rec=iBlock)  &
                  Dt4, Dxyz4_D, Xyz4_D, State4_VC
          end if
       else
          if(IsBfaceSaved)then
             ! Read face centered magnetic field 
             write(iUnitOne, rec=iBlock)  &
                  Dt8, Dxyz8_D, Xyz8_D, State8_VC, &
                  B8_X, B8_Y, B8_Z
          elseif(IsPrevStateSaved)then
             ! Read previous time step for sake of BDF2 scheme
             write(iUnitOne, rec=iBlock)  &
                  Dt8, Dxyz8_D, Xyz8_D, State8_VC, &
                  StatePrev8_CV
          else
             write(iUnitOne, rec=iBlock)  &
                  Dt8, Dxyz8_D, Xyz8_D, State8_VC
          end if
       end if

    end do

    close(iUnitOne)

  end subroutine write_one_restart_file

  !===========================================================================
  subroutine my_stop(String)
    character(len=*), intent(in) :: String
    write(*,*)'ERROR in convert_restart_file_format:'
    write(*,*) String
    stop
  end subroutine my_stop

end program convert_restart_file_format
