!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModBlockData

  use ModSize,   ONLY: MaxBlock
  use ModMain,   ONLY: ProcTest, BlkTest
  use ModProcMH, ONLY: iProc

  implicit none

  private ! except

  public put_block_data    ! store 1 or more values into storage
  interface put_block_data
     module procedure put_point, put_array, put_array2, put_array3, put_array4
  end interface

  public set_block_data    ! indicate that all block data has been set

  public use_block_data    ! function to check if block data has been set

  public get_block_data    ! get 1 or more values from storage
  interface get_block_data
     module procedure get_point, get_array, get_array2, get_array3, get_array4
  end interface

  public n_block_data      ! function for number of data values stored
  interface n_block_data
     module procedure n_data, max_data
  end interface

  public clean_block_data
  interface clean_block_data
     module procedure clean_block, clean_all
  end interface

  public write_block_restart_files
  public read_block_restart_files

  public test_block_data

  ! Maximum number of reals associated with a block. 
  ! This has to be set in ModUser so that load_balance.f90 knows about it.
  integer, public:: MaxBlockData = 0

  ! These arrays can be initialized
  integer :: nData_B(MaxBlock) = -1        ! Number of data elements
  integer :: iData_B(MaxBlock) = -1        ! Current position for put/get
  logical :: UseData_B(MaxBlock) = .false. ! Is the data usable?

  ! Allocatable storage type for block data
  type BlockDataType
     real, pointer :: Array_I(:)
  end type BlockDataType

  ! Array of allocatable storage
  type(BlockDataType) :: Data_B(MaxBlock) 

  character(len=*), parameter :: NameMod = 'ModBlockData'

  logical, parameter :: DoDebug = .false.

contains

  !===========================================================================

  subroutine init_block(iBlock,nValue)
    integer, intent(in) :: iBlock
    integer, intent(in) :: nValue
    character (len=*), parameter :: NameSub = NameMod//'::init_block'
    !------------------------------------------------------------------------
    
    allocate(Data_B(iBlock) % Array_I(nValue))
    iData_B(iBlock)   = 0
    nData_B(iBlock)   = 0
    UseData_B(iBlock) = .false.

    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' finished'
    endif

  end subroutine init_block

  !===========================================================================

  subroutine extend_array(iBlock, nValue)
    integer, intent(in) :: iBlock
    integer, intent(in) :: nValue

    integer :: nSize
    real, pointer :: DataTemp_I(:)
    character (len=*), parameter :: NameSub = NameMod//'::extend_array'
    !------------------------------------------------------------------------
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' called with nValue=',nValue
    endif

    ! Figure out new size
    nSize = max(2*size(Data_B(iBlock) % Array_I), nData_B(iBlock) + nValue)

    ! Allocate new storage
    allocate(DataTemp_I(nSize))

    ! Copy old data into new storage
    DataTemp_I(1:nData_B(iBlock)) = Data_B(iBlock) % Array_I(1:nData_B(iBlock))

    ! Deallocate old storage
    deallocate(Data_B(iBlock) % Array_I)

    ! Set pointer to new storage
    Data_B(iBlock) % Array_I => DataTemp_I

    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' allocated size=',size(DataTemp_I)
    endif

  end subroutine extend_array

  !===========================================================================

  subroutine put_point(iBlock, Value, DoAllowReplace)
    integer, intent(in) :: iBlock
    real,    intent(in) :: Value
    logical, intent(in), optional :: DoAllowReplace

    character (len=*), parameter :: NameSub=NameMod//'::put_point'
    !------------------------------------------------------------------------
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' called'
    endif

    if(UseData_B(iBlock) .and. .not. present(DoAllowReplace)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock
       call stop_mpi('UseData_B=.true. in '//NameSub)
    end if

    if(iData_B(iBlock)+1 > nData_B(iBlock) &
         .or. .not.present(DoAllowReplace))then
       if(nData_B(iBlock) < 0)call init_block(iBlock,1)
       if(nData_B(iBlock)+1 > size(Data_B(iBlock) % Array_I)) &
            call extend_array(iBlock,1)
       nData_B(iBlock) = nData_B(iBlock)+1
    end if

    iData_B(iBlock) = iData_B(iBlock)+1
    Data_B(iBlock) % Array_I(iData_B(iBlock)) = Value

    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' returning Value=',Value
    endif

  end subroutine put_point

  !===========================================================================

  subroutine put_array(iBlock, nValue, Value_I, DoAllowReplace)
    integer, intent(in) :: iBlock
    integer, intent(in) :: nValue
    real,    intent(in) :: Value_I(nValue)
    logical, intent(in), optional :: DoAllowReplace

    integer :: i
    character (len=*), parameter :: NameSub=NameMod//'::put_array'
    !------------------------------------------------------------------------
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' called with nValue=',nValue
    endif

    if(UseData_B(iBlock) .and. .not. present(DoAllowReplace)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock
       call stop_mpi('UseData_B=.true. in '//NameSub)
    end if

    if(iData_B(iBlock) + nValue > nData_B(iBlock) &
         .or. .not.present(DoAllowReplace))then
       if(nData_B(iBlock) < 0)call init_block(iBlock,1)
       if(nData_B(iBlock)+nValue > size(Data_B(iBlock) % Array_I)) &
            call extend_array(iBlock, nValue)
       nData_B(iBlock) = nData_B(iBlock)+nValue
    end if

    i = iData_B(iBlock)
    iData_B(iBlock) = i + nValue
    Data_B(iBlock) % Array_I(i+1:i+nValue) = Value_I

  end subroutine put_array

  !===========================================================================

  subroutine put_array2(iBlock, nI, nJ, Value_II, DoAllowReplace)
    integer, intent(in) :: iBlock
    integer, intent(in) :: nI, nJ
    real,    intent(in) :: Value_II(nI, nJ)
    logical, intent(in), optional :: DoAllowReplace

    integer :: i, j, nValue
    character (len=*), parameter :: NameSub=NameMod//'::put_array2'
    !------------------------------------------------------------------------
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' called with nI, nJ=',nI, nJ
    endif

    if(UseData_B(iBlock) .and. .not. present(DoAllowReplace)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock
       call stop_mpi('UseData_B=.true. in '//NameSub)
    end if

    nValue = nI*nJ

    if(iData_B(iBlock) + nValue > nData_B(iBlock) &
         .or. .not.present(DoAllowReplace))then
       if(nData_B(iBlock) < 0)call init_block(iBlock,1)
       if(nData_B(iBlock)+nValue > size(Data_B(iBlock) % Array_I)) &
            call extend_array(iBlock, nValue)
       nData_B(iBlock) = nData_B(iBlock)+nValue
    end if

    i = iData_B(iBlock)
    do j = 1, nJ
       Data_B(iBlock) % Array_I(i+1:i+nI) = Value_II(:,j)
       i = i+nI
    end do
    iData_B(iBlock) = i

  end subroutine put_array2

  !===========================================================================

  subroutine put_array3(iBlock, nI, nJ, nK, Value_III, DoAllowReplace)
    integer, intent(in) :: iBlock
    integer, intent(in) :: nI, nJ, nK
    real,    intent(in) :: Value_III(nI, nJ, nK)
    logical, intent(in), optional :: DoAllowReplace

    integer :: i, j, k, nValue
    character (len=*), parameter :: NameSub=NameMod//'::put_array3'
    !------------------------------------------------------------------------
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' called with nI, nJ, nK=',nI, nJ, nK
    endif

    if(UseData_B(iBlock) .and. .not. present(DoAllowReplace)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock
       call stop_mpi('UseData_B=.true. in '//NameSub)
    end if

    nValue = nI*nJ*nK

    if(iData_B(iBlock) + nValue > nData_B(iBlock) &
         .or. .not.present(DoAllowReplace))then
       if(nData_B(iBlock) < 0)call init_block(iBlock,1)
       if(nData_B(iBlock)+nValue > size(Data_B(iBlock) % Array_I)) &
            call extend_array(iBlock, nValue)
       nData_B(iBlock) = nData_B(iBlock)+nValue
    end if

    i = iData_B(iBlock)
    do k = 1, nK; do j = 1, nJ
       Data_B(iBlock) % Array_I(i+1:i+nI) = Value_III(:,j,k)
       i = i+nI
    end do; end do
    iData_B(iBlock) = i

  end subroutine put_array3

  !===========================================================================

  subroutine put_array4(iBlock, nI, nJ, nK, nL, Value_IIII, DoAllowReplace)
    integer, intent(in) :: iBlock
    integer, intent(in) :: nI, nJ, nK, nL
    real,    intent(in) :: Value_IIII(nI, nJ, nK, nL)
    logical, intent(in), optional :: DoAllowReplace

    integer :: i, j, k, l, nValue
    character (len=*), parameter :: NameSub=NameMod//'::put_array4'
    !------------------------------------------------------------------------
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' called with nI, nJ, nK=',nI, nJ, nK, nL
    endif

    if(UseData_B(iBlock) .and. .not. present(DoAllowReplace)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock
       call stop_mpi('UseData_B=.true. in '//NameSub)
    end if

    nValue = nI*nJ*nK*nL

    if(iData_B(iBlock) + nValue > nData_B(iBlock) &
         .or. .not.present(DoAllowReplace))then
       if(nData_B(iBlock) < 0)call init_block(iBlock,1)
       if(nData_B(iBlock)+nValue > size(Data_B(iBlock) % Array_I)) &
            call extend_array(iBlock, nValue)
       nData_B(iBlock) = nData_B(iBlock)+nValue
    end if

    i = iData_B(iBlock)
    do l = 1, nL; do k = 1, nK; do j = 1, nJ
       Data_B(iBlock) % Array_I(i+1:i+nI) = Value_IIII(:,j,k,l)
       i = i+nI
    end do; end do; end do
    iData_B(iBlock) = i

  end subroutine put_array4

  !===========================================================================
  subroutine set_block_data(iBlock)
    integer, intent(in) :: iBlock
    character(len=*), parameter :: NameSub = NameMod//'::set_block_data'
    !------------------------------------------------------------------------
    UseData_B(iBlock) = .true.
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) write(*,*)NameSub,' called'
    endif
  end subroutine set_block_data

  !===========================================================================

  logical function use_block_data(iBlock)
    integer, intent(in) :: iBlock
    character (len=*), parameter :: NameSub = NameMod//'::use_block_data'
    !------------------------------------------------------------------------
    use_block_data = UseData_B(iBlock)
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' returning ',UseData_B(iBlock)
    endif
  end function use_block_data

  !===========================================================================

  subroutine get_point(iBlock, Value, DoNotAdvance)
    integer, intent(in) :: iBlock
    real,    intent(out):: Value
    logical, intent(in), optional :: DoNotAdvance

    character(len=*), parameter :: NameSub = NameMod//'::get_point'
    !------------------------------------------------------------------------
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' called'
    endif

    if(nData_B(iBlock) < 1) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock, &
            ' nData_B(iBlock) =',nData_B(iBlock)
       call stop_mpi('nData_B(iBlock) < 1 in '//NameSub)
    end if

    if(.not.UseData_B(iBlock)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock
       call stop_mpi('UseData_B=.false. in '//NameSub)
    end if

    ! wrap around
    if(iData_B(iBlock) >= nData_B(iBlock)) iData_B(iBlock) = 0

    ! Jump to next element and obtain value
    Value = Data_B(iBlock) % Array_I(iData_B(iBlock) + 1)

    if(.not.present(DoNotAdvance)) iData_B(iBlock) = iData_B(iBlock) + 1

    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' returning Value=',Value
    endif
  end subroutine get_point

  !===========================================================================

  subroutine get_array(iBlock, nValue, Value_I, DoNotAdvance)
    integer, intent(in) :: iBlock
    integer, intent(in) :: nValue
    real,    intent(out):: Value_I(nValue)
    logical, intent(in), optional :: DoNotAdvance

    integer :: i

    character(len=*), parameter :: NameSub = NameMod//'::get_array'
    !------------------------------------------------------------------------
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' called with nValue=',nValue
    endif
    ! wrap around
    if(iData_B(iBlock) >= nData_B(iBlock)) iData_B(iBlock) = 0

    i = iData_B(iBlock)
    if(i + nValue > nData_B(iBlock)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock, &
            ' nData_B =',nData_B(iBlock),&
            ' less than iData_B=',i,' + nValue=',nValue
       call stop_mpi('nData_B(iBlock) < iData_B+nValue in '//NameSub)
    end if

    if(.not.UseData_B(iBlock)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock
       call stop_mpi('UseData_B=.false. in '//NameSub)
    end if

    ! Read data
    Value_I = Data_B(iBlock) % Array_I(i+1:i+nValue)

    ! Adjust index
    if(.not.present(DoNotAdvance)) iData_B(iBlock) = i + nValue

  end subroutine get_array

  !===========================================================================

  subroutine get_array2(iBlock, nI, nJ, Value_II, DoNotAdvance)
    integer, intent(in) :: iBlock
    integer, intent(in) :: nI, nJ
    real,    intent(out):: Value_II(nI, nJ)
    logical, intent(in), optional :: DoNotAdvance

    integer :: i, j, nValue

    character(len=*), parameter :: NameSub = NameMod//'::get_array2'
    !------------------------------------------------------------------------
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' called with nI, nJ=',nI, nJ
    endif
    ! wrap around
    if(iData_B(iBlock) >= nData_B(iBlock)) iData_B(iBlock) = 0

    i = iData_B(iBlock)
    nValue = nI*nJ
    if(i + nValue > nData_B(iBlock)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock, &
            ' nData_B =',nData_B(iBlock),&
            ' less than iData_B=',i,' + nValue=',nValue
       call stop_mpi('nData_B(iBlock) < iData_B+nValue in '//NameSub)
    end if

    if(.not.UseData_B(iBlock)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock
       call stop_mpi('UseData_B=.false. in '//NameSub)
    end if

    ! Read data
    do j=1, nJ
       Value_II(:,j) = Data_B(iBlock) % Array_I(i+1:i+nI)
       i = i+nI
    end do

    ! Adjust index
    if(.not.present(DoNotAdvance)) iData_B(iBlock) = i

  end subroutine get_array2

  !===========================================================================

  subroutine get_array3(iBlock, nI, nJ, nK, Value_III, DoNotAdvance)
    integer, intent(in) :: iBlock
    integer, intent(in) :: nI, nJ, nK
    real,    intent(out):: Value_III(nI, nJ, nK)
    logical, intent(in), optional :: DoNotAdvance

    integer :: i, j, k, nValue

    character(len=*), parameter :: NameSub = NameMod//'::get_array3'
    !------------------------------------------------------------------------
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' called with nI, nJ, nK=',nI, nJ, nK
    endif
    ! wrap around
    if(iData_B(iBlock) >= nData_B(iBlock)) iData_B(iBlock) = 0

    i = iData_B(iBlock)
    nValue = nI*nJ*nK
    if(i + nValue > nData_B(iBlock)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock, &
            ' nData_B =',nData_B(iBlock),&
            ' less than iData_B=',i,' + nValue=',nValue
       call stop_mpi('nData_B(iBlock) < iData_B+nValue in '//NameSub)
    end if

    if(.not.UseData_B(iBlock)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock
       call stop_mpi('UseData_B=.false. in '//NameSub)
    end if

    ! Read data
    do k=1, nK; do j=1, nJ
       Value_III(:,j,k) = Data_B(iBlock) % Array_I(i+1:i+nI)
       i = i+nI
    end do; end do

    ! Adjust index
    if(.not.present(DoNotAdvance)) iData_B(iBlock) = i

  end subroutine get_array3
  
  !===========================================================================

  subroutine get_array4(iBlock, nI, nJ, nK, nL, Value_IIII, DoNotAdvance)
    integer, intent(in) :: iBlock
    integer, intent(in) :: nI, nJ, nK, nL
    real,    intent(out):: Value_IIII(nI, nJ, nK, nL)
    logical, intent(in), optional :: DoNotAdvance

    integer :: i, j, k, l, nValue

    character(len=*), parameter :: NameSub = NameMod//'::get_array4'
    !------------------------------------------------------------------------
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' called with nI, nJ, nK, nL=',nI, nJ, nK, nL
    endif
    ! wrap around
    if(iData_B(iBlock) >= nData_B(iBlock)) iData_B(iBlock) = 0

    i = iData_B(iBlock)
    nValue = nI*nJ*nK*nL
    if(i + nValue > nData_B(iBlock)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock, &
            ' nData_B =',nData_B(iBlock),&
            ' less than iData_B=',i,' + nValue=',nValue
       call stop_mpi('nData_B(iBlock) < iData_B+nValue in '//NameSub)
    end if

    if(.not.UseData_B(iBlock)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock
       call stop_mpi('UseData_B=.false. in '//NameSub)
    end if

    ! Read data
    do l=1, nL; do k=1, nK; do j=1, nJ
       Value_IIII(:,j,k,l) = Data_B(iBlock) % Array_I(i+1:i+nI)
       i = i+nI
    end do; end do; end do

    ! Adjust index
    if(.not.present(DoNotAdvance)) iData_B(iBlock) = i

  end subroutine get_array4
  
  !===========================================================================

  subroutine clean_block(iBlock)
    integer, intent(in) :: iBlock
    character (len=*), parameter :: NameSub = NameMod//'::clean_block'
    !------------------------------------------------------------------------
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) write(*,*)NameSub,' called'
    endif
    UseData_B(iBlock) = .false.
    if(nData_B(iBlock) < 0) RETURN
    deallocate(Data_B(iBlock) % Array_I)
    nData_B(iBlock) = -1
    iData_B(iBlock) = -1

  end subroutine clean_block

  !===========================================================================

  subroutine clean_all
    integer :: iBlock
    character (len=*), parameter :: NameSub = NameMod//'::clean_all'
    !------------------------------------------------------------------------
    do iBlock = 1, MaxBlock
       call clean_block(iBlock)
    end do
    if(DoDebug)then
       if(iProc==ProcTest) write(*,*)NameSub,' finished'
    endif
  end subroutine clean_all

  !===========================================================================

  integer function n_data(iBlock)
    integer, intent(in) :: iBlock
    character (len=*), parameter :: NameSub = NameMod//'::n_data'
    !------------------------------------------------------------------------
    n_data = nData_B(iBlock)
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' returning ',nData_B(iBlock)
    endif
  end function n_data

  !===========================================================================

  integer function max_data()
    character (len=*), parameter :: NameSub = NameMod//'::max_data'
    !------------------------------------------------------------------------
    max_data = maxval(nData_B)
    if(DoDebug)then
       if(iProc==ProcTest)write(*,*)NameSub,' returning ',maxval(nData_B)
    endif
  end function max_data

  !===========================================================================

  subroutine test_block_data

    real :: Value, Value_I(3), Value_II(2,3)
    integer :: nData, i
    character (len=*), parameter :: NameSub = NameMod//'::test_block_data'
    !--------------------------------------------------------------------------
    write(*,*)'Testing put_block_data and n_block_data()'

    ! Put 1 value
    call put_block_data(1,1.0)

    nData = n_block_data(1)
    if(nData /= 1)write(*,*)'n_block_data(1) failed, n=',nData,' should be 1'
    nData = n_block_data(2)
    if(nData /= -1)write(*,*)'n_block_data(2) failed, n=',nData,' should be -1'
    nData = n_block_data()
    if(nData /= 1)write(*,*)'n_block_data() failed, n=',nData,' should be 1'
    i = iData_B(1)
    if(i /= 1)write(*,*)'put_block_data failed, iData=',i,' should be 1'

    ! Put 3 values as an array (the optional argument should not matter now)
    call put_block_data(1,3,(/2.0, 3.0, 4.0/)) !, DoAllowReplace=.true.)

    nData = n_block_data(1)
    if(nData /= 4)write(*,*)'n_block_data(1) failed, n=',nData,' should be 4'
    nData = n_block_data()
    if(nData /= 4)write(*,*)'n_block_data() failed, n=',nData,' should be 4'
    i = iData_B(1)
    if(i /= 4)write(*,*)'put_block_data failed, iData=',i,' should be 4'
    i = size(Data_B(1) % Array_I)
    if(i /= 4)write(*,*)'put_block_data failed, size(Array)=',i,' should be 4'

    ! Put another 6 values as a 2d array
    call put_block_data(1,2,3,reshape((/1.0,2.0,3.0,4.0,5.0,6.0/), (/2,3/)))

    nData = n_block_data(1)
    if(nData /= 10)write(*,*)'n_block_data(1) failed, n=',nData,' should be 10'
    nData = n_block_data()
    if(nData /= 10)write(*,*)'n_block_data() failed, n=',nData,' should be 10'
    i = iData_B(1)
    if(i /= 10)write(*,*)'put_block_data failed, iData=',i,' should be 10'
    i = size(Data_B(1) % Array_I)
    if(i /= 10)write(*,*)'put_block_data failed, size(Array)=',i,&
         ' should be 10'

    ! This should fail with error message (get before set)
    ! call get_block_data(1,Value)

    write(*,*)'Testing set_block_data/use_block_data'
    if(use_block_data(1)) &
         write(*,*)'use_block_data(1)failed, =.true., should be .false.'
    call set_block_data(1)
    if(.not.use_block_data(1)) &
         write(*,*)'use_block_data(1)failed, =.false., should be .true.'

    ! This should fail with error message (put after set)
    ! call put_block_data(1,1.0)

    write(*,*)'Testing get_block_data'
    ! Get back first value
    call get_block_data(1,Value)

    if(Value /= 1.0) &
         write(*,*)'get_block_data failed, value =',Value,' should be 1'

    ! Get back next 3 values as an array
    call get_block_data(1,3,Value_I, DoNotAdvance=.true.)
    if(any(Value_I /= (/2.0, 3.0, 4.0/))) &
         write(*,*)'get_block_data do not advance failed, value_I=',&
         Value_I,' should be 2,3,4'

    ! Replace next 3 values as an array
    call put_block_data(1,3,(/2.1, 3.1, 4.1/), DoAllowReplace=.true.)

    ! Get back next 6 values as a 2D array
    call get_block_data(1,2,3,Value_II)
    if(any(Value_II /= reshape( (/1.0, 2.0, 3.0, 4.0, 5.0, 6.0/), (/2,3/)))) &
         write(*,*)'get_block_data failed, value_II=',Value_II, &
         ' should be ((2,3),(4,5),(6,7))'

    write(*,*)'Testing clean_block_data'
    ! Clean storage for block 1
    call clean_block_data(1)

    nData = n_block_data(1)
    if(nData /= -1) &
         write(*,*)'clena_block_daata(1) failed, n=',nData,' should be -1'

    ! Clean storage for all blocks
    call clean_block_data
    nData = n_block_data()
    if(nData /= -1)write(*,*) &
         'clean_block_data failed, n=',nData,' should be -1'

  end subroutine test_block_data

  !===========================================================================
  subroutine write_block_restart_files(NameRestartOutDir, UseRestartOutSeries)

    use ModMain, ONLY: nBlock, Unused_B
    use ModIO,   ONLY: Unit_Tmp

    character(len=*), intent(in) :: NameRestartOutDir
    logical,          intent(in) :: UseRestartOutSeries

    integer            :: iBlock
    character(len=100) :: NameBlockFile

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='write_block_restart_files'
    !--------------------------------------------------------------------
    if(iProc==PROCtest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    do iBlock=1,nBlock
       if (Unused_B(iBlock)) CYCLE

       call get_block_restart_namefile(iBlock, &
            NameRestartOutDir, UseRestartOutSeries, NameBlockFile)

       ! Skip blocks which do not have block data
       if (.not. use_block_data(iBlock)) then
          write(*,*) NameSub, ': not use_block_data, skip ', NameBlockFile
          CYCLE
       end if

       open(unit_tmp, file=NameBlockFile, status="replace", form='UNFORMATTED')
       write(Unit_tmp) nData_B(iBlock)
       write(Unit_tmp) Data_B(iBlock) % Array_I(1:nData_B(iBlock))
       close(unit_tmp)
    end do

    if(DoTestMe)then
       write(*,*)NameSub,': iProc, iBlock  =', iProc, BLKtest
       write(*,*)NameSub,': use_block_data =', use_block_data(BLKtest)
       write(*,*)NameSub,': nData_B(iBlock)=', nData_B(BLKtest)
       write(*,*)NameSub,': Data_B(iBlock)%Array_I(1:5) = ', &
            Data_B(BLKtest)%Array_I(1:5)
       write(*,*)NameSub,' finished'
    end if

  end subroutine write_block_restart_files

  !===========================================================================
  subroutine read_block_restart_files(NameRestartInDir, UseRestartInSeries)

    use ModMain, ONLY: nBlock, Unused_B
    use ModIO,   ONLY: Unit_Tmp

    character(len=*), intent(in) :: NameRestartInDir
    logical,          intent(in) :: UseRestartInSeries

    integer  :: iBlock, iError

    integer :: nData
    real    :: DataTmp_I(1:MaxBlockData)

    character(len=100) :: NameBlockFile

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='read_block_restart_files'
    !--------------------------------------------------------------------
    if(iProc==PROCtest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    do iBlock=1,nBlock
       if (Unused_B(iBlock)) CYCLE

       call get_block_restart_namefile(iBlock, &
            NameRestartInDir, UseRestartInSeries, NameBlockFile)

       open(unit_tmp, file=NameBlockFile, status='old', form='UNFORMATTED',&
            iostat = iError)

       ! Missing block data files (should be blocks without any block data)
       if(iError /= 0) then
          write(*,*) NameSub, ': could not open ', NameBlockFile
          CYCLE
       end if

       ! Read the number of data elements
       read(unit_tmp, iostat = iError) nData
       if(iError /= 0) call stop_mpi(NameSub// &
            ' could not read data from '//trim(NameBlockFile))

       ! Read the data array
       read(unit_tmp, iostat = iError) DataTmp_I(1:nData)
       if(iError /= 0) call stop_mpi(NameSub// &
            ' could not read data from '//trim(NameBlockFile))

       ! Empty the block storage before putting the block data
       if(use_block_data(iBlock)) call clean_block_data(iBlock)

       call put_block_data(iBlock, nData, DataTmp_I(1:nData))
       call set_block_data(iBlock)

       close(unit_tmp)
    end do

    if(DoTestMe)then
       write(*,*)NameSub,': iProc, iBlock  =', iProc, BLKtest
       write(*,*)NameSub,': use_block_data =', use_block_data(BLKtest)
       write(*,*)NameSub,': nData_B(iBlock)=', nData_B(BLKtest)
       write(*,*)NameSub,': Data_B(iBlock)%Array_I(1:5) = ', &
            Data_B(BLKtest)%Array_I(1:5)
       write(*,*)NameSub,' finished'
    end if

  end subroutine read_block_restart_files

  !============================================================================
  subroutine get_block_restart_namefile(iBlock, &
       NameRestartDir, UseRestartSeries, NameBlockFile)

    use BATL_lib, ONLY: iMortonNode_A, iNode_B
    use ModMain,  ONLY: iteration_number

    integer, intent(in) :: iBlock
    logical, intent(in) :: UseRestartSeries
    character(len=*),  intent(in)  :: NameRestartDir
    character(len=100),intent(out) :: NameBlockFile

    integer   :: iBlockRestart

    character :: StringDigit
    character(len=*), parameter :: NameStart        = "blockdata_Blk"
    character(len=*), parameter :: StringRestartExt = ".rst"

    character(len=*), parameter :: NameSub='get_block_restart_namefile'
    !--------------------------------------------------------------------

    iBlockRestart = iMortonNode_A(iNode_B(iBlock))
    write(StringDigit,'(i1)') max(5,int(1+alog10(real(iBlockRestart))))

    if (UseRestartSeries) then
       write(NameBlockFile, &
            '(a,i8.8,a,i'//StringDigit//'.'//StringDigit//',a)') &
            trim(NameRestartDir)//'n', iteration_number,         &
            '_'//NameStart,iBlockRestart,StringRestartExt
    else
       write(NameBlockFile,'(a,i'//StringDigit//'.'//StringDigit//',a)') &
            trim(NameRestartDir)//NameStart,iBlockRestart,StringRestartExt
    end if

  end subroutine get_block_restart_namefile

end module ModBlockData
