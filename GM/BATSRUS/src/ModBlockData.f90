module ModBlockData

  use ModSize,   ONLY: MaxBlock
  use ModMain,   ONLY: ProcTest, BlkTest
  use ModProcMH, ONLY: iProc

  implicit none

  private ! except

  public put_block_data    ! store 1 or more values into storage
  interface put_block_data
     module procedure put_point, put_array
  end interface

  public set_block_data    ! indicate that all block data has been set

  public use_block_data    ! function to check if block data has been set

  public get_block_data    ! get 1 or more values from storage
  interface get_block_data
     module procedure get_point, get_array
  end interface

  public n_block_data      ! function for number of data values stored
  interface n_block_data
     module procedure n_data, max_data
  end interface

  public clean_block_data
  interface clean_block_data
     module procedure clean_block, clean_all
  end interface

  public test_block_data

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

  subroutine put_point(iBlock, Value)
    integer, intent(in) :: iBlock
    real,    intent(in) :: Value
    character (len=*), parameter :: NameSub=NameMod//'::put_point'
    !------------------------------------------------------------------------
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' called'
    endif

    if(UseData_B(iBlock)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock
       call CON_stop('UseData_B=.true. in '//NameSub)
    end if

    if(nData_B(iBlock) < 0)call init_block(iBlock,1)
    if(nData_B(iBlock)+1 > size(Data_B(iBlock) % Array_I)) &
         call extend_array(iBlock,1)
    nData_B(iBlock) = nData_B(iBlock)+1

    iData_B(iBlock) = iData_B(iBlock)+1
    Data_B(iBlock) % Array_I(iData_B(iBlock)) = Value

    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' returning Value=',Value
    endif

  end subroutine put_point

  !===========================================================================

  subroutine put_array(iBlock, nValue, Value_I)
    integer, intent(in) :: iBlock
    integer, intent(in) :: nValue
    real,    intent(in) :: Value_I(nValue)

    integer :: i
    character (len=*), parameter :: NameSub=NameMod//'::put_array'
    !------------------------------------------------------------------------
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' called with nValue=',nValue
    endif

    if(UseData_B(iBlock)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock
       call CON_stop('UseData_B=.true. in '//NameSub)
    end if

    if(nData_B(iBlock) < 0)call init_block(iBlock,1)
    if(nData_B(iBlock)+nValue > size(Data_B(iBlock) % Array_I)) &
         call extend_array(iBlock, nValue)
    nData_B(iBlock) = nData_B(iBlock)+nValue

    i = iData_B(iBlock)
    Data_B(iBlock) % Array_I(i+1:i+nValue) = Value_I

    iData_B(iBlock) = i+nValue

  end subroutine put_array

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

  subroutine get_point(iBlock, Value)
    integer, intent(in) :: iBlock
    real,    intent(out):: Value

    character(len=*), parameter :: NameSub = NameMod//'::get_point'
    !------------------------------------------------------------------------
    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' called'
    endif

    if(nData_B(iBlock) < 1) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock, &
            ' nData_B(iBlock) =',nData_B(iBlock)
       call CON_stop('nData_B(iBlock) < 1 in '//NameSub)
    end if

    if(.not.UseData_B(iBlock)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock
       call CON_stop('UseData_B=.false. in '//NameSub)
    end if

    ! wrap around
    if(iData_B(iBlock) >= nData_B(iBlock)) iData_B(iBlock) = 0

    ! Jump to next element and obtain value
    iData_B(iBlock) = iData_B(iBlock) + 1
    Value = Data_B(iBlock) % Array_I(iData_B(iBlock))

    if(DoDebug)then
       if(iProc==ProcTest .and. iBlock==BlkTest) &
            write(*,*)NameSub,' returning Value=',Value
    endif
  end subroutine get_point

  !===========================================================================

  subroutine get_array(iBlock, nValue, Value_I)
    integer, intent(in) :: iBlock
    integer, intent(in) :: nValue
    real,    intent(out):: Value_I(nValue)

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
       call CON_stop('nData_B(iBlock) < iData_B+nValue in '//NameSub)
    end if

    if(.not.UseData_B(iBlock)) then
       write(*,*)NameSub,' ERROR for iBlock=',iBlock
       call CON_stop('UseData_B=.false. in '//NameSub)
    end if

    ! Read data
    Value_I = Data_B(iBlock) % Array_I(i+1:i+nValue)

    ! Adjust index
    iData_B(iBlock) = i + nValue

  end subroutine get_array

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

    real :: Value, Value_I(3)
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

    ! Put 3 values as an array
    call put_block_data(1,3,(/2.0, 3.0, 4.0/))

    nData = n_block_data(1)
    if(nData /= 4)write(*,*)'n_block_data(1) failed, n=',nData,' should be 4'
    nData = n_block_data()
    if(nData /= 4)write(*,*)'n_block_data() failed, n=',nData,' should be 4'
    i = iData_B(1)
    if(i /= 4)write(*,*)'put_block_data failed, iData=',i,' should be 4'
    i = size(Data_B(1) % Array_I)
    if(i /= 4)write(*,*)'put_block_data failed, size(Array)=',i,' should be 4'

    ! Put another 3 values as an array
    call put_block_data(1,3,(/5.0, 6.0, 7.0/))

    nData = n_block_data(1)
    if(nData /= 7)write(*,*)'n_block_data(1) failed, n=',nData,' should be 7'
    nData = n_block_data()
    if(nData /= 7)write(*,*)'n_block_data() failed, n=',nData,' should be 7'
    i = iData_B(1)
    if(i /= 7)write(*,*)'put_block_data failed, iData=',i,' should be 7'
    i = size(Data_B(1) % Array_I)
    if(i /= 8)write(*,*)'put_block_data failed, size(Array)=',i,' should be 8'

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
    call get_block_data(1,3,Value_I)
    if(any(Value_I /= (/2.0, 3.0, 4.0/))) &
         write(*,*)'get_block_data failed, value_I=',Value_I,' should be 2,3,4'

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

end module ModBlockData
