module BATL_test

  ! Provides basic test functionality

  use ModReadParam, ONLY: lStringLine, read_var
  use BATL_mpi,  ONLY: iProc
  use BATL_size, ONLY: nDim, MaxDim
  use BATL_geometry, ONLY: x_, y_, z_
  use BATL_grid, ONLY: find_grid_block

  implicit none

  private ! except

  public:: read_test_param   ! read parameters for testing
  public:: find_test_cell    ! find test cell
  public:: set_do_test       ! check NameSub for testing
  public:: set_do_test_proc  ! check NameSub and processor for testing
  public:: set_do_test_block ! check NameSub, processor and block for testing
  public:: set_do_test_cell  ! check processor, block and cell indexes to test

  character(lStringLine), public:: StringTest = ' '    ! space separated list
  integer, public:: iTest  = 1, jTest  = 1, kTest  = 1 ! 1st test cell
  integer, public:: iTest2 = 1, jTest2 = 1, kTest2 = 1 ! 2nd test cell
  integer, public:: iBlockTest = 1                     ! 1st test block
  integer, public:: iBlockTest2 = 1                    ! 2nd test block
  integer, public:: iProcTest = 0                      ! 1st test processor
  integer, public:: iProcTest2 = -1                    ! 2nd test processor
  real,    public:: XyzTest_D(MaxDim) = 0.0            ! 1st test position
  real,    public:: xTest = 0.0, yTest = 0.0, zTest = 0.0 
  real,    public:: XyzTest2_D(MaxDim)= 0.0            ! 2nd test position
  real,    public:: xTest2 = 0.0, yTest2 = 0.0, zTest2 = 0.0

  ! verbosity level:
  !   lVerbose=0:   no verbose output
  !   lVerbose=1:   minimal verbose output
  !   lVerbose=10:  verbose output on test processor
  !   lVerbose=100: verbose output on all processors
  integer, public:: lVerbose = 1                       

contains

  subroutine read_test_param(NameCommand)

    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'read_test_param'
    !---------------------------------------------------------
    select case(NameCommand)
    case("#VERBOSE")
       call             read_var('lVerbose', lVerbose)
    case("#TESTXYZ")
       call             read_var('xTest', XyzTest_D(x_))
       if(nDim > 1)call read_var('yTest', XyzTest_D(y_))
       if(nDim > 2)call read_var('zTest', XyzTest_D(z_))
    case("#TESTIJK")
       call read_var('iTest', iTest)
       call read_var('jTest', jTest)
       call read_var('kTest', kTest)
       call read_var('iBlockTest', iBlockTest)
       call read_var('iProcTest', iProcTest)
    case("#TEST")
       call read_var('StringTest', StringTest)
    case default
       call CON_stop(NameSub//': unknown command='//NameCommand)
    end select
  end subroutine read_test_param
  !===========================================================================
  subroutine find_test_cell

    ! Find test cell based on the position

    integer:: iTest_D(MaxDim)
    !------------------------------------------------------------------------
    call find_grid_block(XyzTest_D, iProcTest, iBlockTest, iTest_D)
    iTest = iTest_D(1); jTest = iTest_D(2); kTest = iTest_D(3)

  end subroutine find_test_cell
  !===========================================================================
  subroutine set_do_test(NameSub, DoTestOut)

    ! Report this call on all processors if lVerbose=100
    ! report on the test processor(s) if lVerbose=10 or 
    ! NameSub matches StringTest.
    ! In the latter case set the optional DoTestOut to true.

    character(len=*),  intent(in) :: NameSub
    logical, optional, intent(out):: DoTestOut

    logical:: DoTest
    !------------------------------------------------------------------------
    DoTest = index(' '//StringTest//' ', ' '//NameSub//' ') > 0

    if(present(DoTestOut)) DoTestOut = DoTest

    if(lVerbose == 100 .or. ((DoTest .or. lVerbose == 10) .and. &
         (iProc == iProcTest .or. iProc == iProcTest2))) &
         write(*,*) NameSub,' is starting: iProc=', iProc

  end subroutine set_do_test
  !===========================================================================
  subroutine set_do_test_proc(NameSub, DoTestOut)

    ! Report this call on all processors if lVerbose = 100.
    ! Report this call on the test processor(s) if lVerbose = 10 or
    ! NameSub matches StringTest. In the latter case set the optional
    ! argument DoTestOut=.true. otherwise it is .false.

    character(len=*),  intent(in) :: NameSub
    logical, optional, intent(out):: DoTestOut

    logical:: DoTest
    !------------------------------------------------------------------------
    if(lVerbose >= 100) write(*,*) NameSub,' is starting: iProc=', iProc

    if(iProc /= iProcTest .and. iProc /= iProcTest2)then
       if(present(DoTestOut)) DoTestOut = .false.
       RETURN
    end if

    DoTest = index(' '//StringTest//' ', ' '//NameSub//' ') > 0

    if(present(DoTestOut)) DoTestOut = DoTest

    if(DoTest .or. lVerbose == 10) &
         write(*,*) NameSub,' is starting on processor ',iProc

  end subroutine set_do_test_proc
  !===========================================================================
  subroutine set_do_test_block(NameSub, iBlock, DoTestOut)

    ! Report this call for the tested block(s) if lVerbose > 1 or 
    ! NameSub matches StringTest. In the latter case set the optional
    ! argument DoTestOut=.true. otherwise it is .false.

    character(len=*),  intent(in) :: NameSub
    integer,           intent(in) :: iBlock
    logical, optional, intent(out):: DoTestOut

    logical:: DoTest
    !------------------------------------------------------------------------
    if(  (iProc /= iProcTest  .or. iBlock /= iBlockTest) .and. &
         (iProc /= iProcTest2 .or. iBlock /= iBlockTest2))then
       if(present(DoTestOut)) DoTestOut = .false.
       RETURN
    end if

    DoTest = index(' '//StringTest//' ', ' '//NameSub//' ') > 0

    if(present(DoTestOut)) DoTestOut = DoTest

    if(DoTest .or. lVerbose > 1) &
         write(*,*) NameSub,' is starting: iProc, iBlock=', iProc, iBlock

  end subroutine set_do_test_block
  !===========================================================================
  subroutine set_do_test_cell(iBlock, i, j, k, DoTestCell)

    ! Set DoTestCell to true if the processor, block and cell indexes
    ! match one of the test cells.

    integer, intent(in) :: iBlock, i, j, k
    logical, intent(out):: DoTestCell
    !------------------------------------------------------------------------
    DoTestCell =                                                     &
         (iProc == iProcTest .and. iBlock == iBlockTest .and.        &
         i == iTest .and. j == jTest .and. k == kTest         ) .or. &
         (iProc == iProcTest2 .and. iBlock == iBlockTest2 .and.      &
         i == iTest2 .and. j == jTest2 .and. k == kTest2      )

  end subroutine set_do_test_cell
  !===========================================================================

end module BATL_test
