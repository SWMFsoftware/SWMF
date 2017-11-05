module BATL_test

  ! Provides basic test functionality

  use ModReadParam, ONLY: lStringLine, read_var
  use BATL_size, ONLY: nDim, MaxDim
  use BATL_geometry, ONLY: x_, y_, z_
  use BATL_grid, ONLY: find_grid_block

  implicit none

  private ! except

  public:: read_test_param   ! read parameters for testing
  public:: find_test_cell    ! find test cell
  public:: do_test           ! logical function: test or not

  character(lStringLine), public:: StringTest = ' ' ! space separated list

  integer, public:: iTest=1, jTest=1, kTest=1     ! 1st test cell
  integer, public:: iTest2, jTest2, kTest2  ! 2nd test cell
  integer, public:: iBlockTest = 1          ! 1st test block
  integer, public:: iBlockTest2 = 1         ! 2nd test block
  integer, public:: iProcTest = 0           ! 1st test processor
  integer, public:: iProcTest2 = 0          ! 2nd test processor
  real,    public:: XyzTest_D(MaxDim) = 0.0 ! 1st test position
  real,    public:: xTest=0.0, yTest=0.0, zTest=0.0
  real,    public:: XyzTest2_D(MaxDim)= 0.0 ! 2nd test position
  real,    public:: xTest2=0.0, yTest2=0.0, zTest2=0.0

contains

  subroutine read_test_param(NameCommand)

    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'read_test_param'
    !---------------------------------------------------------
    select case(NameCommand)
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
  logical function do_test(NameSub, iProcIn, iBlockIn, iIn, jIn, kIn)

    character(len=*), intent(in):: NameSub
    integer, intent(in), optional:: iBlockIn
    integer, intent(in), optional:: iProcIn
    integer, intent(in), optional:: iIn, jIn, kIn
    !------------------------------------------------------------------------
    do_test = .false.

    if(present(iProcIn))then
       if(iProcIn /= iProcTest) RETURN
    end if
    if(present(iBlockIn))then
       if(iBlockIn /= iBlockTest) RETURN
    end if
    if(present(iIn))then
       if(iIn /= iTest) RETURN
    end if
    if(nDim > 1 .and. present(jIn))then
       if(jIn /= jTest) RETURN
    end if
    if(nDim > 2 .and. present(kIn))then
       if(kIn /= kTest) RETURN
    end if

    do_test = index(' '//StringTest//' ', ' '//NameSub//' ') > 0

  end function do_test


end module BATL_test
