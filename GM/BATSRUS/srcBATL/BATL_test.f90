module BATL_test

  ! Provides basic test functionality

  use ModReadParam, ONLY: lStringLine, read_var
  use BATL_mpi,  ONLY: iProc, nProc, iComm
  use BATL_size, ONLY: nDim, MaxDim, nBlock, MaxBlock
  use BATL_geometry, ONLY: x_, y_, z_
  use BATL_grid, ONLY: find_grid_block

  implicit none

  private ! except

  public:: read_test_param ! read parameters for testing
  public:: find_test_cell  ! find test cell
  public:: test_start      ! start testing a subroutine/function
  public:: test_stop       ! stop testing a subroutine/function
  public:: test_cell       ! decide if a cell is to be tested

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

  ! Local variables

  ! Test cell based on position or index?
  logical:: UseTestXyz = .false., UseTest2Xyz = .false.

contains

  subroutine read_test_param(NameCommand)

    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'read_test_param'
    !---------------------------------------------------------
    select case(NameCommand)
    case("#VERBOSE")
       call              read_var('lVerbose',    lVerbose)
    case("#TESTXYZ")
       call              read_var('xTest',       XyzTest_D(x_))
       if(nDim > 1)call  read_var('yTest',       XyzTest_D(y_))
       if(nDim > 2)call  read_var('zTest',       XyzTest_D(z_))
       UseTestXyz = .true.
    case("#TESTIJK")
       call              read_var('iTest',       iTest)
       if(nDim > 1) call read_var('jTest',       jTest)
       if(nDim > 2) call read_var('kTest',       kTest)
       call              read_var('iBlockTest',  iBlockTest)
       if(nProc > 1)call read_var('iProcTest',   iProcTest)
       UseTestXyz = .false.
       iBlockTest = min(MaxBlock, max(1, iBlockTest))
       iProcTest  = min(nProc-1, max(0, iProcTest))
    case("#TEST2XYZ")
       call              read_var('xTest2',      XyzTest2_D(x_))
       if(nDim > 1) call read_var('yTest2',      XyzTest2_D(y_))
       if(nDim > 2) call read_var('zTest2',      XyzTest2_D(z_))
       UseTest2Xyz = .true.
    case("#TEST2IJK")
       call              read_var('iTest2',      iTest2)
       if(nDim > 1) call read_var('jTest2',      jTest2)
       if(nDim > 2) call read_var('kTest2',      kTest2)
       call              read_var('iBlockTest2', iBlockTest)
       if(nProc > 1)call read_var('iProcTest2',  iProcTest)
       UseTest2Xyz = .false.
       iBlockTest2 = min(MaxBlock, max(1, iBlockTest2))
       iProcTest2  = min(nProc-1, max(-1, iProcTest2))
    case("#TEST")
       call              read_var('StringTest',  StringTest)
    case default
       call CON_stop(NameSub//': unknown command='//NameCommand)
    end select
  end subroutine read_test_param
  !===========================================================================
  subroutine find_test_cell

    ! Find test cell based on the position

    use BATL_tree, ONLY: Unused_B
    use BATL_grid, ONLY: Xyz_DGB
    use ModMpi, ONLY: MPI_bcast, MPI_REAL

    integer:: iTest_D(MaxDim), iError

    character(len=*), parameter:: NameSub = "BATL_test::find_test_cell"
    !------------------------------------------------------------------------
    if(UseTestXyz)then

       ! Find grid cell based on position
       call find_grid_block(XyzTest_D, iProcTest, iBlockTest, iTest_D)
       if(iProcTest < 0)then
          if(iProc == 0) write(*,*) NameSub,' WARNING test point at ', &
               XyzTest_D,' is outside domain! Setting defaults.'
          iTest = 1; jTest = 1; kTest = 1; iBlockTest = 1; iProcTest = 0
       else
          iTest = iTest_D(1); jTest = iTest_D(2); kTest = iTest_D(3)
       end if

    elseif(iProcTest >= 0 .and. iProcTest < nProc) then
       
       ! Find location of grid cell based on indexes
       if(iProc == iProcTest)then
          if(iBlockTest > nBlock)then
             write(*,*) NameSub,' WARNING iBlockTest=', iBlockTest, &
                  ' is larger than nBlock=', nBlock
          else
             if(Unused_B(iBlockTest))then
                if(lVerbose > 0) write(*,*) &
                     NameSub,' WARNING: test cell is in an unused block'
             else
                XyzTest_D = Xyz_DGB(:,iTest,jTest,kTest,iBlockTest)
             end if
          end if
       end if

       ! Broadcast test cell position to other processors
       call MPI_Bcast(XyzTest_D, 3, MPI_REAL, iProcTest, iComm, iError)

       ! Set the scalars for convenience
       xTest = XyzTest_D(1); yTest = XyzTest_D(2); zTest = XyzTest_D(3)
    end if

    ! Deal with the second test cell

    if(UseTest2Xyz)then

       ! Find grid cell based on position
       call find_grid_block(XyzTest2_D, iProcTest2, iBlockTest2, iTest_D)
       if(iProcTest2 < 0)then
          if(iProc == 0) write(*,*) NameSub,' WARNING test point at ', &
               XyzTest_D,' is outside domain! Setting defaults.'
          iTest2 = 1; jTest2 = 1; kTest2 = 1; iBlockTest2 = 1; iProcTest2 = -1
       else
          iTest2 = iTest_D(1); jTest2 = iTest_D(2); kTest2 = iTest_D(3)
       end if

    elseif(iProcTest2 >= 0 .and. iProcTest2 < nProc) then
       
       ! Find location of grid cell based on indexes
       if(iProc == iProcTest2)then
          if(iBlockTest > nBlock)then
             write(*,*) NameSub,' WARNING iBlockTest2=', iBlockTest2, &
                  ' is larger than nBlock=', nBlock
          else
             if(Unused_B(iBlockTest2))then
                if(lVerbose > 0) write(*,*) &
                     NameSub,' WARNING: 2nd test cell is in an unused block'
             else
                XyzTest2_D = Xyz_DGB(:,iTest2,jTest2,kTest2,iBlockTest2)
             end if
          end if
       end if

       ! Broadcast test cell position to other processors
       call MPI_Bcast(XyzTest2_D, 3, MPI_REAL, iProcTest2, iComm, iError)

       ! Set the scalars for convenience
       xTest2 = XyzTest2_D(1); yTest2 = XyzTest2_D(2); zTest2 = XyzTest2_D(3)
    end if

  end subroutine find_test_cell
  !===========================================================================
  subroutine test_start(NameSub, DoTest, iBlock, DoTestAll)

    ! If optional block index iBlock is present, restrict all actions 
    ! to the test block(s) only.
    !
    ! Report this call on all processors if lVerbose == 100
    ! or DoTestAll is present and true.
    !
    ! Report on the test processor(s) if lVerbose == 10 or 
    ! NameSub matches StringTest and lVerbose /= 0.
    !
    ! In the latter case set the optional DoTestOut to true
    ! on the test processor or possibly on all processors

    character(len=*),  intent(in) :: NameSub   ! method to be tested
    logical,           intent(out):: DoTest    ! return true if testing is on

    integer, optional, intent(in) :: iBlock    ! block index
    logical, optional, intent(in) :: DoTestAll ! test on all processors

    logical:: DoWriteAll
    !------------------------------------------------------------------------

    ! Check block index if present
    if(present(iBlock))then
       if(  (iProc /= iProcTest  .or. iBlock /= iBlockTest) .and. &
            (iProc /= iProcTest2 .or. iBlock /= iBlockTest2))then
          DoTest = .false.
          RETURN
       end if
    end if

    DoTest = index(' '//StringTest//' ', ' '//NameSub//' ') > 0

    if(lVerbose == 0) RETURN

    DoWriteAll = lVerbose == 100
    if(present(DoTestAll)) DoWriteAll = DoTestAll

    if(DoWriteAll .or. ((lVerbose == 10 .or. DoTest) &
         .and. (iProc == iProcTest .or. iProc == iProcTest2)))then
       if(present(iBlock))then
          write(*,*) NameSub,' is starting for iProc, iBlock=', iProc, iBlock
       elseif(nProc > 1 .and. (DoWriteAll .or. iProcTest2 >= 0))then
          write(*,*) NameSub,' is starting on iProc=', iProc
       else
          write(*,*) NameSub,' is starting'
       end if
    end if

  end subroutine 
  !===========================================================================
  subroutine test_stop(NameSub, DoTest, iBlock)

    ! If optional block index iBlock is present, restrict all actions 
    ! to the test block(s) only.
    ! Write out a "finished" message if DoTest is true

    character(len=*),  intent(in):: NameSub
    logical,           intent(in):: DoTest  
    integer, optional, intent(in):: iBlock
    !-----------------------------------------------------------------------

    ! Check block index if present
    if(present(iBlock))then
       if(  (iProc /= iProcTest  .or. iBlock /= iBlockTest) .and. &
            (iProc /= iProcTest2 .or. iBlock /= iBlockTest2)) &
            RETURN
    end if

    if(DoTest .or. lVerbose > 1) then
       if(present(iBlock))then
          write(*,*) NameSub,' is finished for iProc, iBlock=', iProc, iBlock
       elseif(nProc > 1 .and. (lVerbose == 100 .or. iProcTest2 >= 0))then
          write(*,*) NameSub,' is finished on iProc=', iProc
       else
          write(*,*) NameSub,' is finished'
       end if
    end if

  end subroutine test_stop

  !===========================================================================
  subroutine test_cell(iBlock, i, j, k, DoTestCell)

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

  end subroutine test_cell
  !===========================================================================

end module BATL_test