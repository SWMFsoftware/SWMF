!BOP
!MODULE: CON_line_extract - extract field and stream lines in parallel
!INTERFACE:
module CON_line_extract

  !DESCRIPTION:
  ! The ray position and the MHD state along the rays (ray line) 
  ! can be collected into an array and sorted by the length coordinate. 
  ! This class provides the infrastructure for collecting, 
  ! sorting and providing the data extracted along multiple ray lines.

  !USES:
  use ModMpi

  implicit none

  save    ! save all variables

  private ! except

  !PUBLIC MEMBER FUNCTIONS:
  public :: line_init         ! Initialize storage for ray lines
  public :: line_clean        ! Clean up storage for ray lines
  public :: line_put          ! Store state of a single point along a ray line
  public :: line_collect      ! Collect all ray line data onto 1 processor
  public :: line_get          ! Get all (sorted) ray line data from 1 processor
  public :: line_test         ! Unit tester

  !REVISION HISTORY:
  ! 09May04 - Gabor Toth <gtoth@umich.edu> - initial prototype/prolog/code
  !EOP

  ! Private constants
  character(len=*),  parameter :: NameMod='CON_line_extract'

  integer       :: iProc=0       ! Processor rank in MPI_COMM_WORLD
  integer       :: nVar=0        ! Number of variables for each point
  integer       :: nPoint=0      ! Actual number of points
  integer       :: MaxPoint=0    ! Allocated number of points
  real, pointer :: State_VI(:,:) ! State at each point along all rays

contains

  !BOP ========================================================================
  !IROUTINE: line_init - initialize storage for ray lines
  !INTERFACE:
  subroutine line_init(nVarIn)

    !INPUT ARGUMENTS:
    integer, intent(in) :: nVarIn   ! Number of variables to store

    !DESCRIPTION:
    ! Initialize the ray line storage.
    !EOP

    integer :: iError
    character (len=*), parameter :: NameSub = NameMod//'::line_init'
    !-------------------------------------------------------------------------
    
    if(nVar == nVarIn) RETURN          ! nothing to do if the same comm
    if(nVar >  0     ) call line_clean ! clean up previous allocation

    ! Set number of variables for each point
    nVar = nVarIn
    ! Set processor rank
    call MPI_COMM_RANK(MPI_COMM_WORLD, iProc,  iError)
    ! Initialize storaga
    nullify(State_VI)

  end subroutine line_init

  !BOP ========================================================================
  !IROUTINE: line_clean - clean storage for line data
  !INTERFACE:
  subroutine line_clean

    !DESCRIPTION:
    ! Clean the ray line storage.
    !EOP

    character (len=*), parameter :: NameSub = NameMod//'::line_clean'
    !-------------------------------------------------------------------------
    
    if(nVar == 0) RETURN  ! nothing to do, already clean

    nVar     = 0
    nPoint   = 0
    MaxPoint = 0
    if(associated(State_VI)) deallocate(State_VI)

  end subroutine line_clean

  !BOP ========================================================================
  !IROUTINE: line_extend - allocate or extend storage for line data
  !INTERFACE:
  subroutine line_extend(MaxPointIn)

    integer, intent(in) :: MaxPointIn

    real, pointer :: OldState_VI(:,:)
    !------------------------------------------------------------------------
    if(.not.associated(State_VI))then
       allocate(State_VI(0:nVar, MaxPointIn)) ! allocate storage
       MaxPoint = MaxPointIn                  ! set buffer size
    else
       OldState_VI => State_VI                ! store old values
       allocate(State_VI(0:nVar, MaxPointIn)) ! allocate new storage
       State_VI(:,1:nPoint) = &
            OldState_VI(:,1:nPoint)           ! copy old values
       deallocate(OldState_VI)                ! free old storage
       MaxPoint = MaxPointIn                  ! change buffer size
    end if
       
  end subroutine line_extend

  !BOP ========================================================================
  !IROUTINE: line_put - store state for a point along a ray line
  !INTERFACE:
  subroutine line_put(iLine, nVarIn, Line_V)

    !INPUT ARGUMENTS:
    integer, intent(in) :: iLine
    integer, intent(in) :: nVarIn
    real,    intent(in) :: Line_V(nVarIn)

    !DESCRIPTION:
    ! Store the ray index iLine and the state Line_V(1:nVarIn) 
    ! for the line iLine.
    !EOP

    character (len=*), parameter :: NameSub = NameMod//'line_put'
    !------------------------------------------------------------------------
    if(nVarIn /= nVar)then
       write(*,*) NameSub,' ERROR: nVarIn, nVar=',nVarIn, nVar
       call CON_stop(NameSub//' ERROR: incorrect number of variables!')
    end if

    nPoint = nPoint + 1

    if(nPoint > MaxPoint) call line_extend(nPoint + 100)

    State_VI(0       , nPoint) = real(iLine)
    State_VI(1:nVarIn, nPoint) = Line_V

  end subroutine line_put

  !BOP ========================================================================
  !IROUTINE: line_collect - collect ray line data onto 1 processor
  !INTERFACE:
  subroutine line_collect(iProc0From, nProcFrom, DiProcFrom, iProcTo)

    !INPUT ARGUMENTS:
    integer, intent(in) :: iProc0From ! Root PE of sending group
    integer, intent(in) :: nProcFrom  ! Number of sending PE-s
    integer, intent(in) :: DiProcFrom ! Stride of sending PE-s
    integer, intent(in) :: iProcTo    ! Receiving PE rank

    !DESCRIPTION:
    ! Collect Line from the sending processors to the receiving processor.
    ! The sending arrays are emptied, the receiving array is extended.
    ! The receiving PE can be part of the sending group.
    !EOP

    integer, allocatable :: nPoint_P(:)

    integer, parameter :: iTag = 1
    character(len=*), parameter :: NameSub = NameMod//'line_collect'

    integer :: i, iProcFrom, iPoint, Status_I(MPI_STATUS_SIZE), iError
    !-------------------------------------------------------------------------

    ! Collect number of points onto the receiving processors
    
    allocate(nPoint_P(nProcFrom))
    do i = 1, nProcFrom
       iProcFrom = iProc0From + (i-1)*DiProcFrom

       if(iProc == iProcFrom .and. iProc == iProcTo) CYCLE

       if(iProc == iProcFrom) call MPI_SEND(nPoint, &
            1, MPI_INTEGER, iProcTo, iTag, MPI_COMM_WORLD, iError)

       if(iProc == iProcTo) call MPI_RECV(nPoint_P(i), &
            1, MPI_INTEGER, iProcFrom, iTag, MPI_COMM_WORLD, Status_I, iError)
    end do

    ! Extend buffer on receiving processor
    if(iProc == iProcTo) call line_extend(sum(nPoint_P) + 100)

    ! Receive ray line data
    do i = 1, nProcFrom
       iProcFrom = iProc0From + (i-1)*DiProcFrom
       if(iProc == iProcFrom .and. iProc == iProcTo) CYCLE

       ! Send ray line data
       if(iProc == iProcFrom) then
          call MPI_SEND(State_VI(:,1:nPoint), &
               nPoint * (nVar+1), MPI_REAL, &
               iProcTo, iTag, MPI_COMM_WORLD, iError)
          deallocate(State_VI)
          nPoint = 0
       end if

       ! Receive ray line data
       if(iProc == iProcTo) then
          call MPI_RECV(State_VI(:,nPoint+1:nPoint + nPoint_P(i)), &
               nPoint_P(i) * (nVar+1), MPI_REAL, &
               iProcFrom, iTag, MPI_COMM_WORLD, Status_I, iError)
          nPoint = nPoint + nPoint_P(i)
       end if

    end do

    deallocate(nPoint_P)

  end subroutine line_collect

  !BOP ========================================================================
  !IROUTINE: line_get - obtain all the ray states
  !INTERFACE:
  subroutine line_get(nVarOut, nPointOut, Line_VI)

    !INPUT ARGUMENTS:
    integer, intent(inout)         :: nVarOut, nPointOut
    real,    intent(out), optional :: Line_VI(0:nVarOut, nPointOut)

    !DESCRIPTION:
    ! Obtain all the ray states stored. 
    !\begin{verbatim}
    ! First get the number of variables and points:
    !
    !     call line_get(nVarOut, nPointOut)
    !
    ! (Line_VI should not be present!), then 
    !
    !    allocate(Line_VI(0:nVarOut, nPointOut)
    ! 
    ! (the 0 index is needed for the line index) and then obtain the state
    !
    !    call line_get(nVarOut, nPointOut, Line_VI)
    !\end{verbatim}
    !EOP

    character (len=*), parameter :: NameSub = NameMod//'::line_get'

    !------------------------------------------------------------------------
    if(.not. present(Line_VI))then
       ! Provide size
       nVarOut   = nVar
       nPointOUt = nPoint
    else
       ! Check size
       if(nVar /= nVarOut .or. nPoint /= nPointOut)then
          write(*,*) NameSub,' ERROR: nVar  , nVarOut  =', &
               nVar, nVarOut
          write(*,*) NameSub,' ERROR: nPoint, nPointOut=', &
               nPoint, nPointOut
          call CON_stop(NameSub//' ERROR: incorrect array size')
       end if

       Line_VI = State_VI
    end if
       
  end subroutine line_get

  !BOP ========================================================================
  !IROUTINE: line_test - unit tester
  !INTERFACE:
  subroutine line_test

    !DESCRIPTION:
    ! Test the CON_line_extract module. This subroutine should be called from
    ! a small stand alone program.
    !EOP

    integer, parameter :: nVarTest = 4
    integer :: nProc, iError
    integer :: nVarOut, nPointOut, iPoint
    real, allocatable :: Result_VI(:,:)
    !------------------------------------------------------------------------

    call line_init(nVarTest)

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nProc, iError)

    write(*,*)iProc,': nProc = ',nProc

    write(*,*)iProc,': line_init done, nVar =',nVar

    call line_put(1, 4, (/0.0, 10.0+iProc, 20.0+iProc, 30.0+iProc/) )
    call line_put(1, 4, (/1.0, 40.0+iProc, 50.0+iProc, 60.0+iProc/) )
    call line_put(2, 4, (/0.0, 70.0+iProc, 80.0+iProc, 90.0+iProc/) )

    write(*,*)iProc,': line_put done, nPoint, MaxPoint=',nPoint, MaxPoint
    do iPoint = 1, nPoint
       write(*,'(i4,a,5f5.0)')iProc,': State_VI=',State_VI(:,iPoint)
    end do
 
    call line_collect(0,nProc,1,0)

    write(*,*)iProc,': line_collect done, nPoint, MaxPoint=',nPoint, MaxPoint

    if(iProc==0)then
       call line_get(nVarOut, nPointOut)
       write(*,*)'nVarOut, nPointOut=',nVarOut, nPointOut
       allocate(Result_VI(0:nVarOut,nPointOUt))
       call line_get(nVarOut, nPointOut, Result_VI)
       write(*,*)'iPoint Result_VI:'
       do iPoint = 1, nPointOut
          write(*,'(i4,5f5.0)')iPoint, Result_VI(:,iPoint)
       end do
    end if

    call line_clean

    write(*,*)iProc,': line_clean done, nPoint, MaxPoint, State_VI=',&
         nPoint, MaxPoint, associated(State_VI)
    
  end subroutine line_test

  !===========================================================================

end module CON_line_extract
