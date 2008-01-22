!BOP
!MODULE: CON_couple_ie_ps - couple IE and PS components
!
!DESCRIPTION:
! Couple IE and PS components
!
!INTERFACE:
module CON_couple_ie_ps

  !USES:
  use CON_coupler

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_ie_ps_init ! initialize both couplings
  public :: couple_ie_ps      ! couple IE to PS

  !REVISION HISTORY:
  ! 01/21/2008 D. De Zeeuw - initial version
  !
  !EOP

  ! Logicals to simplify message passing and execution
  logical :: IsInitialized = .false., UseMe=.true.

  ! Size of the 2D PS grid
  integer, save :: iSize, jSize
  real, allocatable :: MLTs(:,:), Lats(:,:)

contains

  !BOP =======================================================================
  !IROUTINE: couple_ie_ps_init - initialize IE-PS coupling
  !INTERFACE:
  subroutine couple_ie_ps_init
    !DESCRIPTION:
    ! Set and store PS grid.
    !EOP

    ! MPI status variable
    integer :: iStatus_I(MPI_STATUS_SIZE)

    ! General error code
    integer :: iError

    ! Message size
    integer :: nSize

    !------------------------------------------------------------------------
    if(IsInitialized) RETURN
    IsInitialized = .true.

    UseMe = is_proc(IE_) .or. is_proc(PS_)
    if(.not.UseMe) RETURN

    !\
    ! Set Grid
    ! Assume that all PS processors know grid info,
    !   but it must be passed to IE processors
    !/

    ! Get size
    if(is_proc(PS_)) call PS_get_grid_size(iSize,jSize)

    ! Transfer variables from IE proc0 to PS proc0
    if(i_proc0(IE_) /= i_proc0(PS_))then
       if(is_proc0(IE_)) &
            call MPI_send(iSize,1,MPI_INTEGER,i_proc0(IE_),&
            1,i_comm(),iError)
       if(is_proc(PS_)) &
            call MPI_recv(iSize,1,MPI_INTEGER,i_proc0(PS_),&
            1,i_comm(),iStatus_I,iError)

       if(is_proc0(IE_)) &
            call MPI_send(jSize,1,MPI_INTEGER,i_proc0(IE_),&
            1,i_comm(),iError)
       if(is_proc(PS_)) &
            call MPI_recv(jSize,1,MPI_INTEGER,i_proc0(PS_),&
            1,i_comm(),iStatus_I,iError)
    end if

    ! Allocate grid
    allocate(MLTs(iSize,jSize))
    allocate(Lats(iSize,jSize))

    ! Get grid
    if(is_proc(PS_)) call PS_get_grid(MLTs,Lats)

    ! Transfer variables from IE proc0 to PS proc0
    if(i_proc0(IE_) /= i_proc0(PS_))then
       nSize=iSize*jSize

       if(is_proc0(IE_)) &
            call MPI_send(MLTs,nSize,MPI_REAL,i_proc0(IE_),&
            1,i_comm(),iError)
       if(is_proc(PS_)) &
            call MPI_recv(MLTs,nSize,MPI_REAL,i_proc0(PS_),&
            1,i_comm(),iStatus_I,iError)

       if(is_proc0(IE_)) &
            call MPI_send(Lats,nSize,MPI_REAL,i_proc0(IE_),&
            1,i_comm(),iError)
       if(is_proc(PS_)) &
            call MPI_recv(Lats,nSize,MPI_REAL,i_proc0(PS_),&
            1,i_comm(),iStatus_I,iError)
    end if

    ! Set the coupler grids
    call IE_setnMlts('PS',iSize,iError)
    call IE_setnLats('PS',jSize,iError)
    call IE_setgrid('PS',MLTs,Lats,iError)

    ! Deallocate varibles
    deallocate(MLTs)
    deallocate(Lats)

  end subroutine couple_ie_ps_init

  !BOP =======================================================================
  !IROUTINE: couple_ie_ps - couple components, IE to PS
  !INTERFACE:
  subroutine couple_ie_ps(tSimulation)
    
    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:
    !    Ionosphere    (IE) source
    !    Plasamasphere (PS) target
    !
    ! Send field line volumes, average density and pressure and
    ! geometrical information.
    !EOP

    !\
    ! General coupling variables
    !/

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_ie_ps'

    logical :: DoTest, DoTestMe
    integer :: iProcWorld
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    iProcWorld = i_proc()
    call couple_mpi

    if(DoTest)write(*,*)NameSub,': finished iProc=',iProcWorld

  contains

    !==========================================================================
    subroutine couple_mpi

      character (len=*), parameter :: NameSubSub=NameSub//'.couple_mpi'

      ! Number of variables to pass
      integer, parameter :: nVarPS=6

      ! Buffer for the variables on the 2D PS grid
      real, dimension(:,:,:), allocatable :: Buffer_IIV

      ! MPI related variables

      ! MPI status variable
      integer :: iStatus_I(MPI_STATUS_SIZE)

      ! General error code
      integer :: iError

      ! Message size
      integer :: nSize

      !------------------------------------------------------------------------

      ! After everything is initialized exclude PEs which are not involved
      if(.not.UseMe) RETURN

      if(DoTest)write(*,*)NameSubSub,' starting, iProc=',iProcWorld
      if(DoTest)write(*,*)NameSubSub,', iProc, i_proc0(IE_), i_proc0(PS_)=', &
           iProcWorld,i_proc0(IE_),i_proc0(PS_)

      !\
      ! Allocate buffers in both IE and PS
      !/
      allocate(Buffer_IIV(iSize,jSize,nVarPS), stat=iError)
      call check_allocate(iError,NameSubSub//": Buffer_IIV")

      if(DoTest)write(*,*)NameSubSub,', variables allocated, iProc:',iProcWorld

      !\
      ! Get variables from IE
      !/
      if(is_proc(IE_)) call IE_get_for_ps(Buffer_IIV,iSize,jSize,nVarPS)

      !\
      ! Transfer variables from IE proc0 to PS proc0
      !/
      if(i_proc0(IE_) /= i_proc0(PS_))then
         nSize = iSize*jSize*nVarPS
         if(is_proc0(IE_)) &
              call MPI_send(Buffer_IIV,nSize,MPI_REAL,i_proc0(IE_),&
              1,i_comm(),iError)
         if(is_proc(PS_)) &
              call MPI_recv(Buffer_IIV,nSize,MPI_REAL,i_proc0(PS_),&
              1,i_comm(),iStatus_I,iError)
      end if

      ! Broadcast variables inside PS
      if(n_proc(PS_)>1 .and. is_proc(PS_)) &
           call MPI_bcast(Buffer_IIV,nSize,MPI_REAL,0,i_comm(PS_),iError)

      !\
      ! Put variables into PS
      !/
      if(is_proc(PS_)) call PS_put_from_ie(Buffer_IIV,iSize,jSize,nVarPS)

      !\
      ! Deallocate buffer to save memory
      !/
      deallocate(Buffer_IIV)

      if(DoTest)write(*,*)NameSubSub,', variables deallocated',&
           ', iProc:',iProcWorld

      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',iProcWorld

    end subroutine couple_mpi

  end subroutine couple_ie_ps

end module CON_couple_ie_ps
