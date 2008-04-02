!^CMP COPYRIGHT UM
!^CMP FILE IE
!^CMP FILE PW

!BOP
!MODULE: CON_couple_gm_pw - couple PW and GM components
!
!DESCRIPTION:
! Couple PW and GM components both ways.
! Note that PW -> GM has to be done before GM -> PW 
! so the PW 
!
!INTERFACE:
module CON_couple_gm_pw

  !USES:
  use CON_coupler

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_gm_pw_init ! initialize coupling
  public :: couple_gm_pw      ! couple GM to PW
  public :: couple_pw_gm      ! couple PW to GM

  !REVISION HISTORY:
  ! 01/18/2007 A.Glocer and G.Toth - initial version 
  !EOP

  ! Communicator and logicals to simplify message passing and execution
  integer, save :: iProc0Gm, iProc0Pw
  logical :: UseMe=.true., IsInitialized = .false.

  ! PW grid
  integer, save :: nTotalLine

contains

  !BOP =======================================================================
  !IROUTINE: couple_gm_pw_init - initialize PW-GM couplings
  !INTERFACE:
  subroutine couple_gm_pw_init

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s so that
    ! a union group can be formed. The PW grid size is also stored.
    !EOP

    if(IsInitialized) RETURN
    IsInitialized = .true.

    iProc0Gm = i_proc0(GM_)
    iProc0Pw = i_proc0(PW_)

    UseMe = is_proc(GM_) .or. is_proc(PW_)

    nTotalLine = Grid_C(PW_) % nCoord_D(2)

  end subroutine couple_gm_pw_init

  !BOP =======================================================================
  !IROUTINE: couple_pw_gm - couple PW component to GM component
  !INTERFACE:
  subroutine couple_pw_gm(tSimulation)

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Polar Wind (PW)  source\\
    !    Global Magnetosphere (GM) target
    !
    ! Send electrostatic potential and field aligned current from PW to IE.
    !EOP

    !\
    ! General coupling variables
    !/

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_pw_gm'

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

      ! Variable to pass is potential

      integer, parameter :: nVar = 8

      character (len=*), parameter, dimension(nVar) :: &
           NameVar_V=(/'CoLat    ','Longitude','Density1 ','Density2 ',&
           'Density3 ','Velocity1','Velocity2','Velocity3'/)

      ! Buffer for the potential on the 2D PW grid
      real, dimension(:,:), allocatable :: Buffer_VI

      ! MPI related variables

      ! MPI status variable
      integer :: iStatus_I(MPI_STATUS_SIZE)

      ! General error code
      integer :: iError

      ! Message size
      integer :: nSize

      integer :: iBlock
      integer :: iVar
      !------------------------------------------------------------------------

      ! After everything is initialized exclude PEs which are not involved
      if(.not.UseMe) RETURN

      if(DoTest)write(*,*)NameSubSub,' starting, iProc=',iProcWorld
      if(DoTest)write(*,*)NameSubSub,', iProc, GMi_iProc0, PWi_iProc0=', &
           iProcWorld,iProc0GM,iProc0Pw

      !\
      ! Allocate buffers both in GM and PW
      !/
      allocate(Buffer_VI(nVar, nTotalLine), stat=iError)
      call check_allocate(iError,NameSubSub//": Buffer_VI")

      if(DoTest)write(*,*)NameSubSub,', variables allocated',&
           ', iProc:',iProcWorld

      !\
      ! boundary density from PW
      !/

      if(is_proc(PW_))  &
           call PW_get_for_gm(Buffer_VI, nVar, nTotalLine, &
           NameVar_V, tSimulation)

      !\
      ! Transfer variables from PW to GM
      !/ 
      nSize = nTotalLine*nVar

      if(iProc0Gm /= iProc0Pw)then
         if(is_proc0(PW_)) &
              call MPI_send(Buffer_VI,nSize,MPI_REAL,iProc0Gm,&
              1,i_comm(),iError)
         if(is_proc0(GM_)) &
              call MPI_recv(Buffer_VI,nSize,MPI_REAL,iProc0Pw,&
              1,i_comm(),iStatus_I,iError)
      end if
      if(DoTest)write(*,*)NameSubSub,', variables transferred',&
           ', iProc:',iProcWorld

      !\
      ! Put variables into GM
      !/
      if(is_proc(GM_))then
         ! Broadcast variables inside GM
         if(n_proc(GM_)>1) &
              call MPI_bcast(Buffer_VI,nSize,MPI_REAL,0,i_comm(GM_),iError)

         call GM_put_from_pw(Buffer_VI, nVar, nTotalLine, NameVar_V)
         if(DoTest) &
              write(*,*)NameSubSub//' iProc, Buffer(1,1)=',&
              iProcWorld,Buffer_VI(1,1)
      end if

      !\
      ! Deallocate buffer to save memory
      !/
      deallocate(Buffer_VI)

      if(DoTest)write(*,*)NameSubSub,', variables deallocated',&
           ', iProc:',iProcWorld

      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',iProcWorld

    end subroutine couple_mpi

  end subroutine couple_pw_gm

  !==========================================================================
  subroutine couple_gm_pw(tSimulation)

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Global Magnetosphere (GM) source \\
    !    Polar Wind (PW)  target
    !
    ! Send pressure from GM to PW
    !EOP

    !\
    ! General coupling variables
    !/

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_gm_pw'

    ! Buffer for the coordinates of PW field lines
    real, dimension(:,:), allocatable:: Coord_DI

    ! Names for the coordinates
    character (len=*), parameter:: NameVar_D(2) = (/'CoLat    ','Longitude'/)

    ! Buffer for the pressure on the PW grid
    real, dimension(:), allocatable :: Buffer_I

    ! MPI related variables

    ! MPI status variable
    integer :: iStatus_I(MPI_STATUS_SIZE)

    ! General error code
    integer :: iError

    ! Message size
    integer :: nSize

    integer :: iBlock     
    integer :: iVar       
    integer :: iProcWorld ! global rank of this PE

    logical :: DoTest, DoTestMe
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    ! After everything is initialized exclude PEs which are not involved
    if(.not.UseMe) RETURN

    iProcWorld = i_proc()

    if(DoTest)write(*,*)NameSub,', iProc, GMi_iProc0, PWi_iProc0=', &
         iProcWorld,iProc0GM,iProc0Pw

    !\
    ! Allocate buffers both in GM and PW
    !/
    allocate(Coord_DI(2,nTotalLine), Buffer_I(nTotalLine))

    if(DoTest)write(*,*)NameSub,', variables allocated',&
         ', iProc:',iProcWorld

    !\
    ! coordinates from PW
    !/
    if(is_proc(PW_))  &
         call PW_get_for_gm(Coord_DI, 2, nTotalLine, NameVar_D, tSimulation)

    !\
    ! Transfer coordinates from PW to GM
    !/ 
    nSize = nTotalLine*2

    if(iProc0Gm /= iProc0Pw)then
       if(is_proc0(PW_)) &
            call MPI_send(Coord_DI,nSize,MPI_REAL,iProc0Pw,&
            1,i_comm(),iError)
       if(is_proc0(GM_)) &
            call MPI_recv(Coord_DI,nSize,MPI_REAL,iProc0Gm,&
            1,i_comm(),iStatus_I,iError)
    end if
    if(DoTest)write(*,*)NameSub,', PW coordinates transferred',&
         ', iProc:',iProcWorld

    !\
    ! Put PW field line coordinates into GM, and get pressure from GM
    !/
    if(is_proc(GM_))then
       ! Broadcast variables inside GM
       if(n_proc(GM_)>1) &
            call MPI_bcast(Coord_DI,nSize,MPI_REAL,0,i_comm(GM_),iError)

       call GM_put_from_pw(Coord_DI, 2, nTotalLine, NameVar_D)
       if(DoTest) &
            write(*,*)NameSub//' iProc,Coord_DI(:,1)=',iProcWorld,Coord_DI(:,1)
       
       call GM_get_for_pw(nTotalLine, Buffer_I)
    end if
    
    !\
    ! Send pressure from GM to PW
    !/
    if(iProc0Gm /= iProc0Pw)then
       if(is_proc0(GM_)) &
            call MPI_send(Buffer_I,nTotalLine,MPI_REAL,iProc0Gm,&
            1,i_comm(),iError)
       if(is_proc0(PW_)) &
            call MPI_recv(Buffer_I,nTotalLine,MPI_REAL,iProc0Pw,&
            1,i_comm(),iStatus_I,iError)
    end if

    if(DoTest)write(*,*)NameSub,', variables transferred',&
         ', iProc:',iProcWorld
    
    !\
    ! Put variables into PW
    !/
    if(is_proc(PW_))then
       ! Broadcast variables inside PW
       if(n_proc(PW_)>1) &
            call MPI_bcast(Buffer_I,nTotalLine,MPI_REAL,0,i_comm(PW_),iError)
       call PW_put_from_gm(nTotalLine,Buffer_I)
       if(DoTest) &
            write(*,*)NameSub//' iProc, Buffer_I(1)=',&
            iProcWorld,Buffer_I(1)
    end if
    
    !\
    ! Deallocate buffer to save memory
    !/
    deallocate(Coord_DI, Buffer_I)

    if(DoTest)write(*,*)NameSub,', variables deallocated',&
         ', iProc:',iProcWorld

    if(DoTest)write(*,*)NameSub,' finished, iProc=',iProcWorld

  end subroutine couple_gm_pw

end module CON_couple_gm_pw

