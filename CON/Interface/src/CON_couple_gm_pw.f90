!^CMP COPYRIGHT UM
!^CMP FILE IE
!^CMP FILE PW

!BOP
!MODULE: CON_couple_gm_pw - couple PW and GM components
!
!DESCRIPTION:
! Couple PW and GM components both ways.
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
  !IROUTINE: couple_ie_pw - couple IE component to PW component
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

      integer :: iBlock     ! 1 for northern and 2 for southern hemisphere
      integer :: iVar       ! 1) CoLat, 2)Longitude, 3)Density1, 4)Density2
                            ! 5) Density3, 6)Velocity1, 7)Velocity2, 
                            ! 8) Velocity3
      integer :: iProcFrom  ! PE number sending the potential for current block
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
              call MPI_send(Buffer_VI,nSize,MPI_REAL,iProc0Pw,&
              1,i_comm(),iError)
         if(is_proc0(GM_)) &
              call MPI_recv(Buffer_VI,nSize,MPI_REAL,iProc0Gm,&
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

         !call GM_put_from_pw(Buffer_VI, nTotalLine, nVar, &
         !     NameVar_V, iBlock)
         !if(DoTest) &
         !     write(*,*)NameSubSub//' iProc, Buffer(1,1)=',&
         !     iProcWorld,Buffer_VI(1,1,:)
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
  subroutine couple_gm_pw

    call CON_stop('couple_gm_pw has not been implented yet')

  end subroutine couple_gm_pw

end module CON_couple_gm_pw

