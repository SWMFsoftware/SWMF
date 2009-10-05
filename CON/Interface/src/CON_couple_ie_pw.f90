!^CMP COPYRIGHT UM
!^CMP FILE IE
!^CMP FILE PW

!BOP
!MODULE: CON_couple_ie_pw - couple IE and PW components
!
!DESCRIPTION:
! Couple IE and PW components both ways.
!
!INTERFACE:
module CON_couple_ie_pw

  !USES:
  use CON_coupler

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_ie_pw_init ! initialize coupling
  public :: couple_ie_pw      ! couple IE to PW

  !REVISION HISTORY:
  ! 11/17/2006 A.Glocer and G.Toth - initial version 
  !EOP

  ! Communicator and logicals to simplify message passing and execution
  integer, save :: iCommIePw, iProc0Pw
  logical :: UseMe=.true., IsInitialized = .false.

  ! Size of the 2D spherical structured IE grid
  integer, save :: iSize, jSize, nCells_D(2)

contains

  !BOP =======================================================================
  !IROUTINE: couple_ie_pw_init - initialize IE-PW couplings
  !INTERFACE:
  subroutine couple_ie_pw_init

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s so that
    ! a union group can be formed. The IE grid size is also stored.
    !EOP

    if(IsInitialized) RETURN
    IsInitialized = .true.

    iProc0Pw = i_proc0(PW_)
    call set_router_comm(IE_,PW_,iCommIePw,UseMe,iProc0Pw)

    ! This works for a NODE BASED regular IE grid only
    iSize = Grid_C(IE_) % nCoord_D(1)
    jSize = Grid_C(IE_) % nCoord_D(2)

!    nCells_D = ncells_decomposition_d(IE_) + 1
!    iSize=nCells_D(1); jSize=nCells_D(2)

  end subroutine couple_ie_pw_init

  !BOP =======================================================================
  !IROUTINE: couple_ie_pw - couple IE component to PW component
  !INTERFACE:
  subroutine couple_ie_pw(tSimulation)

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Ionosphere Electrodynamics (IE)  source\\
    !    Polar Wind (PW) target
    !
    ! Send electrostatic potential and field aligned current from IE to PW.
    !EOP

    !\
    ! General coupling variables
    !/

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_ie_pw'

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

      character (len=*), parameter, dimension(2) :: &
           NameHem_B=(/'North','South'/)

      integer, parameter :: nVar = 4
      integer, parameter :: South_ = 2, North_ = 1

      character (len=*), parameter, dimension(nVar) :: &
           NameVar_V=(/'Pot','Jr ','Ave','Tot'/)

      ! Buffer for the potential on the 2D IE grid
      real, dimension(:,:,:), allocatable :: Buffer_IIV

      ! MPI related variables

      ! MPI status variable
      integer :: iStatus_I(MPI_STATUS_SIZE)

      ! General error code
      integer :: iError

      ! Message size
      integer :: nSize

      integer :: iBlock     ! 1 for northern and 2 for southern hemisphere
      integer :: iVar       ! 1 for Pot, 2 for AveE, 3 for EFlux
      integer :: iProcFrom  ! PE number sending the potential for current block
      !------------------------------------------------------------------------

      ! After everything is initialized exclude PEs which are not involved
      if(.not.UseMe) RETURN

      if(DoTest)write(*,*)NameSubSub,' starting, iProc=',iProcWorld
      if(DoTest)write(*,*)NameSubSub,', iProc, PWi_iProc0, IEi_iProc0=', &
           iProcWorld,i_proc0(PW_),i_proc0(IE_)

      !\
      ! Allocate buffers both in PW and IE
      !/
      allocate(Buffer_IIV(iSize,jSize,nVar), stat=iError)
      call check_allocate(iError,NameSubSub//": Buffer_IIV")

      if(DoTest)write(*,*)NameSubSub,', variables allocated',&
           ', iProc:',iProcWorld

!      do iBlock = South_, North_

      iBlock = North_

         !\
         ! Get potential from IE
         !/

         if(is_proc(IE_))  &
              call IE_get_for_pw(Buffer_IIV, iSize, jSize, &
              nVar, NameVar_V, NameHem_B(iBlock), tSimulation)

         !\
         ! Transfer variables from IE to PW
         !/ 

         iProcFrom = i_proc0(IE_)

         nSize = iSize*jSize*nVar

         if(iProcFrom /= i_proc0(PW_))then
            if(i_proc() == iProcFrom) &
                 call MPI_send(Buffer_IIV,nSize,MPI_REAL,i_Proc0(PW_),&
                 1,i_comm(),iError)
            if(is_proc0(PW_)) &
                 call MPI_recv(Buffer_IIV,nSize,MPI_REAL,iProcFrom,&
                 1,i_comm(),iStatus_I,iError)
         end if

         ! Broadcast variables inside PW
         if(n_proc(PW_)>1 .and. is_proc(PW_)) &
              call MPI_bcast(Buffer_IIV,nSize,MPI_REAL,0,i_comm(PW_),iError)

         if(DoTest)write(*,*)NameSubSub,', variables transferred',&
              ', iProc:',iProcWorld

         !\
         ! Put variables into PW
         !/
         if(is_proc(PW_))then
            call PW_put_from_ie(Buffer_IIV, iSize, jSize, nVar, &
                 NameVar_V, iBlock)
            if(DoTest) &
                 write(*,*)NameSubSub//' iProc, Buffer(1,1)=',&
                 iProcWorld,Buffer_IIV(1,1,:)
         end if

!      enddo
      !\
      ! Deallocate buffer to save memory
      !/
      deallocate(Buffer_IIV)

      if(DoTest)write(*,*)NameSubSub,', variables deallocated',&
           ', iProc:',iProcWorld

      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',iProcWorld

    end subroutine couple_mpi

  end subroutine couple_ie_pw

end module CON_couple_ie_pw

