! !  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
! !  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE IE
!^CMP FILE RB

!BOP
!MODULE: CON_couple_ie_rb - couple IE and RB components
!
!DESCRIPTION:
! Couple IE and RB components both ways.
!
!INTERFACE:
module CON_couple_ie_rb

  !USES:
  use CON_coupler

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_ie_rb_init ! initialize coupling
  public :: couple_ie_rb      ! couple IE to RB

  !REVISION HISTORY:
  ! 1/25/2007 A.Glocer and G.Toth - initial version 
  !EOP

  ! Communicator and logicals to simplify message passing and execution
  integer, save :: iCommIeRb, iProc0Rb
  logical :: UseMe=.true., IsInitialized = .false.

  ! Size of the 2D spherical structured IE grid
  integer, save :: iSize, jSize, nCells_D(2)

contains

  !BOP =======================================================================
  !IROUTINE: couple_ie_rb_init - initialize IE-RB couplings
  !INTERFACE:
  subroutine couple_ie_rb_init

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s so that
    ! a union group can be formed. The IE grid size is also stored.
    !EOP

    if(IsInitialized) RETURN
    IsInitialized = .true.

    iProc0Rb = i_proc0(RB_)
    call set_router_comm(IE_,RB_,iCommIeRb,UseMe,iProc0Rb)

    ! This works for a NODE BASED regular IE grid only
    nCells_D = ncells_decomposition_d(IE_) + 1
    iSize=nCells_D(1); jSize=nCells_D(2)

  end subroutine couple_ie_rb_init

  !BOP =======================================================================
  !IROUTINE: couple_ie_rb - couple IE component to RB component
  !INTERFACE:
  subroutine couple_ie_rb(tSimulation)

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Ionosphere Electrodynamics (IE)  source\\
    !    Radiation Belt (RB) target
    !
    ! Send electrostatic potential and field aligned current from IE to RB.
    !EOP

    !\
    ! General coupling variables
    !/

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_ie_rb'

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

      integer, parameter :: nVar = 1
      integer, parameter :: South_ = 1, North_ = 2

      character (len=*), parameter, dimension(nVar) :: &
           NameVar_V=(/'Pot'/)

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
      if(DoTest)write(*,*)NameSubSub,', iProc, RBi_iProc0, IEi_iProc0=', &
           iProcWorld,i_proc0(RB_),i_proc0(IE_)

      !\
      ! Allocate buffers both in RB and IE
      !/
      allocate(Buffer_IIV(iSize,jSize,nVar), stat=iError)
      call check_allocate(iError,NameSubSub//": Buffer_IIV")

      if(DoTest)write(*,*)NameSubSub,', variables allocated',&
           ', iProc:',iProcWorld

      do iBlock = South_, North_

         !\
         ! Get potential from IE
         !/

         if(is_proc(IE_))  &
              call IE_get_for_rb(Buffer_IIV, iSize, jSize, &
              nVar, NameVar_V, NameHem_B(iBlock), tSimulation)

         !\
         ! Transfer variables from IE to RB
         !/ 

         iProcFrom = pe_decomposition(IE_,iBlock)

         nSize = iSize*jSize*nVar

         if(iProcFrom /= i_proc0(RB_))then
            if(i_proc() == iProcFrom) &
                 call MPI_send(Buffer_IIV,nSize,MPI_REAL,i_Proc0(RB_),&
                 1,i_comm(),iError)
            if(is_proc0(RB_)) &
                 call MPI_recv(Buffer_IIV,nSize,MPI_REAL,iProcFrom,&
                 1,i_comm(),iStatus_I,iError)
         end if

         ! Broadcast variables inside RB
         if(n_proc(RB_)>1 .and. is_proc(RB_)) &
              call MPI_bcast(Buffer_IIV,nSize,MPI_REAL,0,i_comm(RB_),iError)

         if(DoTest)write(*,*)NameSubSub,', variables transferred',&
              ', iProc:',iProcWorld

         !\
         ! Put variables into RB
         !/
         if(is_proc(RB_))then
            call RB_put_from_ie(Buffer_IIV, iSize, jSize, nVar, &
                 NameVar_V, iBlock)
            if(DoTest) &
                 write(*,*)NameSubSub//' iProc, Buffer(1,1)=',&
                 iProcWorld,Buffer_IIV(1,1,:)
         end if

      enddo
      !\
      ! Deallocate buffer to save memory
      !/
      deallocate(Buffer_IIV)

      if(DoTest)write(*,*)NameSubSub,', variables deallocated',&
           ', iProc:',iProcWorld

      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',iProcWorld

    end subroutine couple_mpi

  end subroutine couple_ie_rb

end module CON_couple_ie_rb

