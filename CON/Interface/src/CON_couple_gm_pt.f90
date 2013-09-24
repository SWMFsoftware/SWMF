!^CMP COPYRIGHT UM
!^CMP FILE GM
!^CMP FILE PT

!BOP
!MODULE: CON_couple_gm_pt - couple GM and PT components
!
!DESCRIPTION:
! Couple GM and PT components both ways. 
!
!INTERFACE:
module CON_couple_gm_pt

  !USES:
  use CON_coupler

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_gm_pt_init ! initialize both couplings
  public :: couple_gm_pt      ! couple GM to PT

  !REVISION HISTORY:
  ! 09/24/2013 G.Toth <gtoth@umich.edu> - initial version
  !EOP

  ! Communicator and logicals to simplify message passing and execution
  logical       :: UseMe=.true., IsInitialized=.false.

  ! These are useful to store and cannot change
  integer, save:: iProcWorld, iCommWorld

contains

  !BOP =======================================================================
  !IROUTINE: couple_gm_pt_init - initialize GM-PT couplings
  !INTERFACE:
  subroutine couple_gm_pt_init

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='couple_gm_pt_init'

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s
    !EOP
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(IsInitialized) RETURN
    IsInitialized = .true.
    
    iCommWorld = i_comm()
    iProcWorld = i_proc()

    ! Set UseMe to .true. for the participating PE-s
    UseMe = is_proc(IE_) .or. is_proc(GM_)

  end subroutine couple_gm_pt_init

  !BOP =======================================================================
  !IROUTINE: couple_gm_pt - couple GM to PT
  !INTERFACE:
  subroutine couple_gm_pt(tSimulation)

    !INPUT ARGUMENT:
    real, intent(in) :: tSimulation

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Global Magnetosphere       (GM) source\\
    !    Ionosphere Electrodynamics (PT) target
    !
    ! Send information from GM to PT. 
    !EOP

    ! Number of dimensions and number of variables
    integer :: nDim, nVar

    ! List of variables to pass
    character(len=100):: NameVar = 'Bx By Bz'

    logical :: DoTest, DoTestMe

    ! Buffer for the field sent
    real, dimension(:,:), allocatable :: Buffer_VI

    ! MPI related variables

    ! MPI status variable
    integer :: iStatus_I(MPI_STATUS_SIZE)

    ! General error code
    integer :: iError

    ! Message size
    integer :: nSize

    integer :: iBlock, iProcTo

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_gm_pt'
    !-------------------------------------------------------------------------

    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*)NameSub,' starting iProc=',iProcWorld

    ! After everything is initialized exclude PEs which are not involved
    if(.not.UseMe) RETURN

    if(DoTest)write(*,*)NameSub,' starting, iProc=',iProcWorld

    if(DoTest)write(*,*)NameSub,' finished, iProc=',iProcWorld

  end subroutine couple_gm_pt

end module CON_couple_gm_pt
