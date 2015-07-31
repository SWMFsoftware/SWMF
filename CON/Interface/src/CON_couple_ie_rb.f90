!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
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
  use CON_transfer_data, ONLY: transfer_real_array

  use IE_wrapper, ONLY: IE_get_for_rb
  use RB_wrapper, ONLY: RB_put_from_ie

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_ie_rb_init ! initialize coupling
  public :: couple_ie_rb      ! couple IE to RB

  !REVISION HISTORY:
  ! 11/17/2006 A.Glocer and G.Toth - initial version 
  !EOP

  ! Communicator and logicals to simplify message passing and execution
  logical :: IsInitialized = .false.

  ! Size of the 2D spherical structured IE grid
  integer, save :: iSize, jSize

contains

  !BOP =======================================================================
  !IROUTINE: couple_ie_rb_init - initialize IE-RB couplings
  !INTERFACE:
  subroutine couple_ie_rb_init

    !DESCRIPTION:
    ! This subroutin is called from all PE-s. The IE grid size is stored.
    !EOP

    if(IsInitialized) RETURN
    IsInitialized = .true.

    ! Store the IE grid size. 
    ! Only North hemisphere is used. 
    ! The node based grid has 2*nTheta-1 nodes (there is 1 node at equator)
    iSize = (Grid_C(IE_) % nCoord_D(1) + 1)/2
    jSize =  Grid_C(IE_) % nCoord_D(2)

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
    ! Send electrostatic potential from IE to RB.
    !EOP

    ! "block" index for IE model (south = 2)
    integer, parameter :: North_ = 1

    ! Variable to pass is potential
    integer, parameter :: nVar = 1

    character (len=*), parameter, dimension(nVar) :: &
         NameVar_V=(/'Pot'/)

    ! Buffer for the potential on the 2D IE grid
    real, allocatable :: Buffer_IIV(:,:,:)

    logical :: DoTest, DoTestMe
    character (len=*), parameter :: NameSub='couple_ie_rb'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTest)write(*,*)NameSub,' starting, iProc=', i_proc()

    ! Get potential from IE to RB for the North hemisphere
    allocate(Buffer_IIV(iSize,jSize,nVar))
    if(is_proc0(IE_)) call IE_get_for_rb(Buffer_IIV, iSize, jSize, &
         nVar, NameVar_V, 'North', tSimulation)

    ! The potential is only reduced to the IE root
    call transfer_real_array(IE_, RB_, size(Buffer_IIV), Buffer_IIV)

    if(is_proc(RB_)) call RB_put_from_ie(Buffer_IIV, iSize, jSize, nVar, &
         NameVar_V, North_)
    deallocate(Buffer_IIV)

    if(DoTest)write(*,*)NameSub,': finished iProc=', i_proc()

  end subroutine couple_ie_rb

end module CON_couple_ie_rb
