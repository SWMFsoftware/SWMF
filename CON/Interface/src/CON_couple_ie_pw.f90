!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
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
  use CON_transfer_data, ONLY: transfer_real_array

  use IE_wrapper, ONLY: IE_get_for_pw
  use PW_wrapper, ONLY: PW_put_from_ie

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_ie_pw_init ! initialize coupling
  public :: couple_ie_pw      ! couple IE to PW

  !REVISION HISTORY:
  ! 11/17/2006 A.Glocer and G.Toth - initial version 
  !EOP

  ! Communicator and logicals to simplify message passing and execution
  logical :: IsInitialized = .false.

  ! Size of the 2D spherical structured IE grid
  integer, save :: iSize, jSize

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

    ! This works for a NODE BASED regular IE grid only
    iSize = Grid_C(IE_) % nCoord_D(1)
    jSize = Grid_C(IE_) % nCoord_D(2)

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
    ! Send electrostatic potential, field aligned current,
    ! average energy and electron flux from IE to PW.
    !EOP

    ! "block" index for IE model (south = 2)
    integer, parameter :: North_ = 1

    ! Number of variables to pass
    integer, parameter :: nVar = 4

    ! Names of variables to pass
    character (len=*), parameter, dimension(nVar) :: &
         NameVar_V = (/'Pot','Jr ','Ave','Tot'/)

    ! Buffer for the variables on the 2D IE grid
    real, allocatable :: Buffer_IIV(:,:,:)

    logical :: DoTest, DoTestMe
    character (len=*), parameter :: NameSub='couple_ie_pw'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    if(DoTest)write(*,*)NameSub,' starting, iProc=', i_proc()

    ! Get potential from IE to PW for the North hemisphere
    allocate(Buffer_IIV(iSize,jSize,nVar))
    if(is_proc(IE_)) call IE_get_for_pw(Buffer_IIV, iSize, jSize, &
         nVar, NameVar_V, 'North', tSimulation)

    ! The potential is only reduced to the IE root
    call transfer_real_array(IE_, PW_, iSize*jSize*nVar, Buffer_IIV)

    if(is_proc(PW_)) call PW_put_from_ie(Buffer_IIV, iSize, jSize, nVar, &
         NameVar_V, North_)
    deallocate(Buffer_IIV)

    if(DoTest)write(*,*)NameSub,': finished iProc=', i_proc()

  end subroutine couple_ie_pw

end module CON_couple_ie_pw

