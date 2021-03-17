!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE IE
!^CMP FILE PW

!
! Couple PW and GM components both ways.
! Note that PW -> GM has to be done before GM -> PW
! so the PW
!
module CON_couple_gm_pw

  use CON_coupler
  use CON_transfer_data, ONLY: transfer_real_array

  use GM_wrapper, ONLY: GM_get_for_pw, GM_put_from_pw
  use PW_wrapper, ONLY: PW_get_for_gm, PW_put_from_gm

  implicit none

  private ! except

  public :: couple_gm_pw_init ! initialize coupling
  public :: couple_gm_pw      ! couple GM to PW
  public :: couple_pw_gm      ! couple PW to GM

  ! revision history:
  ! 01/18/2007 A.Glocer and G.Toth - initial version

  ! Communicator and logicals to simplify message passing and execution
  logical :: IsInitialized = .false.

  ! PW grid
  integer, save :: nTotalLine

contains
  !============================================================================

  subroutine couple_gm_pw_init

    ! This subroutine should be called from all PE-s so that
    ! a union group can be formed. The PW grid size is also stored.

    !--------------------------------------------------------------------------
    if(IsInitialized) RETURN
    IsInitialized = .true.

    nTotalLine = Grid_C(PW_) % nCoord_D(2)

  end subroutine couple_gm_pw_init
  !============================================================================

  subroutine couple_pw_gm(tSimulation)

    real, intent(in) :: tSimulation     ! simulation time at coupling

    ! Couple between two components:\
    !    Polar Wind (PW)  source\
    !    Global Magnetosphere (GM) target
    !
    ! Send electrostatic potential and field aligned current from PW to IE.

    ! Variables to pass
    integer, parameter :: nVar = 8

    character (len=*), parameter :: NameVar_V(nVar) = [     &
         'CoLat    ', 'Longitude', 'Density1 ', 'Density2 ', &
         'Density3 ', 'Velocity1', 'Velocity2', 'Velocity3' ]

    ! Buffer for the potential on the 2D PW grid
    real, allocatable :: Buffer_VI(:,:)

    logical :: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'couple_pw_gm'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTest)write(*,*)NameSub,' starting, iProc=', i_proc()

    ! Send boundary density from PW to GM
    allocate(Buffer_VI(nVar,nTotalLine))
    if(is_proc(PW_)) call PW_get_for_gm(Buffer_VI, nVar, nTotalLine, &
         NameVar_V, tSimulation)

    ! PW data is only available on the root
    call transfer_real_array(PW_, GM_, nVar*nTotalLine, Buffer_VI)

    if(is_proc(GM_))call GM_put_from_pw(Buffer_VI, nVar, nTotalLine, NameVar_V)
    deallocate(Buffer_VI)

    if(DoTest)write(*,*)NameSub,' finished, iProc=', i_proc()

  end subroutine couple_pw_gm
  !============================================================================

  subroutine couple_gm_pw(tSimulation)

    real, intent(in) :: tSimulation     ! simulation time at coupling

    ! Couple between two components:\
    !    Global Magnetosphere (GM) source \
    !    Polar Wind (PW)  target
    !
    ! Send pressure from GM to PW

    ! Name of variables
    character (len=*), parameter:: NameVar_D(2) = ['CoLat    ','Longitude']

    ! Buffer for the coordinates of PW field lines
    real, allocatable:: Coord_DI(:,:)

    ! Buffer for the pressure on the PW grid
    real, allocatable :: Buffer_I(:)

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'couple_gm_pw'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    if(DoTest)write(*,*)NameSub,' starting, iProc=', i_proc()

    ! Send PW field line coordinates to GM
    allocate(Coord_DI(2,nTotalLine))
    if(is_proc(PW_))  &
         call PW_get_for_gm(Coord_DI, 2, nTotalLine, NameVar_D, tSimulation)
    call transfer_real_array(PW_, GM_, 2*nTotalLine, Coord_DI)
    if(is_proc(GM_))call GM_put_from_pw(Coord_DI, 2, nTotalLine, NameVar_D)
    deallocate(Coord_DI)

    ! Send pressure from GM to PW
    allocate(Buffer_I(nTotalLine))
    if(is_proc(GM_))call GM_get_for_pw(nTotalLine, Buffer_I)
    call transfer_real_array(GM_, PW_, nTotalLine, Buffer_I)
    if(is_proc(PW_)) call PW_put_from_gm(nTotalLine,Buffer_I)
    deallocate(Buffer_I)

    if(DoTest)write(*,*)NameSub,' finished, iProc=',i_proc()

  end subroutine couple_gm_pw
  !============================================================================

end module CON_couple_gm_pw
!==============================================================================

