!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE GM
!^CMP FILE IE

!BOP
!MODULE: CON_couple_gm_ie - couple GM and IE components
!
!DESCRIPTION:
! Couple GM and IE components both ways. 
!
!INTERFACE:
module CON_couple_gm_ie

  !USES:
  use CON_coupler

  use GM_wrapper, ONLY: GM_get_for_ie, GM_put_from_ie, GM_put_mag_from_ie, &
       GM_print_variables
  use IE_wrapper, ONLY: IE_get_for_gm, IE_groundmaginit_for_gm, &
       IE_get_mag_for_gm, IE_put_from_gm

  use CON_transfer_data, ONLY: transfer_integer, transfer_real_array

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_gm_ie_init ! initialize both couplings
  public :: couple_gm_ie      ! couple GM to IE
  public :: couple_ie_gm      ! couple IE to GM

  !REVISION HISTORY:
  ! 07/25/2003 G.Toth <gtoth@umich.edu> - initial version as external 
  !                                       subroutines
  ! 08/27/2003 G.Toth - combined into a module
  ! 12/01/2004 G.Toth - the GM->IE coupling is rewritten for Jr(iSize,jSize)
  !EOP

  ! Communicator and logicals to simplify message passing and execution
  logical       :: IsInitialized=.false.

  ! Size of the 2D spherical structured IE grid
  integer, save :: iSize, jSize, nCells_D(2)

  ! Number of shared ground-based magnetometers
  integer, save :: nShareGroundMag

contains

  !BOP =======================================================================
  !IROUTINE: couple_gm_ie_init - initialize GM-IE couplings
  !INTERFACE:
  subroutine couple_gm_ie_init

    integer :: iCommWorld
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='couple_gm_ie_init'

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s
    ! Store IE grid size and calculate union communicator.
    !EOP
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    if(IsInitialized) RETURN
    IsInitialized = .true.
    
    iCommWorld=i_comm()

    !!!   ! This works for a NODE BASED regular IE grid only
    !!!   nCells_D = ncells_decomposition_d(IE_) + 1
    !!!   iSize = nCells_D(1); jSize = nCells_D(2)

    iSize = Grid_C(IE_) % nCoord_D(1)
    jSize = Grid_C(IE_) % nCoord_D(2)

    ! Get number of shared magnetometers from IE to the GM processors
    if(is_proc(IE_)) call IE_groundmaginit_for_gm(nShareGroundMag)
    call transfer_integer(IE_, GM_, nShareGroundMag, UseSourceRootOnly=.false.)

    if(DoTest) &
         write(*,'(a,i3.3)') 'nShareGroundMag=', nShareGroundMag

  end subroutine couple_gm_ie_init

  !BOP =======================================================================
  !IROUTINE: couple_ie_gm - couple IE component to GM component
  !INTERFACE:
  subroutine couple_ie_gm(tSimulation)

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Inner Magnetosphere (IE)  source\\
    !    Global Magnetosphere (GM) target
    !
    ! The IE component sends the electrostatic potential to GM.
    ! GM can use that to calculate the velocities at the inner boundaries.
    !EOP

    !\
    ! General coupling variables
    !/

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_ie_gm'

    logical :: DoTest, DoTestMe
    integer :: iProcWorld

    ! Buffer for the potential on the 2D IE grid
    real, dimension(:,:,:), allocatable :: Buffer_IIV

    ! Buffer for magnetometer sharing.
    real, dimension(:,:), allocatable   :: Buffer_DI

    ! Number of variables to pass (potential,jouleheating)
    integer, parameter :: nVar = 2
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    iProcWorld = i_proc()

    if(DoTest)write(*,*)NameSub,' starting, iProc=',iProcWorld
    if(DoTest)write(*,*)NameSub,', iProc, GMi_iProc0, IEi_iProc0=', &
         iProcWorld,i_proc0(GM_),i_proc0(IE_)

    ! Transfer variables on IE grid
    allocate(Buffer_IIV(iSize,jSize,nVar))
    if(is_proc(IE_))call IE_get_for_gm(Buffer_IIV, iSize, jSize, tSimulation)

    ! The IE grid data is available from the ROOT IE processor ONLY!
    call transfer_real_array(IE_, GM_, iSize*jSize*nVar, Buffer_IIV)

    if(is_proc(GM_))call GM_put_from_ie(Buffer_IIV, iSize, jSize)
    deallocate(Buffer_IIV)

    if(nShareGroundMag>0)then
       allocate(Buffer_DI(3,nShareGroundMag))
       ! Collect magnetometer values if sharing ground mags.
       if(is_proc(IE_))call IE_get_mag_for_gm(Buffer_DI, nShareGroundMag)

       ! The magnetometer data is allreduced onto all the IE processors
       call transfer_real_array(IE_, GM_, 3*nShareGroundMag, Buffer_DI, &
            UseSourceRootOnly=.false.)
       if(is_proc(GM_))call GM_put_mag_from_ie(Buffer_DI, nShareGroundMag)
       deallocate(Buffer_DI)
    end if

    if(DoTest)write(*,*)NameSub,' finished, iProc=',iProcWorld
    if(DoTest.and.is_proc0(GM_)) call GM_print_variables('IE')

  end subroutine couple_ie_gm

  !BOP =======================================================================
  !IROUTINE: couple_gm_ie - couple GM to IE
  !INTERFACE:
  subroutine couple_gm_ie(tSimulation)

    !INPUT ARGUMENT:
    real, intent(in) :: tSimulation

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Global Magnetosphere       (GM) source\\
    !    Ionosphere Electrodynamics (IE) target
    !
    ! Send field aligned currents from GM to IE. 
    !EOP

    ! Number of variables to pass
    integer, parameter :: nVar=6

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_gm_ie'

    logical :: DoTest, DoTestMe
    integer :: iProcWorld

    ! Buffer for the field aligned current on the 2D IE grid
    real, dimension(:,:,:), allocatable :: Buffer_IIV
    !-------------------------------------------------------------------------

    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    iProcWorld = i_proc()

    if(DoTest)write(*,*)NameSub,' starting iProc=',iProcWorld

    if(DoTest)write(*,*)NameSub,' starting, iProc=',iProcWorld
    if(DoTest)write(*,*)NameSub,', iProc, GMi_iProc0, i_proc0(IE_)=', &
         iProcWorld,i_proc0(GM_),i_proc0(IE_)

    ! Allocate buffers for the variables both in GM and IE
    allocate(Buffer_IIV(iSize, jSize, nVar))

    ! Calculate field aligned currents on GM
    if(is_proc(GM_)) call GM_get_for_ie(Buffer_IIV, iSize, jSize, nVar)

    ! The result is on the root processor of GM only
    call transfer_real_array(GM_, IE_, iSize*jSize*nVar, Buffer_IIV)

    ! Put variables into IE
    if(is_proc(IE_)) call IE_put_from_gm(Buffer_IIV, iSize, jSize, nVar)

    ! Deallocate buffer to save memory
    deallocate(Buffer_IIV)

    if(DoTest)write(*,*)NameSub,' finished, iProc=',iProcWorld

  end subroutine couple_gm_ie

end module CON_couple_gm_ie
