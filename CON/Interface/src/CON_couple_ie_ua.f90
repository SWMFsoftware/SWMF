!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE IE
!^CMP FILE UA

!
! Couple IE and UA components both ways.
!
module CON_couple_ie_ua

  use CON_coupler

  use IE_wrapper, ONLY: IE_get_for_ua, IE_put_from_ua

  implicit none

  private ! except

  public :: couple_ie_ua_init ! initialize both couplings
  public :: couple_ie_ua      ! couple IE to UA
  public :: couple_ua_ie      ! couple UA to IE

  ! revision history:
  ! 08/25/2003 A.Ridley <ridley@umich.edu> - initial version as external
  !                                          subroutines
  ! 08/27/2003 G.Toth <gtoth@umich.edu>    - combined them into a module
  ! 11/10/2021 Burleigh/Welling <dwelling@uta.edu>
  !                                        - Re-implemented & updated coupling.

  ! Communicator and logicals to simplify message passing and execution
  ! Initialization status:
  ! integer, save :: iCommIeUa, iProc0Ua
  logical :: UseMe=.true., IsInitialized = .false.

  ! Information about number and names of variables to share:
  integer, save :: nVarIeUa, nVarUaIe, nUaMagLon, nUaMagLat
  character(len=3), allocatable :: NameVarIeUa_V(:), NameVarUaIe_V(:)

  ! Size of the 2D spherical structured IE grid
  integer, save :: iSize, jSize, nCells_D(2), nRootBlock_d(3)

contains
  !============================================================================
  subroutine couple_ie_ua_init

    ! This subroutine should be called from all PE-s so that
    ! a union group can be formed. The IE grid size is also stored.
    ! Performs handshaking between UA and IE concerning names and number
    ! of variables to transfer.

    use CON_transfer_data, ONLY: transfer_integer, transfer_string_array
    use UA_wrapper, ONLY: UA_get_info_for_ie
    use IE_wrapper, ONLY: IE_get_info_for_ua

    ! General error code
    integer :: iError, i, j

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'couple_ie_ua_init'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTestMe) write(*,*) NameSub//' called;, IsInitialized=', IsInitialized
    if(IsInitialized) RETURN
    IsInitialized = .true.

    ! IE to UA coupling: set names and number of variables
    ! Get number of variables to be passed from UA to IE, pass to IE
    if(is_proc(UA_)) call UA_get_info_for_ie(nVarIeUa)
    call transfer_integer(UA_, IE_, nVarIeUa,  UseSourceRootOnly=.false.)

    ! Allocate the array holding the variable names.
    if(allocated(NameVarIeUa_V)) deallocate(NameVarIeUa_V)
    allocate(NameVarIeUa_V(nVarIeUa))

    ! Get variables names to be passed; transfer to IE
    ! Obtain number of magnetic (not total grid) lats/lons used by UA.
    if(is_proc(UA_)) call UA_get_info_for_ie(nVarIeUa, &
         NameVarIeUa_V, nUaMagLat, nUaMagLon)
    call transfer_integer(UA_, IE_, nUaMagLat, UseSourceRootOnly=.false.)
    call transfer_integer(UA_, IE_, nUaMagLon, UseSourceRootOnly=.false.)
    call transfer_string_array(UA_, IE_, nVarIeUa, NameVarIeUa_V, &
         UseSourceRootOnly=.false.)

    ! UA to IE coupling: set names and number of variables
    ! Get number of variables to be passed from IE to UA, pass to UA
    if(is_proc(IE_)) call IE_get_info_for_ua(nVarUaIe)
    call transfer_integer(IE_, UA_, nVarUaIe, UseSourceRootOnly=.false.)

    ! Allocate the array holding the variable names.
    if(allocated(NameVarUaIe_V)) deallocate(NameVarUaIe_V)
    allocate(NameVarUaIe_V(nVarUaIe))

    ! Get variables names to be passed; transfer to IE
    if(is_proc(IE_)) call IE_get_info_for_ua(nVarUaIe, NameVarUaIe_V)
    call transfer_string_array(IE_, UA_, nVarUaIe, NameVarUaIe_V, &
         UseSourceRootOnly=.false.)

    if(DoTestMe)then
       write(*,*) NameSub//' UA-IE coupling configuration:'
       write(*,*) '   IE requests ', nVarUaIe, ' variables from UA:', &
            NameVarUaIe_V
       write(*,*) '   UA requests ', nVarIeUa, ' variables from IE:', &
            NameVarIeUa_V
    end if

    ! NOT SURE IF NEEDED.
    ! This works for a NODE BASED regular IE grid only
    ! <so then why does IE initialize with SPS, a grid based thing??? MB>
    nCells_D = ncell_id(IE_)
    ! iSize=nCells_D(1); jSize=nCells_D(2)   ! orig. MB
    iSize = nCells_D(1) + 1; jSize = nCells_D(2) + 1 ! Grid size for 1 hemi.

    ! IE should share nVar and varNames to pass.

  end subroutine couple_ie_ua_init
  !============================================================================
  subroutine couple_ie_ua(tSimulation)

    ! Couple between two components:
    !    Ionosphere Electrodynamics (IE)  source
    !    Upper Atmosphere (UA) target
    !
    ! Send electrostatic potential from IE to UA.

    use CON_transfer_data, ONLY: transfer_real_array
    use IE_wrapper, ONLY: IE_get_for_ua
    use UA_wrapper, ONLY: UA_put_from_ie

    real, intent(in) :: tSimulation     ! simulation time at coupling

    ! Buffer for all shared variables on the 2D IE grid
    real, allocatable :: Buffer_IIV(:,:,:)

    ! Variables to assist with coupling
    integer :: nSize, iBlock, iError

    ! Debug variables:
    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'couple_ie_ua'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    ! Allocate buffers both in IE (source) and UA (target):
    allocate(Buffer_IIV(iSize,jSize,nVarIeUa))

    ! Transfer northern then southern hemisphere:
    do iBlock = 1, 2
       ! Get all variables from IE:
       if(is_proc(IE_)) call IE_get_for_ua(Buffer_IIV, iSize, jSize, &
            nVarIeUa, NameVarIeUa_V, iBlock, tSimulation)

       ! Transfer data:
       call transfer_real_array(IE_, UA_, iSize*jSize*nVarIeUa, Buffer_IIV)

       ! UA receives & handles data:
       if(is_proc(UA_)) call UA_put_from_ie(Buffer_IIV, iSize, jSize, &
            nVarIeUa, NameVarIeUa_V, iBlock)
    end do

    ! Deallocate buffer to save memory
    deallocate(Buffer_IIV)

  end subroutine couple_ie_ua
  !============================================================================
  subroutine couple_ua_ie(tSimulation)

    use CON_transfer_data, ONLY: transfer_real_array
    use UA_wrapper, ONLY: UA_get_for_ie
    use IE_wrapper, ONLY: IE_put_from_ua

    real, intent(in) :: tSimulation     ! simulation time at coupling

    ! Couple between two components:
    !    Upper Atmosphere           (UA) source
    !    Ionosphere Electrodynamics (IE) target
    !

    ! General coupling variables

    ! Buffer for the variables on the 2D IE grid: lon, lat, block, vars
    ! Always two blocks, one per hemisphere.
    real, allocatable :: Buffer_IIBV(:,:,:)

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'couple_ua_ie'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Allocate our transfer array:
    allocate(Buffer_IIBV(nUaMagLon, nUaMagLat, 2, nVarUaIe))

    ! Gather values from UA:
    if(is_proc(UA_)) call UA_get_for_ie(Buffer_IIBV, nUaMagLon, &
         nUaMagLat, nVarUaIe, NameVarUaIe_V)

    ! Transfer data:
    call transfer_real_array(UA_, IE_, nUaMagLon*nUaMagLat*2*nVarUaIe, &
         Buffer_IIBV)

    ! Distribute through IE:
    if(is_proc(IE_)) call IE_put_from_ua(Buffer_IIBV,  &
         nUaMagLon, nUaMagLat, nVarUaIe, NameVarUaIe_V)

    deallocate(Buffer_IIbV)

  end subroutine couple_ua_ie
  !============================================================================
end module CON_couple_ie_ua
!==============================================================================

