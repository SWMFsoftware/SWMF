!^CMP COPYRIGHT UM
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
  logical       :: UseMe=.true., IsInitialized=.false.

  ! Size of the 2D spherical structured IE grid
  integer, save :: iSize, jSize, nCells_D(2)

  ! Number of shared ground-based magnetometers
  integer, save :: nShareGroundMag

contains

  !BOP =======================================================================
  !IROUTINE: couple_gm_ie_init - initialize GM-IE couplings
  !INTERFACE:
  subroutine couple_gm_ie_init

    integer :: iError, iCommWorld, iStatus_I(MPI_STATUS_SIZE)
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

    ! Set UseMe to .true. for the participating PE-s
    UseMe = is_proc(IE_) .or. is_proc(GM_)

    ! Get number of shared magnetometers from IE.
    if(is_proc(IE_)) &
         call IE_groundmaginit_for_gm(nShareGroundMag)
    ! Share with GM head node.
    if(i_proc0(IE_) /= i_proc0(GM_))then
       if(is_proc0(IE_)) &
            call MPI_send(nShareGroundMag, 1, MPI_INTEGER, i_proc0(GM_),&
            1, iCommWorld, iError)
       if(is_proc0(GM_))&
            call MPI_recv(nShareGroundMag, 1, MPI_INTEGER, i_proc0(IE_),&
            1, iCommWorld, iStatus_I, iError)
    end if

    ! Broadcast to all GM nodes.
    if(n_proc(GM_)>1 .and. is_proc(GM_)) &
         call MPI_bcast(nShareGroundMag, 1, MPI_INTEGER, 0, i_comm(GM_), iError)
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

    ! Variable to pass is potential
    character (len=*), parameter, dimension(2) :: &
         NameVar_B=(/'PotNorth','PotSouth'/)

    ! Buffer for the potential on the 2D IE grid
    real, dimension(:,:,:), allocatable :: Buffer_IIV
    ! Buffer for magnetometer sharing.
    real, dimension(:,:), allocatable   :: Buffer_DI
    ! Number of variables to pass (potential,jouleheating)
    integer, parameter :: nv=2

    ! MPI related variables

    ! MPI status variable
    integer :: iStatus_I(MPI_STATUS_SIZE)

    ! General error code
    integer :: iError

    ! Message size
    integer :: nSize

    integer :: iBlock     ! 1 for northern and 2 for southern hemisphere
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    iProcWorld = i_proc()

    ! After everything is initialized exclude PEs which are not involved
    if(.not.UseMe) RETURN

    if(DoTest)write(*,*)NameSub,' starting, iProc=',iProcWorld
    if(DoTest)write(*,*)NameSub,', iProc, GMi_iProc0, IEi_iProc0=', &
         iProcWorld,i_proc0(GM_),i_proc0(IE_)

    !\
    ! Allocate buffers both in GM and IE
    !/
    allocate(Buffer_IIV(iSize,jSize,nv), stat=iError)
    call check_allocate(iError,NameSub//": Buffer_IIV")
    if(nShareGroundMag>0) then
       allocate(Buffer_DI(3,nShareGroundMag), stat=iError)
       call check_allocate(iError, NameSub//": Buffer_DI")
    end if

    if(DoTest)write(*,*)NameSub,', variables allocated',&
         ', iProc:',iProcWorld
   
    ! Collect standard coupling variables.
    if(is_proc(IE_)) &
         call IE_get_for_gm(Buffer_IIV,iSize,jSize,tSimulation)
    ! Collect magnetometer values if sharing ground mags.
    if(is_proc(IE_) .AND. (nShareGroundMag>0)) &
         call IE_get_mag_for_gm(Buffer_DI, nShareGroundMag)


    nSize = iSize*jSize*nv
    if(i_proc0(IE_) /= i_proc0(GM_))then
       ! Pass standard variables.
       if(is_proc0(IE_)) &
            call MPI_send(Buffer_IIV,nSize,MPI_REAL,i_Proc0(GM_),&
            1,i_comm(),iError)
       if(is_proc0(GM_)) &
            call MPI_recv(Buffer_IIV,nSize,MPI_REAL,i_Proc0(IE_),&
            1,i_comm(),iStatus_I,iError)
       ! Pass shared magnetometers.
       if(is_proc0(IE_) .AND. (nShareGroundMag>0)) &
            call MPI_send(Buffer_DI,3*nShareGroundMag,MPI_REAL,i_Proc0(GM_),&
            1,i_comm(),iError)
       if(is_proc0(GM_) .AND. (nShareGroundMag>0)) &
            call MPI_recv(Buffer_DI,3*nShareGroundMag,MPI_REAL,i_Proc0(IE_),&
            1,i_comm(),iStatus_I,iError)
    end if

    ! Broadcast variables inside GM
    if(n_proc(GM_)>1 .and. is_proc(GM_)) &
         call MPI_bcast(Buffer_IIV,nSize,MPI_REAL,0,i_comm(GM_),iError)
    if(n_proc(GM_)>1 .and. is_proc(GM_) .and. (nShareGroundMag>0)) &
         call MPI_bcast(Buffer_DI,3*nShareGroundMag,MPI_REAL, &
         0,i_comm(GM_),iError)

    if(DoTest)write(*,*)NameSub,', variables transferred',&
         ', iProc:',iProcWorld

    !\
    ! Put variables into GM
    !/
    if(is_proc(GM_))then
       call GM_put_from_ie(Buffer_IIV,iSize,jSize)
       if(DoTest) &
            write(*,*)NameSub//' iProc, Buffer(1,1)=',&
            iProcWorld,Buffer_IIV(1,1,1:2)
    end if
    if(is_proc(GM_) .and. (nShareGroundMag>0))then
       call GM_put_mag_from_ie(Buffer_DI, nShareGroundMag)
       if(DoTest) &
            write(*,*)NameSub//' iProc, Buffer(1,1)=',&
            iProcWorld,Buffer_DI(1,1)
    end if

    !\
    ! Deallocate buffer to save memory
    !/
    deallocate(Buffer_IIV)
    if(nShareGroundMag>0) deallocate(Buffer_DI)

    if(DoTest)write(*,*)NameSub,', variables deallocated',&
         ', iProc:',iProcWorld

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

    !\
    ! General coupling variables
    !/

    ! Number of variables to pass
    integer, parameter :: nVar=6

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_gm_ie'

    logical :: DoTest, DoTestMe
    integer :: iProcWorld

    ! Names of variables for both blocks
    character (len=*), parameter, dimension(2) :: &
         NameVar_B = (/ 'JrNorth', 'JrSouth' /)


    ! Buffer for the field aligned current on the 2D IE grid
    real, dimension(:,:,:), allocatable :: Buffer_IIV

    ! MPI related variables

    ! MPI status variable
    integer :: iStatus_I(MPI_STATUS_SIZE)

    ! General error code
    integer :: iError

    ! Message size
    integer :: nSize

    integer :: iBlock, iProcTo
    !-------------------------------------------------------------------------

    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    iProcWorld = i_proc()

    if(DoTest)write(*,*)NameSub,' starting iProc=',iProcWorld


    ! After everything is initialized exclude PEs which are not involved
    if(.not.UseMe) RETURN

    if(DoTest)write(*,*)NameSub,' starting, iProc=',iProcWorld
    if(DoTest)write(*,*)NameSub,', iProc, GMi_iProc0, i_proc0(IE_)=', &
         iProcWorld,i_proc0(GM_),i_proc0(IE_)

    !\
    ! Allocate buffers for the variables both in GM and IE
    !/
    allocate(Buffer_IIV(iSize, jSize, nVar), stat=iError)
    call check_allocate(iError,NameSub)

    !\
    ! Calculate field aligned currents on GM
    ! The result will be on the root processor of GM
    !/
    if(is_proc(GM_)) &
         call GM_get_for_ie(Buffer_IIV, iSize, jSize, nVar)
    !\
    ! Transfer variables from GM to IE
    !/
    if(i_proc0(IE_) /= i_proc0(GM_))then
       nSize = iSize*jSize*nVar
       if(is_proc0(GM_)) &
            call MPI_send(Buffer_IIV, nSize, MPI_REAL, i_proc0(IE_),&
            1, i_comm(), iError)
       if(is_proc0(IE_)) &
            call MPI_recv(Buffer_IIV, nSize, MPI_REAL, i_proc0(GM_),&
            1,i_comm(), iStatus_I, iError)
    end if
    !\
    ! Put variables into IE
    !/
    if(is_proc(IE_)) &
         call IE_put_from_gm(Buffer_IIV, iSize, jSize, nVar)

    !\
    ! Deallocate buffer to save memory
    !/
    deallocate(Buffer_IIV)

    if(DoTest)write(*,*)NameSub,' finished, iProc=',iProcWorld

  end subroutine couple_gm_ie

end module CON_couple_gm_ie
