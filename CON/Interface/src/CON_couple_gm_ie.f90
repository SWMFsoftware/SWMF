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

contains

  !BOP =======================================================================
  !IROUTINE: couple_gm_ie_init - initialize GM-IE couplings
  !INTERFACE:
  subroutine couple_gm_ie_init

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s
    ! Store IE grid size and calculate union communicator.
    !EOP
    !------------------------------------------------------------------------
    if(IsInitialized) RETURN
    IsInitialized = .true.

    ! This works for a NODE BASED regular IE grid only
    nCells_D = ncells_decomposition_d(IE_) + 1
    iSize = nCells_D(1); jSize = nCells_D(2)

    ! Set UseMe to .true. for the participating PE-s
    UseMe = is_proc(IE_) .or. is_proc(GM_)

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

    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    iProcWorld = i_proc()
    call couple_mpi

    if(DoTest)write(*,*)NameSub,': finished iProc=',iProcWorld

    if(DoTest.and.is_proc0(GM_)) call GM_print_variables('IE')

  contains

    !=========================================================================
    subroutine couple_mpi

      character (len=*), parameter :: NameSubSub=NameSub//'.couple_mpi'

      ! Variable to pass is potential
      character (len=*), parameter, dimension(2) :: &
           NameVar_B=(/'PotNorth','PotSouth'/)

      ! Buffer for the potential on the 2D IE grid
      real, dimension(:,:), allocatable :: Buffer_II

      ! MPI related variables

      ! MPI status variable
      integer :: iStatus_I(MPI_STATUS_SIZE)

      ! General error code
      integer :: iError

      ! Message size
      integer :: nSize

      integer :: iBlock     ! 1 for northern and 2 for southern hemisphere
      integer :: iProcFrom  ! PE number sending the potential for current block

      !------------------------------------------------------------------------

      ! After everything is initialized exclude PEs which are not involved
      if(.not.UseMe) RETURN

      if(DoTest)write(*,*)NameSubSub,' starting, iProc=',iProcWorld
      if(DoTest)write(*,*)NameSubSub,', iProc, GMi_iProc0, IEi_iProc0=', &
           iProcWorld,i_proc0(GM_),i_proc0(IE_)

      !\
      ! Allocate buffers both in GM and IE
      !/
      allocate(Buffer_II(iSize,jSize), stat=iError)
      call check_allocate(iError,NameSubSub//": Buffer_II")

      if(DoTest)write(*,*)NameSubSub,', variables allocated',&
           ', iProc:',iProcWorld

      do iBlock = 1,2

         !\
         ! Get potential from IE
         !/
         if(is_proc(IE_))  &
              call IE_get_for_gm(Buffer_II,iSize,jSize,NameVar_B(iBlock),&
              tSimulation)
         !\
         ! Transfer variables from IE to GM
         !/ 
         iProcFrom = pe_decomposition(IE_,iBlock)

         nSize = iSize*jSize
         if(iProcFrom /= i_proc0(GM_))then
            if(i_proc() == iProcFrom) &
                 call MPI_send(Buffer_II,nSize,MPI_REAL,i_Proc0(GM_),&
                 1,i_comm(),iError)
            if(is_proc0(GM_)) &
                 call MPI_recv(Buffer_II,nSize,MPI_REAL,iProcFrom,&
                 1,i_comm(),iStatus_I,iError)
         end if

         ! Broadcast variables inside GM
         if(n_proc(GM_)>1 .and. is_proc(GM_)) &
              call MPI_bcast(Buffer_II,nSize,MPI_REAL,0,i_comm(GM_),iError)

         if(DoTest)write(*,*)NameSubSub,', variables transferred',&
              ', iProc:',iProcWorld

         !\
         ! Put variables into GM
         !/
         if(is_proc(GM_))then
            call GM_put_from_ie(Buffer_II,iSize,jSize,NameVar_B(iBlock))
            if(DoTest) &
                 write(*,*)NameSubSub//' iProc, Buffer(1,1)=',&
                 iProcWorld,Buffer_II(1,1)
         end if
      end do
      !\
      ! Deallocate buffer to save memory
      !/
      deallocate(Buffer_II)

      if(DoTest)write(*,*)NameSubSub,', variables deallocated',&
           ', iProc:',iProcWorld

      !\
      ! Map potential to inner BC of GM and calculate velocities
      !/
      if(is_proc(GM_))call GM_calc_iono_bcs

      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',iProcWorld

    end subroutine couple_mpi

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

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_gm_ie'

    logical :: DoTest, DoTestMe
    integer :: iProcWorld
    !-------------------------------------------------------------------------

    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    iProcWorld = i_proc()

    if(DoTest)write(*,*)NameSub,' starting iProc=',iProcWorld

    call couple_mpi

    if(DoTest)write(*,*)NameSub,': finished iProc=',iProcWorld

  contains

    !=========================================================================
    subroutine couple_mpi

      character (len=*), parameter :: NameSubSub = NameSub//'::couple_mpi'

      ! Names of variables for both blocks
      character (len=*), parameter, dimension(2) :: &
           NameVar_B = (/ 'JrNorth', 'JrSouth' /)


      ! Buffer for the field aligned current on the 2D IE grid
      real, dimension(:,:), allocatable :: Buffer_II

      ! MPI related variables

      ! MPI status variable
      integer :: iStatus_I(MPI_STATUS_SIZE)

      ! General error code
      integer :: iError

      ! Message size
      integer :: nSize

      integer :: iBlock, iProcTo
      !-----------------------------------------------------------------------

      ! After everything is initialized exclude PEs which are not involved
      if(.not.UseMe) RETURN

      if(DoTest)write(*,*)NameSubSub,' starting, iProc=',iProcWorld
      if(DoTest)write(*,*)NameSubSub,', iProc, GMi_iProc0, i_proc0(IE_)=', &
           iProcWorld,i_proc0(GM_),i_proc0(IE_)

      ! Do Northern and then Southern hemispheres in this order!
      do iBlock = 1, 2

         if(DoTest)write(*,*)NameSubSub,': iBlock=',iBlock

         ! The IE processor for this block
         iProcTo = pe_decomposition(IE_, iBlock)

         if(DoTest)write(*,*)NameSubSub,': iProcTo=',iProcTo

         !\
         ! Allocate buffers for the variables both in GM and IE
         !/
         allocate(Buffer_II(iSize, jSize), stat=iError)
         call check_allocate(iError,NameSubSub//": "//NameVar_B(iBlock))

         !\
         ! Calculate field aligned currents on GM
         ! The result will be on the root processor of GM
         !/
         if(is_proc(GM_)) &
              call GM_get_for_ie(Buffer_II, iSize, jSize, NameVar_B(iBlock))
         !\
         ! Transfer variables from GM to IE
         !/
         if(iProcTo /= i_proc0(GM_))then
            nSize = iSize*jSize
            if(is_proc0(GM_)) &
                 call MPI_send(Buffer_II, nSize, MPI_REAL, iProcTo,&
                 1, i_comm(), iError)
            if(i_proc() == iProcTo) &
                 call MPI_recv(Buffer_II, nSize, MPI_REAL, i_proc0(GM_),&
                 1,i_comm(), iStatus_I, iError)
         end if
         !\
         ! Put variables into IE
         !/
         if(i_proc() == iProcTo )&
              call IE_put_from_gm(Buffer_II, iSize, jSize, NameVar_B(iBlock))
         !\
         ! Deallocate buffer to save memory
         !/
         deallocate(Buffer_II)
      end do

      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',iProcWorld

    end subroutine couple_mpi

  end subroutine couple_gm_ie

end module CON_couple_gm_ie
