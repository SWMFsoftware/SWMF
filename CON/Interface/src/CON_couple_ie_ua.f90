!^CMP COPYRIGHT UM
!^CMP FILE IE
!^CMP FILE UA

!BOP
!MODULE: CON_couple_ie_ua - couple IE and UA components
!
!DESCRIPTION:
! Couple IE and UA components both ways.
!
!INTERFACE:
module CON_couple_ie_ua

  !USES:
  use CON_coupler

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_ie_ua_init ! initialize both couplings
  public :: couple_ie_ua      ! couple IE to UA
  public :: couple_ua_ie      ! couple UA to IE

  !REVISION HISTORY:
  ! 08/25/2003 A.Ridley <ridley@umich.edu> - initial version as external
  !                                          subroutines
  ! 08/27/2003 G.Toth <gtoth@umich.edu>    - combined them into a module
  !EOP

  ! Communicator and logicals to simplify message passing and execution
  integer, save :: iCommIeUa, iProc0Ua
  logical :: UseMe=.true., IsInitialized = .false.

  ! Size of the 2D spherical structured IE grid
  integer, save :: iSize, jSize, nCells_D(2)

contains

  !BOP =======================================================================
  !IROUTINE: couple_ie_ua_init - initialize IE-UA couplings
  !INTERFACE:
  subroutine couple_ie_ua_init

    use ModNumConst, only:cPi

    ! General error code
    integer :: iError, i, j

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s so that
    ! a union group can be formed. The IE grid size is also stored.
    !EOP

    if(IsInitialized) RETURN
    IsInitialized = .true.

    iProc0Ua = i_proc0(UA_)
    call set_router_comm(IE_,UA_,iCommIeUa,UseMe,iProc0Ua)

    ! This works for a NODE BASED regular IE grid only
    nCells_D = ncells_decomposition_d(IE_) 
    iSize=nCells_D(1); jSize=nCells_D(2)

    write(*,*) "Initializing COUPLER!!", iSize, jSize

    if (is_proc(UA_)) then
       ! This won't work if we are using AMIE in RIM!
       call EIE_InitGrid(iSize, jSize+1, 1, iError)
       call EIE_FillLats(90.0-(Grid_C(IE_) % Coord1_I)*180.0/cPi,iError)
       call EIE_FillMltsOffset((Grid_C(IE_) % Coord2_I)*180.0/cPi,iError)
    endif

  end subroutine couple_ie_ua_init

  !BOP =======================================================================
  !IROUTINE: couple_ie_ua - couple IE component to UA component
  !INTERFACE:
  subroutine couple_ie_ua(tSimulation)

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Ionosphere Electrodynamics (IE)  source\\
    !    Upper Atmosphere (UA) target
    !
    ! Send electrostatic potential from IE to UA.
    !EOP

    !\
    ! General coupling variables
    !/

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_ie_ua'

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

      character (len=*), parameter, dimension(3) :: &
           NameVar_V=(/'Pot','Ave','Tot'/)

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
      integer :: iVar       ! 1 for Pot, 2 for AveE, 3 for EFlux
      integer :: iProcFrom  ! PE number sending the potential for current block
      !------------------------------------------------------------------------

      ! After everything is initialized exclude PEs which are not involved
      if(.not.UseMe) RETURN

      if(DoTest)write(*,*)NameSubSub,' starting, iProc=',iProcWorld
      if(DoTest)write(*,*)NameSubSub,', iProc, UAi_iProc0, IEi_iProc0=', &
           iProcWorld,i_proc0(UA_),i_proc0(IE_)

      !\
      ! Allocate buffers both in UA and IE
      !/
      allocate(Buffer_II(iSize,jSize), stat=iError)
      call check_allocate(iError,NameSubSub//": Buffer_II")

      if(DoTest)write(*,*)NameSubSub,', variables allocated',&
           ', iProc:',iProcWorld

      do iVar = 1, 3

         !\
         ! Get potential from IE
         !/

         if(is_proc(IE_))  &
              call IE_get_for_ua(Buffer_II, iSize, jSize, &
              NameVar_V(iVar),tSimulation)

         !\
         ! Transfer variables from IE to UA
         !/ 
         
         iProcFrom = i_proc0(IE_)

         nSize = iSize*jSize

         if(iProcFrom /= i_proc0(UA_))then
            if(i_proc() == iProcFrom) &
                 call MPI_send(Buffer_II,nSize,MPI_REAL,i_Proc0(UA_),&
                 1,i_comm(),iError)
            if(is_proc0(UA_)) &
                 call MPI_recv(Buffer_II,nSize,MPI_REAL,iProcFrom,&
                 1,i_comm(),iStatus_I,iError)
         end if

         ! Broadcast variables inside UA
         if(n_proc(UA_)>1 .and. is_proc(UA_)) &
              call MPI_bcast(Buffer_II,nSize,MPI_REAL,0,i_comm(UA_),iError)

         if(DoTest)write(*,*)NameSubSub,', variables transferred',&
              ', iProc:',iProcWorld

         !\
         ! Put variables into UA
         !/
         if(is_proc(UA_))then
            call SPS_put_into_ie(Buffer_II, iSize, jSize, &
                 NameVar_V(iVar))
            if(DoTest) &
                 write(*,*)NameSubSub//' iProc, Buffer(1,1)=',&
                 iProcWorld,Buffer_II(1,1)
         end if

      enddo

      !\
      ! Deallocate buffer to save memory
      !/
      deallocate(Buffer_II)

      if(DoTest)write(*,*)NameSubSub,', variables deallocated',&
           ', iProc:',iProcWorld

      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',iProcWorld

    end subroutine couple_mpi

  end subroutine couple_ie_ua

  !BOP =======================================================================
  !IROUTINE: couple_ua_ie - couple UA to IE component
  !INTERFACE:
  subroutine couple_ua_ie(tSimulation)

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Upper Atmosphere           (UA) source\\
    !    Ionosphere Electrodynamics (IE) target
    !
    !EOP

    !\
    ! General coupling variables
    !/

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_ua_ie'

    logical :: DoTest, DoTestMe
    integer :: iProcWorld, iError
    !-------------------------------------------------------------------------

    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    iProcWorld = i_proc()

    if(DoTest)write(*,*)NameSub,' starting iProc=',iProcWorld

    if(is_proc(UA_)) call initialize_ie_ua_buffers(iError)

    call couple_mpi

    if(DoTest)write(*,*)NameSub,': finished iProc=',iProcWorld

  contains

    !==========================================================================
    subroutine couple_mpi

      character (len=*), parameter :: NameSubSub=NameSub//'::couple_mpi'

      ! Number of coordinates and number of variables to pass
      integer, parameter :: nDimLoc=3, nVarIeUa=8

      ! Names of locations for both blocks
      character (len=*), parameter, dimension(2) :: &
           NameLoc_B=(/'North','South'/)

      ! Names of variables for both blocks
      character (len=*), parameter :: NameVar='Lat:MLT:FAC:SigmaP:SigmaH'

      ! Buffer for the variables on the 2D IE grid
      real, dimension(:,:,:), allocatable :: Buffer_IIV

      real, dimension(:,:), allocatable :: UAr2_Mlts, UAr2_Lats
      real, dimension(:,:), allocatable :: UAr2_Hal, UAr2_Ped, UAr2_Fac

      integer :: iStart, iEnd, UAi_nLats, UAi_nMlts, nVarsToPass

      ! MPI related variables

      ! MPI status variable
      integer :: iStatus_I(MPI_STATUS_SIZE)

      ! General error code
      integer :: iError

      ! Message size
      integer :: nSize

      integer :: iBlock, iProcTo

      logical :: IsFirstTime = .true.

      !------------------------------------------------------------------------

      ! After everything is initialized exclude PEs which are not involved
      if(.not.UseMe) RETURN

      if(DoTest)write(*,*)NameSubSub,' starting, iProc=',iProcWorld
      if(DoTest)write(*,*)NameSubSub,', iProc, UAi_iProc0, i_proc0(IE_)=', &
           iProcWorld,i_proc0(UA_),i_proc0(IE_)

!      write(*,*) "Goal #1 is to just get the codes running!"
      write(*,*) "CON_couple_ie_ua -> couple_ua_ie is disabled"
!      write(*,*) "UAi_nMlts, UAi_nLats", UAi_nMlts, UAi_nLats

!      stop

!!!      !\
!!!      ! Calculate field aligned currents in UA
!!!      !/
!!!      if(is_proc(UA_))then
!!!
!!!         if(DoTest)write(*,*)NameSubSub,': call UA_calc_fac'
!!!
!!!         call UA_calc_electrodynamics(UAi_nMlts, UAi_nLats)
!!!
!!!         allocate(UAr2_Fac(UAi_nMlts, UAi_nLats), &
!!!              UAr2_Ped(UAi_nMlts, UAi_nLats), &
!!!              UAr2_Hal(UAi_nMlts, UAi_nLats), &
!!!              UAr2_Lats(UAi_nMlts, UAi_nLats), &
!!!              UAr2_Mlts(UAi_nMlts, UAi_nLats), stat=iError)
!!!         call check_allocate(iError,NameSubSub//': '//NameLoc_B(1))
!!!
!!!         call UA_fill_electrodynamics(UAr2_Fac, UAr2_Ped, UAr2_Hal, &
!!!              UAr2_Lats, UAr2_Mlts)
!!!
!!!         if(DoTest)write(*,*)i_proc(),':',NameSubSub,&
!!!              ' nMlts,nLats=',UAi_nMlts, UAi_nLats
!!!
!!!      end if
!!!      if(DoTest)write(*,*)NameSubSub,': bcast'
!!!
!!!      !\
!!!      ! Broadcast information about mapping points
!!!      !/
!!!
!!!      call MPI_bcast(UAi_nMlts,1,MPI_INTEGER,iProc0Ua,iCommIeUa,iError)
!!!      call MPI_bcast(UAi_nLats,1,MPI_INTEGER,iProc0Ua,iCommIeUa,iError)
!!!
!!!      nVarsToPass = 3
!!!      if (IsFirstTime) nVarsToPass = 5 
!!!
!!!      allocate(Buffer_IIV(UAi_nMlts, UAi_nLats/2, nVarsToPass), stat=iError)
!!!      call check_allocate(iError,NameSubSub//': '//'Buffer_IIV')
!!!
!!!      ! Do Northern and then Southern hemispheres
!!!      do iBlock = 1, 2
!!!
!!!         if(DoTest)write(*,*)NameSubSub,': iBlock=',iBlock
!!!
!!!         ! The IE processor for this block
!!!         iProcTo = pe_decomposition(IE_,iBlock)
!!!
!!!         if(DoTest)write(*,*)NameSubSub,': iProcTo=',iProcTo
!!!
!!!         if (is_proc0(UA_)) then
!!!            !\
!!!            ! UA goes from the South pole to the north pole, while IE goes
!!!            ! from the north pole to the south pole, so the blocks have to
!!!            ! be reversed, basically.
!!!            !/
!!!            if (iBlock == 1) then
!!!               iStart = UAi_nLats/2 + 1
!!!               iEnd   = UAi_nLats
!!!            else
!!!               iStart = 1
!!!               iEnd   = UAi_nLats/2
!!!            endif
!!!            Buffer_IIV(:,1:UAi_nLats/2,1) = UAr2_Fac(:,iStart:iEnd)
!!!            Buffer_IIV(:,1:UAi_nLats/2,2) = UAr2_Ped(:,iStart:iEnd)
!!!            Buffer_IIV(:,1:UAi_nLats/2,3) = UAr2_Hal(:,iStart:iEnd)
!!!            if (IsFirstTime) then
!!!               Buffer_IIV(:,1:UAi_nLats/2,4) = UAr2_Lats(:,iStart:iEnd)
!!!               Buffer_IIV(:,1:UAi_nLats/2,5) = UAr2_Mlts(:,iStart:iEnd)
!!!            endif
!!!
!!!            if(DoTest)write(*,*)i_proc(),':',NameSubSub,': transfer ',NameVar,&
!!!                 ' ',NameLoc_B(iBlock),&
!!!                 ' maxval(Hall)=',maxval(Buffer_IIV(:,40:UAi_nLats/2,3)),&
!!!                 ' maxloc(Hall)=',maxloc(Buffer_IIV(:,40:UAi_nLats/2,3))
!!!            
!!!         endif
!!!
!!!         !\
!!!         ! Transfer values from UA to IE
!!!         !/
!!!
!!!         if(iProcTo /= i_proc0(UA_))then
!!!            nSize = UAi_nMlts * (UAi_nLats/2) * nVarsToPass
!!!            if(is_proc0(UA_)) &
!!!                 call MPI_send(Buffer_IIV,nSize,MPI_REAL,iProcTo,&
!!!                 1,i_comm(),iError)
!!!            if(i_proc() == iProcTo) &
!!!                 call MPI_recv(Buffer_IIV,nSize,MPI_REAL,i_proc0(UA_),&
!!!                 1,i_comm(),iStatus_I,iError)
!!!         end if
!!!
!!!         if(DoTest .and. i_proc() == iProcTo) &
!!!              write(*,*)i_proc(),':',NameSubSub,': put ',NameVar, &
!!!              ' ',NameLoc_B(iBlock), &
!!!              ' maxval(Hall)=',maxval(Buffer_IIV(:,40:UAi_nLats/2,3)),&
!!!              ' maxloc(Hall)=',maxloc(Buffer_IIV(:,40:UAi_nLats/2,3))
!!!
!!!         !\
!!!         ! Put values into IE
!!!         !/
!!!         if(i_proc() == iProcTo )then
!!!            call IE_put_from_ua(Buffer_IIV, iBlock, &
!!!                 UAi_nMlts, UAi_nLats/2, nVarsToPass)
!!!         end if
!!!
!!!      end do
!!!
!!!      !\
!!!      ! Deallocate buffer to save memory
!!!      !/
!!!      deallocate(Buffer_IIV)
!!!
!!!      if(is_proc(UA_)) &
!!!           deallocate(UAr2_Fac, UAr2_Ped, UAr2_Hal, UAr2_Lats, UAr2_Mlts)

      IsFirstTime = .false.

      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',iProcWorld

    end subroutine couple_mpi

  end subroutine couple_ua_ie

end module CON_couple_ie_ua

