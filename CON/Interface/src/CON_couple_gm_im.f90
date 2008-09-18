!^CMP COPYRIGHT UM
!^CMP FILE GM
!^CMP FILE IM

!BOP
!MODULE: CON_couple_gm_im - couple GM and IM components
!
!DESCRIPTION:
! Couple GM and IM components both ways. 
!
!INTERFACE:
module CON_couple_gm_im

  !USES:
  use CON_coupler

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_gm_im_init ! initialize both couplings
  public :: couple_gm_im      ! couple GM to IM
  public :: couple_im_gm      ! couple IM to GM

  !REVISION HISTORY:
  ! 07/25/2003 G.Toth <gtoth@umich.edu> - initial version
  !            O.Volberg and D.DeZeeuw
  !
  ! 08/27/2003 G.Toth - external subroutines combined into a module
  !EOP

  ! Communicator and logicals to simplify message passing and execution
  integer, save :: iCommGmIm
  logical :: IsInitialized = .false., UseMe=.true.

  ! Size of the 2D spherical structured (possibly non-uniform) IM grid
  integer, save :: iSize, jSize, nCells_D(2)

  ! Number of satellites in GM that will also be traced in IM
  integer, save :: nShareSats       !!!DTW 2007

contains

  !BOP =======================================================================
  !IROUTINE: couple_gm_im_init - initialize GM-IM coupling
  !INTERFACE:
  subroutine couple_gm_im_init
    !DESCRIPTION:
    ! Store IM grid size.
    !EOP
    !------------------------------------------------------------------------

    ! MPI status variable
    integer :: iStatus_I(MPI_STATUS_SIZE)

    ! General error code
    integer :: iError

    if(IsInitialized) RETURN
    IsInitialized = .true.

    UseMe = is_proc(IM_) .or. is_proc(GM_)
    if(.not.UseMe) RETURN

    ! This works for a regular IM grid only
    nCells_D=ncells_decomposition_d(IM_)
    iSize=nCells_D(1); jSize=nCells_D(2)

    ! Set number of satellites shared between GM and IM for tracing.
    call GM_satinit_for_im(nShareSats)       !!!!DTW 2007


    if(i_proc0(IM_) /= i_proc0(GM_))then
       if(is_proc0(GM_)) &
            call MPI_send(nShareSats,1,MPI_INTEGER,i_proc0(IM_),&
            1,i_comm(),iError)
       if(is_proc0(IM_)) &
            call MPI_recv(nShareSats,1,MPI_INTEGER,i_proc0(GM_),&
            1,i_comm(),iStatus_I,iError)
    end if

  end subroutine couple_gm_im_init

  !BOP =======================================================================
  !IROUTINE: couple_gm_im - couple GM to IM component
  !INTERFACE:
  subroutine couple_gm_im(tSimulation)
    
    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Global Magnetosphere (GM) source\\
    !    Inner Magnetosphere  (IM) target
    !
    ! Send field line volumes, average density and pressure and
    ! geometrical information.
    !EOP

    !\
    ! General coupling variables
    !/

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_gm_im'

    logical :: DoTest, DoTestMe
    integer :: iProcWorld
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    iProcWorld = i_proc()
    call couple_mpi

    if(DoTest)write(*,*)NameSub,': finished iProc=',iProcWorld

    if(DoTest.and.is_proc0(IM_)) call IM_print_variables('GM')

  contains

    !==========================================================================
    subroutine couple_mpi

      character (len=*), parameter :: NameSubSub=NameSub//'.couple_mpi'

      ! Number of variables to pass
      integer, parameter :: nVarGmIm=6

      character (len=*), parameter :: NameVar='vol:z0x:z0y:bmin:rho:p'

      ! Buffer for the variables on the 2D IM grid
      real, dimension(:,:,:), allocatable :: Buffer_IIV

      ! Buffer for satellite locations   !!! DTW 2007
      real, dimension(:,:,:), allocatable :: SatPos_DII

      ! Buffer for satellite names   !!! DTW 2007
      character (len=100), dimension(:), allocatable:: NameSat_I

      ! MPI related variables
      integer :: iProc0Im, iComm

      ! MPI status variable
      integer :: iStatus_I(MPI_STATUS_SIZE)

      ! General error code
      integer :: iError

      ! Message size
      integer :: nSize

      !------------------------------------------------------------------------

      ! After everything is initialized exclude PEs which are not involved
      if(.not.UseMe) RETURN

      if(DoTest)write(*,*)NameSubSub,' starting, iProc=',iProcWorld
      if(DoTest)write(*,*)NameSubSub,', iProc, GMi_iProc0, i_proc0(IM_)=', &
           iProcWorld,i_proc0(GM_),i_proc0(IM_)

      !\
      ! Allocate buffers both in GM and IM
      !/

      allocate(Buffer_IIV(iSize,jSize,nVarGmIm), stat=iError)
      call check_allocate(iError,NameSubSub//": Buffer_IIV")
      !!! DTW 2007
      if (nShareSats > 0) then
         allocate(SatPos_DII(3,2,nShareSats), stat=iError)
         call check_allocate(iError,NameSubSub//": SatPos_DII")
         allocate(NameSat_I(nShareSats),       stat=iError)
         call check_allocate(iError,NameSubSub//": NameSat_I")
      end if

      if(DoTest)write(*,*)NameSubSub,', variables allocated',&
           ', iProc:',iProcWorld


      iProc0Im = i_proc0(IM_)
      iComm = i_comm()

      !\
      ! Get field line integrals from GM
      !/
      if(is_proc(GM_)) &
           call GM_get_for_im(Buffer_IIV,iSize,jSize,nVarGmIm,NameVar)

      !\
      ! If IM sat tracing is enabled, get sat locations from GM
      !/
      if(is_proc(GM_).AND.(nShareSats > 0)) &
           call GM_get_sat_for_im(SatPos_DII, NameSat_I, nShareSats)
           

      !\
      ! Transfer physical variables from GM to IM
      !/
      if(i_proc0(IM_) /= i_proc0(GM_))then
         nSize = iSize*jSize*nVarGmIm
         if(is_proc0(GM_)) then 
            call MPI_send(Buffer_IIV,nSize,MPI_REAL,i_proc0(IM_),&
                 1,i_comm(),iError)
         endif
         if(is_proc0(IM_)) then
            call MPI_recv(Buffer_IIV,nSize,MPI_REAL,i_proc0(GM_),&
                 1,i_comm(),iStatus_I,iError)
         endif
      end if

      !\
      ! Transfer satellite names from GM to IM   !!!DTW 2007
      !/   

      if(nShareSats > 0 .and. i_proc0(IM_) /= i_proc0(GM_))then
         nSize = nShareSats*100
         if(is_proc0(GM_)) then
            call MPI_send(NameSat_I,nSize,MPI_BYTE,iProc0Im,&
                 1,iComm,iError)
         endif
         if(is_proc0(IM_)) then
            call MPI_recv(NameSat_I,nSize,MPI_BYTE,i_proc0(GM_),&
                 1,iComm,iStatus_I,iError)
         endif

         ! Transfer satellite locations from GM to IM   !!!DTW 2007

         nSize = 3*2*nShareSats
         if(is_proc0(GM_)) then
            call MPI_send(SatPos_DII,nSize,MPI_REAL,iProc0Im,&
                 1,iComm,iError)
         endif
         if(is_proc0(IM_)) then
            call MPI_recv(SatPos_DII,nSize,MPI_REAL,i_proc0(GM_),&
                 1,iComm,iStatus_I,iError)
         endif
      end if

      if(DoTest)write(*,*)NameSubSub,', variables transferred',&
           ', iProc:',iProcWorld

      !\
      ! Put variables into IM
      !/
      if(is_proc0(IM_))then
         call IM_put_from_gm(Buffer_IIV,iSize,jSize,nVarGmIm,NameVar)
         if(nShareSats > 0) &
              call IM_put_sat_from_gm(nShareSats, NameSat_I, SatPos_DII)
         if(DoTest) &
              write(*,*)'get_fieldline_volume: IM iProc, Buffer(1,1)=',&
              iProcWorld,Buffer_IIV(1,1,:)
      end if

      !\
      ! Deallocate buffer to save memory
      !/
      deallocate(Buffer_IIV)

      !!! DTW 2007
      if (nShareSats > 0) then
         deallocate(NameSat_I)
         deallocate(SatPos_DII)
      end if

      if(DoTest)write(*,*)NameSubSub,', variables deallocated',&
           ', iProc:',iProcWorld

      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',iProcWorld

    end subroutine couple_mpi

  end subroutine couple_gm_im

  !BOP =======================================================================
  !IROUTINE: couple_im_gm - couple IM to GM component
  !INTERFACE:
  subroutine couple_im_gm(tSimulation)

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Inner Magnetosphere  (IM) source\\
    !    Global Magnetosphere (GM) target
    !
    ! Send pressure from IM to GM.
    !EOP

    !\
    ! General coupling variables
    !/

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_im_gm'

    logical :: DoTest, DoTestMe
    integer :: iProcWorld

    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    iProcWorld = i_proc()
    call couple_mpi

    if(DoTest)write(*,*)NameSub,': finished iProc=',iProcWorld

    if(DoTest.and.is_proc0(GM_)) call GM_print_variables('IM')

  contains

    !=========================================================================
    subroutine couple_mpi

      character (len=*), parameter :: NameSubSub=NameSub//'.couple_mpi'

      ! Number of variables to pass
      integer, parameter :: nVarImGm=2

      ! Variable to pass is pressure
      character (len=*), parameter :: NameVar='p:rho'

      ! Buffer for the variables on the 2D IM grid
      real, dimension(:,:,:), allocatable :: Buffer_IIV

      ! MPI related variables

      ! MPI status variable
      integer :: iStatus_I(MPI_STATUS_SIZE)

      ! General error code
      integer :: iError

      ! Message size
      integer :: nSize

      ! Communicator and logicals to simplify message passing and execution
      integer, save :: iCommGmIm
      logical :: IsUninitialized = .true., UseMe=.true.
      !-----------------------------------------------------------------------

      ! After everything is initialized exclude PEs which are not involved
      if(.not.UseMe) RETURN

      if(DoTest)write(*,*)NameSubSub,' starting, iProc=',iProcWorld
      if(DoTest)write(*,*)NameSubSub,', iProc, GMi_iProc0, IMi_iProc0=', &
           iProcWorld,i_proc0(GM_),i_proc0(IM_)

      !\
      ! Allocate buffers both in GM and IM
      !/
      allocate(Buffer_IIV(iSize,jSize,nVarImGm), stat=iError)
      call check_allocate(iError,NameSubSub//": Buffer_IIV")

      if(DoTest)write(*,*)NameSubSub,', variables allocated',&
           ', iProc:',iProcWorld

      !\
      ! Get pressure from IM
      !/
      if(is_proc(IM_)) &
           call IM_get_for_gm(Buffer_IIV,iSize,jSize,nVarImGm,NameVar)

      !\
      ! Transfer variables from IM to GM
      !/ 
      nSize = iSize*jSize*nVarImGm
      if(i_proc0(IM_) /= i_proc0(GM_))then
         if(is_proc0(IM_)) &
              call MPI_send(Buffer_IIV,nSize,MPI_REAL,i_Proc0(GM_),&
              1,i_comm(),iError)
         if(is_proc0(GM_)) &
              call MPI_recv(Buffer_IIV,nSize,MPI_REAL,i_proc0(IM_),&
              1,i_comm(),iStatus_I,iError)
      end if

      ! Broadcast variables inside GM
      if(n_proc(GM_)>1 .and. is_proc(GM_)) &
           call MPI_bcast(Buffer_IIV,nSize,MPI_REAL,0,i_comm(GM_),iError)

      if(DoTest)write(*,*)NameSubSub,', variables transferred',&
           ', iProc:',iProcWorld

      !\
      ! Put variables into GM
      !/
      if(is_proc(GM_))then
         call GM_put_from_im(Buffer_IIV,iSize,jSize,nVarImGm,NameVar)
         if(DoTest) &
              write(*,*)NameSubSub//' iProc, Buffer(1,1,1)=',&
              iProcWorld,Buffer_IIV(1,1,1)
      end if

      !\
      ! Deallocate buffer to save memory
      !/
      deallocate(Buffer_IIV)

      if(DoTest)write(*,*)NameSubSub,', variables deallocated',&
           ', iProc:',iProcWorld

      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',iProcWorld

    end subroutine couple_mpi

  end subroutine couple_im_gm


end module CON_couple_gm_im
