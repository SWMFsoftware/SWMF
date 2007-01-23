!^CMP COPYRIGHT UM
!^CMP FILE IE
!^CMP FILE PW

!BOP
!MODULE: CON_couple_pw_gm - couple PW and GM components
!
!DESCRIPTION:
! Couple PW and GM components both ways.
!
!INTERFACE:
module CON_couple_pw_gm

  !USES:
  use CON_coupler
  use ModPWOM, only: nTotalLine,nSpecies

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_pw_gm_init ! initialize coupling
  public :: couple_pw_gm      ! couple PW to GM

  !REVISION HISTORY:
  ! 01/18/2007 A.Glocer and G.Toth - initial version 
  !EOP

  ! Communicator and logicals to simplify message passing and execution
  integer, save :: iCommPwGm, iProc0Gm
  logical :: UseMe=.true., IsInitialized = .false.

  ! PW grid
  !integer, save :: iSize, jSize, nCells_D(2)

contains

  !BOP =======================================================================
  !IROUTINE: couple_pw_gm_init - initialize PW-GM couplings
  !INTERFACE:
  subroutine couple_pw_gm_init

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s so that
    ! a union group can be formed. The PW grid size is also stored.
    !EOP

    if(IsInitialized) RETURN
    IsInitialized = .true.

    iProc0Gm = i_proc0(GM_)
    call set_router_comm(PW_,GM_,iCommPwGm,UseMe,iProc0Gm)

    ! This works for a NODE BASED regular PW grid only
    !nCells_D = ncells_decomposition_d(PW_) + 1
    !iSize=nCells_D(1); jSize=nCells_D(2)

  end subroutine couple_pw_gm_init

  !BOP =======================================================================
  !IROUTINE: couple_ie_pw - couple IE component to PW component
  !INTERFACE:
  subroutine couple_pw_gm(tSimulation)

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Polar Wind (PW)  source\\
    !    Global Magnetosphere (GM) target
    !
    ! Send electrostatic potential and field aligned current from PW to IE.
    !EOP

    !\
    ! General coupling variables
    !/

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_pw_gm'

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

      character (len=*), parameter, dimension(2) :: &
           NameHem_B=(/'North','South'/)

      integer, parameter :: nVar = 8
      integer, parameter :: South_ = 1, North_ = 2

      character (len=*), parameter, dimension(nVar) :: &
           NameVar_V=(/'CoLat    ','Longitude','Density1 ','Density2 ',&
           'Density3 ','Velocity1','Velocity2','Velocity3'/)

      ! Buffer for the potential on the 2D PW grid
      real, dimension(:,:), allocatable :: Buffer_IIV

      ! MPI related variables

      ! MPI status variable
      integer :: iStatus_I(MPI_STATUS_SIZE)

      ! General error code
      integer :: iError

      ! Message size
      integer :: nSize

      integer :: iBlock     ! 1 for northern and 2 for southern hemisphere
      integer :: iVar       ! 1) CoLat, 2)Longitude, 3)Density1, 4)Density2
                            ! 5) Density3, 6)Velocity1, 7)Velocity2, 
                            ! 8) Velocity3
      integer :: iProcFrom  ! PE number sending the potential for current block
      !------------------------------------------------------------------------

      ! After everything is initialized exclude PEs which are not involved
      if(.not.UseMe) RETURN

      if(DoTest)write(*,*)NameSubSub,' starting, iProc=',iProcWorld
      if(DoTest)write(*,*)NameSubSub,', iProc, GMi_iProc0, PWi_iProc0=', &
           iProcWorld,i_proc0(GM_),i_proc0(PW_)

      !\
      ! Allocate buffers both in GM and PW
      !/
      allocate(Buffer_IIV(nTotalLine,nVar), stat=iError)
      call check_allocate(iError,NameSubSub//": Buffer_IIV")

      if(DoTest)write(*,*)NameSubSub,', variables allocated',&
           ', iProc:',iProcWorld

      do iBlock = South_, North_

         !\
         ! boundary density from PW
         !/

         if(is_proc(PW_))  &
              call PW_get_for_gm(Buffer_IIV, nTotalLine, &
              nVar, NameVar_V, NameHem_B(iBlock), tSimulation)

         !\
         ! Transfer variables from PW to GM
         !/ 

         iProcFrom = pe_decomposition(PW_,iBlock)

         nSize = nTotalLine*nVar

         if(iProcFrom /= i_proc0(GM_))then
            if(i_proc() == iProcFrom) &
                 call MPI_send(Buffer_IIV,nSize,MPI_REAL,i_Proc0(GM_),&
                 1,i_comm(),iError)
            if(is_proc0(GM_)) &
                 call MPI_recv(Buffer_IIV,nSize,MPI_REAL,iProcFrom,&
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
            !call GM_put_from_pw(Buffer_IIV, nTotalLine, nVar, &
            !     NameVar_V, iBlock)
            !if(DoTest) &
            !     write(*,*)NameSubSub//' iProc, Buffer(1,1)=',&
            !     iProcWorld,Buffer_IIV(1,1,:)
         end if

      enddo
      !\
      ! Deallocate buffer to save memory
      !/
      deallocate(Buffer_IIV)

      if(DoTest)write(*,*)NameSubSub,', variables deallocated',&
           ', iProc:',iProcWorld

      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',iProcWorld

    end subroutine couple_mpi

  end subroutine couple_pw_gm

end module CON_couple_pw_gm

