!^CMP COPYRIGHT UM
!^CMP FILE GM
!^CMP FILE RB

!BOP
!MODULE: CON_couple_gm_rb - couple GM and RB components
!
!DESCRIPTION:
! Couple GM and RB components one way for now. 
!
!INTERFACE:
module CON_couple_gm_rb

  !USES:
  use CON_coupler

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_gm_rb_init ! initialize both couplings
  public :: couple_gm_rb      ! couple GM to RB

  !REVISION HISTORY:
  ! 05/21/2004 O.Volberg - initial version
  !
  !EOP

  ! Communicator and logicals to simplify message passing and execution
  integer, save :: iCommGMRB
  logical :: IsInitialized = .false., UseMe=.true.

  ! Size of the 2D spherical structured (possibly non-uniform) RB grid
  integer, save :: iSize, jSize, nCells_D(2)

  ! Number of satellites in GM that will also be traced in RB
  integer, save :: nShareSats
contains

  !BOP =======================================================================
  !IROUTINE: couple_gm_rb_init - initialize GM-RB coupling
  !INTERFACE:
  subroutine couple_gm_rb_init
    integer :: iError
    ! MPI status variable
    integer :: iStatus_I(MPI_STATUS_SIZE)
    !DESCRIPTION:
    ! Store RB grid size.
    !EOP
    !------------------------------------------------------------------------
    if(IsInitialized) return
    IsInitialized = .true.

    UseMe = is_proc(RB_) .or. is_proc(GM_)
    if(.not.UseMe) return

    nCells_D=ncells_decomposition_d(RB_)
    iSize=nCells_D(1); jSize=nCells_D(2)

    ! Set number of satellites shared between GM and RB for tracing.
    call GM_satinit_for_rb(nShareSats)

    if(i_proc0(RB_) /= i_proc0(GM_))then
       if(is_proc0(GM_)) &
            call MPI_send(nShareSats,1,MPI_INTEGER,i_proc0(RB_),&
            1,i_comm(),iError)
       if(is_proc0(RB_)) &
            call MPI_recv(nShareSats,1,MPI_INTEGER,i_proc0(GM_),&
            1,i_comm(),iStatus_I,iError)
    end if

  end subroutine couple_gm_rb_init

  !BOP =======================================================================
  !IROUTINE: couple_gm_rb - couple GM to RB component
  !INTERFACE: couple_gm_rb(tSimulation)
  subroutine couple_gm_rb(tSimulation)

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Global Magnetosphere (GM) source\\
    !    Radiation Belt       (RB) target
    !
    ! Send field line volumes, average density and pressure and
    ! geometrical information.
    !EOP

    !\
    ! Coupling variables
    !/

    ! Number of integrals to pass
    integer, parameter :: nIntegral=4

    ! Names of variables to pass
    character (len=*), parameter :: NameVar='Z0x:Z0y:Z0b:I_I:S_I:R_I:B_I:IMF'

    ! Number of variables and points saved into the line data
    integer :: nVarLine, nPointLine

    ! Buffer for the variables on the 2D RB grid and line data
    real, allocatable :: Integral_IIV(:,:,:), BufferLine_VI(:,:)

    ! Buffer for satellite locations   
    real, dimension(:,:,:), allocatable :: SatPos_DII
    
    ! Buffer for satellite names   
    character (len=100), dimension(:), allocatable:: NameSat_I
    
    ! MPI related variables

    ! MPI status variable
    integer :: iStatus_I(MPI_STATUS_SIZE)

    ! General error code
    integer :: iError

    ! Message size
    integer :: nSize

    ! Name of this method
    character (len=*), parameter :: NameSub='couple_gm_rb'

    logical :: DoTest, DoTestMe
    integer :: iProcWorld

    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    iProcWorld = i_proc()

    ! After everything is initialized exclude PEs which are not involved
    if(.not.UseMe) RETURN

    if(DoTest)write(*,*)NameSub,' starting, iProc=',iProcWorld
    if(DoTest)write(*,*)NameSub,', iProc, GMi_iProc0, i_proc0(RB_)=', &
         iProcWorld,i_proc0(GM_),i_proc0(RB_)

    !\
    ! Allocate buffers both in GM and RB
    !/
    allocate(Integral_IIV(iSize,jSize,nIntegral), stat=iError)
    call check_allocate(iError,NameSub//": Integral_IIV")

    if (nShareSats > 0) then
       allocate(SatPos_DII(4,2,nShareSats), stat=iError)
       call check_allocate(iError,NameSub//": SatPos_DII")
       allocate(NameSat_I(nShareSats),       stat=iError)
       call check_allocate(iError,NameSub//": NameSat_I")
    end if

    if(DoTest)write(*,*)NameSub,', variables allocated',&
         ', iProc:',iProcWorld

    !\
    ! Get field line integrals from GM
    !/
    if(is_proc(GM_)) then
       call GM_get_for_rb_trace(iSize, jSize, NameVar, nVarLine, nPointLine)
       allocate(BufferLine_VI(nVarLine, nPointLine))
       call GM_get_for_rb(Integral_IIV, iSize, jSize, nIntegral, &
            BufferLine_VI, nVarLine, nPointLine, NameVar)
    end if
    !\
    ! If RB sat tracing is enabled, get sat locations from GM
    !/
    if(is_proc(GM_).AND.(nShareSats > 0)) &
         call GM_get_sat_for_rb(SatPos_DII, NameSat_I, nShareSats)

    !\
    ! Transfer variables from GM to RB
    !/
    if(i_proc0(RB_) /= i_proc0(GM_))then

       
       nSize = iSize*jSize*nIntegral
       if(is_proc0(GM_)) then
          call MPI_send(nVarLine,1,MPI_INTEGER,i_proc0(RB_),&
               1,i_comm(),iError)
          call MPI_send(nPointLine,1,MPI_INTEGER,i_proc0(RB_),&
               1,i_comm(),iError)
          call MPI_send(Integral_IIV,nSize,MPI_REAL,&
               i_proc0(RB_),1,i_comm(),iError)
          call MPI_send(BufferLine_VI,nVarLine*nPointLine,MPI_REAL,&
               i_proc0(RB_),2,i_comm(),iError)
       end if
       if(is_proc0(RB_))then
          !setup BufferLine in RB when not sharing proc with GM
          call MPI_recv(nVarLine,1,MPI_INTEGER,i_proc0(GM_),&
               1,i_comm(),iStatus_I,iError)
          call MPI_recv(nPointLine,1,MPI_INTEGER,i_proc0(GM_),&
               1,i_comm(),iStatus_I,iError)
          allocate(BufferLine_VI(nVarLine, nPointLine))
          !recieve variables from GM
          call MPI_recv(Integral_IIV,nSize,MPI_REAL,&
               i_proc0(GM_),1,i_comm(),iStatus_I,iError)
          call MPI_recv(BufferLine_VI,nVarLine*nPointLine,MPI_REAL,&
               i_proc0(GM_),2,i_comm(),iStatus_I,iError)
       end if
    end if

    !\
    ! Transfer satellite names from GM to RB
    !/   
    
    if(nShareSats > 0 .and. i_proc0(RB_) /= i_proc0(GM_))then
       nSize = nShareSats*100
       if(is_proc0(GM_)) then
          call MPI_send(NameSat_I,nSize,MPI_BYTE,i_proc0(RB_),&
               1,i_comm(),iError)
       endif
       if(is_proc0(RB_)) then
          call MPI_recv(NameSat_I,nSize,MPI_BYTE,i_proc0(GM_),&
               1,i_comm(),iStatus_I,iError)
       endif
    
   ! Transfer satellite locations from GM to RB
       
       nSize = 3*2*nShareSats
       if(is_proc0(GM_)) then
          call MPI_send(SatPos_DII,nSize,MPI_REAL,i_proc0(RB_),&
               1,i_comm(),iError)
       endif
       if(is_proc0(RB_)) then
          call MPI_recv(SatPos_DII,nSize,MPI_REAL,i_proc0(GM_),&
               1,i_comm(),iStatus_I,iError)
       endif
    end if
    
    if(DoTest)write(*,*)NameSub,', variables transferred',&
         ', iProc:',iProcWorld
       
    !\
    ! Put variables into RB
    !/
    if(is_proc0(RB_))then
       call RB_put_from_gm(Integral_IIV,iSize,jSize,nIntegral,&
            BufferLine_VI,nVarLine,nPointLine,NameVar,tSimulation)
       if(nShareSats > 0) &
            call RB_put_sat_from_gm(nShareSats, NameSat_I, SatPos_DII)
       if(DoTest) &
            write(*,*)'RB got from GM: RB iProc, Buffer(1,1)=',&
            iProcWorld,Integral_IIV(1,1,:)
    end if

    !\
    ! Deallocate buffer to save memory
    !/
    deallocate(Integral_IIV, BufferLine_VI)

    if (nShareSats > 0) then
       deallocate(NameSat_I)
       deallocate(SatPos_DII)
    end if

    if(DoTest)write(*,*)NameSub,', variables deallocated',&
         ', iProc:',iProcWorld

    if(DoTest)write(*,*)NameSub,' finished, iProc=',iProcWorld

    ! if(DoTest.and.is_proc0(RB_)) call RB_print_variables('GM',iError)

  end subroutine couple_gm_rb

end module CON_couple_gm_rb
