!^CMP COPYRIGHT UM
!^CMP FILE GM
!^CMP FILE PT

!BOP
!MODULE: CON_couple_gm_pt - couple GM and PT components
!
!DESCRIPTION:
! Couple GM and PT components both ways. 
!
!INTERFACE:
module CON_couple_gm_pt

  !USES:
  use CON_coupler

  implicit none
  save

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_gm_pt_init ! initialize both couplings
  public :: couple_gm_pt      ! couple GM to PT

  !REVISION HISTORY:
  ! 09/24/2013 G.Toth <gtoth@umich.edu> - initial version
  !EOP

  ! Communicator and logicals to simplify message passing and execution
  logical       :: UseMe = .true., IsInitialized = .false.

  ! Check if the GM and PT processors coincide
  logical:: IsSameLayout

  ! Proc index inside SWMF
  integer:: iProcWorld

  ! Router communicator info
  integer:: iCommGmPt, nProcGmPt, iProcGmPt, iProc0Gm, iProc0Pt

contains

  !BOP =======================================================================
  !IROUTINE: couple_gm_pt_init - initialize GM-PT couplings
  !INTERFACE:
  subroutine couple_gm_pt_init

    integer:: iError

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='couple_gm_pt_init'

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s
    !EOP
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(IsInitialized) RETURN
    IsInitialized = .true.
    
    iProcWorld = i_proc()

    ! Get the union communicator, the UseMe logical, 
    ! and the root proc indexes of the two components in the union comm.
    iProc0Gm = i_proc0(GM_)
    iProc0Pt = i_proc0(PT_)
    call set_router_comm(GM_, PT_, iCommGmPt, UseMe, iProc0Gm, iProc0Pt)

    call MPI_comm_size(iCommGmPt, nProcGmPt, iError)
    call MPI_comm_rank(iCommGmPt, iProcGmPt, iError)

    IsSameLayout = n_proc(GM_)      == n_proc(PT_) &
         .and.     i_proc0(GM_)     == i_proc0(PT_) &
         .and.     i_proc_last(GM_) == i_proc_last(PT_)

    if(.not.IsSameLayout)call CON_stop(NameSub// &
         ': GM and PT must use the same processors for now!')

  end subroutine couple_gm_pt_init

  !BOP =======================================================================
  !IROUTINE: couple_gm_pt - couple GM to PT
  !INTERFACE:
  subroutine couple_gm_pt(tSimulation)

    !INPUT ARGUMENT:
    real, intent(in) :: tSimulation

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Global Magnetosphere       (GM) source\\
    !    Particle Tracker           (PT) target
    !
    ! Send information from GM to PT. 
    !EOP

    ! Stored variables for GM->PT coupler
    !------------------------------------
    ! Number of dimensions and number of variables
    integer, save:: nDim, nVar

    ! List of variables to pass
    character(len=100):: NameVar = 'Bx By Bz'

    ! Grid index
    integer:: iDecompLastGm = -1, iDecompLastPt = -1
    integer:: iGridGm, iDecompGm, iGridPt, iDecompPt

    ! Number of local points for PT and GM cores
    integer:: nPointPt=0, nPointGm=0

    ! Number of points that belong to a given processor of the OTHER component
    integer, allocatable, save:: nPointGm_P(:), nPointPt_P(:)

    ! Permutation of PT points after data is returned
    integer, allocatable, save:: iPointPt_I(:)
    
    ! Point positions local on a GM processor
    real, allocatable, save:: PosGm_DI(:,:)
    !---------------------------------------------------------------
    ! Temporary variables

    ! Is there a need to recalculate the data transfer route?
    logical:: IsNewRoute

    ! Storage for original PT point positions. 
    real, pointer:: PosPt_DI(:,:)

    ! Positions sorted according to the correspongin GM processors
    real, allocatable:: PosSortPt_DI(:,:)

    ! GM processor index for PT points
    integer, allocatable:: iProcPt_I(:)

    ! Buffers for data on GM and PT
    real, allocatable:: DataGm_VI(:,:), DataPt_VI(:,:)

    integer:: iPoint

    logical :: DoTest, DoTestMe

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_gm_pt'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*)NameSub,' starting iProc=',iProcWorld

    ! After everything is initialized exclude PEs which are not involved
    if(.not.UseMe) RETURN

    if(DoTest)write(*,*)NameSub,' starting, iProc=',iProcWorld

    ! Allocate arrays for router
    if(.not.allocated(nPointGm_P)) &
         allocate(nPointGm_P(nProcGmPt), nPointPt_P(nProcGmPt))

    ! Check if there is a need to redo the mapping
    if(IsSameLayout)then
       call GM_get_grid_info(nDim, iGridGm, iDecompGm)
       call PT_get_grid_info(nDim, iGridPt, iDecompPt)
    ! else
    !    if(is_proc0(GM_)) call GM_get_grid_id(nDim, iGridGm, iDecompGm)
    !    if(is_proc0(PT_)) call PT_get_grid_id(nDim, iGridGm, iDecompGm)
    !    call MPI_bcast(iDecompGm from GM root to union group)
    !    call MPI_bcast(iDecompPt from PT root to union group)
    endif

    IsNewRoute = iDecompGm /= iDecompLastGm .or. iDecompLastPt /= iDecompPt

    iDecompLastGm = iDecompGm
    iDecompLastPt = iDecompPt

    if(IsNewRoute)then
       nPointPt = 0
       ! Get positions where info is needed from PT. PT will allocate array.
       if(is_proc(PT_)) call PT_put_from_gm(.false., &
            NameVar, nVar, nPointPt, PosPt_DI, DataPt_VI, iPointPt_I)

       if(IsSameLayout)then
          ! Find processors that own the PT positions in GM
          allocate(iProcPt_I(nPointPt))
          call GM_find_points(nDim, nPointPt, PosPt_DI, iProcPt_I)
       else
          ! Calculate total number of points
          !call MPI_allreduce(nPointPt, nPointAll, 1, MPI_INTEGER, MPI_SUM, &
          !   i_comm(PT_), iError)
          ! Calculate how many points will go to the GM processor
          !nPointGm = nPointAll / n_proc(GM_)
          !if(is_proc0(GM_)) nPointGm = nPointAll - (n_proc(GM_)-1)*nPointGm
          ! allocate(PosGm_DI(nDim,nPointGm))
          ! Transfer PosPt_DI to GM processors using default/previous router
          !call transfer_buffer(iCommGmPt, nProcGmPt, iProcGmPt, nDim, &
          !     nPointPt, nPointPt_P, PosPt_DI, &
          !     nPointGm, nPointGm_P, PosGm_DI)
          !
          ! Find processors that own the PT positions in GM
          !if(is_proc(GM_))then
          !  allocate(iProcGm_I(nPointGm))
          !  call GM_find_pos(nPointGm, PosGm_DI, iProcGm_I)
          !  deallocate(PosGm_DI)
          !endif
          ! Send back the GM proc index info to the PT processors
          !allocate(iProcPt_I(nPointPt))
          !call transfer_buffer(iCommGmPt, nProcGmPt, iProcGmPt, nDim, &
          !     nPointGm, nPointGm_P, iProcGm_I, &
          !     nPointPt, nPointPt_P, iProcPt_I)
          !deallocate(iProcGm_I)
          ! Shift from local GM proc index to router proc index
          !iProcPt_I = iProcGm_I + iProc0Gm
       end if

       call get_buffer_order(nProcGmPt, nPointPt, iProcPt_I, &
            nPointPt_P, iPointPt_I)
       
       deallocate(iProcPt_I)

       ! Could be done in place !!! Ask Valeriy
       allocate(PosSortPt_DI(nDim,nPointPt))
       do iPoint = 1, nPointPt
          PosSortPt_DI(:,iPointPt_I(iPoint)) = PosPt_DI(:,iPoint)
       end do
       deallocate(PosPt_DI)

       call set_recv_info(iCommGmPt, nProcGmPt, iProcGmPt, &
            nPointPt_P, nPointGm, nPointGm_P)

       ! Transfer PT positions to the GM processors that own them
       allocate(PosGm_DI(nDim,nPointGm))
       call transfer_buffer(iCommGmPt, nProcGmPt, iProcGmPt, nDim, &
            nPointPt, nPointPt_P, PosSortPt_DI, &
            nPointGm, nPointGm_P, PosGm_DI)

       deallocate(PosSortPt_DI)

    end if

    ! Get the data from GM
    if(is_proc(GM_)) call GM_get_for_pt( IsNewRoute, &
         NameVar, nVar, nDim, nPointGm, PosGm_DI, DataGm_VI)

    call transfer_buffer(iCommGmPt, nProcGmPt, iProcGmPt, nVar, &
         nPointGm, nPointGm_P, DataGm_VI, &
         nPointPt, nPointPt_P, DataPt_VI) 

    ! Give the data to PT
    if(is_proc(PT_)) call PT_put_from_gm(.true., &
         NameVar, nVar, nPointPt, PosPt_DI, DataPt_VI, iPointPt_I)

    if(DoTest) write(*,*) NameSub,' finished, iProc=',iProcWorld

  end subroutine couple_gm_pt

  !==========================================================================

  subroutine get_buffer_order(nProc, nBuffer, iProc_I, nBuffer_P, iBuffer_I)

    ! nBuffer_P returns the number of points in the buffer belong 
    !    to each processor based on the iProc_I information.
    ! iBuffer_I returns the order of the points with increasing 
    !    processor index.

    integer, intent(in)::  nProc, nBuffer    ! number of procs and points
    integer, intent(in)::  iProc_I(nBuffer)  ! proc index for each point
    integer, intent(out):: nBuffer_P(nProc)  ! number of points per procs
    integer, intent(out):: iBuffer_I(nBuffer)! index for reordering

    integer:: iBuffer, iProc
    integer, allocatable:: iBuffer_P(:)
    !----------------------------------------------------------------------

    ! Buffer chunk sizes for each processor
    nBuffer_P = 0
    do iBuffer = 1, nBuffer
       iProc = iProc_I(iBuffer)
       nBuffer_P(iProc) = nBuffer_P(iProc) + 1
    end do

    ! Set starting point of chunks in the reordered buffer
    allocate(iBuffer_P(0:nProc-1))
    iBuffer_P(0) = 1
    do iProc = 1, nProc-1
       iBuffer_P(iProc) = iBuffer_P(iProc-1) + nBuffer_P(iProc)
    end do

    ! Get the order that is consecutive per processor index
    do iBuffer = 1, nBuffer
       ! This buffer belongs to iProc
       iProc = iProc_I(iBuffer)
       
       ! It will occupy iBufferSort position
       iBuffer_I(iBuffer) = iBuffer_P(iProc)
       
       ! Jump to next position
       iBuffer_P(iProc) = iBuffer_P(iProc) + 1
    end do
    deallocate(iBuffer_P)

  end subroutine get_buffer_order

  !==========================================================================

  subroutine set_recv_info(iComm, nProc, iProc, &
       nBufferS_P, nBufferR, nBufferR_P)

    ! Based on nBufferS_P (number of points sent to each recv proc)
    ! set nBufferR and nBufferR_P (number of points recv from each send proc)
    ! nBufferR = sum(nBufferR_P)

    integer, intent(in)::  iComm, nProc, iProc ! Router communicator
    integer, intent(in)::  nBufferS_P(nProc)   ! Send buffer chunks
    integer, intent(out):: nBufferR            ! Recv buffer size
    integer, intent(out):: nBufferR_P(nProc)   ! Recv buffer chunks

    integer, parameter:: iTag = 76

    integer, allocatable:: nPointRS_PP(:,:), nPoint_P(:)
    integer:: iProcR
    integer:: iError, iStatus_I(MPI_STATUS_SIZE)
    !-----------------------------------------------------------------------

    ! Gather the nBufferS_P information onto the root.
    ! The resulting nPointRS_PP table contains the
    ! the number of points owned by iProcR and iProcS.
    if(iProc == 0)then
       allocate(nPointRS_PP(nProc,nProc), nPoint_P(nProc))
    else
       allocate(nPointRS_PP(1,1))
    end if

    call MPI_gather(nBufferS_P, nProc, MPI_INTEGER, &
         nPointRS_PP, nProc, MPI_INTEGER, 0, iComm, iError)

    ! Local copy for root
    if(iProc==0) nBufferR_P = nPointRS_PP(0,:)

    ! Now send the columns of the table to the recv processors
    do iProcR = 1, nProc-1
       if(iProc==0)then
          nPoint_P = nPointRS_PP(iProcR,:)
          call MPI_send(nPoint_P, nProc, MPI_INTEGER, iProcR, &
               iTag, iComm, iError)
       end if
       if(iProc == iProcR)then
          call MPI_recv(nBufferR_P, nProc, MPI_integer, 0, &
               iTag, iComm, iStatus_I, iError)
       end if
    end do
    if(iProc==0) deallocate(nPoint_P)
    deallocate(nPointRS_PP)

    ! Total number of points owned by this recv proc
    nBufferR = sum(nBufferR_P)

  end subroutine set_recv_info

  !===========================================================================

  subroutine transfer_buffer(iComm, nProc, iProc, nData, &
       nBufferS, nBufferS_P, BufferS_I, &
       nBufferR, nBufferR_P, BufferR_I)

    use ModMpi

    integer, intent(in):: iComm                    ! MPI communicator
    integer, intent(in):: nProc                    ! number of processors
    integer, intent(in):: iProc                    ! local proc index
    integer, intent(in):: nData                    ! number of reals per point
    integer:: nBufferS, nBufferR                   ! total buffer sizes

    integer, intent(in):: nBufferS_P(nProc)        ! send buffer chunks 
    integer, intent(in):: nBufferR_P(nProc)        ! recv buffer chunks

    real, intent(in) :: BufferS_I(nData*nBufferS)  ! send buffer
    real, intent(out):: BufferR_I(nData*nBufferR)  ! recv buffer

    ! index of recv and send processors
    integer:: iProcR, iProcS

    integer:: iBufferS, iBufferR

    ! MPI stuff
    integer, parameter:: iTag = 77
    integer:: iRequestR, iRequestS, iError
    integer, allocatable:: iRequestS_I(:), iRequestR_I(:), iStatus_II(:,:)

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'transfer_buffer'
    !-----------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*)'NameSub called with nProc, iProc, nData=', &
         nProc, iProc, nData

    allocate(iRequestS_I(nProc-1), iRequestR_I(nProc-1), &
         iStatus_II(MPI_STATUS_SIZE,nProc))

    ! Possibly optimize for local copies?

    ! Send the appropriate parts of the send buffer
    iRequestS = 0
    iBufferS  = 1
    do iProcR = 0, nProc - 1
       if(nBufferS_P(iProcR) == 0) CYCLE        ! Skip empty chunks
       iRequestS = iRequestS + 1                ! Count send requests

       ! Post send
       call MPI_isend(BufferS_I(iBufferS), nData*nBufferS_P(iProcR), &
            MPI_REAL, iProcR, iTag, iComm, iRequestS_I(iRequestS), iError)

       ! Jump to the starting point of the next chunk
       iBufferS  = iBufferS + nData*nBufferS_P(iProcR) 
    end do

    ! Recv the appropriate parts from the send processors
    iRequestR = 0
    iBufferR  = 1
    do iProcS = 1, nProc
       if(nBufferR_P(iProcS) == 0) CYCLE            ! Skip empty chunks
       iRequestR = iRequestR + 1                    ! Count recv requests

       ! Post recv
       call MPI_irecv(BufferR_I(iBufferR), nData*nBufferR_P(iProcS), &
            MPI_REAL, iProcS, iTag, iComm, iRequestR_I(iRequestR), iError)

       ! Jump to the starting point of the next chunk
       iBufferR  = iBufferR + nData*nBufferR_P(iProcS)  

    end do

    ! wait for all receives to be completed
    if(iRequestR > 0) &
         call MPI_waitall(iRequestR, iRequestR_I, iStatus_II, iError)

    ! wait for all sends to be completed
    if(iRequestS > 0) &
         call MPI_waitall(iRequestS, iRequestS_I, iStatus_II, iError)

    deallocate(iRequestS_I, iRequestR_I, iStatus_II)

  end subroutine transfer_buffer

end module CON_couple_gm_pt
