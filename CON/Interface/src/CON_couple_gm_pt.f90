! !  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
! !  For more information, see http://csem.engin.umich.edu/tools/swmf
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

    !    if(.not.IsSameLayout)call CON_stop(NameSub// &
    !         ': GM and PT must use the same processors for now!')

    if( .not.(IsSameLayout).and. &
    (nProcGmPt /= n_proc(GM_) + n_proc(PT_)) )call CON_stop(NameSub// &
         ': GM and PT must either share all processors or not share any for now!')

  end subroutine couple_gm_pt_init

  !BOP =======================================================================
  !IROUTINE: couple_gm_pt - couple GM to PT
  !INTERFACE:
  subroutine couple_gm_pt(tSimulation)


    interface
       subroutine PT_put_from_gm(UseData, &
            NameVar, nVar, nPoint, Pos_DI, Data_VI, iPoint_I)

         implicit none

         logical,          intent(in)   :: UseData ! true when data is transferred
         ! false if positions are asked
         character(len=*), intent(inout):: NameVar ! List of variables
         integer,          intent(inout):: nVar    ! Number of variables in Data_VI
         integer,          intent(inout):: nPoint  ! Number of points in Pos_DI

         real, pointer:: Pos_DI(:,:)               ! Position vectors

         real,    intent(in), optional:: Data_VI(nVar,nPoint)! Recv data array
         integer, intent(in), optional:: iPoint_I(nPoint)    ! Order of data
       end subroutine PT_put_from_gm
       subroutine GM_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

         implicit none

         integer, intent(in) :: nDimIn                ! dimension of position vectors
         integer, intent(in) :: nPoint                ! number of positions
         real,    intent(in) :: Xyz_DI(nDimIn,nPoint) ! positions
         integer, intent(out):: iProc_I(nPoint)       ! processor owning position

       end subroutine GM_find_points
    end interface
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
    character(len=100):: NameVar = 'rho ux uy uz bx by bz p'

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
    integer, allocatable:: iProcGm_I(:)

    ! Buffers for data on GM and PT
    real, allocatable:: DataGm_VI(:,:), DataPt_VI(:,:)

    integer:: iPoint, i, iError

    logical :: DoTest, DoTestMe

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_gm_pt'
    !-------------------------------------------------------------------------
    ! Variables for general-case coupling

    ! "Rendezvous": 

    ! number of processors of the OTHER component to communicate with
    integer:: nRndvouGm, nRndvouPt

    ! processors of the OTHER component to communicate with
    integer, allocatable, save:: nRndvouProcGm_I(:), nRndvouProcPt_I(:)

    ! number of entries received/sent by a processor during rendezvous
    integer, allocatable, save:: nRndvouPointGm_I(:), nRndvouPointPt_I(:)
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*)NameSub,' starting iProc=',iProcWorld

    ! After everything is initialized exclude PEs which are not involved
    if(.not.UseMe) RETURN

    if(DoTest)write(*,*)NameSub,' starting, iProc=',iProcWorld

    ! Allocate arrays for router
    if(.not.allocated(nPointGm_P) .and. IsSameLayout) then
       allocate(nPointGm_P(nProcGmPt), nPointPt_P(nProcGmPt))
    elseif(.not.allocated(nPointGm_P)) then
       allocate(nPointGm_P(n_proc(PT_)), nPointPt_P(n_proc(GM_)))
    end if

    ! Check if there is a need to redo the mapping
    if(IsSameLayout)then
       call GM_get_grid_info(nDim, iGridGm, iDecompGm)
       call PT_get_grid_info(nDim, iGridPt, iDecompPt)
    else
       ! added by Dmitry
       if(is_proc0(GM_)) call GM_get_grid_info(nDim, iGridGm, iDecompGm)
       if(is_proc0(PT_)) call PT_get_grid_info(nDim, iGridGm, iDecompGm)
       call MPI_bcast(&
            iDecompGm, 1, MPI_INTEGER, i_proc0(GM_), iCommGmPt, iError)
       call MPI_bcast(&
            iDecompPt, 1, MPI_INTEGER, i_proc0(PT_), iCommGmPt, iError)
       ! added by Dmitry
       call MPI_bcast(&
            nDim, 1, MPI_INTEGER, i_proc0(PT_), iCommGmPt, iError)
    endif

    IsNewRoute = iDecompGm /= iDecompLastGm .or. iDecompLastPt /= iDecompPt

    iDecompLastGm = iDecompGm
    iDecompLastPt = iDecompPt


    if(IsNewRoute)then
       nPointPt = 0
       ! Get positions where info is needed from PT. PT will allocate array.
       if(is_proc(PT_)) call PT_put_from_gm(.false., &
            NameVar, nVar, nPointPt, PosPt_DI)
       call MPI_bcast(&
            nVar, 1, MPI_INTEGER, i_proc0(PT_), iCommGmPt, iError)

       if(IsSameLayout)then

          ! Find processors that own the PT positions in GM
          allocate(iProcPt_I(nPointPt))

          call GM_find_points(nDim, nPointPt, PosPt_DI, iProcPt_I)

          allocate(iPointPt_I(nPointPt))
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

       else
          if(is_proc(PT_))then
             ! \
             ! processor is on PT component
             ! /
             if(allocated(nRndvouProcPt_I))&
                  deallocate(nRndvouProcPt_I)
             if(allocated(nRndvouPointPt_I))&
                  deallocate(nRndvouPointPt_I)
             if(allocated(iPointPt_I))&
                  deallocate(iPointPt_I)

             call set_inquiry_rndvou(GM_, PT_, &
                  nRndvouPt, nRndvouProcPt_I, nRndvouPointPt_I)
             ! nRndvouProcPt_I, nRndvouPointPt_I are allocated

             nRndvouPointPt_I(1:nRndvouPt-1) = nPointPt / nRndvouPt 
             nRndvouPointPt_I(nRndvouPt) = &
                  nPointPt - sum(nRndvouPointPt_I(1:nRndvoupT-1))

             ! send number of points
             call get_recv_buffer_size(iCommGmPt, nProcGmPt, iProcGmPt,&
                  PT_, GM_,&
                  nRndvouPt, nRndvouProcPt_I, nRndvouPointPt_I)
             call transfer_buffer_direct(iCommGmPt, nProcGmPt, iProcGmPt,&
                  PT_, GM_,&
                  nRndvouPt, nRndvouProcPt_I, nRndvouPointPt_I,&
                  nDim, nPointPt, PosPt_DI)

             ! receive processors on GM owning positions of data
             allocate(iProcPt_I(nPointPt))
             call transfer_buffer_direct(iCommGmPt, nProcGmPt, iProcGmPt,&
                  GM_, PT_,&
                  nRndvouPt, nRndvouProcPt_I, nRndvouPointPt_I,&
                  1, nPointPt, iBuffer_I=iProcPt_I)!, .true.)
             ! Shift from local GM proc index to router proc index
             iProcPt_I = iProcPt_I + iProc0Gm

             ! set rendezvous to transfer data
             call set_transfer_rndvou(GM_, PT_, nPointPt, iProcPt_I, &
                  nRndvouPt, nRndvouProcPt_I, nRndvouPointPt_I)

             ! sort iProcPt_I according to order of procs in nRndvouProcPt_I
             allocate(iPointPt_I(nPointPt))

             call get_rndvou_buffer_order(GM_, PT_, &
                  nPointPt, iProcPt_I, &
                  nRndvouPt, nRndvouProcPt_I, nRndvouPointPt_I,&
                  iPointPt_I)

             deallocate(iProcPt_I)

             allocate(PosSortPt_DI(nDim,nPointPt))
             do iPoint = 1, nPointPt
                PosSortPt_DI(:,iPointPt_I(iPoint)) = PosPt_DI(:,iPoint)
             end do
             deallocate(PosPt_DI)


             call transfer_buffer_direct(iCommGmPt, nProcGmPt, iProcGmPt,&
                  PT_, GM_,&
                  nRndvouPt, nRndvouProcPt_I, nRndvouPointPt_I,&
                  nDim, nPointPt, PosSortPt_DI)

             deallocate(PosSortPt_DI)



          elseif(is_proc(GM_))then
             ! \
             ! processor is on GM component
             ! /
             if(allocated(nRndvouProcGm_I))&
                  deallocate(nRndvouProcGm_I)
             if(allocated(nRndvouPointGm_I))&
                  deallocate(nRndvouPointGm_I)
             if(allocated(PosGm_DI))&
                  deallocate(PosGm_DI)

             call set_inquiry_rndvou(GM_, PT_, &
                  nRndvouGm, nRndvouProcGm_I, nRndvouPointGm_I)
             ! nRndvouProcGm_I, nRndvouPointGm_I are allocated

             ! find number of points to be received
             call get_recv_buffer_size(iCommGmPt, nProcGmPt, iProcGmPt,&
                  PT_, GM_,&
                  nRndvouGm, nRndvouProcGm_I, nRndvouPointGm_I)
             nPointGm = sum(nRndvouPointGm_I(:))
             allocate(PosGm_DI(nDim,nPointGm))
             call transfer_buffer_direct(iCommGmPt, nProcGmPt, iProcGmPt,&
                  PT_, GM_,&
                  nRndvouGm, nRndvouProcGm_I, nRndvouPointGm_I,&
                  nDim, nPointGm, PosGm_DI)

             ! Find processors that own the PT positions in GM
             allocate(iProcGm_I(nPointGm))
             call GM_find_points(nDim, nPointGm, PosGm_DI, iProcGm_I)
             deallocate(PosGm_DI)
             ! send processors on GM owning positions of data
             call transfer_buffer_direct(iCommGmPt, nProcGmPt, iProcGmPt,&
                  GM_, PT_,&
                  nRndvouGm, nRndvouProcGm_I, nRndvouPointGm_I,&
                  1, nPointGm, iBuffer_I=iProcGm_I)!, .true.)

             ! set rendezvous to transfer data
             call set_transfer_rndvou(GM_, PT_, nPointGm, iProcGm_I, &
                  nRndvouGm, nRndvouProcGm_I, nRndvouPointGm_I)

             deallocate(iProcGm_I)

             nPointGm = sum(nRndvouPointGm_I(:))
             allocate(PosGm_DI(nDim, nPointGm))

             call transfer_buffer_direct(iCommGmPt, nProcGmPt, iProcGmPt,&
                  PT_, GM_,&
                  nRndvouGm, nRndvouProcGm_I, nRndvouPointGm_I,&
                  nDim, nPointGm, PosGm_DI)
          end if
       end if
    end if

    if(IsSameLayout)then
       ! Get the data from GM
       allocate(DataGm_VI(nVar,nPointGm))
       if(is_proc(GM_)) call GM_get_for_pt( IsNewRoute, &
            NameVar, nVar, nDim, nPointGm, PosGm_DI, DataGm_VI)

       allocate(DataPt_VI(nVar,nPointPt))
       call transfer_buffer(iCommGmPt, nProcGmPt, iProcGmPt, nVar, &
            nPointGm, nPointGm_P, DataGm_VI, &
            nPointPt, nPointPt_P, DataPt_VI) 

       ! Give the data to PT
       if(is_proc(PT_)) call PT_put_from_gm(.true., &
            NameVar, nVar, nPointPt, PosPt_DI, DataPt_VI, iPointPt_I)

       deallocate(DataPt_VI)
       deallocate(DataGm_VI)
    else
       !\
       ! Different layouts for components
       !/
       if(is_proc(PT_))then

          allocate(DataPt_VI(nVar,nPointPt))
          call transfer_buffer_direct(iCommGmPt, nProcGmPt, iProcGmPt,&
               GM_, PT_,&
               nRndvouPt, nRndvouProcPt_I, nRndvouPointPt_I,&
               nVar, nPointPt, DataPt_VI)

          ! Give the data to PT
          call PT_put_from_gm(.true., &
               NameVar, nVar, nPointPt, PosPt_DI, DataPt_VI, iPointPt_I)

          deallocate(DataPt_VI)

       elseif(is_proc(GM_))then

          nPointGm = sum(nRndvouPointGm_I)
          allocate(DataGm_VI(nVar,nPointGm))
          call GM_get_for_pt( IsNewRoute, &
               NameVar, nVar, nDim, nPointGm, PosGm_DI, DataGm_VI)

          call transfer_buffer_direct(iCommGmPt, nProcGmPt, iProcGmPt,&
               GM_, PT_,&
               nRndvouGm, nRndvouProcGm_I, nRndvouPointGm_I,&
               nVar, nPointGm, DataGm_VI)

          deallocate(DataGm_VI)
       end if
    end if

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
    integer, intent(out):: nBuffer_P(0:nProc-1)  ! number of points per procs
    integer, intent(out):: iBuffer_I(nBuffer)! index for reordering

    integer:: iBuffer, iProc
    integer, allocatable:: iBuffer_P(:)
    character(len=*), parameter:: NameSub = 'get_buffer_order'
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
       iBuffer_P(iProc) = iBuffer_P(iProc-1) + nBuffer_P(iProc-1)
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
       allocate(nPointRS_PP(0:nProc-1,0:nProc-1), nPoint_P(0:nProc-1))
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

    integer, intent(in):: nBufferS_P(0:nProc-1)        ! send buffer chunks 
    integer, intent(in):: nBufferR_P(0:nProc-1)        ! recv buffer chunks

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
    if(DoTestMe)write(*,*)NameSub, 'called with nProc, iProc, nData=', &
         nProc, iProc, nData

    allocate(iRequestS_I(nProc), iRequestR_I(nProc), iStatus_II(MPI_STATUS_SIZE,nProc))

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
    do iProcS = 0, nProc - 1
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

  !===========================================================================

  subroutine set_inquiry_rndvou(nCompSource, nCompTarget,&
       nRndvou, nRndvouProc_I, nRndvouPoint_I)
    ! the subroutine returns number of processors to be communicated with
    !     nRndvou
    ! and stores numbers of these processors in 
    !     nRndvouProc_I
    ! this array is ALLOCATED INSIDE this subroutine
    ! nCompSource is the component which Data is stored on
    ! nCompTarget is the component which Data is to be sent to
    integer, intent(in):: nCompSource, nCompTarget
    integer, intent(out):: nRndvou
    integer, allocatable, intent(out):: nRndvouProc_I(:)
    integer, allocatable, intent(out):: nRndvouPoint_I(:)

    integer:: dProc ! remainder of division (nCompSource/nCompTarget)**(+/-1)
    integer:: iProcLocal
    integer:: iRndvou ! loop variable
    integer:: nProcTarget, nProcSource
    integer:: nRndvouOther, nProc0Other, nProcStrideOther
    !------------------------------------------------------------------
    nProcSource = n_proc(nCompSource)
    nProcTarget = n_proc(nCompTarget)

    !\
    ! if nProcSource/nProcTarget is larger and not an exact multiple of the 
    ! other, there is a discrepancy of dProc processors
    !/
    if(nProcSource > nProcTarget)then
       dProc = MOD(nProcSource, nProcTarget)
    else
       dProc = MOD(nProcTarget, nProcSource)
    end if
    if(    is_proc(nCompTarget))then
       nRndvou = nProcSource / nProcTarget
       nRndvouOther   = nProcTarget / nProcSource
       iProcLocal = i_proc(nCompTarget)
       nProc0Other = i_proc0(nCompSource)
       nProcStrideOther = i_proc_stride(nCompSource)
    elseif(is_proc(nCompSource))then
       nRndvou = nProcTarget / nProcSource
       nRndvouOther   = nProcSource / nProcTarget
       iProcLocal = i_proc(nCompSource)
       nProc0Other = i_proc0(nCompTarget)
       nProcStrideOther = i_proc_stride(nCompTarget)
    end if

    !\
    ! account for discrepancy: increase nRndvou of first dProc by 1
    !/
    if(nRndvou > 0 .and. dProc > 0 .and. iProcLocal < dProc)&
         nRndvou = nRndvou + 1

    allocate(nRndvouProc_I(max(1,nRndvou)))

    if(nRndvou == 0)then
       !\
       ! there are more procs on this component than on the other
       ! => only 1 rendezvous
       !/
       nRndvou = 1
       if(iProcLocal < dProc * (nRndvouOther+1))then
          !\
          ! account for discrepancy:: increase nRndvouOther of first dProc by 1
          !/
          nRndvouProc_I(nRndvou) = nProc0Other + &
               (iProcLocal / (nRndvouOther + 1)) * nProcStrideOther
       else
          nRndvouProc_I(nRndvou) = nProc0Other + &
               (dProc + (iProcLocal - dProc * (nRndvouOther+1)) / nRndvouOther)&
               * nProcStrideOther
       endif
    else
       if(iProcLocal < dProc)then
          do iRndvou = 1, nRndvou
             nRndvouProc_I(iRndvou) = nProc0Other + &
                  (nRndvou * iProcLocal + &
                  iRndvou - 1) * nProcStrideOther
          end do
       else
          do iRndvou = 1, nRndvou
             nRndvouProc_I(iRndvou) = nProc0Other + &
                  ((nRndvou+1) * dProc + &
                  nRndvou * (iProcLocal - dProc) + &
                  iRndvou - 1) * nProcStrideOther
          end do
       end if
    end if

    allocate(nRndvouPoint_I(nRndvou))

!write(*,*)'iProc = ', i_proc(), 'Rndvou = ', nRndvouProc_I

  end subroutine set_inquiry_rndvou

  !===========================================================================

  subroutine get_recv_buffer_size(iComm, nProc, iProc,&
       nCompSend, nCompRecv,&
       nRndvou, nRndvouProc_I, nRndvouPoint_I)
    ! This subroutine determines how large is the buffer
    ! to be transfered to recv component processor
    integer, intent(in)   :: iComm                ! MPI communicator
    integer, intent(in)   :: nProc                ! number of processors
    integer, intent(in)   :: iProc                ! local proc index
    integer, intent(in)   :: nCompSend, nCompRecv ! send, recv components
    integer, intent(in)   :: nRndvou              ! # of Rndvou procs
    integer, intent(in)   :: nRndvouProc_I( nRndvou)! Rndvou procs
    integer, intent(inout):: nRndvouPoint_I(nRndvou)! # of points to send/recv

    ! index of recv and send processors
    integer:: iProcRecv, iProcSend 
    integer:: iRndvou

    ! MPI stuff
    integer:: iStatus(MPI_STATUS_SIZE)
    integer, parameter:: iTag = 77
    integer:: iError
    integer, allocatable:: iRequest_I(:), iStatus_II(:,:)
    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'get_recv_buffer_size'
    !-----------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*)NameSub, 'called with nProc, iProc', &
         nProc, iProc

    allocate(iRequest_I(nRndvou), iStatus_II(MPI_STATUS_SIZE,nRndvou))
    if(is_proc(nCompSend)) then
       do iRndvou = 1, nRndvou
          iProcRecv = nRndvouProc_I(iRndvou)
          ! Send size of buffer to be sent
          call MPI_Isend(&
               nRndvouPoint_I(iRndvou), 1, MPI_INTEGER, &
               iProcRecv, iTag, iComm, iRequest_I(iRndvou), iError)
       end do
    elseif(is_proc(nCompRecv))then
       do iRndvou = 1, nRndvou
          iProcSend = nRndvouProc_I(iRndvou)
          ! Recv size of buffer to be recved
          call MPI_Irecv(&
               nRndvouPoint_I(iRndvou), 1, MPI_INTEGER, &
               iProcSend, iTag, iComm, iRequest_I(iRndvou), iError)
       end do
    end if

    call MPI_waitall(nRndvou, iRequest_I, iStatus_II, iError)
    deallocate(iRequest_I, iStatus_II)

  end subroutine get_recv_buffer_size

  !==========================================================================

  subroutine set_transfer_rndvou(&
       nCompSource, nCompTarget,&       
       nBuffer, iProc_I, &
       nRndvou, nRndvouProc_I, nRndvouPoint_I)
    ! subroutine defines rendezvous based on processrs 
    ! of the OTHER components, which own data items
    ! iProc_I
    integer, intent(in)::  nCompSource, nCompTarget ! send, recv components
    integer, intent(in)::  nBuffer    ! number of procs and points
    integer, intent(in)::  iProc_I(nBuffer)  ! proc index for each point
    integer, intent(out):: nRndvou
    integer, intent(inout), allocatable:: nRndvouProc_I(:)
    integer, intent(inout), allocatable:: nRndvouPoint_I(:)

    integer:: nProcSource, nProcTarget
    integer:: iBuffer, iProc, iProcSourceLocal, iProcTargetLocal, iRndvou, iRndvouBuffer
    integer:: iError
    integer, allocatable:: nPoint_PP(:,:), nPointRecv_P(:), nPointSend_P(:) 
    integer, allocatable:: iRndvouProc_I(:), iRndvouPoint_I(:)
    character(len=*), parameter:: NameSub = 'set_rndvou_target'
    !----------------------------------------------------------------------
    nProcTarget = n_proc(nCompTarget)
    nProcSource = n_proc(nCompSource)

    if(is_proc(nCompTarget))then
       if(allocated(nRndvouProc_I))&
            deallocate(nRndvouProc_I)
       if(allocated(nRndvouPoint_I))&
            deallocate(nRndvouPoint_I)
       ! get nRndvou and number of points per proc
       allocate(iRndvouProc_I(0:nProcSource-1),iRndvouPoint_I(0:nProcSource-1))
       iRndvouPoint_I = 0

       do iBuffer = 1, nBuffer
          iProc = iProc_I(iBuffer)
          iRndvouPoint_I(iProc) = iRndvouPoint_I(iProc) + 1
       end do

       nRndvou = 0
       do iProc = 0, nProcSource-1
          if(iRndvouPoint_I(iProc) > 0)then
             nRndvou = nRndvou + 1 
             iRndvouProc_I( nRndvou-1) = iProc
             iRndvouPoint_I(nRndvou-1) = iRndvouPoint_I(iProc)
          end if
       end do

    elseif(is_proc(nCompSource))then
       ! array which contains number of pointst to be sent from 
       ! source component to target
       allocate(nPoint_PP(0:nProcSource-1, 0:nProcTarget-1))
       nPoint_PP = 0

       iRndvou = 1
       iRndvouBuffer = nRndvouPoint_I(iRndvou)
       do iBuffer = 1, nBuffer
          if(iBuffer > iRndvouBuffer)then
             iRndvou = iRndvou + 1
             iRndvouBuffer = iRndvouBuffer + nRndvouPoint_I(iRndvou)
          end if
          iProcSourceLocal = (iProc_I(iBuffer) - i_proc0(nCompSource)) / &
               i_proc_stride(nCompSource)
          iProcTargetLocal = (nRndvouProc_I(iRndvou) - i_proc0(nCompTarget)) / &
               i_proc_stride(nCompTarget)
          nPoint_PP(iProcSourceLocal,iProcTargetLocal) = &
               nPoint_PP(iProcSourceLocal,iProcTargetLocal) + 1
       end do

       deallocate(nRndvouProc_I,nRndvouPoint_I)

       allocate(nPointSend_P(0:nProcTarget-1),nPointRecv_P(0:nProcTarget-1))

       do iProcSourceLocal = 0, nProcSource-1
          nPointSend_P(:) = nPoint_PP(iProcSourceLocal,:)
          call MPI_Reduce(nPointSend_P, nPointRecv_P, nProcTarget,&
               MPI_INTEGER, MPI_SUM, iProcSourceLocal, i_comm(nCompSource), iError)
       end do

       deallocate(nPoint_PP, nPointSend_P)
       allocate(iRndvouProc_I(0:nProcTarget-1),iRndvouPoint_I(0:nProcTarget-1))

       nRndvou = 0
       do iProcTargetLocal = 0, nProcTarget-1
          if(nPointRecv_P(iProcTargetLocal) > 0)then
             nRndvou = nRndvou + 1
             iRndvouProc_I( nRndvou-1) = i_proc0(nCompTarget) + &
                  iProcTargetLocal * i_proc_stride(nCompTarget)
             iRndvouPoint_I(nRndvou-1) = nPointRecv_P(iProcTargetLocal)
          end if
       end do

       deallocate(nPointRecv_P)

    end if

       allocate(nRndvouProc_I(nRndvou), nRndvouPoint_I(nRndvou))
       nRndvouProc_I( 1:nRndvou) = iRndvouProc_I( 0:nRndvou-1)
       nRndvouPoint_I(1:nRndvou) = iRndvouPoint_I(0:nRndvou-1)
       deallocate(iRndvouProc_I, iRndvouPoint_I)

  end subroutine set_transfer_rndvou
  
  !==========================================================================

  subroutine get_rndvou_buffer_order(&
       nCompSource, nCompTarget, &
       nBuffer, iProc_I, &
       nRndvou, nRndvouProc_I, nRndvouPoint_I,&
       iBuffer_I)
    ! iBuffer_I returns the order of the points with processor index
    ! corresponding to nRndvouProc_I
    integer, intent(in)::  nCompSource, nCompTarget ! send, recv components
    integer, intent(in)::  nBuffer    ! number of points
    integer, intent(in)::  iProc_I(nBuffer)  ! proc index for each point
    integer, intent(in)::  nRndvou
    integer, intent(in)::  nRndvouProc_I(nRndvou), nRndvouPoint_I(nRndvou)
    integer, intent(out):: iBuffer_I(nBuffer)! index for reordering

    integer:: iBuffer, iRndvou, iProc, nProc
    integer, allocatable:: iOrder_I(:)
    integer, allocatable:: iRndvouProcInv_P(:)
    character(len=*), parameter:: NameSub = 'get_rndvou_buffer_order'
    !----------------------------------------------------------------------
    nProc = n_proc(nCompSource)
    ! Set starting point of chunks in the reordered buffer
    allocate(iOrder_I(nRndvou))
    iOrder_I(1) = 1
    do iRndvou = 2, nRndvou
       iOrder_I(iRndvou) = iOrder_I(iRndvou-1) + nRndvouPoint_I(iRndvou-1)
    end do

    ! inversion of nRndvouProc_I, i.e. iRndvou corresponding to iProc
    allocate(iRndvouProcInv_P(0:nProc-1))
    iRndvouProcInv_P = -1
    do iRndvou = 1, nRndvou
       iRndvouProcInv_P(nRndvouProc_I(iRndvou)) = iRndvou
    end do

    ! Get the order that is consecutive per processor index
    do iBuffer = 1, nBuffer
       ! This buffer belongs to iProc
       iProc = iProc_I(iBuffer)
       ! It will occupy iBuffer_I position
       iBuffer_I(iBuffer) = iOrder_I(iRndvouProcInv_P(iProc))
       ! Jump to next position
       iOrder_I(iRndvouProcInv_P(iProc)) = &
            iOrder_I(iRndvouProcInv_P(iProc)) + 1
    end do

    deallocate(iRndvouProcInv_P, iOrder_I)

  end subroutine get_rndvou_buffer_order

  !===========================================================================

  subroutine transfer_buffer_direct(iComm, nProc, iProc,&
       nCompSend, nCompRecv,&
       nRndvou, nRndvouProc_I, nRndvouPoint_I,&
       nData, nBuffer, Buffer_I, iBuffer_I)
    ! This subroutine transfers data from send component to recv component
    use ModMpi

    integer, intent(in):: iComm                ! MPI communicator
    integer, intent(in):: nProc                ! number of processors
    integer, intent(in):: iProc                ! local proc index
    integer, intent(in):: nCompSend, nCompRecv ! send and recv components
    integer, intent(in):: nRndvou              ! # of Rndvou procs
    integer, intent(in):: nRndvouProc_I( nRndvou)! Rndvou procs
    integer, intent(in):: nRndvouPoint_I(nRndvou)! # of points to send/recv    

    integer, intent(in):: nData, nBuffer       ! number of reals per point
    real,    intent(inout), optional:: Buffer_I(nData*nBuffer)
    integer, intent(inout), optional:: iBuffer_I(nData*nBuffer)

    integer:: iProcRecv, iProcSend ! index of recv and send processors
    integer:: iRndvou
    integer:: iBuffer

    ! MPI stuff
    integer, parameter:: iTag = 77
    integer:: iError
    integer, allocatable:: iRequest_I(:), iStatus_II(:,:)
    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'transfer_buffer_direct'
    !-----------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*)NameSub, 'called with nProc, iProc, nData=', &
         nProc, iProc, nData
    allocate(iRequest_I(nRndvou), iStatus_II(MPI_STATUS_SIZE,nRndvou))
    iBuffer = 1
    if(present(Buffer_I))then
       !\
       ! Transfer array of reals
       !/
       if(is_proc(nCompSend)) then
          do iRndvou = 1, nRndvou
             iProcRecv = nRndvouProc_I(iRndvou)
             call MPI_Isend(Buffer_I(iBuffer), nData*nRndvouPoint_I(iRndvou), &
                  MPI_REAL, iProcRecv, iTag, iComm, iRequest_I(iRndvou),iError)       
             iBuffer = iBuffer + nData*nRndvouPoint_I(iRndvou)
          end do
       elseif(is_proc(nCompRecv)) then
          do iRndvou = 1, nRndvou
             iProcSend = nRndvouProc_I(iRndvou)
             call MPI_Irecv(Buffer_I(iBuffer), nData*nRndvouPoint_I(iRndvou),&
                  MPI_REAL, iProcSend, iTag, iComm, iRequest_I(iRndvou),iError)
             iBuffer = iBuffer + nData*nRndvouPoint_I(iRndvou)
          end do
       end if
    elseif(present(iBuffer_I))then
       !\
       ! Transfer array of integers
       !/
       if(is_proc(nCompSend)) then
          do iRndvou = 1, nRndvou
             iProcRecv = nRndvouProc_I(iRndvou)
             call MPI_Isend(iBuffer_I(iBuffer), nData*nRndvouPoint_I(iRndvou), &
                  MPI_INTEGER, iProcRecv, iTag, iComm, iRequest_I(iRndvou),iError)       
             iBuffer = iBuffer + nData*nRndvouPoint_I(iRndvou)
          end do
       elseif(is_proc(nCompRecv)) then
          do iRndvou = 1, nRndvou
             iProcSend = nRndvouProc_I(iRndvou)
             call MPI_Irecv(iBuffer_I(iBuffer), nData*nRndvouPoint_I(iRndvou),&
                  MPI_INTEGER, iProcSend, iTag, iComm, iRequest_I(iRndvou), iError)
             iBuffer = iBuffer + nData*nRndvouPoint_I(iRndvou)
          end do
       end if
    end if
    !\
    ! Finalize transfer
    !/
    call MPI_waitall(nRndvou, iRequest_I,iStatus_II,iError)
    deallocate(iRequest_I, iStatus_II)
  end subroutine transfer_buffer_direct

end module CON_couple_gm_pt
