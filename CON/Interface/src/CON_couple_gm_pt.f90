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
  integer:: iGroupCommon,iCompCommon, nProcCommon

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

    ! create group on intersection of GM_ and PT_
    call MPI_Group_intersection(i_group(GM_), i_group(PT_), iGroupCommon, iError)
    call MPI_Group_size(iGroupCommon,nProcCommon,iError)
    call MPI_Group_free(iGroupCommon, iError)

    IsSameLayout = n_proc(GM_)      == n_proc(PT_) &
         .and.     i_proc0(GM_)     == i_proc0(PT_) &
         .and.     i_proc_last(GM_) == i_proc_last(PT_)

    !    if(.not.IsSameLayout)call CON_stop(NameSub// &
    !         ': GM and PT must use the same processors for now!')

    !    if( .not.(IsSameLayout).and.(nProcCommon > 0) )call CON_stop(NameSub// &
    !         ': GM and PT must either share all processors or not share any for now!')

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
    integer, save:: nCoupleGm, nCouplePt

    ! processors of the OTHER component to communicate with
    integer, allocatable, save:: iCoupleProcGm_I(:), iCoupleProcPt_I(:)

    ! number of entries received/sent by a processor during rendezvous
    integer, allocatable, save:: nCouplePointGm_I(:), nCouplePointPt_I(:)
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
            iDecompGm, 1, MPI_INTEGER, iProc0Gm, iCommGmPt, iError)
       call MPI_bcast(&
            iDecompPt, 1, MPI_INTEGER, iProc0Pt, iCommGmPt, iError)
       ! added by Dmitry
       call MPI_bcast(&
            nDim, 1, MPI_INTEGER, iProc0Pt, iCommGmPt, iError)
    endif

    IsNewRoute = iDecompGm /= iDecompLastGm .or. iDecompLastPt /= iDecompPt

    iDecompLastGm = iDecompGm
    iDecompLastPt = iDecompPt

    if(IsNewRoute)then
       nPointPt = 0
       nPointGm = 0
       ! Get positions where info is needed from PT. PT will allocate array.
       if(is_proc(PT_))then
          call PT_put_from_gm(.false., NameVar, nVar, nPointPt, PosPt_DI)
       else
          nPointPt = 0
          allocate(PosPt_DI(nDim,0))
       end if
       call MPI_bcast(&
            nVar, 1, MPI_INTEGER, iProc0Pt, iCommGmPt, iError)

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
          !\
          ! Layouts ovelap
          !/
          if(allocated(iCoupleProcGm_I ))deallocate(iCoupleProcGm_I )
          if(allocated(iCoupleProcPt_I ))deallocate(iCoupleProcPt_I )
          if(allocated(nCouplePointGm_I))deallocate(nCouplePointGm_I)
          if(allocated(nCouplePointPt_I))deallocate(nCouplePointPt_I)
          if(allocated(iPointPt_I      ))deallocate(iPointPt_I      )
          if(allocated(PosGm_DI        ))deallocate(PosGm_DI        )

          call set_inquiry(iCommGmPt,nProcGmPt,iProcGmPt,&
               GM_, PT_, nProcCommon,&
               nCouplePt, iCoupleProcPt_I, nCouplePointPt_I,&
               nCoupleGm, iCoupleProcGm_I, nCouplePointGm_I)
          ! iCoupleProc_I, nCouplePoint_I are allocated

          if(is_proc(PT_))then
             ! \
             ! processor is on PT component
             ! /
             nCouplePointPt_I(1:nCouplePt-1) = nPointPt / nCouplePt 
             nCouplePointPt_I(nCouplePt) = &
                  nPointPt - sum(nCouplePointPt_I(1:nCouplePt-1))
          end if

          ! send number of points
          call get_recv_buffer_size(iCommGmPt, nProcGmPt, iProcGmPt,&
               nCouplePt, iCoupleProcPt_I, nCouplePointPt_I,&
               nCoupleGm, iCoupleProcGm_I, nCouplePointGm_I)
          if(is_proc(GM_))&
               nPointGm = sum(nCouplePointGm_I(:))

          allocate(PosGm_DI(nDim,nPointGm))

          call transfer_buffer_direct(iCommGmPt, nProcGmPt, iProcGmPt,&
               nCouplePt, iCoupleProcPt_I, nCouplePointPt_I,&
               nCoupleGm, iCoupleProcGm_I, nCouplePointGm_I,&
               nDim,&
               nBufferS=nPointPt, BufferS_I=PosPt_DI,&
               nBufferR=nPointGm, BufferR_I=PosGm_DI)

          allocate(iProcGm_I(nPointGm))
          allocate(iProcPt_I(nPointPt))
          ! Find processors that own the PT positions in GM
          if(is_proc(GM_))&
               call GM_find_points(nDim, nPointGm, PosGm_DI, iProcGm_I)
          deallocate(PosGm_DI)
          ! send processors on GM owning positions of data
          call transfer_buffer_direct(iCommGmPt, nProcGmPt, iProcGmPt,&
               nCoupleGm, iCoupleProcGm_I, nCouplePointGm_I,&
               nCouplePt, iCoupleProcPt_I, nCouplePointPt_I,&
               1,&
               nBufferS=nPointGm, iBufferS_I=iProcGm_I,&
               nBufferR=nPointPt, iBufferR_I=iProcPt_I)

          !set_transfer_rndvou
          
          ! Shift from local GM proc index to router proc index
          !if(is_proc(PT_))&
          !     iProcPt_I = iProcPt_I + iProc0Gm               
          call  set_data_transfer(iCommGmPt,&
               GM_, PT_, &
               nPointPt, iProcPt_I, &
               nPointGm, iProcGm_I, &
               nCouplePt, iCoupleProcPt_I, nCouplePointPt_I,&
               nCoupleGm, iCoupleProcGm_I, nCouplePointGm_I)
          allocate(PosSortPt_DI(nDim,nPointPt))

          if(is_proc(PT_))then
             ! sort iProcPt_I according to order of procs in iCoupleProcPt_I
             allocate(iPointPt_I(nPointPt))

             call get_transfer_buffer_order(GM_, PT_, &
                  nPointPt, iProcPt_I, &
                  nCouplePt, iCoupleProcPt_I, nCouplePointPt_I,&
                  iPointPt_I)

             do iPoint = 1, nPointPt
                PosSortPt_DI(:,iPointPt_I(iPoint)) = PosPt_DI(:,iPoint)
             end do
          end if

          deallocate(iProcPt_I)
          deallocate(PosPt_DI)

          deallocate(iProcGm_I)

          if(is_proc(GM_))&
             nPointGm = sum(nCouplePointGm_I(:))

          allocate(PosGm_DI(nDim, nPointGm))

          call transfer_buffer_direct(iCommGmPt, nProcGmPt, iProcGmPt,&
               nCouplePt, iCoupleProcPt_I, nCouplePointPt_I,&
               nCoupleGm, iCoupleProcGm_I, nCouplePointGm_I,&
               nDim,&
               nBufferS=nPointPt, BufferS_I=PosSortPt_DI,&
               nBufferR=nPointGm, BufferR_I=PosGm_DI)

          deallocate(PosSortPt_DI)

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
       ! Layouts overlap
       !/
       nPointGm = sum(nCouplePointGm_I)

       allocate(DataPt_VI(nVar,nPointPt))
       allocate(DataGm_VI(nVar,nPointGm))

       if(is_proc(GM_))&
            call GM_get_for_pt( IsNewRoute, &
            NameVar, nVar, nDim, nPointGm, PosGm_DI, DataGm_VI)

       call transfer_buffer_direct(iCommGmPt, nProcGmPt, iProcGmPt,&
            nCoupleGm, iCoupleProcGm_I, nCouplePointGm_I,&
            nCouplePt, iCoupleProcPt_I, nCouplePointPt_I,&
            nVar,&
            nBufferS=nPointGm, BufferS_I=DataGm_VI,&
            nBufferR=nPointPt, BufferR_I=DataPt_VI)

       deallocate(DataGm_VI)
       
       ! Give the data to PT
       if(is_proc(PT_))&
            call PT_put_from_gm(.true., &
            NameVar, nVar, nPointPt, PosPt_DI, DataPt_VI, iPointPt_I)
       
       deallocate(DataPt_VI)
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

  subroutine set_inquiry(iComm,nProc,iProc,&
       iCompSource, iCompTarget, nProcCommon,&
       nCoupleTarget, iCoupleProcTarget_I, nCouplePointTarget_I,&
       nCoupleSource, iCoupleProcSource_I, nCouplePointSource_I)
    ! the subroutine returns number of processors to be communicated with
    !     nCoupleSource & nCouple Target 
    ! and stores ids of these processors in 
    !     iCoupleProcSource_I & iCoupleProcTarget_I
    ! iCompSource is the component which Data is stored on
    ! iCompTarget is the component which Data is to be sent to
    ! in the Union communicator procs on Source are assumed to come first
    integer,              intent(in ):: iComm, nProc, iProc
    integer,              intent(in ):: iCompSource, iCompTarget
    integer,              intent(in ):: nProcCommon
    !------------------------------------------------------------------
    integer,              intent(out):: nCoupleTarget
    integer, allocatable, intent(out):: iCoupleProcTarget_I(:)
    integer, allocatable, intent(out):: nCouplePointTarget_I(:)
    !------------------------------------------------------------------
    integer,              intent(out):: nCoupleSource
    integer, allocatable, intent(out):: iCoupleProcSource_I(:)
    integer, allocatable, intent(out):: nCouplePointSource_I(:)
    !------------------------------------------------------------------
    integer:: dProc ! remainder of division (nCompSource/nCompTarget)**(+/-1)
    integer:: iProcLocal, iCouple
    integer:: nProcTarget, nProcSource
    integer:: nCouple, nCoupleOther, nProc0Other
    integer, allocatable:: iCoupleProc_I(:)
    integer:: iError
    !------------------------------------------------------------------
    ! Here shared procs are considered to be on Source and not on Target
    ! they are processed separately in order to send inquiry to itself
    nProcSource = n_proc(iCompSource)
    nProcTarget = n_proc(iCompTarget) - nProcCommon

    if(nProcTarget < 1)then
       ! there are no procs on Target only
       if(    is_proc(iCompTarget))then
          nCoupleTarget = 1
          nCoupleSource = 1
       elseif(is_proc(iCompSource))then
          nCoupleTarget = 0
          nCoupleSource = 0
       end if
       allocate(iCoupleProcTarget_I( nCoupleTarget))
       allocate(iCoupleProcSource_I( nCoupleSource))

       iCoupleProcTarget_I(nCoupleTarget) = iProc
       iCoupleProcSource_I(nCoupleSource) = iProc

       allocate(nCouplePointTarget_I(nCoupleTarget))
       allocate(nCouplePointSource_I(nCoupleSource))
       RETURN
    end if

    !\
    ! if nProcSource/nProcTarget is larger and not an exact multiple of the 
    ! other, there is a discrepancy of dProc processors
    !/
    if(nProcSource > nProcTarget)then
       dProc = MOD(nProcSource, nProcTarget)
    else
       dProc = MOD(nProcTarget, nProcSource)
    end if

    if(    is_proc(iCompTarget) .and. .not.(is_proc(iCompSource)))then
       nCouple = nProcSource / nProcTarget
       nCoupleOther   = nProcTarget / nProcSource
       call MPI_Comm_rank(iComm, iProcLocal, iError)
       iProcLocal = iProcLocal - nProcSource
       nProc0Other = 0
    elseif(is_proc(iCompSource))then
       nCouple = nProcTarget / nProcSource
       nCoupleOther   = nProcSource / nProcTarget
       call MPI_Comm_rank(iComm, iProcLocal, iError)
       nProc0Other = nProcSource
    end if

    if(iProcLocal < dProc * (nCoupleOther+1))then
       !\
       ! account for discrepancy:: 
       ! increase nCoupleOther and nCouple by 1
       !/
       dProc = 0
       nCoupleOther = nCoupleOther + 1
       nCouple      = nCouple + 1
    else
       ! change sign of dProc if nCouple==0
       dProc = dProc * (2 * MIN(1,nCouple) - 1)
       nCouple      = MAX(1, nCouple)
       nCoupleOther = MAX(1, nCoupleOther)
    end if

    allocate(iCoupleProc_I( nCouple))

    do iCouple = 1, nCouple
       iCoupleProc_I(iCouple) = &
            nProc0Other + (iCouple - 1) + &
            (nCouple * iProcLocal + dProc)/nCoupleOther
    end do

    if(  is_proc(iCompSource).and.&
         is_proc(iCompTarget))then
       nCoupleTarget = 1
       nCoupleSource = nCouple + 1
       allocate(iCoupleProcTarget_I(nCoupleTarget))
       allocate(iCoupleProcSource_I(nCoupleSource))
       iCoupleProcTarget_I(1            ) = iProc
       iCoupleProcSource_I(nCoupleSource) = iProc
       iCoupleProcSource_I(1:nCoupleSource-1) = iCoupleProc_I(:)
    elseif(is_proc(iCompSource))then
       nCoupleTarget = 0
       nCoupleSource = nCouple
       allocate(iCoupleProcTarget_I(nCoupleTarget))
       allocate(iCoupleProcSource_I(nCoupleSource))
       iCoupleProcSource_I(:) = iCoupleProc_I(:)
    elseif(is_proc(iCompTarget))then
       nCoupleTarget = nCouple
       nCoupleSource = 0
       allocate(iCoupleProcTarget_I(nCoupleTarget))
       allocate(iCoupleProcSource_I(nCoupleSource))
       iCoupleProcTarget_I(:) = iCoupleProc_I(:)
    end if

    deallocate(iCoupleProc_I)

    allocate(nCouplePointTarget_I(nCoupleTarget))
    allocate(nCouplePointSource_I(nCoupleSource))

  end subroutine set_inquiry

  !===========================================================================

  subroutine get_recv_buffer_size(iComm, nProc, iProc,&
       nCoupleS, iCoupleProcS_I, nCouplePointS_I,&
       nCoupleR, iCoupleProcR_I, nCouplePointR_I)
    ! This subroutine determines how large is the buffer
    ! to be transfered to recv component processor
    integer, intent(in)   :: iComm                  ! MPI communicator
    integer, intent(in)   :: nProc                  ! number of processors
    integer, intent(in)   :: iProc                  ! local proc index
    !-----------------------------------------------------------------------
    integer, intent(in):: nCoupleS                 ! 
    integer, intent(in):: iCoupleProcS_I( nCoupleS)! Couple procs
    integer, intent(in):: nCouplePointS_I(nCoupleS)! # of points to send/recv
    !-----------------------------------------------------------------------
    integer, intent(in ):: nCoupleR                 !
    integer, intent(in ):: iCoupleProcR_I( nCoupleR)!
    integer, intent(out):: nCouplePointR_I(nCoupleR)!
    !-----------------------------------------------------------------------
    ! index of recv and send processors
    integer:: iProcR, iProcS 
    integer:: iCouple

    ! MPI stuff
    integer:: iStatus(MPI_STATUS_SIZE)
    integer, parameter:: iTag = 77
    integer:: iError
    integer, allocatable:: iRequest_I(:), iStatus_II(:,:)
    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'get_recv_buffer_size'
    !-----------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*)NameSub, 'called with nProc, iProc=', &
         nProc, iProc
    if(nCoupleS + nCoupleR == 0) RETURN

    allocate(iRequest_I(nCoupleS+nCoupleR))
    allocate(iStatus_II(MPI_STATUS_SIZE,nCoupleS+nCoupleR))
    if(nCoupleS > 0)then
       !\
       ! Send buffer size
       !/
       do iCouple = 1, nCoupleS
          iProcS = iCoupleProcS_I(iCouple)
          call MPI_Isend(&
               nCouplePointS_I(iCouple), 1, MPI_INTEGER, &
               iProcS, iTag, iComm, iRequest_I(iCouple), iError)
       end do
    end if
    if(nCoupleR > 0)then
       !\
       ! Recv buffer size
       !/
       do iCouple = 1, nCoupleR
          iProcR = iCoupleProcR_I(iCouple)
          call MPI_Irecv(&
               nCouplePointR_I(iCouple), 1, MPI_INTEGER, &
               iProcR, iTag, iComm, iRequest_I(nCoupleS+iCouple), iError)
       end do
    end if
    !\
    ! Finalize transfer
    !/
    call MPI_waitall(nCoupleS + nCoupleR, iRequest_I,iStatus_II,iError)
    deallocate(iRequest_I, iStatus_II)
  end subroutine get_recv_buffer_size

  !==========================================================================
  
  subroutine set_data_transfer(iComm, iCompSource, iCompTarget,&
       nBufferTarget, iProcTarget_I, &
       nBufferSource, iProcSource_I, &
       nCoupleTarget, iCoupleProcTarget_I, nCouplePointTarget_I,&
       nCoupleSource, iCoupleProcSource_I, nCouplePointSource_I)
    ! subroutine defines transfer based on processors 
    ! of the OTHER components, which own data items
    ! iProc_I
    integer, intent(in):: iComm
    integer, intent(in):: iCompSource, iCompTarget     ! send, recv components
    integer, intent(in):: nBufferTarget                ! number of points
    integer, intent(in):: iProcTarget_I(nBufferTarget) ! proc indexes
    integer, intent(in):: nBufferSource                ! number of points
    integer, intent(in):: iProcSource_I(nBufferSource) ! proc indexes
    !----------------------------------------------------------------------
    integer, intent(  out)             :: nCoupleTarget
    integer, intent(inout), allocatable:: iCoupleProcTarget_I( :)
    integer, intent(inout), allocatable:: nCouplePointTarget_I(:)
    !----------------------------------------------------------------------
    integer, intent(  out)             :: nCoupleSource
    integer, intent(inout), allocatable:: iCoupleProcSource_I( :)
    integer, intent(inout), allocatable:: nCouplePointSource_I(:)
    !----------------------------------------------------------------------
    integer:: nProcSource, nProcTarget, nProcIntersect
    integer:: iBuffer, iProc, iCouple, iCoupleBuffer 
    integer:: iProcSourceLocal, iProcTargetLocal
    integer:: iError
    integer, allocatable:: nPoint_PP(:,:), nPointRecv_P(:), nPointSend_P(:) 
    integer, allocatable:: iCoupleProc_P(:), nCouplePoint_P(:)
    character(len=*), parameter:: NameSub = 'set_data_transfer'
    !----------------------------------------------------------------------
    nProcTarget    = n_proc(iCompTarget)
    nProcSource    = n_proc(iCompSource)

    if(is_proc(iCompTarget))then
       if(allocated(iCoupleProcTarget_I )) deallocate(iCoupleProcTarget_I)
       if(allocated(nCouplePointTarget_I)) deallocate(nCouplePointTarget_I)
       !\
       ! get nCouple and number of points per proc
       !/
       allocate(iCoupleProc_P(0:nProcSource-1),nCouplePoint_P(0:nProcSource-1))
       nCouplePoint_P = 0

       ! count points per proc on Source
       do iBuffer = 1, nBufferTarget
          iProc = iProcTarget_I(iBuffer)
          nCouplePoint_P(iProc) = nCouplePoint_P(iProc) + 1
       end do

       ! count procs to communicate with
       nCoupleTarget = 0
       do iProc = 0, nProcSource-1
          if(nCouplePoint_P(iProc) > 0)then
             nCoupleTarget = nCoupleTarget + 1 
             iCoupleProc_P( nCoupleTarget-1) = iProc
             nCouplePoint_P(nCoupleTarget-1) = nCouplePoint_P(iProc)
          end if
       end do

       ! return the result
       allocate(iCoupleProcTarget_I( nCoupleTarget))
       allocate(nCouplePointTarget_I(nCoupleTarget))
       iCoupleProcTarget_I( :) = iCoupleProc_P( :)
       nCouplePointTarget_I(:) = nCouplePoint_P(:)
       deallocate(iCoupleProc_P, nCouplePoint_P)
    end if

    if(is_proc(iCompSource))then
       ! nPoint_PP is number of points to sent from Source to Target
       allocate(nPoint_PP(0:nProcSource-1, 0:nProcTarget-1))
       !\
       ! build nPoint_PP
       !/

       ! each proc on Source has only part of the info, store it into nPoint_PP
       nPoint_PP = 0
       iCouple = 1
       iCoupleBuffer = nCouplePointSource_I(iCouple)
       do iBuffer = 1, nBufferSource
          if(iBuffer > iCoupleBuffer)then
             iCouple = iCouple + 1
             iCoupleBuffer = iCoupleBuffer + nCouplePointSource_I(iCouple)
          end if
          ! Source procs come first in the communicator, no need for shift
          iProcSourceLocal = iProcSource_I(iBuffer)
          call get_local_rank(iComm, iCompTarget, &
               iCoupleProcSource_I(iCouple), iProcTargetLocal)
          nPoint_PP(iProcSourceLocal,iProcTargetLocal) = &
               nPoint_PP(iProcSourceLocal,iProcTargetLocal) + 1
       end do

       deallocate(iCoupleProcSource_I,nCouplePointSource_I)

       !\
       ! Reduce the row nPoint_PP(iProcSourceLocal,:)
       !/
       allocate(nPointSend_P(0:nProcTarget-1),nPointRecv_P(0:nProcTarget-1))

       do iProcSourceLocal = 0, nProcSource-1
          nPointSend_P(:) = nPoint_PP(iProcSourceLocal,:)
          call MPI_Reduce(&
               nPointSend_P, nPointRecv_P, nProcTarget, MPI_INTEGER, MPI_SUM, &
               iProcSourceLocal, i_comm(iCompSource), iError)
       end do

       deallocate(nPoint_PP, nPointSend_P)
       allocate(iCoupleProc_P(0:nProcTarget-1),nCouplePoint_P(0:nProcTarget-1))
       !\
       ! get nCouple and number of points per proc
       !/
       nCoupleSource = 0
       do iProcTargetLocal = 0, nProcTarget-1
          if(nPointRecv_P(iProcTargetLocal) > 0)then
             nCoupleSource = nCoupleSource + 1
             call get_union_rank(iComm, iCompTarget, &
                  iProcTargetLocal, iCoupleProc_P(nCoupleSource-1))

             nCouplePoint_P(nCoupleSource-1) = nPointRecv_P(iProcTargetLocal)
          end if
       end do

       allocate(iCoupleProcSource_I( nCoupleSource))
       allocate(nCouplePointSource_I(nCoupleSource))
       iCoupleProcSource_I( :) = iCoupleProc_P( :)
       nCouplePointSource_I(:) = nCouplePoint_P(:)
       deallocate(iCoupleProc_P, nCouplePoint_P,nPointRecv_P)
    end if

  end subroutine set_data_transfer
  
  !==========================================================================

  subroutine get_transfer_buffer_order(&
       iCompSource, iCompTarget, &
       nBuffer, iProc_I, &
       nCouple, iCoupleProc_I, nCouplePoint_I,&
       iBuffer_I)
    ! iBuffer_I returns the order of the points with processor index
    ! corresponding to iCoupleProc_I
    integer, intent(in)::  iCompSource, iCompTarget ! send, recv components
    integer, intent(in)::  nBuffer    ! number of points
    integer, intent(in)::  iProc_I(nBuffer)  ! proc index for each point
    integer, intent(in)::  nCouple
    integer, intent(in)::  iCoupleProc_I(nCouple), nCouplePoint_I(nCouple)
    integer, intent(out):: iBuffer_I(nBuffer)! index for reordering

    integer:: iBuffer, iCouple, iProc, nProc
    integer, allocatable:: iOrder_I(:)
    integer, allocatable:: iCoupleProcInv_P(:)
    character(len=*), parameter:: NameSub = 'get_transfer_buffer_order'
    !----------------------------------------------------------------------
    nProc = n_proc(iCompSource)
    ! Set starting point of chunks in the reordered buffer
    allocate(iOrder_I(nCouple))
    iOrder_I(1) = 1
    do iCouple = 2, nCouple
       iOrder_I(iCouple) = iOrder_I(iCouple-1) + nCouplePoint_I(iCouple-1)
    end do

    ! inversion of iCoupleProc_I, i.e. iCouple corresponding to iProc
    allocate(iCoupleProcInv_P(0:nProc-1))
    iCoupleProcInv_P = -1
    do iCouple = 1, nCouple
       iCoupleProcInv_P(iCoupleProc_I(iCouple)) = iCouple
    end do

    ! Get the order that is consecutive per processor index
    do iBuffer = 1, nBuffer
       ! This buffer belongs to iProc
       iProc = iProc_I(iBuffer)
       ! It will occupy iBuffer_I position
       iBuffer_I(iBuffer) = iOrder_I(iCoupleProcInv_P(iProc))
       ! Jump to next position
       iOrder_I(iCoupleProcInv_P(iProc)) = &
            iOrder_I(iCoupleProcInv_P(iProc)) + 1
    end do

    deallocate(iCoupleProcInv_P, iOrder_I)

  end subroutine get_transfer_buffer_order

  !===========================================================================

  subroutine transfer_buffer_direct(&
       iComm, nProc, iProc,&
       nCoupleS, iCoupleProcS_I, nCouplePointS_I,&
       nCoupleR, iCoupleProcR_I, nCouplePointR_I,&
       nData, &
       nBufferS, BufferS_I, iBufferS_I,&
       nBufferR, BufferR_I, iBufferR_I)
    ! This subroutine transfers data
    ! A processor can both send and recv
    ! Correctness: only buffers of the same type should be present 
    use ModMpi

    integer, intent(in):: iComm  ! MPI communicator
    integer, intent(in):: nProc  ! number of processors
    integer, intent(in):: iProc  ! local proc index
    !------------------------------------------------------------------------
    integer, intent(in):: nCoupleS                 ! #   of procs to send to
    integer, intent(in):: iCoupleProcS_I( nCoupleS)! ids of procs to send to
    integer, intent(in):: nCouplePointS_I(nCoupleS)! # of points  to send
    !------------------------------------------------------------------------
    integer, intent(in):: nCoupleR                 ! #   of procs to recv from
    integer, intent(in):: iCoupleProcR_I( nCoupleR)! ids of procs to recv from
    integer, intent(in):: nCouplePointR_I(nCoupleR)! # of points  to recv
    !------------------------------------------------------------------------
    integer, intent(in   )          :: nData       ! number of items per point
    integer, intent(in   )          :: nBufferS    ! send buffer size
    real,    intent(inout), optional::  BufferS_I(nData*nBufferS) ! send buffer
    integer, intent(inout), optional:: iBufferS_I(nData*nBufferS) ! send buffer
    integer, intent(in   )          :: nBufferR    ! recv buffer size
    real,    intent(inout), optional::  BufferR_I(nData*nBufferR) ! recv buffer
    integer, intent(inout), optional:: iBufferR_I(nData*nBufferR) ! recv buffer

    integer:: iProcR, iProcS ! R/S is proc to Recv from/Send to
    integer:: iCouple, iBuffer
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

    if(nCoupleS + nCoupleR == 0) RETURN

    allocate(iRequest_I(nCoupleS+nCoupleR))
    allocate(iStatus_II(MPI_STATUS_SIZE,nCoupleS+nCoupleR))

    !\
    ! Check correctness
    !/
    if(present(BufferS_I) .and. present(iBufferR_I)&
         .or.&
       present(BufferR_I) .and. present(iBufferS_I))then
       call CON_stop(NameSub// &
         ': Send and Recv are attempted for different variable types')
    end if

    !\
    ! Send buffer
    !/
    if(present(BufferS_I))then
       !\
       ! Transfer array of reals
       !/
       if(nCoupleS > 0)then
          !\
          ! Send array of reals
          !/
          iBuffer = 1
          do iCouple = 1, nCoupleS
             iProcS = iCoupleProcS_I(iCouple)
             call MPI_Isend(&
                  BufferS_I(iBuffer), nData*nCouplePointS_I(iCouple), &
                  MPI_REAL, iProcS, iTag, iComm, &
                  iRequest_I(iCouple),iError)
             iBuffer = iBuffer + nData*nCouplePointS_I(iCouple)
          end do
       end if
       if(nCoupleR > 0)then
          !\
          ! Receive array of reals
          !/
          iBuffer = 1
          do iCouple = 1, nCoupleR
             iProcR = iCoupleProcR_I(iCouple)
             call MPI_Irecv(&
                  BufferR_I(iBuffer), nData*nCouplePointR_I(iCouple),&
                  MPI_REAL, iProcR, iTag, iComm, &
                  iRequest_I(nCoupleS+iCouple),iError)
             iBuffer = iBuffer + nData*nCouplePointR_I(iCouple)
          end do
       end if
    elseif(present(iBufferS_I))then
       !\
       ! Transfer array of integers
       !/
       if(nCoupleS > 0)then
          !\
          ! Send array of integers
          !/
          iBuffer = 1
          do iCouple = 1, nCoupleS
             iProcS = iCoupleProcS_I(iCouple)
             call MPI_Isend(&
                  iBufferS_I(iBuffer), nData*nCouplePointS_I(iCouple), &
                  MPI_INTEGER, iProcS, iTag, iComm, &
                  iRequest_I(iCouple),iError)       
             iBuffer = iBuffer + nData*nCouplePointS_I(iCouple)
          end do
       end if
       if(nCoupleR > 0)then
          !\
          ! Receive array of integers
          !/
          iBuffer = 1
          do iCouple = 1, nCoupleR
             iProcR = iCoupleProcR_I(iCouple)
             call MPI_Irecv(&
                  iBufferR_I(iBuffer), nData*nCouplePointR_I(iCouple),&
                  MPI_INTEGER, iProcR, iTag, iComm, &
                  iRequest_I(nCoupleS+iCouple),iError)
             iBuffer = iBuffer + nData*nCouplePointR_I(iCouple)
          end do
       end if
    end if
    !\
    ! Finalize transfer
    !/
    call MPI_waitall(nCoupleS + nCoupleR, iRequest_I,iStatus_II,iError)
    deallocate(iRequest_I, iStatus_II)
  end subroutine transfer_buffer_direct

  !===========================================================================

  subroutine get_local_rank(iCommUnion, iComp, iProc, iProcLocal)
    ! for processor with rank iProc in the union communicator iCommUnion
    ! returns rank iProcLocal in the group associated with component iComp
    integer, intent(in ) :: iCommUnion, iComp
    integer, intent(in ) :: iProc
    integer, intent(out) :: iProcLocal

    integer:: iGroupUnion
    integer:: iProc_I(1), iProcLocal_I(1)
    integer:: iError
    !----------------------------------------------------------------------
    !\
    ! Create group associated with iCommUnion
    call MPI_Comm_Group(iCommUnion, iGroupUnion, iError)
    !\
    ! Translate rank
    iProc_I = iProc
    call MPI_Group_translate_ranks(&
         iGroupUnion, 1, iProc_I, i_group(iComp), iProcLocal_I, iError)
    iProcLocal = iProcLocal_I(1)
    !\
    ! free union group
    call MPI_Group_free(iGroupUnion, iError)

  end subroutine get_local_rank
  !===========================================================================

  subroutine get_union_rank(iCommUnion, iComp, iProcLocal, iProc)
    ! for processor with rank iProc in the union communicator iCommUnion
    ! returns rank iProcLocal in the group associated with component iComp
    integer, intent(in ) :: iCommUnion, iComp
    integer, intent(in ) :: iProcLocal
    integer, intent(out) :: iProc


    integer:: iGroupUnion
    integer:: iProc_I(1), iProcLocal_I(1)
    integer:: iError
    !----------------------------------------------------------------------
    !\
    ! Create group associated with iCommUnion
    call MPI_Comm_Group(iCommUnion, iGroupUnion, iError)
    !\
    ! Translate rank
    iProcLocal_I = iProcLocal
    call MPI_Group_translate_ranks(&
         i_group(iComp), 1, iProcLocal_I, iGroupUnion, iProc_I, iError)
    iProc = iProc_I(1)
    !\
    ! free union group
    call MPI_Group_free(iGroupUnion, iError)

  end subroutine get_union_rank

end module CON_couple_gm_pt
