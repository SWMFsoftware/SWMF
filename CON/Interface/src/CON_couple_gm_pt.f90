!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
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
  integer:: iCommGmPt, nProcGmPt, iProcGmPt, iProc0Gm, iProc0Pt, nProcCommon
  
contains

  !BOP =======================================================================
  !IROUTINE: couple_gm_pt_init - initialize GM-PT couplings
  !INTERFACE:
  subroutine couple_gm_pt_init

    integer:: iError
    !integer:: nDimPt

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
    nProcCommon = n_proc(GM_) + n_proc(PT_) - nProcGmPt

    IsSameLayout = n_proc(GM_)      == n_proc(PT_) &
         .and.     i_proc0(GM_)     == i_proc0(PT_) &
         .and.     i_proc_last(GM_) == i_proc_last(PT_)

    !! Check if nDim in GM is the same as nDimPt
    !! nDim info could be stored in Grid_C!
    !if(is_proc(GM_)) call GM_get_grid_info(nDim, iGridGm, iDecompGm)
    !if(.not.IsSameLayout) call MPI_bcast(&
    !     nDim, 1, MPI_INTEGER, iProc0Pt, iCommGmPt, iError)
    !if(is_proc0(PT_))then
    !   call PT_get_grid_info(nDimPt, iGridPt, iDecompPt)
    !   if(nDim /= nDimPt)then
    !      write(*,*) NameSub,' ERROR: nDim, nDimPt=', nDim, nDimPt
    !      call CON_stop(NameSub// &
    !          ': GM and PT have different number of spatial dimensions')
    !   endif
    !end if

  end subroutine couple_gm_pt_init

  !BOP =======================================================================
  !IROUTINE: couple_gm_pt - couple GM to PT
  !INTERFACE:
  subroutine couple_gm_pt(tSimulation)
    
    interface

       subroutine PT_put_from_gm( &
            NameVar, nVar, nPoint, Pos_DI, Data_VI, iPoint_I)

         implicit none
         ! List of variables
         character(len=*), intent(inout):: NameVar 

         ! Number of variables in Data_VI
         integer,          intent(inout):: nVar    

         ! Number of points in Pos_DI
         integer,          intent(inout):: nPoint  

         ! Position vectors
         real, pointer:: Pos_DI(:,:)               

         ! Recv data array
         real,    intent(in), optional:: Data_VI(:,:)

         ! Order of data
         integer, intent(in), optional:: iPoint_I(nPoint)    
       end subroutine PT_put_from_gm
       subroutine GM_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

         implicit none

         ! dimension of position vectors
         integer, intent(in) :: nDimIn             

         ! number of positions
         integer, intent(in) :: nPoint                

         ! positions
         real,    intent(in) :: Xyz_DI(nDimIn,nPoint) 

         ! processor owning position
         integer, intent(out):: iProc_I(nPoint)       

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
    character(len=lNameVar):: NameVar

    ! Grid index
    integer:: iDecompLastGm = -1, iDecompLastPt = -1
    integer:: iGridGm, iDecompGm, iGridPt, iDecompPt

    ! Number of local points for PT and GM cores
    integer:: nPointPt = 0, nPointGm = 0
    
    ! Number of data points found on GM component
    integer:: nData = 0

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

    ! Variables for general-case coupling

    ! number of processors of the OTHER component to communicate with
    integer, save:: nCoupleGm, nCouplePt

    ! processors of the OTHER component to communicate with
    integer, allocatable, save:: iCoupleProcGm_I(:), iCoupleProcPt_I(:)

    ! number of entries received/sent by a processor during rendezvous
    integer, allocatable, save:: nCouplePointGm_I(:), nCouplePointPt_I(:)
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! GM sends all its variables to PT. Take information from Grid_C
    NameVar = Grid_C(GM_)%NameVar
    nVar    = Grid_C(GM_)%nVar

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
    call GM_get_grid_info(nDim, iGridGm, iDecompGm)
    call PT_get_grid_info(nDim, iGridPt, iDecompPt)

    if(.not.IsSameLayout)then
       call MPI_bcast(&
            iDecompGm, 1, MPI_INTEGER, iProc0Gm, iCommGmPt, iError)
       call MPI_bcast(&
            iDecompPt, 1, MPI_INTEGER, iProc0Pt, iCommGmPt, iError)
    endif

    IsNewRoute = iDecompGm /= iDecompLastGm .or. iDecompLastPt /= iDecompPt

    iDecompLastGm = iDecompGm
    iDecompLastPt = iDecompPt

    if(IsNewRoute)then
       nPointPt = 0
       nPointGm = 0
       nData    = 0
       ! Get positions where info is needed from PT. PT will allocate array.
       if(is_proc(PT_))then
          call PT_put_from_gm(NameVar, nVar, nPointPt, PosPt_DI)
       else
          nPointPt = 0
          allocate(PosPt_DI(nDim,0))
       end if

       if(IsSameLayout)then

          ! Find processors that own the PT positions in GM
          allocate(iProcPt_I(nPointPt))

          call GM_find_points(nDim, nPointPt, PosPt_DI, iProcPt_I)

          ! Order points according to the owner processor indexes
          allocate(iPointPt_I(nPointPt))
          call get_buffer_order(nProcGmPt, nPointPt, iProcPt_I, &
               nPointPt_P, iPointPt_I, nData)

          deallocate(iProcPt_I)

          ! Rearrange coordinate array according to processor order
          allocate(PosSortPt_DI(nDim,nData))
          do iPoint = 1, nPointPt
             if(iPointPT_I(iPoint) < 0)CYCLE
             PosSortPt_DI(:,iPointPt_I(iPoint)) = PosPt_DI(:,iPoint)
          end do
          deallocate(PosPt_DI)

          ! Set number of points to be received on the GM component
          call set_recv_info(iCommGmPt, nProcGmPt, iProcGmPt, &
               nPointPt_P, nPointGm, nPointGm_P)

          ! Transfer PT positions to the GM processors that own them
          allocate(PosGm_DI(nDim,nPointGm))
          call transfer_buffer(iCommGmPt, nProcGmPt, iProcGmPt, nDim, &
               nData,    nPointPt_P, PosSortPt_DI, &
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

          ! Setup communication pattern for finding points on GM
          ! This allocates and sets iCoupleProc*_I, nCouplePoint*_I arrays
          call set_inquiry(iCommGmPt, nProcGmPt, iProcGmPt,  &
               GM_, PT_, nProcCommon,                        &
               nCouplePt, iCoupleProcPt_I, nCouplePointPt_I, &
               nCoupleGm, iCoupleProcGm_I, nCouplePointGm_I)

          if(is_proc(PT_))then
             ! Number of points to send from PT to the selected GM processors
             nCouplePointPt_I(1:nCouplePt-1) = nPointPt / nCouplePt 
             ! Last processor gets the rest of points
             nCouplePointPt_I(nCouplePt) = &
                  nPointPt - sum(nCouplePointPt_I(1:nCouplePt-1))
          end if

          ! send number of points from PT to GM that will be sent to 
          ! "GM_find_points." result is in nCouplePointGm_I (number of 
          ! points to be recieved)
          call get_recv_buffer_size(iCommGmPt, nProcGmPt, iProcGmPt, &
               nCouplePt, iCoupleProcPt_I, nCouplePointPt_I,         &
               nCoupleGm, iCoupleProcGm_I, nCouplePointGm_I)

          ! Allocate buffer for positions on GM (zero size on PT)
          if(is_proc(GM_)) nPointGm = sum(nCouplePointGm_I)
          allocate(PosGm_DI(nDim,nPointGm))

          ! Send positions from PT to GM, so GM can find the owners
          call transfer_buffer_real(iCommGmPt, nProcGmPt, iProcGmPt,&
               nCouplePt, iCoupleProcPt_I, nCouplePointPt_I,          &
               nCoupleGm, iCoupleProcGm_I, nCouplePointGm_I,          &
               nDim,                                                  &
               nPointPt, PosPt_DI,&
               nPointGm, PosGm_DI)

          ! Find processors that own the PT positions in GM
          allocate(iProcGm_I(nPointGm))
          if(is_proc(GM_))&
               call GM_find_points(nDim, nPointGm, PosGm_DI, iProcGm_I)
          deallocate(PosGm_DI)

          ! send owner processor indexes from GM to PT
          allocate(iProcPt_I(nPointPt))
          call transfer_buffer_int(iCommGmPt, nProcGmPt, iProcGmPt,&
               nCoupleGm, iCoupleProcGm_I, nCouplePointGm_I,&
               nCouplePt, iCoupleProcPt_I, nCouplePointPt_I,&
               1, &
               nPointGm, iProcGm_I,&
               nPointPt, iProcPt_I)

          ! based on owner information set up the final communication pattern
          call  set_data_transfer(iCommGmPt,                 &
               GM_, PT_,                                     &
               nPointGm, iProcGm_I,                          &
               nPointPt, iProcPt_I,                          &
               nCoupleGm, iCoupleProcGm_I, nCouplePointGm_I, &
               nCouplePt, iCoupleProcPt_I, nCouplePointPt_I)


          if(is_proc(PT_))then
             ! Order points according to the owner processor indexes
             ! sort iProcPt_I according to order of procs in iCoupleProcPt_I
             allocate(iPointPt_I(nPointPt))
             call get_transfer_buffer_order(GM_, PT_,           &
                  nPointPt, iProcPt_I,                          &
                  nCouplePt, iCoupleProcPt_I, nCouplePointPt_I, &
                  iPointPt_I, nData)

             allocate(PosSortPt_DI(nDim, nData))

             ! Rearrange coordinate array according to processor order
             do iPoint = 1, nPointPt
                if(iPointPt_I(iPoint) < 0) CYCLE
                PosSortPt_DI(:,iPointPt_I(iPoint)) = PosPt_DI(:,iPoint)
             end do
          else
             allocate(PosSortPt_DI(nDim, nData))
          end if

          deallocate(iProcPt_I, PosPt_DI, iProcGm_I)

          ! Allocate buffer for positions on GM (zero size on PT) 
          if(is_proc(GM_)) nPointGm = sum(nCouplePointGm_I)
          allocate(PosGm_DI(nDim, nPointGm))

          ! Send PT point positions to the owner GM processors
          call transfer_buffer_real(iCommGmPt, nProcGmPt, iProcGmPt, &
               nCouplePt, iCoupleProcPt_I, nCouplePointPt_I,           &
               nCoupleGm, iCoupleProcGm_I, nCouplePointGm_I,           &
               nDim,                                                   &
               nData,    PosSortPt_DI,              &
               nPointGm, PosGm_DI)

          deallocate(PosSortPt_DI)

       end if
    end if

    ! Get the data from GM
    allocate(DataGm_VI(nVar,nPointGm))
    if(is_proc(GM_)) call GM_get_for_pt( IsNewRoute, &
         NameVar, nVar, nDim, nPointGm, PosGm_DI, DataGm_VI)

    ! Send data from GM to PT into DataPt_VI
    allocate(DataPt_VI(nVar, nData))
    if(IsSameLayout)then
       call transfer_buffer(iCommGmPt, nProcGmPt, iProcGmPt, nVar, &
            nPointGm, nPointGm_P, DataGm_VI, &
            nData,    nPointPt_P, DataPt_VI) 
    else
       call transfer_buffer_real(iCommGmPt, nProcGmPt, iProcGmPt, &
            nCoupleGm, iCoupleProcGm_I, nCouplePointGm_I,           &
            nCouplePt, iCoupleProcPt_I, nCouplePointPt_I,           &
            nVar, nPointGm, DataGm_VI, nData, DataPt_VI)
    end if
    deallocate(DataGm_VI)

    ! Give the data to PT
    if(is_proc(PT_)) call PT_put_from_gm( &
         NameVar, nVar, nPointPt, PosPt_DI, DataPt_VI, iPointPt_I)

    deallocate(DataPt_VI)

    if(DoTest) write(*,*) NameSub,' finished, iProc=',iProcWorld
  end subroutine couple_gm_pt

  !==========================================================================

  subroutine get_buffer_order(&
       nProc, nBuffer, iProc_I, nBuffer_P, iBuffer_I, nData)

    ! nBuffer_P returns the number of points in the buffer belong 
    !    to each processor based on the iProc_I information.
    ! iBuffer_I returns the order of the points with increasing 
    !    processor index.

    integer, intent(in)::  nProc, nBuffer    ! number of procs and points
    integer, intent(in)::  iProc_I(nBuffer)  ! proc index for each point
    integer, intent(out):: nBuffer_P(0:nProc-1)  ! number of points per procs
    integer, intent(out):: iBuffer_I(nBuffer)! index for reordering
    integer, intent(out):: nData ! number of points found on Source

    integer:: iBuffer, iProc
    integer, allocatable:: iBuffer_P(:)
    character(len=*), parameter:: NameSub = 'get_buffer_order'
    !----------------------------------------------------------------------
    nData = nBuffer
    ! Buffer chunk sizes for each processor
    nBuffer_P = 0
    do iBuffer = 1, nBuffer
       iProc = iProc_I(iBuffer)
       if(iProc < 0) CYCLE
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
       if(iProc < 0)then
          ! the point hasn't been found
          iBuffer_I(iBuffer) = - 1
          nData = nData - 1
       else
          ! It will occupy iBufferSort position
          iBuffer_I(iBuffer) = iBuffer_P(iProc)
          ! Jump to next position
          iBuffer_P(iProc) = iBuffer_P(iProc) + 1
       end if
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

    allocate(iRequestS_I(nProc), iRequestR_I(nProc), &
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

  subroutine set_inquiry(iComm, nProc, iProc,                    &
       iCompSource, iCompTarget, nProcCommon,                    &
       nCoupleTarget, iCoupleProcTarget_I, nCouplePointTarget_I, &
       nCoupleSource, iCoupleProcSource_I, nCouplePointSource_I)

    ! the subroutine returns number of processors to be communicated with
    !     nCoupleSource, nCoupleTarget 
    ! and stores indexes of these processors in 
    !     iCoupleProcSource_I, iCoupleProcTarget_I
    ! iCompSource is the component which Data is stored on
    ! iCompTarget is the component which Data is to be sent to
    ! in the Union communicator procs on Source are assumed to come first

    integer,              intent(in ):: iComm, nProc, iProc
    integer,              intent(in ):: iCompSource, iCompTarget
    integer,              intent(in ):: nProcCommon

    integer,              intent(out):: nCoupleTarget
    integer, allocatable, intent(out):: iCoupleProcTarget_I(:)
    integer, allocatable, intent(out):: nCouplePointTarget_I(:)

    integer,              intent(out):: nCoupleSource
    integer, allocatable, intent(out):: iCoupleProcSource_I(:)
    integer, allocatable, intent(out):: nCouplePointSource_I(:)

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
       if(is_proc(iCompTarget))then
          ! Target and source procs coincide, communicates with itself
          nCoupleTarget = 1
          nCoupleSource = 1
       elseif(is_proc(iCompSource))then
          ! A source proc that is not involved
          nCoupleTarget = 0
          nCoupleSource = 0
       end if
       allocate(iCoupleProcTarget_I(nCoupleTarget))
       allocate(iCoupleProcSource_I(nCoupleSource))

       if(nCoupleSource*nCoupleTarget > 0)then
          iCoupleProcTarget_I(nCoupleTarget) = iProc
          iCoupleProcSource_I(nCoupleSource) = iProc
       end if

       allocate(nCouplePointTarget_I(nCoupleTarget))
       allocate(nCouplePointSource_I(nCoupleSource))

       RETURN
    end if


    ! if nProcSource/nProcTarget is larger and not an exact multiple of the 
    ! other, there is a discrepancy of dProc processors
    if(nProcSource > nProcTarget)then
       dProc = MOD(nProcSource, nProcTarget)
    else
       dProc = MOD(nProcTarget, nProcSource)
    end if

    ! in the communicator iComm processor are ordered as follows
    ! [ SOURCE | TARGET ]
    ! find the local index of iProc on its component,
    ! index of 0-processor of the other component,
    ! # of procs to couple with: nCouple, 
    ! # of procs other component is coupled with: nCoupleOther
    if(is_proc(iCompSource))then
       nCouple      = nProcTarget / nProcSource
       nCoupleOther = nProcSource / nProcTarget
       iProcLocal   = iProc
       nProc0Other  = nProcSource
    elseif(is_proc(iCompTarget))then
       nCouple      = nProcSource / nProcTarget
       nCoupleOther = nProcTarget / nProcSource
       iProcLocal   = iProc - nProcSource
       nProc0Other  = 0
    end if

    ! account for discrepancy
    if(iProcLocal < dProc * (nCoupleOther+1))then
       dProc = 0
       nCoupleOther = nCoupleOther + 1
       nCouple = nCouple + 1
    elseif(nCouple == 0)then
       dProc =  - dProc
    end if
    
    nCouple      = MAX(1, nCouple)
    nCoupleOther = MAX(1, nCoupleOther)

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
       nBufferSource, iProcSource_I, &
       nBufferTarget, iProcTarget_I, &
       nCoupleSource, iCoupleProcSource_I, nCouplePointSource_I,&
       nCoupleTarget, iCoupleProcTarget_I, nCouplePointTarget_I)

    ! 

    ! subroutine defines transfer based on processors 
    ! of the OTHER components, which own data items
    ! iProc_I

    ! Union communicator
    integer, intent(in):: iComm                        

    ! Source and target component indexes
    integer, intent(in):: iCompSource, iCompTarget

    ! Number of point positions asked from this source proc
    integer, intent(in):: nBufferSource
    ! Index of source processors owning these points
    integer, intent(in):: iProcSource_I(nBufferSource)

    ! Number of points required by the target processor
    integer, intent(in):: nBufferTarget
    ! Index of source processor owning the points of this target proc.
    integer, intent(in):: iProcTarget_I(nBufferTarget)

    ! Number of couplings performed by source on this processor
    integer, intent(out)               :: nCoupleSource
    ! Processor indexes in iComm that source is coupled with
    integer, intent(inout), allocatable:: iCoupleProcSource_I( :)
    ! Number of points per processor source is coupled with
    integer, intent(inout), allocatable:: nCouplePointSource_I(:)

    ! Number of couplings performed by target on this processor
    integer, intent(out)               :: nCoupleTarget
    ! Processor indexes in iComm that target is coupled with
    integer, intent(inout), allocatable:: iCoupleProcTarget_I( :)
    ! Number of points per processor target is coupled with
    integer, intent(inout), allocatable:: nCouplePointTarget_I(:)

    integer:: nProcSource, nProcTarget, nProcIntersect
    integer:: iBuffer, iProc, iCouple, iCoupleBuffer 
    integer:: iProcSourceLocal, iProcTargetLocal
    integer:: iError
    integer, allocatable:: iProcSourceLocal_P(:), iProcTargetLocal_P(:)
    integer, allocatable:: iProcSourceUnion_P(:), iProcTargetUnion_P(:)
    integer, allocatable:: nPoint_PP(:,:), nPointRecv_P(:), nPointSend_P(:) 
    integer, allocatable:: iCoupleProc_P(:), nCouplePoint_P(:)
    character(len=*), parameter:: NameSub = 'set_data_transfer'
    integer:: iTmp
    !----------------------------------------------------------------------
    nProcTarget    = n_proc(iCompTarget)
    nProcSource    = n_proc(iCompSource)

    if(is_proc(iCompTarget))then
       if(allocated(iCoupleProcTarget_I )) deallocate(iCoupleProcTarget_I)
       if(allocated(nCouplePointTarget_I)) deallocate(nCouplePointTarget_I)
       !\
       ! get nCouple and number of points per proc
       !/

       ! Allocate temporary processor index and buffer size arrays.
       ! nProcSource is the maximum possible number of processors to 
       ! communicate with. The actual size will be reduced later.
       allocate(iCoupleProc_P(0:nProcSource-1), &
            nCouplePoint_P(0:nProcSource-1))
       nCouplePoint_P = 0

       ! count points received from each source processors
       do iBuffer = 1, nBufferTarget
          iProc = iProcTarget_I(iBuffer)
          if(iProc < 0) CYCLE
          nCouplePoint_P(iProc) = nCouplePoint_P(iProc) + 1
       end do

       ! Skip source processors that have no data to send.
       ! Count the source processors that have data (nCoupleTarget).
       ! Compact the arrays to the source processors that have data.
       nCoupleTarget = 0
       do iProc = 0, nProcSource-1
          if(nCouplePoint_P(iProc) == 0) CYCLE
          ! One more source proc with data
          nCoupleTarget = nCoupleTarget + 1
          ! Move up the proc index and number of points 
          ! from iProc to nCoupleTarget-1
          iCoupleProc_P( nCoupleTarget-1) = iProc
          nCouplePoint_P(nCoupleTarget-1) = nCouplePoint_P(iProc)
       end do

       ! Copy compacted arrays into the output arguments
       allocate(&
            iCoupleProcTarget_I( nCoupleTarget), &
            nCouplePointTarget_I(nCoupleTarget))

       iCoupleProcTarget_I  = iCoupleProc_P(0:nCoupleTarget-1)
       nCouplePointTarget_I = nCouplePoint_P(0:nCoupleTarget-1)
       deallocate(iCoupleProc_P, nCouplePoint_P)
    end if

    if(is_proc(iCompSource))then
       allocate(iProcSourceUnion_P(0:nProcSource-1))
       allocate(iProcTargetUnion_P(0:nProcTarget-1))
       allocate(iProcSourceLocal_P(0:nProcSource+nProcTarget-1))
       allocate(iProcTargetLocal_P(0:nProcSource+nProcTarget-1))
       call get_rank_translation(&
            iComm, iCompSource, iCompTarget, nProcSource, nProcTarget, &
            iProcSourceLocal_P, iProcTargetLocal_P,                    &
            iProcSourceUnion_P, iProcTargetUnion_P)

       ! nPoint_PP is number of points to sent from Source to Target
       allocate(nPoint_PP(0:nProcSource-1, 0:nProcTarget-1))

       ! build nPoint_PP
       ! each proc on Source has only part of the info, store it into nPoint_PP
       nPoint_PP = 0
       if(nCoupleSource > 0)then
          iCouple = 1
          iProcTargetLocal = iProcTargetLocal_P(iCoupleProcSource_I(iCouple))
          iCoupleBuffer = nCouplePointSource_I(iCouple)
          do iBuffer = 1, nBufferSource
             iProc = iProcSource_I(iBuffer)
             if(iProc < 0) CYCLE
             if(iBuffer > iCoupleBuffer)then
                iCouple = iCouple + 1
                iCoupleBuffer = iCoupleBuffer + nCouplePointSource_I(iCouple)
                ! iCoupleProcSource_I(iCouple) is the index of target proc for
                ! union communicator that sent the request to this source 
                ! processor. iProcTargetLocal is the index of taget processor 
                ! on the target component
                iProcTargetLocal = iProcTargetLocal_P(&
                     iCoupleProcSource_I(iCouple))
             end if
             ! iProcSource_I is the source processor index on source component
             iProcSourceLocal = iProcSourceLocal_P(iProc)
             
             ! Count this point in the communication matrix
             nPoint_PP(iProcSourceLocal,iProcTargetLocal) = &
                  nPoint_PP(iProcSourceLocal,iProcTargetLocal) + 1
          end do
       end if

       deallocate(iCoupleProcSource_I, nCouplePointSource_I)

       !\
       ! Reduce the row nPoint_PP(iProcSourceLocal,:)
       !/
       allocate(nPointSend_P(0:nProcTarget-1), nPointRecv_P(0:nProcTarget-1))

       ! Each source processor collects the row of the matrix that it needs.
       ! Each row is added up from all source processors, because they
       ! were asked about point positions randomly.
       do iProcSourceLocal = 0, nProcSource-1
          ! Extract the particular row (partial info)
          nPointSend_P = nPoint_PP(iProcSourceLocal,:)
          ! Add up partial rows and collect it on iProcSourceLocal
          call MPI_Reduce(&
               nPointSend_P, nPointRecv_P, nProcTarget, MPI_INTEGER, MPI_SUM, &
               iProcSourceLocal, i_comm(iCompSource), iError)
       end do

       deallocate(nPoint_PP, nPointSend_P)

       ! Skip target processors that have no data to receive.
       ! Count the target processors that have data (nCoupleSource).
       ! Compact the arrays to the target processors that have data.
       ! 
       allocate(iCoupleProc_P(0:nProcTarget-1), &
            nCouplePoint_P(0:nProcTarget-1))
       nCoupleSource = 0
       do iProcTargetLocal = 0, nProcTarget-1
          if(nPointRecv_P(iProcTargetLocal) == 0) CYCLE
          nCoupleSource = nCoupleSource + 1

          ! Translate local rank iProcTargetLocal to union rank
          iCoupleProc_P(nCOupleSource-1) = iProcTargetUnion_P(iProcTargetLocal)
          nCouplePoint_P(nCoupleSource-1) = nPointRecv_P(iProcTargetLocal)
       end do

       ! Copy compacted arrays into the output arguments
       allocate(iCoupleProcSource_I( nCoupleSource))
       allocate(nCouplePointSource_I(nCoupleSource))
       iCoupleProcSource_I  = iCoupleProc_P(0:nCoupleSource-1)
       nCouplePointSource_I = nCouplePoint_P(0:nCoupleSource-1)

       deallocate(iCoupleProc_P, nCouplePoint_P, nPointRecv_P)
       deallocate(iProcSourceLocal_P, iProcSourceUnion_P)
       deallocate(iProcTargetLocal_P, iProcTargetUnion_P)
    end if

  end subroutine set_data_transfer

  !==========================================================================

  subroutine get_transfer_buffer_order(&
       iCompSource, iCompTarget, &
       nBuffer, iProc_I, &
       nCouple, iCoupleProc_I, nCouplePoint_I,&
       iBuffer_I, nData)
    ! iBuffer_I returns the order of the points with processor index
    ! corresponding to iCoupleProc_I
    ! nData is the number of points found on Source
    integer, intent(in) :: iCompSource, iCompTarget ! send, recv components
    integer, intent(in) :: nBuffer    ! number of points
    integer, intent(in) :: iProc_I(nBuffer)  ! proc index for each point
    integer, intent(in) :: nCouple
    integer, intent(in) :: iCoupleProc_I(nCouple), nCouplePoint_I(nCouple)
    integer, intent(out):: iBuffer_I(nBuffer)! index for reordering
    integer, intent(out):: nData ! number of points found on Source

    integer:: iBuffer, iCouple, iProc, nProc, nPointNotFound
    integer, allocatable:: iOrder_I(:)
    integer, allocatable:: iCoupleProcInv_P(:)
    character(len=*), parameter:: NameSub = 'get_transfer_buffer_order'
    !----------------------------------------------------------------------
    nProc = n_proc(iCompSource)
    nData = nBuffer
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
       if(iProc < 0)then
          ! the point hasn't been found
          iBuffer_I(iBuffer) = - 1
          nData = nData - 1
       else
          ! It will occupy iBuffer_I position
          iBuffer_I(iBuffer) = iOrder_I(iCoupleProcInv_P(iProc))
          ! Jump to next position
          iOrder_I(iCoupleProcInv_P(iProc)) = &
               iOrder_I(iCoupleProcInv_P(iProc)) + 1
       end if
    end do

    deallocate(iCoupleProcInv_P, iOrder_I)

  end subroutine get_transfer_buffer_order

  !===========================================================================

  subroutine transfer_buffer_real(&
       iComm, nProc, iProc,&
       nCoupleS, iCoupleProcS_I, nCouplePointS_I,&
       nCoupleR, iCoupleProcR_I, nCouplePointR_I,&
       nData, &
       nBufferS, BufferS_I, nBufferR, BufferR_I)

    ! This subroutine transfers array of reals
    ! A processor can both send and recv
    ! Correctness: only buffers of the same type should be present 

    use ModMpi


    integer, intent(in):: iComm  ! MPI communicator
    integer, intent(in):: nProc  ! number of processors
    integer, intent(in):: iProc  ! local proc index

    integer, intent(in):: nCoupleS                 ! #   of procs to send to
    integer, intent(in):: iCoupleProcS_I( nCoupleS)! ids of procs to send to
    integer, intent(in):: nCouplePointS_I(nCoupleS)! # of points  to send

    integer, intent(in):: nCoupleR                 ! #   of procs to recv from
    integer, intent(in):: iCoupleProcR_I( nCoupleR)! ids of procs to recv from
    integer, intent(in):: nCouplePointR_I(nCoupleR)! # of points  to recv

    integer, intent(in)   :: nData                 ! number of items per point
    integer, intent(in)   :: nBufferS                  ! send buffer size
    real,    intent(inout):: BufferS_I(nData*nBufferS) ! send buffer
    integer, intent(in)   :: nBufferR                  ! recv buffer size
    real,    intent(inout):: BufferR_I(nData*nBufferR) ! recv buffer

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
    ! Send buffer
    !/
    if(nCoupleS > 0)then
       !\
       ! Send array of reals
       !/
       iBuffer = 1
       do iCouple = 1, nCoupleS
          iProcS = iCoupleProcS_I(iCouple)
          call MPI_isend(&
               BufferS_I(iBuffer), nData*nCouplePointS_I(iCouple), &
               MPI_REAL, iProcS, iTag, iComm, iRequest_I(iCouple),iError)
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
          call MPI_irecv(&
               BufferR_I(iBuffer), nData*nCouplePointR_I(iCouple),&
               MPI_REAL, iProcR, iTag, iComm, &
               iRequest_I(nCoupleS+iCouple),iError)
          iBuffer = iBuffer + nData*nCouplePointR_I(iCouple)
       end do
    end if
    !\
    ! Finalize transfer
    !/
    call MPI_waitall(nCoupleS + nCoupleR, iRequest_I, iStatus_II, iError)
    deallocate(iRequest_I, iStatus_II)
  end subroutine transfer_buffer_real

  !===========================================================================

  subroutine transfer_buffer_int(&
       iComm, nProc, iProc,&
       nCoupleS, iCoupleProcS_I, nCouplePointS_I,&
       nCoupleR, iCoupleProcR_I, nCouplePointR_I,&
       nData, &
       nBufferS, iBufferS_I, nBufferR, iBufferR_I)

    ! This subroutine transfers array of integers
    ! A processor can both send and recv
    ! Correctness: only buffers of the same type should be present 

    use ModMpi

    integer, intent(in):: iComm  ! MPI communicator
    integer, intent(in):: nProc  ! number of processors
    integer, intent(in):: iProc  ! local proc index

    integer, intent(in):: nCoupleS                 ! #   of procs to send to
    integer, intent(in):: iCoupleProcS_I( nCoupleS)! ids of procs to send to
    integer, intent(in):: nCouplePointS_I(nCoupleS)! # of points  to send

    integer, intent(in):: nCoupleR                 ! #   of procs to recv from
    integer, intent(in):: iCoupleProcR_I( nCoupleR)! ids of procs to recv from
    integer, intent(in):: nCouplePointR_I(nCoupleR)! # of points  to recv

    integer, intent(in   ):: nData                 ! number of items per point
    integer, intent(in   ):: nBufferS                   ! send buffer size
    integer, intent(inout):: iBufferS_I(nData*nBufferS) ! send buffer
    integer, intent(in)   :: nBufferR                   ! recv buffer size
    integer, intent(inout):: iBufferR_I(nData*nBufferR) ! recv buffer

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
    ! Send buffer
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

    !\
    ! Finalize transfer
    !/
    call MPI_waitall(nCoupleS + nCoupleR, iRequest_I, iStatus_II, iError)
    deallocate(iRequest_I, iStatus_II)
  end subroutine transfer_buffer_int

  !===========================================================================

  subroutine get_rank_translation(                                      &
       iCommUnion, iCompSource, iCompTarget, nProcSource, nProcTarget,  &
       iProcSourceLocal_I, iProcTargetLocal_I,                          &
       iProcSourceUnion_I, iProcTargetUnion_I)
    ! for processor with rank iProc in the union communicator iCommUnion
    ! returns rank iProcLocal in the group associated with component iComp
    integer, intent(in) :: iCommUnion, iCompSource, iCompTarget
    integer, intent(in) :: nProcSource, nProcTarget
    integer, intent(out):: iProcSourceLocal_I(nProcSource+nProcTarget)
    integer, intent(out):: iProcTargetLocal_I(nProcSource+nProcTarget)
    integer, intent(out):: iProcSourceUnion_I(nProcSource)
    integer, intent(out):: iProcTargetUnion_I(nProcTarget)

    integer:: iGroupUnion
    integer:: iProc
    integer, allocatable:: iProc_I(:)
    integer:: iError
    !----------------------------------------------------------------------
    ! Create group associated with iCommUnion
    call MPI_Comm_Group(iCommUnion, iGroupUnion, iError)

    ! Translate ranks for Source Local to Union
    allocate(iProc_I(nProcSource))
    do iProc = 0, nProcSource-1
       iProc_I(iProc+1) = iProc
    end do
    call MPI_Group_translate_ranks(&
         i_group(iCompSource), nProcSource, iProc_I, &
         iGroupUnion, iProcSourceUnion_I, iError)
    deallocate(iProc_I)

    ! Translate ranks for Target Local to Union
    allocate(iProc_I(nProcTarget))
    do iProc = 0, nProcTarget-1
       iProc_I(iProc+1) = iProc
    end do
    call MPI_Group_translate_ranks(&
         i_group(iCompTarget), nProcTarget, iProc_I, &
         iGroupUnion, iProcTargetUnion_I, iError)
    deallocate(iProc_I)

    ! Translate ranks Union to Local
    allocate(iProc_I(nProcSource+nProcTarget))
    do iProc = 0, nProcSource+nProcTarget-1
       iProc_I(iProc+1) = iProc
    end do
    call MPI_Group_translate_ranks(&
         iGroupUnion, nProcSource+nProcTarget, iProc_I, &
         i_group(iCompSource), iProcSourceLocal_I, iError)
    call MPI_Group_translate_ranks(&
         iGroupUnion, nProcSource+nProcTarget, iProc_I, &
         i_group(iCompTarget), iProcTargetLocal_I, iError)
    deallocate(iProc_I)

    ! free union group
    call MPI_Group_free(iGroupUnion, iError)
  end subroutine get_rank_translation

end module CON_couple_gm_pt
