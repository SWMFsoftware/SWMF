!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_global_message_pass

  use CON_router
  use ModUtilities, ONLY: CON_stop

  implicit none

  ! The toolkit for coupling the different data sets within a
  ! single component code, or within a framework, as it is now (see
  ! the date above) consists of:
  !
  ! ModMpi (iclude for mpif90.h)
  !
  ! ModNumConst
  !
  ! CON\_registry
  !
  ! CON\_domain\_decomposition
  !
  ! CON\_grid\_descriptor
  !
  ! CON\_router
  !
  ! CON\_global\_message\_pass
  !
  ! Sokolov I.V.
  ! 7.20.03-7.21.03
  ! igorsok@umich.edu
  ! phone(734)647-4705
  !
  ! The final file of the toolkit. Using the router, prepared by
  ! CON\_router, the messages are sent and received by different PEs
  ! within a single component and/or between the components
  !
  !
  ! prototypes:
  !
  ! BATSRUS/src/message\_pass\_cells.f90 (D. DeZeeuw,2002)
  ! BATSRUS/src/message\_pass\_nodes.f90 (D. DeZeeuw,2002-2003)
  ! BATSRUS/src/message\_pass\_faces.f90 (D. DeZeeuw,2002)
  !
  ! New features as compared to the prototypes:
  ! The use of a router of a standard format.
  ! Reusable single-index buffers.
  ! Generic interface.
  !
  ! At any call of the CON\_global\_message\_pass the same buffer
  ! arrays are used
  !\begin{verbatim}

  real, allocatable, save,private :: VSend_I(:), VRecv_I(:)
  !\end{verbatim}
  ! These parameters can be readjusted for a better performance
  ! with a given computer.

  ! This parameter does not allow to use too large amount of non-
  ! blocking send-receives
  integer, parameter,private :: maxMessages=10000

  ! Ususally there is not more than one message sent within any
  ! pair of PEs. If this appears to be enefficient or is not
  ! stable, the messages can be broken up...
  logical, parameter,private :: DoBreakUpMessages=.true.

  !... and in this case the following parameter defines the size
  ! for the message pieces.
  integer, parameter,private  :: MessageSize=99999

  ! NOTE: with DoRSend=.false. there is no barrier at all in the
  ! global\_message\_pass procedure.

contains
  !============================================================================
  subroutine global_message_pass(Router, nVar, fill_buffer, apply_buffer)

    ! Router: should be preset. See CON_router for the instructions.
    ! nVar is the number of REAL variables which are sent-received.
    ! nVar>1 means that more than one scalar field is involved. In
    ! 3D Magnetohydrodynamics nVar may amount to 8 (density, three
    ! components of the velocity/momentum, three components of the
    ! magnetic field, temperature or pressure).
    ! You should use two (optionally three) external subroutines. It
    ! is important to mention that in the procedure itself there is
    ! absolutely no names of the data sets to be sent-received. You
    ! can send whatever you want but you should form the procedure
    ! which GET the data you want from the points defined by router
    ! (the cell and block indexes for these points are submitted by
    ! global_message_pass itself, and you should place
    ! point-by-point, whatever you  want to a short (of the nVar
    ! length) piece of the buffer which then will be placed to a
    ! particular place of the global buffer and proberly
    ! sent-received. A name of your  procedure which gets your data
    ! and place them to a small buffer should be used here. Use
    ! fill_buffer=your_procedure_name best of all...
    ! ... and then you should apply the data from the RECV buffer.
    ! Just analogously you will obtain nVar real variables with the
    ! "address" in the form of cell and block index, at the proper
    ! PE and should process them as desired and (usually) put to some
    ! global arrays with the obtained values for the cell and block
    ! indexes (use apply_buffer=yet_another_procedure)
    ! The last procedure is optional and can be only used for
    ! somewhat better efficiency. The global_message_pass applies
    ! this procedure only if the data to be sent and the points they
    ! are to be received appear to be located at the same PE. Of
    ! course the data can be put for a while to some temporary buffer
    ! using fill_buffer and then invoke apply_buffer for the content
    ! of this temporary buffer (this is what is done if the optional
    ! procedure "copy" is not present), but we leave the opportunity
    ! for a user to organise this operation more efficiently, if
    ! desired. In the exchange between the components with
    ! non-clashing PE sets, this procedure is not applied at all and
    ! it is pointless to have it present.

    type(RouterType),intent(in) :: Router
    integer,         intent(in) :: nVar
    interface
       ! Should get nVar values from some arrays using indexes
       !
       ! Get%Index_II(:,iGetStart:iGetStart+nPartial-1), take the sum(s)
       ! with weight coefficients.
       ! Weight_I(iGetStart:iGetStart+nPartial-1), and put the nVar
       ! results to the buffer array VSend_I, starting from iBuffStart
       ! value of indexes.
       !
       !
       ! and the vector V can be the set of the conserved variables in
       ! the thus numbered control volumes

       subroutine fill_buffer(nPartial, iGetStart, Get, Weight, Buff_I, nVar)
         use CON_router
         implicit none
         integer,intent(in)::nPartial,iGetStart,nVar
         type(IndexPtrType),intent(in)::Get
         type(WeightPtrType),intent(in)::Weight
         real,dimension(nVar),intent(out)::Buff_I
       end subroutine fill_buffer
       ! Inverse operation: gets nVar variables from
       ! iRecv_I(iBuff:iBuff+nVar-1) and puts it to the state vector
       ! with the index vector being equal to Put%iCB_II(:,iPut)
       subroutine apply_buffer(nPartial, &
            iPutStart, Put, Weight, DoAdd, Buff_I, nVar)
         use CON_router
         implicit none
         integer,intent(in)::nPartial,iPutStart,nVar
         type(IndexPtrType),intent(in)::Put
         type(WeightPtrType),intent(in)::Weight
         logical,intent(in)::DoAdd
         real,dimension(nVar),intent(in)::Buff_I
       end subroutine apply_buffer
    end interface

    integer :: iProc,nProc,iComm,iCommUnion
    integer :: iV,iPE, iError, nError, lBuff, lGet, lPut
    integer :: nPartialGet, nPartialPut
    integer :: iTag, nSends, i
    integer :: nSendAll, nRecvAll
    integer :: nSENDrequests, SENDrequests(MaxMessages)
    integer :: nRECVrequests, RECVrequests(MaxMessages)
    integer :: MESGstatus(MPI_STATUS_SIZE, MaxMessages)
    integer,dimension(0:Router%nProc-1) :: iSendStart_P,iRecvStart_P
    real    :: TempBuff_I(nVar)
    ! The communication set is stored in the Router
    !--------------------------------------------------------------------------
    if(.not.Router%IsProc)RETURN
    iProc=Router%iProc
    nProc=Router%nProc
    iComm=Router%iComm
    iCommUnion=Router%iCommUnion

    iSendStart_P=0
    iRecvStart_P=0
    do iPE=0,nProc-2
       if(iPE==iProc)then
          iSendStart_P(iPE+1) = iSendStart_P(iPE)
          iRecvStart_P(iPE+1) = iRecvStart_P(iPE)
       else
          iSendStart_P(iPE+1) = iSendStart_P(iPE) + Router%nSend_P(iPE)
          iRecvStart_P(iPE+1) = iRecvStart_P(iPE) + Router%nRecv_P(iPE)
       end if
    end do
    if(nVar < 1)call CON_stop('In GlobalMessagePass you can not use nVar<1')
    ! Check send/recv buffers
    nSendAll = (iSendStart_P(nProc-1) + Router%nSend_P(nProc-1))*nVar

    if(nSendAll>0)call check_send_buffer(nSendAll)

    nRecvAll=(iRecvStart_P(nProc-1) + Router%nRecv_P(nProc-1))*nVar
    if(nRecvAll>0)call check_recv_buffer(nRecvAll)

    iPE = iProc
    lGet=0; lPut=0
    if(Router%nSend_P(iPE)>0)then
       do iV=1,Router%nSend_P(iPE)
          nPartialGet = Router%iGet_P(iPE)%iCB_II(0,lGet+1)
          nPartialPut = Router%iPut_P(iPE)%iCB_II(0,lPut+1)
          call fill_buffer(nPartialGet, lGet+1, Router%iGet_P(iPE),   &
               Router%Get_P(iPE), TempBuff_I,nVar)
          call apply_buffer(nPartialPut, lPut+1, Router%iPut_P(iPE),  &
               Router%Put_P(iPE), Router%DoAdd_P(iPE)%DoAdd_I(lPut+1),&
               TempBuff_I,nVar)
          lGet = lGet + nPartialGet
          lPut = lPut + nPartialPut
       end do
    end if

    ! Collect values into VSend that need to be passed to other processors
    lBuff=0
    do iPE = 0, nProc - 1
       if(iPE==iProc) CYCLE
       if(Router%nSend_P(iPE)==0) CYCLE
       lGet=0
       do iV=1,Router%nSend_P(iPE)
          nPartialGet=Router%iGet_P(iPE)%iCB_II(0,lGet+1)
          call fill_buffer(nPartialGet, lGet+1, Router%iGet_P(iPE),&
               Router%Get_P(iPE), VSend_I(lBuff+1:lBuff+nVar), nVar)
          lGet=lGet+nPartialGet
          lBuff=lBuff+nVar
       end do
    end do

    ! Post receives first so that they are ready
    nRECVrequests = 0
    do iPE = 0, nProc - 1
       if(iPE == iProc) CYCLE
       if(Router%nRecv_P(iPE)==0) CYCLE
       if(DoBreakUpMessages)then
          nSends=1+((Router%nRecv_P(iPE)-1)/MessageSize)
          do i=1,nSends
             itag = nProc*(i-1)+iPE
             nRECVrequests = nRECVrequests + 1
             if(nRECVrequests > maxMessages)then
                write(*,*)'Too many messages!!'
                call MPI_Abort(iComm,iError,nError)
             end if
             call MPI_irecv(VRecv_I(1 + nVar*                   &
                  (iRecvStart_P(iPE) + (i - 1)*MessageSize)),   &
                  nVar*min(MessageSize,                         &
                  Router%nRecv_P(iPE) - ((i - 1)*MessageSize)), &
                  MPI_REAL, iPE, itag, iComm,                   &
                  RECVrequests(nRECVrequests),iError)
          end do
       else
          itag = iPE
          nRECVrequests = nRECVrequests + 1
          if(nRECVrequests > maxMessages)&
               call CON_stop("Too many RECVs in mp_SendValues")
          call MPI_irecv(VRecv_I(iRecvStart_P(iPE)*nVar+1),         &
               nVar*Router%nRecv_P(iPE), MPI_REAL, iPE, itag, iComm,&
               RECVrequests(nRECVrequests), iError)
       end if
    end do

    if(DoRSend)then
       ! Make sure all recv's are posted before using an rsend
       call MPI_BARRIER(iCommUnion,iError)
    end if

    ! VSend array sent to VRecv array on other PEs
    nSENDrequests = 0
    do iPE = 0, nProc - 1
       if(iPE == iProc) CYCLE
       if(Router%nSend_P(iPE)==0) CYCLE
       if(DoBreakUpMessages)then
          nSends=1 + ((Router%nSend_P(iPE)-1)/MessageSize)
          do i=1,nSends
             itag = nProc*(i - 1) + iProc
             if(DoRSend)then
                call MPI_rsend(VSend_I(1+nVar*&
                     (iSendStart_P(iPE) + (i-1)*MessageSize)),&
                     nVar*min(MessageSize,&
                     Router%nSend_P(iPE) - ((i - 1)*MessageSize)), &
                     MPI_REAL, iPE, itag, iComm, iError)
             else
                nSENDrequests = nSENDrequests + 1
                call MPI_isend(VSend_I(1+nVar*&
                     (iSendStart_P(iPE) + (i-1)*MessageSize)),&
                     nVar*min(MessageSize,&
                     Router%nSend_P(iPE) - ((i - 1)*MessageSize)), &
                     MPI_REAL, iPE, itag, iComm, &
                     SENDrequests(nSENDrequests), iError)
             end if
          end do
       else
          itag = iProc
          if(DoRSend)then
             call MPI_rsend(VSend_I(1 + nVar*iSendStart_P(iPE)),&
                  nVar*Router%nSend_P(iPE), MPI_REAL, iPE, itag, iComm, iError)
          else
             nSENDrequests = nSENDrequests + 1
             call MPI_isend(VSend_I(1 + nVar*iSendStart_P(iPE)),        &
                  nVar*Router%nSend_P(iPE), MPI_REAL, iPE, itag, iComm, &
                  SENDrequests(nSENDrequests), iError)
          end if
       end if
    end do

    ! Wait for messages to be received before continuing.
    call MPI_waitall(nRECVrequests, RECVrequests(1), MESGstatus(1,1), iError)
    lBuff = 0
    do iPE = 0, nProc - 1
       if(iPE==iProc) CYCLE
       if(Router%nRecv_P(iPE)==0) CYCLE
       lPut=0
       do iV = 1,Router%nRecv_P(iPE)
          nPartialPut=Router%iPut_P(iPE)%iCB_II(0,lPut+1)
          call apply_buffer(nPartialPut, lPut+1, Router%iPut_P(iPE),  &
               Router%Put_P(iPE), Router%DoAdd_P(iPE)%DoAdd_I(lPut+1),&
               VRecv_I(lBuff+1:lBuff+nVar),nVar)
          lPut  = lPut  + nPartialPut
          lBuff = lBuff + nVar
       end do
    end do
  contains
    !==========================================================================
    subroutine check_recv_buffer(nSize)
      integer,intent(in)::nSize
      integer::iError,nError
      !------------------------------------------------------------------------
      if(allocated(VRecv_I))then
         if(ubound(VRecv_I,1) < nSize)then
            deallocate(VRecv_I)
         else
            RETURN
         end if
      end if
      allocate(VRecv_I(nSize), stat = iError)
      call check_allocate(iError, 'VRecv')
    end subroutine check_recv_buffer
    !==========================================================================
    subroutine check_send_buffer(nSize)
      integer,intent(in)::nSize
      integer::iError,nError
      !------------------------------------------------------------------------
      if(allocated(VSend_I))then
         if(ubound(VSend_I,1)<nSize)then
            deallocate(VSend_I)
         else
            RETURN
         end if
      end if
      allocate(VSend_I(nSize), stat = iError)
      call check_allocate(iError, 'VSend_I')
    end subroutine check_send_buffer
    !==========================================================================
  end subroutine global_message_pass
  !============================================================================
end module CON_global_message_pass
!==============================================================================

