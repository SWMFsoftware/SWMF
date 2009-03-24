!^CFG COPYRIGHT UM
module ModMPCells

  integer :: numSendRecv, numSend, numRecv, numCopy
  integer :: numSendRecvMax=0, numSendMax=0, numRecvMax=0, &
       numCopyMax=0,  numSendMax8=0, numRecvMax8=0

  integer, dimension(:),   allocatable, save :: &
       nSend, nRecv, nSendStart, nRecvStart
  integer, dimension(:,:), allocatable, save :: &
       VSendI, VRecvI, VSendIlocal, VRecvIlocal

  real,    dimension(:),   allocatable, save :: VSend, VRecv
  real,    dimension(:,:), allocatable, save :: VSend8, VRecv8

  integer :: iLastGrid = -1, iLastDecomposition = -1
  integer :: itag, lS(0:7), lR(0:7), nSends
  integer, parameter :: maxMessages=10000
  integer :: MessageSize=99999

  ! Coarse to Fine Exchange options
  ! iCFExchangeType=1  Old style: send 8 identical values to 8 fine cells
  ! iCFExchangeType=2  New style: send 1 value to 8 fine cells
  ! iCFExchangeType=3  Prolongation: send 8 different values to 8 fine cells
  !                    using simple prolongation
  integer, parameter :: iCFExchangeType=2
  real, parameter    :: Inv12 = 1.0/12.0

  !With DoOneCoarserLayer=.false., the finer cells of the "second" layer are filled in 
  !with the values from the "second" layer of the coarser block
  logical::DoOneCoarserLayer=.true.
  Logical::DoOneLayerCorner=.false.
  
  !With DoOneCoarserLayer=.true., the following nDuplicateIJK=2 and
  !the values from only one layer of the physical cells are sent to the finer neighbor
  
  integer:: nDuplicateI=2,nDuplicateJ=2,nDuplicateK=2
  integer:: iTwoOrOneForTwoCoarserLayers=2,iZeroOrOneForTwoCoarserLayers=0


  logical :: DoFacesOnlyLast, DoOneLayerLast, DoOneLayer_D(26)=.false.
  logical, parameter :: DoLimitCornerMemory=.false.
  logical, parameter :: DoRSend=.true.
  logical, parameter :: DoBreakUpMessages=.false.
  logical, parameter :: DoImplicitUnusedBlock=.true.

  ! For UseGhostsInProlongation=.false. no ghost cells are used in prolongation for sending
  !   values from coarse to fine blocks.  When there is a load balance, amr, or starting up,
  !   the ghost cells may have bad or unset values.  In this case, correct values for the
  !   prolongation need two message passes.  Even when the values are OK, they may be one
  !   timestep old.  Two message passes will also fix this.
  logical, parameter :: UseGhostsInProlongation=.false.

end module ModMPCells

!==========================================================================
!==========================================================================
!==========================================================================
subroutine message_pass_cells(DoOneLayer,DoFacesOnly,UseMonotoneRestrict,V)
  !
  ! This routine will complete a messagepass of passed in variable V.
  !   It will fill in all ghost cell values where a neighbor in that
  !   direction exists.  This includes refinement level changes and
  !   fills face values and optionally edge and corner values.
  !
  ! NOTE: restriction is the average of 8 finer values, prolongation gives
  !       all 8 new values the same value as coarser cell.
  !
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK, nBLK, iNewGrid, iNewDecomposition
  use ModMPCells
  use ModMpi
  implicit none

  !Subroutine arguements
  real, intent(inout), dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK) :: V
  logical, intent(in) :: DoOneLayer,DoFacesOnly,UseMonotoneRestrict

  !Local variables
  real :: Counter
  integer :: i,j,k, iV,iBLK,iPE, iError
  integer :: nSENDrequests, SENDrequests(maxMessages)
  integer :: nRECVrequests, RECVrequests(maxMessages)
  integer :: MESGstatus(MPI_STATUS_SIZE, maxMessages)

  !------------------------------------------

  ! Check that indices are up to date and expected values of DoOneLayer,DoFacesOnly used
  if(iNewGrid/=iLastGrid .or. iNewDecomposition/=iLastDecomposition &
       .or. (DoOneLayer .neqv. DoOneLayerLast) &
       .or. (DoFacesOnly .neqv. DoFacesOnlyLast)) &
       call mp_cells_set_indices(DoOneLayer,DoFacesOnly)

  ! When neighbor is on the same processor, Collect/Send/Assign are all
  !    done in one step, without intermediate memory use.
  iPE=iProc
  do iV=1,nSend(iPE)
     lS(:)=VSendIlocal(:,iV)
     lR(:)=VRecvIlocal(:,iV)
     if(lS(0)==0)then
        do k=0,lR(6)-lR(5); do j=0,lR(4)-lR(3); do i=0,lR(2)-lR(1)
           V(     lR(1)+i,lR(3)+j,lR(5)+k,lR(7)) = &
                V(lS(1)+i,lS(3)+j,lS(5)+k,lS(7))
        end do; end do; end do
     elseif(lS(0)==1)then
        V(lR(1),lR(3),lR(5),lR(7)) = V(lS(1),lS(3),lS(5),lS(7))
     elseif(lS(0)==2)then
        do k=lR(5),lR(6); do j=lR(3),lR(4); do i=lR(1),lR(2)
           V(i,j,k,lR(7)) = V(lS(1),lS(3),lS(5),lS(7))
        end do; end do; end do
     elseif(lS(0)==3)then
        V(lR(1),lR(3),lR(5),lR(7)) = Inv12*(9.0* &
             V(lS(1),lS(3),lS(5),lS(7))+ &
             V(lS(2),lS(3),lS(5),lS(7))+ &
             V(lS(1),lS(4),lS(5),lS(7))+ &
             V(lS(1),lS(3),lS(6),lS(7)))
     else
        if(UseMonotoneRestrict) call monotone_restrict_indexes(lS)
        Counter=1.0/real( (1+lS(2)-lS(1)) * (1+lS(4)-lS(3)) * (1+lS(6)-ls(5)) )
        V(lR(1),lR(3),lR(5),lR(7)) = Counter*sum(V(lS(1):lS(2),lS(3):lS(4),lS(5):lS(6),lS(7)))
     end if
  end do

  ! Collect values into VSend that need to be passed to other processors
  do iPE=0,nProc-1
     if(iPE==iProc) CYCLE
     if(nSend(iPE)==0) CYCLE
     do iV=nSendStart(iPE)+1,nSendStart(iPE)+nSend(iPE)
        lS(:)=VSendI(:,iV)
        if(lS(0)==1 .or. lS(0)==2)then
           VSend(iV) = V(lS(1),lS(3),lS(5),lS(7))
        elseif(lS(0)==3)then
           VSend(iV) = Inv12*(9.0* &
                V(lS(1),lS(3),lS(5),lS(7))+ &
                V(lS(2),lS(3),lS(5),lS(7))+ &
                V(lS(1),lS(4),lS(5),lS(7))+ &
                V(lS(1),lS(3),lS(6),lS(7)))
        else
           if(UseMonotoneRestrict) call monotone_restrict_indexes(lS)
           Counter=1.0/real( (1+lS(2)-lS(1)) * (1+lS(4)-lS(3)) * (1+lS(6)-ls(5)) )
           VSend(iV) = Counter*sum(V(lS(1):lS(2),lS(3):lS(4),lS(5):lS(6),lS(7)))
        end if
     end do
  end do

  ! Post receives first so that they are ready
  nRECVrequests = 0  
  do iPE=0,nProc-1
     if(iPE == iProc) CYCLE
     if(nRecv(iPE)==0) CYCLE
     if(DoBreakUpMessages)then
        nSends=1+((nRecv(iPE)-1)/MessageSize)
        do i=1,nSends
           itag = nProc*(i-1)+iPE
           nRECVrequests = nRECVrequests + 1
           if(nRECVrequests>maxMessages) call stop_mpi("Too many RECVs in mp_SendValues")
           call MPI_irecv(VRecv(nRecvStart(iPE)+1+((i-1)*MessageSize)), &
                min(MessageSize,nRecv(iPE)-((i-1)*MessageSize)), &
                MPI_REAL,iPE,itag,iComm,RECVrequests(nRECVrequests),iError)
        end do
     else
        itag = iPE
        nRECVrequests = nRECVrequests + 1
        if(nRECVrequests>maxMessages) call stop_mpi("Too many RECVs in mp_SendValues")
        call MPI_irecv(VRecv(nRecvStart(iPE)+1),nRecv(iPE), &
             MPI_REAL,iPE,itag,iComm,RECVrequests(nRECVrequests),iError)
     end if
  end do

  ! Make sure all recv's are posted before using an rsend
  if(DoRSend)then
     call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
  end if

  ! VSend array sent to VRecv array on other PEs
  nSENDrequests = 0
  do iPE=0,nProc-1
     if(iPE == iProc) CYCLE
     if(nSend(iPE)==0) CYCLE
     if(DoBreakUpMessages)then
        nSends=1+((nSend(iPE)-1)/MessageSize)
        do i=1,nSends
           itag = nProc*(i-1)+iProc
           if(DoRSend)then
              call MPI_rsend(VSend(nSendStart(iPE)+1+((i-1)*MessageSize)), &
                   min(MessageSize,nSend(iPE)-((i-1)*MessageSize)), &
                   MPI_REAL,iPE,itag,iComm,iError)
           else
              nSENDrequests = nSENDrequests + 1
              call MPI_isend(VSend(nSendStart(iPE)+1+((i-1)*MessageSize)), &
                   min(MessageSize,nSend(iPE)-((i-1)*MessageSize)), &
                   MPI_REAL,iPE,itag,iComm,SENDrequests(nSENDrequests),iError)
           end if
        end do
     else
        itag = iProc
        if(DoRSend)then
           call MPI_rsend(VSend(nSendStart(iPE)+1),nSend(iPE), &
                MPI_REAL,iPE,itag,iComm,iError)
        else
           nSENDrequests = nSENDrequests + 1
           call MPI_isend(VSend(nSendStart(iPE)+1),nSend(iPE), &
                MPI_REAL,iPE,itag,iComm,SENDrequests(nSENDrequests),iError)
        end if
     end if
  end do

  ! Wait for messages to be received before continuing.
  call MPI_waitall(nRECVrequests, RECVrequests(1), MESGstatus(1,1), iError)

  ! This loop copies the values from VRecv to their destination
  do iPE=0,nProc-1
     if(iPE==iProc) CYCLE
     if(nRecv(iPE)==0) CYCLE
     do iV=nRecvStart(iPE)+1,nRecvStart(iPE)+nRecv(iPE)
        lR(:)=VRecvI(:,iV)
        if(lR(0)==2)then
           do k=lR(5),lR(6); do j=lR(3),lR(4); do i=lR(1),lR(2)
              V(i,j,k,lR(7)) = VRecv(iV)
           end do; end do; end do
        else
           V(lR(1),lR(3),lR(5),lR(7)) = VRecv(iV)
        end if
     end do
  end do

  ! Wait for sent messages to be received before exiting
  if(.not.DoRSend)then
     call MPI_waitall(nSENDrequests, SENDrequests(1), MESGstatus(1,1), iError)
  end if

end subroutine message_pass_cells

!==========================================================================
!==========================================================================
!==========================================================================
subroutine message_pass_cells8(DoOneLayer,DoFacesOnly,UseMonotoneRestrict,nVar,State_VGB)
  !
  ! This routine will complete a messagepass of passed in variables
  !   State_VGB
  !   It will fill in all ghost cell values where a neighbor in that
  !   direction exists.  This includes refinement level changes and
  !   fills face values and optionally edge and corner values.
  !
  ! NOTE: restriction is the average of 8 finer values, prolongation gives
  !       all 8 new values the same value as coarser cell.
  !
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK, nBLK, iNewGrid, iNewDecomposition
  use ModMPCells
  use ModMpi
  use ModVarIndexes,ONLY:nVarBuff=>nVar
  implicit none

  !Subroutine arguements
  integer,intent(in)::nVar
  real, intent(inout), dimension(nVar,-1:nI+2,-1:nJ+2,-1:nK+2,nBLK) :: State_VGB
  logical, intent(in) :: DoOneLayer,DoFacesOnly,UseMonotoneRestrict

  !Local variables
  real :: Counter
  integer :: i,j,k, iV,iBLK,iPE, iError, iVar
  integer :: nSENDrequests, SENDrequests(maxMessages)
  integer :: nRECVrequests, RECVrequests(maxMessages)
  integer :: MESGstatus(MPI_STATUS_SIZE, maxMessages)

  !------------------------------------------

  ! If sending corners and you want to limit the amount of memory needed,
  !    then send values one at a time to limit memory use
  if(DoLimitCornerMemory .and. .not.DoFacesOnly) then
     do iVar=1,nVar
        call message_pass_cells(&
             DoOneLayer,DoFacesOnly,UseMonotoneRestrict,State_VGB(iVar,:,:,:,:))
     end do
     return
  end if

  ! Check that indices are up to date and expected values of DoOneLayer,DoFacesOnly used
  if(iNewGrid/=iLastGrid .or. iNewDecomposition/=iLastDecomposition &
       .or. (DoOneLayer .neqv. DoOneLayerLast) &
       .or. (DoFacesOnly .neqv. DoFacesOnlyLast)) &
       call mp_cells_set_indices(DoOneLayer,DoFacesOnly)

  ! When neighbor is on the same processor, Collect/Send/Assign are all
  !    done in one step, without intermediate memory use.
  iPE=iProc
  do iV=1,nSend(iPE)
     lS(:)=VSendIlocal(:,iV)
     lR(:)=VRecvIlocal(:,iV)
     if(lS(0)==0)then
        do k=0,lR(6)-lR(5); do j=0,lR(4)-lR(3); do i=0,lR(2)-lR(1)
           State_VGB(     :,lR(1)+i,lR(3)+j,lR(5)+k,lR(7)) = &
                State_VGB(:,lS(1)+i,lS(3)+j,lS(5)+k,lS(7))
        end do; end do; end do
     elseif(lS(0)==1)then
        State_VGB(     :,lR(1),lR(3),lR(5),lR(7)) = &
             State_VGB(:,lS(1),lS(3),lS(5),lS(7))
     elseif(lS(0)==2)then
        do k=lR(5),lR(6); do j=lR(3),lR(4); do i=lR(1),lR(2)
           State_VGB(:,i,j,k,lR(7)) = State_VGB(:,lS(1),lS(3),lS(5),lS(7))
        end do; end do; end do
     elseif(lS(0)==3)then
        State_VGB(:,lR(1),lR(3),lR(5),lR(7)) = Inv12*(9.0* &
             State_VGB(:,lS(1),lS(3),lS(5),lS(7))+ &
             State_VGB(:,lS(2),lS(3),lS(5),lS(7))+ &
             State_VGB(:,lS(1),lS(4),lS(5),lS(7))+ &
             State_VGB(:,lS(1),lS(3),lS(6),lS(7)))
     else
        if(UseMonotoneRestrict) call monotone_restrict_indexes(lS)
        Counter=1.0/real( (1+lS(2)-lS(1)) * (1+lS(4)-lS(3)) * (1+lS(6)-ls(5)) )
        do iVar=1,nVar
           State_VGB(iVar,lR(1),lR(3),lR(5),lR(7)) =  Counter*sum(&
               State_VGB(iVar,lS(1):lS(2),lS(3):lS(4),lS(5):lS(6),lS(7)))
        end do
     end if
  end do

  ! Collect values into VSend that need to be passed to other processors
  do iPE=0,nProc-1
     if(iPE==iProc) CYCLE
     if(nSend(iPE)==0) CYCLE
     do iV=nSendStart(iPE)+1,nSendStart(iPE)+nSend(iPE)
        lS(:)=VSendI(:,iV)
        if(lS(0)==0 .or. lS(0)==1)then
           VSend8(1:nVar,iV) =  State_VGB(:,lS(1),lS(3),lS(5),lS(7))
        elseif(lS(0)==3)then
           VSend8(1:nVar,iV) = Inv12*(9.0* &
                State_VGB(:,lS(1),lS(3),lS(5),lS(7))+ &
                State_VGB(:,lS(2),lS(3),lS(5),lS(7))+ &
                State_VGB(:,lS(1),lS(4),lS(5),lS(7))+ &
                State_VGB(:,lS(1),lS(3),lS(6),lS(7)))
        else
           if(UseMonotoneRestrict) call monotone_restrict_indexes(lS)
           Counter=1.0/real( (1+lS(2)-lS(1)) * (1+lS(4)-lS(3)) * (1+lS(6)-ls(5)) )
           do iVar=1,nVar
              VSend8(iVar,iV) = Counter*sum(&
                   State_VGB(iVar,lS(1):lS(2),lS(3):lS(4),lS(5):lS(6),lS(7)))
           end do
        end if
     end do
  end do

  ! Post receives first so that they are ready
  nRECVrequests = 0  
  do iPE=0,nProc-1
     if(iPE == iProc) CYCLE
     if(nRecv(iPE)==0) CYCLE
     if(DoBreakUpMessages)then
        nSends=1+((nRecv(iPE)-1)/MessageSize)
        do i=1,nSends
           itag = nProc*(i-1)+iPE
           nRECVrequests = nRECVrequests + 1
           if(nRECVrequests>maxMessages) call stop_mpi("Too many RECVs in mp_SendValues")
           call MPI_irecv(VRecv8(1,nRecvStart(iPE)+1+((i-1)*MessageSize)), &
                nVarBuff * min(MessageSize,nRecv(iPE)-((i-1)*MessageSize)), &
                MPI_REAL,iPE,itag,iComm,RECVrequests(nRECVrequests),iError)
        end do
     else
        itag = iPE
        nRECVrequests = nRECVrequests + 1
        if(nRECVrequests>maxMessages) call stop_mpi("Too many RECVs in mp_SendValues")
        call MPI_irecv(VRecv8(1,nRecvStart(iPE)+1),nVarBuff * nRecv(iPE), &
             MPI_REAL,iPE,itag,iComm,RECVrequests(nRECVrequests),iError)
     end if
  end do

  ! Make sure all recv's are posted before using an rsend
  if(DoRSend)then
     call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
  end if

  ! VSend array sent to VRecv array on other PEs
  nSENDrequests = 0
  do iPE=0,nProc-1
     if(iPE == iProc) CYCLE
     if(nSend(iPE)==0) CYCLE
     if(DoBreakUpMessages)then
        nSends=1+((nSend(iPE)-1)/MessageSize)
        do i=1,nSends
           itag = nProc*(i-1)+iProc
           if(DoRSend)then
              call MPI_rsend(VSend8(1,nSendStart(iPE)+1+((i-1)*MessageSize)), &
                   nVarBuff * min(MessageSize,nSend(iPE)-((i-1)*MessageSize)), &
                   MPI_REAL,iPE,itag,iComm,iError)
           else
              nSENDrequests = nSENDrequests + 1
              call MPI_isend(VSend8(1,nSendStart(iPE)+1+((i-1)*MessageSize)), &
                   nVarBuff * min(MessageSize,nSend(iPE)-((i-1)*MessageSize)), &
                   MPI_REAL,iPE,itag,iComm,SENDrequests(nSENDrequests),iError)
           end if
        end do
     else
        itag = iProc
        if(DoRSend)then
           call MPI_rsend(VSend8(1,nSendStart(iPE)+1), &
                nVarBuff * nSend(iPE), &
                MPI_REAL,iPE,itag,iComm,iError)
        else
           nSENDrequests = nSENDrequests + 1
           call MPI_isend(VSend8(1,nSendStart(iPE)+1), &
                nVarBuff * nSend(iPE), &
                MPI_REAL,iPE,itag,iComm,SENDrequests(nSENDrequests),iError)
        end if
     end if
  end do

  ! Wait for messages to be received before continuing.
  call MPI_waitall(nRECVrequests, RECVrequests(1), MESGstatus(1,1), iError)

  ! This loop copies the values from VRecv to their destination
  do iPE=0,nProc-1
     if(iPE==iProc) CYCLE
     if(nRecv(iPE)==0) CYCLE
     do iV=nRecvStart(iPE)+1,nRecvStart(iPE)+nRecv(iPE)
        lR(:)=VRecvI(:,iV)
        if(lR(0)==2)then
           do k=lR(5),lR(6); do j=lR(3),lR(4); do i=lR(1),lR(2)
              State_VGB(:,i,j,k,lR(7)) = VRecv8(1:nVar,iV)
           end do; end do; end do
        else
           State_VGB(:,lR(1),lR(3),lR(5),lR(7)) = VRecv8(1:nVar,iV)
        end if
     end do
  end do

  ! Wait for sent messages to be received before exiting
  if(.not.DoRSend)then
     call MPI_waitall(nSENDrequests, SENDrequests(1), MESGstatus(1,1), iError)
  end if

end subroutine message_pass_cells8
!==========================================================================
subroutine message_pass_boundary_cells(UseMonotoneRestrict)
  !
  ! This routine will complete a messagepass of boundary cell flags.
  !
  use ModProcMH
  use ModBoundaryCells
  use ModMPCells
  use ModNumConst
  use ModMpi
  implicit none

  !Subroutine arguements
  logical, intent(in) :: UseMonotoneRestrict

  !Local variables
  integer :: iBoundary,i,j,k, iV,iPE, iError
  integer :: nSENDrequests, SENDrequests(maxMessages)
  integer :: nRECVrequests, RECVrequests(maxMessages)
  integer :: MESGstatus(MPI_STATUS_SIZE, maxMessages)

  !logical buffers
  logical, dimension(:,:), allocatable, save :: IsBuffSend_II,IsBuffRecv_II
  !------------------------------------------

  if(.not.SaveBoundaryCells)return
  if(numsend>0)then
     allocate( IsBuffSend_II(MinBoundarySaved:MaxBoundarySaved,numSend),&
          stat=iError )
     call alloc_check(iError,"IsBuffSend_II")
  end if

  if(numrecv>0)then
     allocate( IsBuffRecv_II(MinBoundarySaved:MaxBoundarySaved,numRecv),&
          stat=iError );   
     call alloc_check(iError,"IsBuffRecv_II")
  end if

  ! When neighbor is on the same processor, Collect/Send/Assign are all
  !    done in one step, without intermediate memory use.
  iPE=iProc
  do iV=1,nSend(iPE)
     lS(:)=VSendIlocal(:,iV)
     lR(:)=VRecvIlocal(:,iV)
     if(lS(0)==0)then
        do k=0,lR(6)-lR(5); do j=0,lR(4)-lR(3); do i=0,lR(2)-lR(1)
           IsBoundaryCell_IGB(     :,lR(1)+i,lR(3)+j,lR(5)+k,lR(7)) = &
                IsBoundaryCell_IGB(:,lS(1)+i,lS(3)+j,lS(5)+k,lS(7))
        end do; end do; end do
     elseif(lS(0)==1)then
        IsBoundaryCell_IGB(     :,lR(1),lR(3),lR(5),lR(7)) = &
             IsBoundaryCell_IGB(:,lS(1),lS(3),lS(5),lS(7))
     elseif(lS(0)==2)then
        do k=lR(5),lR(6); do j=lR(3),lR(4); do i=lR(1),lR(2)
           IsBoundaryCell_IGB(:,i,j,k,lR(7)) = IsBoundaryCell_IGB(:,lS(1),lS(3),lS(5),lS(7))
        end do; end do; end do
     elseif(lS(0)==3)then
        IsBoundaryCell_IGB(:,lR(1),lR(3),lR(5),lR(7)) = &
             IsBoundaryCell_IGB(:,lS(1),lS(3),lS(5),lS(7)).or. &
             IsBoundaryCell_IGB(:,lS(2),lS(3),lS(5),lS(7)).or. &
             IsBoundaryCell_IGB(:,lS(1),lS(4),lS(5),lS(7)).or. &
             IsBoundaryCell_IGB(:,lS(1),lS(3),lS(6),lS(7))
     else
        if(UseMonotoneRestrict) call monotone_restrict_indexes(lS)
        do iBoundary=MinBoundarySaved,MaxBoundarySaved
           IsBoundaryCell_IGB(iBoundary,lR(1),lR(3),lR(5),lR(7)) = any( &
                IsBoundaryCell_IGB(iBoundary,lS(1):lS(2),lS(3):lS(4),lS(5):lS(6),lS(7)))
        end do
     end if
  end do

  ! Collect values into VSend that need to be passed to other processors
  do iPE=0,nProc-1
     if(iPE==iProc) CYCLE
     if(nSend(iPE)==0) CYCLE
     do iV=nSendStart(iPE)+1,nSendStart(iPE)+nSend(iPE)
        lS(:)=VSendI(:,iV)
        if(lS(0)==0 .or. lS(0)==1)then
           IsBuffSend_II(:,iV) = IsBoundaryCell_IGB(:,lS(1),lS(3),lS(5),lS(7))
        elseif(lS(0)==3)then
           IsBuffSend_II(:,iV) =                      &
                IsBoundaryCell_IGB(:,lS(1),lS(3),lS(5),lS(7)).or. &
                IsBoundaryCell_IGB(:,lS(2),lS(3),lS(5),lS(7)).or. &
                IsBoundaryCell_IGB(:,lS(1),lS(4),lS(5),lS(7)).or. &
                IsBoundaryCell_IGB(:,lS(1),lS(3),lS(6),lS(7))
        else
           if(UseMonotoneRestrict) call monotone_restrict_indexes(lS)
           do iBoundary=MinBoundarySaved,MaxBoundarySaved
              IsBuffSend_II(iBoundary,iV) = any( &
                   IsBoundaryCell_IGB(iBoundary,lS(1):lS(2),lS(3):lS(4),lS(5):lS(6),lS(7)))
           end do
        end if
     end do
  end do

  ! Post receives first so that they are ready
  nRECVrequests = 0  
  do iPE=0,nProc-1
     if(iPE == iProc) CYCLE
     if(nRecv(iPE)==0) CYCLE
     if(DoBreakUpMessages)then
        nSends=1+((nRecv(iPE)-1)/MessageSize)
        do i=1,nSends
           itag = nProc*(i-1)+iPE
           nRECVrequests = nRECVrequests + 1
           if(nRECVrequests>maxMessages) call stop_mpi("Too many RECVs in mp_SendValues")
           call MPI_irecv(IsBuffRecv_II(MinBoundarySaved,nRecvStart(iPE)+1+((i-1)*MessageSize)),&
                nBoundarySaved*min(MessageSize,nRecv(iPE)-((i-1)*MessageSize)), &
                MPI_LOGICAL,iPE,itag,iComm,RECVrequests(nRECVrequests),iError)
        end do
     else
        itag = iPE
        nRECVrequests = nRECVrequests + 1
        if(nRECVrequests>maxMessages) call stop_mpi("Too many RECVs in mp_SendValues")
        call MPI_irecv(IsBuffRecv_II(MinBoundarySaved,nRecvStart(iPE)+1),nBoundarySaved*nRecv(iPE), &
             MPI_LOGICAL,iPE,itag,iComm,RECVrequests(nRECVrequests),iError)
     end if
  end do

  ! Make sure all recv's are posted before using an rsend
  if(DoRSend)then
     call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
  end if

  ! VSend array sent to VRecv array on other PEs
  nSENDrequests = 0
  do iPE=0,nProc-1
     if(iPE == iProc) CYCLE
     if(nSend(iPE)==0) CYCLE
     if(DoBreakUpMessages)then
        nSends=1+((nSend(iPE)-1)/MessageSize)
        do i=1,nSends
           itag = nProc*(i-1)+iProc
           if(DoRSend)then
              call MPI_rsend(IsBuffSend_II(MinBoundarySaved,nSendStart(iPE)+1+((i-1)*MessageSize)),&
                   nBoundarySaved*min(MessageSize,nSend(iPE)-((i-1)*MessageSize)), &
                   MPI_LOGICAL,iPE,itag,iComm,iError)
           else
              nSENDrequests = nSENDrequests + 1
              call MPI_isend(IsBuffSend_II(MinBoundarySaved,nSendStart(iPE)+1+((i-1)*MessageSize)),&
                   nBoundarySaved*min(MessageSize,nSend(iPE)-((i-1)*MessageSize)), &
                   MPI_LOGICAL,iPE,itag,iComm,SENDrequests(nSENDrequests),iError)
           end if
        end do
     else
        itag = iProc
        nSENDrequests = nSENDrequests + 1
        if(DoRSend)then
           call MPI_rsend(IsBuffSend_II(MinBoundarySaved ,nSendStart(iPE)+1),&
                nBoundarySaved*nSend(iPE), &
                MPI_LOGICAL,iPE,itag,iComm,iError)
        else
           call MPI_isend(IsBuffSend_II(MinBoundarySaved,nSendStart(iPE)+1),&
                nBoundarySaved*nSend(iPE), &
                MPI_LOGICAL,iPE,itag,iComm,SENDrequests(nSENDrequests),iError)
        end if
     end if
  end do

  ! Wait for messages to be received before continuing.
  call MPI_waitall(nRECVrequests, RECVrequests(1), MESGstatus(1,1), iError)

  ! This loop copies the values from VRecv to their destination
  do iPE=0,nProc-1
     if(iPE==iProc) CYCLE
     if(nRecv(iPE)==0) CYCLE
     do iV=nRecvStart(iPE)+1,nRecvStart(iPE)+nRecv(iPE)
        lR(:)=VRecvI(:,iV)
        if(lR(0)==2)then
           do k=lR(5),lR(6); do j=lR(3),lR(4); do i=lR(1),lR(2)
              IsBoundaryCell_IGB(:,i,j,k,lR(7)) = IsBuffRecv_II(:,iV)
           end do; end do; end do
        else
           IsBoundaryCell_IGB(:,lR(1),lR(3),lR(5),lR(7)) = IsBuffRecv_II(:,iV)
        end if
     end do
  end do

  ! Wait for sent messages to be received before exiting
  if(.not.DoRSend)then
     call MPI_waitall(nSENDrequests, SENDrequests(1), MESGstatus(1,1), iError)
  end if
  if(allocated(IsBuffSend_II)) deallocate(IsBuffSend_II)
  if(allocated(IsBuffRecv_II)) deallocate(IsBuffRecv_II)
end subroutine message_pass_boundary_cells

!========================================================================
!========================================================================
!========================================================================
subroutine monotone_restrict_indexes(lS)
  use ModMain, ONLY : nI,nJ,nK
  integer, intent(inout),dimension(0:7)::lS
  !\
  ! This subroutine insures "monotone" restriction. The index range is changed near the
  ! refinement resolution so that only 4 (rather than 8) correspondent values got from the FINER block
  ! are used to form the ghostcell value in the coarsen block 


  select case(lS(0))
  case(-1) ! send from finer block to the coarser WEST block
     if(lS(2)==nI) lS(1)=lS(2)
  case(-2) ! send from finer block to the coarser EAST block
     if(lS(1)==1) lS(2)=lS(1)
  case(-3) ! send from finer block to the coarser NORTH block
     if(lS(4)==nJ) lS(3)=lS(4)
  case(-4) ! send from finer block to the coarser SOUTH block
     if(lS(3)==1) lS(4)=lS(3)
  case(-5) ! send from finer block to the coarser TOP block
     if(lS(6)==nK) lS(5)=lS(6)
  case(-6) ! send from finer block to the coarser BOT block
     if(lS(5)==1) lS(6)=lS(5)
  end select
end subroutine monotone_restrict_indexes
!==========================================================================
!==========================================================================
!==========================================================================
subroutine mp_cells_set_indices(DoOneLayer,DoFacesOnly)
  !
  !
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK, nBLK, iNewGrid, iNewDecomposition
  use ModVarIndexes,ONLY:nVar
  use ModMPCells
  use ModMpi
  implicit none

  !Subroutine arguements
  logical, intent(in) :: DoOneLayer,DoFacesOnly

  !Local variables
  integer :: iPE, iError
  !------------------------------------------

  DoOneLayer_D = DoOneLayer 
  
  !If corners are to be send and DoOneLayerCorner is set to .true.,
  !then only one layer of ghostcells is send to fill in corners,
 
  DoOneLayer_D(7:26) = DoOneLayer_D(7:26).or.DoOneLayerCorner

  DoOneLayerLast     = DoOneLayer
  DoFacesOnlyLast    = DoFacesOnly
  iLastGrid          = iNewGrid
  iLastDecomposition = iNewDecomposition

  call mp_allocate_cell_arrays1
  call mp_build_cell_indices(.true.)

  call mp_allocate_cell_arrays2
  call mp_build_cell_indices(.false.)

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------


contains

  !==========================================================================
  subroutine mp_allocate_cell_arrays1
    !
    ! Allocate initial arrays for counting the number of sends and receives
    !

    !------------------------------------------

    ! Allocate new memory only if not already allocated, then initialize

    if(.not.allocated(nSend))then
       allocate( nSend(0:nProc-1), stat=iError )
       call alloc_check(iError,"nSend")
    end if
    nSend=0

    if(.not.allocated(nRecv))then
       allocate( nRecv(0:nProc-1), stat=iError )
       call alloc_check(iError,"nRecv")
    end if
    nRecv=0

    if(.not.allocated(nSendStart))then
       allocate( nSendStart(0:nProc-1), stat=iError )
       call alloc_check(iError,"nSendStart")
    end if
    nSendStart=0

    if(.not.allocated(nRecvStart))then
       allocate( nRecvStart(0:nProc-1), stat=iError )
       call alloc_check(iError,"nRecvStart")
    end if
    nRecvStart=0

  end subroutine mp_allocate_cell_arrays1

  !==========================================================================
  subroutine mp_allocate_cell_arrays2
    !
    ! Allocate arrays to hold the index mapping and send/recv values based on
    !    how many sends/recvs are needed
    !

    !------------------------------------------

    numCopy=nSend(iProc)
    numSend=sum(nSend)-nSend(iProc)
    numRecv=sum(nRecv)-nRecv(iProc)
    numSendRecv = max(1,max(numSend,numRecv))

    nSendStart=0
    nRecvStart=0
    do iPE=0,nProc-2
       if(iPE==iProc)then
          nSendStart(iPE+1)=nSendStart(iPE)
          nRecvStart(iPE+1)=nRecvStart(iPE)
       else
          nSendStart(iPE+1)=nSendStart(iPE)+nSend(iPE)
          nRecvStart(iPE+1)=nRecvStart(iPE)+nRecv(iPE)
       end if
    end do

!!$    write(*,*)'MPC: PE=',iProc,' SEND:',nSend,' RECV:',nRecv,' SUMMARY:',numSend,numRecv,numCopy

    ! Re-Allocate memory if past use was not large enough

    if(  numSendRecv > numSendRecvMax .or. &
         numSend     > numSendMax     .or. &
         numRecv     > numRecvMax     .or. &
         numCopy     > numCopyMax&
         .or.( &
         (.not.(DoLimitCornerMemory .and. .not.DoFacesOnly))&
         .and.(numSend>numSendMax8.or.numRecv>numRecvMax8)))then

       ! Deallocate old memory and allocate new memory

       if(allocated(VSendI)) deallocate(VSendI)
       allocate( VSendI(0:7,numSendRecv), stat=iError )
       call alloc_check(iError,"VSendI")

       if(allocated(VRecvI)) deallocate(VRecvI)
       allocate( VRecvI(0:7,numSendRecv), stat=iError )
       call alloc_check(iError,"VRecvI")

       if(allocated(VSendIlocal)) deallocate(VSendIlocal)
       allocate( VSendIlocal(0:7,numCopy), stat=iError )
       call alloc_check(iError,"VSendIlocal")

       if(allocated(VRecvIlocal)) deallocate(VRecvIlocal)
       allocate( VRecvIlocal(0:7,numCopy), stat=iError )
       call alloc_check(iError,"VRecvIlocal")

       if(allocated(VSend )) deallocate(VSend )
       allocate( VSend(numSend), stat=iError )
       call alloc_check(iError,"VSend")

       if(allocated(VRecv )) deallocate(VRecv )
       allocate( VRecv(numRecv), stat=iError )
       call alloc_check(iError,"VRecv")

       ! If sending corners and you want to limit the amount of memory needed,
       !    then send values one at a time to limit memory use
       !    and don't allocate this memory.
       if(.not.(DoLimitCornerMemory .and. .not.DoFacesOnly)) then
          if(allocated(VSend8)) deallocate(VSend8)
          allocate( VSend8(nVar,numSend), stat=iError )
          call alloc_check(iError,"VSend8")

          if(allocated(VRecv8)) deallocate(VRecv8)
          allocate( VRecv8(nVar,numRecv), stat=iError )
          call alloc_check(iError,"VRecv8")

          numSendMax8=numSend
          numRecvMax8=numRecv
       end if

       ! Update new max values

       numCopyMax = numCopy
       numSendMax = numSend
       numRecvMax = numRecv
       numSendRecvMax = numSendRecv

    end if

    ! Initialize memory

    VSendI=-99
    VRecvI=-99
    VSendIlocal=-99
    VRecvIlocal=-99
    VSend=0
    VRecv=0

    ! If sending corners and you want to limit the amount of memory needed,
    !    then send values one at a time to limit memory use
    !    and don't allocate this memory.
    ! We need to check if allocated because "previous use" test above may not have
    !    allocated the memory but may be enough for other message passing.
    if(.not.(DoLimitCornerMemory .and. .not.DoFacesOnly)) then
       if(.not.allocated(VSend8))then
          allocate( VSend8(nVar,numSend), stat=iError )
          call alloc_check(iError,"VSend8")
          numSendMax8=numSend
       end if
       VSend8=0
       if(.not.allocated(VRecv8))then
          allocate( VRecv8(nVar,numRecv), stat=iError )
          call alloc_check(iError,"VRecv8")
          numRecvMax8=numRecv
       end if
       VRecv8=0
    end if

  end subroutine mp_allocate_cell_arrays2

end subroutine mp_cells_set_indices

!==========================================================================
!==========================================================================
!==========================================================================
subroutine mp_build_cell_indices(JustCount)
  !
  !
  use ModProcMH
  use ModMain
  use ModOctree
  use ModMPCells
  use ModAMR, ONLY : unusedBlock_BP
  use ModGeometry,ONLY:TypeGeometry                
  use ModParallel,ONLY:NOBLK         
  use ModPolarNeighbor
  implicit none

  !Subroutine arguements
  logical, intent(in) :: JustCount

  !Local variables
  integer :: iBLK,iPE, iChild, idir, sSubF
  integer :: i1S,i2S, j1S,j2S, k1S,k2S, i1R,i2R, j1R,j2R, k1R,k2R
  integer :: sS
  integer,dimension(3)::iMinS_D,iMaxS_D,iMinR_D,iMaxR_D,nDuplicate_D
  integer::nLayerS,nLayerR
 
  
  !For passing message across the pole some indexes should be flipped
  !Typical example:
  !Usual situation: if in the loop for sending the limits for k are 1,2
  !SEND:do k=1,2;........
  !then in the loop for receiving, the limits are nK+1,nK+2
  !RECEIVE:do k=nK+1,nK+2;...
  !However,
  !for message pass across the pole of the spherical grid the RECEIVE loop becomes:
  !RECEIVE:do k=0,-1,-1;...
  !So we need a logical to figure out if the message pass occurs across the pole:

  logical::IsPolarBlockS=.false.,DoMPassAcrossPole=.false.
  integer::iDirPole,iLoopPole

  !as well as the integers for strides which are negative is the indexes are flipped


  !j index always relates to a cycling coordinate around the pole, so that this index
  !never flips and never requires negative stride

  integer::iDirMax=26

  integer, dimension(26) :: nsubF
  integer, dimension(26,3) :: dLOOP
 


  integer :: iLevelR    !level of the receiving block
  integer, dimension(4) :: neighborPE,neighborBLK,neighborCHILD

  !------------------------------------------

  ! face=1-6, edge=7-18, corner=19-26
  !    1   2   3   4   5   6    7   8   9  10  11  12  13  14  15  16  17  18   19  20  21  22  23  24  25  26
  !    W   E   N   S   T   B   WN  ES  WS  EN  NT  SB  NB  ST  TW  BE  TE  BW  WNT ESB WNB EST WST ENB WSB ENT
  data nsubF / &
       4,  4,  4,  4,  4,  4,   2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,   1,  1,  1,  1,  1,  1,  1,  1 /

  data dLOOP / &
       1, -1,  0,  0,  0,  0,   1, -1,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,   1, -1,  1, -1,  1, -1,  1, -1, &
       0,  0,  1, -1,  0,  0,   1, -1, -1,  1,  1, -1,  1, -1,  0,  0,  0,  0,   1, -1,  1, -1, -1,  1, -1,  1, &
       0,  0,  0,  0,  1, -1,   0,  0,  0,  0,  1, -1, -1,  1,  1, -1,  1, -1,   1, -1, -1,  1,  1, -1, -1,  1 /

  nSend=0
  nRecv=0


  if(DoFacesOnlyLast)then
     sS=0 ;  iDirMax=6    
  else
     sS=1 ;  iDirMax=26
  end if

  if(.not.DoOneCoarserLayer) then
     iTwoOrOneForTwoCoarserLayers=1
     iZeroOrOneForTwoCoarserLayers=1
  else
     iTwoOrOneForTwoCoarserLayers=2
     iZeroOrOneForTwoCoarserLayers=0
  end if


  do iPE = 0,nProc-1
     do iBLK = 1,nBLK
        if (.not.associated(global_block_ptrs(iBLK, iPE+1) % ptr)) CYCLE

        if (.not.global_block_ptrs(iBLK, iPE+1) % ptr % used) CYCLE

        if (DoImplicitUnusedBlock)then
           if (unusedBlock_BP(iBLK,iPE)) CYCLE
        end if

        !valid used block found, setup indices
        iCHILD = global_block_ptrs(iBLK, iPE+1) % ptr % child_number

        call check_pole_inside_blkpe(iBlk,iPE,&
             IsPolarBlockS,iDirPole,iLoopPole)

        do idir=1,iDirMax
           DoMPassAcrossPole=&
                IsPolarBlockS.and.dLoop(iDir,iDirPole)==iLoopPole
           if(DoMPassAcrossPole)then
              call tree_neighbor_across_pole(iPE,iBLK,&
                   dLOOP(idir,1),dLOOP(idir,2),dLOOP(idir,3),iDirPole, &
                   neighborPE,neighborBLK,neighborCHILD,iLevelR)
           else
              call treeNeighbor(iPE,iBLK,&
                   dLOOP(idir,1),dLOOP(idir,2),dLOOP(idir,3), &
                   neighborPE,neighborBLK,neighborCHILD,iLevelR)
           end if
           nDuplicate_D=2
           select case(iLevelR)
           case(0)
              !Build indices for send to same block level
              sSubF=0
              if(DoOneLayer_D(iDir))then
                 nLayerS = 1 ; nLayerR = 1
              else
                 nLayerS = 2 ; nLayerR = 2
              end if
              call set_indices(&
                   iDirS2R_D=dLoop(iDir,:), &
                   nLayerS=nLayerS,nLayerR=nLayerR,&
                   iMinS_D=iMinS_D,iMaxS_D=iMaxS_D,&
                   iMinR_D=iMinR_D,iMaxR_D=iMaxR_D,&
                   iLevelR=iLevelR)
              call build_i
           case(1)
              !Build indices for send to coarser block level
              sSubF=-1
              if(is_not_at_face(&
                   iDirC2F_D=-dLoop(iDir,:),iChild=iChild))CYCLE
              if(DoOneLayer_D(iDir))then
                 nLayerS = 2 ; nLayerR = 1
              else
                 nLayerS = 4 ; nLayerR = 2
              end if
              call set_indices(&
                   iDirS2R_D=dLoop(iDir,:), &
                   nLayerS=nLayerS,nLayerR=nLayerR,&
                   iMinS_D=iMinS_D,iMaxS_D=iMaxS_D,&
                   iMinR_D=iMinR_D,iMaxR_D=iMaxR_D,&
                   iLevelR=iLevelR,iChild=iChild)
              call build_i
           case(-1)
              !Build indices for send to finer block level
              nLayerS = 1 ; nLayerR = 1 
              if((.not.DoOneLayer_D(iDir)))then
                 nLayerS = 1 + iZeroOrOneForTwoCoarserLayers
                 nLayerR = 2
                 where(dLoop(iDir,:)/=0)
                    nDuplicate_D=iTwoOrOneForTwoCoarserLayers
                 end where
              end if
              do sSubF=1,nsubF(idir)
                 call set_indices(&
                      iDirS2R_D=dLoop(iDir,:), &
                      nLayerS=nLayerS,nLayerR=nLayerR,&
                      iMinS_D=iMinS_D,iMaxS_D=iMaxS_D,&
                      iMinR_D=iMinR_D,iMaxR_D=iMaxR_D,&
                      iLevelR=iLevelR,&
                      iChild =neighborCHILD(sSubf),&
                      iChild1=neighborCHILD(  1  ),&
                      nExtraCell=sS)
                 call build_i
              end do
           case(NOBLK)                   
           end select
        end do
     end do
  end do


contains
  subroutine build_i
    !
    !
    !Local variables
    integer :: i,j,k, n, iAdd,jAdd,kAdd
    integer :: nborPE, nborBLK
    integer::iStrideR=1,kStrideR=1
    !------------------------------------------
    iStrideR=1; kStrideR=1
   
    nborPE  = neighborPE (max(1,sSubF))
    nborBLK = neighborBLK(max(1,sSubF))

    if (DoImplicitUnusedBlock)then
       if (unusedBlock_BP(nborBLK,nborPE)) return
    end if

    i1S=iMinS_D(1) ;j1S=iMinS_D(2) ; k1S=iMinS_D(3) 
    i2S=iMaxS_D(1) ;j2S=iMaxS_D(2) ; k2S=iMaxS_D(3) 
    i1R=iMinR_D(1) ;j1R=iMinR_D(2) ; k1R=iMinR_D(3) 
    i2R=iMaxR_D(1) ;j2R=iMaxR_D(2) ; k2R=iMaxR_D(3) 
    nDuplicateI=nDuplicate_D(1)
    nDuplicateJ=nDuplicate_D(2)
    nDuplicateK=nDuplicate_D(3)

    if(DoMPassAcrossPole)then
       select case(iDirPole)
       case(1)
          !flip i idexes for recv:
          i1R=nI+1-i1R; i2R=nI+1-i2R; iStrideR=-1
       case(3)
          !flip k indexes for recv
          k1R=nK+1-k1R; k2R=nK+1-k2R; kStrideR=-1
       case default
          call stop_mpi(&
               'Message_pass_cells: in flipping indexes, unknown idir')
       end select
    end if
    
    if(sSubF == 0)then !Send to same level

       if(iProc == iPE)then  ! Setup sends
          SENDLOOP: do i=i1S,i2S; do j=j1S,j2S; do k=k1S,k2S
             nSend(nborPE)=nSend(nborPE)+1
             if(JustCount .and. iPE==nborPE.and.&
             iStrideR==1.and.kStrideR==1) EXIT SENDLOOP  ! Use for array syntax copy
             if(JustCount) CYCLE

             if(iPE==nborPE)then
                if(iStrideR==-1.or.kStrideR==-1)then
                   ! Old point by point copy
                   n=nSend(nborPE)
                   VSendIlocal(0,n)=1
                   VSendIlocal(1,n)=i; VSendIlocal(3,n)=j; VSendIlocal(5,n)=k
                   VSendIlocal(2,n)=i; VSendIlocal(4,n)=j; VSendIlocal(6,n)=k
                   VSendIlocal(7,n)=iBLK
                else
                   ! For local copies of blocks at same level, we can do array syntax copy
                   n=nSend(nborPE)
                   VSendIlocal(0,n)=0
                   VSendIlocal(1,n)=i1S; VSendIlocal(3,n)=j1S; VSendIlocal(5,n)=k1S
                   VSendIlocal(2,n)=i2S; VSendIlocal(4,n)=j2S; VSendIlocal(6,n)=k2S
                   VSendIlocal(7,n)=iBLK
                   EXIT SENDLOOP
                end if
             else
                n=nSendStart(nborPE)+nSend(nborPE)
                VSendI(0,n)=1
                VSendI(1,n)=i; VSendI(3,n)=j; VSendI(5,n)=k
                VSendI(2,n)=i; VSendI(4,n)=j; VSendI(6,n)=k
                VSendI(7,n)=iBLK
             end if
          end do; end do; end do SENDLOOP
       end if
       if(iProc == nborPE)then  ! Setup recvs
          RECVLOOP: do i=i1R,i2R,iStrideR; do j=j1R,j2R; do k=k1R,k2R,kStrideR
             nRecv(iPE)=nRecv(iPE)+1
             if(JustCount .and. iPE==nborPE.and.iStrideR==1.and.kStrideR==1) EXIT RECVLOOP  ! Use for array syntax copy
             if(JustCount) CYCLE

             if(iPE==nborPE)then
                if(iStrideR==-1.or.kStrideR==-1)then
                   ! Old point by point copy
                   n=nRecv(iPE)
                   VRecvIlocal(1,n)=i; VRecvIlocal(3,n)=j; VRecvIlocal(5,n)=k
                   VRecvIlocal(2,n)=i; VRecvIlocal(4,n)=j; VRecvIlocal(6,n)=k
                   VRecvIlocal(7,n)=nborBLK
                else
                   ! For local copies of blocks at same level, we can do array syntax copy
                   n=nRecv(iPE)
                   VRecvIlocal(1,n)=i1R; VRecvIlocal(3,n)=j1R; VRecvIlocal(5,n)=k1R
                   VRecvIlocal(2,n)=i2R; VRecvIlocal(4,n)=j2R; VRecvIlocal(6,n)=k2R
                   VRecvIlocal(7,n)=nborBLK
                   EXIT RECVLOOP
                end if
             else
                n=nRecvStart(iPE)+nRecv(iPE)
                VRecvI(1,n)=i; VRecvI(3,n)=j; VRecvI(5,n)=k
                VRecvI(2,n)=i; VRecvI(4,n)=j; VRecvI(6,n)=k
                VRecvI(7,n)=nborBLK
             end if
          end do; end do; end do RECVLOOP
       end if

    elseif(sSubF > 0)then !Send to finer level

       if(iProc == iPE)then  ! Setup sends
          if(iCFExchangeType==1)then
             ! Old Style: send 8 separate identical values for 8 fine cells to use
             iAdd=1; jAdd=1; kAdd=1
             if(i1R==i2R) iAdd=0
             if(j1R==j2R) jAdd=0
             if(k1R==k2R) kAdd=0
             do i=2*i1S,2*i2S+iAdd; do j=2*j1S,2*j2S+jAdd; do k=2*k1S,2*k2S+kAdd
                nSend(nborPE)=nSend(nborPE)+1
                if(JustCount) CYCLE

                if(iPE==nborPE)then
                   n=nSend(nborPE)
                   VSendIlocal(0,n)=1
                   VSendIlocal(1,n)=i/2; VSendIlocal(3,n)=j/2; VSendIlocal(5,n)=k/2
                   VSendIlocal(2,n)=i/2; VSendIlocal(4,n)=j/2; VSendIlocal(6,n)=k/2
                   VSendIlocal(7,n)=iBLK
                else
                   n=nSendStart(nborPE)+nSend(nborPE)
                   VSendI(0,n)=1
                   VSendI(1,n)=i/2; VSendI(3,n)=j/2; VSendI(5,n)=k/2
                   VSendI(2,n)=i/2; VSendI(4,n)=j/2; VSendI(6,n)=k/2
                   VSendI(7,n)=iBLK
                end if
             end do; end do; end do

          elseif(iCFExchangeType==2)then
             ! New Style: send 1 value for 8 fine cells to use
             do i=i1S,i2S; do j=j1S,j2S; do k=k1S,k2S
                if(iPE==nborPE)then
                   nSend(nborPE)=nSend(nborPE)+1
                   if(JustCount) CYCLE

                   n=nSend(nborPE)
                   VSendIlocal(0,n)=2
                   VSendIlocal(1,n)=i; VSendIlocal(3,n)=j; VSendIlocal(5,n)=k
                   VSendIlocal(2,n)=i; VSendIlocal(4,n)=j; VSendIlocal(6,n)=k
                   VSendIlocal(7,n)=iBLK
                else
                   nSend(nborPE)=nSend(nborPE)+1
                   if(JustCount) CYCLE

                   n=nSendStart(nborPE)+nSend(nborPE)
                   VSendI(0,n)=2
                   VSendI(1,n)=i; VSendI(3,n)=j; VSendI(5,n)=k
                   VSendI(2,n)=i; VSendI(4,n)=j; VSendI(6,n)=k
                   VSendI(7,n)=iBLK
                end if
             end do; end do; end do

          elseif(iCFExchangeType==3)then
             ! TEST
             iAdd=1; jAdd=1; kAdd=1
             if(i1R==i2R) iAdd=0
             if(j1R==j2R) jAdd=0
             if(k1R==k2R) kAdd=0
             do i=2*i1S,2*i2S+iAdd; do j=2*j1S,2*j2S+jAdd; do k=2*k1S,2*k2S+kAdd
                nSend(nborPE)=nSend(nborPE)+1
                if(JustCount) CYCLE

                if(iPE==nborPE)then
                   n=nSend(nborPE)
                   VSendIlocal(0,n)=3
                   VSendIlocal(1,n)=i/2; VSendIlocal(3,n)=j/2; VSendIlocal(5,n)=k/2
                   VSendIlocal(2,n)=VSendIlocal(1,n)-1
                   if(((VSendIlocal(1,n)/2)*2)==VSendIlocal(1,n)) VSendIlocal(2,n)=VSendIlocal(1,n)+1
                   if(.not.UseGhostsInProlongation) VSendIlocal(2,n)=max(1,min(nI,VSendIlocal(2,n)))
                   VSendIlocal(4,n)=VSendIlocal(1,n)-1
                   if(((VSendIlocal(3,n)/2)*2)==VSendIlocal(3,n)) VSendIlocal(4,n)=VSendIlocal(3,n)+1
                   if(.not.UseGhostsInProlongation) VSendIlocal(4,n)=max(1,min(nJ,VSendIlocal(4,n)))
                   VSendIlocal(6,n)=VSendIlocal(1,n)-1
                   if(((VSendIlocal(5,n)/2)*2)==VSendIlocal(5,n)) VSendIlocal(6,n)=VSendIlocal(5,n)+1
                   if(.not.UseGhostsInProlongation) VSendIlocal(6,n)=max(1,min(nK,VSendIlocal(6,n)))
                   VSendIlocal(7,n)=iBLK
                else
                   n=nSendStart(nborPE)+nSend(nborPE)
                   VSendI(0,n)=3
                   VSendI(1,n)=i/2; VSendI(3,n)=j/2; VSendI(5,n)=k/2
                   VSendI(2,n)=VSendI(1,n)-1
                   if(((VSendI(1,n)/2)*2)==VSendI(1,n)) VSendI(2,n)=VSendI(1,n)+1
                   if(.not.UseGhostsInProlongation) VSendI(2,n)=max(1,min(nI,VSendI(2,n)))
                   VSendI(4,n)=VSendI(1,n)-1
                   if(((VSendI(3,n)/2)*2)==VSendI(3,n)) VSendI(4,n)=VSendI(3,n)+1
                   if(.not.UseGhostsInProlongation) VSendI(4,n)=max(1,min(nJ,VSendI(4,n)))
                   VSendI(6,n)=VSendI(1,n)-1
                   if(((VSendI(5,n)/2)*2)==VSendI(5,n)) VSendI(6,n)=VSendI(5,n)+1
                   if(.not.UseGhostsInProlongation) VSendI(6,n)=max(1,min(nK,VSendI(6,n)))
                   VSendI(7,n)=iBLK
                end if
             end do; end do; end do

          end if
       end if
       if(iProc == nborPE)then  ! Setup recvs
          if(iCFExchangeType==1 .or. iCFExchangeType==3)then
             ! Old Style: receive 8 separate values for 8 fine cells to use
             do i=i1R,i2R,iStrideR; do j=j1R,j2R; do k=k1R,k2R,kStrideR
                nRecv(iPE)=nRecv(iPE)+1
                if(JustCount) CYCLE

                if(iPE==nborPE)then
                   n=nRecv(iPE)
                   VRecvIlocal(1,n)=i; VRecvIlocal(3,n)=j; VRecvIlocal(5,n)=k
                   VRecvIlocal(2,n)=i; VRecvIlocal(4,n)=j; VRecvIlocal(6,n)=k
                   VRecvIlocal(7,n)=nborBLK
                else
                   n=nRecvStart(iPE)+nRecv(iPE)
                   VRecvI(1,n)=i; VRecvI(3,n)=j; VRecvI(5,n)=k
                   VRecvI(2,n)=i; VRecvI(4,n)=j; VRecvI(6,n)=k
                   VRecvI(7,n)=nborBLK
                end if
             end do; end do; end do

          elseif(iCFExchangeType==2)then
             ! New Style: receive 1 value for 8 fine cells to use
             iAdd=iStrideR; jAdd=1; kAdd=kStrideR
             if(i1R==i2R.or.nDuplicateI==1) iAdd=0
             if(j1R==j2R.or.nDuplicateJ==1) jAdd=0
             if(k1R==k2R.or.nDuplicateK==1) kAdd=0
             do i=i1R,i2R,nDuplicateI*iStrideR;do j=j1R,j2R,nDuplicateJ;do k=k1R,k2R,nDuplicateK*kStrideR
                if(iPE==nborPE)then
                   nRecv(iPE)=nRecv(iPE)+1
                   if(JustCount) CYCLE
                   !Consider case kStride=-1. In this case k2R=k1R-1
                   !For nDuplicate=1 automatically kAdd=0, so that we should first take
                   !VRecvLocal(5,n)=k, VRecvLocal(5,n)=k for k=k1R and then for k=k2R
                   !For nDuplicate=2 k=k1R. So we should take
                   !VRecvLocal(5,n)=k2R=k1R-1,VRecvLocal(6,n)=k1R,
                   n=nRecv(iPE)
                   VRecvIlocal(0,n)=2
                   VRecvIlocal(1,n)=min(i,i+iAdd) ; VRecvIlocal(3,n)=j     ; VRecvIlocal(5,n)=min(k,k+kAdd)
                   VRecvIlocal(2,n)=max(i,i+iAdd) ; VRecvIlocal(4,n)=j+jAdd; VRecvIlocal(6,n)=max(k,k+kAdd)
                   VRecvIlocal(7,n)=nborBLK
                else
                   nRecv(iPE)=nRecv(iPE)+1
                   if(JustCount) CYCLE

                   n=nRecvStart(iPE)+nRecv(iPE)
                   VRecvI(0,n)=2
                   VRecvI(1,n)=min(i,i+iAdd)     ; VRecvI(3,n)=j     ; VRecvI(5,n)=min(k,k+kAdd)
                   VRecvI(2,n)=max(i,i+iAdd)     ; VRecvI(4,n)=j+jAdd; VRecvI(6,n)=max(k,k+kAdd)
                   VRecvI(7,n)=nborBLK
                end if
             end do;end do;end do
          end if
       end if

    elseif(sSubF == -1)then !Send to coarser level

       if(iProc == iPE)then  ! Setup sends
          do i=i1S,i2S,2; do j=j1S,j2S,2; do k=k1S,k2S,2
             nSend(nborPE)=nSend(nborPE)+1
             if(JustCount) CYCLE

             if(iPE==nborPE)then
                n=nSend(nborPE)
                VSendIlocal(0,n)=-idir
                VSendIlocal(1,n)=i   ; VSendIlocal(3,n)=j   ; VSendIlocal(5,n)=k
                VSendIlocal(2,n)=i+1 ; VSendIlocal(4,n)=j+1 ; VSendIlocal(6,n)=k+1
                VSendIlocal(7,n)=iBLK
             else
                n=nSendStart(nborPE)+nSend(nborPE)
                VSendI(0,n)=-idir
                VSendI(1,n)=i   ; VSendI(3,n)=j   ; VSendI(5,n)=k
                VSendI(2,n)=i+1 ; VSendI(4,n)=j+1 ; VSendI(6,n)=k+1
                VSendI(7,n)=iBLK
             end if
          end do; end do; end do
       end if
       if(iProc == nborPE)then  ! Setup recvs
          do i=i1R,i2R,iStrideR; do j=j1R,j2R; do k=k1R,k2R,kStrideR
             nRecv(iPE)=nRecv(iPE)+1
             if(JustCount) CYCLE

             if(iPE==nborPE)then
                n=nRecv(iPE)
                VRecvIlocal(1,n)=i; VRecvIlocal(3,n)=j; VRecvIlocal(5,n)=k
                VRecvIlocal(2,n)=i; VRecvIlocal(4,n)=j; VRecvIlocal(6,n)=k
                VRecvIlocal(7,n)=nborBLK
             else
                n=nRecvStart(iPE)+nRecv(iPE)
                VRecvI(1,n)=i; VRecvI(3,n)=j; VRecvI(5,n)=k
                VRecvI(2,n)=i; VRecvI(4,n)=j; VRecvI(6,n)=k
                VRecvI(7,n)=nborBLK
             end if
          end do; end do; end do
       end if

    end if

  end subroutine build_i
end subroutine mp_build_cell_indices

!==========================================================================
!==========================================================================
!==========================================================================
!!$!------------------------------------------------------------------------------------!
!!$!
!!$!    North                      8--------------7
!!$!    y,j    Top               / |            / |
!!$!     ^     z,k             /   |          /   |
!!$!     |   /               /     |        /     |
!!$!     | /     West      5-------+------6       |
!!$!     |-----> x,i       |       |      |       |
!!$!                       |       1------+-------2
!!$!                       |     /        |     /
!!$!                       |   /          |   /
!!$!                       | /            | /
!!$!                       4--------------3 --(front, right, and bottom 3D box corner)
!!$!
!!$!  Point 7 is:  +x, +y, +z, West, North, Top
!!$!
!!$!  Point 4 is:  -x, -y, -z, East, South, Bottom
!!$!
!!$!------------------------------------------------------------------------------------!

!^CFG IF DEBUGGING BEGIN
!==========================================================================
!==========================================================================
!==========================================================================
subroutine testmessage_pass_cells
  !
  ! This routine tests the new message passing for filling ghost cell values
  !   It builds a desired result state in V0 made up of a linear combinations
  !   of X, Y, and Z.  Values are then put into V1(to V8), but only the inner
  !   computational cells.  After message pass, V0 should equal V1(to V8).
  !
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK, nBLK, unusedBLK, iNewGrid
  use ModMPCells
  use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK
  use ModParallel, ONLY : BLKneighborLEV
  use ModMpi
  implicit none

  !Local variables
  integer :: i,j,k, iVar, n, iBLK, d1,d2,d3, i0,i1,i2, j0,j1,j2, k0,k1,k2, iError, iPE
  real, parameter :: f1=0.13, f2=0.71
  real, dimension(:,:,:,:), allocatable :: V0,V1,V2
  real, dimension(:,:,:,:,:), allocatable :: State_VGB

  !------------------------------------------

  write(*,*)' '
  write(*,*)'Testing message_pass_cells, PE=',iProc,'  Starting tests ...'

  allocate( V0(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), stat=iError )
  call alloc_check(iError,"V0")
  allocate( V1(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), stat=iError )
  call alloc_check(iError,"V1")
  allocate( V2(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), stat=iError )
  call alloc_check(iError,"V2")
  allocate(State_VGB(8,-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), stat=iError )
  call alloc_check(iError,"State_VGB")

  V0=-999999.
  V1=V0; V2=V0
  do iBLK=1,nBLK
     if(unusedBLK(iBLK)) CYCLE
     do i=-1,nI+2; do j=-1,nJ+2; do k=-1,nK+2
        V0(i,j,k,iBLK)= &
             x_BLK(i,j,k,iBLK)+f1*y_BLK(i,j,k,iBLK)+f2*z_BLK(i,j,k,iBLK)
     end do; end do; end do
     do d1=-1,1; do d2=-1,1; do d3=-1,1
        select case(d1)
        case(-1)
           i1=-1   ; i2= 0
        case(0)
           i1= 1   ; i2=nI
        case(1)
           i1=nI+1 ; i2=nI+2
        end select

        select case(d2)
        case(-1)
           j1=-1   ; j2= 0
        case(0)
           j1= 1   ; j2=nJ
        case(1)
           j1=nJ+1 ; j2=nJ+2
        end select

        select case(d3)
        case(-1)
           k1=-1   ; k2= 0
        case(0)
           k1= 1   ; k2=nK
        case(1)
           k1=nK+1 ; k2=nK+2
        end select

        if(BLKneighborLEV(d1,d2,d3,iBLK)>0)then
           do i=i1,i2; do j=j1,j2; do k=k1,k2
              i0=2*((i+1)/2) ; j0=2*((j+1)/2) ; k0=2*((k+1)/2)
              V0(i,j,k,iBLK)= 0.125*( &
                   (x_BLK(i0  ,j0  ,k0  ,iBLK)+f1*y_BLK(i0  ,j0  ,k0  ,iBLK)+f2*z_BLK(i0  ,j0  ,k0  ,iBLK)) + &
                   (x_BLK(i0-1,j0  ,k0  ,iBLK)+f1*y_BLK(i0-1,j0  ,k0  ,iBLK)+f2*z_BLK(i0-1,j0  ,k0  ,iBLK)) + &
                   (x_BLK(i0  ,j0-1,k0  ,iBLK)+f1*y_BLK(i0  ,j0-1,k0  ,iBLK)+f2*z_BLK(i0  ,j0-1,k0  ,iBLK)) + &
                   (x_BLK(i0-1,j0-1,k0  ,iBLK)+f1*y_BLK(i0-1,j0-1,k0  ,iBLK)+f2*z_BLK(i0-1,j0-1,k0  ,iBLK)) + &
                   (x_BLK(i0  ,j0  ,k0-1,iBLK)+f1*y_BLK(i0  ,j0  ,k0-1,iBLK)+f2*z_BLK(i0  ,j0  ,k0-1,iBLK)) + &
                   (x_BLK(i0-1,j0  ,k0-1,iBLK)+f1*y_BLK(i0-1,j0  ,k0-1,iBLK)+f2*z_BLK(i0-1,j0  ,k0-1,iBLK)) + &
                   (x_BLK(i0  ,j0-1,k0-1,iBLK)+f1*y_BLK(i0  ,j0-1,k0-1,iBLK)+f2*z_BLK(i0  ,j0-1,k0-1,iBLK)) + &
                   (x_BLK(i0-1,j0-1,k0-1,iBLK)+f1*y_BLK(i0-1,j0-1,k0-1,iBLK)+f2*z_BLK(i0-1,j0-1,k0-1,iBLK)))
           end do; end do; end do
        end if

     end do; end do; end do

     if (all(BLKneighborLEV(:,:,:,iBLK) > -2)) then
        do i=1,nI; do j=1,nJ; do k=1,nK
           V1(i,j,k,iBLK)= &
                x_BLK(i,j,k,iBLK)+f1*y_BLK(i,j,k,iBLK)+f2*z_BLK(i,j,k,iBLK)
        end do; end do; end do
     else
        do i=-1,nI+2; do j=-1,nJ+2; do k=-1,nK+2
           V1(i,j,k,iBLK)= &
                x_BLK(i,j,k,iBLK)+f1*y_BLK(i,j,k,iBLK)+f2*z_BLK(i,j,k,iBLK)
        end do; end do; end do
     end if
  end do

  V2=V1
  do iVar=1,8; State_VGB(iVar,:,:,:,:)=V1;end do

  iNewGrid=iNewGrid+1

  call message_pass_cells(.true.,.true.,.false.,V1)
  write(*,*)'MP SUMMARY: l=1 corner=F PE=',iProc, &
       ' numSendMPI=',numSend,' (numSendCopy=',nSend(iProc),')', &
       ' numRecvMPI=',numRecv,' (numRecvCopy=',nRecv(iProc),')'
  do iPE=0,nProc-1
     write(*,*)' PE=',iProc,' Send/Recv with PE=',iPE,' sending:',nSend(iPE),' receiving:',nRecv(iPE)
  end do

  call message_pass_cells(.true.,.false.,.false.,V1)
  write(*,*)'MP SUMMARY: l=1 corner=T PE=',iProc, &
       ' numSendMPI=',numSend,' (numSendCopy=',nSend(iProc),')', &
       ' numRecvMPI=',numRecv,' (numRecvCopy=',nRecv(iProc),')'
  do iPE=0,nProc-1
     write(*,*)' PE=',iProc,' Send/Recv with PE=',iPE,' sending:',nSend(iPE),' receiving:',nRecv(iPE)
  end do

  call message_pass_cells(.false.,.true.,.false.,V1)
  write(*,*)'MP SUMMARY: l=2 corner=F PE=',iProc, &
       ' numSendMPI=',numSend,' (numSendCopy=',nSend(iProc),')', &
       ' numRecvMPI=',numRecv,' (numRecvCopy=',nRecv(iProc),')'
  do iPE=0,nProc-1
     write(*,*)' PE=',iProc,' Send/Recv with PE=',iPE,' sending:',nSend(iPE),' receiving:',nRecv(iPE)
  end do

  V1=V2
  call message_pass_cells(.false.,.false.,.false.,V1)
  write(*,*)'MP SUMMARY: l=2 corner=T PE=',iProc, &
       ' numSendMPI=',numSend,' (numSendCopy=',nSend(iProc),')', &
       ' numRecvMPI=',numRecv,' (numRecvCopy=',nRecv(iProc),')'
  do iPE=0,nProc-1
     write(*,*)' PE=',iProc,' Send/Recv with PE=',iPE,' sending:',nSend(iPE),' receiving:',nRecv(iPE)
  end do

  if(max(maxval(V1-V0),maxval(V0-V1))>.01)then
     write(*,*)'Testing message_pass_cells, PE=',iProc,&
          ' max difference=', &
          max(maxval(V1-V0),maxval(V0-V1)), &
          ' printing out values and exiting.'
     write(*,*)' '
     do iBLK=1,nBLK
        if(unusedBLK(iBLK)) CYCLE
        n=0
        do i=-1,nI+2; do j=-1,nJ+2; do k=-1,nK+2
           if(abs(V0(i,j,k,iBLK)-V1(i,j,k,iBLK))>.01) n=n+1
        end do; end do; end do
        if(n/=0) write(*,*)' PE=',iProc,' BLK=',iBLK,' has',n,'bad values ...'
!!$        do i=-1,nI+2; do j=-1,nJ+2; do k=-1,nK+2
!!$           if(abs(V0(i,j,k,iBLK)-V1(i,j,k,iBLK))>.01)then
!!$              write(*,*)' PE=',iProc,' BLK=',iBLK,' IJK=',i,j,k,&
!!$                   ' Diff=',abs(V0(i,j,k,iBLK)-V1(i,j,k,iBLK)), &
!!$                   ' Values=',V0(i,j,k,iBLK),V1(i,j,k,iBLK)
!!$           end if
!!$        end do; end do; end do
     end do
     call stop_mpi("testmessage_pass_cells error")
  end if

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------


  call message_pass_cells8(.false.,.false.,.false.,8,State_VGB)

  if(  max(maxval(State_VGB(1,:,:,:,:)-V0),maxval(V0-State_VGB(1,:,:,:,:)),&
           maxval(State_VGB(2,:,:,:,:)-V0),maxval(V0-State_VGB(2,:,:,:,:)),&
           maxval(State_VGB(3,:,:,:,:)-V0),maxval(V0-State_VGB(3,:,:,:,:)),&
           maxval(State_VGB(4,:,:,:,:)-V0),maxval(V0-State_VGB(4,:,:,:,:)),&
           maxval(State_VGB(5,:,:,:,:)-V0),maxval(V0-State_VGB(5,:,:,:,:)),&
           maxval(State_VGB(6,:,:,:,:)-V0),maxval(V0-State_VGB(6,:,:,:,:)),&
           maxval(State_VGB(7,:,:,:,:)-V0),maxval(V0-State_VGB(7,:,:,:,:)),&
           maxval(State_VGB(8,:,:,:,:)-V0),maxval(V0-State_VGB(8,:,:,:,:)))>.01 )then
     write(*,*)'Testing message_pass_cells8, PE=',iProc,&
          ' max difference=', &
          max(maxval(State_VGB(1,:,:,:,:)-V0),maxval(V0-State_VGB(1,:,:,:,:))),&
          max(maxval(State_VGB(2,:,:,:,:)-V0),maxval(V0-State_VGB(2,:,:,:,:))),&
          max(maxval(State_VGB(3,:,:,:,:)-V0),maxval(V0-State_VGB(3,:,:,:,:))),&
          max(maxval(State_VGB(4,:,:,:,:)-V0),maxval(V0-State_VGB(4,:,:,:,:))),&
          max(maxval(State_VGB(5,:,:,:,:)-V0),maxval(V0-State_VGB(5,:,:,:,:))),&
          max(maxval(State_VGB(6,:,:,:,:)-V0),maxval(V0-State_VGB(6,:,:,:,:))),&
          max(maxval(State_VGB(7,:,:,:,:)-V0),maxval(V0-State_VGB(7,:,:,:,:))),&
          max(maxval(State_VGB(8,:,:,:,:)-V0),maxval(V0-State_VGB(8,:,:,:,:))),&
          ' printing out values and exiting.'
     write(*,*)' '
     do iBLK=1,nBLK
        if(unusedBLK(iBLK)) CYCLE
        n=0
        do i=-1,nI+2; do j=-1,nJ+2; do k=-1,nK+2
           if(abs(V0(i,j,k,iBLK)-State_VGB(1,i,j,k,iBLK))>.01) n=n+1
        end do; end do; end do
        if(n/=0) write(*,*)' PE=',iProc,' BLK=',iBLK,' has',n,'bad values ...'
!!$        do i=-1,nI+2; do j=-1,nJ+2; do k=-1,nK+2
!!$           if(abs(V0(i,j,k,iBLK)-State_VGB(1,i,j,k,iBLK))>.01)then
!!$              write(*,*)' PE=',iProc,' BLK=',iBLK,' IJK=',i,j,k,&
!!$                   ' Diff=',abs(V0(i,j,k,iBLK)-State_VGB(1,i,j,k,iBLK)), &
!!$                   ' Values=',V0(i,j,k,iBLK),State_VGB(1,i,j,k,iBLK)
!!$           end if
!!$        end do; end do; end do
     end do
     call stop_mpi("testmessage_pass_cells error")
  end if

  deallocate(V0)
  deallocate(V1)
  deallocate(V2)
  deallocate(State_VGB)

  write(*,*)'Testing message_pass_cells and message_pass_cells8, PE=',iProc, &
       '  All tests passed.'
  write(*,*)' '

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
!!$  stop

end subroutine testmessage_pass_cells
!==========================================================================
!==========================================================================
!==========================================================================
subroutine timemessage_pass_cells
  use ModProcMH
  use ModMPCells, ONLY : MessageSize
  use ModMpi
  use ModAdvance, ONLY: nVar, State_VGB
  implicit none

  integer :: i,j, iError

  real*8 :: time_this,timeSUM

  !---------------------------------------------------------------------------

  call message_pass_cells8(.false.,.true.,.true.,nVar, State_VGB)

  timeSUM=0.
  do j=1,10
     call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------  
     time_this=MPI_WTIME()  
     call message_pass_cells8(.false.,.true.,.true.,nVar, State_VGB)
     call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
     timeSUM=timeSUM+(MPI_WTIME()-time_this)
  end do
  if(iProc==0) &
       write(*,'(a,f8.5,a,i7)')' Pass STATE   faceonly  ',timeSUM/10.,' sec, MessageSize=',MessageSize

  MessageSize=1000000
  timeSUM=0.
  do j=1,10
     call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------  
     time_this=MPI_WTIME()  
     call message_pass_cells8(.false.,.true.,.true.,nVar, State_VGB)
     call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
     timeSUM=timeSUM+(MPI_WTIME()-time_this)
  end do
  if(iProc==0) &
       write(*,'(a,f8.5,a,i7)')' Pass STATE   faceonly  ',timeSUM/10.,' sec, MessageSize=',MessageSize

  MessageSize=7
  do i=1,50
     MessageSize=MessageSize+1
     if(MessageSize>1000000)EXIT

     timeSUM=0.
     do j=1,10
        call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------  
        time_this=MPI_WTIME()  
        call message_pass_cells8(.false.,.true.,.true.,nVar, State_VGB)
        call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
        timeSUM=timeSUM+(MPI_WTIME()-time_this)
     end do
     if(iProc==0) &
          write(*,'(a,f8.5,a,i7)')' Pass STATE   faceonly  ',timeSUM/10.,' sec, MessageSize=',MessageSize

  end do

end subroutine timemessage_pass_cells

!^CFG END DEBUGGING
