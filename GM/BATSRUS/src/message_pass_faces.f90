!^CFG COPYRIGHT UM
!==========================================================================
!==========================================================================
!==========================================================================
module ModMPFaces

  integer :: numSendRecv, numSend, numRecv, numCopy
  integer :: numSendRecvMax=0, numSendMax=0, numRecvMax=0, numCopyMax=0

  integer, dimension(:),   allocatable, save :: &
       nSend, nRecv, nSendStart, nRecvStart
  integer, dimension(:,:), allocatable, save :: &
       VSendI, VRecvI, VSendIlocal,VRecvIlocal

  real,    dimension(:),   allocatable, save :: VSend, VRecv
  real,    dimension(:,:), allocatable, save :: VSend9, VRecv9

  integer :: iLastGrid = -1, iLastDecomposition = -1
  integer :: itag, lS(0:7), lR(0:4), nSends
  integer, parameter :: maxMessages=10000
  integer :: nVar=99999

  logical, parameter :: DoRSend=.false.
  logical, parameter :: DoBreakUpMessages=.false.
  logical, parameter :: DoImplicitUnusedBlock=.true.

end module ModMPFaces

!==========================================================================
!==========================================================================
!==========================================================================
subroutine message_pass_faces(Vx,Vy,Vz)
  !
  ! This routine will complete a messagepass of passed in face variables
  !   Vx, Vy, and Vz.  It will update the face values for every cell with
  !   finer neighbors in the base six directions.
  !
  ! NOTE: restriction is the sum of 4 finer face values
  !
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK, nBLK, iNewGrid, iNewDecomposition
  use ModMPFaces
  use ModNumConst
  use ModMpi
  implicit none

  !Subroutine arguements
  real, intent(inout), dimension(1:nJ,1:nK,1:2,nBLK) :: Vx
  real, intent(inout), dimension(1:nI,1:nK,1:2,nBLK) :: Vy
  real, intent(inout), dimension(1:nI,1:nJ,1:2,nBLK) :: Vz

  !Local variables
  integer :: i,j,k, iV,iBLK,iPE, iError
  integer :: nSENDrequests, SENDrequests(maxMessages)
  integer :: nRECVrequests, RECVrequests(maxMessages)
  integer :: MESGstatus(MPI_STATUS_SIZE, maxMessages)

  !------------------------------------------

  ! Check that indices are up to date
  if(iNewGrid/=iLastGrid .or. iNewDecomposition/=iLastDecomposition) &
       call mp_faces_set_indices

  ! When neighbor is on the same processor, Collect/Send/Assign are all
  !    done in one step, without intermediate memory use.
  iPE=iProc
  do iV=1,nSend(iPE)
     lS(:)=VSendIlocal(:,iV)
     lR(:)=VRecvIlocal(:,iV)
     select case(-lS(0))
     case(1,2)
        Vx(lR(1),lR(2),lR(3),lR(4)) = sum(Vx(lS(1):lS(2),lS(3):lS(4),lS(5),lS(7)))
     case(3,4)
        Vy(lR(1),lR(2),lR(3),lR(4)) = sum(Vy(lS(1):lS(2),lS(3):lS(4),lS(5),lS(7)))
     case(5,6)
        Vz(lR(1),lR(2),lR(3),lR(4)) = sum(Vz(lS(1):lS(2),lS(3):lS(4),lS(5),lS(7)))
     end select
  end do

  ! Collect values into VSend that need to be passed to other processors
  do iPE=0,nProc-1
     if(iPE==iProc) CYCLE
     if(nSend(iPE)==0) CYCLE
     do iV=nSendStart(iPE)+1,nSendStart(iPE)+nSend(iPE)
        lS(:)=VSendI(:,iV)
        select case(-lS(0))
        case(1,2)
           VSend(iV) = sum(Vx(lS(1):lS(2),lS(3):lS(4),lS(5),lS(7)))
        case(3,4)
           VSend(iV) = sum(Vy(lS(1):lS(2),lS(3):lS(4),lS(5),lS(7)))
        case(5,6)
           VSend(iV) = sum(Vz(lS(1):lS(2),lS(3):lS(4),lS(5),lS(7)))
        end select
     end do
  end do

  ! Post receives first so that they are ready
  nRECVrequests = 0  
  do iPE=0,nProc-1
     if(iPE == iProc) CYCLE
     if(nRecv(iPE)==0) CYCLE
     if(DoBreakUpMessages)then
        nSends=1+((nRecv(iPE)-1)/nVar)
        do i=1,nSends
           itag = nProc*(i-1)+iPE
           nRECVrequests = nRECVrequests + 1
           if(nRECVrequests>maxMessages) call stop_mpi("Too many RECVs in mp_SendValues")
           call MPI_irecv(VRecv(nRecvStart(iPE)+1+((i-1)*nVar)),min(nVar,nRecv(iPE)-((i-1)*nVar)), &
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
        nSends=1+((nSend(iPE)-1)/nVar)
        do i=1,nSends
           itag = nProc*(i-1)+iProc
           if(DoRSend)then
              call MPI_rsend(VSend(nSendStart(iPE)+1+((i-1)*nVar)),min(nVar,nSend(iPE)-((i-1)*nVar)), &
                   MPI_REAL,iPE,itag,iComm,iError)
           else
              nSENDrequests = nSENDrequests + 1
              call MPI_isend(VSend(nSendStart(iPE)+1+((i-1)*nVar)),min(nVar,nSend(iPE)-((i-1)*nVar)), &
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
        select case(-lR(0))
        case(1,2)
           Vx(lR(1),lR(2),lR(3),lR(4)) = VRecv(iV)
        case(3,4)
           Vy(lR(1),lR(2),lR(3),lR(4)) = VRecv(iV)
        case(5,6)
           Vz(lR(1),lR(2),lR(3),lR(4)) = VRecv(iV)
        end select
     end do
  end do

  ! Wait for sent messages to be received before exiting
  if(.not.DoRSend)then
     call MPI_waitall(nSENDrequests, SENDrequests(1), MESGstatus(1,1), iError)
  end if

end subroutine message_pass_faces

!==========================================================================
!==========================================================================
!==========================================================================
subroutine message_pass_faces_9conserve
  !
  ! This routine will complete a messagepass of conserved face variables
  !   rho, rhoUx, rhoUy, rhoUz,Bx, By, Bz, E, and Vdt.  It will update the
  !   face values for every cell with finer neighbors in the base six
  !   directions.
  !
  ! NOTE: restriction is the sum of 4 finer face values
  !
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK, nBLK, iNewGrid, iNewDecomposition
  use ModConserveFlux, ONLY : nCorrectedFaceValues,&
       CorrectedFlux_VXB, CorrectedFlux_VYB, CorrectedFlux_VZB
  use ModMPFaces
  use ModNumConst
  use ModMpi
  implicit none

  !Local variables
  integer :: i,j,k, iV,iBLK,iPE, iError, iVar
  integer :: nSENDrequests, SENDrequests(maxMessages)
  integer :: nRECVrequests, RECVrequests(maxMessages)
  integer :: MESGstatus(MPI_STATUS_SIZE, maxMessages)

  !------------------------------------------

  ! Check that indices are up to date
  if(iNewGrid/=iLastGrid .or. iNewDecomposition/=iLastDecomposition) &
       call mp_faces_set_indices

  ! When neighbor is on the same processor, Collect/Send/Assign are all
  !    done in one step, without intermediate memory use.
  iPE=iProc
  do iV=1,nSend(iPE)
     lS(:)=VSendIlocal(:,iV)
     lR(:)=VRecvIlocal(:,iV)
     select case(-lS(0))
     case(1,2)
        do iVar=1,nCorrectedFaceValues
           CorrectedFlux_VXB(iVar,lR(1),lR(2),lR(3),lR(4)) = &
                sum( CorrectedFlux_VXB(iVar,lS(1):lS(2),lS(3):lS(4),lS(5),lS(7)))
        end do
     case(3,4)
        do iVar=1,nCorrectedFaceValues
           CorrectedFlux_VYB(iVar,lR(1),lR(2),lR(3),lR(4)) = &
                sum( CorrectedFlux_VYB(iVar,lS(1):lS(2),lS(3):lS(4),lS(5),lS(7)))
        end do
     case(5,6)
        do iVar=1,nCorrectedFaceValues
           CorrectedFlux_VZB(iVar,lR(1),lR(2),lR(3),lR(4)) = &
                sum( CorrectedFlux_VZB(iVar,lS(1):lS(2),lS(3):lS(4),lS(5),lS(7)))
        end do
     end select
  end do

  ! Collect values into VSend9 that need to be passed to other processors
  do iPE=0,nProc-1
     if(iPE==iProc) CYCLE
     if(nSend(iPE)==0) CYCLE
     do iV=nSendStart(iPE)+1,nSendStart(iPE)+nSend(iPE)
        lS(:)=VSendI(:,iV)
        select case(-lS(0))
        case(1,2)
           do iVar=1,nCorrectedFaceValues
               VSend9(iVar,iV) = &
                   sum( CorrectedFlux_VXB(iVar,lS(1):lS(2),lS(3):lS(4),lS(5),lS(7)))
           end do
        case(3,4)
           do iVar=1,nCorrectedFaceValues
               VSend9(iVar,iV) = &
                   sum( CorrectedFlux_VYB(iVar,lS(1):lS(2),lS(3):lS(4),lS(5),lS(7)))
           end do
        case(5,6)
           do iVar=1,nCorrectedFaceValues
               VSend9(iVar,iV) = &
                   sum( CorrectedFlux_VZB(iVar,lS(1):lS(2),lS(3):lS(4),lS(5),lS(7)))
           end do
        end select
     end do
  end do

  ! Post receives first so that they are ready
  nRECVrequests = 0  
  do iPE=0,nProc-1
     if(iPE == iProc) CYCLE
     if(nRecv(iPE)==0) CYCLE
     if(DoBreakUpMessages)then
        nSends=1+((nRecv(iPE)-1)/nVar)
        do i=1,nSends
           itag = nProc*(i-1)+iPE
           nRECVrequests = nRECVrequests + 1
           if(nRECVrequests>maxMessages) call stop_mpi("Too many RECVs in mp_SendValues")
           call MPI_irecv(VRecv9(1,nRecvStart(iPE)+1+((i-1)*nVar)),&
                nCorrectedFaceValues*min(nVar,nRecv(iPE)-((i-1)*nVar)), &
                MPI_REAL,iPE,itag,iComm,RECVrequests(nRECVrequests),iError)
        end do
     else
        itag = iPE
        nRECVrequests = nRECVrequests + 1
        if(nRECVrequests>maxMessages) call stop_mpi("Too many RECVs in mp_SendValues")
        call MPI_irecv(VRecv9(1,nRecvStart(iPE)+1),&
             nCorrectedFaceValues*nRecv(iPE), &
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
        nSends=1+((nSend(iPE)-1)/nVar)
        do i=1,nSends
           itag = nProc*(i-1)+iProc
           if(DoRSend)then
              call MPI_rsend(VSend9(1,nSendStart(iPE)+1+((i-1)*nVar)),&
                   nCorrectedFaceValues*min(nVar,nSend(iPE)-((i-1)*nVar)), &
                   MPI_REAL,iPE,itag,iComm,iError)
           else
              nSENDrequests = nSENDrequests + 1
              call MPI_isend(VSend9(1,nSendStart(iPE)+1+((i-1)*nVar)),&
                   nCorrectedFaceValues*min(nVar,nSend(iPE)-((i-1)*nVar)), &
                   MPI_REAL,iPE,itag,iComm,SENDrequests(nSENDrequests),iError)
           end if
        end do
     else
        itag = iProc
        if(DoRSend)then
           call MPI_rsend(VSend9(1,nSendStart(iPE)+1),&
                nCorrectedFaceValues*nSend(iPE), &
                MPI_REAL,iPE,itag,iComm,iError)
        else
           nSENDrequests = nSENDrequests + 1
           call MPI_isend(VSend9(1,nSendStart(iPE)+1),&
                nCorrectedFaceValues*nSend(iPE), &
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
        select case(-lR(0))
        case(1,2)
           do iVar=1,nCorrectedFaceValues
                 CorrectedFlux_VXB(iVar,lR(1),lR(2),lR(3),lR(4))= VRecv9(iVar,iV) 
           end do
        case(3,4)
           do iVar=1,nCorrectedFaceValues
                 CorrectedFlux_VYB(iVar,lR(1),lR(2),lR(3),lR(4))= VRecv9(iVar,iV) 
           end do
        case(5,6)
           do iVar=1,nCorrectedFaceValues
                 CorrectedFlux_VZB(iVar,lR(1),lR(2),lR(3),lR(4))= VRecv9(iVar,iV) 
           end do
        end select
     end do
  end do

  ! Wait for sent messages to be received before exiting
  if(.not.DoRSend)then
     call MPI_waitall(nSENDrequests, SENDrequests(1), MESGstatus(1,1), iError)
  end if

end subroutine message_pass_faces_9conserve

!==========================================================================
!==========================================================================
!==========================================================================
subroutine mp_faces_set_indices
  !
  !
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK, nBLK, iNewGrid, iNewDecomposition
  use ModConserveFlux, ONLY: nCorrectedFaceValues
  use ModMPFaces
  use ModMpi
  implicit none

  !Local variables
  integer :: iPE, iError
  !------------------------------------------

  iLastGrid          = iNewGrid
  iLastDecomposition = iNewDecomposition

  call mp_allocate_face_arrays1
  call mp_build_face_indices(.true.)

  call mp_allocate_face_arrays2
  call mp_build_face_indices(.false.)

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------


contains

  !==========================================================================
  subroutine mp_allocate_face_arrays1
    !
    ! Allocate initial arrays for counting the number of sends and receives
    !

    !------------------------------------------

    ! Allocate new memory only if not already allocated, then initialize

    if(.not.allocated(nSend))then
       allocate( nSend(0:nProc-1), stat=iError ); call alloc_check(iError,"nSend")
    end if
    nSend=0

    if(.not.allocated(nRecv))then
       allocate( nRecv(0:nProc-1), stat=iError ); call alloc_check(iError,"nRecv")
    end if
    nRecv=0

    if(.not.allocated(nSendStart))then
       allocate( nSendStart(0:nProc-1), stat=iError ); call alloc_check(iError,"nSendStart")
    end if
    nSendStart=0

    if(.not.allocated(nRecvStart))then
       allocate( nRecvStart(0:nProc-1), stat=iError ); call alloc_check(iError,"nRecvStart")
    end if
    nRecvStart=0

  end subroutine mp_allocate_face_arrays1

  !==========================================================================
  subroutine mp_allocate_face_arrays2
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

!!$    write(*,*)'MPF: PE=',iProc,' SEND:',nSend,' RECV:',nRecv,' SUMMARY:',numSend,numRecv,numCopy

    ! Re-Allocate memory if past use was not large enough

    if(  numSendRecv > numSendRecvMax .or. &
         numSend     > numSendMax     .or. &
         numRecv     > numRecvMax     .or. &
         numCopy     > numCopyMax           )then

       ! Deallocate old memory and allocate new memory

       if(allocated(VSendI)) deallocate(VSendI)
       allocate( VSendI(0:7,numSendRecv), stat=iError ); call alloc_check(iError,"VSendI")

       if(allocated(VRecvI)) deallocate(VRecvI)
       allocate( VRecvI(0:4,numSendRecv), stat=iError ); call alloc_check(iError,"VRecvI")

       if(allocated(VSendIlocal)) deallocate(VSendIlocal)
       allocate( VSendIlocal(0:7,numCopy), stat=iError ); call alloc_check(iError,"VSendIlocal")

       if(allocated(VRecvIlocal)) deallocate(VRecvIlocal)
       allocate( VRecvIlocal(0:4,numCopy), stat=iError ); call alloc_check(iError,"VRecvIlocal")

       if(allocated(VSend )) deallocate(VSend )
       allocate( VSend(numSend), stat=iError ); call alloc_check(iError,"VSend")

       if(allocated(VRecv )) deallocate(VRecv )
       allocate( VRecv(numRecv), stat=iError ); call alloc_check(iError,"VRecv")

       if(allocated(VSend9)) deallocate(VSend9)
       allocate( VSend9(nCorrectedFaceValues,numSend), stat=iError ); call alloc_check(iError,"VSend9")

       if(allocated(VRecv9)) deallocate(VRecv9)
       allocate( VRecv9(nCorrectedFaceValues,numRecv), stat=iError ); call alloc_check(iError,"VRecv9")

       ! Update new max values

       numCopyMax = numCopy
       numSendMax = numSend
       numRecvMax = numRecv
       numSendRecvMax = numSendRecv

    end if

    ! Initialize memory
    VSendI=-9
    VRecvI=-9
    VSendIlocal=-9
    VRecvIlocal=-9
    VSend=0
    VRecv=0
    VSend9=0
    VRecv9=0

  end subroutine mp_allocate_face_arrays2

end subroutine mp_faces_set_indices

!==========================================================================
!==========================================================================
!==========================================================================
subroutine mp_build_face_indices(JustCount)
  !
  !
  use ModProcMH
  use ModMain
  use ModOctree
  use ModMPFaces
  use ModAMR, ONLY : unusedBlock_BP
  implicit none

  !Subroutine arguements
  logical, intent(in) :: JustCount

  !Local variables
  integer :: iBLK,iPE, iCHILD, idir
  

  integer, dimension(6,3) :: dLOOP

  integer :: neighborLEV
  integer, dimension(4) :: neighborPE,neighborBLK,neighborCHILD
  !------------------------------------------

  ! face=1-6
  !    1   2   3   4   5   6  
  !    W   E   N   S   T   B   

  data dLOOP / &
       1, -1,  0,  0,  0,  0, &
       0,  0,  1, -1,  0,  0, &
       0,  0,  0,  0,  1, -1  /


  nSend=0
  nRecv=0

  do iPE = 0,nProc-1
     do iBLK = 1,nBLK
        if (.not.associated(global_block_ptrs(iBLK, iPE+1) % ptr)) CYCLE

        if (.not.global_block_ptrs(iBLK, iPE+1) % ptr % used) CYCLE

        if (DoImplicitUnusedBlock)then
           if (unusedBlock_BP(iBLK,iPE)) CYCLE
        end if

        !valid used block found, setup indices
        iCHILD = global_block_ptrs(iBLK, iPE+1) % ptr % child_number

        do idir=1,6
           call treeNeighbor(iPE,iBLK,dLOOP(idir,1),dLOOP(idir,2),dLOOP(idir,3), &
                neighborPE,neighborBLK,neighborCHILD,neighborLEV)

           if(neighborLEV==1)call build_i
              !Build indices for send to coarser block level
        end do
     end do
  end do


contains

  !==========================================================================
  subroutine build_i
    !
    !
    Use ModCube
    !Local variables
    integer :: i,j,k, n
    integer,dimension(3)::iMinS_D,iMaxS_D,iMinR_D,iMaxR_D
    integer :: i1S,i2S, j1S,j2S, k1S,k2S, i1R,i2R, j1R,j2R, k1R,k2R
    integer :: nborPE, nborBLK

    !------------------------------------------

    if(idir>6)then
       write(*,*)'ERROR: message_pass_faces can only work with directions 1-6 (face direction).'
       write(*,*)'       PE=',iProc,'  idir=',idir
       call stop_mpi('bad idir value in message_pass_faces')
    end if

    nborPE  = neighborPE (1)
    nborBLK = neighborBLK(1)

    if (DoImplicitUnusedBlock)then
       if (unusedBlock_BP(nborBLK,nborPE)) return
    end if

    call set_indices(&
                   iDirS2R_D=dLoop(iDir,:), &
                   nLayerS=1,nLayerR=1,&
                   iMinS_D=iMinS_D,iMaxS_D=iMaxS_D,&
                   iMinR_D=iMinR_D,iMaxR_D=iMaxR_D,&
                   iLevelR=1,iChild=iChild)

    
    select case(iDir)
    case(1)
       i1S=iMinS_D(2); j1S=iMinS_D(3)
       i2S=iMaxS_D(2); j2S=iMaxS_D(3)      
       i1R=iMinR_D(2); j1R=iMinR_D(3)
       i2R=iMaxR_D(2); j2R=iMaxR_D(3)       
       k1S=2         ; k1R=1
    case(2)
       i1S=iMinS_D(2); j1S=iMinS_D(3)
       i2S=iMaxS_D(2); j2S=iMaxS_D(3)      
       i1R=iMinR_D(2); j1R=iMinR_D(3)
       i2R=iMaxR_D(2); j2R=iMaxR_D(3)       
       k1S=1         ; k1R=2
    case(3)
       i1S=iMinS_D(1); j1S=iMinS_D(3)
       i2S=iMaxS_D(1); j2S=iMaxS_D(3)      
       i1R=iMinR_D(1); j1R=iMinR_D(3)
       i2R=iMaxR_D(1); j2R=iMaxR_D(3)       
       k1S=2         ; k1R=1
    case(4)
       i1S=iMinS_D(1); j1S=iMinS_D(3)
       i2S=iMaxS_D(1); j2S=iMaxS_D(3)      
       i1R=iMinR_D(1); j1R=iMinR_D(3)
       i2R=iMaxR_D(1); j2R=iMaxR_D(3)       
       k1S=1         ; k1R=2
    case(5)
       i1S=iMinS_D(1); j1S=iMinS_D(2)
       i2S=iMaxS_D(1); j2S=iMaxS_D(2)      
       i1R=iMinR_D(1); j1R=iMinR_D(2)
       i2R=iMaxR_D(1); j2R=iMaxR_D(2)       
       k1S=2         ; k1R=1
    case(6)
       i1S=iMinS_D(1); j1S=iMinS_D(2)
       i2S=iMaxS_D(1); j2S=iMaxS_D(2)      
       i1R=iMinR_D(1); j1R=iMinR_D(2)
       i2R=iMaxR_D(1); j2R=iMaxR_D(2)       
       k1S=1         ; k1R=2
    end select
    k2S=k1S  ;  k2R=k1R
       if(iProc == iPE)then
          do i=i1S,i2S,2; do j=j1S,j2S,2; do k=k1S,k2S
             nSend(nborPE)=nSend(nborPE)+1
             if(JustCount) CYCLE

             if(iPE==nborPE)then
                n=nSend(nborPE)
                VSendIlocal(0,n)=-idir
                VSendIlocal(1,n)=i   ; VSendIlocal(3,n)=j   ; VSendIlocal(5,n)=k
                VSendIlocal(2,n)=i+1 ; VSendIlocal(4,n)=j+1 ; VSendIlocal(6,n)=k
                VSendIlocal(7,n)=iBLK
             else
                n=nSendStart(nborPE)+nSend(nborPE)
                VSendI(0,n)=-idir
                VSendI(1,n)=i   ; VSendI(3,n)=j   ; VSendI(5,n)=k
                VSendI(2,n)=i+1 ; VSendI(4,n)=j+1 ; VSendI(6,n)=k
                VSendI(7,n)=iBLK
             end if
          end do; end do; end do
       end if
       if(iProc == nborPE)then
          do i=i1R,i2R; do j=j1R,j2R; do k=k1R,k2R
             nRecv(iPE)=nRecv(iPE)+1
             if(JustCount) CYCLE

             if(iPE==nborPE)then
                n=nRecv(iPE)
                VRecvIlocal(0,n)=-idir
                VRecvIlocal(1,n)=i; VRecvIlocal(2,n)=j; VRecvIlocal(3,n)=k
                VRecvIlocal(4,n)=nborBLK
             else
                n=nRecvStart(iPE)+nRecv(iPE)
                VRecvI(0,n)=-idir
                VRecvI(1,n)=i; VRecvI(2,n)=j; VRecvI(3,n)=k
                VRecvI(4,n)=nborBLK
             end if
          end do; end do; end do
       end if


  end subroutine build_i
end subroutine mp_build_face_indices

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

