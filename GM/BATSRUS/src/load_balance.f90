!^CFG COPYRIGHT UM
subroutine load_balance(DoMoveCoord, DoMoveData, IsNewBlock)
  use ModProcMH
  use ModMain
  use ModImplicit, ONLY : UsePartImplicit !^CFG IF IMPLICIT
  use ModAdvance, ONLY: iTypeAdvance_B, iTypeAdvance_BP,&
       SkippedBlock_, SteadyBlock_, SteadyBoundBlock_, ExplBlock_, ImplBlock_
  use ModGeometry, ONLY: True_Blk
  use ModGeometry, ONLY: UseCovariant     
  use ModPartSteady, ONLY: UsePartSteady
  use ModAMR, ONLY : availableBLKs, UnusedBlock_BP
  use ModParallel
  use ModIO, ONLY: write_prefix, iUnitOut
  use ModMpi
  use ModEnergy, ONLY: calc_energy_ghost
  use ModConserveFlux, ONLY: init_cons_flux
  implicit none

  ! Load balance grid using Peano-Hilbert ordering of blocks
  ! Coordinates are moved if DoMoveCoord is true.
  ! Data is moved with the blocks if DoMoveData is true.
  ! There are new blocks (due to initial refinement, restart or AMR)
  ! when IsNewBlock is new (so update neighbors etc).

  logical, intent(in) :: DoMoveCoord, DoMoveData, IsNewBlock

  ! Maximum number of attempts to accomplish the load balancing
  ! The algorithm needs multiple tries if the actual number of blocks
  ! is very close to the maximum number of blocks and many blocks are moved
  integer, parameter :: MaxTry = 100
  
  ! Set this logical to .false. to return to the previous version,
  ! which did not move the B0 and body force variables.
  ! There is another declaration in subroutine move_block! Change together!!!
  logical, parameter :: DoMoveExtraData = .true.

  ! Maximum number of block types to be load balanced separately
  integer, parameter :: MaxType = 10

  ! Actual number of block types
  integer :: nType

  ! Index for block type
  integer :: iType

  ! Global block index for the various block types
  integer :: iBlockALL_I(MaxType)

  ! Number of blocks for each type
  integer :: nBlockALL_I(MaxType)

  ! Conversion from iTypeAdvance (including body block info) to iType
  integer :: iType_I(-ImplBlock_:ImplBlock_)

  ! load balance distribute each type
  integer :: nBlock_PI(0:nProc-1,MaxType), iProcTo_I(MaxType)
  integer :: iProcStart, iProcStop, iProcExtraBlock

  ! Number of blocks moved around
  integer :: nBlockMoved

  integer :: iError
  integer :: iBlockALL, iBlock

  integer :: iBlockFrom, iProcFrom, iBlockTo, iProcTo, iTry

  logical :: SkippedAnyBlock

  logical :: DoTest, DoTestMe

  logical :: DoFixVar_B(MaxBlock)

  !---------------------------------------------------------------------------
  call set_oktest('load_balance',DoTest,DoTestMe)

  if (DoTestMe) write(*,*)'load_balance: ',&
       'iProc, DoMoveCoord, DoMoveData, IsNewBlock=',&
       iProc, DoMoveCoord, DoMoveData, IsNewBlock

  ! starting value of number of blocks moved between processors
  nBlockMoved = 0

  ! Find the last used block on the processor
  do iBlock = nBlockMax,1,-1
     nBlock = iBlock
     if (.not.unusedBLK(iBlock)) EXIT
  end do

  if(DoMoveData .and. .not.DoMoveCoord)call stop_mpi(&
       'ERROR in load_balance: DoMoveData=T and DoMoveCoord=F !!!')

  !^CFG IF IMPLICIT BEGIN
  !\
  ! Select blocks for implicit/local time stepping
  !/
  ! Find number of implicit and explicit blocks.   
  ! Partial selection is only possible if dt_BLK and coords are known
  ! i.e. if DoMoveCoord is true.
  call select_stepping(DoMoveCoord)
  !^CFG END IMPLICIT

  if (nProc==1 .or. index(test_string,'NOLOADBALANCE')>0) then
     call finish_load_balance
     RETURN
  end if

  ! If the coordinates are known then include the body block info into
  ! iTypeAdvance_B and _BP by changing the sign to negative for body blocks
  if(DoMoveCoord)then
!!! If there was a IsTrueBlock_BP array there was no need for MPI_ALLGATHER
     where(.not. True_BLK(1:nBlock)) &
          iTypeAdvance_B(1:nBlock) = -abs(iTypeAdvance_B(1:nBlock))

     ! Update iTypeAdvance_BP
     call MPI_ALLGATHER(iTypeAdvance_B, MaxBlock, MPI_INTEGER, &
          iTypeAdvance_BP, MaxBlock, MPI_INTEGER, iComm, iError)
  end if

  ! Set the transformation from iTypeAdvance to iType
  iType_I = 1
  iType_I(SkippedBlock_) = 0
  if(UsePartSteady)then
     iType_I(SteadyBoundBlock_) = 2
     iType_I(ExplBlock_)        = 2

     iType_I(-SteadyBoundBlock_) = 3
     iType_I(-ExplBlock_)        = 3
  elseif(UsePartImplicit)then                !^CFG IF IMPLICIT
     iType_I( ImplBlock_) = 2                !^CFG IF IMPLICIT
     iType_I(-ImplBlock_) = 2                !^CFG IF IMPLICIT
  else
     iType_I(-ExplBlock_) = 2
  endif
  nType = maxval(iType_I)

  ! Count the number of blocks of various types
  nBlockALL_I = 0
  do iProcTo = 0, nProc-1; do iBlock = 1, nBlockMax
     iType = iType_I(iTypeAdvance_BP(iBlock, iProcTo))
     if(iType > 0) nBlockALL_I(iType) = nBlockALL_I(iType) + 1
  end do; end do

  ! Construct load balance table for various types
  do iType=1,nType
     ! minimum load for each processor
     nBlock_PI(0:nProc-1,iType) = nBlockALL_I(iType)/nProc

     ! The processors with extra blocks are filled in from nProc-1 backwards
     iProcStart = nProc - modulo(sum(nBlockALL_I(1:iType)),nProc)
     iProcStop = iProcStart + modulo(nBlockALL_I(iType),nProc) - 1
     do iProcExtraBlock=iProcStart,iProcStop
        iProcTo = modulo(iProcExtraBlock,nProc)
        nBlock_PI(iProcTo,iType) = nBlock_PI(iProcTo,iType) + 1
     end do

     ! convert to accumilative load table
     do iProcTo=1,nProc-1
        nBlock_PI(iProcTo,iType) = nBlock_PI(iProcTo,iType) &
             + nBlock_PI(iProcTo-1,iType)
     end do
  end do

  if(DoTestMe)then
     write(*,'(a,i6)') &
          'load_balance starting: nBlockMax=',nBlockMax
     write(*,'(a,i4,2i6)') &
          'load_balance starting: me, nBlock, nBlockUsed=',&
          iProc, nBlock, count(.not.unusedBLK(1:nBlock))
     if(nType > 1) then
        write(*,'(a,3i4)') &
             'load_balance starting: iProc, max, min(iTypeAdvance_BP)=',&
             iProc,&
             maxval(iTypeAdvance_BP(1:nBlockMax,:)),&
             minval(iTypeAdvance_BP(1:nBlockMax,:))
        write(*,'(a,i4,10i7)')'load_balance starting: iProc, nBlockALL_I=',&
             iProc, nBlockALL_I(1:nType)
        do iType = 1, nType
           write(*,'(a,i4,i2,i6)')'load_balance starting: '// &
                'iProc, iType, count=',&
                iProc, iType,&
                count(iType_I(iTypeAdvance_B(1:nBlockMax))==iType)
        end do
     end if
  end if

  call timing_start('load_balance')

  if(DoMoveData) DoFixVar_B = .false. ! initialize variable fixing flags

  TRY: do iTry=1,MaxTry

     skippedAnyBlock=.false.
     iBlockALL_I = 0
     iProcTo_I = 0

     TREE: do iBlockALL  = 1, nBlockALL

        iBlockFrom = iBlock_A(iBlockALL)
        iProcFrom  = iProc_A(iBlockALL)

        ! Figure out the block type
        iType = iType_I(iTypeAdvance_BP(iBlockFrom,iProcFrom))

        ! Increase the index for this block type and select target processor
        iBlockALL_I(iType) = iBlockALL_I(iType) + 1

        PROC: do iProcTo=iProcTo_I(iType),nProc-1
           if(iBlockALL_I(iType) <= nBlock_PI(iProcTo,iType))then
              iProcTo_I(iType) = iProcTo
              EXIT PROC
           end if
        end do PROC
        if(iProcTo>=nProc) call stop_mpi("load_balance error: iProcTo==nProc")
        if(iProcTo/=iProcTo_I(iType)) call stop_mpi("load_balance error: F90 failed")

        !DEBUG
        !if(iProc==0)write(*,*)'iBlockALL,iBlockFrom,iProcFrom,iProcTo=',&
        !     iBlockALL,iBlockFrom,iProcFrom,iProcTo

        if(iProcTo == iProcFrom) CYCLE TREE

        ! Check if ProcTo would get too many blocks
        if(availableBLKs(0,iProcTo) > MaxBlock) then
           skippedAnyBlock = .true.
           CYCLE TREE
        end if

        ! Select next available block
        iBlockTo = availableBLKs(availableBLKs(0,iProcTo),iProcTo)
        availableBLKs(0,iProcTo)   = availableBLKs(0,iProcTo)+1
        availableBLKs(0,iProcFrom) = availableBLKs(0,iProcFrom)-1
        availableBLKs(availableBLKs(0,iProcFrom),iProcFrom) = iBlockFrom

        ! Move the block from ProcFrom to Procto
        call move_block(DoMoveCoord, DoMoveData, iBlockALL,&
             iBlockFrom, iProcFrom, iBlockTo, iProcTo)

        nBlockMoved=nBlocKMoved+1

        if(DoMoveData .and. iProc==iProcTo) DoFixVar_B(iBlockTo)=.true.

        iBlock_A(iBlockALL) = iBlockTo
        iProc_A(iBlockALL)  = iProcTo

     end do TREE
     if(.not.skippedAnyBlock) EXIT TRY
  end do TRY

  call MPI_ALLREDUCE(nBlock, nBlockMax, 1, MPI_INTEGER, MPI_MAX, &
       iComm, iError)

  ! Change decomposition ID if any blocks were moved
  if(nBlockMoved > 0)iNewDecomposition = mod(iNewDecomposition+1,10000)

  ! Fix variables if any data was moved
  if(DoMoveData .and. nBlockMoved > 0)then
     !call timing_start('load_fix_var')
     do iBlock = 1,nBlock
        if(.not.DoFixVar_B(iBlock)) CYCLE
        if(useConstrainB) call Bface2Bcenter(iBlock)  !^CFG IF CONSTRAINB
        call calc_energy_ghost(iBlock)

        if(DoMoveExtraData)then
           if(UseB0)then
              if(UseCovariant)then                      
                 call calc_b0source_covar(iBlock)  
              else                                       
                 call set_b0_matrix(iBlock)              
              end if
           end if
           call init_cons_flux(iBlock)
        else
           call calc_other_soln_vars(iBlock)
        end if
     end do
     !call timing_stop('load_fix_var')
  end if

  if(DoTestMe)then
     write(*,'(a,i4,2i6)') &
          'load_balance finished: nTry, nBlockMax, nBlockMoved=',&
          iTry, nBlockMax, nBlockMoved

     write(*,'(a,i4,3i6)')&
       'load_balance finished: me, nBlock, nBlockUsed=',&
       iProc, nBlock, count(.not.unusedBLK(1:nBlock)), &
       count(iTypeAdvance_B /= SkippedBlock_)

     if(nType > 1) then
        do iType = 1, nType
           write(*,'(a,i4,i2,i6)') &
                'load_balance finished: iProc, iType, count=',&
                iProc, iType,count(iType_I(iTypeAdvance_B(1:nBlockMax))==iType)
        end do
     end if
  end if

  ! restore iTypeAdvance_B and _BP to be positive
  iTypeAdvance_B  = abs(iTypeAdvance_B)
  iTypeAdvance_BP = abs(iTypeAdvance_BP)

  call finish_load_balance

  call timing_stop('load_balance')

contains
  !===========================================================================
  subroutine finish_load_balance

    use ModCovariant,ONLY: UseVertexBasedGrid, &       
         do_fix_geometry_at_reschange                  

    integer :: iBlock

    ! Update neighbor and test cell info if there are new and/or moved blocks
    if(IsNewBlock .or. nBlockMoved > 0)then

       call find_neighbors

       if(DoMoveCoord .and. UseVertexBasedGrid)then
          do iBlock=1, nBlock
             if(do_fix_geometry_at_reschange(iBlock)) &       
                  call fix_geometry_at_reschange(iBlock) 
          end do
       end if

       ! Analyze only when new blocks are created   !^CFG IF DEBUGGING
       if(IsNewBlock)call analyze_neighbors         !^CFG IF DEBUGGING 

       ! If coordinates are known, find the test cell
       if(DoMoveCoord)call find_test_cell

       ! When load balancing is done Skipped and Unused blocks coincide
       UnusedBlock_BP = iTypeAdvance_BP == SkippedBlock_

       ! Report number of moved blocks due to redistribution,
       ! i.e. when there are no new blocks due to AMR or initial startup.
       if(.not.IsNewBlock .and. iProc==0 .and. lVerbose > 0)then
          call write_prefix; write(iUnitOut,*)&
               'load_balance finished: nBlockMoved=',nBlockMoved
       end if

    end if

  end subroutine finish_load_balance

end subroutine load_balance

!=============================================================================

subroutine move_block(DoMoveCoord, DoMoveData, iBlockALL, &
     iBlockFrom, iProcFrom, iBlockTo,iProcTo)

  ! Move block with global block number iBlockALL 
  ! from block iBlockFrom from processor iProcFrom 
  ! to   block iBlockTo   on   processor iProcTo.
  ! Move flow variables if DoMoveData is true.

  use ModProcMH
  use ModMain
  use ModVarIndexes
  use ModAdvance, ONLY : State_VGB, &
       fbody_x_BLK, fbody_y_BLK, fbody_z_BLK, &
       B0_DGB, B0ResChange_DXSB,B0ResChange_DYSB,B0ResChange_DZSB, &
       iTypeAdvance_B, iTypeAdvance_BP, SkippedBlock_
  use ModGeometry, ONLY : dx_BLK,dy_BLK,dz_BLK,xyzStart_BLK,UseCovariant
  use ModParallel
  use ModBlockData, ONLY: get_block_data, put_block_data, &
       n_block_data, use_block_data, set_block_data, clean_block_data
  use ModImplicit                                         !^CFG IF IMPLICIT
  use ModCT, ONLY : Bxface_BLK,Byface_BLK,Bzface_BLK      !^CFG IF CONSTRAINB
  use ModRaytrace, ONLY : ray                             !^CFG IF RCM
  use ModIo,       ONLY : log_vars                        !^CFG IF RCM
  use ModMpi
  implicit none

  logical, intent(in) :: DoMoveCoord, DoMoveData
  integer, intent(in) :: iBlockALL, iBlockFrom, iProcFrom, iBlockTo,iProcTo

  ! Set this logical to .false. to return to the previous version,
  ! which did not move the B0 and body force variables.
  ! There is another declaration in subroutine load_balance! Change together!!!
  logical, parameter :: DoMoveExtraData = .true.

  integer, parameter :: nScalarBLK=14, &
       nCellGhostBLK=(nI+4)*(nJ+4)*(nK+4)
  integer, parameter :: nExtraData = &
       3*nCellGhostBLK +                 & ! B0*Cell
       6*(nJ*nK                          & ! B0*Face_x
       +  nI*nK                          & ! B0*Face_y
       +  nI*nJ)                +        & ! B0*Face_z
       3*nIJK                              ! fbody_
       
  integer, parameter :: nDataBLK= &
       nVar*nIJK +                       & !^CFG IF IMPLICIT
       3*2*nIJK +                        & !^CFG IF RCM
       nScalarBLK +                      & ! scalars
       nVar*nCellGhostBLK +              & ! State_VGB
       nExtraData +                      & ! B0, fbody
       3*nVar*nIJK                         ! max for data in ModBlockData

  ! Buffer for send and recieve
  real, dimension(nDataBLK) :: BlockData_I

  integer :: iData, itag, i, j, k, i1,i2, iVar, iw, iSize,iDim
  integer :: iError
  integer :: Status_I(MPI_STATUS_SIZE)

  ! Number of data values stored in ModBlockData
  integer :: nDynamicData
  
  logical :: DoTest, DoTestMe, DoSendRay
  !----------------------------------------------------------------------------

  if((iProcFrom==PROCtest .and. iBlockFrom==BLKtest).or. &
       (iProcTo  ==PROCtest .and. iBlockTo  ==BLKtest))then
     call set_oktest('move_block',DoTest,DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  end if
  if(DoTestMe)write(*,*)'iBlockALL,iBlockFrom,iProcFrom, iBlockTo,iProcTo=',&
       iBlockALL,iBlockFrom,iProcFrom, iBlockTo,iProcTo

  DoSendRay = UseIM .or. index(log_vars, 'status') > 0     !^CFG IF RCM

  if (iProc == iProcFrom) then
     if (DoMoveCoord) call send_block_data
     unusedBLK(iBlockFrom)   = .true.
     iTypeAdvance_B(iBlockFrom) = SkippedBlock_
     call clean_block_data(iBlockFrom)
     do
        if (nBlock==0) EXIT
        if (.not.unusedBLK(nBlock)) EXIT
        nBlock = nBlock-1
     end do
  end if

  if (iProc == iProcTo) then
     unusedBLK(iBlockTo) = .false.
     iTypeAdvance_B(iBlockTo) = iTypeAdvance_BP(iBlockFrom,iProcFrom)

     if (DoMoveCoord) call recv_block_data
     nBlock = max(nBlock, iBlockTo)
     global_block_number(iBlockTo) = iBlockALL
  end if

  ! Update global advance type info
  iTypeAdvance_BP(iBlockTo, iProcTo) = iTypeAdvance_BP(iBlockFrom,iProcFrom)
  iTypeAdvance_BP(iBlockFrom,iProcFrom) = SkippedBlock_

  ! Fix pointers
  call move_octree_block(iBlockFrom,iProcFrom,iBlockTo,iProcTo)

contains
  !============================================================================
  subroutine send_block_data
    globalBLK = iBlockFrom
    BlockData_I(1)   = dx_BLK(iBlockFrom)
    BlockData_I(2)   = dy_BLK(iBlockFrom)
    BlockData_I(3)   = dz_BLK(iBlockFrom)
    BlockData_I(4:6) = xyzStart_BLK(:,iBlockFrom)
    BlockData_I(7)   = dt_BLK(iBlockFrom)
    BlockData_I(8:13)= real(neiLev(:,iBlockFrom))

    nDynamicData = 0
    if(use_block_data(iBlockFrom)) nDynamicData = n_block_data(iBlockFrom)
    BlockData_I(14)  = real(nDynamicData)

    iData = nScalarBLK

    if(DoMoveData)then
       if (UseConstrainB) then                      !^CFG IF CONSTRAINB BEGIN
          do iVar=1,Bx_-1
             do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
                iData = iData+1
                BlockData_I(iData) = State_VGB(iVar,i,j,k,iBlockFrom)
             end do; end do; end do
          end do
          do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
             iData = iData+1
             BlockData_I(iData) = Bxface_BLK(i,j,k,iBlockFrom)
          end do; end do; end do
          do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
             iData = iData+1
             BlockData_I(iData) = Byface_BLK(i,j,k,iBlockFrom)
          end do; end do; end do
          do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
             iData = iData+1
             BlockData_I(iData) = Bzface_BLK(i,j,k,iBlockFrom)
          end do; end do; end do
          do iVar=Bz_+1,nVar
             do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
                iData = iData+1
                BlockData_I(iData) = State_VGB(iVar,i,j,k,iBlockFrom)
             end do; end do; end do
          end do
       else                                         !^CFG END CONSTRAINB
          do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn;do i=1-gcn,nI+gcn;do iVar=1,nVar
             iData = iData+1
             BlockData_I(iData) = State_VGB(iVar,i,j,k,iBlockFrom)
          end do; end do; end do; end do
       end if                                       !^CFG IF CONSTRAINB

       if(DoMoveExtraData)then
          if(UseB0)then
             ! B0*Cell
             do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn;do i=1-gcn,nI+gcn
                BlockData_I(iData+1:iData+3) = B0_DGB(:,i,j,k,iBlockFrom)
                iData = iData+3
             end do; end do; end do
             if(.not.UseCovariant)then
                ! B0*Face_x
                do i=1,2; do k=1,nK; do j=1,nJ 
                   BlockData_I(iData+1:iData+3) = &
                        B0ResChange_DXSB(:,j,k,i,iBlockFrom)
                   iData = iData+3
                end do; end do; end do
                
                ! B0*Face_y
                do j=3,4; do k=1,nK; do i=1,nI
                   BlockData_I(iData+1:iData+3) = &
                        B0ResChange_DYSB(:,i,k,j,iBlockFrom)
                   iData = iData+3
                end do; end do; end do
                
                ! B0*Face_z
                do k=5,6; do j=1,nJ; do i=1,nI
                   BlockData_I(iData+1:iData+3) = &
                        B0ResChange_DZSB(:,i,j,k,iBlockFrom)
                   iData = iData+3
                end do; end do; end do
             end if
          end if
          ! fbody*
          if(UseGravity.or.UseRotatingFrame)then
             do k=1,nK; do j=1,nJ; do i=1,nI
                iData = iData+1
                BlockData_I(iData) = fBody_x_BLK(i,j,k,iBlockFrom)
                iData = iData+1
                BlockData_I(iData) = fBody_y_BLK(i,j,k,iBlockFrom)
                iData = iData+1
                BlockData_I(iData) = fBody_z_BLK(i,j,k,iBlockFrom)
             end do; end do; end do
          end if

       end if ! DoMoveExtraData

       if(UseBDF2 .and. n_prev > 0)then             !^CFG IF IMPLICIT BEGIN
          do iw=1,nw; do k=1,nK; do j=1,nJ; do i=1,nI; iData = iData+1
             BlockData_I(iData) = w_prev(i,j,k,iw,iBlockFrom)
          end do; end do; end do; end do
       end if                                       !^CFG END IMPLICIT

       if(DoSendRay)then                            !^CFG IF RCM BEGIN
          do k=1,nK; do j=1,nJ; do i=1,nI; do i2=1,2; do i1=1,3
             iData = iData+1
             BlockData_I(iData) = ray(i1,i2,i,j,k,iBlockFrom)
          end do; end do; end do; end do; end do
       end if                                       !^CFG END RCM

       if(nDynamicData > 0)then
          call get_block_data(iBlockFrom, nDynamicData, &
               BlockData_I(iData+1:iData+nDynamicData))
          iData = iData + nDynamicData
       endif

       if(iData > nDataBLK)then
          write(*,*)'ERROR in load_balance: iData=',iData,&
               ' > nDataBLK=',nDataBLK
          call CON_stop('load_balnce: increase nDataBLK')
       end if


    end if

    if(DoTest)write(*,*)'sending BlockData_I: iData=',iData,' from',&
         iProc,' to',iProcTo
    itag=1
    call MPI_SEND(BlockData_I, iData, MPI_REAL, iProcTo, &
         itag, iComm, iError)

    if(DoTest)write(*,*)'send done, me,iBlockFrom,nDynamic=', &
         iProc,iBlockFrom,nDynamicData

  end subroutine send_block_data

  !============================================================================

  subroutine recv_block_data

    globalBLK = iBlockTo

    if(DoTest)write(*,*)'recv BlockData_I: iData=',nDataBLK,' from',&
         iProcFrom,' to',iProc

    if(DoMoveData)then
       ! This is only an upper estimate of the number of reals received
       iSize = nDataBLK
       if(.not.DoMoveExtraData) &
            iSize = iSize - nExtraData
       if(.not.(UseBDF2 .and. n_prev > 0)) &    !^CFG IF IMPLICIT
            iSize = iSize - nWIJK               !^CFG IF IMPLICIT
       if(.not.(DoSendRay)) &                   !^CFG IF RCM
            iSize = iSize - 3*2*nIJK            !^CFG IF RCM
    else
       iSize= nScalarBLK
    end if
    itag=1
    call MPI_RECV(BlockData_I, iSize, MPI_REAL, iProcFrom, &
         itag, iComm, Status_I, iError)

    dx_BLK(iBlockTo)           = BlockData_I(1)
    dy_BLK(iBlockTo)           = BlockData_I(2)
    dz_BLK(iBlockTo)           = BlockData_I(3)
    xyzStart_BLK(1:3,iBlockTo) = BlockData_I(4:6)
    dt_BLK(iBlockTo)           = BlockData_I(7)
    neiLev(:,iBlockTo)         = nint(BlockData_I(8:13))
    nDynamicData               = nint(BlockData_I(14))

    if(DoTest)write(*,*)'recv done, me,iBlockTo,nDynamic=', &
         iProc,iBlockTo,nDynamicData

    ! Put neighbor info into other arrays 
    ! (used for B0 face restriction)
    neiLeast(iBlockTo)         = neiLev(east_,iBlockTo)
    neiLwest(iBlockTo)         = neiLev(west_,iBlockTo)
    neiLsouth(iBlockTo)        = neiLev(south_,iBlockTo)
    neiLnorth(iBlockTo)        = neiLev(north_,iBlockTo)
    neiLbot(iBlockTo)          = neiLev(bot_,iBlockTo)
    neiLtop(iBlockTo)          = neiLev(top_,iBlockTo)

    ! Fix geometry
    call fix_block_geometry(iBlockTo)

    if (.not.DoMoveData) RETURN

    ! Read rest of the blockData buffer
    iData = nScalarBLK
    if (UseConstrainB) then                      !^CFG IF CONSTRAINB BEGIN
       do iVar=1,Bx_-1
          do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
             iData = iData+1
             State_VGB(iVar,i,j,k,iBlockTo) = BlockData_I(iData)
          end do; end do; end do
       end do
       do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
          iData = iData+1
          Bxface_BLK(i,j,k,iBlockTo) = BlockData_I(iData) 
       end do; end do; end do
       do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
          iData = iData+1
          Byface_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
       end do; end do; end do
       do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
          iData = iData+1
          Bzface_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
       end do; end do; end do
       do iVar=Bz_+1,nVar
          do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
             iData = iData+1
             State_VGB(iVar,i,j,k,iBlockTo) = BlockData_I(iData)
          end do; end do; end do
       end do
    else                                         !^CFG END CONSTRAINB
       do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn; do iVar=1,nVar
          iData = iData+1
          State_VGB(iVar,i,j,k,iBlockTo) = BlockData_I(iData)
       end do; end do; end do; end do
    end if                                      !^CFG IF CONSTRAINB

    if(DoMoveExtraData)then
       if(UseB0)then
          ! B0*Cell
          do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn;do i=1-gcn,nI+gcn
             B0_DGB(:,i,j,k,iBlockTo) = BlockData_I(iData+1:iData+3)
             iData = iData+3
          end do; end do; end do
          if(.not.UseCovariant)then
             ! B0*Face_x
             do i=1,2; do k=1,nK; do j=1,nJ 
                B0ResChange_DXSB(:,j,k,i,iBlockTo) = &
                     BlockData_I(iData+1:iData+3)
                iData = iData+3
             end do; end do; end do
             
             ! B0*Face_y
             do j=3,4; do k=1,nK; do i=1,nI
                B0ResChange_DYSB(:,i,k,j,iBlockTo) = &
                     BlockData_I(iData+1:iData+3)
                iData = iData+3
             end do; end do; end do
             
             ! B0*Face_z
             do k=5,6; do j=1,nJ; do i=1,nI
                B0ResChange_DZSB(:,i,j,k,iBlockTo) = &
                     BlockData_I(iData+1:iData+3)
                iData = iData+3
             end do; end do; end do
          end if
       end if
       ! fbody*
       if(UseGravity.or.UseRotatingFrame)then
          do k=1,nK; do j=1,nJ; do i=1,nI
             iData = iData+1
             fBody_x_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
             iData = iData+1
             fBody_y_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
             iData = iData+1
             fBody_z_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
          end do; end do; end do
       else
          fbody_x_BLK(:,:,:,iBlockTo) = 0.0
          fbody_y_BLK(:,:,:,iBlockTo) = 0.0
          fbody_z_BLK(:,:,:,iBlockTo) = 0.0
       end if

    end if ! DoMoveExtraData

    if(UseBDF2 .and. n_prev > 0)then            !^CFG IF IMPLICIT BEGIN
       do iw=1,nw; do k=1,nK; do j=1,nJ; do i=1,nI
          iData = iData+1
          w_prev(i,j,k,iw,iBlockTo) = BlockData_I(iData)
       end do; end do; end do; end do
    end if                                      !^CFG END IMPLICIT

    if(DoSendRay)then                           !^CFG IF RCM BEGIN
       do k=1,nK; do j=1,nJ; do i=1,nI; do i2=1,2; do i1=1,3
          iData = iData+1
          ray(i1,i2,i,j,k,iBlockTo) = BlockData_I(iData)
       end do; end do; end do; end do; end do
    end if                                      !^CFG END RCM

    if(nDynamicData > 0)then
       call put_block_data(iBlockTo, nDynamicData, &
            BlockData_I(iData+1:iData+nDynamicData))
       call set_block_data(iBlockTo)
    end if

  end subroutine recv_block_data
  !============================================================================
end subroutine move_block

!^CFG IF IMPLICIT BEGIN
!=============================================================================
subroutine select_stepping(DoPartSelect)
  ! Set logical arrays for implicit blocks, 
  ! set number of implicit and explicit blocks,
  ! and if DoPartSelect is true then select explicit and implicit blocks
  ! based on the stepping selection criteria.

  use ModProcMH
  use ModMain
  use ModFaceFlux, ONLY : calc_face_flux
  use ModFaceValue,ONLY : calc_face_value
  use ModAdvance,  ONLY : iTypeAdvance_B, iTypeAdvance_BP, &
       SkippedBlock_, ExplBlock_, ImplBlock_
  use ModGeometry, ONLY : Rmin_BLK
  use ModImplicit, ONLY : UseImplicit, UseFullImplicit, UsePartImplicit, &
       UseSemiImplicit, ImplCritType, ExplCFL, rImplicit
  use ModIO,       ONLY : write_prefix, iUnitOut
  use ModB0,ONLY:set_b0_face
  use ModMpi
  implicit none

  logical, intent(in) :: DoPartSelect

  integer :: nBlockExpl, nBlockImpl
  integer :: iError
  logical :: oktest, oktest_me
  !---------------------------------------------------------------------------
  call set_oktest('select_stepping',oktest,oktest_me)

  if(oktest)then
     write(*,*) 'select_stepping starting with iProc ',&
          'UseFullImplicit, UsePartImplicit, DoPartSelect=',&
          iProc, UseFullImplicit, UsePartImplicit, DoPartSelect
     write(*,*)'select_stepping starting with ',&
          'iProc,nBlockExpl, nBlockImpl, nBlockSkipped=', iProc, &
          count(iTypeAdvance_B == ExplBlock_),&
          count(iTypeAdvance_B == ImplBlock_),&
          count(iTypeAdvance_B == SkippedBlock_)
     write(*,*)'select_stepping starting with ',&
          'iProc,nBlockExplALL, nBlockImplALL, nBlockSkippedALL=', iProc, &
          count(iTypeAdvance_BP == ExplBlock_),&
          count(iTypeAdvance_BP == ImplBlock_),&
          count(iTypeAdvance_BP == SkippedBlock_)
     write(*,*)'select_stepping starting with ',&
          'iProc,min(advance),max(advance)=',iProc, &
          minval(iTypeAdvance_B),maxval(iTypeAdvance_B)
     write(*,*)'select_stepping starting with ',&
          'iProc,min(advanceall),max(advanceall)=',iProc, &
          minval(iTypeAdvance_BP),maxval(iTypeAdvance_BP)
  end if

  if(UseFullImplicit .or. UseSemiImplicit)then
     nBlockExplALL = 0
     nBlockImplALL = nBlockALL
     where(iTypeAdvance_BP(1:nBlockMax,:) /= SkippedBlock_) &
          iTypeAdvance_BP(1:nBlockMax,:) = ImplBlock_
     iTypeAdvance_B(1:nBlockMax) = iTypeAdvance_BP(1:nBlockMax,iProc)

  elseif((UsePartImplicit .and. .not. DoPartSelect) .or. .not. UseImplicit)then
     nBlockExplALL    = nBlockALL
     nBlockImplALL    = 0
     where(iTypeAdvance_BP(1:nBlockMax,:) == ImplBlock_) &
          iTypeAdvance_BP(1:nBlockMax,:) = ExplBlock_
     iTypeAdvance_B(1:nBlockMax) = iTypeAdvance_BP(1:nBlockMax,iProc)
          
  else
     ! First set all blocks to be explicit
     where(iTypeAdvance_B(1:nBlockMax) /= SkippedBlock_) &
          iTypeAdvance_B(1:nBlockMax) = ExplBlock_

     ! Select implicitly treated blocks
     select case(ImplCritType)
     case('dt')
        ! Just checking
        if(.not.time_accurate)call stop_mpi(&
             'ImplCritType=dt is only valid in time_accurate mode')

        ! Set implicitBLK based on the time step.
        do globalBLK=1,nBlockMax
           if(unusedBLK(globalBLK)) CYCLE

           ! Obtain the time step based on CFL condition

           ! For first iteration calculate dt_BLK when inside time loop,
           ! otherwise use the available dt_BLK from previous time step,
           ! or from the restart file, or simply 0 set in read_inputs.
           ! The latter two choices will be overruled later anyways.
           if(n_step==0 .and. time_loop)then
              ! For first iteration in the time loop
              ! calculate stable time step
              call set_b0_face(globalBLK)
              call calc_face_value(.false., GlobalBlk)
              call calc_face_flux(.false., GlobalBlk)
              call calc_timestep
           end if

           ! If the smallest allowed timestep is below the fixed DtFixed
           ! then only implicit scheme will work
           if(dt_BLK(globalBLK)*explCFL <= DtFixed) &
                iTypeAdvance_B(globalBLK) = ImplBlock_
        end do

        if(oktest_me)write(*,*)&
             'SELECT: advancetype,dt_BLK,explCFL,dt=',&
             iTypeAdvance_B(BLKtest),dt_BLK(BLKtest),explCFL,dt

     case('r','R')
        ! implicitly treated blocks are within rImplicit and not Unused
        where(rMin_BLK(1:nBlockMax) <= rImplicit .and. &
             .not.UnusedBLK(1:nBlockMax)) &
             iTypeAdvance_B(1:nBlockMax) = ImplBlock_
     case('test')
        if(iProc==PROCtest) iTypeAdvance_B(BLKtest) = ImplBlock_
     end select

     ! Gather global information
     call MPI_ALLGATHER(iTypeAdvance_B, MaxBlock, MPI_INTEGER, &
          iTypeAdvance_BP, MaxBlock, MPI_INTEGER, iComm, iError)

     nBlockImplALL = count(iTypeAdvance_BP == ImplBlock_)
     nBlockExplALL = count(iTypeAdvance_BP == ExplBlock_)

     if(iProc==0.and.lVerbose>0)then
        call write_prefix; 
        write(iUnitOut,*)'select_stepping: nBlockExplALL, nBlockImplALL=',&
             nBlockExplALL,nBlockImplALL
     end if

  end if
  if(oktest)then
     write(*,*)'select_stepping finished with ',&
          'iProc,nBlockExpl, nBlockImpl, nBlockSkipped=', iProc, &
          count(iTypeAdvance_B == ExplBlock_),&
          count(iTypeAdvance_B == ImplBlock_),&
          count(iTypeAdvance_B == SkippedBlock_)
     write(*,*)'select_stepping finished with ',&
          'iProc,nBlockExplALL, nBlockImplALL, nBlockSkippedALL=', iProc, &
          count(iTypeAdvance_BP == ExplBlock_),&
          count(iTypeAdvance_BP == ImplBlock_),&
          count(iTypeAdvance_BP == SkippedBlock_)
     write(*,*)'select_stepping finished with ',&
          'iProc,min(advance),max(advance)=',iProc, &
          minval(iTypeAdvance_B),maxval(iTypeAdvance_B)
     write(*,*)'select_stepping finished with ',&
          'iProc,min(advanceall),max(advanceall)=',iProc, &
          minval(iTypeAdvance_BP),maxval(iTypeAdvance_BP)
  end if
end subroutine select_stepping
!^CFG END IMPLICIT
