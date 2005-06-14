!^CFG COPYRIGHT UM
subroutine load_balance(DoMoveCoord, DoMoveData, nBlockMoved)
  use ModProcMH
  use ModMain
  use ModImplicit, ONLY : implicitBLK, UsePartImplicit !^CFG IF IMPLICIT
  use ModAMR, ONLY : availableBLKs
  use ModParallel
  use ModIO
  use ModMpi
  implicit none

  ! Load balance grid using Peano-Hilbert ordering of blocks
  ! Coordinates are moved if DoMoveCoord is true.
  ! Data is moved with the blocks if DoMoveData is true.

  logical, intent(in) :: DoMoveCoord, DoMoveData
  integer, intent(out):: nBlockMoved

  ! Maximum number of attempts to accomplish the load balancing
  ! The algorithm needs multiple tries if the actual number of blocks
  ! is very close to the maximum number of blocks and many blocks are moved
  integer, parameter :: MaxTry = 100
  
  ! Set this logical to .false. to return to the previous version,
  ! which did not move the B0, body force and heating variables.
  ! There is another declaration in subroutine move_block! Change together!!!
  logical, parameter :: DoMoveExtraData = .true.

  integer :: iError
  integer :: iBlockALL, iBlock
  integer :: iBlockExplALL, iBlockImplALL     !^CFG IF IMPLICIT

  integer :: iBlockFrom, iProcFrom, iBlockTo, iProcTo, iTry

  logical :: SkippedAnyBlock

  logical :: DoTest, DoTestMe

  logical :: DoFixVar_B(MaxBlock)

  !---------------------------------------------------------------------------
  call set_oktest('load_balance',DoTest,DoTestMe)

  if (DoTestMe) write(*,*)'load_balance: DoMoveCoord, DoMoveData=',&
       DoMoveCoord, DoMoveData

  ! starting value of number of blocks moved between processors
  nBlockMoved = 0

  ! Find the last used block on the processor
  do iBlock = nBlockMax,1,-1
     nBlock = iBlock
     if (.not.unusedBLK(iBlock)) EXIT
  end do

  if(DoMoveData.and. .not.DoMoveCoord)call stop_mpi(&
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

  if(DoTestMe)write(*,*)'load_balance starting: nBlockMax=',nBlockMax
  if(DoTest)write(*,*)'load_balance starting: me, nBlock, nBlockUsed=',&
       iProc, nBlock, count(.not.unusedBLK(1:nBlock))
  if(UseImplicit)then                              !^CFG IF IMPLICIT BEGIN
     if(DoTestMe)write(*,*)'load_balance starting: ',&
          'me, nBlockExplALL, nBlockImplALL=',&
          iProc,nBlockExplALL, nBlockImplALL
     if(DoTest)write(*,*)&
          'load_balance starting: me, nBlockImpl=',&
          iProc, count(implicitBLK(1:nBlock))
  end if                                           !^CFG END IMPLICIT

  !DEBUG
  !write(*,*)'START me, nBlock, iBlockALL_B=',iProc, &
  !     nBlock, (global_block_number(iBlock),unusedBLK(iBlock),&
  !     iBlock=1,nBlock)

  if (nProc==1) RETURN

  if (index(test_string,'NOLOADBALANCE')>0) RETURN

  call timing_start('load_balance')

  if(DoMoveData) DoFixVar_B = .false. ! initialize variable fixing flags

  TRY: do iTry=1,MaxTry

     skippedAnyBlock=.false.
     iBlockExplALL = 0; iBlockImplALL = 0           !^CFG IF IMPLICIT

     TREE: do iBlockALL  = 1, nBlockALL

        iBlockFrom = iBlock_A(iBlockALL)
        iProcFrom  = iProc_A(iBlockALL)

        iProcTo = (nProc*iBlockALL-1)/nBlockALL

        if(implicitBlock_BP(iBlockFrom,iProcFrom))then  !^CFG IF IMPLICIT BEGIN
           iBlockImplALL = iBlockImplALL + 1
           iProcTo = (nProc*iBlockImplALL-1) / nBlockImplALL
        else
           iBlockExplALL = iBlockExplALL + 1
           iProcTo = (nProc*iBlockExplALL-1) / nBlockExplALL
        end if                                          !^CFG END IMPLICIT

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
        globalBLK = iBlock
        if(useConstrainB) call Bface2Bcenter          !^CFG IF CONSTRAINB
        call correctE

        if(DoMoveExtraData)then
           call set_b0_matrix(iBlock)                 !^CFG IF CARTESIAN
           !call calc_b0source_covar(iBlock)          !^CFG IF NOT CARTESIAN
           
           if(UsePartImplicit)&                            !^CFG IF IMPLICIT
                call init_conservative_facefluxes(iBlock)  !^CFG IF IMPLICIT
        else
           call calc_other_soln_vars(iBlock)
        end if
     end do
     !call timing_stop('load_fix_var')
  end if

  call timing_stop('load_balance')

  if(DoTestMe)write(*,*)'load_balance finished: ',&
       'nTry, nBlockMax, nBlockMoved=',&
       iTry, nBlockMax, nBlockMoved

  if(DoTest)write(*,*)&
       'load_balance finished: me, nBlock, nBlockUsed=',&
       iProc, nBlock, count(.not.unusedBLK(1:nBlock))

  if(DoTest.and.UseImplicit)write(*,*)&               !^CFG IF IMPLICIT
       'load_balance finished: me, nBlockImpl=',&     !^CFG IF IMPLICIT
       iProc, count(implicitBLK(1:nBlock))         !^CFG IF IMPLICIT

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
       fbody_x_BLK, fbody_y_BLK, fbody_z_BLK, qheat_BLK, &
       B0xCell_BLK, B0yCell_BLK, B0zCell_BLK, &
       B0xFace_x_BLK, B0yFace_x_BLK, B0zFace_x_BLK, &
       B0xFace_y_BLK, B0yFace_y_BLK, B0zFace_y_BLK, &
       B0xFace_z_BLK, B0yFace_z_BLK, B0zFace_z_BLK
  use ModGeometry, ONLY : dx_BLK,dy_BLK,dz_BLK,xyzStart_BLK
  use ModParallel
  use ModImplicit                                         !^CFG IF IMPLICIT
  use ModCT, ONLY : Bxface_BLK,Byface_BLK,Bzface_BLK      !^CFG IF CONSTRAINB
  use ModRaytrace, ONLY : ray                             !^CFG IF RCM
  use ModMpi
  implicit none

  logical, intent(in) :: DoMoveCoord, DoMoveData
  integer, intent(in) :: iBlockALL, iBlockFrom, iProcFrom, iBlockTo,iProcTo

  ! Set this logical to .false. to return to the previous version,
  ! which did not move the B0, body force and heating variables.
  ! There is another declaration in subroutine load_balance! Change together!!!
  logical, parameter :: DoMoveExtraData = .true.

  integer, parameter :: nScalarBLK=13, &
       nCellGhostBLK=(nI+4)*(nJ+4)*(nK+4)
  integer, parameter :: nExtraData = &
       3*nCellGhostBLK +                 & ! B0*Cell
       3*((nI+3)*(nJ+2)*(nK+2)           & ! B0*Face_x
       +  (nI+2)*(nJ+3)*(nK+2)           & ! B0*Face_y
       +  (nI+2)*(nJ+2)*(nK+3)) +        & ! B0*Face_x
       4*nIJK                              ! fbody_* and qheat
       
  integer, parameter :: nDataBLK= &
       nwIJK +                           & !^CFG IF IMPLICIT
       3*2*nIJK +                        & !^CFG IF RCM
       nScalarBLK +                      & ! scalars
       nVar*nCellGhostBLK +              & ! State_VGB
       nExtraData                          ! B0, fbody, qheat

  real, dimension(nDataBLK) :: BlockData_I
  integer :: iData, itag, i, j, k, i1,i2, iVar, iw, iSize
  integer :: iError
  integer :: status(MPI_STATUS_SIZE,1)

  logical :: DoTest, DoTestMe
  !----------------------------------------------------------------------------

  if((iProcFrom==PROCtest .and. iBlockFrom==BLKtest).or. &
       (iProcTo  ==PROCtest .and. iBlockTo  ==BLKtest))then
     call set_oktest('move_block',DoTest,DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  end if
  if(DoTestMe)write(*,*)'iBlockALL,iBlockFrom,iProcFrom, iBlockTo,iProcTo=',&
       iBlockALL,iBlockFrom,iProcFrom, iBlockTo,iProcTo

  if (iProc == iProcFrom) then
     if (DoMoveCoord) call send_block_data
     unusedBLK(iBlockFrom)   = .true.
     implicitBLK(iBlockFrom) = .false.      !^CFG IF IMPLICIT
     do
        if (nBlock==0) EXIT
        if (.not.unusedBLK(nBlock)) EXIT
        nBlock = nBlock-1
     end do
  end if

  if (iProc == iProcTo) then
     unusedBLK(iBlockTo) = .false.
     implicitBLK(iBlockTo) = implicitBlock_BP(iBlockFrom,iProcFrom) !^CFG IF IMPLICIT
     if (DoMoveCoord) call recv_block_data
     nBlock = max(nBlock, iBlockTo)
     global_block_number(iBlockTo) = iBlockALL
  end if

  !^CFG IF IMPLICIT BEGIN
  ! Update global implicit block information
  implicitBlock_BP(iBlockTo, iProcTo) = implicitBlock_BP(iBlockFrom,iProcFrom)
  implicitBlock_BP(iBlockFrom, iProcFrom) = .false.
  !^CFG END IMPLICIT

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
          ! B0*Cell
          do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn;do i=1-gcn,nI+gcn
             iData = iData+1
             BlockData_I(iData) = B0xCell_BLK(i,j,k,iBlockFrom)
             iData = iData+1
             BlockData_I(iData) = B0yCell_BLK(i,j,k,iBlockFrom)
             iData = iData+1
             BlockData_I(iData) = B0zCell_BLK(i,j,k,iBlockFrom)
          end do; end do; end do

          ! B0*Face_x
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+2
             iData = iData+1
             BlockData_I(iData) = B0xFace_x_BLK(i,j,k,iBlockFrom)
             iData = iData+1
             BlockData_I(iData) = B0yFace_x_BLK(i,j,k,iBlockFrom)
             iData = iData+1
             BlockData_I(iData) = B0zFace_x_BLK(i,j,k,iBlockFrom)
          end do; end do; end do

          ! B0*Face_y
          do k=0,nK+1; do j=0,nJ+2; do i=0,nI+1
             iData = iData+1
             BlockData_I(iData) = B0xFace_y_BLK(i,j,k,iBlockFrom)
             iData = iData+1
             BlockData_I(iData) = B0yFace_y_BLK(i,j,k,iBlockFrom)
             iData = iData+1
             BlockData_I(iData) = B0zFace_y_BLK(i,j,k,iBlockFrom)
          end do; end do; end do

          ! B0*Face_z
          do k=0,nK+2; do j=0,nJ+1; do i=0,nI+1
             iData = iData+1
             BlockData_I(iData) = B0xFace_z_BLK(i,j,k,iBlockFrom)
             iData = iData+1
             BlockData_I(iData) = B0yFace_z_BLK(i,j,k,iBlockFrom)
             iData = iData+1
             BlockData_I(iData) = B0zFace_z_BLK(i,j,k,iBlockFrom)
          end do; end do; end do

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

          ! heating
          if(UseUserHeating)then
             do k=1,nK; do j=1,nJ; do i=1,nI
                iData = iData+1
                BlockData_I(iData) = qHeat_BLK(i,j,k,iBlockFrom)
             end do; end do; end do
          end if
       end if ! DoMoveExtraData

       if(UseBDF2 .and. n_prev > 0)then             !^CFG IF IMPLICIT BEGIN
          do iw=1,nw; do k=1,nK; do j=1,nJ; do i=1,nI; iData = iData+1
             BlockData_I(iData) = w_prev(i,j,k,iw,iBlockFrom)
          end do; end do; end do; end do
       end if                                       !^CFG END IMPLICIT

       if(UseIM)then                                !^CFG IF RCM BEGIN
          do k=1,nK; do j=1,nJ; do i=1,nI; do i2=1,2; do i1=1,3
             iData = iData+1
             BlockData_I(iData) = ray(i1,i2,i,j,k,iBlockFrom)
          end do; end do; end do; end do; end do
       end if                                       !^CFG END RCM
    end if

    if(DoTest)write(*,*)'sending BlockData_I: iData=',iData,' from',&
         iProc,' to',iProcTo
    itag=1
    call MPI_SEND(BlockData_I, iData, MPI_REAL, iProcTo, &
         itag, iComm, iError)

    if(DoTest)write(*,*)'send done, me=',iProc

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
       if(.not.(UseIM)) &                       !^CFG IF RCM
            iSize = iSize - 3*2*nIJK            !^CFG IF RCM
    else
       iSize= nScalarBLK
    end if
    itag=1
    call MPI_RECV(BlockData_I, iSize, MPI_REAL, iProcFrom, &
         itag, iComm, status, iError)

    if(DoTest)write(*,*)'recv done, me=',iProc

    dx_BLK(iBlockTo)           = BlockData_I(1)
    dy_BLK(iBlockTo)           = BlockData_I(2)
    dz_BLK(iBlockTo)           = BlockData_I(3)
    xyzStart_BLK(1:3,iBlockTo) = BlockData_I(4:6)
    dt_BLK(iBlockTo)           = BlockData_I(7)
    neiLev(:,iBlockTo)         = nint(BlockData_I(8:13))
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
       ! B0*Cell
       do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn;do i=1-gcn,nI+gcn
          iData = iData+1
          B0xCell_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
          iData = iData+1
          B0yCell_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
          iData = iData+1
          B0zCell_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
       end do; end do; end do

       ! B0*Face_x
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+2
          iData = iData+1
          B0xFace_x_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
          iData = iData+1
          B0yFace_x_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
          iData = iData+1
          B0zFace_x_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
       end do; end do; end do
       
       ! B0*Face_y
       do k=0,nK+1; do j=0,nJ+2; do i=0,nI+1
          iData = iData+1
          B0xFace_y_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
          iData = iData+1
          B0yFace_y_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
          iData = iData+1
          B0zFace_y_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
       end do; end do; end do

       ! B0*Face_z
       do k=0,nK+2; do j=0,nJ+1; do i=0,nI+1
          iData = iData+1
          B0xFace_z_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
          iData = iData+1
          B0yFace_z_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
          iData = iData+1
          B0zFace_z_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
       end do; end do; end do

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

       ! heating
       if(UseUserHeating)then
          do k=1,nK; do j=1,nJ; do i=1,nI
             iData = iData+1
             qHeat_BLK(i,j,k,iBlockTo) = BlockData_I(iData)
          end do; end do; end do
       else
          qheat_BLK(:,:,:,iBlockTo) = 0.0
       end if
    end if ! DoMoveExtraData

    if(UseBDF2 .and. n_prev > 0)then            !^CFG IF IMPLICIT BEGIN
       do iw=1,nw; do k=1,nK; do j=1,nJ; do i=1,nI
          iData = iData+1
          w_prev(i,j,k,iw,iBlockTo) = BlockData_I(iData)
       end do; end do; end do; end do
    end if                                      !^CFG END IMPLICIT

    if(UseIM)then                               !^CFG IF RCM BEGIN
       do k=1,nK; do j=1,nJ; do i=1,nI; do i2=1,2; do i1=1,3
          iData = iData+1
          ray(i1,i2,i,j,k,iBlockTo) = BlockData_I(iData)
       end do; end do; end do; end do; end do
    end if                                      !^CFG END RCM

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
  use ModGeometry, ONLY : Rmin_BLK
  use ModImplicit, ONLY : UseFullImplicit,UsePartImplicit, &
       implicitBLK,ImplCritType,explCFL,Rimplicit
  use ModParallel, ONLY : implicitBlock_BP
  use ModMpi
  implicit none

  logical, intent(in) :: DoPartSelect

  integer :: nBlockExpl, nBlockImpl
  integer :: iError
  logical :: oktest, oktest_me
  !---------------------------------------------------------------------------
  call set_oktest('select_stepping',oktest,oktest_me)

  if(oktest_me)write(*,*) 'select_stepping starting with ',&
       'UseFullImplicit, UsePartImplicit, UsePartLocal, DoPartSelect=',&
       UseFullImplicit, UsePartImplicit, UsePartLocal, DoPartSelect

  if(((UsePartLocal .or. UsePartImplicit) .and. .not. DoPartSelect) &
       .or. .not. (UseImplicit .or. UsePartLocal))then
     nBlockExplALL    = nBlockALL
     nBlockImplALL    = 0
     implicitBLK(1:nBlock)           = .false.
     implicitBlock_BP(1:nBlockMax,:) = .false.

  elseif(UseFullImplicit)then
     nBlockExplALL = 0
     nBlockImplALL = nBlockALL
     implicitBLK(1:nBlock)           = .not.unusedBLK(1:nBlock)

     call MPI_ALLGATHER(implicitBLK,      nBLK, MPI_LOGICAL, &
          implicitBlock_BP, nBLK, MPI_LOGICAL, &
          iComm, iError)
  else

     if(iProc==0.and.lVerbose>0)&
          write(*,*)'select_stepping: ImplCritType=',ImplCritType

     ! Select implicitly treated blocks
     select case(ImplCritType)
     case('dt')
        ! Just checking
        if(.not.time_accurate)call stop_mpi(&
             'ImplCritType=dt is only valid in time_accurate mode')

        ! Set implicitBLK based on the time step.
        do globalBLK=1,nBlockMax
           if(unusedBLK(globalBLK))then
              implicitBLK(globalBLK)=.false.
           else
              ! Obtain the time step based on CFL condition

              ! For first iteration calculate dt_BLK when inside time loop,
              ! otherwise use the available dt_BLK from previous time step,
              ! or from the restart file, or simply 0 set in read_inputs.
              ! The latter two choices will be overruled later anyways.
              if(iteration_number==0 .and. time_loop)then
                 ! For first iteration in the time loop
                 ! calculate stable time step
                 call calc_facevalues(.false.)
                 call calc_facefluxes(.false.)
                 call calc_timestep
              end if
              ! If the smallest allowed timestep is below the fixed DtFixed
              ! then only implicit scheme will work
              implicitBLK(globalBLK) = dt_BLK(globalBLK)*explCFL <= DtFixed
           end if
        end do

        if(oktest_me)write(*,*)&
             'SELECT: unused,implicit,dt_BLK,explCFL,dt=',&
             unusedBLK(BLKtest),implicitBLK(BLKtest),dt_BLK(BLKtest),&
             explCFL,dt

     case('r','R')
        ! implicitly treated blocks are within Rimplicit and not unused
        implicitBLK(1:nBlockMax) = &
             Rmin_BLK(1:nBlockMax) <= Rimplicit .and. .not.unusedBLK(1:nBlockMax)
     case('test')
        implicitBLK(1:nBlockMax) = .false.
        if(iProc==PROCtest) implicitBLK(BLKtest) = .true.
     end select

     nBlockImpl = count(implicitBLK(1:nBlock))
     nBlockExpl = count(.not.(unusedBLK(1:nBlock).or.implicitBLK(1:nBlock)))

     call MPI_allreduce(nBlockImpl, nBlockImplALL, 1, MPI_INTEGER, MPI_SUM, &
          iComm, iError)

     call MPI_allreduce(nBlockExpl, nBlockExplALL, 1, MPI_INTEGER, MPI_SUM, &
          iComm, iError)

     call MPI_ALLGATHER(implicitBLK,      nBLK, MPI_LOGICAL, &
          implicitBlock_BP, nBLK, MPI_LOGICAL, &
          iComm, iError)
     if(iProc==0.and.lVerbose>0)oktest_me=.true. ! report for part implicit
  end if
  if(oktest_me)write(*,*)'select_stepping finished with ',&
       'nBlockExplALL, nBlockImplALL=',nBlockExplALL, nBlockImplALL

end subroutine select_stepping
!^CFG END IMPLICIT
