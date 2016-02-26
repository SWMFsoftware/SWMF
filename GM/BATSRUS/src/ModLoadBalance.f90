!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModLoadBalance

  use ModMain, ONLY: UseConstrainB, UseB0, UseIM
  use BATL_size, ONLY: nI, nJ, nK, nIJK, &
       MinI, MaxI, MinJ, MaxJ, MinK, MaxK, MaxIJK
  use ModBlockData, ONLY: MaxBlockData, get_block_data, put_block_data, &
       n_block_data, use_block_data, set_block_data, clean_block_data
  use ModImplicit, ONLY: UseImplicit, UseBDF2, n_prev, ImplOld_VCB
  use ModCT, ONLY: Bxface_BLK,Byface_BLK,Bzface_BLK
  use ModRaytrace, ONLY: ray
  use ModAdvance, ONLY: nVar
  use ModB0, ONLY: B0_DGB
  use ModIo, ONLY: log_vars

  implicit none

  private ! except

  public:: init_load_balance
  public:: pack_load_balance
  public:: unpack_load_balance

  ! Upper estimate on the size of buffer to be sent or received
  integer, public:: nBuffer

  ! This could be public, so one can tell if we need to recalculate stuff
  logical, parameter :: DoMoveExtraData = .false.

  ! Local variables

  ! Number of scalars passed in the buffer
  integer, parameter :: nScalarData = 1

  logical:: DoSendRay

contains
  !============================================================================
  subroutine init_load_balance

    DoSendRay = UseIM .or. index(log_vars, 'status') > 0  !!! to be improved

    nBuffer = nScalarData

    if(UseConstrainB) nBuffer = nBuffer + 3*MaxIJK
    if(DoSendRay) &
         nBuffer = nBuffer + 6*nIJK
    if(UseImplicit .and. UseBDF2 .and. n_prev > 0) &
         nBuffer = nBuffer + nVar*nIJK
    if(DoMoveExtraData)then
       if(UseB0) nBuffer = nBuffer + 3*MaxIJK
    end if
    if(MaxBlockData > 0) &
         nBuffer = nBuffer + MaxBlockData
    
  end subroutine init_load_balance

  !============================================================================
  subroutine pack_load_balance(iBlock, nBuffer, Buffer_I)

    integer, intent(in) :: iBlock
    integer, intent(in) :: nBuffer
    real,    intent(out):: Buffer_I(nBuffer)

    integer:: i, j, k, i1, i2, iVar, iData, nDynamicData

    character(len=*), parameter:: NameSub = 'pack_load_balance'
    !------------------------------------------------------------------------

    ! Amount of user defined data for this block
    nDynamicData = 0
    if(use_block_data(iBlock)) nDynamicData = n_block_data(iBlock)
    if(nDynamicData > MaxBlockData)then
       write(*,*)NameSub,' ERROR: iBlock, nDynamicData=', iBlock, nDynamicData
       write(*,*)NameSub,' MaxBlockData=',MaxBlockData,' is too small!'
       call stop_mpi(NameSub//': MaxBlockData has to be set in ModUser')
    end if

    Buffer_I(1)  = real(nDynamicData)

    iData = nScalarData

    if (UseConstrainB) then

       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          iData = iData+1
          Buffer_I(iData) = Bxface_BLK(i,j,k,iBlock)
          iData = iData+1
          Buffer_I(iData) = Byface_BLK(i,j,k,iBlock)
          iData = iData+1
          Buffer_I(iData) = Bzface_BLK(i,j,k,iBlock)
       end do; end do; end do

    endif

    if(DoMoveExtraData)then
       if(UseB0)then
          ! Cell centered B0_DGB
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Buffer_I(iData+1:iData+3) = B0_DGB(:,i,j,k,iBlock)
             iData = iData + 3
          end do; end do; end do
       end if

    end if ! DoMoveExtraData

    if(UseImplicit .and. UseBDF2 .and. n_prev > 0)then
       do k=1,nK; do j=1,nJ; do i=1,nI; do iVar=1,nVar; iData = iData+1
          Buffer_I(iData) = ImplOld_VCB(iVar,i,j,k,iBlock)
       end do; end do; end do; end do
    end if

    if(DoSendRay)then
       do k=1,nK; do j=1,nJ; do i=1,nI; do i2=1,2; do i1=1,3
          iData = iData+1
          Buffer_I(iData) = ray(i1,i2,i,j,k,iBlock)
       end do; end do; end do; end do; end do
    end if

    if(nDynamicData > 0)then
       call get_block_data(iBlock, nDynamicData, &
            Buffer_I(iData+1:iData+nDynamicData))
       iData = iData + nDynamicData
       call clean_block_data(iBlock)
    endif

    if(iData > nBuffer)then
       write(*,*)'ERROR in load_balance: iData=',iData,&
            ' > nBuffer=',nBuffer
       call CON_stop(NameSub//'load_balnce: increase nBuffer')
    end if

  end subroutine pack_load_balance

  !==========================================================================

  subroutine unpack_load_balance(iBlock, nBuffer, Buffer_I)

    integer, intent(in):: iBlock
    integer, intent(in):: nBuffer
    real,    intent(in):: Buffer_I(nBuffer)

    integer:: i, j, k, i1, i2, iVar, iData, nDynamicData

    character(len=*), parameter:: NameSub = 'unpack_load_balance'
    !------------------------------------------------------------------------

    ! Amount of user defined data for this block
    nDynamicData = nint(Buffer_I(1))

    ! Read rest of the blockData buffer
    iData = nScalarData

    if (UseConstrainB) then
       do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
          iData = iData+1
          Bxface_BLK(i,j,k,iBlock) = Buffer_I(iData) 
          iData = iData+1
          Byface_BLK(i,j,k,iBlock) = Buffer_I(iData)
          iData = iData+1
          Bzface_BLK(i,j,k,iBlock) = Buffer_I(iData)
       end do; end do; end do
    end if

    if(DoMoveExtraData)then
       if(UseB0)then
          ! Cell centered B0_DGB
          do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
             B0_DGB(:,i,j,k,iBlock) = Buffer_I(iData+1:iData+3)
             iData = iData+3
          end do; end do; end do
       end if
    end if ! DoMoveExtraData

    if(UseImplicit .and. UseBDF2 .and. n_prev > 0)then
       do k=1,nK; do j=1,nJ; do i=1,nI; do iVar=1,nVar
          iData = iData+1
          ImplOld_VCB(iVar,i,j,k,iBlock) = Buffer_I(iData)
       end do; end do; end do; end do
    end if

    if(DoSendRay)then
       do k=1,nK; do j=1,nJ; do i=1,nI; do i2=1,2; do i1=1,3
          iData = iData+1
          ray(i1,i2,i,j,k,iBlock) = Buffer_I(iData)
       end do; end do; end do; end do; end do
    end if

    if(nDynamicData > 0)then
       call put_block_data(iBlock, nDynamicData, &
            Buffer_I(iData+1:iData+nDynamicData))
       call set_block_data(iBlock)
    end if

  end subroutine unpack_load_balance
  !===========================================================================

end module ModLoadBalance
!=============================================================================

subroutine load_balance(DoMoveCoord, DoMoveData, IsNewBlock)

  use ModLoadBalance, ONLY: nBuffer, &
       init_load_balance, pack_load_balance, unpack_load_balance

  use ModProcMH
  use ModMain
  use ModImplicit, ONLY : UsePartImplicit, &
       TypeSemiImplicit, iBlockFromSemi_B, nBlockSemi
  use ModAdvance, ONLY: iTypeAdvance_B, iTypeAdvance_BP,&
       SkippedBlock_, SteadyBoundBlock_, ExplBlock_, ImplBlock_, SteadyBlock_,&
       State_VGB, iNotSkippedBlock, iTrueBlock, iImplBlock, iSemiImplBlock,&
       iFieldLineThreadsBlock, iPicBlock, iSteadyBlock, iSteadyBoundBlock,&
       iHighOrderBlock, IsLowOrderOnly_B
  use ModGeometry,   ONLY: True_Blk, true_cell, far_field_BCs_Blk
  use ModPartSteady, ONLY: UsePartSteady
  use BATL_lib,      ONLY: Unused_BP
  use ModParallel
  use ModMpi
  use ModPIC, ONLY: UsePic, pic_find_node, IsPicNode_A, DoBalancePicBlock
  
  use BATL_lib, ONLY: MaxNode, nNode, iTree_IA, Status_, Proc_, Block_, Used_,&
       regrid_batl, IsCartesianGrid, iNode_B
  use ModBatlInterface, ONLY: set_batsrus_grid, set_batsrus_state
  use ModLocalTimeStep, ONLY: UseLocalTimeStep
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
  ! There is another declaration in subroutine move_block! Change together!
  logical, parameter :: DoMoveExtraData = .true.

  ! Maximum number of block types to be load balanced separately
  integer, parameter :: MaxType = 10

  ! Actual number of block types
  integer :: nType

  ! Index for block type

  ! Global block index for the various block types

  ! Number of blocks for each type

  ! load balance distribute each type

  ! Number of blocks moved around

  integer :: iError
  integer:: iNode, iBlock, ncount

  integer:: iType, iCrit, nCritBalance

  ! Conversion from iTypeBalance to iType
  integer, allocatable:: iType_I(:)

  logical, allocatable:: IsTypeExist_I(:)
  
  ! We should switch to these variables instead of _BP indexes !!!
  integer, allocatable:: iTypeAdvance_A(:), iTypeBalance_A(:), &
       iTypeBalanceLocal_A(:)

  integer, allocatable:: iTypeConvert_I(:)
  
  logical :: IsSemiImplBlock_B(MaxBlock)
  integer :: i, j
  integer, allocatable:: nTypeProc_PI(:,:)
  
  logical:: DoTest, DoTestMe
  character(len=*), parameter:: NameSub = 'load_balance'
  !---------------------------------------------------------------------------
  call set_oktest(NameSub, DoTest, DoTestMe)
  call timing_start(NameSub)
  if(DoTestMe)write(*,*) NameSub, &
       ' starting with DoMoveCoord,DoMoveData,IsNewBlock=', &
       DoMoveCoord, DoMoveData, IsNewBlock
 
  call select_stepping(DoMoveCoord)

  nCritBalance = 9
  if (nProc>1 .and. index(test_string,'NOLOADBALANCE') < 1) then
     allocate(iTypeBalance_A(MaxNode),iTypeBalanceLocal_A(MaxNode))
     iTypeBalance_A = 0
     iTypeBalanceLocal_A=0
     
     if(DoMoveCoord)then        
        IsSemiImplBlock_B = .false.
        if(TypeSemiImplicit(1:6) == 'resist' .and. nBlockSemi >= 0) &
             IsSemiImplBlock_B(iBlockFromSemi_B(1:nBlockSemi)) = .true.        

        if(UsePic .and. DoBalancePicBlock .and. IsNewBlock) call pic_find_node()
        do iBlock = 1, nBlock
           iType = 0
           ! N criterions are used to decide the type of a block.
           ! Each criterion determins the value of one bit (0 or 1). For N
           ! criterions, the iType of a block could be -1 or a value
           ! between 1 and 2^N-1 (the first bit SkippedBlock is special).
           ! Only part of the total possibilities exist.
           ! 1st bit: NotSkippedBlock-> 1, otherwise -> iType will be -1.
           ! 2nd bit: true_blk       -> 1, otherwise -> 0
           ! 3rd bit: implicit       -> 1, explicit  -> 0
           ! 4th bit: semi-implicit  -> 1, otherwise -> 0
           ! 5th bit: threaded field BC -> 1, otherwise -> 0
           ! 6th bit: pic-block      -> 1, otherwise -> 0
           ! 7th bit: steady block   -> 1, otherwise -> 0
           ! 8th bit: high-order scheme -> 1, otherwise -> 0
           
           if(iTypeAdvance_B(iBlock) ==SkippedBlock_) then              
              ! If it is skipped block, then iType = -1.
              iType = -1
              CYCLE
           else
              iCrit = 1
              iType = iType + iCrit
              
              ! true block 
              iCrit = 2*iCrit
              if(True_Blk(iBlock)) iType = iType + iCrit

              ! implicit
              iCrit = 2*iCrit
              if(iTypeAdvance_B(iBlock)==ImplBlock_) iType = iType + iCrit

              ! semi-implicit
              iCrit = 2*iCrit
              if(IsSemiImplBlock_B(iBlock)) iType = iType + iCrit

              ! threaded field line BC
              iCrit = 2*iCrit
              if(UseFieldLineThreads .and. far_field_BCs_Blk(iBlock) &
                   .and. NeiLev(1,iBlock)==NOBLK) iType = iType + iCrit

              ! PIC region
              iCrit = 2*iCrit
               if(UsePic .and. DoBalancePicBlock) then
                  if(IsPicNode_A(iNode_B(iBlock))) iType = iType + iCrit
               endif

              ! Steady block
               iCrit = 2*iCrit
              if(UsePartSteady .and. iTypeAdvance_B(iBlock)==SteadyBlock_)&
                   iType = iType + iCrit

              ! High-order scheme block: at least part of the faces use
              ! high-order face values. 
              iCrit = 2*iCrit
              if(.not.IsLowOrderOnly_B(iBlock)) iType = iType + iCrit  
           endif
           iTypeBalanceLocal_A(iNode_B(iBlock)) = iType
        enddo
        ! Update iTypeAdvance_BP
        ! CHEATING: only indicate first index to circumvent ModMpiInterfaces
        ! Set displacement equal to MaxBlock so we get same behavior 
        ! as MPI_allgather. Use nBlockMax for maximum receive data for speed
        nBlockMax_P = nBlockMax
        call MPI_allgatherv(iTypeAdvance_B(1), nBlockMax, MPI_INTEGER, &
             iTypeAdvance_BP(1,0), nBlockMax_P, MaxBlockDisp_P,&
             MPI_INTEGER, iComm, iError)

        ! In order to do MPI_allreduce, the default value of iTypeBalanceLocal_A
        ! is set to 0. While the skipped blocks are also 0.         
        call MPI_allreduce(iTypeBalanceLocal_A, iTypeBalance_A, MaxNode, &
             MPI_INTEGER, MPI_SUM, iComm, iError)

     end if
     if(.not.allocated(iType_I)) allocate(iType_I(0:2**nCritBalance-1))     
     if(.not.allocated(IsTypeExist_I))&
          allocate(IsTypeExist_I(0:2**nCritBalance-1))     
     IsTypeExist_I = .false.
     do iNode = 1, nNode
        if(iTypeBalance_A(iNode)>=0)&
             IsTypeExist_I(iTypeBalance_A(iNode)) = .true.        
     enddo

     iType_I = -1
     ncount = 0
     ! Ignore iType=0, the skipped blocks. 
     do iType = 1, 2**nCritBalance-1
        if(IsTypeExist_I(iType)) then
           ncount = ncount + 1
           iType_I(iType) = ncount
        endif
     enddo
     nType = ncount


     ! Convert from block/proc index to node index.
     allocate(iTypeAdvance_A(MaxNode))
     if(DoMoveCoord) then
        do iNode = 1, nNode       
           if(iTree_IA(Status_,iNode) /= Used_) CYCLE
           iTypeAdvance_A(iNode) = &
                iTypeAdvance_BP(iTree_IA(Block_,iNode),iTree_IA(Proc_,iNode))
           iTypeBalance_A(iNode) = iType_I(iTypeBalance_A(iNode))
        end do
     else         
        iTypeBalance_A = 1
        iTypeAdvance_A = 1
     endif
     
     ! load balance depending on block types
     if(DoMoveData)then
        call init_load_balance

        call regrid_batl(nVar, State_VGB, Dt_BLK,  &
             DoBalanceEachLevelIn=UseLocalTimeStep,&
             iTypeBalance_A=iTypeBalance_A,        &
             iTypeNode_A=iTypeAdvance_A,           &
             DoBalanceOnlyIn=.true.,               &
             nExtraData=nBuffer,                   &
             pack_extra_data=pack_load_balance,    &
             unpack_extra_data=unpack_load_balance )

        call set_batsrus_grid
        call set_batsrus_state
     else        
        call regrid_batl(nVar, State_VGB, Dt_BLK, Used_GB=true_cell, &
             iTypeBalance_A=iTypeBalance_A, iTypeNode_A=iTypeAdvance_A)
        call set_batsrus_grid
     end if

     ! Go up to nBlockMax instead of nBlock (!) to clean up beyond nBlock
     iTypeAdvance_BP(1:nBlockMax,:) = SkippedBlock_
     do iNode = 1, nNode
        if(iTree_IA(Status_,iNode) /= Used_) CYCLE
        iTypeAdvance_BP(iTree_IA(Block_,iNode),iTree_IA(Proc_,iNode)) = &
             iTypeAdvance_A(iNode)
     end do
     iTypeAdvance_B(1:nBlockMax)  = iTypeAdvance_BP(1:nBlockMax,iProc)
     if(DoTestMe .and. nType>0) then
        ! Create the table so that we can check the status of a block.
        if(.not.allocated(iTypeConvert_I)) allocate(iTypeConvert_I(nType))
        ncount=1
        do iType = 1, 2**nCritBalance-1
           if(iType_I(iType)>0) then
              iTypeConvert_I(ncount) = iType
              ncount = ncount + 1           
           endif
        enddo
        allocate(nTypeProc_PI(0:nProc-1,nType))
        nTypeProc_PI = 0
        do iNode = 1, nNode
           if(iTree_IA(Status_,iNode) /= Used_) CYCLE
           nTypeProc_PI(iTree_IA(Proc_,iNode),iTypeBalance_A(iNode)) = &
                nTypeProc_PI(iTree_IA(Proc_,iNode),iTypeBalance_A(iNode))+1
        end do
        do i=1, nType
           write(*,*) 'iType= ',iTypeConvert_I(i),' nCount= ',sum(nTypeProc_PI(:,i))
        enddo       
        deallocate(nTypeProc_PI)
     endif ! DoTestMe
     
     deallocate(iTypeAdvance_A, iTypeBalance_A, iTypeBalanceLocal_A)
  end if
  
  ! When load balancing is done Skipped and Unused blocks coincide
  Unused_BP(1:nBlockMax,:) = &
       iTypeAdvance_BP(1:nBlockMax,:) == SkippedBlock_

  call find_test_cell

  if(allocated(iType_I))       deallocate(iType_I)
  if(allocated(IsTypeExist_I)) deallocate(IsTypeExist_I)
  call timing_stop(NameSub)
end subroutine load_balance


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
       ImplCritType, ExplCFL, rImplicit
  use ModIO,       ONLY: write_prefix, iUnitOut
  use ModB0,       ONLY: set_b0_face
  use ModMpi
  use ModParallel, ONLY: nBlockMax_P, MaxBlockDisp_P
  use ModTimeStepControl, ONLY: calc_timestep
  implicit none

  logical, intent(in) :: DoPartSelect

  integer :: iBlock, iError
  logical :: oktest, oktest_me
  !---------------------------------------------------------------------------
  call set_oktest('select_stepping',oktest,oktest_me)

  if(oktest_me)then
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

  if(UseFullImplicit)then
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
        do iBlock=1,nBlockMax
           if(Unused_B(iBlock)) CYCLE

           ! Obtain the time step based on CFL condition

           ! For first iteration calculate dt_BLK when inside time loop,
           ! otherwise use the available dt_BLK from previous time step,
           ! or from the restart file, or simply 0 set in read_inputs.
           ! The latter two choices will be overruled later anyways.
           if(n_step==0 .and. time_loop)then
              ! For first iteration in the time loop
              ! calculate stable time step
              call set_b0_face(iBlock)
              call calc_face_value(.false., iBlock)
              call calc_face_flux(.false., iBlock)
              call calc_timestep(iBlock)
           end if

           ! If the smallest allowed timestep is below the fixed DtFixed
           ! then only implicit scheme will work
           if(dt_BLK(iBlock)*explCFL <= DtFixed) &
                iTypeAdvance_B(iBlock) = ImplBlock_
        end do

        if(oktest_me)write(*,*)&
             'SELECT: advancetype,dt_BLK,explCFL,dt=',&
             iTypeAdvance_B(BLKtest),dt_BLK(BLKtest),explCFL,dt

     case('r','R')
        ! implicitly treated blocks are within rImplicit and not Unused
        where(rMin_BLK(1:nBlockMax) <= rImplicit .and. &
             .not.Unused_B(1:nBlockMax)) &
             iTypeAdvance_B(1:nBlockMax) = ImplBlock_
     case('test')
        if(iProc==PROCtest) iTypeAdvance_B(BLKtest) = ImplBlock_
     end select

     ! Gather global information
     ! CHEATING: only indicate first index to circumvent ModMpiInterfaces.
     ! Set displacement equal to MaxBlock so we get same behavior 
     ! as MPI_allgather. Use nBlockMax for maximum receive data for speed!
     nBlockMax_P = nBlockMax
     call MPI_allgatherv(iTypeAdvance_B(1), nBlockMax, MPI_INTEGER, &
          iTypeAdvance_BP(1,0), nBlockMax_P, MaxBlockDisp_P,&
          MPI_INTEGER, iComm, iError)

     nBlockImplALL = count(iTypeAdvance_BP(1:nBlockMax,:) == ImplBlock_)
     nBlockExplALL = count(iTypeAdvance_BP(1:nBlockMax,:) == ExplBlock_)

     if(iProc==0.and.lVerbose>0)then
        call write_prefix; 
        write(iUnitOut,*)'select_stepping: nBlockExplALL, nBlockImplALL=',&
             nBlockExplALL,nBlockImplALL
     end if

  end if
  if(oktest_me)then
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

