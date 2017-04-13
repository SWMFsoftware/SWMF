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

  SAVE

  private ! except

  public:: init_load_balance
  public:: load_balance
  public:: load_balance_blocks
  public:: select_stepping

  ! Upper estimate on the size of buffer to be sent or received
  integer, public:: nBuffer

  ! We should switch to these variables instead of _BP indexes !!!
  integer, allocatable, public:: iTypeBalance_A(:)


  integer, allocatable:: iTypeAdvance_A(:)


  ! Local variables

  ! This could be public, so one can tell if we need to recalculate stuff
  logical, parameter :: DoMoveExtraData = .false.

  ! Number of scalars passed in the buffer
  integer, parameter :: nScalarData = 1

  logical:: DoSendRay

  ! The index of block type features. 
  integer, parameter:: &
       iNotSkippedBlock     =  1, & ! Skipped block or not.
       iTrueBlock           =  2    ! True block or not.
  integer :: &
       iPartImplBlock       = -1, & ! Implicit or explicit.
       iSemiImplBlock       = -1, & ! Semi-implicit or not.
       iFieldLineThreadBlock= -1, & ! Threaded file line BC or not.
       iPicBlock            = -1, & ! Overlaped with PIC or not.
       iSteadyBlock         = -1, & ! Steady block or not.
       iHighOrderBlock      = -1, & ! Use high-order scheme or not. 
       iSubCycleBlock       = -1    ! First bit of time level info

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

  subroutine load_balance(DoMoveCoord, DoMoveData, IsNewBlock)

    use ModProcMH
    use ModMain
    use ModImplicit, ONLY : UsePartImplicit, UseSemiImplicit, &
         TypeSemiImplicit, iBlockFromSemi_B, nBlockSemi
    use ModAdvance, ONLY: &
         State_VGB, iTypeAdvance_B, iTypeAdvance_BP,                 &
         SkippedBlock_, ImplBlock_, SteadyBlock_, &
         UseLowOrderRegion, IsLowOrderOnly_B
    use ModGeometry,   ONLY: True_Blk, true_cell, far_field_BCs_Blk
    use ModPartSteady, ONLY: UsePartSteady
    use BATL_lib,      ONLY: Unused_BP
    use ModParallel
    use ModMpi
    use ModPIC, ONLY: UsePic, pic_find_node, IsPicNode_A, DoBalancePicBlock

    use BATL_lib, ONLY: MaxNode, nNode, iTree_IA, Status_, Proc_, Block_, &
         Used_, iNode_B, regrid_batl
    use ModBatlInterface, ONLY: set_batsrus_grid, set_batsrus_state
    use ModTimeStepControl, ONLY: UseMaxTimeStep, DtMax, DtMin
    use ModUserInterface ! user_action

    ! Load balance grid using space filling (Morton) ordering of blocks
    ! Coordinates are moved if DoMoveCoord is true.
    ! Data is moved with the blocks if DoMoveData is true.
    ! There are new blocks (due to initial refinement, restart or AMR)
    ! when IsNewBlock is true (so update neighbors etc).

    logical, optional, intent(in) :: DoMoveCoord, DoMoveData, IsNewBlock

    ! Number of different block types
    integer :: nType

    integer :: iError
    integer:: iNode, iBlock

    ! Index for block types
    integer:: iType, iCrit, iTypeMax, nCount

    ! Conversion from iTypeBalance to iType
    integer, allocatable:: iType_I(:)

    logical, allocatable:: IsTypeExist_I(:)

    integer, allocatable:: iTypeConvert_I(:)

    logical :: IsSemiImplBlock_B(MaxBlock)
    integer :: i
    integer, allocatable:: nTypeProc_PI(:,:)

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'load_balance'
    !--------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)
    call timing_start(NameSub)
    if(DoTestMe)write(*,*) NameSub, &
         ' starting with DoMoveCoord,DoMoveData,IsNewBlock=', &
         DoMoveCoord, DoMoveData, IsNewBlock

    call select_stepping(DoMoveCoord)

    if (nProc > 1 .and. index(test_string,'NOLOADBALANCE') < 1) then
       if(.not.allocated(iTypeBalance_A)) allocate(iTypeBalance_A(MaxNode))
       iTypeBalance_A = 0

       if(DoMoveCoord)then        
          IsSemiImplBlock_B = .false.
          if(TypeSemiImplicit(1:6) == 'resist' .and. nBlockSemi >= 0) &
               IsSemiImplBlock_B(iBlockFromSemi_B(1:nBlockSemi)) = .true.

          if(UsePic .and. DoBalancePicBlock .and. IsNewBlock) &
               call pic_find_node

          ! N criteria are used to decide the type of a block.
          ! Each criterion determines the value of one bit (0 or 1). 
          ! For N criteria, the iType of a block could be -1 or a value
          ! between 1 and 2^N-1 (the first bit SkippedBlock is special).
          ! In a typical case only some of the possible combinations occur.
          ! 

          ! 1st  bit: NotSkippedBlock   -> 1, otherwise -> iType will be -1.
          ! 2nd  bit: true_blk          -> 1, otherwise -> 0
          iCrit = 2

          ! next bit: implicit          -> 1, explicit  -> 0
          if(UsePartImplicit)then
             iCrit = 2*iCrit
             iPartImplBlock = iCrit
          end if

          ! next bit: semi-implicit     -> 1, otherwise -> 0
          if(UseSemiImplicit)then
             iCrit = 2*iCrit
             iSemiImplBlock = iCrit
          end if

          ! next bit: threaded field BC -> 1, otherwise -> 0
          if(UseFieldLineThreads)then
             iCrit = 2*iCrit
             iFieldLineThreadBlock = iCrit
          end if

          ! next bit: PIC-block         -> 1, otherwise -> 0
          if(UsePic .and. DoBalancePicBlock)then
             iCrit = 2*iCrit
             iPicBlock = iCrit
          endif

          ! next bit: steady block      -> 1, otherwise -> 0
          if(UsePartSteady)then
             iCrit = 2*iCrit
             iSteadyBlock = iCrit
          end if

          ! next bit: high-order scheme -> 1, otherwise -> 0
          if(UseLowOrderRegion)then
             ! High-order scheme block: some faces use high-order scheme
             iCrit = 2*iCrit
             iHighOrderBlock = iCrit
          end if

          ! Maximum possible value is when all the bits are 1 from 0 to iCrit.
          iTypeMax = max(1, 2*iCrit - 1)

          ! Time level uses the last nTimeLevel bits
          if(UseMaxTimeStep .and. DtMax > DtMin .and. DtMin > 0)then

             ! First bit used by time levels
             iSubCycleBlock = 2*iCrit
             ! Number of time levels is log_2(DtMax/DtMin)
             iTypeMax = iTypeMax &
                  + iSubCycleBlock*nint(alog(DtMax/DtMin) / alog(2.0))

          end if

          do iBlock = 1, nBlock

             if(iTypeAdvance_B(iBlock) == SkippedBlock_) CYCLE

             ! used block has 1st bit set
             iType = iNotSkippedBlock

             ! true block has second bit set
             if(True_Blk(iBlock)) iType = iType + iTrueBlock

             if(UsePartImplicit)then
                if(iTypeAdvance_B(iBlock)==ImplBlock_) &
                     iType = iType + iPartImplBlock
             end if

             if(UseSemiImplicit)then
                if(IsSemiImplBlock_B(iBlock)) iType = iType + iSemiImplBlock
             end if

             if(UseFieldLineThreads)then
                ! Field line threads are at the 1st boundary (minimum r)
                if(far_field_BCs_Blk(iBlock) .and. NeiLev(1,iBlock)==NOBLK)&
                     iType = iType + iFieldLineThreadBlock
             end if

             if(UsePic .and. DoBalancePicBlock)then
                if(IsPicNode_A(iNode_B(iBlock))) iType = iType + iPicBlock
             endif

             if(UsePartSteady)then
                if(iTypeAdvance_B(iBlock)==SteadyBlock_) &
                     iType = iType + iSteadyBlock
             end if

             if(UseLowOrderRegion)then
                ! High-order scheme block: some faces use high-order scheme
                if(.not.IsLowOrderOnly_B(iBlock)) &
                     iType = iType + iHighOrderBlock
             end if

             if(UseMaxTimeStep .and. DtMax > DtMin .and. DtMin > 0)then
                iType = iType + iSubCycleBlock &
                     *nint(alog(dt_BLK(iBlock)/DtMin) / alog(2.0))

                if(iType > iTypeMax)then
                   write(*,*) NameSub,' ERROR for iBlock, iProc=', iBlock, iProc
                   write(*,*) NameSub,'iType, iTypeMax =', iType, iTypeMax
                   write(*,*) NameSub,' iSubCycleBlock =', iSubCycleBlock
                   write(*,*) NameSub,' DtMin, DtMax   =', DtMin, DtMax
                   write(*,*) NameSub,' dt_BLK         =', dt_BLK(iBlock)
                   write(*,*) NameSub,' time level     =', &
                        nint(alog(dt_BLK(iBlock)/DtMin) / alog(2.0))
                end if
             end if

             iTypeBalance_A(iNode_B(iBlock)) = iType
          enddo
          ! Update iTypeAdvance_BP
          ! CHEATING: only indicate first index to circumvent ModMpiInterfaces
          ! Set displacement equal to MaxBlock so we get same behavior 
          ! as MPI_allgather. Use nBlockMax for maximum receive data for speed
          nBlockMax_P = nBlockMax
          call MPI_allgatherv(iTypeAdvance_B(1), nBlockMax, MPI_INTEGER, &
               iTypeAdvance_BP(1,0), nBlockMax_P, MaxBlockDisp_P,&
               MPI_INTEGER, iComm, iError)

          ! In order to do MPI_allreduce with MPI_SUM, the default value of 
          ! iTypeBalance_A=0. The skipped blocks also have 0 type value.
          call MPI_allreduce(MPI_IN_PLACE, iTypeBalance_A, MaxNode, &
               MPI_INTEGER, MPI_SUM, iComm, iError)

          ! Find all different block types that occur (ignore skipped blocks of type 0)
          if(.not.allocated(iType_I))       allocate(iType_I(0:iTypeMax))
          if(.not.allocated(IsTypeExist_I)) allocate(IsTypeExist_I(0:iTypeMax))     
          IsTypeExist_I = .false.
          do iNode = 1, nNode
             if(iTypeBalance_A(iNode) >= 0)&
                  IsTypeExist_I(iTypeBalance_A(iNode)) = .true.        
          enddo

          ! Count the number of different block types that occur and 
          ! create an indirect index from the full range of possible types.
          iType_I = -1
          nType = 0
          do iType = 1, iTypeMax
             if(IsTypeExist_I(iType)) then
                nType = nType + 1
                iType_I(iType) = nType
             endif
          enddo
       end if

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
               DoBalanceEachLevelIn= &
               (UseLocalTimeStep .and. .not.UseMaxTimeStep), &
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

       if(DoTestMe .and. nType > 0) then
          ! Create the table so that we can check the status of a block.
          if(.not.allocated(iTypeConvert_I)) allocate(iTypeConvert_I(nType))
          nCount = 1
          do iType = 1, iTypeMax
             if(iType_I(iType) > 0) then
                iTypeConvert_I(nCount) = iType
                nCount = nCount + 1           
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
             write(*,*) 'iType= ',iTypeConvert_I(i), &
                  ' nCount= ',sum(nTypeProc_PI(:,i))
          enddo
          deallocate(nTypeProc_PI)
       endif ! DoTestMe

       deallocate(iTypeAdvance_A)
    end if

    ! When load balancing is done Skipped and Unused blocks coincide
    Unused_BP(1:nBlockMax,:) = &
         iTypeAdvance_BP(1:nBlockMax,:) == SkippedBlock_

    call find_test_cell

    if(allocated(iType_I))       deallocate(iType_I)
    if(allocated(IsTypeExist_I)) deallocate(IsTypeExist_I)
    call timing_stop(NameSub)

    ! Allow the user to do something after load balancing was done
    if (nProc>1 .and. index(test_string,'NOLOADBALANCE') < 1 .and. DoMoveData)&
         call user_action('load balance done')

  end subroutine load_balance
  !============================================================================
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

    elseif((UsePartImplicit .and. .not. DoPartSelect) &
         .or. .not. UseImplicit)then
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
             if(n_step==1 .and. time_loop)then
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
  !==========================================================================
  subroutine load_balance_blocks

    use ModProcMH
    use ModImplicit, ONLY : UsePartImplicit, nBlockSemi, IsDynamicSemiImpl
    use ModPartSteady, ONLY: UsePartSteady, IsNewSteadySelect
    use ModTimeStepControl, ONLY: UseMaxTimeStep

    !LOCAL VARIABLES:
    logical:: DoBalanceSemiImpl = .true.

    character(len=*), parameter :: NameSub = 'load_balance_blocks'
    !--------------------------------------------------------------------------

    ! Select and load balance blocks
    if(  UseMaxTimeStep .or. &                         ! subcycling scheme
         UsePartImplicit .or. &                        ! part implicit scheme
         UsePartSteady .and. IsNewSteadySelect .or. &  ! part steady scheme
         nBlockSemi >= 0 .and. DoBalanceSemiImpl) then ! semi-implicit scheme

       ! Redo load balancing
       call load_balance(DoMoveCoord=.true., DoMoveData=.true., &
            IsNewBlock=.false.)

       IsNewSteadySelect = .false.

       ! Repeated semi implicit load balancing is only needed if the
       ! semi-implicit condition is changing dynamically. 
       DoBalanceSemiImpl = IsDynamicSemiImpl
    end if

  end subroutine load_balance_blocks

end module ModLoadBalance
!=============================================================================

