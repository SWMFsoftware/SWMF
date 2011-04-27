!^CFG COPYRIGHT UM

module BATL_amr

  implicit none

  SAVE

  private ! except

  public init_amr
  public do_amr
  public test_amr

  ! Parameter of slope limiter used by prolongation
  real, public:: BetaProlong = 1.5

  ! Status change due to AMR is registered in this 
  integer, public, allocatable:: iAmrChange_B(:)
  integer, public, parameter  :: AmrRemoved_ = -1, &
       AmrUnchanged_ = 0, AmrMoved_ = 1, AmrRefined_ = 2, AmrCoarsened_ = 3

contains

  !===========================================================================

  subroutine init_amr()

    use BATL_size, ONLY: MaxBlock
    !-------------------------------------------------------------------------

    if(.not.allocated(iAmrChange_B)) &
         allocate(iAmrChange_B(MaxBlock))

    ! Initialize iAmrChange_B
    iAmrChange_B = AmrUnchanged_

  end subroutine init_amr

  !===========================================================================

  subroutine do_amr(nVar, State_VGB, Dt_B, Used_GB, DoTestIn)

    use BATL_size, ONLY: MaxBlock, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         nI, nJ, nK, nIJK, iRatio, jRatio, kRatio
    use BATL_mpi,  ONLY: iComm, nProc, iProc

    use BATL_tree, ONLY: nNode, Unused_BP, &
         iTree_IA, iProcNew_A, Proc_, Block_, Coord1_, Coord2_, Coord3_, &
         Status_, Child1_, ChildLast_, &
         Used_, Unused_, Refine_, Refined_, CoarsenNew_, Coarsened_

    use BATL_geometry, ONLY: IsCartesian
    use BATL_grid, ONLY: create_grid_block, CellVolume_GB

    use ModMpi

    ! Arguments

    ! Number of variables
    integer, intent(in) :: nVar

    ! State variable to be refined/coarsened
    real, intent(inout) :: &
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    ! Time step limit for each block
    real, intent(inout), optional :: Dt_B(MaxBlock)

    ! Cells that can be used for AMR 
    ! For example cells inside internal boundaries cannot be used in BATSRUS.
    logical, intent(in), optional:: &
         Used_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    ! Write debug info if true
    logical, optional:: DoTestIn

    ! Dynamic arrays
    real,    allocatable :: Buffer_I(:), StateP_VG(:,:,:,:)
    real,    allocatable :: SlopeL_V(:), SlopeR_V(:), Slope_V(:)

    ! Permanently allocated array
    integer, save, allocatable :: iBlockAvailable_P(:)

    integer :: iNodeSend, iNodeRecv
    integer :: iProcSend, iProcRecv, iBlockSend, iBlockRecv
    integer :: iChild

    integer:: iStatus_I(MPI_STATUS_SIZE), iError

    integer, parameter:: iMinP = 2-iRatio, iMaxP = nI/iRatio + iRatio - 1
    integer, parameter:: jMinP = 2-jRatio, jMaxP = nJ/jRatio + jRatio - 1
    integer, parameter:: kMinP = 2-kRatio, kMaxP = nK/kRatio + kRatio - 1
    integer, parameter:: nSizeP = &
         (iMaxP-iMinP+1)*(jMaxP-jMinP+1)*(kMaxP-kMinP+1)

    integer, parameter:: MaxTry=100
    integer:: iTry
    logical:: DoTryAgain
    logical:: UseMask
    integer:: nVarBuffer

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'BATL_AMR::do_amr'
    !-------------------------------------------------------------------------
    DoTest = .false.
    if(present(DoTestIn)) DoTest = DoTestIn

    if(DoTest)write(*,*) NameSub,' starting'

    ! Simple logical
    UseMask = present(Used_GB)    

    ! If UseMask is true we will have one extra variable in StatP_VG
    ! containing a "boolean" represented as .true. = 1.0 and .false. = 0.0
    ! The masked cells cannot be used for 2nd order prolongation
    if(UseMask) then
       nVarBuffer = nVar+1
    else
       nVarBuffer = nVar
    end if

    ! Small arrays are allocated once 
    if(.not.allocated(iBlockAvailable_P)) &
         allocate(iBlockAvailable_P(0:nProc-1))
    !Moved to init_amr
    !if(.not.allocated(iAmrChange_B)) &
    !      allocate(iAmrChange_B(MaxBlock))


    ! Initialize iAmrChange_B
    !iAmrChange_B = AmrUnchanged_

    ! nVar dependent arrays are allocated and deallocated every call
    ! Buffer_I +1 for Dt_B and +1 for DoCheckMask header information
    allocate(Buffer_I(nVarBuffer*nIJK+2), &
         StateP_VG(nVarBuffer,iMinP:iMaxP,jMinP:jMaxP,kMinP:kMaxP), &
         SlopeL_V(nVar), SlopeR_V(nVar), Slope_V(nVar))

    ! Set iBlockAvailable_P to first available block
    iBlockAvailable_P = -1
    do iProcRecv = 0, nProc-1
       do iBlockRecv = 1, MaxBlock
          if(Unused_BP(iBlockRecv, iProcRecv))then
             iBlockAvailable_P(iProcRecv) = iBlockRecv
             EXIT
          end if
       end do
    end do

    LOOPTRY: do iTry = 1, MaxTry
       DoTryAgain = .false.

       ! Coarsen and move blocks
       do iNodeRecv = 1, nNode
          if(iTree_IA(Status_,iNodeRecv) /= CoarsenNew_) CYCLE

          if(DoTest) write(*,*)NameSub,' CoarsenNew iNode=',iNodeRecv

          iProcRecv  = iProcNew_A(iNodeRecv)
          if(iBlockAvailable_P(iProcRecv) > MaxBlock)then
             ! Continue with other nodes and then try again
             DoTryAgain = .true.
             CYCLE
          end if

          iBlockRecv = i_block_available(iProcRecv, iNodeRecv, AmrCoarsened_)

          do iChild = Child1_, ChildLast_
             iNodeSend = iTree_IA(iChild,iNodeRecv)

             iProcSend  = iTree_IA(Proc_,iNodeSend)
             iBlockSend = iTree_IA(Block_,iNodeSend)

             if(iProc == iProcSend) call send_coarsened_block
             if(iProc == iProcRecv) call recv_coarsened_block

             call make_block_available(iNodeSend, iBlockSend, iProcSend)
          end do

          ! This parent block was successfully coarsened
          iTree_IA(Status_,iNodeRecv) = Coarsened_
       end do

       ! Move blocks
       do iNodeSend = 1, nNode

          if(iTree_IA(Status_,iNodeSend) /= Used_) CYCLE

          iProcSend = iTree_IA(Proc_,iNodeSend)
          iProcRecv = iProcNew_A(iNodeSend)

          if(iProcRecv == iProcSend) CYCLE

          iBlockSend = iTree_IA(Block_,iNodeSend)
          if(iBlockAvailable_P(iProcRecv) > MaxBlock)then
             ! Continue with other nodes and then try again
             DoTryAgain = .true.
             CYCLE
          end if

          iBlockRecv = i_block_available(iProcRecv, iNodeSend, AmrMoved_)

          if(DoTest) write(*,*)NameSub, &
               ' node to move iNode,iProcS/R,iBlockS/R=',&
               iNodeSend, iProcSend, iProcRecv, iBlockSend, iBlockRecv

          if(iProc == iProcSend) call send_block
          if(iProc == iProcRecv) call recv_block

          call make_block_available(iNodeSend, iBlockSend, iProcSend)
       end do

       ! Prolong and move blocks
       LOOPNODE: do iNodeSend = 1, nNode

          if(iTree_IA(Status_,iNodeSend) /= Refine_) CYCLE

          iProcSend  = iTree_IA(Proc_,iNodeSend)
          iBlockSend = iTree_IA(Block_,iNodeSend)

          if(DoTest) write(*,*)NameSub,' Refine iNode=',iNodeSend

          do iChild = Child1_, ChildLast_
             iNodeRecv = iTree_IA(iChild,iNodeSend)

             ! Check if this child has been created already
             if(iTree_IA(Status_,iNodeRecv) == Refined_) CYCLE

             ! Check if there is a free block available
             iProcRecv  = iProcNew_A(iNodeRecv)

             if(iBlockAvailable_P(iProcRecv) > MaxBlock)then
                ! Continue with other nodes and then try again
                DoTryAgain = .true.
                CYCLE LOOPNODE
             end if

             iBlockRecv = i_block_available(iProcRecv, iNodeRecv, AmrRefined_)

             if(iProc == iProcSend) call send_refined_block
             if(iProc == iProcRecv) call recv_refined_block

             ! This child block was successfully refined
             iTree_IA(Status_,iNodeRecv) = Refined_
          end do
          call make_block_available(iNodeSend, iBlockSend, iProcSend)

       end do LOOPNODE

       !if(iProc==0)then
       !   write(*,*)'!!! Number of tries=', iTry
       !   write(*,*)'!!! iBlockAvailable_P=', iBlockAvailable_P
       !   do jProc = 0, nProc-1
       !      write(*,*)'jProc, count(Unused_BP)=', &
       !           jProc, count(Unused_BP(:,jProc))
       !   end do
       !end if

       if(.not.DoTryAgain) EXIT LOOPTRY

    end do LOOPTRY

    if(DoTryAgain) call CON_stop(NameSub//': could not fit blocks')

    deallocate(Buffer_I, StateP_VG, SlopeL_V, SlopeR_V, Slope_V)

  contains

    !==========================================================================
    integer function i_block_available(iProcRecv, iNodeRecv, iAmrChange)

      integer, intent(in):: iProcRecv, iNodeRecv, iAmrChange
      integer :: iBlock
      character(len=*), parameter:: NameSub = 'BATL_amr::i_block_available'
      !-----------------------------------------------------------------------
      ! Assign the processor index
      iTree_IA(Proc_,iNodeRecv) = iProcRecv

      ! Find and assign the block index
      iBlock = iBlockAvailable_P(iProcRecv)
      iTree_IA(Block_,iNodeRecv) = iBlock

      ! Return the block index
      i_block_available = iBlock

      ! Make the new block used
      Unused_BP(iBlock,iProcRecv) = .false.

      ! Create the grid info (x,y,z,volume,faces)
      if(iProc == iProcRecv)then
         iAmrChange_B(iBlock) = iAmrChange
         call create_grid_block(iBlock, iNodeRecv)
         State_VGB(:,:,:,:,iBlock) = 777.0
         if(present(Dt_B))  Dt_B(iBlock) = Huge(1.0)
      end if

      ! Find next available block
      do
         iBlock = iBlock + 1
         iBlockAvailable_P(iProcRecv) = iBlock
         if(iBlock > MaxBlock) RETURN
         if(Unused_BP(iBlock,iProcRecv)) RETURN
      end do

    end function i_block_available
    !==========================================================================
    subroutine make_block_available(iNodeSend, iBlockSend, iProcSend)

      integer, intent(in):: iNodeSend, iBlockSend, iProcSend
      !-----------------------------------------------------------------------

      iTree_IA(Status_,iNodeSend) = Unused_
      Unused_BP(iBlockSend,iProcSend) = .true.
      iBlockAvailable_P(iProcSend) = &
           min(iBlockAvailable_P(iProcSend), iBlockSend)

      if(iProc == iProcSend) iAmrChange_B(iBlockSend) = AmrRemoved_

    end subroutine make_block_available
    !==========================================================================
    subroutine send_block

      ! Copy buffer into recv block of State_VGB

      integer:: iBuffer, i, j, k
      !------------------------------------------------------------------------

      iBuffer = 0
      do k = 1, nK; do j = 1, nJ; do i = 1, nI
         Buffer_I(iBuffer+1:iBuffer+nVar) = State_VGB(:,i,j,k,iBlockSend)
         iBuffer = iBuffer + nVar
      end do; end do; end do

      if(present(Dt_B))then
         iBuffer = iBuffer + 1
         Buffer_I(iBuffer) = Dt_B(iBlockSend)
      end if

      ! Put more things into buffer here if necessary !!!
      call MPI_send(Buffer_I, iBuffer, MPI_REAL, iProcRecv, 1, iComm, iError)

    end subroutine send_block
    !==========================================================================
    subroutine recv_block

      ! Copy buffer into recv block of State_VGB

      integer:: iBuffer, i, j, k
      !------------------------------------------------------------------------

      iBuffer = nIJK*nVar
      if(present(Dt_B))  iBuffer = iBuffer + 1
      call MPI_recv(Buffer_I, iBuffer, MPI_REAL, iProcSend, 1, iComm, &
           iStatus_I, iError)

      iBuffer = 0
      do k = 1, nK; do j = 1, nJ; do i = 1, nI
         State_VGB(:,i,j,k,iBlockRecv) = Buffer_I(iBuffer+1:iBuffer+nVar)
         iBuffer = iBuffer + nVar
      end do; end do; end do

      if(present(Dt_B)) Dt_B(iBlockRecv) = Buffer_I(iBuffer+1)

    end subroutine recv_block

    !=========================================================================
    subroutine send_coarsened_block

      use BATL_size, ONLY:  InvIjkRatio

      integer :: i, j, k, i2, j2, k2, iVar, iBuffer
      logical:: DoCheckMask
      real :: CellVolume_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK), Volume, InvVolume
      !-----------------------------------------------------------------------
      iBuffer = 0

      ! Do we need to check the mask?
      if(UseMask)then
         DoCheckMask = .not. all(Used_GB(:,:,:,iBlockSend))
      else
         DoCheckMask = .false.
      end if


      ! Averaging all used fine cells inside the coarse cell. 
      if(DoCheckMask) then
         ! Set the volume of unused cells to zero
         where(Used_GB(:,:,:,iBlockSend)) 
            CellVolume_G = CellVolume_GB(:,:,:,iBlockSend)
         elsewhere
            CellVolume_G = 0.0
         endwhere

         do k = 1, nK, kRatio; do j = 1, nJ, jRatio; do i=1, nI, iRatio
            i2 = i+iRatio-1; j2 = j+jRatio-1; k2 = k+kRatio-1
            Volume = sum(CellVolume_G(i:i2,j:j2,k:k2))
            ! If there are any used fine cells, the total volume is positive
            if(Volume > 0.0)then
               ! Take volume weighted average excluding the unused cells
               InvVolume = 1.0/Volume
               do iVar = 1, nVar
                  Buffer_I(iBuffer+iVar) = InvVolume * sum( &
                       CellVolume_G(i:i2,j:j2,k:k2)* &
                       State_VGB(iVar,i:i2,j:j2,k:k2,iBlockSend))
               end do
            else
               ! None of the cells are used, use true cell volumes
               InvVolume = 1.0/sum(CellVolume_GB(i:i2,j:j2,k:k2,iBlockSend))
               do iVar = 1, nVar
                  Buffer_I(iBuffer+iVar) = InvVolume * sum( &
                       CellVolume_GB(i:i2,j:j2,k:k2,iBlockSend)* &
                       State_VGB(iVar,i:i2,j:j2,k:k2,iBlockSend))
               end do
            end if
            iBuffer = iBuffer + nVar
         end do; end do; end do
      else
         if(IsCartesian)then
            do k = 1, nK, kRatio; do j = 1, nJ, jRatio; do i=1, nI, iRatio
               do iVar = 1, nVar
                  Buffer_I(iBuffer+iVar) = InvIjkRatio * &
                       sum(State_VGB(iVar,i:i+iRatio-1,j:j+jRatio-1,&
                       k:k+kRatio-1, iBlockSend))
               end do
               iBuffer = iBuffer + nVar
            end do; end do; end do
         else
            do k = 1, nK, kRatio; do j = 1, nJ, jRatio; do i=1, nI, iRatio
               i2 = i+iRatio-1; j2 = j+jRatio-1; k2 = k+kRatio-1
               InvVolume = 1.0/sum(CellVolume_GB(i:i2,j:j2,k:k2,iBlockSend))
               do iVar = 1, nVar
                  Buffer_I(iBuffer+iVar) = InvVolume * sum( &
                       CellVolume_GB(i:i2,j:j2,k:k2,iBlockSend)* &
                       State_VGB(iVar,i:i2,j:j2,k:k2,iBlockSend))
               end do
               iBuffer = iBuffer + nVar
            end do; end do; end do
         end if
      end if

      if(present(Dt_B))then
         iBuffer = iBuffer + 1
         Buffer_I(iBuffer) = Dt_B(iBlockSend)
      end if

      if(iProcRecv /= iProcSend) call MPI_send(Buffer_I, iBuffer, &
           MPI_REAL, iProcRecv, 1, iComm, iError)

    end subroutine send_coarsened_block
    !==========================================================================
    subroutine recv_coarsened_block

      use BATL_size, ONLY: IjkRatio

      integer:: iBuffer, iSide, jSide, kSide
      integer:: iMin, jMin, kMin, iMax, jMax, kMax
      integer:: i, j, k
      !----------------------------------------------------------------------

      if(iProcRecv /= iProcSend)then
         iBuffer = nIJK*nVar/IjkRatio 
         if(present(Dt_B))  iBuffer = iBuffer + 1
         call MPI_recv(Buffer_I, iBuffer, MPI_REAL, iProcSend, 1, iComm, &
              iStatus_I, iError)
      end if

      ! Find the part of the block to be written into
      iSide = modulo(iTree_IA(Coord1_,iNodeSend)-1, iRatio)
      jSide = modulo(iTree_IA(Coord2_,iNodeSend)-1, jRatio)
      kSide = modulo(iTree_IA(Coord3_,iNodeSend)-1, kRatio)

      iMin = 1 + iSide*nI/2; iMax = iMin + nI/iRatio - 1
      jMin = 1 + jSide*nJ/2; jMax = jMin + nJ/jRatio - 1
      kMin = 1 + kSide*nK/2; kMax = kMin + nK/kRatio - 1

      iBuffer = 0
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         State_VGB(:,i,j,k,iBlockRecv) = Buffer_I(iBuffer+1:iBuffer+nVar)
         iBuffer = iBuffer + nVar
      end do; end do; end do


      ! Take the smallest of the doubled (valid for full AMR only) time steps 
      ! of the children blocks
      if(present(Dt_B)) &
           Dt_B(iBlockRecv) = min(Dt_B(iBlockRecv), 2*Buffer_I(iBuffer+1))

    end subroutine recv_coarsened_block
    !==========================================================================
    subroutine send_refined_block

      use BATL_size, ONLY: nDim

      integer:: iSide, jSide, kSide
      integer:: iMin, jMin, kMin, iMax, jMax, kMax
      integer:: i, j, k, iBuffer, iVar
      real   :: CellVolume_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK), Volume
      logical:: DoCheckMask
      !----------------------------------------------------------------------

      ! Find the part of the block to be prolonged
      iSide = modulo(iTree_IA(Coord1_,iNodeRecv)-1, iRatio)
      jSide = modulo(iTree_IA(Coord2_,iNodeRecv)-1, jRatio)
      kSide = modulo(iTree_IA(Coord3_,iNodeRecv)-1, kRatio)

      ! Send parent part of the block with one ghost cell
      if(iRatio == 2)then
         iMin = iSide*nI/2; iMax = iMin + nI/2 + 1
      else
         iMin = 1; iMax = nI
      endif
      if(jRatio == 2)then
         jMin = jSide*nJ/2; jMax = jMin + nJ/2 + 1
      else
         jMin = 1; jMax = nJ
      end if
      if(kRatio == 2)then
         kMin = kSide*nK/2; kMax = kMin + nK/2 + 1
      else
         kMin = 1; kMax = nK
      end if

      if(UseMask)then
         DoCheckMask = &
              .not. all(Used_GB(iMin:iMax,jMin:jMax,kMin:kMax,iBlockSend))
      else
         DoCheckMask = .false.
      end if

      iBuffer = 1
      Buffer_I(iBuffer) = logical_to_real(DoCheckMask)

      if(DoCheckMask) then
         ! Set the volume of unused cells to zero
         where(Used_GB(:,:,:,iBlockSend))
            CellVolume_G = CellVolume_GB(:,:,:,iBlockSend)
         elsewhere
            CellVolume_G = 0.0
         endwhere

         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax

            if(.not.Used_GB(i,j,k,iBlockSend))then

               ! Find the total volume of used neighbor cells 
               Volume = sum(CellVolume_G(i-1:i+1,j,k))
               if(nDim >1) Volume = Volume + sum(CellVolume_G(i,j-1:j+1,k))
               if(nDim >2) Volume = Volume + sum(CellVolume_G(i,j,k-1:k+1))

               if(Volume > 0.0) then
                  do iVar=1,nVar
                     Buffer_I(iBuffer+iVar) = &
                          sum(State_VGB(iVar,i-1:i+1,j,k,iBlockSend)* &
                          CellVolume_G(i-1:i+1,j,k))
                     if(nDim>1) Buffer_I(iBuffer+iVar) = &
                          Buffer_I(iBuffer+iVar) + &
                          sum(State_VGB(iVar,i,j-1:j+1,k,iBlockSend)*&
                          CellVolume_G(i,j-1:j+1,k))
                     if(nDim>2) Buffer_I(iBuffer+iVar) = &
                          Buffer_I(iBuffer+iVar) + &
                          sum(State_VGB(iVar,i,j,k-1:k+1,iBlockSend)*&
                          CellVolume_G(i,j,k-1:k+1))

                     Buffer_I(iBuffer+iVar) = Buffer_I(iBuffer+iVar)/Volume

                     if(.not.IsCartesian) &
                          Buffer_I(iBuffer+iVar) = Buffer_I(iBuffer+iVar) &
                          *CellVolume_GB(i,j,k,iBlockSend)

                  end do
                  ! Send the information about the masked cells to
                  ! the prolongation operation converting 
                  ! .true. = 1.0 and .false. = 0.0
                  Buffer_I(iBuffer+nVar+1) = &
                       logical_to_real(Used_GB(i,j,k,iBlockSend))
                  iBuffer = iBuffer + nVar+1
                  CYCLE
               end if
            end if
            if(IsCartesian)then
               Buffer_I(iBuffer+1:iBuffer+nVar) = &
                    State_VGB(:,i,j,k,iBlockSend)
            else
               Buffer_I(iBuffer+1:iBuffer+nVar) = &
                    State_VGB(:,i,j,k,iBlockSend) * &
                    CellVolume_GB(i,j,k,iBlockSend)
            end if
            Buffer_I(iBuffer+nVar+1) = &
                 logical_to_real(Used_GB(i,j,k,iBlockSend))
            iBuffer = iBuffer + nVar+1
         end do; end do; end  do
      else
         iBuffer = 1
         if(IsCartesian)then
            do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
               Buffer_I(iBuffer+1:iBuffer+nVar) = &
                    State_VGB(:,i,j,k,iBlockSend)
               iBuffer = iBuffer + nVar
            enddo; end do; end do
         else
            do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
               Buffer_I(iBuffer+1:iBuffer+nVar) = &
                    State_VGB(:,i,j,k,iBlockSend) &
                    *CellVolume_GB(i,j,k,iBlockSend)
               iBuffer = iBuffer + nVar
            enddo; end do; end do
         end if
      end if

      if(present(Dt_B))then
         iBuffer = iBuffer + 1
         Buffer_I(iBuffer) = Dt_B(iBlockSend)
      end if

      if(iProcRecv /= iProcSend) &
           call MPI_send(Buffer_I, iBuffer, MPI_REAL, iProcRecv, 1, &
           iComm, iError)

    end subroutine send_refined_block
    !==========================================================================
    real function logical_to_real(IsTrue)

      logical, intent(in):: IsTrue
      !-----------------------------------------------------------------------
      if(IsTrue)then
         logical_to_real = 1.0
      else
         logical_to_real = 0.0
      end if

    end function logical_to_real
    !==========================================================================
    subroutine recv_refined_block

      ! Copy buffer into recv block of State_VGB

      use BATL_size, ONLY: InvIjkRatio

      integer:: iBuffer, i, j, k
      integer:: nVarUsed
      integer:: iP, jP, kP, iR, jR, kR
      integer, parameter:: Di = iRatio-1, Dj = jRatio-1, Dk = kRatio-1

      logical:: DoCheckMask
      logical:: UseSlopeI, UseSlopeJ, UseSlopeK
      !------------------------------------------------------------------------

      UseSlopeI = .true.
      UseSlopeJ = .true.
      UseSlopeK = .true.

      if(iProcRecv /= iProcSend)then
         iBuffer = nSizeP*nVarBuffer+1
         if(present(Dt_B)) iBuffer = iBuffer + 1
         call MPI_recv(Buffer_I, iBuffer, MPI_REAL, iProcSend, 1, iComm, &
              iStatus_I, iError)
      end if

      DoCheckMask = Buffer_I(1) > 0.5
      !      print *,"DoCheckMask",Buffer_I(1),DoCheckMask, UseMask
      if(DoCheckMask) then
         nVarUsed = nVarBuffer 
      else 
         nVarUsed = nVar
      end if

      ! StateP_VG(nVar+1,:,:,:) is the Used_GB for the parent block
      iBuffer = 1
      do kP = kMinP, kMaxP; do jP = jMinP, jMaxP; do iP = iMinP, iMaxP
         StateP_VG(1:nVarUsed,iP,jP,kP) = Buffer_I(iBuffer+1:iBuffer+nVarUsed)
         iBuffer = iBuffer + nVarUsed
      end do; end do; end do
      ! Set time step to half of the parent block
      if(present(Dt_B)) Dt_B(iBlockRecv) = 0.5*Buffer_I(iBuffer+1)

      ! 1st order prolongation
      do kR = 1, nK
         kP = (kR + Dk)/kRatio
         do jR = 1, nJ
            jP = (jR + Dj)/jRatio
            do iR = 1, nI
               iP = (iR + Di)/iRatio
               State_VGB(:,iR,jR,kR,iBlockRecv) = StateP_VG(1:nVar,iP,jP,kP)
            end do
         end do
      end do


      if(BetaProlong > 0.0)then
         ! Add corection for 2nd order prolongation
         do kP = 1, nK/kRatio
            kR = kRatio*(kP - 1) + 1
            do jP = 1, nJ/jRatio
               jR = jRatio*(jP - 1) + 1
               do iP = 1, nI/iRatio
                  iR = iRatio*(iP - 1) + 1

                  ! If one of the neighboring cells or the cell are masked we
                  ! will only be able to use 1st order prolongation in that
                  ! dimension 

                  if(DoCheckMask) then
                     ! If any of the neighbor cells are masked the product
                     ! will be zero and no 2nd order prolongation can be done
                     if(iRatio == 2) UseSlopeI = &
                          product(StateP_VG(nVar+1,iP-1:iP+1,jP,kP)) > 0.5
                     if(jRatio == 2) UseSlopeJ = &
                          product(StateP_VG(nVar+1,iP,jP-1:jP+1,kP)) > 0.5
                     if(kRatio == 2) UseSlopeK = &
                          product(StateP_VG(nVar+1,iP,jP,kP-1:kP+1)) > 0.5
                  end if

                  if(iRatio == 2 .and. UseSlopeI)then

                     SlopeL_V = StateP_VG(1:nVar,iP  ,jP,kP)-&
                          StateP_VG(1:nVar,iP-1,jP,kP)
                     SlopeR_V = StateP_VG(1:nVar,iP+1,jP,kP)-&
                          StateP_VG(1:nVar,iP  ,jP,kP)
                     Slope_V  = &
                          (sign(0.125,SlopeL_V)+sign(0.125,SlopeR_V))*min( &
                          BetaProlong*abs(SlopeL_V), &
                          BetaProlong*abs(SlopeR_V), &
                          0.5*abs(SlopeL_V+SlopeR_V))
                     do k = kR, kR+Dk; do j = jR, jR+Dj
                        State_VGB(:,iR,j,k,iBlockRecv) = &
                             State_VGB(:,iR,j,k,iBlockRecv) - Slope_V
                        State_VGB(:,iR+1,j,k,iBlockRecv) = &
                             State_VGB(:,iR+1,j,k,iBlockRecv) + Slope_V
                     end do; end do

                  end if

                  if(jRatio == 2 .and. UseSlopeJ)then

                     SlopeL_V = StateP_VG(1:nVar,iP,jP  ,kP)-&
                          StateP_VG(1:nVar,iP,jP-1,kP)
                     SlopeR_V = StateP_VG(1:nVar,iP,jP+1,kP)-&
                          StateP_VG(1:nVar,iP,jP  ,kP)
                     Slope_V  = &
                          (sign(0.125,SlopeL_V)+sign(0.125,SlopeR_V))*min( &
                          BetaProlong*abs(SlopeL_V), &
                          BetaProlong*abs(SlopeR_V), &
                          0.5*abs(SlopeL_V+SlopeR_V))
                     do k = kR, kR+Dk; do i = iR, iR+Di
                        State_VGB(:,i,jR,k,iBlockRecv) = &
                             State_VGB(:,i,jR,k,iBlockRecv) - Slope_V
                        State_VGB(:,i,jR+1,k,iBlockRecv) = &
                             State_VGB(:,i,jR+1,k,iBlockRecv) + Slope_V
                     end do; end do

                  end if

                  if(kRatio == 2 .and. UseSlopeK)then

                     SlopeL_V = StateP_VG(1:nVar,iP,jP,kP  )-&
                          StateP_VG(1:nVar,iP,jP,kP-1)
                     SlopeR_V = StateP_VG(1:nVar,iP,jP,kP+1)-&
                          StateP_VG(1:nVar,iP,jP,kP)
                     Slope_V  = &
                          (sign(0.125,SlopeL_V)+sign(0.125,SlopeR_V))*min( &
                          BetaProlong*abs(SlopeL_V), &
                          BetaProlong*abs(SlopeR_V), &
                          0.5*abs(SlopeL_V+SlopeR_V))
                     do j = jR, jR+Dj; do i = iR, iR+Di
                        State_VGB(:,i,j,kR,iBlockRecv) = &
                             State_VGB(:,i,j,kR,iBlockRecv) - Slope_V
                        State_VGB(:,i,j,kR+1,iBlockRecv) = &
                             State_VGB(:,i,j,kR+1,iBlockRecv) + Slope_V
                     end do; end do

                  end if

               end do
            end do
         end do
      endif

      if(IsCartesian) RETURN

      ! Divide back by volume and IjkRatio
      do kR = 1, nK; do jR = 1, nJ; do iR = 1, nI
         State_VGB(:,iR,jR,kR,iBlockRecv) = State_VGB(:,iR,jR,kR,iBlockRecv) &
              * InvIjkRatio / CellVolume_GB(iR,jR,kR,iBlockRecv)
      end do; end do; end do

    end subroutine recv_refined_block

  end subroutine do_amr

  !============================================================================

  subroutine test_amr

    use BATL_mpi,  ONLY: iProc, nProc, barrier_mpi
    use BATL_size, ONLY: MaxDim, nDim, nIJK_D, iDimAmr_D, &
         MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nI, nJ, nK, nBlock, &
         iRatio, jRatio, kRatio
    use BATL_tree, ONLY: init_tree, set_tree_root, refine_tree_node, &
         coarsen_tree_node, distribute_tree, move_tree, show_tree, clean_tree,&
         iProcNew_A, Unused_B, nNode, iNode_B, iTree_IA, Proc_, Block_, &
         Child1_, unset_, nNode
    use BATL_grid, ONLY: init_grid, create_grid, clean_grid, Xyz_DGB, &
         CellSize_DB
    use BATL_geometry, ONLY: init_geometry

    integer, parameter:: MaxBlockTest            = 100
    integer, parameter:: nRootTest_D(MaxDim)     = (/3,3,3/)
    logical, parameter:: IsPeriodicTest_D(MaxDim)= .true.
    real, parameter:: DomainMin_D(MaxDim) = (/ 1.0, 10.0, 100.0 /)
    real, parameter:: DomainMax_D(MaxDim) = (/10.0,100.0,1000.0 /)
    real, parameter:: DomainSize_D(MaxDim) = DomainMax_D - DomainMin_D

    integer, parameter:: nVar = nDim
    real, allocatable:: State_VGB(:,:,:,:,:), Dt_B(:)
    real, allocatable:: TestState_VC(:,:,:,:)
    logical, allocatable:: Used_GB(:,:,:,:)
    integer, allocatable:: iEffectedNode_A(:)
    integer:: iBlock, iDim, iNode, iVar,i,j,k
    integer:: iChild

    logical:: DoTestMe
    character(len=*), parameter :: NameSub = 'test_amr'
    !-----------------------------------------------------------------------
    DoTestMe = iProc == 0

    if(DoTestMe) write(*,*) 'Starting ',NameSub

    call init_tree(MaxBlockTest)
    call init_grid( DomainMin_D(1:nDim), DomainMax_D(1:nDim) )
    call init_geometry( IsPeriodicIn_D = IsPeriodicTest_D(1:nDim) )
    call set_tree_root( nRootTest_D(1:nDim))
    call distribute_tree(.true.)
    call create_grid
    call init_amr

    if(DoTestMe) call show_tree('after create_grid')

    allocate(State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlockTest), &
         Dt_B(MaxBlockTest))

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       State_VGB(:,:,:,:,iBlock) = Xyz_DGB(1:nDim,:,:,:,iBlock)
       ! set the time step to the cells size in the first AMR direction
       iDim = iDimAmr_D(1)
       Dt_B(iBlock) = DomainSize_D(iDim) / (nIjk_D(iDim)*nRootTest_D(iDim))
    end do

    if(DoTestMe) write(*,*)'test prolong and balance'
    call refine_tree_node(1)
    if(DoTestMe) call show_tree('after refine_tree_node')
    call distribute_tree(.false.)
    if(DoTestMe)then
       call show_tree('after distribute_tree(.false.)')
       write(*,*)'iProc, iProcNew_A=',iProc, iProcNew_A(1:nNode)
    end if

    call do_amr(nVar, State_VGB, Dt_B)
    if(DoTestMe) call show_tree('after do_amr')
    call move_tree
    if(DoTestMe) call show_tree('after move_tree')
    if(DoTestMe) write(*,*)'iAmrChange_B=', iAmrChange_B(1:nBlock)

    call check_state

    if(DoTestMe) write(*,*)'test restrict and balance'
    call coarsen_tree_node(1)
    if(DoTestMe) call show_tree('after coarsen_tree_node')
    call distribute_tree(.false.)
    if(DoTestMe)then
       call show_tree('after distribute_tree(.false.)')
       write(*,*)'iProc, iProcNew_A=',iProc, iProcNew_A(1:nNode)
    end if
    call do_amr(nVar, State_VGB, Dt_B)
    if(DoTestMe) call show_tree('after do_amr')
    call move_tree
    if(DoTestMe) call show_tree('after move_tree')
    if(DoTestMe) write(*,*)'iAmrChange_B=', iAmrChange_B(1:nBlock)

    call check_state

    ! tests with mask
    if(DoTestMe) write(*,*) 'test masked cells'

    allocate(Used_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlockTest),&
         TestState_VC(nVar,nI,nJ,nK))

    Used_GB = .true.

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       State_VGB(:,:,:,:,iBlock)    = Xyz_DGB(1:nDim,:,:,:,iBlock)
       ! set the time step to the cells size in the first AMR direction
       iDim = iDimAmr_D(1)
       Dt_B(iBlock) = DomainSize_D(iDim) / (nIjk_D(iDim)*nRootTest_D(iDim))
    end do

    call refine_tree_node(1)
    call distribute_tree(.false.)

    ! Mask out the middle cell of node 1
    if(iTree_IA(Proc_,1) == iProc) then
       iBlock = iTree_IA(Block_,1)
       State_VGB(:,(MinI+MaxI)/2,(MinJ+MaxJ)/2,(MinK+MaxK)/2,iBlock) = -7777
       Used_GB((MinI+MaxI)/2,(MinJ+MaxJ)/2,(MinK+MaxK)/2,iBlock)     = .false.
    end if

    call do_amr(nVar, State_VGB, Dt_B, Used_GB=Used_GB)
    call move_tree

    ! Nodes that needs special care in the testing
    allocate(iEffectedNode_A(nNode)) 
    iEffectedNode_A = unset_

    ! 0: the masked cell, 1: x+ neighbor, 2: y+ neighbor
    ! 3: z+ neighbor
    ! The first child of Node of the refined node
    iChild = Child1_
    iEffectedNode_A(iTree_IA(iChild,1))                 = 0
    iChild = iChild + 1
    ! As masked cell is on the corner of the childe nodes
    ! it will also effect the neighboring cell in i-direction
    ! The running "iChild" cunter is using the morton indexing
    ! to give the right neighbor independent of which directions
    ! are refined.
    if(iRatio == 2) then
       iEffectedNode_A(iTree_IA(iChild,1))              = 1
       iChild = iChild + 1
    end if
    ! Neoghbor in j direction
    if(jRatio == 2) then
       iEffectedNode_A(iTree_IA(iChild,1))              = 2
       iChild = iChild + iRatio
    end if
    ! and k direction
    if(kRatio == 2) iEffectedNode_A(iTree_IA(iChild,1)) = 3

    call check_state_mask

    !--------------- END refine --------------------
    iEffectedNode_A = unset_
    Used_GB = .true.

    if(iTree_IA(Proc_,iTree_IA(Child1_,1)) == iProc) then
       iBlock = iTree_IA(Block_,iTree_IA(Child1_,1))
       if(iRatio == 2) then
          i = nI
       else
          i = (MinI+MaxI)/2
       end if
       if(jRatio == 2) then
          j = nJ
       else
          j = (MinJ+MaxJ)/2
       end if
       if(kRatio == 2) then
          k = nK
       else
          k = (MinK+MaxK)/2
       end if

       State_VGB(:,i,j,k,iBlock) = -7777
       Used_GB(i,j,k,iBlock)     = .false.
    end if

    call coarsen_tree_node(1)
    call distribute_tree(.false.)

    call do_amr(nVar, State_VGB, Dt_B,Used_GB=Used_GB)
    call move_tree

    call check_state_mask


    deallocate(State_VGB, Dt_B,Used_GB, iEffectedNode_A,TestState_VC)

    call clean_grid
    call clean_tree

  contains
    !==========================================================================
    subroutine check_state

      do iBlock = 1, nBlock
         if(Unused_B(iBlock)) CYCLE

         if(any(abs(State_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
              -     Xyz_DGB(1:nDim,1:nI,1:nJ,1:nK,iBlock)) > 1e-6))then
            write(*,*)NameSub,' error for iProc,iBlock,maxloc=',iProc,iBlock,&
                 maxloc(abs(State_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
                 -    Xyz_DGB(1:nDim,1:nI,1:nJ,1:nK,iBlock)))
         end if

         iDim = iDimAmr_D(1)
         if(abs(Dt_B(iBlock) - CellSize_DB(iDim,iBlock)) > 1e-6) &
              write(*,*)NameSub,' error for iProc,iBlock,dt,dx=', &
              iProc, iBlock, Dt_B(iBlock), CellSize_DB(iDim,iBlock)

      end do

    end subroutine check_state
    !==========================================================================
    subroutine check_state_mask

      integer :: iDn, iUp, jDn, jUp, kDn, kUp
      integer :: Di, Dj, Dk


      do iBlock = 1, nBlock
         if(Unused_B(iBlock)) CYCLE
         iNode = iNode_B(iBlock)
         select case(iEffectedNode_A(iNode))
         case(unset_)
            if(any(abs(State_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
                 -     Xyz_DGB(1:nDim,1:nI,1:nJ,1:nK,iBlock)) > 1e-6))then
               write(*,*)NameSub,' error for iProc,iBlock,maxloc=', &
                    iProc,iBlock,&
                    maxloc(abs(State_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
                    - Xyz_DGB(1:nDim,1:nI,1:nJ,1:nK,iBlock)))
            end if

            iDim = iDimAmr_D(1)
            if(abs(Dt_B(iBlock) - CellSize_DB(iDim,iBlock)) > 1e-6) &
                 write(*,*)NameSub,' error for iProc,iBlock,dt,dx=', &
                 iProc, iBlock, Dt_B(iBlock), CellSize_DB(iDim,iBlock)
         case(0) ! the masked block

            if(iRatio == 2) then
               iDn = nI-iRatio+1
               iUp = nI
               Di = iRatio
            else
               iDn = (MinI+MaxI)/2
               iUp = iDn
               Di = 0
            end if

            if(jRatio == 2) then
               jDn = nJ-jRatio+1 
               jUp = nJ
               Dj = jRatio
            else
               jDn = (MinJ+MaxJ)/2
               jUp = jDn
               Dj = 0
            end if

            if(kRatio == 2) then
               kDn = nK-kRatio+1
               kUp = nK
               Dk = kRatio
            else
               kDn = (MinK+MaxK)/2
               kUp = kDn
               Dk = 0
            end if

            TestState_VC = Xyz_DGB(1:nDim,1:nI,1:nJ,1:nK,iBlock)

            ! the average of the masked Cell
            do iVar=1,nVar
               TestState_VC(iVar,iDn:iUp,jDn:jUp,kDn:kUp) = &
                    sum(Xyz_DGB(iVar,iDn:iUp,jDn:jUp,kDn:kUp,iBlock))/&
                    (iRatio*jRatio*kRatio)
            end do

            ! The neighbor cell in i direction will only have 
            ! 1st order prolongation 
            TestState_VC(1,iDn-Di:iUp-Di,jDn:jUp,kDn:kUp) = &
                 sum(Xyz_DGB(1,iDn-Di:iUp-Di,jDn:jUp,kDn:kUp,iBlock))/&
                 (iRatio*jRatio*kRatio)                 

            ! The neighbor cell in j direction will only have 
            ! 1st order prolongation 
            if(nDim >1) then
               TestState_VC(2,iDn:iUp,jDn-Dj:jUp-Dj,kDn:kUp) = &
                    sum(Xyz_DGB(2,iDn:iUp,jDn-Dj:jUp-Dj,kDn:kUp,iBlock))/&
                    (iRatio*jRatio*kRatio)
            end if


            ! The neighbor cell in k direction will only have 
            ! 1st order prolongation 
            if(nDim >2) then
               TestState_VC(3,iDn:iUp,jDn:jUp,kDn-Dk:kUp-Dk) = &
                    sum(Xyz_DGB(3,iDn:iUp,jDn:jUp,kDn-Dk:kUp-Dk,iBlock))/&
                    (iRatio*jRatio*kRatio)
            end if

            if(any(abs(State_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
                 -     TestState_VC) > 1e-6))then
               write(*,*)NameSub,' case 0 error for iProc,iBlock,maxloc=',&
                    iProc,iBlock,&
                    maxloc(abs(State_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
                    -    TestState_VC))
            end if

            iDim = iDimAmr_D(1)
            if(abs(Dt_B(iBlock) - CellSize_DB(iDim,iBlock)) > 1e-6) &
                 write(*,*)NameSub,' error for iProc,iBlock,dt,dx=', &
                 iProc, iBlock, Dt_B(iBlock), CellSize_DB(iDim,iBlock)

         case(1) ! i neigbor of mask block

            if(iRatio == 2) then
               iDn = 1
               iUp = iRatio
               Di  = iRatio
            else
               iDn = (MinI+MaxI)/2
               iUp = iDn
               Di = 0
            end if

            if(jRatio == 2) then
               jDn = nJ-jRatio+1 
               jUp = nJ
               Dj = jRatio
            else
               jDn = (MinJ+MaxJ)/2
               jUp = jDn
               Dj = 0
            end if

            if(kRatio == 2) then
               kDn = nK-kRatio+1
               kUp = nK
               Dk = kRatio
            else
               kDn = (MinK+MaxK)/2
               kUp = kDn
               Dk = 0
            end if

            TestState_VC = Xyz_DGB(1:nDim,1:nI,1:nJ,1:nK,iBlock)

            ! The neighbor cell in i direction will only have 
            ! 1st order prolongation 
            TestState_VC(1,iDn:iUp,jDn:jUp,kDn:kUp) = &
                 sum(Xyz_DGB(1,iDn:iUp,jDn:jUp,kDn:kUp,iBlock))/&
                 (iRatio*jRatio*kRatio)

            if(any(abs(State_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
                 -     TestState_VC) > 1e-6))then
               write(*,*)NameSub,' case 1 error for iProc,iBlock,maxloc=',&
                    iProc,iBlock,&
                    maxloc(abs(State_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
                    -    TestState_VC))
            end if

            iDim = iDimAmr_D(1)
            if(abs(Dt_B(iBlock) - CellSize_DB(iDim,iBlock)) > 1e-6) &
                 write(*,*)NameSub,' error for iProc,iBlock,dt,dx=', &
                 iProc, iBlock, Dt_B(iBlock), CellSize_DB(iDim,iBlock)

         case(2) ! j neighbore of masked block

            if(jRatio == 2) then
               jDn = 1 
               jUp = jRatio
               Dj  = jRatio
            else
               jDn = (MinJ+MaxJ)/2
               jUp = jDn
               Dj  = 0
            end if

            if(iRatio == 2) then
               iDn = nI-iRatio+1
               iUp = nI
               Di = iRatio
            else
               iDn = (MinI+MaxI)/2
               iUp = iDn
               Di = 0
            end if

            if(kRatio == 2) then
               kDn = nK-kRatio+1
               kUp = nK
               Dk = kRatio
            else
               kDn = (MinK+MaxK)/2
               kUp = kDn
               Dk = 0
            end if

            TestState_VC = Xyz_DGB(1:nDim,1:nI,1:nJ,1:nK,iBlock)

            ! The neighbor cell in j direction will only have 
            ! 1st order prolongation 
            TestState_VC(2,iDn:iUp,jDn:jUp,kDn:kUp) = &
                 sum(Xyz_DGB(2,iDn:iUp,jDn:jUp,kDn:kUp,iBlock))/&
                 (iRatio*jRatio*kRatio)

            if(any(abs(State_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
                 -     TestState_VC) > 1e-6))then
               write(*,*)NameSub,' case 2 error for iProc,iBlock,maxloc=',&
                    iProc,iBlock,&
                    maxloc(abs(State_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
                    -    TestState_VC))
            end if

            iDim = iDimAmr_D(1)
            if(abs(Dt_B(iBlock) - CellSize_DB(iDim,iBlock)) > 1e-6) &
                 write(*,*)NameSub,' error for iProc,iBlock,dt,dx=', &
                 iProc, iBlock, Dt_B(iBlock), CellSize_DB(iDim,iBlock)

         case(3) ! k neighbore of masked block

            if(kRatio == 2) then
               kDn = 1
               kUp = kRatio
               Dk = kRatio
            else
               kDn = (MinK+MaxK)/2
               kUp = kDn
               Dk = 0
            end if

            if(iRatio == 2) then
               iDn = nI-iRatio+1
               iUp = nI
               Di = iRatio
            else
               iDn = (MinI+MaxI)/2
               iUp = iDn
               Di = 0
            end if

            if(jRatio == 2) then
               jDn = nJ-jRatio+1 
               jUp = nJ
               Dj = jRatio
            else
               jDn = (MinJ+MaxJ)/2
               jUp = jDn
               Dj = 0
            end if

            TestState_VC = Xyz_DGB(1:nDim,1:nI,1:nJ,1:nK,iBlock)

            ! The neighbor cell in k direction will only have 
            ! 1st order prolongation 
            TestState_VC(3,iDn:iUp,jDn:jUp,kDn:kUp) = &
                 sum(Xyz_DGB(3,iDn:iUp,jDn:jUp,kDn:kUp,iBlock))/&
                 (iRatio*jRatio*kRatio)

            if(any(abs(State_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
                 -     TestState_VC) > 1e-6))then
               write(*,*)NameSub,' case 3 error for iProc,iBlock,maxloc=',&
                    iProc,iBlock,&
                    maxloc(abs(State_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
                    -    TestState_VC))
            end if

            iDim = iDimAmr_D(1)
            if(abs(Dt_B(iBlock) - CellSize_DB(iDim,iBlock)) > 1e-6) &
                 write(*,*)NameSub,' error for iProc,iBlock,dt,dx=', &
                 iProc, iBlock, Dt_B(iBlock), CellSize_DB(iDim,iBlock)

         end select
      end do

    end subroutine check_state_mask
    !==========================================================================
    subroutine show_state

      integer :: iProcShow
      !-----------------------------------------------------------------------
      call barrier_mpi

      do iProcShow = 0, nProc - 1
         if(iProc == iProcShow)then
            do iBlock = 1, nBlock
               if(Unused_B(iBlock)) CYCLE
               write(*,'(a,2i4,100f8.4)') &
                    'iProc, iBlock, State(1,1:nI,1,1,iBlock)=', &
                    iProc, iBlock, State_VGB(1,1:nI,1,1,iBlock)

               if(nDim > 1) write(*,'(a,2i4,100f8.3)') &
                    'iProc, iBlock, State(2,1,1:nJ,1,iBlock)=', &
                    iProc, iBlock, State_VGB(2,1,1:nJ,1,iBlock)

               if(nDim > 2) write(*,'(a,2i4,100f8.2)') &
                    'iProc, iBlock, State(3,1,1,1:nK,iBlock)=', &
                    iProc, iBlock, State_VGB(3,1,1,1:nK,iBlock)
            end do
         end if
         call barrier_mpi
      end do

    end subroutine show_state
    !========================================================================
    subroutine show_dt

      integer:: iProcShow
      !---------------------------------------------------------------------
      call barrier_mpi

      do iProcShow = 0, nProc - 1
         if(iProc == iProcShow)then
            do iBlock = 1, nBlock
               if(Unused_B(iBlock)) CYCLE
               write(*,*)'iProc, iBlock, Dx, Dt =', &
                    iProc, iBlock, CellSize_DB(1,iBlock), Dt_B(iBlock)
            end do
         end if
         call barrier_mpi
      end do

    end subroutine show_dt
    !========================================================================

  end subroutine test_amr

end module BATL_amr

