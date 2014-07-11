!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_pass_cell

  use BATL_geometry, ONLY: IsCartesianGrid, IsRotatedCartesian, IsRoundCube, &
       IsCylindricalAxis, IsSphericalAxis, IsLatitudeAxis, Lat_, Theta_
  use ModNumConst, ONLY: cPi, cHalfPi, cTwoPi
  use BATL_high_order, ONLY: get_ghost_for_coarse_blk
  ! Possible improvements:
  ! (1) Instead of sending the receiving block number
  !     and the 2*nDim range limits, we can send only the tag which
  !     we would use in a block to block communication:
  !        iTag = 100*iBlockRecv + iRecv + 4*(jRecv + 4*kRecv)
  !     There are 2 advantages:
  !     (a) The amount of info reduces: 1+2*nDim numbers --> 1 number (iTag)
  !     (b) This procedure allows to send for 1st order prolongation only
  !         one copy per 2**nDim cells
  ! (2) Instead of waiting for receiving buffers from ALL processors, we
  !     we can wait for ANY receiving and already start unpacking
  ! (3) Instead of determining the receive (and send) buffer size during
  !     the message_pass_cell, we can determine the sizes a priori:
  !     (a) We can then allocate a small known buffer size
  !     (b) we do at least two times message_pass_cell per time iteration,
  !         each time determining the buffer size. This would be reduced to
  !         only once (there is a small complication with operator split
  !         schemes)

  implicit none

  SAVE

  private ! except

  public message_pass_cell, message_pass_ng_int1
  public test_pass_cell

  interface message_pass_cell
     module procedure            &
          message_pass_ng_int1,  &  ! Integer scalar with nG ghost cells
          message_pass_ng_real1, &  ! Real scalar with nG ghost cells
          message_pass_real1,    &  ! Real scalar with arbitrary ghost cells
          message_pass_ng_real,  &  ! Real array with nG ghost cells
          message_pass_real         ! Real array with arbitrary ghost cells
  end interface

contains
  !===========================================================================
  subroutine message_pass_ng_real(nVar, State_VGB, &
       nWidthIn, nProlongOrderIn, nCoarseLayerIn, DoSendCornerIn, &
       DoRestrictFaceIn, TimeOld_B, Time_B, DoTestIn, NameOperatorIn,&
       DoResChangeOnlyIn, UseHighResChangeIn)

    ! Message pass real array with nVar variables and BATL_size::nG ghost cells

    use BATL_size, ONLY: MaxBlock, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nG

    ! Arguments
    integer, intent(in)   :: nVar
    real,    intent(inout):: &
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    ! Optional arguments
    integer, optional, intent(in) :: nWidthIn
    integer, optional, intent(in) :: nProlongOrderIn
    integer, optional, intent(in) :: nCoarseLayerIn
    logical, optional, intent(in) :: DoSendCornerIn
    logical, optional, intent(in) :: DoRestrictFaceIn
    logical, optional, intent(in) :: DoTestIn
    logical, optional, intent(in) :: DoResChangeOnlyIn
    logical, optional, intent(in) :: UseHighResChangeIn
    real,    optional, intent(in) :: TimeOld_B(MaxBlock)
    real,    optional, intent(in) :: Time_B(MaxBlock)
    character(len=*), optional,intent(in) :: NameOperatorIn 

    character(len=*), parameter:: NameSub = 'message_pass_ng_real'
    !--------------------------------------------------------------------------

    call message_pass_real(nVar, nG, State_VGB, nWidthIn=nWidthIn, &
         nProlongOrderIn=nProlongOrderIn, nCoarseLayerIn=nCoarseLayerIn, &
         DoSendCornerIn=DoSendCornerIn, DoRestrictFaceIn=DoRestrictFaceIn, &
         TimeOld_B=TimeOld_B, Time_B=Time_B, DoTestIn=DoTestIn, &
         NameOperatorIn=NameOperatorIn, DoResChangeOnlyIn=DoResChangeOnlyIn, &
         UseHighResChangeIn=UseHighResChangeIn)

  end subroutine message_pass_ng_real

  !===========================================================================
  subroutine message_pass_ng_real1(State_GB, &
       nWidthIn, nProlongOrderIn, nCoarseLayerIn, DoSendCornerIn, &
       DoRestrictFaceIn, TimeOld_B, Time_B, DoTestIn, NameOperatorIn,&
       DoResChangeOnlyIn)

    ! Message pass real scalar with BATL_size::nG ghost cells

    use BATL_size, ONLY: MaxBlock, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nG

    ! Arguments
    real, intent(inout):: &
         State_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    ! Optional arguments
    integer, optional, intent(in) :: nWidthIn
    integer, optional, intent(in) :: nProlongOrderIn
    integer, optional, intent(in) :: nCoarseLayerIn
    logical, optional, intent(in) :: DoSendCornerIn
    logical, optional, intent(in) :: DoRestrictFaceIn
    logical, optional, intent(in) :: DoTestIn
    logical, optional, intent(in) :: DoResChangeOnlyIn
    real,    optional, intent(in) :: TimeOld_B(MaxBlock)
    real,    optional, intent(in) :: Time_B(MaxBlock)
    character(len=*), optional,intent(in) :: NameOperatorIn 

    character(len=*), parameter:: NameSub = 'message_pass_ng_real'
    !--------------------------------------------------------------------------

    call message_pass_real(1, nG, State_GB, nWidthIn=nWidthIn, &
         nProlongOrderIn=nProlongOrderIn, nCoarseLayerIn=nCoarseLayerIn, &
         DoSendCornerIn=DoSendCornerIn, DoRestrictFaceIn=DoRestrictFaceIn, &
         TimeOld_B=TimeOld_B, Time_B=Time_B, DoTestIn=DoTestIn, &
         NameOperatorIn=NameOperatorIn, DoResChangeOnlyIn=DoResChangeOnlyIn)

  end subroutine message_pass_ng_real1

  !===========================================================================
  subroutine message_pass_real1(nG, State_GB, &
       nWidthIn, nProlongOrderIn, nCoarseLayerIn, DoSendCornerIn, &
       DoRestrictFaceIn, TimeOld_B, Time_B, DoTestIn, NameOperatorIn,&
       DoResChangeOnlyIn)

    ! Message pass real scalar with BATL_size::nG ghost cells

    use BATL_size, ONLY: MaxBlock, nI, nJ, nK, jDim_, kDim_

    ! Arguments
    integer, intent(in):: nG
    real, intent(inout):: State_GB(1-nG:nI+nG,&
         1-nG*jDim_:nJ+nG*jDim_,1-nG*kDim_:nK+nG*kDim_,MaxBlock)

    ! Optional arguments
    integer, optional, intent(in) :: nWidthIn
    integer, optional, intent(in) :: nProlongOrderIn
    integer, optional, intent(in) :: nCoarseLayerIn
    logical, optional, intent(in) :: DoSendCornerIn
    logical, optional, intent(in) :: DoRestrictFaceIn
    logical, optional, intent(in) :: DoTestIn
    logical, optional, intent(in) :: DoResChangeOnlyIn
    real,    optional, intent(in) :: TimeOld_B(MaxBlock)
    real,    optional, intent(in) :: Time_B(MaxBlock)
    character(len=*), optional,intent(in) :: NameOperatorIn 

    character(len=*), parameter:: NameSub = 'message_pass_real1'
    !--------------------------------------------------------------------------

    call message_pass_real(1, nG, State_GB, nWidthIn=nWidthIn, &
         nProlongOrderIn=nProlongOrderIn, nCoarseLayerIn=nCoarseLayerIn, &
         DoSendCornerIn=DoSendCornerIn, DoRestrictFaceIn=DoRestrictFaceIn, &
         TimeOld_B=TimeOld_B, Time_B=Time_B, DoTestIn=DoTestIn, &
         NameOperatorIn=NameOperatorIn, DoResChangeOnlyIn=DoResChangeOnlyIn)

  end subroutine message_pass_real1

  !===========================================================================
  subroutine message_pass_ng_int1(Int_GB, &
       nWidthIn, nProlongOrderIn, nCoarseLayerIn, DoSendCornerIn, &
       DoRestrictFaceIn, TimeOld_B, Time_B, DoTestIn, NameOperatorIn,&
       DoResChangeOnlyIn)

    ! Message pass scalar integer data with BATL_size::nG ghost cells

    use BATL_size, ONLY: MaxBlock, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nG, &
         nBlock

    ! Arguments
    integer, intent(inout) :: Int_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    ! Optional arguments
    integer, optional, intent(in) :: nWidthIn
    integer, optional, intent(in) :: nProlongOrderIn
    integer, optional, intent(in) :: nCoarseLayerIn
    logical, optional, intent(in) :: DoSendCornerIn
    logical, optional, intent(in) :: DoRestrictFaceIn
    logical, optional, intent(in) :: DoTestIn
    logical, optional, intent(in) :: DoResChangeOnlyIn
    real,    optional, intent(in) :: TimeOld_B(MaxBlock)
    real,    optional, intent(in) :: Time_B(MaxBlock)
    character(len=*), optional,intent(in) :: NameOperatorIn 


    ! help array for converting between Scalar_GB and State_VGB
    ! used by message_pass_cell
    real, allocatable, save:: Scalar_VGB(:,:,:,:,:)

    character(len=*), parameter:: NameSub = 'message_pass_ng_int1'
    !--------------------------------------------------------------------------

    if(.not.allocated(Scalar_VGB)) &
         allocate(Scalar_VGB(1,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))

    Scalar_VGB(1,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,1:nBlock) = &
         Int_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,1:nBlock)

    call message_pass_cell(1, nG, Scalar_VGB, nWidthIn=nWidthIn, &
         nProlongOrderIn=nProlongOrderIn, nCoarseLayerIn=nCoarseLayerIn, &
         DoSendCornerIn=DoSendCornerIn, DoRestrictFaceIn=DoRestrictFaceIn, &
         TimeOld_B=TimeOld_B, Time_B=Time_B, DoTestIn=DoTestIn, &
         NameOperatorIn=NameOperatorIn, DoResChangeOnlyIn=DoResChangeOnlyIn)

    Int_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,1:nBlock) = &
         nint(Scalar_VGB(1,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,1:nBlock))

  end subroutine message_pass_ng_int1

  !============================================================================

  subroutine message_pass_real(nVar, nG, State_VGB, &
       nWidthIn, nProlongOrderIn, nCoarseLayerIn, DoSendCornerIn, &
       DoRestrictFaceIn, TimeOld_B, Time_B, DoTestIn, NameOperatorIn,&
       DoResChangeOnlyIn, UseHighResChangeIn)

    use BATL_size, ONLY: MaxBlock, nBlock, nI, nJ, nK, nIjk_D, &
         MaxDim, nDim, jDim_, kDim_, &
         iRatio, jRatio, kRatio, iRatio_D, InvIjkRatio, &
         MinI, MinJ, MinK, MaxI, MaxJ, MaxK

    use BATL_mpi, ONLY: iComm, nProc, iProc, barrier_mpi

    use BATL_tree, ONLY: &
         iNodeNei_IIIB, DiLevelNei_IIIB, Unused_B, Unused_BP, iNode_B, &
         iTree_IA, Proc_, Block_, Coord1_, Coord2_, Coord3_

    use BATL_grid, ONLY: CoordMin_DB, CoordMax_DB

    use ModMpi

    use ModUtilities, ONLY: lower_case

    ! Arguments
    integer, intent(in) :: nVar  ! number of variables
    integer, intent(in) :: nG    ! number of ghost cells for 1..nDim
    real, intent(inout) :: State_VGB(nVar,&
         1-nG:nI+nG,1-nG*jDim_:nJ+nG*jDim_,1-nG*kDim_:nK+nG*kDim_,MaxBlock)

    ! Optional arguments
    integer, optional, intent(in) :: nWidthIn
    integer, optional, intent(in) :: nProlongOrderIn
    integer, optional, intent(in) :: nCoarseLayerIn
    logical, optional, intent(in) :: DoSendCornerIn
    logical, optional, intent(in) :: DoRestrictFaceIn
    logical, optional, intent(in) :: DoResChangeOnlyIn
    character(len=*), optional,intent(in) :: NameOperatorIn
    real,    optional, intent(in) :: TimeOld_B(MaxBlock)
    real,    optional, intent(in) :: Time_B(MaxBlock)
    logical, optional, intent(in) :: DoTestIn
    logical, optional, intent(in) :: UseHighResChangeIn
    ! Fill in the nVar variables in the ghost cells of State_VGB.
    !
    ! nWidthIn is the number of ghost cell layers to be set. Default is all.
    ! nProlongOrderIn is the order of accuracy of prolongation. Default is 2.
    ! nCoarseLayerIn is the number of coarse layers sent during first order
    !     prolongation. Default is 1, ie all fine cells are equal.
    !     If it is set to 2, the 2 (or more) coarse layers are copied into 
    !     the fine cell layers one by one.
    ! DoSendCornerIn determines if edges/corners are filled. Default is true.
    ! DoRestrictFaceIn determines if restriction is applied to a single layer
    !     of ghost cells instead of two layers. Default is false.
    !     Only works with first order prolongation.
    ! DoResChangeOnlyIn determines if only ghost cells next to resolution
    !    changes are filled in.
    ! NameOperatorIn is used for the minimum or the maximum at the fine
    !    Grid to the course grid cell. If not given the average will be used
    ! TimeOld_B and Time_B are the simulation times associated with the
    !    ghost cells and physical cells of State_VGB, respectively. 
    !    If these arguments are present, the ghost cells are interpolated 
    !    in time. Default is a simple update with no temporal interpolation.
    ! DoTestIn determines if verbose information should be printed

    ! Local variables

    logical, parameter :: UseRSend = .false.

    ! local variables corresponding to optional arguments
    integer :: nWidth
    integer :: nCoarseLayer
    integer :: nProlongOrder
    logical :: DoSendCorner
    logical :: DoRestrictFace
    logical :: DoResChangeOnly
    character(len=4) :: NameOperator
    logical:: UseMin, UseMax  ! logicals for min and max operators
    logical :: UseTime        ! true if Time_B and TimeOld_B are present
    logical :: DoTest
    logical :: UseHighResChange

    ! Various indexes
    integer :: iSendStage  ! index for 2 stage scheme for 2nd order prolong or 
                           ! high order restriction. 
    integer :: iCountOnly     ! index for 2 stage scheme for count, sendrecv
    logical :: DoCountOnly    ! logical for count vs. sendrecv stages

    integer :: iSend, jSend, kSend, iRecv, jRecv, kRecv, iSide, jSide, kSide
    integer :: iDir, jDir, kDir
    integer :: iNodeRecv, iNodeSend
    integer :: iBlockRecv, iProcRecv, iBlockSend, iProcSend, DiLevel

    ! Is the sending node next to the symmetry axis?
    logical :: IsAxisNode

    ! Fast lookup tables for index ranges per dimension
    integer, parameter:: Min_=1, Max_=2
    integer:: iEqualS_DII(MaxDim,-1:1,Min_:Max_)
    integer:: iEqualR_DII(MaxDim,-1:1,Min_:Max_)
    integer:: iRestrictS_DII(MaxDim,-1:1,Min_:Max_)
    integer:: iRestrictR_DII(MaxDim,0:3,Min_:Max_)
    integer:: iProlongS_DII(MaxDim,0:3,Min_:Max_)
    integer:: iProlongR_DII(MaxDim,0:3,Min_:Max_)

    ! Index range for recv and send segments of the blocks
    integer :: iRMin, iRMax, jRMin, jRMax, kRMin, kRMax
    integer :: iSMin, iSMax, jSMin, jSMax, kSMin, kSMax

    ! Message passing across the pole can reverse the recv. index range
    integer :: DiR = 1, DjR = 1, DkR = 1

    ! Variables related to recv and send buffers
    integer, allocatable, save:: iBufferS_P(:), nBufferS_P(:), nBufferR_P(:)

    integer :: iBufferS, iBufferR
    integer :: MaxBufferS = -1, MaxBufferR = -1
    real, allocatable, save:: BufferR_I(:), BufferS_I(:)

    integer:: iRequestR, iRequestS, iError
    integer, allocatable, save:: iRequestR_I(:), iRequestS_I(:), &
         iStatus_II(:,:)

    ! Slopes for 2nd order prolongation
    real, allocatable:: Slope_VG(:,:,:,:)

    character(len=*), parameter:: NameSub = 'BATL_pass_cell::message_pass_cell'

    integer:: nSendStage

    integer:: iBlock
    !--------------------------------------------------------------------------
    DoTest = .false.; if(present(DoTestIn)) DoTest = DoTestIn
    if(DoTest)write(*,*)NameSub,' starting with nVar=',nVar

    call timing_start('batl_pass')

    call timing_start('init_pass')

    ! Set values or defaults for optional arguments
    nWidth = nG
    if(present(nWidthIn)) nWidth = nWidthIn

    nProlongOrder = 2
    if(present(nProlongOrderIn)) nProlongOrder = nProlongOrderIn

    nCoarseLayer = 1
    if(present(nCoarseLayerIn)) nCoarseLayer = nCoarseLayerIn

    DoSendCorner = .true.
    if(present(DoSendCornerIn)) DoSendCorner = DoSendCornerIn

    DoRestrictFace = .false.
    if(present(DoRestrictFaceIn)) DoRestrictFace = DoRestrictFaceIn

    DoResChangeOnly =.false.
    if(present(DoResChangeOnlyIn)) DoResChangeOnly = DoResChangeOnlyIn


    UseHighResChange = .false. 
    if(present(UseHighResChangeIn)) UseHighResChange = UseHighResChangeIn

    ! Check arguments for consistency
    if(nProlongOrder == 2 .and. DoRestrictFace) call CON_stop(NameSub// &
         ' cannot use 2nd order prolongation with face restriction')

    if(nProlongOrder == 2 .and. nCoarseLayer>1) call CON_stop(NameSub// &
         ' cannot use 2nd order prolongation nCoarseLayer > 1')

    if(nProlongOrder < 1 .or. nProlongOrder > 2) call CON_stop(NameSub// &
         ' only nProlongOrder = 1 or 2 is implemented ')

    if(nWidth < 1 .or. nWidth > nG) call CON_stop(NameSub// &
         ' nWidth do not contain the ghost cells or too many')

    if(nCoarseLayer < 1 .or.  nCoarseLayer > 2 ) call CON_stop(NameSub// &
         ' nCoarseLayer are only defined for value or 1 or 2 ')

    UseMin =.false.
    UseMax =.false.

    if(present(NameOperatorIn)) then
       NameOperator = adjustl(NameOperatorIn)
       call lower_case(NameOperator)
       select case(NameOperator)
       case("min")
          UseMin=.true.
       case("max")
          UseMax=.true.
       end select
    end if

    if(present(Time_B) .and. present(NameOperatorIn)) then
       call CON_stop(NameSub// &
            ': Time_B can not be used with '//trim(NameOperator))
    end if

    if(DoTest)write(*,*) NameSub, &
         ': Width, Prolong, Coarse, Corner, RestrictFace, ResChangeOnly=', &
         nWidth, nProlongOrder, nCoarseLayer, DoSendCorner, &
         DoRestrictFace, DoResChangeOnly

    ! Initialize logical for time interpolation/extrapolation
    UseTime = .false.

    ! Set index ranges based on arguments
    call set_range

    if(nProc > 1)then
       ! Allocate fixed size communication arrays.
       if(.not.allocated(iBufferS_P))then
          allocate(iBufferS_P(0:nProc-1), nBufferS_P(0:nProc-1), &
               nBufferR_P(0:nProc-1))
          allocate(iRequestR_I(nProc-1), iRequestS_I(nProc-1))
          allocate(iStatus_II(MPI_STATUS_SIZE,nProc-1))
       end if
    end if

    ! Allocate slope for prolongation. Size depends on nVar and nWidth
    allocate(Slope_VG(nVar,1-nWidth:nI+nWidth,&
         1-nWidth*jDim_:nJ+nWidth*jDim_,1-nWidth*kDim_:nK+nWidth*kDim_))
    ! Set to zero so we can add it for first order prolongation too
    Slope_VG = 0.0

    call timing_stop('init_pass')

    nSendStage = nProlongOrder
    if(UseHighResChange) nSendStage = 2
    do iSendStage = 1, nSendStage
       do iCountOnly = 1, 2
          DoCountOnly = iCountOnly == 1

          ! No need to count data for send/recv in serial runs
          if(nProc == 1 .and. DoCountOnly) CYCLE

          call timing_start('local_pass')

          ! Second order prolongation needs two stages: 
          ! first stage fills in equal and coarser ghost cells
          ! second stage uses these to prolong and fill in finer ghost cells

          if(nProc>1)then
             if(DoCountOnly)then
                ! initialize buffer size
                nBufferR_P = 0
                nBufferS_P = 0
             else
                ! Make buffers large enough
                if(sum(nBufferR_P) > MaxBufferR) then
                   if(allocated(BufferR_I)) deallocate(BufferR_I)
                   MaxBufferR = sum(nBufferR_P) 
                   allocate(BufferR_I(MaxBufferR))
                end if

                if(sum(nBufferS_P) > MaxBufferS) then
                   if(allocated(BufferS_I)) deallocate(BufferS_I)
                   MaxBufferS = sum(nBufferS_P)
                   allocate(BufferS_I(MaxBufferS))
                end if

                ! Initialize buffer indexes for storing data into BufferS_I
                iBufferS = 0
                do iProcRecv = 0, nProc-1
                   iBufferS_P(iProcRecv) = iBufferS
                   iBufferS = iBufferS + nBufferS_P(iProcRecv)
                end do
             end if
          end if

          ! Loop through all blocks that may send a message
          do iBlockSend = 1, nBlock

             if(Unused_B(iBlockSend)) CYCLE

             iNodeSend = iNode_B(iBlockSend)

             IsAxisNode = .false.

             do kDir = -1, 1
                ! Do not message pass in ignored dimensions
                if(nDim < 3 .and. kDir /= 0) CYCLE

                if(nDim > 2 .and. IsLatitudeAxis) IsAxisNode = &
                     kDir == -1 .and. &
                     CoordMin_DB(Lat_,iBlockSend) < -cHalfPi + 1e-8 .or. &
                     kDir == +1 .and. &
                     CoordMax_DB(Lat_,iBlockSend) > +cHalfPi - 1e-8

                do jDir = -1, 1
                   if(nDim < 2 .and. jDir /= 0) CYCLE
                   ! Skip edges
                   if(.not.DoSendCorner .and. jDir /= 0 .and. kDir /= 0) CYCLE

                   if(nDim > 2 .and. IsSphericalAxis) IsAxisNode = &
                        jDir == -1 .and. &
                        CoordMin_DB(Theta_,iBlockSend) < 1e-8 .or. &
                        jDir == +1 .and. &
                        CoordMax_DB(Theta_,iBlockSend) > cPi-1e-8

                   do iDir = -1,1
                      ! Ignore inner parts of the sending block
                      if(iDir == 0 .and. jDir == 0 .and. kDir == 0) CYCLE

                      ! Exclude corners where i and j or k is at the edge
                      if(.not.DoSendCorner .and. iDir /= 0 .and. &
                           (jDir /= 0 .or.  kDir /= 0)) CYCLE

                      if(nDim > 1 .and. IsCylindricalAxis) IsAxisNode = &
                           iDir == -1 .and. iTree_IA(Coord1_,iNodeSend) == 1

                      DiLevel = DiLevelNei_IIIB(iDir,jDir,kDir,iBlockSend)

                      ! Do prolongation in the second stage if nProlongOrder=2
                      ! We still need to call restriction and prolongation in
                      ! both stages to calculate the amount of received data
                      if(iSendStage == 2 .and. DiLevel == 0) CYCLE

                      if(DiLevel == 0)then
                         if(.not.DoResChangeOnly) call do_equal
                      elseif(DiLevel == 1)then
                         call do_restrict
                      elseif(DiLevel == -1)then
                         call do_prolong
                      endif
                   end do ! iDir
                end do ! jDir
             end do ! kDir
          end do ! iBlockSend

          call timing_stop('local_pass')

       end do ! iCountOnly

       ! Done for serial run
       if(nProc == 1) CYCLE
       
       call timing_start('recv_pass')

       ! post requests
       iRequestR = 0
       iBufferR  = 1
       do iProcSend = 0, nProc - 1
          if(nBufferR_P(iProcSend) == 0) CYCLE
          iRequestR = iRequestR + 1

          call MPI_irecv(BufferR_I(iBufferR), nBufferR_P(iProcSend), &
               MPI_REAL, iProcSend, 10, iComm, iRequestR_I(iRequestR), &
               iError)

          iBufferR  = iBufferR  + nBufferR_P(iProcSend)
       end do
       
       call timing_stop('recv_pass')

       if(UseRSend) then
          call timing_start('barrier_pass')
          call barrier_mpi
          call timing_stop('barrier_pass')
       end if

       call timing_start('send_pass')

       ! post sends
       iRequestS = 0
       iBufferS  = 1
       do iProcRecv = 0, nProc-1
          if(nBufferS_P(iProcRecv) == 0) CYCLE
          iRequestS = iRequestS + 1

          if(UseRSend)then
             call MPI_rsend(BufferS_I(iBufferS), nBufferS_P(iProcRecv), &
                  MPI_REAL, iProcRecv, 10, iComm, iError)
          else
             call MPI_isend(BufferS_I(iBufferS), nBufferS_P(iProcRecv), &
                  MPI_REAL, iProcRecv, 10, iComm, iRequestS_I(iRequestS), &
                  iError)
          end if

          iBufferS  = iBufferS  + nBufferS_P(iProcRecv)
       end do
       call timing_stop('send_pass')

       call timing_start('wait_pass')
       ! wait for all requests to be completed
       if(iRequestR > 0) &
            call MPI_waitall(iRequestR, iRequestR_I, iStatus_II, iError)
       
       ! wait for all sends to be completed
       if(.not.UseRSend .and. iRequestS > 0) &
            call MPI_waitall(iRequestS, iRequestS_I, iStatus_II, iError)
       call timing_stop('wait_pass')

       call timing_start('buffer_to_state')
       call buffer_to_state
       call timing_stop('buffer_to_state')
       
    end do ! iSendStage

     if(UseHighResChange) then
        do iBlock = 1, nBlock
           if (Unused_B(iBlock)) CYCLE
           call calc_high_ghost(iBlock)
        enddo
     endif

    deallocate(Slope_VG)

    call timing_stop('batl_pass')

  contains
    !======================================================================
    subroutine calc_high_ghost(iBlock) 
      ! The ghost cells of coarse blocks have been got in do_restrict
      ! Ghost cells of fine blocks are calculated here. 

      use BATL_high_order, ONLY: calc_high_ghost_for_fine_blk, &
           correct_ghost_for_fine_blk, correct_ghost_for_coarse_blk

      integer, intent(in):: iBlock
      real, allocatable:: Field1_VG(:,:,:,:)
      integer:: neiLev_I(6)
      !----------------------------------------------------------------------

      if(.not. allocated(Field1_VG)) &
           allocate(Field1_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK))            

      call calc_high_ghost_for_fine_blk(&
           iBlock, nVar, Field1_VG, State_VGB(:,:,:,:,iBlock))

      neiLev_I(1) = DiLevelNei_IIIB(-1,0,0,iBlock)
      neiLev_I(2) = DiLevelNei_IIIB(+1,0,0,iBlock)
      neiLev_I(3) = DiLevelNei_IIIB(0,-1,0,iBlock)
      neiLev_I(4) = DiLevelNei_IIIB(0,+1,0,iBlock)
      neiLev_I(5) = DiLevelNei_IIIB(0,0,-1,iBlock)
      neiLev_I(6) = DiLevelNei_IIIB(0,0,+1,iBlock)

      ! If the corner block is not a coarse block, the ghost values for 
      ! fine block need to be corrected. 
      if(.not. all(neiLev_I /=1)) call correct_ghost_for_fine_blk(&
           iBlock, nVar, State_VGB(:,:,:,:,iBlock))

      ! If the corner block is not a fine block, the ghost values for 
      ! coarse block need to be corrected.     
      if(.not. all(neiLev_I /= -1)) call correct_ghost_for_coarse_blk(&
           iBlock, nVar, State_VGB(:,:,:,:,iBlock))

    end subroutine calc_high_ghost

    !==========================================================================
    subroutine buffer_to_state

      ! Copy buffer into recv block of State_VGB

      integer:: iBufferR, i, j, k
      real :: TimeSend, WeightOld, WeightNew
      !------------------------------------------------------------------------

      jRMin = 1; jRMax = 1
      kRMin = 1; kRMax = 1

      iBufferR = 0
      do iProcSend = 0, nProc - 1
         if(nBufferR_P(iProcSend) == 0) CYCLE

         do
            iBlockRecv = nint(BufferR_I(iBufferR+1))
            iRMin      = nint(BufferR_I(iBufferR+2))
            iRMax      = nint(BufferR_I(iBufferR+3))
            if(nDim > 1) DiR = sign(1,iRMax - iRMin)
            if(nDim > 1) jRMin = nint(BufferR_I(iBufferR+4))
            if(nDim > 1) jRMax = nint(BufferR_I(iBufferR+5))
            if(nDim > 2) DjR   = sign(1, jRmax - jRMin)
            if(nDim > 2) kRMin = nint(BufferR_I(iBufferR+6))
            if(nDim > 2) kRMax = nint(BufferR_I(iBufferR+7))
            if(nDim > 2) DkR   = sign(1, kRmax - kRMin)

            iBufferR = iBufferR + 1 + 2*nDim
            if(present(Time_B))then
               ! Get time of neighbor and interpolate/extrapolate ghost cells
               iBufferR = iBufferR + 1
               TimeSend  = BufferR_I(iBufferR)
               UseTime = abs(TimeSend - Time_B(iBlockRecv)) > 1e-30
            end if

            if(UseTime)then
               WeightOld = (TimeSend - Time_B(iBlockRecv)) &
                    /      (TimeSend - TimeOld_B(iBlockRecv))
               WeightNew = 1 - WeightOld
               do k=kRMin,kRmax,DkR; do j=jRMin,jRMax,DjR; do i=iRMin,iRmax,DiR
                  State_VGB(:,i,j,k,iBlockRecv) = &
                       WeightOld*State_VGB(:,i,j,k,iBlockRecv) + &
                       WeightNew*BufferR_I(iBufferR+1:iBufferR+nVar)

                  iBufferR = iBufferR + nVar
               end do; end do; end do
            else
               do k=kRMin,kRmax,DkR; do j=jRMin,jRMax,DjR; do i=iRMin,iRmax,DiR
                  State_VGB(:,i,j,k,iBlockRecv) = &
                       BufferR_I(iBufferR+1:iBufferR+nVar)

                  iBufferR = iBufferR + nVar
               end do; end do; end do
            end if
            if(iBufferR >= sum(nBufferR_P(0:iProcSend))) EXIT
         end do
      end do

    end subroutine buffer_to_state

    !==========================================================================

    subroutine do_equal

      integer :: iBufferS, i, j, k, nSize
      real    :: WeightOld, WeightNew
      !------------------------------------------------------------------------

      iSend = (3*iDir + 3)/2
      jSend = (3*jDir + 3)/2
      kSend = (3*kDir + 3)/2

      iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)
      iProcRecv  = iTree_IA(Proc_,iNodeRecv)
      iBlockRecv = iTree_IA(Block_,iNodeRecv)

      ! For part implicit and part steady schemes
      if(Unused_BP(iBlockRecv,iProcRecv)) RETURN

      ! No need to count data for local copy
      if(DoCountOnly .and. iProc == iProcRecv) RETURN

      iRMin = iEqualR_DII(1,iDir,Min_)
      iRMax = iEqualR_DII(1,iDir,Max_)
      jRMin = iEqualR_DII(2,jDir,Min_)
      jRMax = iEqualR_DII(2,jDir,Max_)
      kRMin = iEqualR_DII(3,kDir,Min_)
      kRMax = iEqualR_DII(3,kDir,Max_)

      if(DoCountOnly)then
         ! Number of reals to send to and received from the other processor
         nSize = nVar*(iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1) &
              + 1 + 2*nDim
         if(present(Time_B)) nSize = nSize + 1
         nBufferR_P(iProcRecv) = nBufferR_P(iProcRecv) + nSize
         nBufferS_P(iProcRecv) = nBufferS_P(iProcRecv) + nSize
         RETURN
      end if

      if(IsAxisNode)then
         if(IsLatitudeAxis)then
            kRMin = iEqualR_DII(3,-kDir,Max_)
            kRMax = iEqualR_DII(3,-kDir,Min_)
         elseif(IsSphericalAxis)then
            jRMin = iEqualR_DII(2,-jDir,Max_)
            jRMax = iEqualR_DII(2,-jDir,Min_)
         elseif(IsCylindricalAxis)then
            iRMin = iEqualR_DII(1,1,Max_)
            iRMax = iEqualR_DII(1,1,Min_)
         end if
      end if

      iSMin = iEqualS_DII(1,iDir,Min_)
      iSMax = iEqualS_DII(1,iDir,Max_)
      jSMin = iEqualS_DII(2,jDir,Min_)
      jSMax = iEqualS_DII(2,jDir,Max_)
      kSMin = iEqualS_DII(3,kDir,Min_)
      kSMax = iEqualS_DII(3,kDir,Max_)

      !if(DoTest)then
      !   write(*,*)'!!! iSMin, iSMax, jSMin, jSMax, kSMin, kSMax=', &
      !        iSMin, iSMax, jSMin, jSMax, kSMin, kSMax
      !   write(*,*)'!!! iRMin, iRMax, jRMin, jRMax, kRMin, kRMax=', &
      !        iRMin, iRMax, jRMin, jRMax, kRMin, kRMax
      !   write(*,*)'!!! State(SendMin)=', &
      !        State_VGB(:,iSMin,jSMin,kSMin,iBlockSend)
      !end if

      if(iProc == iProcRecv)then
         ! Local copy
         if(nDim > 1) DiR = sign(1, iRMax - iRMin)
         if(nDim > 2) DjR = sign(1, jRMax - jRMin)
         if(nDim > 2) DkR = sign(1, kRMax - kRMin)

         if(present(Time_B)) UseTime = &
              abs(Time_B(iBlockSend) - Time_B(iBlockRecv)) > 1e-30
         if(UseTime)then
            WeightOld = (Time_B(iBlockSend) - Time_B(iBlockRecv)) &
                 /      (Time_B(iBlockSend) - TimeOld_B(iBlockRecv))
            WeightNew = 1 - WeightOld
            State_VGB(:,iRMin:iRMax:DiR,jRMin:jRMax:DjR,kRMin:kRMax:DkR,&
                 iBlockRecv)= WeightOld* &
                 State_VGB(:,iRMin:iRMax:DiR,jRMin:jRMax:DjR,kRMin:kRMax:DkR, &
                 iBlockRecv) + WeightNew* &
                 State_VGB(:,iSMin:iSMax,jSMin:jSMax,kSMin:kSMax,iBlockSend)
         else
            State_VGB(:,iRMin:iRMax:DiR,jRMin:jRMax:DjR,kRMin:kRMax:DkR, &
                 iBlockRecv)= &
                 State_VGB(:,iSMin:iSMax,jSMin:jSMax,kSMin:kSMax,iBlockSend)
         end if
      else
         ! Put data into the send buffer
         iBufferS = iBufferS_P(iProcRecv)

         BufferS_I(            iBufferS+1) = iBlockRecv
         BufferS_I(            iBufferS+2) = iRMin
         BufferS_I(            iBufferS+3) = iRMax
         if(nDim > 1)BufferS_I(iBufferS+4) = jRMin
         if(nDim > 1)BufferS_I(iBufferS+5) = jRMax
         if(nDim > 2)BufferS_I(iBufferS+6) = kRMin
         if(nDim > 2)BufferS_I(iBufferS+7) = kRMax

         iBufferS = iBufferS + 1 + 2*nDim

         if(present(Time_B))then
            iBufferS = iBufferS + 1
            BufferS_I(iBufferS) = Time_B(iBlockSend)
         end if

         do k = kSMin,kSmax; do j = jSMin,jSMax; do i = iSMin,iSmax
            BufferS_I(iBufferS+1:iBufferS+nVar) = State_VGB(:,i,j,k,iBlockSend)
            iBufferS = iBufferS + nVar
         end do; end do; end do

         iBufferS_P(iProcRecv) = iBufferS

      end if

    end subroutine do_equal

    !==========================================================================

    subroutine do_restrict

      integer :: iR, jR, kR, iS1, jS1, kS1, iS2, jS2, kS2, iVar
      integer :: iRatioRestr, jRatioRestr, kRatioRestr
      real    :: InvIjkRatioRestr
      integer :: iBufferS, nSize, i, j, k
      real    :: WeightOld, WeightNew

      real, allocatable:: Primitive_VIII(:,:,:,:)
      real:: Ghost_I(3)

      integer:: iS0,iS, jS0, jS, kS
      integer:: iBegin, iEnd, jBegin, jEnd, kBegin, kEnd, Di, Dj, Dk
      real :: CoarseCell
      real, allocatable:: State_VG(:,:,:,:)
      !------------------------------------------------------------------------

      ! For sideways communication from a fine to a coarser block
      ! the coordinate parity of the sender block tells 
      ! if the receiver block fills into the 
      ! lower (D*Recv = 0) or upper (D*Rev=1) half of the block
      iSide = 0; if(iRatio==2) iSide = modulo(iTree_IA(Coord1_,iNodeSend)-1, 2)
      jSide = 0; if(jRatio==2) jSide = modulo(iTree_IA(Coord2_,iNodeSend)-1, 2)
      kSide = 0; if(kRatio==2) kSide = modulo(iTree_IA(Coord3_,iNodeSend)-1, 2)

      ! Do not restrict diagonally in the direction of the sibling.
      if(iDir == -1 .and. iSide==1 .and. iRatio == 2) RETURN
      if(iDir == +1 .and. iSide==0 .and. iRatio == 2) RETURN
      if(jDir == -1 .and. jSide==1 .and. jRatio == 2) RETURN
      if(jDir == +1 .and. jSide==0 .and. jRatio == 2) RETURN
      if(kDir == -1 .and. kSide==1 .and. kRatio == 2) RETURN
      if(kDir == +1 .and. kSide==0 .and. kRatio == 2) RETURN

      iSend = (3*iDir + 3 + iSide)/2
      jSend = (3*jDir + 3 + jSide)/2
      kSend = (3*kDir + 3 + kSide)/2

      iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)
      iProcRecv  = iTree_IA(Proc_,iNodeRecv)
      iBlockRecv = iTree_IA(Block_,iNodeRecv)

      ! For part implicit and part steady schemes
      if(Unused_BP(iBlockRecv,iProcRecv)) RETURN

      ! No need to count data for local copy
      if(DoCountOnly .and. iProc == iProcRecv) RETURN

      if(DoCountOnly .and. (&
           (.not. UseHighResChange .and. iSendStage == nProlongOrder) .or. &
           (UseHighResChange .and. iSendStage == 1)))then 
         ! For high resolution change, finer block only receives data 
         ! when iSendStage = 1. 

         ! This processor will receive a prolonged buffer from
         ! the other processor and the "recv" direction of the prolongations
         ! will be the same as the "send" direction for this restriction:
         ! iSend,kSend,jSend = 0..3
         iRMin = iProlongR_DII(1,iSend,Min_)
         iRMax = iProlongR_DII(1,iSend,Max_)
         jRMin = iProlongR_DII(2,jSend,Min_)
         jRMax = iProlongR_DII(2,jSend,Max_)
         kRMin = iProlongR_DII(3,kSend,Min_)
         kRMax = iProlongR_DII(3,kSend,Max_)

         nSize = nVar*(iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1) &
              + 1 + 2*nDim
         if(present(Time_B)) nSize = nSize + 1
         nBufferR_P(iProcRecv) = nBufferR_P(iProcRecv) + nSize
      end if

      ! If this is the pure prolongation stage, all we did was counting
      if(iSendStage == 2 .and. .not. UseHighResChange) RETURN

      ! For high resolution change, the finer block receives data from 
      ! coarser or equal blocks when iSendStage = 1. Restriction will 
      ! be done when iSendStage = 2. 
      if(UseHighResChange .and. iSendStage == 1) RETURN

      iRecv = iSend - 3*iDir
      jRecv = jSend - 3*jDir
      kRecv = kSend - 3*kDir

      ! Receiving range depends on iRecv,kRecv,jRecv = 0..3
      iRMin = iRestrictR_DII(1,iRecv,Min_)
      iRMax = iRestrictR_DII(1,iRecv,Max_)
      jRMin = iRestrictR_DII(2,jRecv,Min_)
      jRMax = iRestrictR_DII(2,jRecv,Max_)
      kRMin = iRestrictR_DII(3,kRecv,Min_)
      kRMax = iRestrictR_DII(3,kRecv,Max_)

      if(DoCountOnly)then
         ! Number of reals to send to the other processor
         nSize = nVar*(iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1) &
              + 1 + 2*nDim 
         if(present(Time_B)) nSize = nSize + 1
         nBufferS_P(iProcRecv) = nBufferS_P(iProcRecv) + nSize
         RETURN
      end if

      if(IsAxisNode)then
         if(IsLatitudeAxis)then
            kRMin = iRestrictR_DII(3,kSend,Max_)
            kRMax = iRestrictR_DII(3,kSend,Min_)
         elseif(IsSphericalAxis)then
            jRMin = iRestrictR_DII(2,jSend,Max_)
            jRMax = iRestrictR_DII(2,jSend,Min_)
         elseif(IsCylindricalAxis)then
            iRMin = iRestrictR_DII(1,0,Max_)
            iRMax = iRestrictR_DII(1,0,Min_)
         end if
      end if

      if(nDim > 1) DiR = sign(1, iRMax - iRMin)
      if(nDim > 2) DjR = sign(1, jRMax - jRMin)
      if(nDim > 2) DkR = sign(1, kRMax - kRMin)

      ! Index range that gets restricted depends on iDir,jDir,kDir only
      iSMin = iRestrictS_DII(1,iDir,Min_)
      iSMax = iRestrictS_DII(1,iDir,Max_)
      jSMin = iRestrictS_DII(2,jDir,Min_)
      jSMax = iRestrictS_DII(2,jDir,Max_)
      kSMin = iRestrictS_DII(3,kDir,Min_)
      kSMax = iRestrictS_DII(3,kDir,Max_)

      iRatioRestr = iRatio; jRatioRestr = jRatio; kRatioRestr = kRatio
      InvIjkRatioRestr = InvIjkRatio
      if(DoRestrictFace)then
         if(iDir /= 0) iRatioRestr = 1
         if(jDir /= 0) jRatioRestr = 1
         if(kDir /= 0) kRatioRestr = 1
         InvIjkRatioRestr = 1.0/(iRatioRestr*jRatioRestr*kRatioRestr)
      end if

      !write(*,*)'iDir, jDir, kDir =',iDir, jDir, kDir
      !write(*,*)'iRecv,jRecv,kRecv=',iRecv,jRecv,kRecv
      !
      !write(*,*)'iSMin,iSmax,jSMin,jSMax,kSMin,kSmax=',&
      !     iSMin,iSmax,jSMin,jSMax,kSMin,kSmax
      !
      !write(*,*)'iRMin,iRmax,jRMin,jRMax,kRMin,kRmax=',&
      !     iRMin,iRmax,jRMin,jRMax,kRMin,kRmax
      !
      !write(*,*)'iRatioRestr,InvIjkRatioRestr=',iRatioRestr,InvIjkRatioRestr

      if(UseHighResChange .and. .not. allocated(Primitive_VIII)) &
           allocate(Primitive_VIII(nVar,8,6,min(6,nK)))

      if(iProc == iProcRecv)then

         if(present(Time_B)) UseTime = &
              abs(Time_B(iBlockSend) - Time_B(iBlockRecv)) > 1e-30
         if(UseTime)then

            ! Get time of neighbor and interpolate/extrapolate ghost cells
            WeightOld = (Time_B(iBlockSend) - Time_B(iBlockRecv)) &
                 /      (Time_B(iBlockSend) - TimeOld_B(iBlockRecv))
            WeightNew = 1 - WeightOld

            do kR = kRMin, kRMax, DkR
               kS1 = kSMin + kRatioRestr*abs(kR-kRMin)
               kS2 = kS1 + kRatioRestr - 1
               do jR = jRMin, jRMax, DjR
                  jS1 = jSMin + jRatioRestr*abs(jR-jRMin)
                  jS2 = jS1 + jRatioRestr - 1
                  do iR = iRMin, iRMax, DiR
                     iS1 = iSMin + iRatioRestr*abs(iR-iRMin)
                     iS2 = iS1 + iRatioRestr - 1
                     if(UseMin) then
                        do iVar = 1, nVar
                           State_VGB(iVar,iR,jR,kR,iBlockRecv) = &
                                minval(State_VGB(iVar,iS1:iS2,jS1:jS2,kS1:kS2,&
                                iBlockSend))
                        end do
                     else if(UseMax) then
                        do iVar = 1, nVar
                           State_VGB(iVar,iR,jR,kR,iBlockRecv) = &
                                maxval(State_VGB(iVar,iS1:iS2,jS1:jS2,kS1:kS2,&
                                iBlockSend))
                        end do
                     else
                        do iVar = 1, nVar
                           State_VGB(iVar,iR,jR,kR,iBlockRecv) = &
                                WeightOld*State_VGB(iVar,iR,jR,kR,iBlockRecv)+&
                                WeightNew*InvIjkRatioRestr * &
                                sum(State_VGB(iVar,iS1:iS2,jS1:jS2,kS1:kS2,&
                                iBlockSend))
                        end do
                     end if
                  end do
               end do
            end do

         else
            ! No time interpolation/extrapolation is needed

            if(UseHighResChange .and. iDir /= 0 .and. jDir == 0 .and. kDir == 0) then
               ! Resolution change in x-dir. 
               ! Only work for 2D x-y plane case.
               ! The case iDir = 1 and iDir = -1 are symmetric.

               iBegin = 1; iEnd = 8
               jBegin = 1; jEnd = 6
               if(nK == 1) kS = 1 

               if(iDir == 1) then 
                  iS0 = iSMax
               else
                  iS0 = iSMin
               endif

               do kR = kRMin, kRMax, DkR ! This loop is useless for 2D.  
                  do jR = jRMin, jRMax, DjR
                     jS1 = jSMin + jRatioRestr*abs(jR-jRMin)

                     iS = iS0
                     do i = iBegin, iEnd   
                        do j = jBegin, jEnd
                           jS = j - 3 + jS1

                           ! i is the index along the resolution change direction.
                           ! j is the index parallel to the resolution change direction.
                           Primitive_VIII(:,i,j,1) = &
                                State_VGB(:,iS,jS,kS,iBlockSend)
                        enddo
                        iS = iS - iDir
                     enddo

                     do iVar = 1, nVar
                        CoarseCell = State_VGB(iVar,iS0+iDir,jS1,kS,iBlockSend)

                        call get_ghost_for_coarse_blk(CoarseCell,&
                             Primitive_VIII(iVar,:,:,1), Ghost_I)

                        if(iDir == 1) then
                           iS = 3
                           do iR = iRmax-2, iRmax
                              State_VGB(iVar,iR,jR,kR,iBlockRecv) = Ghost_I(iS)
                              iS = iS -iDir
                           enddo
                        else
                           iS = 1
                           do iR = iRmin, iRmin+2
                              State_VGB(iVar,iR,jR,kR,iBlockRecv) = Ghost_I(iS)
                              iS = iS -iDir
                           enddo
                        endif
                     enddo ! iVar
                  enddo ! jR
               enddo ! kR

            else if(UseHighResChange .and. iDir == 0 .and. jDir /= 0 .and. kDir ==0) then

               iBegin = 1; iEnd = 6
               jBegin = 1; jEnd = 8
               if(nK == 1) kS = 1

               if(jDir == 1) then 
                  jS0 = jSMax
               else
                  jS0 = jSMin
               endif

               do kR = kRMin, kRMax, DkR
                  do iR = iRMin, iRMax, DiR
                     iS1 = iSMin + iRatioRestr*abs(iR-iRMin)

                     jS = jS0
                     do j = jBegin, jEnd
                        do i = iBegin, iEnd
                           iS = i - 3 + iS1

                           ! j is the index along the resolution change direction.
                           ! i is the index parallel to the resolution change direction.
                           Primitive_VIII(:,j,i,1) = &
                                State_VGB(:,iS,jS,kS,iBlockSend) 
                        enddo
                        jS = jS - jDir
                     enddo

                     do iVar = 1, nVar
                        CoarseCell = State_VGB(iVar,iS1,jS0+jDir,kS,iBlockSend)

                        call get_ghost_for_coarse_blk(CoarseCell,&
                             Primitive_VIII(iVar,:,:,1), Ghost_I)

                        if(jDir == 1) then
                           jS = 3
                           do jR = jRmax-2, jRmax
                              State_VGB(iVar,iR,jR,kR,iBlockRecv) = Ghost_I(jS)
                              jS = jS -jDir
                           enddo
                        else
                           jS = 1
                           do jR = jRmin, jRmin+2
                              State_VGB(iVar,iR,jR,kR,iBlockRecv) = Ghost_I(jS)
                              jS = jS -jDir
                           enddo
                        endif

                     enddo ! iVar
                  enddo ! jR
               enddo ! kR
            else
               do kR = kRMin, kRMax, DkR
                  kS1 = kSMin + kRatioRestr*abs(kR-kRMin)
                  kS2 = kS1 + kRatioRestr - 1
                  do jR = jRMin, jRMax, DjR
                     jS1 = jSMin + jRatioRestr*abs(jR-jRMin)
                     jS2 = jS1 + jRatioRestr - 1
                     do iR = iRMin, iRMax, DiR
                        iS1 = iSMin + iRatioRestr*abs(iR-iRMin)
                        iS2 = iS1 + iRatioRestr - 1
                        if(UseMin)then
                           do iVar = 1, nVar
                              State_VGB(iVar,iR,jR,kR,iBlockRecv) = &
                                   minval(State_VGB(iVar,iS1:iS2,jS1:jS2,kS1:kS2,&
                                   iBlockSend))
                           end do
                        else if(UseMax) then
                           do iVar = 1, nVar
                              State_VGB(iVar,iR,jR,kR,iBlockRecv) = &
                                   maxval(State_VGB(iVar,iS1:iS2,jS1:jS2,kS1:kS2,&
                                   iBlockSend))
                           end do
                        else
                           do iVar = 1, nVar
                              State_VGB(iVar,iR,jR,kR,iBlockRecv) = &
                                   InvIjkRatioRestr * &
                                   sum(State_VGB(iVar,iS1:iS2,jS1:jS2,kS1:kS2, &
                                   iBlockSend))                              
                           end do
                        end if
                     end do ! iR
                  end do ! jR
               end do ! kR

            endif ! UseHighResChange
         end if ! UseTime

      else ! iProc /= iProcRecv

         if(UseHighResChange.and. iDir /= 0 .and. jDir == 0 .and. kDir == 0) then
            if(.not. allocated(State_VG)) then
               allocate(State_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK))
               State_VG = 1
            endif

            iBegin = 1; iEnd = 8
            jBegin = 1; jEnd = 6
            if(nK ==1) kS = 1

            if(iDir == 1) then 
               iS0 = iSMax
            else
               iS0 = iSMin
            endif

            do kR = kRMin, kRMax, DkR
               do jR = jRMin, jRMax, DjR
                  jS1 = jSMin + jRatioRestr*abs(jR-jRMin)

                  iS = iS0
                  do i = iBegin,iEnd
                     do j = jBegin,jEnd
                        jS = j - 3 + jS1                        
                        Primitive_VIII(:,i,j,1) = &
                             State_VGB(:,iS,jS,kS,iBlockSend)       
                     enddo
                     iS = iS - iDir
                  enddo

                  do iVar = 1, nVar
                     CoarseCell = State_VGB(iVar,iS0+iDir,jS1,kS,iBlockSend)

                     call get_ghost_for_coarse_blk(CoarseCell,&
                          Primitive_VIII(iVar,:,:,1), Ghost_I)

                     if(iDir == 1) then
                        iS = 3
                        do iR = iRmax-2, iRmax
                           State_VG(iVar,iR,jR,kR) = Ghost_I(iS)
                           iS = iS -iDir
                        enddo
                     else
                        iS = 1
                        do iR = iRmin, iRmin+2
                           State_VG(iVar,iR,jR,kR) = Ghost_I(iS)
                           iS = iS -iDir
                        enddo
                     endif
                  enddo ! iVar
               enddo ! jR
            enddo ! kR

         else if(UseHighResChange .and. iDir == 0 .and. jDir /= 0 .and. kDir == 0) then
            if(.not. allocated(State_VG)) then
               allocate(State_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK))
               State_VG = 1
            endif

            iBegin = 1; iEnd = 6
            jBegin = 1; jEnd = 8
            if(nK == 1) kS = 1

            if(jDir == 1) then 
               jS0 = jSMax
            else
               jS0 = jSMin
            endif

            do kR = kRMin, kRMax, DkR
               do iR = iRMin, iRMax, DiR
                  iS1 = iSMin + iRatioRestr*abs(iR-iRMin)

                  jS = jS0
                  do j = jBegin,jEnd
                     do i = iBegin,iEnd
                        iS = i - 3 + iS1
                        Primitive_VIII(:,j,i,1) = &
                             State_VGB(:,iS,jS,kS,iBlockSend)   
                     enddo
                     jS = jS - jDir
                  enddo

                  do iVar = 1, nVar
                     CoarseCell = State_VGB(iVar,iS1,jS0+jDir,kS,iBlockSend)

                     call get_ghost_for_coarse_blk(CoarseCell,&
                          Primitive_VIII(iVar,:,:,1), Ghost_I)

                     if(jDir == 1) then
                        jS = 3
                        do jR = jRmax-2, jRmax
                           State_VG(iVar,iR,jR,kR) = Ghost_I(jS)
                           jS = jS -jDir
                        enddo
                     else
                        jS = 1
                        do jR = jRmin, jRmin+2
                           State_VG(iVar,iR,jR,kR) = Ghost_I(jS)
                           jS = jS -jDir
                        enddo
                     endif
                  enddo ! iVar
               enddo ! jR
            enddo ! kR

         endif ! UseHighResChange

         !write(*,*) '2 iproc,iprocrecv:', iProc, iProcRecv
         iBufferS = iBufferS_P(iProcRecv)

         BufferS_I(            iBufferS+1) = iBlockRecv
         BufferS_I(            iBufferS+2) = iRMin
         BufferS_I(            iBufferS+3) = iRMax
         if(nDim > 1)BufferS_I(iBufferS+4) = jRMin
         if(nDim > 1)BufferS_I(iBufferS+5) = jRMax
         if(nDim > 2)BufferS_I(iBufferS+6) = kRMin
         if(nDim > 2)BufferS_I(iBufferS+7) = kRMax

         iBufferS = iBufferS + 1 + 2*nDim

         if(present(Time_B))then
            iBufferS = iBufferS + 1
            BufferS_I(iBufferS) = Time_B(iBlockSend)
         end if

         do kR = kRMin, kRMax, DkR
            kS1 = kSMin + kRatioRestr*abs(kR-kRMin)
            kS2 = kS1 + kRatioRestr - 1
            do jR = jRMin, jRMax, DjR
               jS1 = jSMin + jRatioRestr*abs(jR-jRMin)
               jS2 = jS1 + jRatioRestr - 1
               do iR = iRMin, iRMax, DiR
                  iS1 = iSMin + iRatioRestr*abs(iR-iRMin)
                  iS2 = iS1 + iRatioRestr - 1
                  if(UseHighResChange .and. &
                       ((iDir /= 0 .and. jDir == 0 .and. kDir ==0) &
                       .or. (iDir == 0 .and. jDir /=0 .and. kDir ==0))) then
                     do iVar = 1, nVar
                        BufferS_I(iBufferS+iVar) = State_VG(iVar,iR,jR,kR)
                     end do
                  else if(UseMin) then
                     do iVar = 1, nVar
                        BufferS_I(iBufferS+iVar) = &
                             minval(State_VGB(iVar,iS1:iS2,jS1:jS2,kS1:kS2,&
                             iBlockSend))
                     end do
                  else if(UseMax) then
                     do iVar = 1, nVar
                        BufferS_I(iBufferS+iVar) = &
                             maxval(State_VGB(iVar,iS1:iS2,jS1:jS2,kS1:kS2,&
                             iBlockSend))
                     end do
                  else
                     do iVar = 1, nVar
                        BufferS_I(iBufferS+iVar) = &
                             InvIjkRatioRestr * &
                             sum(State_VGB(iVar,iS1:iS2,jS1:jS2,kS1:kS2,&
                             iBlockSend))
                     end do
                  end if
                  iBufferS = iBufferS + nVar
               end do
            end do
         end do
         iBufferS_P(iProcRecv) = iBufferS

      end if

    end subroutine do_restrict

    !==========================================================================

    subroutine do_prolong

      use ModCoordTransform, ONLY: cross_product
      use BATL_grid, ONLY: Xyz_DGB, CoordMin_D, CoordMax_D
      use BATL_tree, ONLY: get_tree_position
      use BATL_geometry, ONLY: coord_to_xyz
      use BATL_size,     ONLY: nDimAmr

      integer :: iR, jR, kR, iS, jS, kS, iS1, jS1, kS1
      integer :: iRatioRestr, jRatioRestr, kRatioRestr
      integer :: iBufferS, nSize
      integer, parameter:: Di=iRatio-1, Dj=jRatio-1, Dk=kRatio-1
      real    :: WeightOld, WeightNew, Weight, WeightI, WeightJ, WeightK, InvV
      real, dimension(MaxDim):: Xyz_D, dI_D, dJ_D, dK_D, dR_D, &
           PositionMinR_D, PositionMaxR_D, CoordMinR_D, CoordMaxR_D, &
           CellSizeR_D, CoordR_D

      logical :: UseSimpleWeights
      !------------------------------------------------------------------------

      UseSimpleWeights = nDim == 1 .or. nDimAmr < nDim &
           .or. IsCartesianGrid .or. IsRotatedCartesian .or. IsRoundCube

      ! Loop through the subfaces or subedges
      do kSide = (1-kDir)/2, 1-(1+kDir)/2, 3-kRatio
         kSend = (3*kDir + 3 + kSide)/2
         kRecv = kSend - 3*kDir
         do jSide = (1-jDir)/2, 1-(1+jDir)/2, 3-jRatio
            jSend = (3*jDir + 3 + jSide)/2
            jRecv = jSend - 3*jDir
            do iSide = (1-iDir)/2, 1-(1+iDir)/2, 3-iRatio
               iSend = (3*iDir + 3 + iSide)/2
               iRecv = iSend - 3*iDir

               iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)
               iProcRecv  = iTree_IA(Proc_,iNodeRecv)
               iBlockRecv = iTree_IA(Block_,iNodeRecv)

               ! For part implicit and part steady schemes
               if(Unused_BP(iBlockRecv,iProcRecv)) CYCLE

               ! No need to count data for local copy
               if(DoCountOnly .and. iProc == iProcRecv) CYCLE

               if(DoCountOnly .and. (.not. UseHighResChange .and. &
                    iSendStage == 1 .or. &
                    (UseHighResChange .and. iSendStage == 2)))then
                  ! This processor will receive a restricted buffer from
                  ! the other processor and the "recv" direction of the
                  ! restriction will be the same as the "send" direction for
                  ! this prolongation: iSend,kSend,jSend = 0..3
                  iRMin = iRestrictR_DII(1,iSend,Min_)
                  iRMax = iRestrictR_DII(1,iSend,Max_)
                  jRMin = iRestrictR_DII(2,jSend,Min_)
                  jRMax = iRestrictR_DII(2,jSend,Max_)
                  kRMin = iRestrictR_DII(3,kSend,Min_)
                  kRMax = iRestrictR_DII(3,kSend,Max_)

                  nSize = nVar*(iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1)&
                       + 1 + 2*nDim
                  if(present(Time_B)) nSize = nSize + 1
                  nBufferR_P(iProcRecv) = nBufferR_P(iProcRecv) + nSize

               end if

               ! For 2nd order prolongation no prolongation is done in stage 1
               if(.not. UseHighResChange .and. iSendStage < nProlongOrder) CYCLE

               ! For HighResChange, only do prolongation in stage 1. 
               if(UseHighResChange .and. iSendStage == 2) CYCLE

               ! Receiving range depends on iRecv,kRecv,jRecv = 0..3
               iRMin = iProlongR_DII(1,iRecv,Min_)
               iRMax = iProlongR_DII(1,iRecv,Max_)
               jRMin = iProlongR_DII(2,jRecv,Min_)
               jRMax = iProlongR_DII(2,jRecv,Max_)
               kRMin = iProlongR_DII(3,kRecv,Min_)
               kRMax = iProlongR_DII(3,kRecv,Max_)

               if(DoCountOnly)then
                  ! Number of reals to send to the other processor
                  nSize = nVar*(iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1)&
                       + 1 + 2*nDim
                  if(present(Time_B)) nSize = nSize + 1
                  nBufferS_P(iProcRecv) = nBufferS_P(iProcRecv) + nSize
                  CYCLE
               end if

               if(IsAxisNode)then
                  if(IsLatitudeAxis)then
                     kRMin = iProlongR_DII(3,kSend,Max_)
                     kRMax = iProlongR_DII(3,kSend,Min_)
                  elseif(IsSphericalAxis)then
                     jRMin = iProlongR_DII(2,jSend,Max_)
                     jRMax = iProlongR_DII(2,jSend,Min_)
                  elseif(IsCylindricalAxis)then
                     iRMin = iProlongR_DII(1,0,Max_)
                     iRMax = iProlongR_DII(1,0,Min_)
                  end if
               end if

               if(nDim > 1) DiR = sign(1, iRMax - iRMin)
               if(nDim > 2) DjR = sign(1, jRMax - jRMin)
               if(nDim > 2) DkR = sign(1, kRMax - kRMin)

               ! Sending range depends on iSend,jSend,kSend = 0..3
               iSMin = iProlongS_DII(1,iSend,Min_)
               iSMax = iProlongS_DII(1,iSend,Max_)
               jSMin = iProlongS_DII(2,jSend,Min_)
               jSMax = iProlongS_DII(2,jSend,Max_)
               kSMin = iProlongS_DII(3,kSend,Min_)
               kSMax = iProlongS_DII(3,kSend,Max_)

               iRatioRestr = iRatio; jRatioRestr = jRatio; kRatioRestr = kRatio
               if(nCoarseLayer > 1)then
                  if(iDir /= 0) iRatioRestr = 1
                  if(jDir /= 0) jRatioRestr = 1
                  if(kDir /= 0) kRatioRestr = 1
               end if

               ! if(DoTest)then
               !    write(*,*)'iNodeSend, iProc    , iBlockSend=', &
               !         iNodeSend, iProc, iBlockSend
               !    write(*,*)'iNodeRecv, iProcRecv, iBlockRecv=', &
               !         iNodeRecv, iProcRecv, iBlockRecv
               !    write(*,*)'iSide,jSide,kSide=',iSide,jSide,kSide
               !    write(*,*)'iSend,jSend,kSend=',iSend,jSend,kSend
               !    write(*,*)'iRecv,jRecv,kRecv=',iRecv,jRecv,kRecv
               !    write(*,*)'iSMin,iSmax,jSMin,jSMax,kSMin,kSmax=',&
               !         iSMin,iSmax,jSMin,jSMax,kSMin,kSmax
               !    write(*,*)'iRMin,iRmax,jRMin,jRMax,kRMin,kRmax=',&
               !         iRMin,iRmax,jRMin,jRMax,kRMin,kRmax
               ! end if

               if(nProlongOrder == 2)then
                  ! Add up 2nd order corrections for all AMR dimensions
                  ! Use simple interpolation, should be OK for ghost cells
                  Slope_VG(:,iRMin:iRmax:DiR,jRMin:jRMax:DjR,kRMin:kRMax:DkR)&
                       = 0.0

                  if(.not.UseSimpleWeights .and. iProcRecv /= iProc)then
                     call get_tree_position(iNodeRecv, &
                          PositionMinR_D, PositionMaxR_D)
                     CoordMinR_D = CoordMin_D &
                          + (CoordMax_D - CoordMin_D)*PositionMinR_D
                     CoordMaxR_D = CoordMin_D &
                          + (CoordMax_D - CoordMin_D)*PositionMaxR_D
                     CellSizeR_D = (CoordMaxR_D - CoordMinR_D)/nIjk_D
                  end if

                  do kR = kRMin, kRMax, DkR
                     ! For kRatio = 1 simple shift: kS = kSMin + |kR - kRMin|
                     ! For kRatio = 2 coarsen both kR and kRMin before shift
                     ! We add 9 both to kR and kRMin before dividing by kRatio
                     ! so that all values remain positive and get rounded down.
                     ! This works up to nG=10 ghost cells: likely to be enough.
                     kS = kSMin + abs((kR+9)/kRatio - (kRMin+9)/kRatio)

                     ! DkR=+1: interpolate left for odd kR, right for even kR 
                     ! DkR=-1: interpolate left for even kR, right for odd kR 
                     if(kRatio == 1) kS1 = kS
                     if(kRatio == 2) kS1 = kS + DkR*(1 - 2*modulo(kR,2))

                     do jR = jRMin, jRMax, DjR
                        jS = jSMin + abs((jR+9)/jRatio - (jRMin+9)/jRatio)
                        if(jRatio == 1) jS1 = jS
                        if(jRatio == 2) jS1 = jS + DjR*(1 - 2*modulo(jR,2))

                        do iR = iRMin, iRMax, DiR
                           iS = iSMin + abs((iR+9)/iRatio - (iRMin+9)/iRatio)

                           if(iRatio == 1) iS1 = iS
                           if(iRatio == 2) iS1 = iS + DiR*(1 - 2*modulo(iR,2))

                           if(UseSimpleWeights)then
                              ! For Cartesian-like grids the weights are 0.25
                              if(iRatio == 2) WeightI = 0.25
                              if(jRatio == 2) WeightJ = 0.25
                              if(kRatio == 2) WeightK = 0.25
                           else
                              ! The weights are area/volume fractions
                              Xyz_D= Xyz_DGB(:,iS,jS,kS,iBlockSend)
                              dI_D = Xyz_DGB(:,iS1,jS,kS,iBlockSend) - Xyz_D
                              dJ_D = Xyz_DGB(:,iS,jS1,kS,iBlockSend) - Xyz_D
                              dK_D = Xyz_DGB(:,iS,jS,kS1,iBlockSend) - Xyz_D

                              if(iProcRecv == iProc)then
                                 dR_D = Xyz_DGB(:,iR,jR,kR,iBlockRecv) - Xyz_D
                              else
                                 CoordR_D = CoordMinR_D + &
                                      ((/iR,jR,kR/) - 0.5)*CellSizeR_D
                                 call coord_to_xyz(CoordR_D, dR_D)
                                 dR_D = dR_D - Xyz_D
                              end if

                              ! The max(0.0, and the normalization to 1
                              ! avoids extrapolation when the
                              ! receiving point is outside the sending 
                              ! polyhedron. Remove these for exact
                              ! second order test.
                              if(nDim == 2)then
                                 InvV = 1/ &
                                      (dI_D(1)*dJ_D(2)-dI_D(2)*dJ_D(1))
                                 WeightI = max(0.0, InvV* &
                                      (dR_D(1)*dJ_D(2)-dR_D(2)*dJ_D(1)))
                                 WeightJ = max(0.0, InvV* &
                                      (dI_D(1)*dR_D(2)-dI_D(2)*dR_D(1)))
                                 Weight = WeightI + WeightJ
                                 if(Weight > 1)then
                                    WeightI = WeightI / Weight
                                    WeightJ = WeightJ / Weight
                                 end if
                              else
                                 InvV = 1/ &
                                      sum(dI_D*cross_product(dJ_D,dK_D))
                                 WeightI = max(0.0, InvV* &
                                      sum(dR_D*cross_product(dJ_D,dK_D)))
                                 WeightJ = max(0.0, InvV* &
                                      sum(dI_D*cross_product(dR_D,dK_D)))
                                 WeightK = max(0.0, InvV* &
                                      sum(dI_D*cross_product(dJ_D,dR_D)))

                                 Weight  = WeightI + WeightJ + WeightK
                                 if(Weight > 1.0)then
                                    WeightI = WeightI / Weight
                                    WeightJ = WeightJ / Weight
                                    WeightK = WeightK / Weight
                                 end if


                                 !if(IsSphericalAxis .and. &
                                 !     iNodeSend==20 .and. iNodeRecv==26 .and.&
                                 !     iR==7 .and. jR==7 .and. kR==-1)then
                                 !   write(*,*)'!!! iNodeSend, iBlockSend=', &
                                 !        iNodeSend, iBlockSend
                                 !   write(*,*)'!!! iS,jS,kS,iS1,jS1,kS1=', &
                                 !        iS,jS,kS,iS1,jS1,kS1
                                 !   write(*,*)'!!! Xyz_D=',Xyz_D
                                 !   write(*,*)'!!! dI_D =',dI_D
                                 !   write(*,*)'!!! dJ_D =',dJ_D
                                 !   write(*,*)'!!! dK_D =',dK_D
                                 !   write(*,*)'!!! dR_D =',dR_D
                                 !   write(*,*)'InvV,WeightI, J, K=',&
                                 !        InvV, WeightI, WeightJ, WeightK
                                 !end if


                              end if
                           end if

                           if(iRatio == 2) Slope_VG(:,iR,jR,kR) = &
                                Slope_VG(:,iR,jR,kR) + WeightI* &
                                ( State_VGB(:,iS1,jS,kS,iBlockSend) &
                                - State_VGB(:,iS ,jS,kS,iBlockSend) )

                           if(jRatio == 2) Slope_VG(:,iR,jR,kR) = &
                                Slope_VG(:,iR,jR,kR) + WeightJ* &
                                ( State_VGB(:,iS,jS1,kS,iBlockSend) &
                                - State_VGB(:,iS,jS ,kS,iBlockSend) )

                           if(kRatio == 2) Slope_VG(:,iR,jR,kR) = &
                                Slope_VG(:,iR,jR,kR) + WeightK* &
                                ( State_VGB(:,iS,jS,kS1,iBlockSend) &
                                - State_VGB(:,iS,jS,kS ,iBlockSend) )

                        end do
                     end do
                  end do
               end if ! nProlongOrder = 2

               if(iProc == iProcRecv)then

                  if(present(Time_B)) UseTime = &
                       abs(Time_B(iBlockSend) - Time_B(iBlockRecv)) > 1e-30
                  if(UseTime)then
                     ! Interpolate/extrapolate ghost cells in time
                     WeightOld = (Time_B(iBlockSend) - Time_B(iBlockRecv)) &
                          /      (Time_B(iBlockSend) - TimeOld_B(iBlockRecv))
                     WeightNew = 1 - WeightOld

                     do kR = kRMin, kRMax, DkR
                        ! For kRatio = 1 simple shift: kS = kSMin + kR - kRMin 
                        ! For kRatio = 2 coarsen both kR and kRMin before shift
                        kS = kSMin + abs((kR+9)/kRatioRestr &
                             -           (kRMin+9)/kRatioRestr)
                        do jR = jRMin, jRMax, DjR
                           jS = jSMin + abs((jR+9)/jRatioRestr &
                                -           (jRMin+9)/jRatioRestr)
                           do iR = iRMin, iRMax, DiR
                              iS = iSMin + abs((iR+9)/iRatioRestr &
                                   -           (iRMin+9)/iRatioRestr)
                              State_VGB(:,iR,jR,kR,iBlockRecv) = &
                                   WeightOld*State_VGB(:,iR,jR,kR,iBlockRecv)+&
                                   WeightNew*(State_VGB(:,iS,jS,kS,iBlockSend)&
                                   +          Slope_VG(:,iR,jR,kR))
                           end do
                        end do
                     end do
                  else
                     do kR = kRMin, kRMax, DkR
                        kS = kSMin + abs((kR+9)/kRatioRestr &
                             -           (kRMin+9)/kRatioRestr)
                        do jR = jRMin, jRMax, DjR
                           jS = jSMin + abs((jR+9)/jRatioRestr &
                                -           (jRMin+9)/jRatioRestr)
                           do iR = iRMin, iRMax, DiR
                              iS = iSMin + abs((iR+9)/iRatioRestr &
                                   -           (iRMin+9)/iRatioRestr)

                              State_VGB(:,iR,jR,kR,iBlockRecv) = &
                                   State_VGB(:,iS,jS,kS,iBlockSend) &
                                   + Slope_VG(:,iR,jR,kR)
                           end do
                        end do
                     end do
                  end if

               else
                  iBufferS = iBufferS_P(iProcRecv)

                  BufferS_I(            iBufferS+1) = iBlockRecv
                  BufferS_I(            iBufferS+2) = iRMin
                  BufferS_I(            iBufferS+3) = iRMax
                  if(nDim > 1)BufferS_I(iBufferS+4) = jRMin
                  if(nDim > 1)BufferS_I(iBufferS+5) = jRMax
                  if(nDim > 2)BufferS_I(iBufferS+6) = kRMin
                  if(nDim > 2)BufferS_I(iBufferS+7) = kRMax

                  iBufferS = iBufferS + 1 + 2*nDim
                  if(present(Time_B))then
                     iBufferS = iBufferS + 1
                     BufferS_I(iBufferS) = Time_B(iBlockSend)
                  end if

                  do kR = kRMin, kRMax, DkR
                     kS = kSMin + abs((kR+9)/kRatioRestr &
                          -           (kRMin+9)/kRatioRestr)
                     do jR=jRMin, jRMax, DjR
                        jS = jSMin + abs((jR+9)/jRatioRestr &
                             -           (jRMin+9)/jRatioRestr)
                        do iR = iRMin, iRMax, DiR
                           iS = iSMin + abs((iR+9)/iRatioRestr &
                                -           (iRMin+9)/iRatioRestr)
                           BufferS_I(iBufferS+1:iBufferS+nVar)= &
                                State_VGB(:,iS,jS,kS,iBlockSend) &
                                + Slope_VG(:,iR,jR,kR)
                           iBufferS = iBufferS + nVar
                        end do
                     end do
                  end do

                  iBufferS_P(iProcRecv) = iBufferS

               end if
            end do
         end do
      end do

    end subroutine do_prolong

    !==========================================================================

    subroutine set_range

      integer:: nWidthProlongS_D(MaxDim), iDim
      !------------------------------------------------------------------------

      ! Indexed by iDir/jDir/kDir for sender = -1,0,1
      iEqualS_DII(:,-1,Min_) = 1
      iEqualS_DII(:,-1,Max_) = nWidth
      iEqualS_DII(:, 0,Min_) = 1
      iEqualS_DII(:, 0,Max_) = nIjk_D
      iEqualS_DII(:, 1,Min_) = nIjk_D + 1 - nWidth
      iEqualS_DII(:, 1,Max_) = nIjk_D

      ! Indexed by iDir/jDir/kDir for sender = -1,0,1
      iEqualR_DII(:,-1,Min_) = nIjk_D + 1
      iEqualR_DII(:,-1,Max_) = nIjk_D + nWidth
      iEqualR_DII(:, 0,Min_) = 1
      iEqualR_DII(:, 0,Max_) = nIjk_D
      iEqualR_DII(:, 1,Min_) = 1 - nWidth
      iEqualR_DII(:, 1,Max_) = 0

      ! Indexed by iDir/jDir/kDir for sender = -1,0,1
      iRestrictS_DII(:,-1,Min_) = 1
      iRestrictS_DII(:, 0,Min_) = 1
      iRestrictS_DII(:, 0,Max_) = nIjk_D
      iRestrictS_DII(:, 1,Max_) = nIjk_D
      if(DoRestrictFace)then
         iRestrictS_DII(:,-1,Max_) = nWidth
         iRestrictS_DII(:, 1,Min_) = nIjk_D + 1 - nWidth
      else
         iRestrictS_DII(:,-1,Max_) = iRatio_D*nWidth
         iRestrictS_DII(:, 1,Min_) = nIjk_D + 1 - iRatio_D*nWidth
      end if

      ! Indexed by iRecv/jRecv/kRecv = 0..3
      iRestrictR_DII(:,0,Min_) = 1 - nWidth
      iRestrictR_DII(:,0,Max_) = 0
      iRestrictR_DII(:,1,Min_) = 1
      do iDim = 1, MaxDim
         ! This loop is used to avoid the NAG 5.1 (282) bug on nyx
         iRestrictR_DII(iDim,1,Max_) = nIjk_D(iDim)/iRatio_D(iDim)
         iRestrictR_DII(iDim,2,Min_) = nIjk_D(iDim)/iRatio_D(iDim) + 1
      end do
      iRestrictR_DII(:,2,Max_) = nIjk_D
      iRestrictR_DII(:,3,Min_) = nIjk_D + 1
      iRestrictR_DII(:,3,Max_) = nIjk_D + nWidth

      ! Number of ghost cells sent from coarse block.
      ! Divided by resolution ratio and rounded up.
      nWidthProlongS_D         = 0
      if(nCoarseLayer == 1)then
         nWidthProlongS_D(1:nDim) = 1 + (nWidth-1)/iRatio_D(1:nDim)
      else
         nWidthProlongS_D(1:nDim) = nWidth
      end if

      ! Indexed by iSend/jSend,kSend = 0..3
      do iDim = 1, MaxDim
         ! This loop is used to avoid the NAG 5.1 (282) bug on nyx
         iProlongS_DII(iDim,0,Min_) = 1
         iProlongS_DII(iDim,0,Max_) = nWidthProlongS_D(iDim)
         iProlongS_DII(iDim,1,Min_) = 1
         iProlongS_DII(iDim,1,Max_) = nIjk_D(iDim)/iRatio_D(iDim)
         iProlongS_DII(iDim,2,Min_) = nIjk_D(iDim)/iRatio_D(iDim) + 1
         iProlongS_DII(iDim,2,Max_) = nIjk_D(iDim)
         iProlongS_DII(iDim,3,Min_) = nIjk_D(iDim) + 1 - nWidthProlongS_D(iDim)
         iProlongS_DII(iDim,3,Max_) = nIjk_D(iDim)
      end do

      ! Indexed by iRecv/jRecv/kRecv = 0,1,2,3
      iProlongR_DII(:, 0,Min_) = 1 - nWidth
      iProlongR_DII(:, 0,Max_) = 0
      iProlongR_DII(:, 1,Min_) = 1
      iProlongR_DII(:, 1,Max_) = nIjk_D
      iProlongR_DII(:, 2,Min_) = 1
      iProlongR_DII(:, 2,Max_) = nIjk_D
      iProlongR_DII(:, 3,Min_) = nIjk_D + 1
      iProlongR_DII(:, 3,Max_) = nIjk_D + nWidth

      if(DoSendCorner)then
         ! Face + two edges + corner or edge + one corner 
         ! are sent/recv together from fine to coarse block

         do iDim = 1, nDim
            if(iRatio_D(iDim) == 1)CYCLE

            ! The extension is by nWidth/2 rounded upwards independent of
            ! the value of nCoarseLayers. There is no need to send 
            ! two coarse layers into corner/edge ghost cells.

            iProlongS_DII(iDim,1,Max_) = iProlongS_DII(iDim,1,Max_) &
                 + (nWidth+1)/2
            iProlongS_DII(iDim,2,Min_) = iProlongS_DII(iDim,2,Min_) &
                 - (nWidth+1)/2
            iProlongR_DII(iDim,1,Max_) = iProlongR_DII(iDim,1,Max_) + nWidth
            iProlongR_DII(iDim,2,Min_) = iProlongR_DII(iDim,2,Min_) - nWidth
         end do
      end if

    end subroutine set_range

  end subroutine message_pass_real

  !============================================================================

  subroutine test_pass_cell

    use BATL_mpi,  ONLY: iProc, iComm
    use BATL_size, ONLY: MaxDim, nDim, nDimAmr, iRatio, jRatio, kRatio, &
         MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nG, nI, nJ, nK, nBlock,&
         nIJK_D, iRatio_D
    use BATL_tree, ONLY: init_tree, set_tree_root, find_tree_node, &
         refine_tree_node, distribute_tree, show_tree, clean_tree, &
         Unused_B, DiLevelNei_IIIB, iNode_B
    use BATL_grid, ONLY: init_grid, create_grid, clean_grid, &
         Xyz_DGB, CellSize_DB, CoordMin_DB
    use BATL_geometry, ONLY: init_geometry, z_, IsPeriodic_D

    use ModMpi, ONLY: MPI_allreduce, MPI_REAL, MPI_MIN, MPI_MAX

    integer, parameter:: MaxBlockTest            = 50
    logical:: IsPeriodicTest_D(MaxDim)= .true.
    integer:: nRootTest_D(MaxDim) = (/3,3,3/)
    real   :: DomainMin_D(MaxDim) = (/ 1.0, 10.0, 100.0 /)
    real   :: DomainMax_D(MaxDim) = (/ 4.0, 40.0, 400.0 /)
    real   :: DomainSize_D(MaxDim)

    real   :: Tolerance = 1e-6

    integer, parameter:: nVar = nDim
    real, allocatable:: State_VGB(:,:,:,:,:)
    real, allocatable:: Scalar_GB(:,:,:,:)
    real, allocatable:: FineGridLocal_III(:,:,:)
    real, allocatable:: FineGridGlobal_III(:,:,:)
    real, allocatable:: XyzCorn_DGB(:,:,:,:,:)
    real :: CourseGridCell_III(iRatio,jRatio,kRatio)

    integer:: nWidth
    integer:: nProlongOrder
    integer:: nCoarseLayer
    integer:: iSendCorner,  iRestrictFace
    logical:: DoSendCorner, DoRestrictFace

    real:: Xyz_D(MaxDim)
    integer:: iNode, iBlock, i, j, k, iMin, iMax, jMin, jMax, kMin, kMax, iDim
    integer:: iDir, jDir, kDir, Di, Dj, Dk

    integer ::iOp
    integer, parameter :: nOp=2
    character(len=4) :: NameOperator_I(nOp) = (/ "min", "max" /)
    character(len=4) :: NameOperator = "Min"
    real :: FineGridStep_D(MaxDim)
    integer :: iFG, jFG, kFG
    integer :: nFineCell
    integer :: iMpiOperator
    integer :: iError, iTest

    character(len=20):: NameGeometry

    logical:: DoTestMe
    character(len=*), parameter :: NameSub = 'test_pass_cell'
    !-----------------------------------------------------------------------
    DoTestMe = iProc == 0

    if(DoTestMe) write(*,*) 'Starting ',NameSub
    DomainSize_D = DomainMax_D - DomainMin_D
    call init_tree(MaxBlockTest)
    call init_geometry( IsPeriodicIn_D = IsPeriodicTest_D(1:nDim) )
    call init_grid( DomainMin_D(1:nDim), DomainMax_D(1:nDim) )
    call set_tree_root( nRootTest_D(1:nDim))

    call find_tree_node( (/0.5,0.5,0.5/), iNode)
    if(DoTestMe)write(*,*) NameSub,' middle node=',iNode
    call refine_tree_node(iNode)
    call distribute_tree(.true.)
    call create_grid

    if(DoTestMe) call show_tree(NameSub,.true.)

    allocate(State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlockTest))

    do nProlongOrder = 1, 2; do nCoarseLayer = 1, 2; do nWidth = 1, nG

       ! Second order prolongation does not work with sending multiple coarse 
       ! cell layers into the fine cells with their original values. 
       if(nProlongOrder == 2 .and. nCoarseLayer == 2) CYCLE

       ! Cannot send more coarse layers than the number of ghost cell layers
       if(nCoarseLayer > nWidth) CYCLE

       if(DoTestMe)write(*,*) 'testing message_pass_cell with', &
            ' nProlongOrder=',  nProlongOrder, &
            ' nCoarseLayer=',   nCoarseLayer,  &
            ' nWidth=',         nWidth

       ! Set the range of ghost cells that should be set
       iMin =  1 - nWidth
       jMin =  1; if(nDim > 1) jMin = 1 - nWidth
       kMin =  1; if(nDim > 2) kMin = 1 - nWidth
       iMax = nI + nWidth
       jMax = nJ; if(nDim > 1) jMax = nJ + nWidth
       kMax = nK; if(nDim > 2) kMax = nK + nWidth

       do iSendCorner = 1, 2; do iRestrictFace = 1, 2

          DoSendCorner   = iSendCorner   == 2
          DoRestrictFace = iRestrictFace == 2

          ! Second order prolongation does not work with restricting face:
          ! the first order restricted cell cannot be used in the prolongation.
          if(DoRestrictFace .and. nProlongOrder == 2) CYCLE

          if(DoTestMe)write(*,*) 'testing message_pass_cell with', &
               ' DoSendCorner=',   DoSendCorner, &
               ' DoRestrictFace=', DoRestrictFace

          State_VGB = 0.0

          do iBlock = 1, nBlock
             if(Unused_B(iBlock)) CYCLE
             State_VGB(:,1:nI,1:nJ,1:nK,iBlock) = &
                  Xyz_DGB(1:nDim,1:nI,1:nJ,1:nK,iBlock)
          end do

          call message_pass_cell(nVar, State_VGB, &
               nProlongOrderIn =nProlongOrder,    &
               nCoarseLayerIn  =nCoarseLayer,     &
               nWidthIn        =nWidth,           &
               DoSendCornerIn  =DoSendCorner,     &
               DoRestrictFaceIn=DoRestrictFace)

          do iBlock = 1, nBlock
             if(Unused_B(iBlock)) CYCLE

             ! Loop through all cells including ghost cells
             do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI

                ! The filled in second order accurate ghost cell value 
                ! should be the same as the coordinates of the cell center
                Xyz_D = Xyz_DGB(:,i,j,k,iBlock)

                ! Check that no info is sent in the non-used dimensions,
                ! i.e. for all iDim: nDim+1 < iDim < MaxDim
                if(  i < iMin .or. i > iMax .or. &
                     j < jMin .or. j > jMax .or. &
                     k < kMin .or. k > kMax) then
                   do iDim = 1, nDim
                      if(abs(State_VGB(iDim,i,j,k,iBlock)) > 1e-6)then
                         write(*,*)'Face should not be set: ', &
                              'iProc,iBlock,i,j,k,iDim,State,Xyz=', &
                              iProc,iBlock,i,j,k,iDim, &
                              State_VGB(iDim,i,j,k,iBlock), &
                              Xyz_D(iDim)
                      end if
                   end do

                   CYCLE
                end if

                ! Get the direction vector
                iDir = 0; if(i<1) iDir = -1; if(i>nI) iDir = 1
                jDir = 0; if(j<1) jDir = -1; if(j>nJ) jDir = 1
                kDir = 0; if(k<1) kDir = -1; if(k>nK) kDir = 1

                ! If nCoarseLayer==2 and DoSendCorner is true
                ! the second ghost cells in the corner/edges
                ! are not well defined (they may contain
                ! the value coming from the first or second coarse cell).

                if(nCoarseLayer==2 .and. DoSendCorner .and. ( &
                     (i<0 .or. i>nI+1) .and. (jDir /= 0 .or. kDir /= 0) .or. &
                     (j<0 .or. j>nJ+1) .and. (iDir /= 0 .or. kDir /= 0) .or. &
                     (k<0 .or. k>nK+1) .and. (iDir /= 0 .or. jDir /= 0) &
                     )) CYCLE

                ! if we do not send corners and edges, check that the
                ! State_VGB in these cells is still the unset value
                if(.not.DoSendCorner .and. ( &
                     iDir /= 0 .and. jDir /= 0 .or. &
                     iDir /= 0 .and. kDir /= 0 .or. &
                     jDir /= 0 .and. kDir /= 0 ))then

                   do iDim = 1, nDim
                      if(abs(State_VGB(iDim,i,j,k,iBlock)) > 1e-6)then
                         write(*,*)'corner/edge should not be set: ', &
                              'iProc,iBlock,i,j,k,iDim,State,Xyz=', &
                              iProc,iBlock,i,j,k,iDim, &
                              State_VGB(iDim,i,j,k,iBlock), &
                              Xyz_D
                      end if
                   end do

                   CYCLE
                end if

                ! Shift ghost cell coordinate into periodic domain
                Xyz_D = DomainMin_D + modulo(Xyz_D - DomainMin_D, DomainSize_D)

                ! Calculate distance of ghost cell layer
                Di = 0; Dj = 0; Dk = 0
                if(i <  1 .and. iRatio == 2) Di = 2*i-1
                if(i > nI .and. iRatio == 2) Di = 2*(i-nI)-1
                if(j <  1 .and. jRatio == 2) Dj = 2*j-1
                if(j > nJ .and. jRatio == 2) Dj = 2*(j-nJ)-1
                if(k <  1 .and. kRatio == 2) Dk = 2*k-1
                if(k > nK .and. kRatio == 2) Dk = 2*(k-nK)-1

                if(DoRestrictFace .and. &
                     DiLevelNei_IIIB(iDir,jDir,kDir,iBlock) == -1)then

                   ! Shift coordinates if only 1 layer of fine cells
                   ! is averaged in the orthogonal direction
                   Xyz_D(1) = Xyz_D(1) - 0.25*Di*CellSize_DB(1,iBlock)
                   Xyz_D(2) = Xyz_D(2) - 0.25*Dj*CellSize_DB(2,iBlock)
                   Xyz_D(3) = Xyz_D(3) - 0.25*Dk*CellSize_DB(3,iBlock)

                end if

                if(nProlongOrder == 1 .and. &
                     DiLevelNei_IIIB(iDir,jDir,kDir,iBlock) == 1)then

                   ! Shift depends on the parity of the fine ghost cell
                   ! except when there is no AMR or multiple coarse cell 
                   ! layers are sent in that direction
                   if(iRatio == 2 .and. (nCoarseLayer == 1 .or. iDir == 0)) &
                        Di = 2*modulo(i,2) - 1
                   if(jRatio == 2 .and. (nCoarseLayer == 1 .or. jDir == 0)) &
                        Dj = 2*modulo(j,2) - 1
                   if(kRatio == 2 .and. (nCoarseLayer == 1 .or. kDir == 0)) &
                        Dk = 2*modulo(k,2) - 1

                   Xyz_D(1) = Xyz_D(1) + 0.5*Di*CellSize_DB(1,iBlock)
                   Xyz_D(2) = Xyz_D(2) + 0.5*Dj*CellSize_DB(2,iBlock)
                   Xyz_D(3) = Xyz_D(3) + 0.5*Dk*CellSize_DB(3,iBlock)

                end if

                do iDim = 1, nDim
                   if(abs(State_VGB(iDim,i,j,k,iBlock) - Xyz_D(iDim)) &
                        > Tolerance)then
                      write(*,*)'iProc,iBlock,i,j,k,iDim,State,Xyz=', &
                           iProc,iBlock,i,j,k,iDim, &
                           State_VGB(iDim,i,j,k,iBlock), &
                           Xyz_D(iDim)
                   end if
                end do
             end do; end do; end do
          end do

       end do; end do; end do
    end do; end do ! test parameters
    deallocate(State_VGB)

    call clean_grid
    call clean_tree
    !------------------------ Test Scalar -----------------------------

    ! To test the message pass for the cell with min max operators we 
    ! generate a fine uniform grid 
    ! for the whole domain and transfer the cell values from the
    ! block cells to the cells on the fine grid. Then we gather all the
    ! data on the fine grid with the proper operator.
    ! We can then compare the values on the coresponding node after
    ! message_pass_cell_scalar is called with the fine grid values.

    ! rescale the domain to make indexing easyer
    DomainSize_D = iRatio_D*nRootTest_D*nIJK_D
    DomainMin_D = (/ 0.0, 0.0, 0.0 /)
    DomainMax_D = DomainSize_D

    call init_tree(MaxBlockTest)
    call init_grid( DomainMin_D(1:nDim), DomainMax_D(1:nDim) )
    call init_geometry( IsPeriodicIn_D = IsPeriodicTest_D(1:nDim) )
    call set_tree_root( nRootTest_D(1:nDim))

    call find_tree_node( (/0.5,0.5,0.5/), iNode)
    call refine_tree_node(iNode)
    call distribute_tree(.true.)
    call create_grid

    ! Position of cell corners, for solving problems with round-off 
    ! when getting fine grid positions
    allocate(XyzCorn_DGB(MaxDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlockTest))
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI 
          XyzCorn_DGB(:,i,j,k,iBlock) = &
               Xyz_DGB(:,i,j,k,iBlock) - &
               0.5*CellSize_DB(:,iBlock)*&
               (/ min(1,nI-1),min(1,nJ-1),min(1,nK-1) /)
       end do; end do; end do
    end do

    allocate(Scalar_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlockTest))
    Scalar_GB = -7777

    allocate(FineGridLocal_III( &
         nI*iRatio*nRootTest_D(1),&
         nJ*jRatio*nRootTest_D(2),&
         nK*kRatio*nRootTest_D(3)))
    allocate(FineGridGlobal_III( &
         (nI)*iRatio*nRootTest_D(1),&
         (nJ)*jRatio*nRootTest_D(2),&
         (nK)*kRatio*nRootTest_D(3)))


    nFineCell = ((nI)*iRatio*nRootTest_D(1))*&
         ((nJ)*jRatio*nRootTest_D(2))*&
         ((nK)*kRatio*nRootTest_D(3))

    FineGridStep_D = DomainSize_D &
         / (DomainMax_D - DomainMin_D)

    do iOp = 1, nOp

       NameOperator = NameOperator_I(iOp)
       select case(NameOperator)
       case("min")
          FineGridLocal_III(:,:,:)  =  1.0e8
          FineGridGlobal_III(:,:,:) =  1.0e8
          iMpiOperator = MPI_MIN
       case("max")
          FineGridLocal_III(:,:,:)  = -1.0e8
          FineGridGlobal_III(:,:,:) = -1.0e8 
          iMpiOperator=MPI_MAX
       case default
          call CON_stop(NameSub//': incorrect operator name')
       end select

       if(DoTestMe) write(*,*) 'testing message_pass_cell_scalar ', &
            'with operator= ',NameOperator

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             Scalar_GB(i,j,k,iBlock)= iNode_B(iBlock) + &
                  sum(CoordMin_DB(:,iBlock) + &
                  ( (/i, j, k/) ) * CellSize_DB(:,iBlock))
          end do; end do; end do
       end do

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          do k = 1, nK; do j = 1, nJ; do i = 1, nI

             iFG = nint(XyzCorn_DGB(1,i,j,k,iBlock)*FineGridStep_D(1)) + 1 
             jFG = nint(XyzCorn_DGB(2,i,j,k,iBlock)*FineGridStep_D(2)) + 1
             kFG = nint(XyzCorn_DGB(3,i,j,k,iBlock)*FineGridStep_D(3)) + 1
    
             FineGridLocal_III(iFG,jFG,kFG) = Scalar_GB(i,j,k,iBlock)
          end do; end do; end do
       end do

       call message_pass_cell(Scalar_GB, &
            nProlongOrderIn=1, nCoarseLayerIn=2, &
            DoSendCornerIn=.true., DoRestrictFaceIn=.false., &
            NameOperatorIn=NameOperator_I(iOp))

       call MPI_ALLREDUCE(FineGridLocal_III(1,1,1), &
            FineGridGlobal_III(1,1,1),              &
            nFineCell, MPI_REAL, iMpiOperator, iComm, iError)

       ! making sure that we have the center cell along the x=0 side
       ! so the boundary are not tested.
       call find_tree_node( (/0.0,0.5,0.5/), iNode)

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          if(iNode_B(iBlock) == iNode) then
             do k = MinK, MaxK; do j = MinJ, MaxJ; do i = 1, MaxI

                iFG = nint(XyzCorn_DGB(1,i,j,k,iBlock)*FineGridStep_D(1)) + 1 
                jFG = nint(XyzCorn_DGB(2,i,j,k,iBlock)*FineGridStep_D(2)) + 1
                kFG = nint(XyzCorn_DGB(3,i,j,k,iBlock)*FineGridStep_D(3)) + 1

                !copy cells that are inside the course grid cell
                CourseGridCell_III = FineGridGlobal_III(&
                     iFG:iFG+min(1,iRatio-1),&
                     jFG:jFG+min(1,jRAtio-1),&
                     kFG:kFG+min(1,kRatio-1))

                select case(NameOperator)
                case("min") 
                   if(Scalar_GB(i,j,k,iBlock) /= &
                        minval(CourseGridCell_III))then
                      write (*,*) "Error for operator, iNode,  iBlock= ",&
                           NameOperator, iNode_B(iBlock), iBlock, ", value=",&
                           minval(CourseGridCell_III),&
                           " should be ", Scalar_GB(i,j,k,iBlock), "index : " &
                           ,i,j,k, " ::", iFG, jFG,kFG

                   end if
                case("max") 
                   if(Scalar_GB(i,j,k,iBlock) /= &
                        maxval(CourseGridCell_III))then
                      write (*,*) "Error for operator, iNode,  iBlock= ",&
                           NameOperator, iNode_B(iBlock), iBlock, ", value=",&
                           maxval(CourseGridCell_III),&
                           " should be ", Scalar_GB(i,j,k,iBlock), "index : " &
                           ,i,j,k, " ::", iFG, jFG,kFG
                   end if
                end select

             end do; end do; end do
          end if
       end do

    end do
    deallocate(Scalar_GB, FineGridLocal_III, FineGridGlobal_III,XyzCorn_DGB)
    call clean_grid
    call clean_tree

    if(nDim == 1) RETURN !------------------------

    do iTest = 1,6

       ! The code is quite inaccurate for partial AMR across the pole
       if(nDimAmr < nDim .and. iTest > 3) EXIT

       call init_tree(MaxBlockTest)

       ! Do not test ghost cells in the radial direction
       iMin = 1; iMax = nI
       jMin = MinJ; jMax = MaxJ
       kMin = MinK; kMax = MaxK

       select case(iTest)
       case(1,4)
          NameGeometry = 'cylindrical'

          ! 0 < r < 10, 0 < phi < 360deg, -5 < z < 5
          DomainMin_D = (/ 0.0,  0.0, -5.0 /)
          DomainMax_D = (/ 8.0, cTwoPi, +5.0 /)
          IsPeriodicTest_D = (/ .false., .true., .true. /)

          ! There must be an even number of root blocks in the phi direction
          ! There are 3 root blocks in z so that we can refine the middle
          ! and avoid issues of periodicity in the testing
          nRootTest_D = (/ 2, 4, 3 /)

          ! Test ghost cells at rMin
          iMin = MinI

       case(2,5)
          if(nDim < 3)CYCLE
          NameGeometry = 'spherical'

          ! 1 < r < 9, 0 < theta < 180deg, 0 < phi < 360deg
          DomainMin_D = (/ 1.0,  0.0, 0.0 /)
          DomainMax_D = (/ 9.0,  cPi, cTwoPi /)
          IsPeriodicTest_D = (/ .false., .false., .true. /)

          ! There must be an even number of root blocks in the phi direction
          ! There are 3 root blocks in r so that we can refine the middle
          ! and avoid issues at inner and outer radial boundaries
          nRootTest_D = (/ 3, 2, 4 /)

       case(3,6)
          if(nDim < 3)CYCLE
          NameGeometry = 'rlonlat'

          ! 1 < r < 9, 0 < phi < 360deg, -90 < lat < 90
          DomainMin_D = (/ 1.0,    0.0, -cHalfPi /)
          DomainMax_D = (/ 9.0, cTwoPi, cHalfPi /)
          IsPeriodicTest_D = (/ .false., .true., .false. /)

          ! There must be an even number of root blocks in the phi direction
          ! There are 3 root blocks in r so that we can refine the middle
          ! and avoid issues at inner and outer radial boundaries
          nRootTest_D = (/ 3, 4, 2 /)

       end select
       DomainSize_D = DomainMax_D - DomainMin_D

       if(DoTestMe)then
          if(iTest <= 3)write(*,*) &
               'testing message_pass_cell across '//trim(NameGeometry)//' pole'
          if(iTest >= 4)write(*,*) &
               'testing message_pass_cell across '//trim(NameGeometry)// &
               ' pole with resolution change'
       end if

       call init_geometry(NameGeometry, &
            IsPeriodicIn_D=IsPeriodicTest_D(1:nDim))

       call init_grid(DomainMin_D(1:nDim), DomainMax_D(1:nDim), &
            UseDegreeIn=.false.)
       call set_tree_root( nRootTest_D(1:nDim))

       if(any(IsPeriodic_D(1:nDim) .neqv. IsPeriodicTest_D(1:nDim))) &
            write(*,*) NameSub,': IsPeriodic_D=', IsPeriodic_D(1:nDim), &
            ' should agree with ', IsPeriodicTest_D(1:nDim)


       if(iTest > 3)then
          ! Test with refined grid
          if(iTest==4)then
             ! refine node next to r=0 axis but middle in Z direction
             call refine_tree_node(9)
          else
             ! refine root nodes at min and max (theta/lat/phi) coordinates
             ! but middle in the R direction
             call refine_tree_node(2)
             call refine_tree_node(23)
          end if
          ! Restriction is not linear so there is truncation error
          Tolerance = 0.15
       else
          ! For tests with no AMR the error is round-off only
          Tolerance = 1e-6
       end if

       call distribute_tree(.true.)
       call create_grid

       !if(DoTestMe) call show_tree(NameSub,.true.)
       !do iBlock = 3, 3
       !   do j = 0, 3; do i = 0, 3
       !      write(*,*) '!!! iBlock, i, j, iNodeNei_IIIB(i,j,:,iBlock)=', &
       !           iBlock, i, j, iNodeNei_IIIB(i,j,:,iBlock)
       !   end do; end do
       !end do

       allocate(State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlockTest))

       State_VGB = 0.0
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          State_VGB(:,1:nI,1:nJ,1:nK,iBlock) = &
               Xyz_DGB(1:nDim,1:nI,1:nJ,1:nK,iBlock)
       end do

       ! Second order
       call message_pass_cell(nVar, State_VGB)

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          ! Loop through all cells including ghost cells
          do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax

             ! The filled in second order accurate ghost cell value 
             ! should be the same as the coordinates of the cell center
             Xyz_D = Xyz_DGB(:,i,j,k,iBlock)

             ! For 3D cylindrical Z coordinate is periodic
             if( (iTest==1 .or. iTest==4) .and. nDim == 3) &
                  Xyz_D(z_) = DomainMin_D(z_) &
                  + modulo(Xyz_D(z_) - DomainMin_D(z_), DomainSize_D(z_))

             do iDim = 1, nDim
                if(abs(State_VGB(iDim,i,j,k,iBlock) - Xyz_D(iDim)) &
                     /abs( Xyz_D(iDim)) > Tolerance)then
                   write(*,*)'iProc,iBlock,i,j,k,iDim,State,Xyz=', &
                        iProc,iBlock,i,j,k,iDim, &
                        State_VGB(iDim,i,j,k,iBlock), &
                        Xyz_D(iDim)

                end if
             end do
          end do; end do; end do
       end do

       deallocate(State_VGB)

       call clean_grid
       call clean_tree
    end do

  end subroutine test_pass_cell

end module BATL_pass_cell
