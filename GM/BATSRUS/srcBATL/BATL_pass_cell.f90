!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_pass_cell

  use BATL_geometry, ONLY: IsCartesianGrid, IsRotatedCartesian, IsRoundCube, &
       IsCylindricalAxis, IsSphericalAxis, IsLatitudeAxis, Lat_, Theta_
  use ModNumConst, ONLY: cPi, cHalfPi, cTwoPi
  use BATL_high_order, ONLY: restriction_high_order_reschange
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
         iTree_IA, Proc_, Block_, Coord1_, Coord2_, Coord3_, find_neighbor_for_anynode

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
    ! for  high order resolution change. 
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

    ! Variables for coarsened block. 
    real, allocatable:: State_VIIIB(:,:,:,:,:)
    logical, allocatable:: IsAccurate_B(:)

    logical:: Do6thCorrect = .false. 

    ! For high order resolution change, a few face ghost cells need to be 
    ! calculated remotely after the corase block have got accurate 
    ! ghost cells. For this case with block cell size smaller than 
    ! 8, these few inaccurate face ghost cells need to be sent to 
    ! neighbour block before they are overwritten by 5th order values
    ! calculated remotely.  So, nSubStage == 2 for this case. Do_equal 
    ! is called first and then use do_prolong. 
    logical:: DoSendFace = .false., DoRecvFace = .false. 
    logical, allocatable:: IsAccurateFace_GB(:,:,:,:)
    integer:: iSubStage, nSubStage 
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

    if(UseHighResChange .and. nProlongOrder /=1) call CON_stop(NameSub// &
         'nProlongOrder should be 1 for high order resolution change')

    if(UseHighResChange .and. .not. DoSendCorner) call CON_stop(NameSub// &
         'DoSendCorner should be true when UseHighResChange = .true.')

    if(UseHighResChange .and. nCoarseLayer == 1) call CON_stop(NameSub// &
         'nCoarseLayer should be 2 when UseHighResChange = .true.')

    if(UseHighResChange) DoRestrictFace = .false.

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
    if(.not. allocated(IsAccurateFace_GB)) &
         allocate(IsAccurateFace_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK, nBlock))
    if(UseHighResChange) then

       ! stage 1: first order prolongation and do_equal.
       ! stage 2: a) Do high order restriction remotely for all kinds of ghost
       !          cells. b) Do high order prolongation for face ghost cell 
       !          locally.  c) Pass restricted cells from fine to coarse.
       ! stage 3: a) Refine with high order accuracy on coarse blocks for 
       !          edges and corners and for faces that are too complex to do
       !          locally. b)  Pass locally high order prolonged face ghost cells 
       !          to edge/corner ghost cells of neighboring fine block.
       !          c) Pass remotely high order refined ghost cells from coarse
       !          to fine block for other edge and corner ghost cells and 
       !          also for not yet done face ghost cells.
       nSendStage = 3

       ! For 6th order correction, which may be better because of symmetry,
       ! 8 cells are needed in each direction. If it is not satisfied, 
       ! use 5th order correction.        
       Do6thCorrect = nI>7 .and. nJ>7 .and. (nDim==2 .or. nK>7)

       ! Used for stage 2a. 
       if(.not. allocated(State_VIIIB))&
            allocate(&
            State_VIIIB(nVar,1:max(nI/2,1),1:max(nJ/2,1),1:max(nK/2,1),&
            nBlock), IsAccurate_B(nBlock))

    endif

    do iSendStage = 1, nSendStage
       if(UseHighResChange) then
          State_VIIIB = 0
          IsAccurate_B = .false. 
       endif

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

          nSubStage = 1
          if(iSendStage == 3) nSubStage = 2

          do iSubStage = 1, nSubStage
             ! For the last stage of high order resolution change, 
             ! do_equal first, then do_prolong. 

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
                      if(.not.DoSendCorner .and. jDir /= 0 .and. kDir /= 0) &
                           CYCLE

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

                         ! Do prolongation in the second stage if 
                         ! nProlongOrder=2. We still need to call restriction 
                         ! and prolongation in both stages to calculate the
                         ! amount of received data
                         if(iSendStage == 2 .and. DiLevel == 0) CYCLE

                         ! Fill in the edge/corner ghost cells with values from
                         ! face ghost cells first (do_equal). Replace the 
                         ! inaccurate face ghost cells with the values from 
                         ! remote restriction. 
                         if(iSendStage == 3 .and. DiLevel /= 0 &
                              .and. iSubStage == 1) CYCLE
                         if(iSendStage == 3 .and. DiLevel ==0 &
                              .and. iSubStage == 2) CYCLE

                         if(DiLevel == 0)then
                            if(iSendStage == 3) then
                               call corrected_do_equal
                            else
                               if(.not.DoResChangeOnly) call do_equal
                            endif
                         elseif(DiLevel == 1)then
                            call do_restrict
                         elseif(DiLevel == -1)then
                            call do_prolong
                         endif
                      end do ! iDir
                   end do ! jDir
                end do ! kDir
             end do ! iBlockSend
          end do ! iSubStage

          call timing_stop('local_pass')

       end do ! iCountOnly

       ! Done for serial run
       if(nProc == 1) then
          if(UseHighResChange .and. iSendStage == 2) then
             do iBlock = 1, nBlock
                if (Unused_B(iBlock)) CYCLE
                call high_prolong_for_face_ghost(iBlock)                
             enddo
          endif
          ! Done for serial run
          CYCLE
       endif

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

       if(UseHighResChange .and. iSendStage == 2) then
          do iBlock = 1, nBlock
             if (Unused_B(iBlock)) CYCLE
             call high_prolong_for_face_ghost(iBlock)
          enddo
       endif
    end do ! iSendStage

    deallocate(Slope_VG)

    call timing_stop('batl_pass')

  contains
    !======================================================================

    logical function only_corner_fine(iNode, iDir0, jDir0, kDir0, &
         iDirCorner, jDirCorner, kDirCorner)
      integer, intent(in):: iNode, iDir0, jDir0, kDir0
      integer, optional, intent(inout):: iDirCorner,jDirCorner,kDirCorner
      integer:: DiLevelNei_III(-1:1,-1:1,-1:1)
      integer:: iDir1, jDir1, kDir1, iDir2, jDir2, kDir2
      integer:: iDirBegin, iDirEnd, jDirBegin, jDirEnd, kDirBegin, kDirEnd
      integer:: DiDir, DjDir, DkDir
      integer:: nRefinedEdge
      logical:: IsOnlyCornerFine = .false. 
      !----------------------------------------------------------------------

      call find_neighbor_for_anynode(iNode,DiLevelNei_III)

      IsOnlyCornerFine = .false. 
      if(abs(iDir0) + abs(jDir0) + abs(kDir0) == 1) then        

         ! Loop through 4 corners block corresponding to this face, 
         ! check whether it is finer. 
         iDirBegin = -1; iDirEnd = 1; DiDir = 2
         jDirBegin = -1; jDirEnd = 1; DjDir = 2
         kDirBegin = -1; kDirEnd = 1; DkDir = 2

         if(iDir0 /= 0) then
            iDirBegin = iDir0; iDirEnd = iDirBegin; DiDir = 1
         elseif(jDir0 /= 0) then
            jDirBegin = jDir0; jDirEnd = jDirBegin; DjDir = 1
         elseif(kDir0 /= 0) then
            kDirBegin = kDir0; kDirEnd = kDirBegin; DkDir = 1
         endif

         ! 4 corners
         do kDir1=kDirBegin,kDirEnd,DkDir; do jDir1=jDirBegin,jDirEnd,DjDir; &
              do iDir1=iDirBegin,iDirEnd,DiDir            
            if(DiLevelNei_III(iDir1,jDir1,kDir1) == 0) then
               nRefinedEdge = 0

               ! Check first edge block. 
               if(iDir0 /= 0) then
                  iDir2 = iDir1
                  jDir2 = jDir1
                  kDir2 = 0
               endif
               if(jDir0 /= 0) then
                  jDir2 = jDir1
                  kDir2 = kDir1
                  iDir2 = 0
               endif
               if(kDir0 /= 0) then
                  kDir2 = kDir1
                  iDir2 = iDir1
                  jDir2 = 0
               endif
               if(DiLevelNei_III(iDir2,jDir2,kDir2) == 0) &
                    nRefinedEdge = nRefinedEdge + 1

               ! Check second edge block. 
               if(iDir0 /= 0) then
                  iDir2 = iDir1
                  jDir2 = 0
                  kDir2 = kDir1
               endif
               if(jDir0 /= 0) then
                  jDir2 = jDir1
                  kDir2 = 0
                  iDir2 = iDir1
               endif
               if(kDir0 /= 0) then
                  kDir2 = kDir1
                  iDir2 = 0
                  jDir2 = jDir1
               endif
               if(DiLevelNei_III(iDir2,jDir2,kDir2) == 0) &
                    nRefinedEdge = nRefinedEdge + 1


               if(nRefinedEdge == 0) then
                  IsOnlyCornerFine = .true.
                  if(present(iDirCorner)) iDirCorner = iDir1
                  if(present(jDirCorner)) jDirCorner = jDir1
                  if(present(kDirCorner)) kDirCorner = kDir1
               endif
            endif
         enddo; enddo; enddo
      endif
      only_corner_fine = IsOnlyCornerFine
    end function only_corner_fine

    !======================================================================

    subroutine corrected_do_equal
      integer:: iEqualSOrig_DII(MaxDim,-1:1,Min_:Max_)
      integer:: iEqualROrig_DII(MaxDim,-1:1,Min_:Max_)
      integer:: iDir1,jDir1, kDir1, nDir
      integer:: iDir2, jDir2, kDir2, iDir3,jDir3,kDir3
      !----------------------------------------------------------------------

      nDir = abs(iDir)+abs(jDir)+abs(kDir)

      if(nDir > nDim-1) RETURN

      if(nDim == 2) then
         kDir1 = 0
         if(iDir /=0) then
            iDir1 = 0
            do jDir1 = -1, 1, 2
               ! Some information passed here is useless. Advantage: do not need 
               ! to change do_equal.
               if(DiLevelNei_IIIB(iDir1,jDir1,kDir1,iBlockSend) == 1 .or.&
                    DiLevelNei_IIIB(iDir, jDir1,kDir1,iBlockSend) == 1 ) then
                  iEqualSOrig_DII = iEqualS_DII
                  iEqualROrig_DII = iEqualR_DII

                  if(jDir1 == -1) then
                     iEqualS_DII(2,0,Min_) = 1-nWidth
                     iEqualS_DII(2,0,Max_) = 0

                     iEqualR_DII(2,0,Min_) = 1-nWidth
                     iEqualR_DII(2,0,Max_) = 0
                  elseif(jDir1 == 1) then
                     iEqualS_DII(2,0,Min_) = nJ + 1
                     iEqualS_DII(2,0,Max_) = nJ + nWidth

                     iEqualR_DII(2,0,Min_) = nJ + 1
                     iEqualR_DII(2,0,Max_) = nJ + nWidth
                  endif
                  call do_equal               
                  iEqualS_DII = iEqualSOrig_DII
                  iEqualR_DII = iEqualROrig_DII
               endif
            enddo ! jDir1

         elseif(jDir /=0) then
            jDir1 = 0
            do iDir1 = -1, 1, 2
               if(DiLevelNei_IIIB(iDir1,jDir1,kDir1,iBlockSend) == 1 .or. &
                    DiLevelNei_IIIB(iDir1,jDir, kDir1,iBlockSend) == 1) then
                  iEqualSOrig_DII = iEqualS_DII
                  iEqualROrig_DII = iEqualR_DII
                  if(iDir1 == -1) then
                     iEqualS_DII(1,0,Min_) = 1-nWidth
                     iEqualS_DII(1,0,Max_) = 0

                     iEqualR_DII(1,0,Min_) = 1-nWidth
                     iEqualR_DII(1,0,Max_) = 0
                  elseif(iDir1 == 1) then
                     iEqualS_DII(1,0,Min_) = nI + 1
                     iEqualS_DII(1,0,Max_) = nI + nWidth

                     iEqualR_DII(1,0,Min_) = nI + 1
                     iEqualR_DII(1,0,Max_) = nI + nWidth
                  endif
                  call do_equal               
                  iEqualS_DII = iEqualSOrig_DII
                  iEqualR_DII = iEqualROrig_DII
               endif
            enddo ! iDir1
         endif

      elseif(nDim == 3) then
         if(nDir == 2) then
            if(iDir == 0) then
               jDir1 = 0; kDir1 = 0
               do iDir1 = -1, 1, 2
                  if(DiLevelNei_IIIB(iDir1,jDir1,kDir1,iBlockSend) == 1 .or.&
                       DiLevelNei_IIIB(iDir1,jDir,kdir,iBlockSend) == 1) then 


                     !-------------------
                     ! The face values of iBlockSend is not accurate. Do not 
                     ! send to iBlockRecv. The edge/corner cells of iBlockRecv
                     ! will be filled by do_prolongation.
                     iSend = (3*iDir + 3)/2
                     jSend = (3*jDir + 3)/2
                     kSend = (3*kDir + 3)/2
                     iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)

                     iDir2 = -100; jDir2 = -100; kDir2 = -100
                     iDir3 = -100; jDir3 = -100; kDir3 = -100
                     if(only_corner_fine(&
                          iNodeRecv,iDir1,0,0,iDir2,jDir2,kDir2) .or. &
                          only_corner_fine(&
                          iNode_B(iBlockSend),iDir1,0,0,iDir3,jDir3,kDir3))then

                        if((jDir == -jDir2 .and. kDir == -kDir2) .or. &
                             (jDir == jDir3 .and. kDir == kDir3))&
                             CYCLE                                             
                     endif
                     !-------------------
                     
                     iEqualSOrig_DII = iEqualS_DII
                     iEqualROrig_DII = iEqualR_DII

                     if(iDir1 == -1) then
                        iEqualS_DII(1,0,Min_) = 1-nWidth
                        iEqualS_DII(1,0,Max_) = 0

                        iEqualR_DII(1,0,Min_) = 1-nWidth
                        iEqualR_DII(1,0,Max_) = 0
                     elseif(iDir1 == 1) then
                        iEqualS_DII(1,0,Min_) = nI + 1
                        iEqualS_DII(1,0,Max_) = nI + nWidth

                        iEqualR_DII(1,0,Min_) = nI + 1
                        iEqualR_DII(1,0,Max_) = nI + nWidth
                     endif

                     call do_equal               
                     iEqualS_DII = iEqualSOrig_DII
                     iEqualR_DII = iEqualROrig_DII
                  endif

               enddo

            elseif(jDir == 0) then
               iDir1 = 0; kDir1 = 0
               do jDir1 = -1, 1, 2
                  if(DiLevelNei_IIIB(iDir1,jDir1,kDir1,iBlockSend) == 1 .or.&
                       DiLevelNei_IIIB(iDir,jDir1,kDir,iBlockSend) == 1) then

                     iSend = (3*iDir + 3)/2
                     jSend = (3*jDir + 3)/2
                     kSend = (3*kDir + 3)/2
                     iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)

                     iDir2 = -100; jDir2 = -100; kDir2 = -100
                     iDir3 = -100; jDir3 = -100; kDir3 = -100
                     if(only_corner_fine(&
                          iNodeRecv,0,jDir1,0,iDir2,jDir2,kDir2) .or. &
                          only_corner_fine(&
                          iNode_B(iBlockSend),0,jDir1,0,iDir3,jDir3,kDir3))then

                        if((iDir == -iDir2 .and. kDir == -kDir2) .or. &
                             (iDir == iDir3 .and. kDir == kDir3))&
                             CYCLE                                             
                     endif

                     iEqualSOrig_DII = iEqualS_DII
                     iEqualROrig_DII = iEqualR_DII

                     if(jDir1 == -1) then
                        iEqualS_DII(2,0,Min_) = 1-nWidth
                        iEqualS_DII(2,0,Max_) = 0

                        iEqualR_DII(2,0,Min_) = 1-nWidth
                        iEqualR_DII(2,0,Max_) = 0
                     elseif(jDir1 == 1) then
                        iEqualS_DII(2,0,Min_) = nJ + 1
                        iEqualS_DII(2,0,Max_) = nJ + nWidth

                        iEqualR_DII(2,0,Min_) = nJ + 1
                        iEqualR_DII(2,0,Max_) = nJ + nWidth
                     endif

                     call do_equal               
                     iEqualS_DII = iEqualSOrig_DII
                     iEqualR_DII = iEqualROrig_DII
                  endif
               enddo ! jDir1
            elseif(kDir == 0) then
               iDir1 = 0; jDir1 = 0

               do kDir1 = -1, 1, 2
                  if(DiLevelNei_IIIB(iDir1,jDir1,kDir1,iBlockSend) == 1 .or.&
                       DiLevelNei_IIIB(iDir,jDir,kDir1,iBlockSend) == 1) then

                     iSend = (3*iDir + 3)/2
                     jSend = (3*jDir + 3)/2
                     kSend = (3*kDir + 3)/2
                     iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)

                     iDir2 = -100; jDir2 = -100; kDir2 = -100
                     iDir3 = -100; jDir3 = -100; kDir3 = -100
                     if(only_corner_fine(&
                          iNodeRecv,0,0,kDir1,iDir2,jDir2,kDir2) .or. &
                          only_corner_fine(&
                          iNode_B(iBlockSend),0,0,kDir1,iDir3,jDir3,kDir3))then

                        if((jDir == -jDir2 .and. iDir == -iDir2) .or. &
                             (jDir == jDir3 .and. iDir == iDir3))&
                             CYCLE                                             
                     endif

                     iEqualSOrig_DII = iEqualS_DII
                     iEqualROrig_DII = iEqualR_DII

                     if(kDir1 == -1) then
                        iEqualS_DII(3,0,Min_) = 1-nWidth
                        iEqualS_DII(3,0,Max_) = 0

                        iEqualR_DII(3,0,Min_) = 1-nWidth
                        iEqualR_DII(3,0,Max_) = 0
                     elseif(kDir1 == 1) then
                        iEqualS_DII(3,0,Min_) = nK + 1
                        iEqualS_DII(3,0,Max_) = nK + nWidth

                        iEqualR_DII(3,0,Min_) = nK + 1
                        iEqualR_DII(3,0,Max_) = nK + nWidth
                     endif

                     call do_equal               

                     iEqualS_DII = iEqualSOrig_DII
                     iEqualR_DII = iEqualROrig_DII
                  endif
               enddo ! kDir1
            endif
         elseif(nDir == 1) then
            if(iDir /= 0) then
               iDir1 = 0
               do kDir1 = -1, 1; do jDir1 = -1, 1
                  if(kDir1 == 0 .and. jDir1 == 0) CYCLE
                  if(kDir1 /= 0 .and. jDir1 /= 0) CYCLE
                  if(DiLevelNei_IIIB(iDir1,jDir1,kDir1,iBlockSend) == 1 .or. &
                       DiLevelNei_IIIB(iDir,jDir1,kDir1,iBlockSend) == 1) then


                     iSend = (3*iDir + 3)/2
                     jSend = (3*jDir + 3)/2
                     kSend = (3*kDir + 3)/2
                     iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)

                     iDir2 = -100; jDir2 = -100; kDir2 = -100
                     iDir3 = -100; jDir3 = -100; kDir3 = -100

                     if(only_corner_fine(&
                          iNodeRecv,iDir1,jDir1,kDir1,iDir2,jDir2,kDir2) .or. &
                          only_corner_fine(&
                          iNode_B(iBlockSend),iDir1,jDir1,kDir1,&
                          iDir3,jDir3,kDir3))then

                        if((iDir == iDir3 .and. &
                             (jDir1 == jDir3 .or. kDir1 == kDir3))&
                             .or. &
                             (iDir == -iDir2 .and. &
                             (jDir1 == jDir2 .or. kDir1 == kDir2))) then
                           CYCLE                                             
                        endif

                     endif


                     iEqualSOrig_DII = iEqualS_DII
                     iEqualROrig_DII = iEqualR_DII

                     if(kDir1 == -1) then
                        iEqualS_DII(3,0,Min_) = 1-nWidth
                        iEqualS_DII(3,0,Max_) = 0

                        iEqualR_DII(3,0,Min_) = 1-nWidth
                        iEqualR_DII(3,0,Max_) = 0
                     elseif(kDir1 == 1) then
                        iEqualS_DII(3,0,Min_) = nK + 1
                        iEqualS_DII(3,0,Max_) = nK + nWidth

                        iEqualR_DII(3,0,Min_) = nK + 1
                        iEqualR_DII(3,0,Max_) = nK + nWidth
                     endif

                     if(jDir1 == -1) then
                        iEqualS_DII(2,0,Min_) = 1-nWidth
                        iEqualS_DII(2,0,Max_) = 0

                        iEqualR_DII(2,0,Min_) = 1-nWidth
                        iEqualR_DII(2,0,Max_) = 0
                     elseif(jDir1 == 1) then
                        iEqualS_DII(2,0,Min_) = nJ + 1
                        iEqualS_DII(2,0,Max_) = nJ + nWidth

                        iEqualR_DII(2,0,Min_) = nJ + 1
                        iEqualR_DII(2,0,Max_) = nJ + nWidth
                     endif
                     call do_equal               
                     iEqualS_DII = iEqualSOrig_DII
                     iEqualR_DII = iEqualROrig_DII
                  endif
               enddo; enddo
            elseif(jDir /= 0) then
               jDir1 = 0
               do kDir1 = -1, 1; do iDir1 = -1, 1
                  if(kDir1 == 0 .and. iDir1 == 0) CYCLE
                  if(kDir1 /= 0 .and. iDir1 /= 0) CYCLE
                  if(DiLevelNei_IIIB(iDir1,jDir1,kDir1,iBlockSend) == 1 .or. &
                       DiLevelNei_IIIB(iDir1,jDir,kDir1,iBlockSend) == 1) then
                     
                     iSend = (3*iDir + 3)/2
                     jSend = (3*jDir + 3)/2
                     kSend = (3*kDir + 3)/2
                     iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)

                     iDir2 = -100; jDir2 = -100; kDir2 = -100
                     iDir3 = -100; jDir3 = -100; kDir3 = -100

                     if(only_corner_fine(&
                          iNodeRecv,iDir1,jDir1,kDir1,iDir2,jDir2,kDir2) .or. &
                          only_corner_fine(&
                          iNode_B(iBlockSend),iDir1,jDir1,kDir1,&
                          iDir3,jDir3,kDir3))then

                        if((jDir == jDir3 .and. &
                             (iDir1 == iDir3 .or. kDir1 == kDir3))&
                             .or. &
                             (jDir == -jDir2 .and. &
                             (iDir1 == iDir2 .or. kDir1 == kDir2))) then
                           CYCLE                                             
                        endif
                     endif
                     
                     iEqualSOrig_DII = iEqualS_DII
                     iEqualROrig_DII = iEqualR_DII

                     if(kDir1 == -1) then
                        iEqualS_DII(3,0,Min_) = 1-nWidth
                        iEqualS_DII(3,0,Max_) = 0

                        iEqualR_DII(3,0,Min_) = 1-nWidth
                        iEqualR_DII(3,0,Max_) = 0
                     elseif(kDir1 == 1) then
                        iEqualS_DII(3,0,Min_) = nK + 1
                        iEqualS_DII(3,0,Max_) = nK + nWidth

                        iEqualR_DII(3,0,Min_) = nK + 1
                        iEqualR_DII(3,0,Max_) = nK + nWidth
                     endif

                     if(iDir1 == -1) then
                        iEqualS_DII(1,0,Min_) = 1-nWidth
                        iEqualS_DII(1,0,Max_) = 0

                        iEqualR_DII(1,0,Min_) = 1-nWidth
                        iEqualR_DII(1,0,Max_) = 0
                     elseif(iDir1 == 1) then
                        iEqualS_DII(1,0,Min_) = nI + 1
                        iEqualS_DII(1,0,Max_) = nI + nWidth

                        iEqualR_DII(1,0,Min_) = nI + 1
                        iEqualR_DII(1,0,Max_) = nI + nWidth
                     endif

                     call do_equal               
                     iEqualS_DII = iEqualSOrig_DII
                     iEqualR_DII = iEqualROrig_DII
                  endif
               enddo; enddo

            elseif(kDir /= 0) then
               kDir1 = 0
               do jDir1 = -1, 1; do iDir1 = -1, 1
                  if(jDir1 == 0 .and. iDir1 == 0) CYCLE
                  if(jDir1 /= 0 .and. iDir1 /= 0) CYCLE
                  if(DiLevelNei_IIIB(iDir1,jDir1,kDir1,iBlockSend) == 1 .or. &
                       DiLevelNei_IIIB(iDir1,jDir1,kDir,iBlockSend) == 1) then
                     
                     iSend = (3*iDir + 3)/2
                     jSend = (3*jDir + 3)/2
                     kSend = (3*kDir + 3)/2
                     iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)

                     iDir2 = -100; jDir2 = -100; kDir2 = -100
                     iDir3 = -100; jDir3 = -100; kDir3 = -100

                     if(only_corner_fine(&
                          iNodeRecv,iDir1,jDir1,kDir1,iDir2,jDir2,kDir2) .or. &
                          only_corner_fine(&
                          iNode_B(iBlockSend),iDir1,jDir1,kDir1,&
                          iDir3,jDir3,kDir3))then

                        if((kDir == kDir3 .and. &
                             (iDir1 == iDir3 .or. jDir1 == jDir3))&
                             .or. &
                             (kDir == -kDir2 .and. &
                             (iDir1 == iDir2 .or. jDir1 == jDir2))) then
                           CYCLE                                             
                        endif
                     endif

                     iEqualSOrig_DII = iEqualS_DII
                     iEqualROrig_DII = iEqualR_DII

                     if(jDir1 == -1) then
                        iEqualS_DII(2,0,Min_) = 1-nWidth
                        iEqualS_DII(2,0,Max_) = 0

                        iEqualR_DII(2,0,Min_) = 1-nWidth
                        iEqualR_DII(2,0,Max_) = 0
                     elseif(jDir1 == 1) then
                        iEqualS_DII(2,0,Min_) = nJ + 1
                        iEqualS_DII(2,0,Max_) = nJ + nWidth

                        iEqualR_DII(2,0,Min_) = nJ + 1
                        iEqualR_DII(2,0,Max_) = nJ + nWidth
                     endif

                     if(iDir1 == -1) then
                        iEqualS_DII(1,0,Min_) = 1-nWidth
                        iEqualS_DII(1,0,Max_) = 0

                        iEqualR_DII(1,0,Min_) = 1-nWidth
                        iEqualR_DII(1,0,Max_) = 0
                     elseif(iDir1 == 1) then
                        iEqualS_DII(1,0,Min_) = nI + 1
                        iEqualS_DII(1,0,Max_) = nI + nWidth

                        iEqualR_DII(1,0,Min_) = nI + 1
                        iEqualR_DII(1,0,Max_) = nI + nWidth
                     endif

                     call do_equal               
                     iEqualS_DII = iEqualSOrig_DII
                     iEqualR_DII = iEqualROrig_DII
                  endif
               enddo; enddo
            endif
         endif
      endif ! nDim
    end subroutine corrected_do_equal

    !======================================================================
    subroutine calc_accurate_coarsened_block(iBlock)

      use BATL_high_order, ONLY: limit_interpolation, restriction_high_order_amr
      ! For nI*nJ*nK fine block, calculate its coarsened nI/2 * nJ/2 * nK/2
      ! overlaped block. 

      integer, intent(in):: iBlock
      integer:: iPerp
      integer:: iDir1, jDir1, kDir1, iVar
      integer:: iDir2, jDir2, kDir2
      integer:: iDirMin, iDirMax, jDirMin, jDirMax, kDirMin, kDirMax
      integer:: DiDir, DjDir, DkDir
      integer:: i, j, k, Di, Dj, Dk
      integer:: i0, j0, k0, ic0, jc0, kc0
      integer:: i1, i2, j1, j2, k1, k2
      integer:: iBegin, iEnd, jBegin, jEnd, kBegin, kEnd

      real, allocatable:: Fine_VIII(:,:,:,:)
      real:: CoarseCell, Coarse_I(3),Cell_I(5)
      real:: Cell1_I(6)=0, Cell2_I(6)=0, Cell3_I(6)=0
      real:: Distance_I(4) = (/-2,-1,1,2/)
      real:: Cell_III(6,6,6)
      logical, allocatable:: IsAccurate_III(:,:,:) 
      integer, allocatable:: nCorrected_III(:,:,:)
      logical:: DoResChange_D(3), IsAccurateGhost, DoSymInterp
      real, parameter:: c1=0.05, c2=-0.3, c3=0.75, c1over10=1./10, c1over2=1./2
      real:: Coef_I(6) 
      real:: Orig, Orig1, Orig2, Orig3, Res1, Res2, Res3
      integer:: nResChange, nEdge, nCorrect 
      character(len=*), parameter :: NameSub = 'calc_accurate_coarsened_block'
      !----------------------------------------------------------------------

      if(Do6thCorrect) then
         Coef_I = (/0.05, -0.3, 0.75, 0.75, -0.3, 0.05/)
      else
         Coef_I = (/0.1, -0.5, 1.0, 0.5, -0.1, 0.0/)
      endif

      if(.not. allocated(Fine_VIII))&
           allocate(Fine_VIII(nVar,8,6,min(6,nK)))

      if(.not. allocated(IsAccurate_III)) &
           allocate(IsAccurate_III(max(nI/2,1), max(nJ/2,1),max(nK/2,1)))
      IsAccurate_III = .false.

      DoResChange_D = .false.
      DoSymInterp = .true.

      ! Resolution change in x-dir
      jDir1 = 0; kDir1 = 0
      do iDir1 = -1, 1, 2
         if(DiLevelNei_IIIB(iDir1,jDir1,kDir1,iBlock) /=1) CYCLE
         ! Resolution change can happen in at most three directions.
         DoResChange_D(1) = .true.

         ! Are the ghost cells close to the opposite face accurate?
         ! For the more than 2 levels refinement case. 
         IsAccurateGhost = all(DiLevelNei_IIIB(-iDir1,:,:,iBlock) /= -1)
         DoSymInterp = nI .ge. 8 .or. IsAccurateGhost

         if(iDir1 == -1) then
            ! i0 is the index of the origin block. 
            ! ic0 is the index of the coarsened block. 
            i0 = 1; ic0 = 1
         elseif(iDir1 == 1) then
            i0 = nI; ic0 = nI/2
         endif

         jBegin = 1; jEnd = max(nJ/2,1)
         if(DiLevelNei_IIIB(0,1,0,iBlock) == 1) then
            ! For this situation, the calculation of cells j = nJ/2 
            ! involves values of coarsened neighbour block. These 
            ! cells need to be treated in a special way. 
            jBegin = 1; jEnd = max(nJ/2-1,1)
         elseif(DiLevelNei_IIIB(0,-1,0,iBlock) == 1) then
            jBegin = max(nJ/2,1); jEnd = min(2,nJ)
         endif
         Dj = sign(1,jEnd - jBegin)

         ! For 3D.
         kBegin = 1; kEnd = max(nK/2,1)
         if(DiLevelNei_IIIB(0,0,1,iBlock) == 1) then
            kBegin = 1; kEnd = max(nK/2-1,1)
         elseif(DiLevelNei_IIIB(0,0,-1,iBlock) == 1) then
            kBegin = max(nK/2,1); kEnd = min(2,nK)
         endif
         Dk = sign(1,kEnd - kBegin)

         do k = kBegin, kEnd, Dk; do j = jBegin, jEnd, Dj
            if(nK == 1) then ! 2D
               Fine_VIII(:,:,:,k) = &
                    State_VGB(:,&
                    i0:i0-iDir1*7:-iDir1,&
                    2*j-3:2*j+2,&
                    k,iBlock)
            else  ! 3D
               Fine_VIII = &
                    State_VGB(:,&
                    i0:i0-iDir1*7:-iDir1,&
                    2*j-3:2*j+2,&
                    2*k-3:2*k+2,iBlock)
            endif

            do iVar = 1, nVar
               CoarseCell = State_VGB(iVar,i0+iDir1,2*j,min(2*k,nK),iBlock)
               call restriction_high_order_reschange(CoarseCell, &
                    Fine_VIII(iVar,:,:,:), Coarse_I, DoSymInterp)
               State_VIIIB(iVar,ic0:ic0-2*iDir1:-iDir1,j,k,iBlock) = Coarse_I
            enddo ! iVar

            IsAccurate_III(ic0:ic0-2*iDir1:-iDir1,j,k) = .true.
         enddo; enddo
      enddo

      ! Resolution change in y-dir
      iDir1 = 0; kDir1 = 0
      do jDir1 = -1, 1, 2
         if(DiLevelNei_IIIB(iDir1,jDir1,kDir1,iBlock) /=1) CYCLE
         DoResChange_D(2) = .true.

         IsAccurateGhost = all(DiLevelNei_IIIB(:,-jDir1,:,iBlock) /= -1)
         DoSymInterp = nJ .ge. 8 .or. IsAccurateGhost

         if(jDir1 == -1) then
            j0 = 1; jc0 = 1
         elseif(jDir1 == 1) then
            j0 = nJ; jc0 = nJ/2
         endif

         iBegin = 1; iEnd = max(nI/2,1)
         if(DiLevelNei_IIIB(1,0,0,iBlock) == 1) then
            iBegin = 1; iEnd = max(nI/2-1,1)
         elseif(DiLevelNei_IIIB(-1,0,0,iBlock) == 1) then
            iBegin = max(nI/2,1); iEnd = min(2,nI)
         endif
         Di = sign(1,iEnd - iBegin)

         ! For 3D.
         kBegin = 1; kEnd = max(nK/2,1)
         if(DiLevelNei_IIIB(0,0,1,iBlock) == 1) then
            kBegin = 1; kEnd = max(nK/2-1,1)
         elseif(DiLevelNei_IIIB(0,0,-1,iBlock) == 1) then
            kBegin = max(nK/2,1); kEnd = min(2,nK)
         endif
         Dk = sign(1,kEnd - kBegin)

         do k = kBegin, kEnd, Dk; do i = iBegin, iEnd, Di
            ! The coarsened cell has been calculated. 
            if(IsAccurate_III(i,jc0,k)) CYCLE

            if(nK == 1) then
               do iPerp = 1, 8
                  Fine_VIII(:,iPerp,:,k) = &
                       State_VGB(:,&
                       2*i-3:2*i+2, &
                       j0-jDir1*(iPerp-1),&
                       k,iBlock)
               enddo
            else ! 3D
               do iPerp = 1, 8
                  Fine_VIII(:,iPerp,:,:) = &
                       State_VGB(:,&
                       2*i-3:2*i+2, &
                       j0-jDir1*(iPerp-1),&
                       2*k-3:2*k+2,iBlock)
               enddo
            endif

            do iVar = 1, nVar
               CoarseCell = State_VGB(iVar,2*i,j0+jDir1,min(2*k,nK),iBlock)

               call restriction_high_order_reschange(CoarseCell, &
                    Fine_VIII(iVar,:,:,:), Coarse_I, DoSymInterp)

               State_VIIIB(iVar,i,jc0:jc0-2*jDir1:-jDir1,k,iBlock) = Coarse_I

            enddo ! iVar
            IsAccurate_III(i,jc0:jc0-2*jDir1:-jDir1,k) = .true.
         enddo; enddo
      enddo

      ! Resolution change in z-dir.
      iDir1 = 0; jDir1 = 0
      do kDir1 = -1, 1, 2
         if(DiLevelNei_IIIB(iDir1,jDir1,kDir1,iBlock) /=1) CYCLE
         DoResChange_D(3) = .true.

         IsAccurateGhost = all(DiLevelNei_IIIB(:,:,-kDir1,iBlock) /= -1)
         DoSymInterp = nK .ge. 8 .or. IsAccurateGhost

         if(kDir1 == -1) then
            k0 = 1; kc0 = 1
         else ! kDir1 == 1
            k0 = nK; kc0 = nK/2
         endif

         iBegin = 1; iEnd = max(nI/2,1)
         if(DiLevelNei_IIIB(1,0,0,iBlock) == 1) then
            iBegin = 1; iEnd = max(nI/2-1,1)
         elseif(DiLevelNei_IIIB(-1,0,0,iBlock) == 1) then
            iBegin = max(nI/2,1); iEnd = min(2,nI)
         endif
         Di = sign(1,iEnd - iBegin)

         jBegin = 1; jEnd = max(nJ/2,1)
         if(DiLevelNei_IIIB(0,1,0,iBlock) == 1) then
            jBegin = 1; jEnd = max(nJ/2-1,1)
         elseif(DiLevelNei_IIIB(0,-1,0,iBlock) == 1) then
            jBegin = max(nJ/2,1); jEnd = min(2,nJ)
         endif
         Dj = sign(1,jEnd - jBegin)

         do j = jBegin, jEnd, Dj; do i = iBegin, iEnd, Di
            if(IsAccurate_III(i,j,kc0)) CYCLE

            do iPerp = 1, 8
               Fine_VIII(:,iPerp,:,:) = &
                    State_VGB(:,&
                    2*i-3:2*i+2,&
                    2*j-3:2*j+2,&
                    k0-kDir1*(iPerp-1),&
                    iBlock)
            enddo

            do iVar = 1, nVar
               CoarseCell = State_VGB(iVar,2*i,2*j,k0+kDir1,iBlock)
               call restriction_high_order_reschange(CoarseCell, &
                    Fine_VIII(iVar,:,:,:), Coarse_I, DoSymInterp)

               State_VIIIB(iVar,i,j,kc0:kc0-2*kDir1:-kDir1,iBlock) = Coarse_I
            enddo
            IsAccurate_III(i,j,kc0:kc0-2*kDir1:-kDir1) = .true.
         enddo; enddo
      enddo ! kDir1

      ! At least one neighbour block is coarse when this subroutine called. 
      ! If it is not a face block, it will be a edge/corner block. 
      if(nK == 1) then ! 2D. Combine 2D and 3D part??
         ! Resolution change in the edge direction.
         if(.not. DoResChange_D(1) .and. .not.DoResChange_D(2)) then

            ! ___________________________________________
            ! |                    |         |          |                     
            ! |                    |         |          |                     
            ! |                    |         |          |                     
            ! |        3           |_________|__________|
            ! |                    |         |          |                      
            ! |                    |   14    |          |                      
            ! |                    |         |          |        y             
            ! |____________________|_________|__________|        | 
            ! |          |         |*        |          |        |             
            ! |          |   13    |   11    |   12     |        -----> x     
            ! |          |         |         |          |                     
            ! |__________|_________|_________|__________|                     
            ! |          |         |         |          |                     
            ! |          |         |   9     |   10     |                     
            ! |          |         |         |          |                      
            ! |__________|_________|_________|__________|
            !  
            ! 9-14 are the top layer of their parent block. 

            ! Coarsened cell * is interpolated diagonally. 

            ! Edge ghost cells.
            kDir1 = 0
            do iDir1 = -1, 1, 2; do jDir1 = -1, 1, 2
               if(DiLevelNei_IIIB(iDir1,jDir1,kDir1,iBlock) /=1) CYCLE

               if(iDir1 == 1) then
                  ic0 = nI/2; i0 = nI; Di  = -1
               else ! iDir1 = -1
                  ic0 = 1; i0 = 1; Di = 1
               endif

               if(jDir1 == 1) then
                  jc0 = nJ/2; j0 = nJ; Dj = -1
               else ! jDir1 = -1
                  jc0 = 1; j0 = 1; Dj = 1
               endif

               do iVar = 1, nVar
                  k = 1
                  do j = jc0, jc0+2*Dj, Dj; do i = ic0, ic0+2*Di, Di
                     if(j == jc0 .and. i == ic0) CYCLE
                     Cell_III(:,:,k) = &
                          State_VGB(iVar,2*i-3:2*i+2,2*j-3:2*j+2,k,iBlock)
                     State_VIIIB(iVar,i,j,k,iBlock) = &
                          restriction_high_order_amr(Cell_III)
                  enddo; enddo

                  ! Interpolate in diagonal direction.
                  j = jc0; i = ic0; k = 1

                  Cell_I(1) = State_VIIIB(iVar,i+2*Di,j+2*Dj,k,iBlock)
                  Cell_I(2) = State_VIIIB(iVar,i+  Di,j+  Dj,k,iBlock)
                  Cell_I(3) = State_VGB  (iVar,i0-  Di,j0-  Dj,k,iBlock)
                  Cell_I(4) = State_VGB  (iVar,i0-2*Di,j0-2*Dj,k,iBlock)
                  Cell_I(5) = State_VGB  (iVar,i0-3*Di,j0-3*Dj,k,iBlock)

                  Orig = -c1over10*Cell_I(1) + c1over2*Cell_I(2) + Cell_I(3) &
                       -c1over2*Cell_I(4) + c1over10*Cell_I(5)
                  State_VIIIB(iVar,i,j,k,iBlock) = &
                       limit_interpolation(Orig,Cell_I(1:4),Distance_I)
               enddo ! iVar
            enddo; enddo
         elseif(DoResChange_D(1) .and. DoResChange_D(2)) then

            ! ___________________________________________
            ! |                    |                    |
            ! |                    |                    |
            ! |                    |                    |
            ! |        3           |                    |
            ! |                    |                    |
            ! |                    |   4                |
            ! |                    |                    |        y 
            ! |____________________|_________ __________|        | 
            ! |                    |*        |          |        |            
            ! |                    |   8     |   7      |        -----> x     
            ! |                    |         |          |                    
            ! |        1           |_________|__________|                    
            ! |                    |         |          |                   
            ! |                    |   5     |   6      |                    
            ! |                    |         |          |                  
            ! |____________________|_________|__________|
            !       
            ! 5-8 are the top layer of their parent block. 

            ! Coarsened cell * can be corrected in x or y direction. Use
            ! the average.

            if(DiLevelNei_IIIB(-1,0,0,iBlock) == 1) then
               ic0 = 1; i0 = 1; Di = 1
            elseif(DiLevelNei_IIIB(1,0,0,iBlock) == 1) then
               ic0 = nI/2; i0 = nI; Di = -1
            else
               call CON_stop(NameSub//': This case should not happen! - case1')
            endif

            if(DiLevelNei_IIIB(0,-1,0,iBlock) == 1) then
               jc0 = 1; j0 = 1; Dj = 1
            elseif(DiLevelNei_IIIB(0,1,0,iBlock) == 1) then
               jc0 = nJ/2; j0 = nJ; Dj = -1
            else
               call CON_stop(NameSub//': This case should not happen! - case2')
            endif

            k = 1
            do iVar = 1,nVar
               ! Use the neighbour coarsened cells to correct the corner cell.
               Cell1_I(1:3) = State_VGB(iVar,i0-3*Di:i0-Di:Di, j0, k,iBlock)
               Cell2_I(1:3) = State_VGB(iVar,i0,j0-3*Dj:j0-Dj:Dj,k,iBlock)
               if(Do6thCorrect) then
                  Cell1_I(4:6) = &
                       State_VIIIB(iVar,ic0+Di:ic0+3*Di:Di,jc0,k,iBlock)
                  Cell2_I(4:6) = &
                       State_VIIIB(iVar,ic0, jc0+Dj:jc0+3*Dj:Dj,k,iBlock)
               else
                  Cell1_I(4:5) = &
                       State_VIIIB(iVar,ic0+Di:ic0+2*Di:Di,jc0,k,iBlock)
                  Cell2_I(4:5) = &
                       State_VIIIB(iVar,ic0, jc0+Dj:jc0+2*Dj:Dj,k,iBlock)
               endif

               Orig1 = sum(Coef_I*Cell1_I)
               Orig2 = sum(Coef_I*Cell2_I)

               Res1 = limit_interpolation(Orig1, Cell1_I(2:5), Distance_I)
               Res2 = limit_interpolation(Orig2, Cell2_I(2:5), Distance_I)

               State_VIIIB(iVar,ic0,jc0,k,iBlock) = 0.5*(Res1 + Res2)
            enddo
         endif

      else ! 3D
         nResChange = 0
         do i = 1, 3
            if(DoResChange_D(i)) nResChange = nResChange + 1
         enddo

         if(nResChange > 1) then ! nResChange is 2 or 3. 

            ! Example: coarsen block 11, nResChange == 2
            ! ___________________________________________
            ! |                    |                    |                     
            ! |                    |                    |                     
            ! |                    |                    |                    
            ! |         7          |         8          |                    
            ! |                    |                    |                   
            ! |                    |                    |                  
            ! |                    |                    |                 
            ! |____________________|____________________|
            ! |                    |        |           |                  
            ! |                    |        |           |                 
            ! |                    |        |           |               
            ! |          5         |________|___________|
            ! |                    |        |           |              
            ! |                    |        |           |               
            ! |                    |        |           |          
            ! |____________________|________|___________|
            !               TOP LAYER

            ! ___________________________________________
            ! |                    |                    |
            ! |                    |                    |
            ! |                    |                    |
            ! |        3           |                    |
            ! |                    |                    |
            ! |                    |   4                |
            ! |                    |                    |        y 
            ! |____________________|_________ __________|        | 
            ! |                    |*        |          |        |         
            ! |                    |   11    |   12     |        -----> x    
            ! |                    |         |          |                 
            ! |        1           |_________|__________|               
            ! |                    |         |          |               
            ! |                    |   9     |   10     |               
            ! |                    |         |          |               
            ! |____________________|_________|__________|
            !                 BOTTOM LAYER
            ! 9-12 are the top layer of their parent block. 

            ! For block 11, resolution change happens in x and y direcitons. 
            ! Only cells ic == 1 .and. jc == nJ/2 .and. 1 =< kc =< nK/2 need 
            ! to be corrected with coarse cell values. Interpolations can be 
            ! done in x direction or y direction. Use the average of both. 


            ! Example: coarsen block 11, nResChange == 3
            ! ___________________________________________
            ! |                    |                    |               
            ! |                    |                    |                 
            ! |                    |                    |                 
            ! |         7          |         8          |                  
            ! |                    |                    |                    
            ! |                    |                    |                  
            ! |                    |                    |               
            ! |____________________|____________________|
            ! |                    |                    |                 
            ! |                    |                    |               
            ! |                    |                    |           
            ! |          5         |         6          |
            ! |                    |                    |           
            ! |                    |                    |           
            ! |                    |                    |       
            ! |____________________|____________________|
            !               TOP LAYER

            ! ___________________________________________
            ! |                    |                    |
            ! |                    |                    |
            ! |                    |                    |
            ! |        3           |                    |
            ! |                    |                    |
            ! |                    |   4                |
            ! |                    |                    |        y 
            ! |____________________|_________ __________|        | 
            ! |                    |         |          |        |       
            ! |                    |   11    |   12     |        -----> x  
            ! |                    |         |          |                  
            ! |        1           |_________|__________|                 
            ! |                    |         |          |               
            ! |                    |   9     |   10     |            
            ! |                    |         |          |         
            ! |____________________|_________|__________|
            !                 BOTTOM LAYER
            ! 9-12 are the top layer of their parent block. 

            ! For block 11, resolution change happens in x, y and z directions.
            ! Correction similar with the nResChange==2 case, but one cell : 
            ! ic == 1 .and. jc == nJ/2 .and. kc == nK/2 can be corrected in 
            ! three directions. Also use the average. 

            if(DiLevelNei_IIIB(-1,0,0,iBlock) == 1) then
               ! i0: index of origin block. 
               ! ic0 & iBegin & iEnd: index of coarsened block. 
               ic0 = 1; i0 = 1; Di = 1
               iBegin = nI/2; iEnd = 2
            elseif(DiLevelNei_IIIB(1,0,0,iBlock) == 1) then
               ic0 = nI/2; i0 = nI; Di = -1
               iBegin = 1; iEnd = nI/2 - 1
            endif

            if(DiLevelNei_IIIB(0,-1,0,iBlock) == 1) then
               jc0 = 1; j0 = 1; Dj = 1
               jBegin = nJ/2; jEnd = 2
            elseif(DiLevelNei_IIIB(0,1,0,iBlock) == 1) then
               jc0 = nJ/2; j0 = nJ; Dj = -1
               jBegin = 1; jEnd = nJ/2 - 1
            endif

            if(DiLevelNei_IIIB(0,0,-1,iBlock) == 1) then
               kc0 = 1; k0 = 1; Dk = 1
               kBegin = nK/2; kEnd = 2
            elseif(DiLevelNei_IIIB(0,0,1,iBlock) == 1) then
               kc0 = nK/2; k0 = nK; Dk = -1
               kBegin = 1; kEnd = nK/2 - 1
            endif

            if(.not. DoResChange_D(3) .or. nResChange == 3) then
               if(.not. DoResChange_D(3)) then
                  kBegin = 1; kEnd = nK/2
               endif
               do k = kBegin, kEnd, sign(1,kEnd-kBegin)
                  do iVar = 1,nVar
                     Cell1_I(1:3) = &
                          State_VGB(iVar,i0-3*Di:i0-Di:Di, j0, 2*k,iBlock)
                     Cell2_I(1:3) = &
                          State_VGB(iVar,i0,j0-3*Dj:j0-Dj:Dj,2*k,iBlock)

                     if(Do6thCorrect) then
                        Cell1_I(4:6) = &
                             State_VIIIB(iVar,ic0+Di:ic0+3*Di:Di,jc0,k,iBlock)
                        Cell2_I(4:6) = &
                             State_VIIIB(iVar,ic0, jc0+Dj:jc0+3*Dj:Dj,k,iBlock)
                     else
                        Cell1_I(4:5) = &
                             State_VIIIB(iVar,ic0+Di:ic0+2*Di:Di,jc0,k,iBlock)
                        Cell2_I(4:5) = &
                             State_VIIIB(iVar,ic0, jc0+Dj:jc0+2*Dj:Dj,k,iBlock)
                     endif

                     Orig1 = sum(Coef_I*Cell1_I)
                     Orig2 = sum(Coef_I*Cell2_I)

                     Res1 = limit_interpolation(Orig1, Cell1_I(2:5),Distance_I)
                     Res2 = limit_interpolation(Orig2, Cell2_I(2:5),Distance_I)

                     State_VIIIB(iVar,ic0,jc0,k,iBlock) = 0.5*(Res1 + Res2)
                  enddo
               enddo
            endif

            if(.not. DoResChange_D(2) .or. nResChange == 3) then
               if(.not. DoResChange_D(2)) then
                  jBegin = 1; jEnd = nJ/2
               endif
               do j = jBegin, jEnd, sign(1,jEnd - jBegin)
                  do iVar = 1, nVar
                     Cell1_I(1:3) = &
                          State_VGB(iVar,i0-3*Di:i0-Di:Di,2*j,k0,iBlock)
                     Cell2_I(1:3) = &
                          State_VGB(iVar,i0,2*j,k0-3*Dk:k0-Dk:Dk,iBlock)

                     if(Do6thCorrect) then
                        Cell1_I(4:6) = &
                             State_VIIIB(iVar,ic0+Di:ic0+3*Di:Di,j,kc0,iBlock)
                        Cell2_I(4:6) = &
                             State_VIIIB(iVar,ic0,j,kc0+Dk:kc0+3*Dk:Dk,iBlock)
                     else
                        Cell1_I(4:5) = &
                             State_VIIIB(iVar,ic0+Di:ic0+2*Di:Di,j,kc0,iBlock)
                        Cell2_I(4:5) = &
                             State_VIIIB(iVar,ic0,j,kc0+Dk:kc0+2*Dk:Dk,iBlock)
                     endif

                     Orig1 = sum(Coef_I*Cell1_I)
                     Orig2 = sum(Coef_I*Cell2_I)

                     Res1 = limit_interpolation(Orig1, Cell1_I(2:5),Distance_I)
                     Res2 = limit_interpolation(Orig2, Cell2_I(2:5),Distance_I)

                     State_VIIIB(iVar,ic0,j,kc0,iBlock) = 0.5*(Res1 + Res2)
                  enddo
               enddo
            endif

            if(.not. DoResChange_D(1) .or. nResChange == 3) then
               if(.not. DoResChange_D(1)) then
                  iBegin = 1; iEnd = nI/2
               endif
               do i = iBegin, iEnd, sign(1,iEnd - iBegin)
                  do iVar = 1, nVar
                     Cell1_I(1:3) = &
                          State_VGB(iVar,2*i,j0-3*Dj:j0-Dj:Dj,k0,iBlock)
                     Cell2_I(1:3) = &
                          State_VGB(iVar,2*i,j0,k0-3*Dk:k0-Dk:Dk,iBlock)


                     if(Do6thCorrect) then
                        Cell1_I(4:6) = &
                             State_VIIIB(iVar,i,jc0+Dj:jc0+3*Dj:Dj,kc0,iBlock)
                        Cell2_I(4:6) = &
                             State_VIIIB(iVar,i,jc0,kc0+Dk:kc0+3*Dk:Dk,iBlock)
                     else
                        Cell1_I(4:5) = &
                             State_VIIIB(iVar,i,jc0+Dj:jc0+2*Dj:Dj,kc0,iBlock)
                        Cell2_I(4:5) = &
                             State_VIIIB(iVar,i,jc0,kc0+Dk:kc0+2*Dk:Dk,iBlock)
                     endif


                     Orig1 = sum(Coef_I*Cell1_I)
                     Orig2 = sum(Coef_I*Cell2_I)

                     Res1 = limit_interpolation(Orig1, Cell1_I(2:5),Distance_I)
                     Res2 = limit_interpolation(Orig2, Cell2_I(2:5),Distance_I)

                     State_VIIIB(iVar,i,jc0,kc0,iBlock) = 0.5*(Res1 + Res2)
                  enddo
               enddo
            endif

            if(nResChange == 3) then
               ! One corner cell. 
               do iVar = 1, nVar
                  Cell1_I(1:3) = &
                       State_VGB(iVar,i0-3*Di:i0-Di:Di, j0, k0,iBlock)
                  Cell2_I(1:3) = &
                       State_VGB(iVar,i0,j0-3*Dj:j0-Dj:Dj,k0,iBlock)
                  Cell3_I(1:3) = &
                       State_VGB(iVar,i0,j0,k0-3*Dk:k0-Dk:Dk,iBlock)
                  if(Do6thCorrect) then
                     Cell1_I(4:6) = &
                          State_VIIIB(iVar,ic0+Di:ic0+3*Di:Di,jc0,kc0,iBlock)
                     Cell2_I(4:6) = &
                          State_VIIIB(iVar,ic0, jc0+Dj:jc0+3*Dj:Dj,kc0,iBlock)
                     Cell3_I(4:6) = &
                          State_VIIIB(iVar,ic0,jc0,kc0+Dk:kc0+3*Dk:Dk,iBlock)
                  else
                     Cell1_I(4:5) = &
                          State_VIIIB(iVar,ic0+Di:ic0+2*Di:Di,jc0,kc0,iBlock)
                     Cell2_I(4:5) = &
                          State_VIIIB(iVar,ic0, jc0+Dj:jc0+2*Dj:Dj,kc0,iBlock)
                     Cell3_I(4:5) = &
                          State_VIIIB(iVar,ic0,jc0,kc0+Dk:kc0+2*Dk:Dk,iBlock)
                  endif



                  Orig1 = sum(Coef_I*Cell1_I)
                  Orig2 = sum(Coef_I*Cell2_I)
                  Orig3 = sum(Coef_I*Cell3_I)

                  Res1 = limit_interpolation(Orig1, Cell1_I(2:5), Distance_I)
                  Res2 = limit_interpolation(Orig2, Cell2_I(2:5), Distance_I)
                  Res3 = limit_interpolation(Orig3, Cell3_I(2:5), Distance_I)

                  State_VIIIB(iVar,ic0,jc0,kc0,iBlock) = &
                       (Res1 + Res2 + Res3)/3.0
               enddo
            endif

         elseif(nResChange == 1) then
            ! Example: coarsen block 11
            ! ___________________________________________
            ! |                    |                    |                  
            ! |                    |                    |            
            ! |                    |                    |              
            ! |         7          |         8          |         
            ! |                    |                    |         
            ! |                    |                    |           
            ! |                    |                    |            
            ! |____________________|____________________|
            ! |                    |                    |            
            ! |                    |                    |             
            ! |                    |                    |           
            ! |          5         |         6          |           
            ! |                    |                    |           
            ! |                    |                    |              
            ! |                    |                    |          
            ! |____________________|____________________|
            !               TOP LAYER

            ! ___________________________________________
            ! |                    |         |          |          
            ! |                    |         |          |             
            ! |                    |         |          |           
            ! |        3           |_________|__________|
            ! |                    |         |          |         
            ! |                    |   14    |          |        
            ! |                    |         |          |        y     
            ! |____________________|_________|__________|        | 
            ! |          |         |*        |          |        |     
            ! |          |   13    |   11    |   12     |        -----> x 
            ! |          |         |         |          |                
            ! |__________|_________|_________|__________|          
            ! |          |         |         |          |           
            ! |          |         |   9     |   10     |          
            ! |          |         |         |          |     
            ! |__________|_________|_________|__________|
            !                 BOTTOM LAYER
            ! 9-14 are the top layer of their parent block. 


            ! Assume resolution change happens in z direction for block 11.
            ! If block 3 is also refined, only need to coarsen block 11 to 
            ! fill in the face ghost cells of block 6, which is trivial. 
            ! If block 3 is coarse, there will be four kinds cells need 
            ! to be coarsened for block 11: 
            ! ic, jc, kc are the index of coarsened block, not the original
            ! block 11. 
            ! Type 1: 1 =< ic =< nI/2; 1 =< jc =< nJ/2; nk/2 -2 =< kc =< nk/2
            !         except for ic == 1 .and. jc==nJ/2.
            ! Type 2: ic == 1 .and. jc==nJ/2 .and. nk/2 -2 =< kc =< nk/2
            ! Type 3: ic == 1 .and. jc==nJ/2 .and. 1 =< kc =< nk/2 - 3
            ! Type 4: 1 =< ic =< nI/2; 1 =< jc =< nJ/2; 1 =< kc =< nk/2 -3
            !         except those belong to Type 3. 

            ! Type 1 and Type 2 are face ghost cells of block 6. Some of them
            ! are also edge ghost cells for block 3.
            ! Type 3 and Type 4 are only used for edge ghost cells of block 3.
            ! If nK is 6, there are no type 3 and type 4 cells. 

            ! Calculation: 
            ! Type 1: the same as simple resolution change. 
            ! Type 4: simple restriction use 6*6*6 fine cells. 
            ! Type 2&3: interpolation diagionally in x-y plane.

            if(DoResChange_D(1)) then               
               iDirMin =  0; iDirMax = 0; DiDir = 1
               jDirMin = -1; jDirMax = 1; DjDir = 2
               kDirMin = -1; kDirMax = 1; DkDir = 2
               iBegin = 1; iEnd = nI/2; Di = 1
            elseif(DoResChange_D(2)) then
               iDirMin = -1; iDirMax = 1; DiDir = 2
               jDirMin =  0; jDirMax = 0; DjDir = 1
               kDirMin = -1; kDirMax = 1; DkDir = 2
               jBegin = 1; jEnd = nJ/2; Dj = 1
            elseif(DoResChange_D(3)) then
               iDirMin = -1; iDirMax = 1; DiDir = 2
               jDirMin = -1; jDirMax = 1; DjDir = 2
               kDirMin =  0; kDirMax = 0; DkDir = 1
               kBegin = 1; kEnd = nK/2; Dk = 1
            endif

            do kDir2 = kDirMin, kDirMax, DkDir
               do jDir2 = jDirMin, jDirMax, DjDir
                  do iDir2 = iDirMin, iDirMax, DiDir
                     if(DiLevelNei_IIIB(iDir2,jDir2,kDir2,iBlock) /=1) CYCLE

                     if(iDir2 == 1) then
                        iBegin = nI/2; iEnd = nI/2 - 2; Di = -1; i0 = nI
                     elseif(iDir2 == -1) then
                        iBegin = 1; iEnd = 3; Di = 1; i0 = 1
                     endif

                     if(jDir2 == 1) then
                        jBegin = nJ/2; jEnd = nJ/2 - 2; Dj = -1; j0 = nJ
                     elseif(jDir2 == -1) then
                        jBegin = 1; jEnd = 3; Dj = 1; j0 = 1
                     endif

                     if(kDir2 == 1) then
                        kBegin = nK/2; kEnd = nK/2 - 2; Dk = -1; k0 = nK
                     elseif(kDir2 == -1) then
                        kBegin = 1; kEnd = 3; Dk = 1; k0 = 1
                     endif

                     ! Calc Type 4 cells. 
                     do k = kBegin, kEnd, Dk
                        do j = jBegin, jEnd, Dj
                           do i = iBegin, iEnd, Di
                              if(IsAccurate_III(i,j,k)) CYCLE
                              if(DoResChange_D(1) .and. &
                                   j == jBegin .and. k == kBegin) CYCLE
                              if(DoResChange_D(2) .and. &
                                   i == iBegin .and. k == kBegin) CYCLE
                              if(DoResChange_D(3) .and. &
                                   i == iBegin .and. j == jBegin) CYCLE

                              do iVar = 1, nVar
                                 Cell_III = State_VGB(iVar,&
                                      2*i-3:2*i+2,&
                                      2*j-3:2*j+2,&
                                      2*k-3:2*k+2,&
                                      iBlock)
                                 State_VIIIB(iVar,i,j,k,iBlock) = &
                                      restriction_high_order_amr(Cell_III)
                              enddo ! iVar
                              IsAccurate_III(i,j,k) = .true.

                           enddo ! i
                        enddo ! j
                     enddo ! k                                         

                     ! Calc Type 2 and Type 3. 
                     if(DoResChange_D(1)) then
                        j = jBegin; k = kBegin                      
                        do i = iBegin, iEnd, Di
                           do iVar = 1, nVar
                              Cell_I(1) = State_VIIIB&
                                   (iVar,i,j+2*Dj,k+2*Dk,iBlock)
                              Cell_I(2) = State_VIIIB&
                                   (iVar,i,j+  Dj,k+ Dk,iBlock)
                              Cell_I(3) = State_VGB  &
                                   (iVar,2*i,j0-  Dj,k0-  Dk,iBlock)
                              Cell_I(4) = State_VGB  &
                                   (iVar,2*i,j0-2*Dj,k0-2*Dk,iBlock)
                              Cell_I(5) = State_VGB  &
                                   (iVar,2*i,j0-3*Dj,k0-3*Dk,iBlock)

                              Orig = -c1over10*Cell_I(1) + c1over2*Cell_I(2) +&
                                   Cell_I(3) - c1over2*Cell_I(4) + &
                                   c1over10*Cell_I(5)
                              State_VIIIB(iVar,i,j,k,iBlock) =&
                                   limit_interpolation&
                                   (Orig,Cell_I(1:4),Distance_I)
                           enddo ! iVar
                           IsAccurate_III(i,j,k) = .true.
                        enddo ! i 
                     elseif(DoResChange_D(2)) then
                        i = iBegin; k = kBegin                        
                        do j = jBegin, jEnd, Dj
                           do iVar = 1, nVar
                              Cell_I(1) = State_VIIIB&
                                   (iVar,i+2*Di,j,k+2*Dk,iBlock)
                              Cell_I(2) = State_VIIIB&
                                   (iVar,i+  Di,j,k+  Dk,iBlock)
                              Cell_I(3) = State_VGB  &
                                   (iVar,i0-  Di,2*j,k0-  Dk,iBlock)
                              Cell_I(4) = State_VGB  &
                                   (iVar,i0-2*Di,2*j,k0-2*Dk,iBlock)
                              Cell_I(5) = State_VGB  &
                                   (iVar,i0-3*Di,2*j,k0-3*Dk,iBlock)

                              Orig = -c1over10*Cell_I(1) + c1over2*Cell_I(2) +&
                                   Cell_I(3) - c1over2*Cell_I(4) + &
                                   c1over10*Cell_I(5)
                              State_VIIIB(iVar,i,j,k,iBlock) =&
                                   limit_interpolation&
                                   (Orig,Cell_I(1:4),Distance_I)
                           enddo ! iVar
                           IsAccurate_III(i,j,k) = .true.
                        enddo ! j 

                     elseif(DoResChange_D(3)) then
                        i = iBegin; j = jBegin
                        do k = kBegin, kEnd, Dk
                           do iVar = 1, nVar
                              Cell_I(1) = State_VIIIB&
                                   (iVar,i+2*Di,j+2*Dj,k,iBlock)
                              Cell_I(2) = State_VIIIB&
                                   (iVar,i+  Di,j+  Dj,k,iBlock)
                              Cell_I(3) = State_VGB  &
                                   (iVar,i0-  Di,j0-  Dj,2*k,iBlock)
                              Cell_I(4) = State_VGB  &
                                   (iVar,i0-2*Di,j0-2*Dj,2*k,iBlock)
                              Cell_I(5) = State_VGB  &
                                   (iVar,i0-3*Di,j0-3*Dj,2*k,iBlock)

                              Orig = -c1over10*Cell_I(1) + c1over2*Cell_I(2) +&
                                   Cell_I(3) - c1over2*Cell_I(4) + &
                                   c1over10*Cell_I(5)
                              State_VIIIB(iVar,i,j,k,iBlock) =&
                                   limit_interpolation&
                                   (Orig,Cell_I(1:4),Distance_I)
                           enddo ! iVar
                           IsAccurate_III(i,j,k) = .true.
                        enddo ! k
                     endif
                  enddo
               enddo
            enddo
         elseif(nResChange == 0) then
            ! Example: coarsen block 11
            ! ___________________________________________
            ! |                    |                    |             
            ! |                    |                    |             
            ! |                    |                    |                
            ! |         7          |         8          |           
            ! |                    |                    |            
            ! |                    |                    |            
            ! |                    |                    |       
            ! |____________________|____________________|
            ! |                    |        |           |           
            ! |                    |        |           |             
            ! |                    |        |           |           
            ! |          5         _________|___________|
            ! |                    |        |           |          
            ! |                    |        |           |          
            ! |                    |        |           |          
            ! |____________________|________|___________|
            !               TOP LAYER

            ! ___________________________________________
            ! |                    |         |          |        
            ! |                    |         |          |           
            ! |                    |         |          |              
            ! |        3           |_________|__________|
            ! |                    |         |          |         
            ! |                    |   14    |          |       
            ! |                    |         |          |        y          
            ! |____________________|_________|__________|        | 
            ! |          |         |*        |          |        |      
            ! |          |   13    |   11    |   12     |        -----> x  
            ! |          |         |         |          |                
            ! |__________|_________|_________|__________|              
            ! |          |         |         |          |          
            ! |          |         |   9     |   10     |   
            ! |          |         |         |          |  
            ! |__________|_________|_________|__________|
            !                 BOTTOM LAYER
            ! 9-14 are the top layer of their parent block. 

            ! Six faces neighbours of block 11 are fine blocks. 
            ! But at most three edg blocks (block 3, 5, 8) and one 
            ! corner (block 7) may not. 

            ! Case 1: only one edge block (like block 3) is coarse.
            !         ic = 1 .and. jc = nJ/2 .and. 1 =< kc =< nK/2
            !         are interpolated diagonally in x-y plane. Other
            !         needed values are restricted with 6*6*6 fine cells.
            ! Case 2: n ( 1<n<=3 ) edge blocks are coarse. Similar 
            !         with Case 1, but one cell: ic=1,jc=nk/2,kc=nk/2
            !         can be interpolated in n directions. Use the average
            !         of these n interpolations. 
            ! Case 3: only the corner block is coarse. Need to coarsen 11 
            !         to fill in the corner ghost of block 7. 
            !         All needed cells except for ic=1,jc=nk/2,kc=nk/2 can 
            !         be rectricted with 6*6*6 fine cells. 
            !         The corner one (ic=1,jc=nk/2,kc=nk/2) is interpolated 
            !         diagionally in 3D space. 

            State_VIIIB = 0
            nEdge = 0
            if(.not. allocated(nCorrected_III)) allocate(&
                 nCorrected_III(max(nI/2,1), max(nJ/2,1),max(nK/2,1)))       
            nCorrected_III = 0

            do kDir1 = -1, 1; do jDir1 = -1, 1; do iDir1 = -1, 1
               ! Concentrate on the edge neighbour blocks.
               if(abs(iDir1) + abs(jDir1) + abs(kDir1) /= 2) CYCLE
               if(DiLevelNei_IIIB(iDir1,jDir1,kDir1,iBlock) /=1) CYCLE

               ! Record how many edge neighbour blocks are coarse.
               nEdge = nEdge + 1

               if(iDir1 == 1) then
                  iBegin = nI/2; iEnd = nI/2 - 2; Di = -1; i0 = nI
               elseif(iDir1 == -1) then
                  iBegin = 1; iEnd = 3; Di = 1; i0 = 1
               elseif(iDir1 == 0) then
                  iBegin = 1; iEnd = nI/2; Di = 1
               endif

               if(jDir1 == 1) then
                  jBegin = nJ/2; jEnd = nJ/2 - 2; Dj = -1; j0 = nJ
               elseif(jDir1 == -1) then
                  jBegin = 1; jEnd = 3; Dj = 1; j0 = 1
               elseif(jDir1 == 0) then
                  jBegin = 1; jEnd = nJ/2; Dj = 1
               endif


               if(kDir1 == 1) then
                  kBegin = nK/2; kEnd = nK/2 - 2; Dk = -1; k0 = nK
               elseif(kDir1 == -1) then
                  kBegin = 1; kEnd = 3; Dk = 1; k0 = 1
               elseif(kDir1 == 0) then
                  kBegin = 1; kEnd = nK/2; Dk = 1
               endif

               ! Simple restrictioin for 'inner' cells. 
               ! Simple restriction with 6*6*6 fine cells. But some cells
               ! calculated in this way may not accurate and will be 
               ! overwritten in the later part. 
               do k = kBegin, kEnd, Dk
                  do j = jBegin, jEnd, Dj
                     do i = iBegin, iEnd, Di
                        if(IsAccurate_III(i,j,k)) CYCLE
                        if(kDir1 == 0 .and. &
                             (i==iBegin .and. j==jBegin)) CYCLE
                        if(jDir1 == 0 .and. &
                             i==iBegin .and. k==kBegin) CYCLE
                        if(iDir1 == 0 .and. &
                             j==jBegin .and. k==kBegin) CYCLE
                        do iVar = 1, nVar
                           Cell_III = State_VGB(iVar,&
                                2*i-3:2*i+2,2*j-3:2*j+2,2*k-3:2*k+2,iBlock)
                           ! Some value calculated here is not accurate, 
                           ! which will be corrected in the next part. 
                           State_VIIIB(iVar,i,j,k,iBlock) = &
                                restriction_high_order_amr(Cell_III)
                        enddo

                        ! It is somewhat complicated to tell weather it is
                        ! accurate or not. So, do not set IsAccurate_III value.

                     enddo ! i 
                  enddo ! j 
               enddo ! k 

               do k = kBegin, kEnd, Dk
                  do j = jBegin, jEnd, Dj
                     do i = iBegin, iEnd, Di
                        if(kDir1 == 0 .and. &
                             .not. (i==iBegin .and. j==jBegin)) CYCLE
                        if(jDir1 == 0 .and. &
                             .not. (i==iBegin .and. k==kBegin)) CYCLE
                        if(iDir1 == 0 .and. &
                             .not. (j==jBegin .and. k==kBegin)) CYCLE

                        do iVar = 1, nVar
                           if(kDir1 == 0) then
                              Cell_I(1) = &
                                   State_VIIIB(iVar,i+2*Di,j+2*Dj,k,iBlock)
                              Cell_I(2) = &
                                   State_VIIIB(iVar,i+  Di,j+  Dj,k,iBlock)
                              Cell_I(3) = &
                                   State_VGB(iVar,i0-  Di,j0-  Dj,2*k,iBlock)
                              Cell_I(4) = &
                                   State_VGB(iVar,i0-2*Di,j0-2*Dj,2*k,iBlock)
                              Cell_I(5) = &
                                   State_VGB(iVar,i0-3*Di,j0-3*Dj,2*k,iBlock)
                           endif

                           if(jDir1 == 0) then
                              Cell_I(1) = &
                                   State_VIIIB(iVar,i+2*Di,j,k+2*Dk,iBlock)
                              Cell_I(2) = &
                                   State_VIIIB(iVar,i+  Di,j,k+  Dk,iBlock)
                              Cell_I(3) = &
                                   State_VGB(iVar,i0-  Di,2*j,k0-  Dk,iBlock)
                              Cell_I(4) = &
                                   State_VGB(iVar,i0-2*Di,2*j,k0-2*Dk,iBlock)
                              Cell_I(5) = &
                                   State_VGB(iVar,i0-3*Di,2*j,k0-3*Dk,iBlock)
                           endif

                           if(iDir1 == 0) then
                              Cell_I(1) = &
                                   State_VIIIB(iVar,i,j+2*Dj,k+2*Dk,iBlock)
                              Cell_I(2) = &
                                   State_VIIIB(iVar,i,j+  Dj,k+  Dk,iBlock)
                              Cell_I(3) = &
                                   State_VGB(iVar,2*i,j0-  Dj,k0-  Dk,iBlock)
                              Cell_I(4) = &
                                   State_VGB(iVar,2*i,j0-2*Dj,k0-2*Dk,iBlock)
                              Cell_I(5) = &
                                   State_VGB(iVar,2*i,j0-3*Dj,k0-3*Dk,iBlock)

                           endif
                           Orig = -c1over10*Cell_I(1) + c1over2*Cell_I(2) + &
                                Cell_I(3) - c1over2*Cell_I(4) + &
                                c1over10*Cell_I(5)
                           nCorrect = nCorrected_III(i,j,k)
                           if(nCorrect == 0) then
                              ! 2D diagonally interpolation. 
                              State_VIIIB(iVar,i,j,k,iBlock) = &
                                   limit_interpolation&
                                   (Orig,Cell_I(1:4),Distance_I)
                           else
                              ! Corner cell for Case 2. 
                              ! Some cells can be corrected in different 
                              ! direction. Use the average of these correcitons.
                              State_VIIIB(iVar,i,j,k,iBlock) = &
                                   (nCorrect*State_VIIIB(iVar,i,j,k,iBlock) +&
                                   limit_interpolation&
                                   (Orig,Cell_I(1:4),Distance_I)) / &
                                   (nCorrect + 1)
                           endif
                        enddo ! iVar

                        ! If it is accurate, it will not be overwritten by 
                        ! simple restriction for the inner cell code. 
                        IsAccurate_III(i,j,k) = .true.

                        nCorrected_III(i,j,k) = nCorrected_III(i,j,k) + 1
                     enddo ! i 
                  enddo ! j 
               enddo ! k 
            enddo; enddo; enddo

            if(nEdge == 0) then
               ! Case 3. 
               ! Take care corner neighbour block.
               do kDir1 = -1, 1, 2; do jDir1 = -1, 1, 2; do iDir1 = -1, 1,2
                  if(DiLevelNei_IIIB(iDir1,jDir1,kDir1,iBlock) /=1) CYCLE

                  if(iDir1 == 1) then
                     iBegin = nI/2; iEnd = nI/2 - 2; Di = -1; i0 = nI
                  elseif(iDir1 == -1) then
                     iBegin = 1; iEnd = 3; Di = 1; i0 = 1
                  endif

                  if(jDir1 == 1) then
                     jBegin = nJ/2; jEnd = nJ/2 - 2; Dj = -1; j0 = nJ
                  elseif(jDir1 == -1) then
                     jBegin = 1; jEnd = 3; Dj = 1; j0 = 1
                  endif

                  if(kDir1 == 1) then
                     kBegin = nK/2; kEnd = nK/2 - 2; Dk = -1; k0 = nK
                  elseif(kDir1 == -1) then
                     kBegin = 1; kEnd = 3; Dk = 1; k0 = 1
                  endif

                  ! Simple 6*6*6 restriction. 
                  do k = kBegin, kEnd, Dk
                     do j = jBegin, jEnd, Dj
                        do i = iBegin, iEnd, Di
                           if(i == iBegin .and. j == jBegin .and. k == kBegin)&
                                CYCLE
                           do iVar = 1, nVar
                              Cell_III = State_VGB(iVar,&
                                   2*i-3:2*i+2,2*j-3:2*j+2,2*k-3:2*k+2,iBlock)
                              State_VIIIB(iVar,i,j,k,iBlock) = &
                                   restriction_high_order_amr(Cell_III)
                           enddo
                        enddo
                     enddo
                  enddo

                  do iVar = 1, nVar
                     ! 3D diagonal interpolation. 
                     i = iBegin; j = jBegin; k = kBegin
                     Cell_I(1) = &
                          State_VIIIB(iVar,i+2*Di,j+2*Dj,k+2*Dk,iBlock)
                     Cell_I(2) = &
                          State_VIIIB(iVar,i+  Di,j+  Dj,k+  Dk,iBlock)
                     Cell_I(3) = &
                          State_VGB(iVar,i0-  Di,j0-  Dj,k0-  Dk,iBlock)
                     Cell_I(4) = &
                          State_VGB(iVar,i0-2*Di,j0-2*Dj,k0-2*Dk,iBlock)
                     Cell_I(5) = &
                          State_VGB(iVar,i0-3*Di,j0-3*Dj,k0-3*Dk,iBlock)

                     Orig = -c1over10*Cell_I(1) + c1over2*Cell_I(2) + &
                          Cell_I(3) - c1over2*Cell_I(4) + c1over10*Cell_I(5)
                     State_VIIIB(iVar,i,j,k,iBlock) = &
                          limit_interpolation(Orig,Cell_I(1:4),Distance_I)

                  enddo ! iVar

               enddo; enddo; enddo

            endif ! nEdge

         endif ! nResChange
      endif ! nK 

      IsAccurate_B(iBlock) = .true.

    end subroutine calc_accurate_coarsened_block

    !======================================================================

    subroutine high_prolong_for_face_ghost(iBlock) 
      ! High order prolongation for face ghost cells. It is done locally. 
      use BATL_high_order, ONLY: prolongation_high_order_for_face_ghost, &
           correct_face_ghost_for_fine_block

      integer, intent(in):: iBlock
      real, allocatable:: Field1_VG(:,:,:,:)
      integer:: neiLev_I(6)
      !----------------------------------------------------------------------

      if(.not. allocated(Field1_VG)) &
           allocate(Field1_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK))            


      call is_face_accurate(iBlock)

      call prolongation_high_order_for_face_ghost(&
           iBlock, nVar, Field1_VG, State_VGB(:,:,:,:,iBlock), &
           IsAccurateFace_GB(:,:,:,iBlock))

      neiLev_I(1) = DiLevelNei_IIIB(-1,0,0,iBlock)
      neiLev_I(2) = DiLevelNei_IIIB(+1,0,0,iBlock)
      neiLev_I(3) = DiLevelNei_IIIB(0,-1,0,iBlock)
      neiLev_I(4) = DiLevelNei_IIIB(0,+1,0,iBlock)
      neiLev_I(5) = DiLevelNei_IIIB(0,0,-1,iBlock)
      neiLev_I(6) = DiLevelNei_IIIB(0,0,+1,iBlock)

      ! If the corner/edge block is not a coarse block, the ghost values for 
      ! fine block need to be corrected. 
      if(.not. all(neiLev_I /=1)) call correct_face_ghost_for_fine_block(&
           iBlock, nVar, State_VGB(:,:,:,:,iBlock))

    end subroutine high_prolong_for_face_ghost

    !==========================================================================

    subroutine is_face_accurate(iBlock)
      integer, intent(in):: iBlock
      logical:: IsOnlyCornerFine
      integer:: iDirCorner, jDirCorner, kDirCorner
      integer:: iBegin, iEnd, jBegin, jEnd, kBegin, kEnd, Di, Dj, Dk
      character(len=*), parameter :: NameSub = 'is_face_accurate'
      integer:: i,j,k

      logical:: DoTestMe = .false. 
      !----------------------------------------------------------------------
      ! Non-face ghost cells are also set false.
      IsAccurateFace_GB = .false.

      ! Assume face ghost cells are accurate. 
      IsAccurateFace_GB(-2:0,      1:nJ,1:nK,iBlock) = .true. 
      IsAccurateFace_GB(nI+1:nI+3, 1:nJ,1:nK,iBlock) = .true.
      IsAccurateFace_GB(1:nI,   -2:0,   1:nK,iBlock) = .true.
      IsAccurateFace_GB(1:nI,nJ+1:nJ+3, 1:nK,iBlock) = .true.
      if(nK == 1) RETURN
      IsAccurateFace_GB(1:nI,1:nJ,   -2:0,   iBlock) = .true.
      IsAccurateFace_GB(1:nI,1:nJ, nK+1:nK+3,iBlock) = .true.

      do iDir = -1, 1; do jDir = -1, 1; do kDir = -1, 1
         if(abs(iDir)+abs(jDir)+abs(kDir) /= 1) CYCLE
         IsOnlyCornerFine = only_corner_fine(iNode_B(iBlock),iDir,jDir,kDir,&
              iDirCorner,jDirCorner,kDirCorner)

         if(.not. IsOnlyCornerFine) CYCLE

         if(iDirCorner == 1) then 
            iBegin = nI; iEnd = nI -3; Di = -1
         elseif(iDirCorner==-1) then 
            iBegin = 1; iEnd = 4; Di = 1
         else 
            call CON_stop(NameSub//': This case should not happen! - case1')
         endif

         if(jDirCorner == 1) then 
            jBegin = nJ; jEnd = nJ -3; Dj = -1
         elseif(jDirCorner == -1) then 
            jBegin = 1; jEnd = 4; Dj = 1
         else
            call CON_stop(NameSub//': This case should not happen! - case2')
         endif

         if(kDirCorner == 1) then 
            kBegin = nK; kEnd = nK - 3; Dk = -1
         elseif(kDirCorner == -1) then 
            kBegin = 1; kEnd = 4; Dk = 1
         else
            call CON_stop(NameSub//': This case should not happen! - case3')
         endif

         ! This kind of things should be replaced by a matrix finally.  
         if(iDir == 1) then
            iBegin = nI+1; iEnd = nI+3; Di = 1
         elseif(iDir == -1) then
            iBegin = 0; iEnd = -2; Di = -1
         elseif(jDir == 1) then
            jBegin = nJ+1; jEnd = nJ+3; Dj = 1
         elseif(jDir == -1) then
            jBegin = 0; jEnd = -2; Dj = -1
         elseif(kDir == 1) then
            kBegin = nK+1; kEnd = nK+3; Dk = 1
         elseif(kDir == -1) then
            kBegin = 0; kEnd = -2; Dk = -1
         endif

         ! Find out the not accurate face cells. 
         IsAccurateFace_GB(iBegin:iEnd:Di,jBegin:jEnd:Dj,kBegin:kEnd:Dk,iBlock)&
              = .false.
      enddo; enddo; enddo ! iDir jDir kDir

    end subroutine is_face_accurate
    !======================================================================

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
                  if(.not. (nK >1 .and. iSendStage == 3 &
                       .and. IsAccurateFace_GB(i,j,k,iBlockRecv)))then    
                     State_VGB(:,i,j,k,iBlockRecv) = &
                          BufferR_I(iBufferR+1:iBufferR+nVar)
                  endif
                  iBufferR = iBufferR + nVar
               end do; end do; end do
            end if
            if(iBufferR >= sum(nBufferR_P(0:iProcSend))) EXIT
         end do
      end do

    end subroutine buffer_to_state

    !==========================================================================

    subroutine do_equal

      integer :: iBufferS, i, j, k, nSize, nWithin
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

      if(iSendStage == 3) then
         ! Only edge/corner cells need to be overwritten. 
         nWithin = 0
         if(.not.(iRMin .ge. 0 .and. iRMin .le. nI)) nWithin = nWithin + 1
         if(.not.(jRMin .ge. 0 .and. jRMin .le. nJ)) nWithin = nWithin + 1
         if(.not.(kRMin .ge. 0 .and. kRMin .le. nK)) nWithin = nWithin + 1
         if(nWithin < 1) RETURN
      endif

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

      if(iSendStage == 3 .and. nK > 1 .and. &
           abs(iDir)+abs(jDir)+abs(kDir) == 1) then
         DoRecvFace = only_corner_fine(iNode_B(iBlockSend),iDir,jDir,kDir)
         if(.not.DoRecvFace) RETURN
      endif

      ! For part implicit and part steady schemes
      if(Unused_BP(iBlockRecv,iProcRecv)) RETURN

      ! No need to count data for local copy
      if(DoCountOnly .and. iProc == iProcRecv) RETURN

      if(DoCountOnly .and. (&
           (.not. UseHighResChange .and. iSendStage == nProlongOrder) .or. &
           (UseHighResChange .and. (iSendStage == 1 .or. iSendStage == 3))))then 
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

      ! Do prolongation for edge/corner ghost cells remotely.
      if(UseHighResChange .and. iSendStage == 3) RETURN

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
            if(UseHighResChange) then
               if(.not.IsAccurate_B(iBlockSend)) &
                    call calc_accurate_coarsened_block(iBlockSend)
            endif

            if(UseHighResChange) then
               do kR = kRMin, kRMax, DkR
                  kS1 = kSMin + kRatioRestr*abs(kR-kRMin)
                  do jR = jRMin, jRMax, DjR
                     jS1 = jSMin + jRatioRestr*abs(jR-jRMin)
                     do iR = iRMin, iRMax, DiR
                        iS1 = iSMin + iRatioRestr*abs(iR-iRMin)
                        do iVar = 1, nVar
                           State_VGB(iVar,iR,jR,kR,iBlockRecv) = &
                                State_VIIIB(iVar,(iS1+1)/2,(jS1+1)/2,&
                                (kS1+1)/2,iBlockSend)
                        end do
                     enddo
                  enddo
               enddo

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
         if(UseHighResChange) then
            if(.not.IsAccurate_B(iBlockSend)) &
                 call calc_accurate_coarsened_block(iBlockSend)
            if(.not. allocated(State_VG)) then
               allocate(State_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK))
               State_VG = 1
            endif

            do kR = kRMin, kRMax, DkR
               kS1 = kSMin + kRatioRestr*abs(kR-kRMin)
               do jR = jRMin, jRMax, DjR
                  jS1 = jSMin + jRatioRestr*abs(jR-jRMin)
                  do iR = iRMin, iRMax, DiR
                     iS1 = iSMin + iRatioRestr*abs(iR-iRMin)                        
                     do iVar = 1, nVar
                        State_VG(iVar,iR,jR,kR) = &
                             State_VIIIB(iVar,(iS1+1)/2,(jS1+1)/2,(kS1+1)/2,&
                             iBlockSend)
                     end do
                  enddo
               enddo
            enddo
         endif

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
                  if(UseHighResChange) then
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
      use BATL_high_order, ONLY: prolongation_high_order_amr

      integer :: iR, jR, kR, iS, jS, kS, iS1, jS1, kS1
      integer :: iRatioRestr, jRatioRestr, kRatioRestr
      integer :: iBufferS, nSize
      integer, parameter:: Di=iRatio-1, Dj=jRatio-1, Dk=kRatio-1
      real    :: WeightOld, WeightNew, Weight, WeightI, WeightJ, WeightK, InvV
      real, dimension(MaxDim):: Xyz_D, dI_D, dJ_D, dK_D, dR_D, &
           PositionMinR_D, PositionMaxR_D, CoordMinR_D, CoordMaxR_D, &
           CellSizeR_D, CoordR_D

      logical :: UseSimpleWeights

      integer :: iVar
      integer:: nWidthProlongS_D(MaxDim), iDim
      real:: CoarseCell_III(5,5,5)
      integer:: i5_,j5_,k5_, iDir1, jDir1, kDir1
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

               if(iSendStage == 3 .and. nK > 1 .and. &
                    abs(iDir)+abs(jDir)+abs(kDir) .eq. 1 ) then 
                  ! Do_prolongation for edge/corner ghost cells and for 
                  ! some special face cells. 
                  DoSendFace = only_corner_fine(iNodeRecv,-iDir,-jDir,-kDir)
                  if(.not. DoSendFace) CYCLE
               endif

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

               ! For HighResChange, only do restriction in stage 2.
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

               if(UseHighResChange .and. iSendStage == 3) then
                  ! The values set in set_range are used for iSendStage == 1, 
                  ! Which is first order prolongtion. Now, for high order 
                  ! prolongation, some values need to be corrected.
                  nWidthProlongS_D(1:nDim) = 1 + (nWidth-1)/iRatio_D(1:nDim)
                  do iDim = 1, MaxDim
                     ! This loop is used to avoid the NAG 5.1 (282) bug on nyx
                     iProlongS_DII(iDim,0,Min_) = 1
                     iProlongS_DII(iDim,0,Max_) = nWidthProlongS_D(iDim)
                     iProlongS_DII(iDim,1,Min_) = 1
                     iProlongS_DII(iDim,1,Max_) = nIjk_D(iDim)/iRatio_D(iDim)
                     iProlongS_DII(iDim,2,Min_) = nIjk_D(iDim)/iRatio_D(iDim) + 1
                     iProlongS_DII(iDim,2,Max_) = nIjk_D(iDim)
                     iProlongS_DII(iDim,3,Min_) = &
                          nIjk_D(iDim) + 1 - nWidthProlongS_D(iDim)
                     iProlongS_DII(iDim,3,Max_) = nIjk_D(iDim)
                  end do

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
                     end do
                  end if
               endif

               ! Sending range depends on iSend,jSend,kSend = 0..3
               iSMin = iProlongS_DII(1,iSend,Min_)
               iSMax = iProlongS_DII(1,iSend,Max_)
               jSMin = iProlongS_DII(2,jSend,Min_)
               jSMax = iProlongS_DII(2,jSend,Max_)
               kSMin = iProlongS_DII(3,kSend,Min_)
               kSMax = iProlongS_DII(3,kSend,Max_)

               iRatioRestr = iRatio; jRatioRestr = jRatio; kRatioRestr = kRatio      
               if(iSendStage /= 3 .and. nCoarseLayer > 1)then
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
                     if(UseHighResChange .and. iSendStage == 3) then
                        iDir1 = 0; jDir1 = 0; kDir1 = 0
                        i5_ = max(5*Di,1); j5_ = max(5*Dj,1); k5_ =  max(5*Dk,1)
                        do kR = kRMin, kRMax, DkR
                           kS = kSMin + abs((kR+9)/kRatioRestr &
                                -           (kRMin+9)/kRatioRestr)
                           ! kDir = -1 if kR is even; kDir = 1 if kR is odd. 
                           ! kR may be negative (1-nK), kR+2*nK will be positive.
                           if(kRatioRestr == 2) kDir1 = 2*mod(kR+2*nK,2) - 1
                           do jR = jRMin, jRMax, DjR
                              jS = jSMin + abs((jR+9)/jRatioRestr &
                                   -           (jRMin+9)/jRatioRestr)
                              if(jRatioRestr == 2) jDir1 = 2*mod(jR+2*nJ,2) - 1
                              do iR = iRMin, iRMax, DiR
                                 iS = iSMin + abs((iR+9)/iRatioRestr &
                                      -           (iRMin+9)/iRatioRestr)
                                 if(iRatioRestr == 2) iDir1 = 2*mod(iR+2*nI,2) - 1

                                 if(IsAccurateFace_GB(iR,jR,kR,iBlockRecv))&
                                      CYCLE

                                 do iVar = 1, nVar
                                    CoarseCell_III(1:i5_,1:j5_,1:k5_) = &
                                         State_VGB(iVar,&
                                         iS-2*iDir1:iS+2*iDir1:sign(1,iDir1),&
                                         jS-2*jDir1:jS+2*jDir1:sign(1,jDir1),&
                                         kS-2*kDir1:kS+2*kDir1:sign(1,kDir1),&
                                         iBlockSend)

                                    State_VGB(iVar,iR,jR,kR,iBlockRecv) = &
                                         prolongation_high_order_amr(CoarseCell_III)
                                 enddo

                              end do ! iR
                           end do ! jR
                        end do ! kR
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

                     endif
                  end if ! UseTime

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

                  if(UseHighResChange .and. iSendStage == 3) then
                     iDir1 = 0; jDir1 = 0; kDir1 = 0
                     i5_ = max(5*Di,1); j5_ = max(5*Dj,1); k5_ =  max(5*Dk,1)
                     do kR = kRMin, kRMax, DkR
                        kS = kSMin + abs((kR+9)/kRatioRestr &
                             -           (kRMin+9)/kRatioRestr)
                        ! kDir = -1 if kR is even; kDir = 1 if kR is odd. 
                        ! kR may be negative (1-nK), kR+2*nK will be positive.
                        if(kRatioRestr == 2) kDir1 = 2*mod(kR+2*nK,2) - 1
                        do jR = jRMin, jRMax, DjR
                           jS = jSMin + abs((jR+9)/jRatioRestr &
                                -           (jRMin+9)/jRatioRestr)
                           if(jRatioRestr == 2) jDir1 = 2*mod(jR+2*nJ,2) - 1
                           do iR = iRMin, iRMax, DiR
                              iS = iSMin + abs((iR+9)/iRatioRestr &
                                   -           (iRMin+9)/iRatioRestr)
                              if(iRatioRestr == 2) iDir1 = 2*mod(iR+2*nI,2) - 1

                              do iVar = 1, nVar
                                 CoarseCell_III(1:i5_,1:j5_,1:k5_) = &
                                      State_VGB(iVar,&
                                      iS-2*iDir1:iS+2*iDir1:sign(1,iDir1),&
                                      jS-2*jDir1:jS+2*jDir1:sign(1,jDir1),&
                                      kS-2*kDir1:kS+2*kDir1:sign(1,kDir1),&
                                      iBlockSend)

                                 BufferS_I(iBufferS+iVar)=&
                                      prolongation_high_order_amr(CoarseCell_III) 
                                 
                              enddo
                              iBufferS = iBufferS + nVar
                           end do ! iR
                        end do ! jR
                     end do ! kR

                  else
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
                  endif ! UseHighResChange

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
    use BATL_geometry, ONLY: init_geometry, z_, IsPeriodic_D, rot_to_cart, &
         xyz_to_coord

    use ModMpi, ONLY: MPI_allreduce, MPI_REAL, MPI_MIN, MPI_MAX

    integer, parameter:: MaxBlockTest            = 200
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

    call test_switches
    call test_scalar

    if(nDim == 1) RETURN !------------------------
    call test_non_cartesian


    if(nG < 3) RETURN
    call test_high_order_cartesian
    call test_high_order_non_cartesian
    
  contains
    !----------------------------------------------------------------------
    subroutine test_switches
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
    end subroutine test_switches
    !----------------------------------------------------------------------

    subroutine test_scalar
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
    end subroutine test_scalar
    !----------------------------------------------------------------------

    subroutine test_non_cartesian
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
    end subroutine test_non_cartesian
    !----------------------------------------------------------------------

    subroutine test_high_order_cartesian
      real    :: ExactSolution, Error, ErrorTotal
      integer :: iCount, nCount, nRefineNode, iRefinement,nRefinement
      integer :: iNode_I(8)
      integer :: iNode1_I(8)
      logical :: DoTestMeOnly = .false.
      integer :: nPoly = 3
      
      if(nDim == 2) then
         nCount = 16; nRefineNode = 4
      else 
         nCount = 256; nRefineNode = 8
      endif

      if(nDimAmr < nDim) RETURN

      iMin = MinI; iMax = MaxI
      jMin = MinJ; jMax = MaxJ
      kMin = MinK; kMax = MaxK

      NameGeometry = 'cartesian'

      DomainMin_D = (/0.0, 0.0, 0.0/)
      DomainMax_D = (/8.0, 8.0, 8.0/)
      DomainSize_D = DomainMax_D - DomainMin_D

      IsPeriodicTest_D = (/.true., .true., .true./)
      nRootTest_D = (/4,4,4/)

      nRefinement = 2
      do iRefinement = 1, nRefinement
         ! iRefinement = 1: 1 level refine
         ! iRefinement = 2: 2 level refine
         
         if(DoTestMe)then
            write(*,*) &
                 'testing message_pass_cell across '//trim(NameGeometry)// &
                 ' with high resolution change with refinement level =', &
                 iRefinement
         end if

         if(nDim == 2) then
            if(iRefinement == 1) then
               iNode_I = (/6,7,10,11,-1,-1,-1,-1/)
            else
               iNode1_I = (/6,7,10,11,-1,-1,-1,-1/)
               iNode_I  = (/20,23,26,29,-1,-1,-1,-1/)
            endif
         else ! 3D
            if(iRefinement == 1) then
               iNode_I = (/22,23,26,27,38,39,42,43/)
            else
               iNode1_I = (/22,23,26,27,38,39,42,43/)
               iNode_I = (/72,79,86,93,100,107,114,121/)
            endif
         endif

         do iCount = 0, nCount-1
            if(DoTestMeOnly) then
               write(*,*) ''
               write(*,*) 'test_high_order iCount = ', iCount
            endif
            call init_tree(MaxBlockTest)
            call init_geometry(NameGeometry, &
                 IsPeriodicIn_D=IsPeriodicTest_D(1:nDim))

            call init_grid(DomainMin_D(1:nDim), DomainMax_D(1:nDim), &
                 UseDegreeIn=.false.)
            call set_tree_root( nRootTest_D(1:nDim))

            if(any(IsPeriodic_D(1:nDim) .neqv. IsPeriodicTest_D(1:nDim))) &
                 write(*,*) NameSub,': IsPeriodic_D=', IsPeriodic_D(1:nDim), &
                 ' should agree with ', IsPeriodicTest_D(1:nDim)

            if(iRefinement == 2) then
               do iNode = 1,nRefineNode
                  call refine_tree_node(iNode1_I(iNode))
               enddo
            endif

            do iNode = 1,nRefineNode
               if(btest(iCount,iNode-1)) then
                  call refine_tree_node(iNode_I(iNode))
                  if(DoTestMeOnly) &
                       write(*,*) 'iNode IsRefined:', iNode_I(iNode), 'TRUE'
               else
                  if(DoTestMeOnly) &
                       write(*,*) 'iNode IsRefined:', iNode_I(iNode), 'FALSE'
               endif
            enddo

            Tolerance = 5e-15

            call distribute_tree(.true.)
            call create_grid

            allocate(&
                 State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlockTest))
            State_VGB = 0

            do iBlock = 1, nBlock
               if(Unused_B(iBlock)) CYCLE
               do i = 1,nI; do j = 1,nJ; do k = 1,nK
                  State_VGB(1,i,j,k,iBlock) = &
                       exact_solution(Xyz_DGB(:,i,j,k,iBlock),nPolyIn=nPoly)
               enddo; enddo; enddo
            end do

            call message_pass_real(nVar,nG,State_VGB,nProlongOrderIn=1,&
                 nCoarseLayerIn=2,DoResChangeOnlyIn=.false., &
                 UseHighResChangeIn=.true.)
            ErrorTotal = 0
            do iBlock = 1, nBlock
               if(Unused_B(iBlock)) CYCLE

               ! Loop through all cells including ghost cells
               do k = kMin,kMax; do j = jMin,jMax; do i = iMin, iMax
                  Xyz_D = Xyz_DGB(:,i,j,k,iBlock)                  
                  if(.not. (all(Xyz_D(1:nDim) < DomainMax_D(1:nDim)) &
                       .and. all(Xyz_D(1:nDim) > DomainMin_D(1:nDim)))) then
                     CYCLE                
                  endif

                  ExactSolution= exact_solution(Xyz_D, nPolyIn=nPoly)
                  Error = abs(ExactSolution - State_VGB(1,i,j,k,iBlock))
                  ErrorTotal = ErrorTotal + Error
                  if(abs(Error)/abs(ExactSolution)>Tolerance)&
                       then
                     write(*,*)&
                          'iProc,iNode,i,j,k,x,y,z,',&
                          'state,exact-solution,error,relative-error='
                     write(*,'(5I5,7e20.12)')&
                          iProc,iNode_B(iBlock),i,j,k, &
                          Xyz_D,State_VGB(1,i,j,k,iBlock), &
                          ExactSolution, Error, abs(Error)/abs(ExactSolution)

                  end if
               end do; end do; end do
            end do
            if(DoTestMeOnly ) then
               write(*,*) 'Refine level = ', iRefinement
               write(*,*) 'Total error  = ', ErrorTotal
            endif
            deallocate(State_VGB)

            call clean_grid
            call clean_tree
         enddo ! iCount
      enddo ! iRefinement

    end subroutine test_high_order_cartesian
    !----------------------------------------------------------------------

    subroutine test_high_order_non_cartesian
      real    :: ErrorTotal, ExactSolution, Error
      real    :: xyz1_D(3), xyzGeneral_D(3)
      integer :: nPoly
      do iTest = 1,2

         ! The code is quite inaccurate for partial AMR across the pole
         if(nDimAmr < nDim .and. iTest > 3) EXIT

         call init_tree(MaxBlockTest)

         ! Do not test ghost cells in the radial direction
         iMin = 1; iMax = nI
         jMin = MinJ; jMax = MaxJ
         kMin = MinK; kMax = MaxK

         select case(iTest)
         case(1)

            NameGeometry = 'cylindrical'

            ! 2 < r < 8, 0 < phi < 360deg, -5 < z < 5
            DomainMin_D = (/ 2.0,  0.0, -5.0 /)
            DomainMax_D = (/ 8.0, cTwoPi, +5.0 /)
            IsPeriodicTest_D = (/ .false., .true., .true. /)

            ! There must be an even number of root blocks in the phi direction
            ! There are 3 root blocks in z so that we can refine the middle
            ! and avoid issues of periodicity in the testing
            nRootTest_D = (/3,4,3/)

            ! Test ghost cells at rMin
            iMin = MinI

         case(2)
            NameGeometry = 'rotatedcartesian'

            DomainMin_D = (/0.0, 0.0, 0.0/)
            DomainMax_D = (/6.0, 6.0, 6.0/)

            IsPeriodicTest_D = (/.false., .false., .false./)
            nRootTest_D = (/3,3,3/)
         end select
         DomainSize_D = DomainMax_D - DomainMin_D

         if(DoTestMe)then
            write(*,*) &
                 'testing message_pass_cell across '//trim(NameGeometry)// &
                 ' with high resolution change'
         end if

         call init_geometry(NameGeometry, &
              IsPeriodicIn_D=IsPeriodicTest_D(1:nDim))

         call init_grid(DomainMin_D(1:nDim), DomainMax_D(1:nDim), &
              UseDegreeIn=.false.)
         call set_tree_root( nRootTest_D(1:nDim))

         if(any(IsPeriodic_D(1:nDim) .neqv. IsPeriodicTest_D(1:nDim))) &
              write(*,*) NameSub,': IsPeriodic_D=', IsPeriodic_D(1:nDim), &
              ' should agree with ', IsPeriodicTest_D(1:nDim)

         if(iTest==1)then
            ! refine node next to r=0 axis but middle in Z direction
            if(nDim==2) then
               call refine_tree_node(2)
            elseif(nDim==3) then
               call refine_tree_node(17)
            endif
         elseif(iTest == 2) then
            if(nDim==2) then
               call refine_tree_node(5)
            elseif(nDim==3) then
               call refine_tree_node(14)
            endif
         end if

         if(iTest == 2) then
            Tolerance = 1e-14
            nPoly=3
         else 
            Tolerance = 7e-3
            nPoly=1
         endif

         call distribute_tree(.true.)
         call create_grid

         allocate(State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlockTest))
         State_VGB = 0

         do iBlock = 1, nBlock
            if(Unused_B(iBlock)) CYCLE
            do i = 1,nI; do j = 1,nJ; do k = 1,nK
               State_VGB(1,i,j,k,iBlock) = &
                    exact_solution(Xyz_DGB(:,i,j,k,iBlock),nPolyIn=nPoly)
            enddo; enddo; enddo
         end do

         call message_pass_real(nVar,nG,State_VGB,nProlongOrderIn=1,&
              nCoarseLayerIn=2,DoResChangeOnlyIn=.false., &
              UseHighResChangeIn=.true.)

         ! Second order
         ! call message_pass_cell(nVar, State_VGB)

         ErrorTotal = 0
         do iBlock = 1, nBlock
            if(Unused_B(iBlock)) CYCLE
            ! Loop through all cells including ghost cells
            do k = kMin,kMax; do j = jMin,jMax; do i = iMin, iMax

               Xyz_D = Xyz_DGB(:,i,j,k,iBlock)

               if(iTest == 2) then
                  xyz1_D = rot_to_cart(Xyz_D)
               else
                  call xyz_to_coord(Xyz_D, xyz1_D)
               endif

               if(.not. (all(Xyz1_D(1:nDim) < DomainMax_D(1:nDim)) &
                    .and. all(Xyz1_D(1:nDim) > DomainMin_D(1:nDim)))) then
                  CYCLE                
               endif

               ExactSolution= exact_solution(Xyz_D,nPolyIn=nPoly)
               Error = abs(ExactSolution - State_VGB(1,i,j,k,iBlock))
               ErrorTotal = ErrorTotal + Error/abs(ExactSolution)

               if(abs(Error)/abs(ExactSolution)> Tolerance)then
                  write(*,*)&
                       'iProc,iNode,i,j,k,x,y,z,',&
                       'state,exact-solution,error,relative-error='
                  write(*,'(5I5,7e20.12)')&
                       iProc,iNode_B(iBlock),i,j,k, &
                       Xyz_D,State_VGB(1,i,j,k,iBlock), &
                       ExactSolution, Error, abs(Error)/abs(ExactSolution)
                  call xyz_to_coord(Xyz_D, xyzGeneral_D)
                  write(*,*) 'xyz general = ',xyzGeneral_D
                  write(*,*) ''
               end if

            end do; end do; end do
         end do

         deallocate(State_VGB)

         call clean_grid
         call clean_tree
      end do ! iTest

    end subroutine test_high_order_non_cartesian

    !----------------------------------------------------------------------

    real function exact_solution(Xyz_D, nPolyIn)
      real, intent(in):: Xyz_D(3)
      integer, optional,intent(in) :: nPolyIn
      integer:: nPoly
      real:: x, y, z

      nPoly = 4
      if(present(nPolyIn)) nPoly = nPolyIn

      x = Xyz_D(1) 
      y = Xyz_D(2)
      z = Xyz_D(3)

      !exact_solution = 4.0+sin(y*2*cpi/8.0)+cos(z*2*cpi/8.0)+sin(x*2*cpi/8.0)
      select case(nPoly)
      case(4)
         exact_solution = x**4 +y**4+ z**4 + x**3*y &
              + y**3*x + x**3*z + x*z**3 + y**3*z + y*z**3 &
              + x**2*y**2 + x**2*z**2 + y**2*z**2 &
              + y*z*x**2 + x*z*y**2 + x*y*z**2
      case(3)
         exact_solution = x**3 + y**3 + z**3 + x*y*z + x**2*y + x*y**2 + &
              x**2*z + x*z**2 + y**2*z + y*z**2
      case(1)
         exact_solution = x + y + z
      end select
    end function exact_solution
    !======================================================================


  end subroutine test_pass_cell
end module BATL_pass_cell
