!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModImplHypre

  use ModKind,   ONLY: Int8_
  use BATL_size, ONLY: nDim, nI, nJ, nK, nIJK, nIJK_D, iRatio, jRatio, kRatio
  use BATL_lib,  ONLY: iProc, iNode_B, &
       iNodeNei_IIIB, DiLevelNei_IIIB, Unset_, iComm, &
       nRoot_D, MaxCoord_I, IsPeriodic_D, iTree_IA, &
       Proc_, Coord0_, Coord1_, Coord2_, Coord3_, Level_
  use ModImplicit, ONLY: impl2iBlk, nImplBlock => nImplBlk, &
       nStencil, Stencil1_, Stencil2_, Stencil3_, Stencil4_, Stencil5_, &
       Stencil6_, Stencil7_

  implicit none

  SAVE

  private ! except

  public:: hypre_read_param
  public:: hypre_initialize
  public:: hypre_finalize
  public:: hypre_set_matrix_block
  public:: hypre_set_matrix
  public:: hypre_preconditioner

  logical, public, parameter:: IsHypreAvailable = .true. ! true here
  logical, public:: DoInitHypreAmg = .true. ! set to false to reinitialize AMG

  ! local variables

  integer, parameter:: CoordLast_ = Coord0_+nDim
  integer, parameter:: nCell_D(nDim) = nIjk_D(1:nDim)

  ! This is defined in HYPREf.h (Fortran header file)                         
  integer, parameter:: HYPRE_PARCSR = 5555

  ! This is defined in HYPRE_sstruct_mv.h (a C header file)
  integer, parameter:: HYPRE_SSTRUCT_VARIABLE_CELL = 0

  ! Hypre uses C type 8-byte integers for pointers
  integer(Int8_):: i8Grid, i8Graph, i8Stencil, i8Precond
  integer(Int8_):: i8A, i8B, i8X, i8ParA, i8ParB, i8ParX

  ! These indexes are never out of the 1:nDim range
  integer, parameter:: Dim1_ = 1, Dim2_ = min(2,nDim), Dim3_=nDim

  integer:: nPart = 1 ! grid consists of one or more parts
  integer:: nVar  = 1 ! there is only one variable per grid cell
  integer:: iVar  = 0 ! index of the variable
  integer:: iVarType_I(1) = (/ HYPRE_SSTRUCT_VARIABLE_CELL /)
  integer:: iLower_D(nDim), iUpper_D(nDim) ! index limits for each box

  integer:: iLowerBc_D(nDim), iUpperBc_D(nDim) ! index limit for ghost cells
  integer:: jLowerBc_D(nDim), jUpperBc_D(nDim) ! index limit for phys. cells

  ! Mapping between boundaries: same index order and same direction
  integer:: iDim
  integer, parameter:: iIndexMap_D(nDim) = (/ (iDim, iDim=0,nDim-1)/) 
  integer, parameter:: iIndexDir_D(nDim) = 1

  ! Stencil description
  integer:: iStencil
  integer:: iStencil_I(nStencil) = (/ (iStencil, iStencil=0,nStencil-1) /)
  integer:: DiStencil_DI(nDim,nStencil)

  integer:: iObjectType = HYPRE_PARCSR

  real, allocatable:: Value_I(:) ! matrix elements for 1 block

  ! BoomerAMG parameters (the defaults are optimal for high res. 1D CRASH)
  integer:: iVerboseAmg        = 0    ! 0..3
  integer:: MaxRowElementsAmg  = 0    ! 2, 4 or 6 (2*nDim)
  integer:: iCoarsenAmg        =10    ! 0,1,3,6,7,8,9,10,11,21,22
  integer:: iRelaxAmg          = 6    ! 0..6,8,9,15..18
  integer:: iInterpolateAmg    = 6    ! 0..14
  real::    StrongThresholdAmg = 0.0  ! 0.25 for 2D, 0.5-0.6 for 3D
  real::    TruncFactorAmg     = 0.0  ! 0..1 ?

  ! Check if the grid has changed
  integer:: iLastDecomposition = -1

contains
  !==========================================================================
  subroutine hypre_read_param

    use ModReadParam, ONLY: read_var

    call read_var('iVerboseAmg',        iVerboseAmg)
    call read_var('MaxRowElementsAmg',  MaxRowElementsAmg)
    call read_var('iCoarsenAmg',        iCoarsenAmg)
    call read_var('iRelaxAmg',          iRelaxAmg)
    call read_var('iInterpolateAmg',    iInterpolateAmg)
    call read_var('StrongThresholdAmg', StrongThresholdAmg)
    call read_var('TruncFactorAmg',     TruncFactorAmg)

  end subroutine hypre_read_param
  !==========================================================================
  subroutine hypre_initialize

    use ModMain, ONLY: iNewDecomposition
    use BATL_lib, ONLY: nLevel

    integer:: iLevel, DiLevel, iNode, jNode
    integer:: iCoord_D(nDim), jCoord_D(nDim), iCoord0_D(nDim), jCoord0_D(nDim)

    integer:: iImplBlock, iBlock, iError, iPart, jPart

    integer:: iReduce, jReduce, iSide, jSide, kSide, i, j, k
    integer:: iSideMe, jSideMe, kSideMe, iSideMax, jSideMax, kSideMax

    logical:: DoTest, DoTestMe

    logical, parameter:: DoDebug = .false.

    character(len=*), parameter:: NameSub = 'hypre_initialize'
    !-------------------------------------------------------------------------
    ! Nothing to do if the grid has not changed
    if(iNewDecomposition == iLastDecomposition) RETURN

    call set_oktest(NameSub,DoTest,DoTestMe)

    if(DoTestMe)write(*,*) NameSub,' starting'

    if(iLastDecomposition == -1)then
       allocate(Value_I(nStencil*nIJK))
    else
       ! This is not the first call, so destroy previously created HYPRE objects
       call hypre_clean
    end if
    iLastDecomposition = iNewDecomposition

    ! One part per level. We can possibly use MaxLevel (<=30), or count
    ! the number of current levels (1-5 or so).

    ! Level indexes go from 0 to nLevel, hence the +1
    nPart = nLevel + 1

    ! Create an empty 3D grid object
    call HYPRE_SStructGridCreate(iComm, nDim, nPart, i8Grid, iError)

    if(DoTestMe)write(*,*)'HYPRE_SStructGridCreate done with nPart=', nPart

    ! Add each block as a local box in the corresponding part (level)
    do iImplBlock = 1, nImplBlock
       iBlock   = impl2iBlk(iImplBlock)
       iNode    = iNode_B(iBlock)
       iPart    = iTree_IA(Level_,iNode)
       iLower_D = 1 + (iTree_IA(Coord1_:CoordLast_,iNode)-1)*nCell_D
       iUpper_D = iLower_D - 1 + nCell_D

       if(DoDebug)write(*,*)'iPart, iNode, iLower_D, iUpper_D=',&
            iPart, iNode, iLower_D, iUpper_D

       call HYPRE_SStructGridSetExtents(i8Grid, iPart, iLower_D, iUpper_D, &
            iError)
    end do

    if(DoTestMe)write(*,*)'HYPRE_SStructGridSetExtents done'

    ! Single cell centered variable on all parts
    do jPart = 0, nPart - 1
       call HYPRE_SStructGridSetVariables(i8Grid, jPart, nVar, iVarType_I, &
            iError)
    end do

    if(DoTestMe)write(*,*)'HYPRE_SStructGridSetVariables done'

    ! For periodic boundaries the part is connected to itself
    do iDim = 1, nDim
       if(.not.IsPeriodic_D(iDim)) CYCLE

       do iLevel = 0, nLevel
          iPart = iLevel
          iLower_D = 1
          iUpper_D = MaxCoord_I(iLevel)*nRoot_D(1:nDim)*nCell_D

          ! initialize boundary ranges
          iLowerBc_D = iLower_D
          iUpperBc_D = iUpper_D
          jLowerBc_D = iLower_D
          jUpperBc_D = iUpper_D

          ! lower ghost cell connected to last cell
          iLowerBc_D(iDim) = 0
          iUpperBc_D(iDim) = 0
          jLowerBc_D(iDim) = iUpper_D(iDim)
          jUpperBc_D(iDim) = iUpper_D(iDim)

          call HYPRE_SStructGridSetNeighborPart( i8Grid, &
               iPart, iLowerBc_D, iUpperBc_D, &
               iPart, jLowerBc_D, jUpperBc_D, &
               iIndexMap_D, iIndexDir_D, iError)

          ! upper ghost cell connected to first cell
          iLowerBc_D(iDim) = iUpper_D(iDim) + 1
          iUpperBc_D(iDim) = iUpper_D(iDim) + 1
          jLowerBc_D(iDim) = 1
          jUpperBc_D(iDim) = 1

          call HYPRE_SStructGridSetNeighborPart( i8Grid, &
               iPart, iLowerBc_D, iUpperBc_D, &
               iPart, jLowerBc_D, jUpperBc_D, &
               iIndexMap_D, iIndexDir_D, iError)

       end do
    end do

    if(DoTestMe)write(*,*)'HYPRE_SStructGridSetNeighborPart done'

    ! Assemble grid from all processors
    call HYPRE_SStructGridAssemble(i8Grid, iError)
    if(iError/=0)write(*,*)'ERROR: HYPRE_SStructGridAssemble failed'
    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructGridAssemble done'

    ! Define index offsets for the 2*nDim+1-point stencil
    DiStencil_DI = 0
    DiStencil_DI(           Dim1_,Stencil2_) = -1
    DiStencil_DI(           Dim1_,Stencil3_) = +1
    if(nJ > 1) DiStencil_DI(Dim2_,Stencil4_) = -1
    if(nJ > 1) DiStencil_DI(Dim2_,Stencil5_) = +1
    if(nK > 1) DiStencil_DI(Dim3_,Stencil6_) = -1
    if(nK > 1) DiStencil_DI(Dim3_,Stencil7_) = +1

    call HYPRE_SStructStencilCreate(nDim, nStencil, i8Stencil, iError)
    do iStencil = 1, nStencil
       call HYPRE_SStructStencilSetEntry(i8Stencil, &
            iStencil-1, DiStencil_DI(1,iStencil), iVar, iError)
    enddo

    ! Create the graph object
    call HYPRE_SStructGraphCreate(iComm, i8Grid, i8Graph, iError)
    call HYPRE_SStructGraphSetObjectType(i8Graph, iObjectType, iError)

    ! Tell the graph which stencil to use for each variable on each part 
    do jPart = 0, nPart - 1
       call HYPRE_SStructGraphSetStencil(i8Graph, jPart, iVar, i8Stencil, &
            iError)
    end do

    ! Add non-stencil entries to the graph at resolution changes

    do iImplBlock = 1, nImplBlock

       iBlock = impl2iBlk(iImplBlock)
       iNode  = iNode_B(iBlock)
       iPart  = iTree_IA(Level_,iNode)

       ! Which side of the parent block
       iSideMe = modulo(iTree_IA(Coord1_,iNode)-1,2) + 1
       jSideMe = modulo(iTree_IA(Coord2_,iNode)-1,2) + 1
       kSideMe = modulo(iTree_IA(Coord3_,iNode)-1,2) + 1

       ! -I neighbor
       DiLevel = DiLevelNei_IIIB(-1,0,0,iBlock)
       if(abs(DiLevel) == 1)then
          call set_jpart_iratio_jratio

          ! Loop through neighboring node(s)
          do kSide = 1, kSideMax; do jSide = 1, jSideMax
             jNode = iNodeNei_IIIB(0,jSide,kSide,iBlock)

             if(iProc /= iTree_IA(Proc_,iNode) &
                  .and. iProc /= iTree_IA(Proc_,jNode)) CYCLE

             ! Global cell index of the first cell of iNode
             iCoord0_D = (iTree_IA(Coord1_:CoordLast_,iNode) - 1)*nCell_D + 1

             ! Global cell index of the last cell of jNode
             jCoord0_D = iTree_IA(Coord1_:CoordLast_,jNode)*nCell_D

             ! Shift back to the lower corner in the J and K dimensions
             if(nJ > 1) jCoord0_D(Dim2_) = jCoord0_D(Dim2_) - nJ + 1
             if(nK > 1) jCoord0_D(Dim3_) = jCoord0_D(Dim3_) - nK + 1

             if(DoDebug)write(*,*)'Connect -I direction'
             call connect_i_direction

          end do; end do
       end if

       ! +I neighbor
       DiLevel = DiLevelNei_IIIB(+1,0,0,iBlock)
       if(abs(DiLevel) == 1)then
          call set_jpart_iratio_jratio

          ! Loop through neighboring node(s)
          do kSide = 1, kSideMax; do jSide = 1, jSideMax
             jNode = iNodeNei_IIIB(3,jSide,kSide,iBlock)

             if(iProc /= iTree_IA(Proc_,iNode) &
                  .and. iProc /= iTree_IA(Proc_,jNode)) CYCLE

             ! Global cell index of the last cell of iNode
             iCoord0_D = iTree_IA(Coord1_:CoordLast_,iNode)*nCell_D

             ! Global cell index of the first cell of jNode
             jCoord0_D = (iTree_IA(Coord1_:CoordLast_,jNode) - 1)*nCell_D + 1

             ! Shift back to the lower corner in the J and K dimensions
             if(nJ > 1) iCoord0_D(Dim2_) = iCoord0_D(Dim2_) - nJ + 1
             if(nK > 1) iCoord0_D(Dim3_) = iCoord0_D(Dim3_) - nK + 1

             if(DoDebug)write(*,*)'Connect +I direction'
             call connect_i_direction

          end do; end do
       end if

       ! Done with the block for 1D
       if(nJ == 1) CYCLE

       ! -J neighbor
       DiLevel = DiLevelNei_IIIB(0,-1,0,iBlock)
       if(abs(DiLevel) == 1)then
          call set_jpart_iratio_jratio

          ! Loop through neighboring node(s)
          do kSide = 1, kSideMax; do iSide = 1, iSideMax
             jNode = iNodeNei_IIIB(iSide,0,kSide,iBlock)

             if(iProc /= iTree_IA(Proc_,iNode) &
                  .and. iProc /= iTree_IA(Proc_,jNode)) CYCLE

             ! Global cell index of the first cell of iNode
             iCoord0_D = (iTree_IA(Coord1_:CoordLast_,iNode) - 1)*nCell_D + 1

             ! Global cell index of the last cell of jNode
             jCoord0_D = iTree_IA(Coord1_:CoordLast_,jNode)*nCell_D

             ! Shift back to the lower corner in the I and K dimensions
             jCoord0_D(Dim1_)          = jCoord0_D(Dim1_) - nI + 1
             if(nK>1) jCoord0_D(Dim3_) = jCoord0_D(Dim3_) - nK + 1

             if(DoDebug)write(*,*)'Connect -J direction'
             call connect_j_direction

          end do; end do
       end if

       ! +J neighbor
       DiLevel = DiLevelNei_IIIB(0,+1,0,iBlock)
       if(abs(DiLevel) == 1)then
          call set_jpart_iratio_jratio

          ! Loop through neighboring node(s)
          do kSide = 1, kSideMax; do iSide = 1, jSideMax
             jNode = iNodeNei_IIIB(iSide,3,kSide,iBlock)

             if(iProc /= iTree_IA(Proc_,iNode) &
                  .and. iProc /= iTree_IA(Proc_,jNode)) CYCLE

             ! Global cell index of the last cell of iNode
             iCoord0_D = iTree_IA(Coord1_:CoordLast_,iNode)*nCell_D

             ! Global cell index of the first cell of jNode
             jCoord0_D = (iTree_IA(Coord1_:CoordLast_,jNode) - 1)*nCell_D + 1

             ! Shift back to the lower corner in the J and K dimensions
             iCoord0_D(Dim1_)          = iCoord0_D(Dim1_) - nI + 1
             if(nK>1) iCoord0_D(Dim3_) = iCoord0_D(Dim3_) - nK + 1

             if(DoDebug)write(*,*)'Connect +J direction'
             call connect_j_direction

          end do; end do
       end if

       ! Done with the block in 2D
       if(nK == 1) CYCLE

       ! -K neighbor
       DiLevel = DiLevelNei_IIIB(0,0,-1,iBlock)
       if(abs(DiLevel) == 1)then
          call set_jpart_iratio_jratio

          ! Loop through neighboring node(s)
          do jSide = 1, jSideMax; do iSide = 1, iSideMax
             jNode = iNodeNei_IIIB(iSide,jSide,0,iBlock)

             if(iProc /= iTree_IA(Proc_,iNode) &
                  .and. iProc /= iTree_IA(Proc_,jNode)) CYCLE

             ! Global cell index of the first cell of iNode
             iCoord0_D = (iTree_IA(Coord1_:CoordLast_,iNode) - 1)*nCell_D + 1

             ! Global cell index of the last cell of jNode
             jCoord0_D = iTree_IA(Coord1_:CoordLast_,jNode)*nCell_D

             ! Shift back to the lower corner in the J and K dimensions
             jCoord0_D(Dim1_) = jCoord0_D(Dim1_) - nI + 1
             jCoord0_D(Dim2_) = jCoord0_D(Dim2_) - nJ + 1

             if(DoDebug)write(*,*)'Connect -K direction'
             call connect_k_direction

          end do; end do
       end if

       ! +K neighbor
       DiLevel = DiLevelNei_IIIB(0,0,+1,iBlock)
       if(abs(DiLevel) == 1)then
          call set_jpart_iratio_jratio

          ! Loop through neighboring node(s)
          do jSide = 1, jSideMax; do iSide = 1, iSideMax
             jNode = iNodeNei_IIIB(iSide,jSide,3,iBlock)

             if(iProc /= iTree_IA(Proc_,iNode) &
                  .and. iProc /= iTree_IA(Proc_,jNode)) CYCLE

             ! Global cell index of the last cell of iNode
             iCoord0_D = iTree_IA(Coord1_:CoordLast_,iNode)*nCell_D

             ! Global cell index of the first cell of jNode
             jCoord0_D = (iTree_IA(Coord1_:CoordLast_,jNode) - 1)*nCell_D + 1

             ! Shift back to the lower corner in the J and K dimensions
             iCoord0_D(Dim1_) = iCoord0_D(Dim1_) - nI + 1
             iCoord0_D(Dim2_) = iCoord0_D(Dim2_) - nJ + 1

             if(DoDebug)write(*,*)'Connect +K direction'
             call connect_k_direction

          end do; end do
       end if

    end do

    ! Assemble the graph
    call HYPRE_SStructGraphAssemble(i8Graph, iError)

    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructGraphAssemble done'

    ! Create an empty matrix object
    call HYPRE_SStructMatrixCreate(iComm, i8Graph, i8A, iError)

    ! Set storage type
    call HYPRE_SStructMatrixSetObjectTyp(i8A, iObjectType, iError)

    ! Initialize matrix
    call HYPRE_SStructMatrixInitialize(i8A, iError)
    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructMatrixInitialize done'

    ! Create empty vector objects for RHS and solution
    call HYPRE_SStructVectorCreate(iComm, i8Grid, i8B, iError)
    call HYPRE_SStructVectorCreate(iComm, i8Grid, i8X, iError)

    ! Set object type for the vectors
    call HYPRE_SStructVectorSetObjectTyp(i8B, iObjectType, iError)
    call HYPRE_SStructVectorSetObjectTyp(i8X, iObjectType, iError)

    ! Initialize vectors
    call HYPRE_SStructVectorInitialize(i8B, iError)
    call HYPRE_SStructVectorInitialize(i8X, iError)
    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructVectorInitialize done'

    ! Create the BoomerAMG as a preconditioner
    call HYPRE_BoomerAMGCreate(i8Precond, iError)
    if(DoTestMe)write(*,*) NameSub,' HYPRE_BoomerAMGCreate done'

    ! Set BoomerAMG parameters

    ! As a preconditioner always do one sweep
    call HYPRE_BoomerAMGSetMaxIter(i8Precond, 1, iError)
    call HYPRE_BoomerAMGSetTol(i8Precond,   0.0, iError)

    ! These parameters are adjustables with the #HYPRE command
    call HYPRE_BoomerAMGSetPrintLevel(   i8Precond, iVerboseAmg, iError)
    call HYPRE_BoomerAMGSetPMaxElmts(    i8Precond, MaxRowElementsAmg, iError)
    call HYPRE_BoomerAMGSetCoarsenType(  i8Precond, iCoarsenAmg, iError)
    call HYPRE_BoomerAMGSetRelaxType(    i8Precond, iRelaxAmg, iError)
    call HYPRE_BoomerAMGSetInterpType(   i8Precond, iInterpolateAmg, iError)
    call HYPRE_BoomerAMGSetStrongThrshld(i8Precond, StrongThresholdAmg, iError)
    call HYPRE_BoomerAMGSetTruncFactor(  i8Precond, TruncFactorAmg, iError)

    if(DoTestMe)write(*,*) NameSub,' finished'

  contains
    !========================================================================
    subroutine set_jpart_iratio_jratio

      if(DiLevel == 1)then
         ! jNode is coarser than iNode
         jPart   = iPart - 1
         iReduce = 1
         jReduce = 2
         ! There are no subfaces
         iSideMax = 1
         jSideMax = 1
         kSideMax = 1
      else
         ! jNode is finer than iNode
         jPart   = iPart + 1
         iReduce = 2
         jReduce = 1
         ! Set up ranges for the subfaces
         iSideMax = iRatio
         jSideMax = jRatio
         kSideMax = kRatio
      end if

    end subroutine set_jpart_iratio_jratio
    !========================================================================
    subroutine connect_i_direction

      ! Shift index ranges to appropriate subface
      if(iReduce == 2)then
         ! iNode is coarser than jNode: do current subface of iNode
         if(nJ>1) iCoord0_D(Dim2_) = iCoord0_D(Dim2_) + (jSide-1)*nJ/2
         if(nK>1) iCoord0_D(Dim3_) = iCoord0_D(Dim3_) + (kSide-1)*nK/2
      else
         ! iNode is finer than jNode: do the corresponding subface of jNode
         if(nJ>1) jCoord0_D(Dim2_) = jCoord0_D(Dim2_) + (jSideMe-1)*nJ/2
         if(nK>1) jCoord0_D(Dim3_) = jCoord0_D(Dim3_) + (kSideMe-1)*nK/2
      end if

      if(DoDebug)then
         write(*,*)'iPart, jPart   =', iPart, jPart
         write(*,*)'iNode, jNode   =', iNode, jNode
         write(*,*)'iCoord0,jCoord0=', iCoord0_D, jCoord0_D
      end if

      iCoord_D = iCoord0_D
      jCoord_D = jCoord0_D
      do k = 1, nK
         ! On the coarse side k index changes every second time only
         if(nK>1) iCoord_D(Dim3_) = iCoord0_D(Dim3_) + (k - 1)/iReduce
         if(nK>1) jCoord_D(Dim3_) = jCoord0_D(Dim3_) + (k - 1)/jReduce
         do j = 1,nJ
            ! On the coarse side j index changes every second time only
            if(nJ>1) iCoord_D(Dim2_) = iCoord0_D(Dim2_) + (j - 1)/iReduce
            if(nJ>1) jCoord_D(Dim2_) = jCoord0_D(Dim2_) + (j - 1)/jReduce

            if(DoDebug)write(*,*)'iCoord,jCoord=',iCoord_D, jCoord_D

            ! Add the connecting graph entry
            call HYPRE_SStructGraphAddEntries(i8Graph, iPart, &
                 iCoord_D, iVar, jPart, jCoord_D, iVar, iError)

         end do
      end do

    end subroutine connect_i_direction
    !========================================================================
    subroutine connect_j_direction

      ! Shift to appropriate side
      if(iReduce == 2)then
         iCoord0_D(Dim1_)          = iCoord0_D(Dim1_) + (iSide-1)*nI/2
         if(nK>1) iCoord0_D(Dim3_) = iCoord0_D(Dim3_) + (kSide-1)*nK/2
      else
         jCoord0_D(Dim1_)          = jCoord0_D(Dim1_) + (iSideMe-1)*nI/2
         if(nK>1) jCoord0_D(Dim3_) = jCoord0_D(Dim3_) + (kSideMe-1)*nK/2
      end if

      if(DoDebug)then
         write(*,*)'iPart, jPart   =', iPart, jPart
         write(*,*)'iNode, jNode   =', iNode, jNode
         write(*,*)'iCoord0,jCoord0=', iCoord0_D, jCoord0_D
      end if

      iCoord_D = iCoord0_D
      jCoord_D = jCoord0_D
      do k = 1, nK
         if(nK>1) iCoord_D(Dim3_) = iCoord0_D(Dim3_) + (k - 1)/iReduce
         if(nK>1) jCoord_D(Dim3_) = jCoord0_D(Dim3_) + (k - 1)/jReduce
         do i = 1,nI
            iCoord_D(Dim1_) = iCoord0_D(Dim1_) + (i - 1)/iReduce
            jCoord_D(Dim1_) = jCoord0_D(Dim1_) + (i - 1)/jReduce

            if(DoDebug)write(*,*)'iCoord,jCoord=',iCoord_D, jCoord_D

            call HYPRE_SStructGraphAddEntries(i8Graph, iPart, &
                 iCoord_D, iVar, jPart, jCoord_D, iVar, iError)

         end do
      end do

    end subroutine connect_j_direction
    !========================================================================
    subroutine connect_k_direction

      ! Shift to appropriate side
      if(iReduce == 2)then
         iCoord0_D(Dim1_) = iCoord0_D(Dim1_) + (iSide-1)*nI/2
         iCoord0_D(Dim2_) = iCoord0_D(Dim2_) + (jSide-1)*nJ/2
      else
         jCoord0_D(Dim1_) = jCoord0_D(Dim1_) + (iSideMe-1)*nI/2
         jCoord0_D(Dim2_) = jCoord0_D(Dim2_) + (jSideMe-1)*nJ/2
      end if

      if(DoDebug)then
         write(*,*)'iPart, jPart   =', iPart, jPart
         write(*,*)'iNode, jNode   =', iNode, jNode
         write(*,*)'iCoord0,jCoord0=', iCoord0_D, jCoord0_D
      end if

      iCoord_D = iCoord0_D
      jCoord_D = jCoord0_D
      do j = 1,nJ
         iCoord_D(Dim2_) = iCoord0_D(Dim2_) + (j - 1)/iReduce
         jCoord_D(Dim2_) = jCoord0_D(Dim2_) + (j - 1)/jReduce
         do i = 1,nI
            iCoord_D(Dim1_) = iCoord0_D(Dim1_) + (i - 1)/iReduce
            jCoord_D(Dim1_) = jCoord0_D(Dim1_) + (i - 1)/jReduce

            if(DoDebug)write(*,*)'iCoord,jCoord=',iCoord_D, jCoord_D

            call HYPRE_SStructGraphAddEntries(i8Graph, iPart, &
                 iCoord_D, iVar, jPart, jCoord_D, iVar, iError)

         end do
      end do

    end subroutine connect_k_direction

  end subroutine hypre_initialize

  !===========================================================================

  subroutine hypre_set_matrix_block(iImplBlock, Jacobian_CI)

    ! Set the maxtrix elements corresponding to iImplBlock

    use ModMain, ONLY: TypeBc_I

    integer, intent(in):: iImplBlock
    real,    intent(inout):: Jacobian_CI(nI,nJ,nK,nStencil)

    ! Matrix element
    real   :: Jac

    ! Index of next extra connection for every cell in the block
    integer:: iStencilExtra_C(nI,nJ,nK)

    ! Index of first connection and number of extra connection to a single cell
    integer:: iStencilFirst, nStencilExtra

    integer:: iValue, i, j, k, iPart, iBlock, iError
    integer:: DiLevel

    ! Node index and global cell index
    integer:: iNode, iCoord_D(nDim)

    logical, parameter :: DoDebug = .false.
    !------------------------------------------------------------------------
    iBlock = impl2iblk(iImplBlock)

    ! DoDebug = iProc == ProcTest

    if(TypeBc_I(1)=='reflect' .and. &
         DiLevelNei_IIIB(-1,0,0,iBlock) == Unset_)&
         Jacobian_CI(1,:,:,Stencil1_) = Jacobian_CI(1,:,:,Stencil1_) &
         + Jacobian_CI(1,:,:,Stencil2_)

    if(TypeBc_I(2)=='reflect' .and. &
         DiLevelNei_IIIB(+1,0,0,iBlock) == Unset_)&
         Jacobian_CI(nI,:,:,Stencil1_) = Jacobian_CI(nI,:,:,Stencil1_) &
         + Jacobian_CI(nI,:,:,Stencil3_)

    if(nJ > 1)then
       if(TypeBc_I(3)=='reflect' .and. &
            DiLevelNei_IIIB(0,-1,0,iBlock) == Unset_)&
            Jacobian_CI(:,1,:,Stencil1_) = Jacobian_CI(:,1,:,Stencil1_) &
            + Jacobian_CI(:,1,:,Stencil4_)

       if(TypeBc_I(4)=='reflect' .and. &
            DiLevelNei_IIIB(0,+1,0,iBlock) == Unset_)&
            Jacobian_CI(:,nJ,:,Stencil1_) = Jacobian_CI(:,nJ,:,Stencil1_) &
            + Jacobian_CI(:,nJ,:,Stencil5_)
    end if

    if(nK > 1)then
       if(TypeBc_I(5)=='reflect' .and. &
            DiLevelNei_IIIB(0,0,-1,iBlock) == Unset_)&
            Jacobian_CI(:,:,1,Stencil1_) = Jacobian_CI(:,:,1,Stencil1_) &
            + Jacobian_CI(:,:,1,Stencil6_)

       if(TypeBc_I(6)=='reflect' .and. &
            DiLevelNei_IIIB(0,0,+1,iBlock) == Unset_)&
            Jacobian_CI(:,:,1,Stencil1_) = Jacobian_CI(:,:,1,Stencil1_) &
            + Jacobian_CI(:,:,1,Stencil7_)
    end if

    ! Reorder matrix elements so that stencil comes first
    iValue = 0
    do k = 1, nK; do j = 1, nJ; do i = 1, nI; do iStencil = 1, nStencil
       iValue = iValue + 1

       Value_I(iValue) = Jacobian_CI(i,j,k,iStencil)

    end do; end do; end do; end do

    iNode    = iNode_B(iBlock)
    iPart    = iTree_IA(Level_,iNode)
    iLower_D = 1 + (iTree_IA(Coord1_:CoordLast_,iNode)-1)*nCell_D
    iUpper_D = iLower_D - 1 + nCell_D

    if(DoDebug)write(*,*) &
         'iPart, iImplBlk, iBlock, iNode, iLower, iUpper, maxval, maxloc=',&
         iPart, iImplBlock, iBlock, iNode, iLower_D, iUpper_D, &
         maxval(Value_I), maxloc(Value_I)

    call HYPRE_SStructMatrixSetBoxValues(i8A, iPart, iLower_D, iUpper_D, &
         iVar, nStencil, iStencil_I, Value_I, iError)

    ! Initialize index array for extra matrix elements at resolution changes
    iStencilExtra_C = nStencil

    ! check -I neighbor for resolution change
    DiLevel = DiLevelNei_IIIB(-1,0,0,iBlock)
    if( abs(DiLevel) == 1)then

       ! Number of connections per cell for coarser and finer neighbor
       nStencilExtra = 1
       if(DiLevel == -1) nStencilExtra = jRatio*kRatio

       if(DoDebug)write(*,*)'-I iProc, iPart, DiLevel=',iProc, iPart, DiLevel

       iCoord_D = iLower_D
       do k = 1, nK
          if(nK > 1)iCoord_D(Dim3_) = iLower_D(Dim3_) + k - 1
          do j = 1, nJ
             if(nJ > 1)iCoord_D(Dim2_) = iLower_D(Dim2_) + j - 1

             Jac = Jacobian_CI(1,j,k,Stencil2_)/nStencilExtra
             do iStencil = nStencil, nStencil + nStencilExtra - 1
                call HYPRE_SStructMatrixSetValues(i8A, iPart, iCoord_D, &
                     iVar, 1, (/iStencil/), (/Jac/), iError)
             end do
             iStencilExtra_C(1,j,k) = iStencil

             if(DoDebug)write(*,*)'-I iCoord, Jac =', iCoord_D, Jac
             
          end do
       end do
    end if

    ! check +I neighbor for resolution change
    DiLevel = DiLevelNei_IIIB(+1,0,0,iBlock)
    if(abs(DiLevel) == 1)then

       ! Number of connections per cell for coarser and finer neighbor
       nStencilExtra = 1
       if(DiLevel == -1) nStencilExtra = jRatio*kRatio

       if(DoDebug)write(*,*)'+I iProc, iPart, DiLevel=',iProc, iPart, DiLevel

       iCoord_D = iUpper_D
       do k = 1, nK
          if(nK > 1) iCoord_D(Dim3_) = iUpper_D(Dim3_) + k - nK
          do j = 1, nJ
             if(nJ > 1) iCoord_D(Dim2_) = iUpper_D(Dim2_) + j - nJ
                
             Jac = Jacobian_CI(nI,j,k,Stencil3_)/nStencilExtra
             do iStencil = nStencil, nStencil + nStencilExtra - 1
                call HYPRE_SStructMatrixSetValues(i8A, iPart, &
                     iCoord_D, iVar, 1, (/iStencil/), (/Jac/), iError)
             end do
             iStencilExtra_C(nI,j,k) = iStencil

             if(DoDebug)write(*,*)'+I iCoord, Jac =', iCoord_D, Jac

          end do
       end do

    end if

    ! Done with block in 1D
    if(nJ == 1) RETURN

    ! check -J neighbor for resolution change
    DiLevel = DiLevelNei_IIIB(0,-1,0,iBlock)
    if(abs(DiLevel) == 1)then

       ! Number of connections per cell for coarser and finer neighbor
       nStencilExtra = 1
       if(DiLevel == -1) nStencilExtra = iRatio*kRatio

       if(DoDebug)write(*,*)'-J iProc, iPart, DiLevel=',iProc, iPart, DiLevel

       iCoord_D = iLower_D
       do k = 1, nK
          if(nK > 1) iCoord_D(Dim3_) = iLower_D(Dim3_) + k - 1
          do i = 1, nI
             iCoord_D(Dim1_) = iLower_D(Dim1_) + i - 1

             Jac = Jacobian_CI(i,1,k,Stencil4_)/nStencilExtra
             iStencilFirst = iStencilExtra_C(i,1,k)
             do iStencil = iStencilFirst, iStencilFirst + nStencilExtra - 1
                call HYPRE_SStructMatrixSetValues(i8A, iPart, iCoord_D, &
                     iVar, 1, (/iStencil/), (/Jac/), iError)
             end do
             iStencilExtra_C(i,1,k) = iStencil

             if(DoDebug)write(*,*)'-J iCoord, Jac =', iCoord_D, Jac
             
          end do
       end do

    end if

    ! check +J neighbor for resolution change
    DiLevel = DiLevelNei_IIIB(0,+1,0,iBlock)
    if(abs(DiLevel) == 1)then

       ! Number of connections per cell for coarser and finer neighbor
       nStencilExtra = 1
       if(DiLevel == -1) nStencilExtra = iRatio*kRatio

       if(DoDebug)write(*,*)'+J iProc, iPart, DiLevel=',iProc, iPart, DiLevel

       iCoord_D = iUpper_D
       do k = 1, nK
          if(nK > 1) iCoord_D(Dim3_) = iUpper_D(Dim3_) + k - nK
          do i = 1, nI
             iCoord_D(Dim1_) = iUpper_D(Dim1_) + i - nI

             Jac = Jacobian_CI(i,nJ,k,Stencil5_)/nStencilExtra
             iStencilFirst = iStencilExtra_C(i,nJ,k)
             do iStencil = iStencilFirst, iStencilFirst + nStencilExtra - 1
                call HYPRE_SStructMatrixSetValues(i8A, iPart, &
                     iCoord_D, iVar, 1, (/iStencil/), (/Jac/), iError)
             end do
             iStencilExtra_C(i,nJ,k) = iStencil

             if(DoDebug)write(*,*)'+J iCoord, Jac =', iCoord_D, Jac

          end do
       end do

    end if

    ! Done with block in 2D
    if(nK == 1) RETURN

    ! check -K neighbor for resolution change
    DiLevel = DiLevelNei_IIIB(0,0,-1,iBlock)
    if(abs(DiLevel) == 1)then

       ! Number of connections per cell for coarser and finer neighbor
       nStencilExtra = 1
       if(DiLevel == -1) nStencilExtra = iRatio*jRatio

       if(DoDebug)write(*,*)'-K iProc, iPart, DiLevel=',iProc, iPart, DiLevel

       iCoord_D = iLower_D
       do j = 1, nJ
          iCoord_D(Dim2_) = iLower_D(Dim2_) + j - 1
          do i = 1, nI
             iCoord_D(Dim1_) = iLower_D(Dim1_) + i - 1

             Jac = Jacobian_CI(i,j,1,Stencil6_)/nStencilExtra
             iStencilFirst = iStencilExtra_C(i,j,1)
             do iStencil = iStencilFirst, iStencilFirst + nStencilExtra - 1
                call HYPRE_SStructMatrixSetValues(i8A, iPart, iCoord_D, &
                     iVar, 1, (/iStencil/), (/Jac/), iError)
             end do
             iStencilExtra_C(i,j,1) = iStencil

             if(DoDebug)write(*,*)'-K iCoord, Jac =', iCoord_D, Jac
             
          end do
       end do

    end if

    ! check +K neighbor for resolution change
    DiLevel = DiLevelNei_IIIB(0,0,+1,iBlock)
    if(abs(DiLevel) == 1)then

       ! Number of connections per cell for coarser and finer neighbor
       nStencilExtra = 1
       if(DiLevel == -1) nStencilExtra = iRatio*jRatio

       if(DoDebug)write(*,*)'+K iProc, iPart, DiLevel=',iProc, iPart, DiLevel

       iCoord_D = iUpper_D
       do j = 1, nJ
          iCoord_D(Dim2_) = iUpper_D(Dim2_) + j - nJ
          do i = 1, nI
             iCoord_D(Dim1_) = iUpper_D(Dim1_) + i - nI

             Jac = Jacobian_CI(i,j,nK,Stencil7_)/nStencilExtra
             iStencilFirst = iStencilExtra_C(i,j,nK)
             do iStencil = iStencilFirst, iStencilFirst + nStencilExtra - 1
                call HYPRE_SStructMatrixSetValues(i8A, iPart, &
                     iCoord_D, iVar, 1, (/iStencil/), (/Jac/), iError)
             end do
             iStencilExtra_C(i,j,nK) = iStencil

             if(DoDebug)write(*,*)'+K iCoord, Jac =', iCoord_D, Jac

          end do
       end do

    end if

  end subroutine hypre_set_matrix_block

  !============================================================================
  subroutine hypre_set_matrix

    use ModMain, ONLY: test_string

    integer:: iError

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'hypre_set_matrix'
    !-------------------------------------------------------------------------
    ! Assemble matrix
    call set_oktest(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*) NameSub,' starting with DoInitHypreAmg=', &
         DoInitHypreAmg

    call HYPRE_SStructMatrixAssemble(i8A, iError)

    if(index(test_string,'HYPRE_PRINT_MATRIX') > 0)then
       call HYPRE_SStructMatrixPrint("matrix.dat", i8A, 0, iError)
       call stop_mpi('debug')
    end if

    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructMatrixAssemble done'

    if(.not.DoInitHypreAmg) RETURN

    ! Pass matrix to the solvers
    call HYPRE_SStructMatrixGetObject(i8A, i8ParA, iError)
    if(DoTestMe) write(*,*) NameSub,' HYPRE_SStructMatrixGetObject done'

    ! Setup AMG preconditioner for Krylov solver
    call timing_start('BoomerAMGSetup')
    call HYPRE_BoomerAMGSetup(i8Precond, i8ParA, i8ParB, i8ParX, iError)
    call timing_stop('BoomerAMGSetup')

    DoInitHypreAmg = .false.

    if(DoTestMe)write(*,*) NameSub,' finished'

  end subroutine hypre_set_matrix

  !============================================================================

  subroutine hypre_preconditioner(n, y_I)

    integer, intent(in):: n
    real, intent(inout):: y_I(n)

    integer:: i, iImplBlock, iPart, iError, iBlock, iNode
    real, allocatable:: Zero_I(:)

    logical, parameter:: DoDebug = .false.

    character(len=*), parameter:: NameSub = 'hypre_preconditioner'
    !-------------------------------------------------------------------------

    ! DoDebug = iProc == 1
    
    if(DoDebug)write(*,*) NameSub,' starting n, sum(y_I**)=', &
         n, sum(y_I**2)

    ! Preconditioning: y'= AMG.y

    ! Set initial guess value to zero
    allocate(Zero_I(nIJK))
    Zero_I = 0.0

    ! Set y_I as the RHS
    i = 1
    do iImplBlock = 1, nImplBlock
       iBlock   = impl2iblk(iImplBlock)
       iNode    = iNode_B(iBlock)
       iPart    = iTree_IA(Level_,iNode)
       iLower_D = 1 + (iTree_IA(Coord1_:CoordLast_,iNode)-1)*nCell_D
       iUpper_D = iLower_D - 1 + nCell_D

       call HYPRE_SStructVectorSetBoxValues(i8B, iPart, iLower_D, iUpper_D, &
            iVar, y_I(i), iError)

       call HYPRE_SStructVectorSetBoxValues(i8X, iPart, iLower_D, iUpper_D, &
            iVar, Zero_I, iError)

       if(DoDebug)write(*,*)'iImplBlock, iNode, iLower, iUpper, sum(y**2)=', &
            iImplBlock, iNode, iLower_D, iUpper_D, sum(y_I(i:i+nIJK-1)**2)

       i = i + nIJK
    end do
    deallocate(Zero_I)

    if(DoDebug)write(*,*) NameSub,' set X=0'

    ! Assemble vectors
    call HYPRE_SStructVectorAssemble(i8X, iError)
    call HYPRE_SStructVectorAssemble(i8B, iError)

    if(DoDebug)write(*,*) NameSub,' HYPRE_SStructVectorAssemble done'

    ! Pass the vectors to the solvers
    call HYPRE_SStructVectorGetObject(i8B, i8ParB, iError)
    call HYPRE_SStructVectorGetObject(i8X, i8ParX, iError)

    if(DoDebug)write(*,*) NameSub,' passed vectors to AMG'

    call HYPRE_BoomerAMGSolve(i8Precond, i8ParA, i8ParB, i8ParX, iError)

    if(DoDebug)write(*,*) NameSub,' applied AMG preconditioner'

    ! Get back solution
    call HYPRE_SStructVectorGather(i8x, iError);

    i = 1
    do iImplBlock = 1, nImplBlock
       iBlock   = impl2iblk(iImplBlock)
       iNode    = iNode_B(iBlock)
       iPart    = iTree_IA(Level_,iNode)
       iLower_D = 1 + (iTree_IA(Coord1_:CoordLast_,iNode)-1)*nCell_D
       iUpper_D = iLower_D - 1 + nCell_D

       call HYPRE_SStructVectorGetBoxValues(i8X, iPart, &
            iLower_D, iUpper_D, iVar, y_I(i), iError)

       if(DoDebug)write(*,*)'iImplBlock, iNode, iLower, iUpper, sum(y**2)=', &
            iImplBlock, iNode, iLower_D, iUpper_D, sum(y_I(i:i+nIJK-1)**2)

       i = i + nIJK
    end do

    if(DoDebug)write(*,*) NameSub,' finished'

  end subroutine hypre_preconditioner

  !================================================================================
  subroutine hypre_clean

    integer:: iError
    !-----------------------------------------------------------------------------

    call HYPRE_SStructGridDestroy(i8Grid, iError)
    call HYPRE_SStructStencilDestroy(i8Stencil, iError)
    call HYPRE_SStructGraphDestroy(i8Graph, iError)
    call HYPRE_SStructMatrixDestroy(i8A, iError)
    call HYPRE_SStructVectorDestroy(i8B, iError)
    call HYPRE_SStructVectorDestroy(i8X, iError)
    call HYPRE_BoomerAMGDestroy(i8Precond, iError)

  end subroutine hypre_clean
  !================================================================================
  subroutine hypre_finalize

    if(iLastDecomposition == -1) RETURN

    call hypre_clean
    deallocate(Value_I)
    iLastDecomposition = -1

  end subroutine hypre_finalize

end module ModImplHypre

