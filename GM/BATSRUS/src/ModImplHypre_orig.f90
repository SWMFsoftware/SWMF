module ModImplHypre

  use ModKind,   ONLY: Int8_
  use BATL_size, ONLY: nDim, nI, nJ, nK, nIJK, nIJK_D
  use BATL_lib,  ONLY: iProc, nNode, nNodeUsed, iNode_B, iMortonNode_A, &
       iNodeNei_IIIB, DiLevelNei_IIIB, Unset_, iComm, &
       nRoot_D, MaxCoord_I, IsPeriodic_D, iTree_IA, &
       Proc_, Coord0_, Coord1_, Level_
  use ModImplicit, ONLY: impl2iBlk, MAT, nImplBlock => nImplBlk, &
       nStencil, Stencil1_, Stencil2_, Stencil3_, Stencil4_, Stencil5_, &
       Stencil6_, Stencil7_

  implicit none

  SAVE

  private ! except

  public:: hypre_read_param
  public:: hypre_initialize
  public:: hypre_set_matrix_block
  public:: hypre_set_matrix
  public:: hypre_preconditioner

  logical, public, parameter:: IsHypreAvailable = .true. ! true here
  logical, public:: DoInitHypreAmg = .true. ! set to false to reinitialize AMG

  ! local variables

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

    integer:: nLevel = 30
    integer:: iLevel, DiLevel, iNode, jNode, iCoord_D(nDim), jCoord_D(nDim)

    integer:: iImplBlock, iBlock, iError, iPart, jPart, iPrintLevel

    logical:: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'hypre_initialize'
    !-------------------------------------------------------------------------
    if(allocated(Value_I)) RETURN !!! More conditions needed here !!!

    call set_oktest(NameSub,DoTest,DoTestMe)

    if(DoTestMe)write(*,*) NameSub,' starting'

    allocate(Value_I(nStencil*nIJK))

    ! One part per level. We can possibly use MaxLevel (<=30), or count
    ! the number of current levels (1-5 or so).

    ! Level indexes go from 0 to nLevel, hence the +1
    nPart = nLevel + 1

    ! Create an empty 3D grid object
    call HYPRE_SStructGridCreate(iComm, nDim, nPart, i8Grid, iError)

    if(DoTestMe)write(*,*)'HYPRE_SStructGridCreate done'

    ! Add each block as a local box in the corresponding part (level)
    do iImplBlock = 1, nImplBlock
       iBlock   = impl2iBlk(iImplBlock)
       iNode    = iNode_B(iBlock)
       iPart    = iTree_IA(Level_,iNode)
       iLower_D = 1 + (iTree_IA(Coord1_:Coord0_+nDim,iNode)-1)*nIjk_D(1:nDim)
       iUpper_D = iLower_D - 1 + nIjk_D(1:nDim)

       if(DoTestMe)write(*,*)'iPart, iNode, iLower_D, iUpper_D=',&
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
          iUpper_D = MaxCoord_I(iLevel)*nRoot_D(1:nDim)*nIjk_D(1:nDim)

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

       ! -I neighbor
       DiLevel = DiLevelNei_IIIB(-1,0,0,iBlock)
       if(abs(DiLevel) == 1)then
          ! Neighbor node index and level
          jNode = iNodeNei_IIIB(0,1,1,iBlock)

          if(iProc == iTree_IA(Proc_,iNode) &
               .or. iProc == iTree_IA(Proc_,jNode))then

             jPart = iTree_IA(Level_,jNode)

             iCoord_D = (iTree_IA(Coord1_:Coord0_+nDim,iNode) - 1) &
                  *nIjk_D(1:nDim)+1
             jCoord_D = iTree_IA(Coord1_:Coord0_+nDim,jNode)*nIjk_D(1:nDim)

             if(DoTestMe)then
                write(*,*)'-I iPart, jPart =', iPart, jPart
                write(*,*)'-I iNode, jNode =', iNode, jNode
                write(*,*)'-I iCoord,jCoord=', iCoord_D, jCoord_D
             end if

             call HYPRE_SStructGraphAddEntries(i8Graph, iPart, &
                  iCoord_D, iVar, jPart, jCoord_D, iVar, iError)
          end if
       end if

       ! +I neighbor
       DiLevel = DiLevelNei_IIIB(+1,0,0,iBlock)
       if(abs(DiLevel) == 1)then
          ! Neighbor node index and level
          jNode = iNodeNei_IIIB(3,1,1,iBlock)

          if(iProc == iTree_IA(Proc_,iNode) &
               .or. iProc == iTree_IA(Proc_,jNode))then

             jPart = iTree_IA(Level_,jNode)

             iCoord_D = iTree_IA(Coord1_:Coord0_+nDim,iNode)*nIjk_D(1:nDim)
             jCoord_D = (iTree_IA(Coord1_:Coord0_+nDim,jNode)-1) &
                  *nIjk_D(1:nDim)+1

             if(DoTestMe)then
                write(*,*)'+I iPart, jPart =', iPart, jPart
                write(*,*)'+I iNode, jNode =', iNode, jNode
                write(*,*)'+I iCoord,jCoord=', iCoord_D, jCoord_D
             end if

             call HYPRE_SStructGraphAddEntries(i8Graph, iPart, &
                  iCoord_D, iVar, jPart, jCoord_D, iVar, iError)
          end if
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

  end subroutine hypre_initialize

  !===========================================================================

  subroutine hypre_set_matrix_block(iImplBlock, Jacobian_CI)

    ! Set the maxtrix elements corresponding to iImplBlock

    use ModMain, ONLY: TypeBc_I, ProcTest

    integer, intent(in):: iImplBlock
    real,    intent(inout):: Jacobian_CI(nI,nJ,nK,nStencil)

    real   :: Jac
    integer:: iValue, i, j, k, iPart, iBlock, iError
    integer:: DiLevel

    integer:: iNode

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
    iLower_D = 1 + (iTree_IA(Coord1_:Coord0_+nDim,iNode)-1)*nIjk_D(1:nDim)
    iUpper_D = iLower_D - 1 + nIjk_D(1:nDim)

    ! write(*,*)'!!! iPart, iProc, iBlock, iNode, iLower_D, iUpper_D=',&
    !     iPart, iProc, iBlock, iNode, iLower_D, iUpper_D, &
    !     maxval(Value_I), minval(Value_I)

    call HYPRE_SStructMatrixSetBoxValues(i8A, iPart, iLower_D, iUpper_D, &
         iVar, nStencil, iStencil_I, Value_I, iError)

    ! check -I neighbor for resolution change
    DiLevel = DiLevelNei_IIIB(-1,0,0,iBlock)
    if( abs(DiLevel) == 1)then

       Jac = Jacobian_CI(1,1,1,2)

       call HYPRE_SStructMatrixSetValues(i8A, iPart, &
            iLower_D, iVar, 1, (/nStencil/), (/Jac/), iError)

       if(DoDebug)write(*,*)'-I iProc, iPart, DiLevel, iCoord_D, Jac =', &
            iProc, iPart, DiLevel, iLower_D, Jac
    end if

    ! check +I neighbor for resolution change
    DiLevel = DiLevelNei_IIIB(+1,0,0,iBlock)
    if(abs(DiLevel) == 1)then

       Jac = Jacobian_CI(nI,1,1,3)

       call HYPRE_SStructMatrixSetValues(i8A, iPart, &
            iUpper_D, iVar, 1, (/nStencil/), (/Jac/), iError)

       if(DoDebug)write(*,*)'+I iProc, iPart, DiLevel, iCoord_D, Jac =', &
            iProc, iPart, DiLevel, iUpper_D, Jac
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
    real, allocatable:: Value_I(:)

    logical, parameter:: DoDebug = .false.

    character(len=*), parameter:: NameSub = 'hypre_preconditioner'
    !-------------------------------------------------------------------------

    ! DoDebug = iProc == 1
    
    if(DoDebug)write(*,*) NameSub,' starting n, maxval, minval, y_I(1)=', &
         n, maxval(y_I), minval(y_I), y_I(1)

    ! Preconditioning: y'= AMG.y

    ! Set initial guess value to zero
    allocate(Value_I(nIJK))
    Value_I = 0.0

    ! Set y_I as the RHS
    i = 1
    do iImplBlock = 1, nImplBlock
       iBlock   = impl2iblk(iImplBlock)
       iNode    = iNode_B(iBlock)
       iPart    = iTree_IA(Level_,iNode)
       iLower_D = 1 + (iTree_IA(Coord1_:Coord0_+nDim,iNode)-1)*nIjk_D(1:nDim)
       iUpper_D = iLower_D - 1 + nIjk_D(1:nDim)

       call HYPRE_SStructVectorSetBoxValues(i8B, iPart, iLower_D, iUpper_D, &
            iVar, y_I(i), iError)

       call HYPRE_SStructVectorSetBoxValues(i8X, iPart, iLower_D, iUpper_D, &
            iVar, Value_I, iError)

       i = i + nIJK
    end do
    deallocate(Value_I)

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
       iNode    = iNode_B(iImplBlock)
       iPart    = iTree_IA(Level_,iNode)
       iLower_D = 1 + (iTree_IA(Coord1_:Coord0_+nDim,iNode)-1)*nIjk_D(1:nDim)
       iUpper_D = iLower_D - 1 + nIjk_D(1:nDim)

       call HYPRE_SStructVectorGetBoxValues(i8X, iPart, &
            iLower_D, iUpper_D, iVar, y_I(i), iError)

       i = i + nIJK
    end do

    if(DoDebug)write(*,*) NameSub,' finished'

  end subroutine hypre_preconditioner

end module ModImplHypre

