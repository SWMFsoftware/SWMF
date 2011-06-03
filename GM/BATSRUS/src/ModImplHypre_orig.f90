module ModImplHypre

  use ModKind,   ONLY: Int8_
  use BATL_size, ONLY: nDim, nI, nJ, nK, nIJK, nIJK_D
  use BATL_lib,  ONLY: nNode, nNodeUsed, iNode_B, iMortonNode_A, &
       iNodeNei_IIIB, DiLevelNei_IIIB, Unset_, iComm
  use ModImplicit, ONLY: impl2iBlk, MAT, nImplBlock => nImplBlk, &
       nStencil, Stencil1_, Stencil2_, Stencil3_, Stencil4_, Stencil5_, &
       Stencil6_, Stencil7_

  implicit none

  SAVE

  private
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

  integer, allocatable:: iPart_I(:), iPart_A(:)
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

contains

  !==========================================================================
  subroutine hypre_initialize

    integer:: iImplBlock, iBlock, iError, iPart, jPart, iPrintLevel

    logical:: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'hypre_initialize'
    !-------------------------------------------------------------------------
    if(allocated(iPart_I)) RETURN !!! More conditions needed here !!!

    call set_oktest(NameSub,DoTest,DoTestMe)

    if(DoTestMe)write(*,*) NameSub,' starting'

    ! Set the part index array indexed by implicit block index
    allocate(iPart_I(nImplBlock), iPart_A(nNode), Value_I(nStencil*nIJK))

    ! For now we assume that all used blocks are implicit...
    ! Otherwise we could loop through iTypeBlock_A or something...
    ! Note that parts are indexed from 0, not from 1 !

    nPart = nNodeUsed
    iPart_A = iMortonNode_A(1:nNode) - 1
    do iImplBlock = 1, nImplBlock
       iPart_I(iImplBlock) = iPart_A(iNode_B(impl2iBlk(iImplBlock)))
    end do

    ! part = local domain with local index space
    iLower_D = 1
    iUpper_D = nIJK_D(1:nDim)

    if(DoTestMe)write(*,*)'nPart, iLower_D, iUpper_D=',&
         nPart, iLower_D, iUpper_D

    ! Create an empty 3D grid object
    call HYPRE_SStructGridCreate(iComm, nDim, nPart, i8Grid, iError)

    if(DoTestMe)write(*,*)'HYPRE_SStructGridCreate done'

    ! Add local parts
    do iImplBlock = 1, nImplBlock
       iPart = iPart_I(iImplBlock)
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

    ! Setup connection between parts
    do iImplBlock = 1, nImplBlock
       iBlock = impl2iBlk(iImplBlock)
       iPart  = iPart_I(iImplBlock)

       ! -I neighbor
       if(DiLevelNei_IIIB(-1,0,0,iBlock) == 0)then
          jPart = iPart_A(iNodeNei_IIIB(0,1,1,iBlock))
          if(DoTestMe) write(*,*)'-I iPart, jPart=',iPart, jPart

          iLowerBc_D = 1;              iLowerBc_D(1) = 0
          iUpperBc_D = nIJK_D(1:nDim); iUpperBc_D(1) = 0
          jLowerBc_D = 1;              jLowerBc_D(1) = nI
          jUpperBc_D = nIJK_D(1:nDim)

          call HYPRE_SStructGridSetNeighborPart( i8Grid, &
               iPart, iLowerBc_D, iUpperBc_D, &
               jPart, jLowerBc_D, jUpperBc_D, &
               iIndexMap_D, iIndexDir_D, iError)
       end if

       ! +I neighbor
       if(DiLevelNei_IIIB(+1,0,0,iBlock) == 0)then
          jPart = iPart_A(iNodeNei_IIIB(3,1,1,iBlock))
          if(DoTestMe) write(*,*)'+I iPart, jPart=', iPart, jPart

          iLowerBc_D(1) = nI+1
          iUpperBc_D(1) = nI+1
          jLowerBc_D(1) = 1
          jUpperBc_D(1) = 1

          call HYPRE_SStructGridSetNeighborPart( i8Grid, &
               iPart, iLowerBc_D, iUpperBc_D, &
               jPart, jLowerBc_D, jUpperBc_D, &
               iIndexMap_D, iIndexDir_D, iError)
       end if

       if(nJ > 1)then
          ! -J neighbor
          if(DiLevelNei_IIIB(0,-1,0,iBlock) == 0)then
             jPart = iPart_A(iNodeNei_IIIB(1,0,1,iBlock))
             if(DoTestMe) write(*,*)'-J iPart, jPart=', iPart, jPart

             iLowerBc_D = 1;              iLowerBc_D(Dim2_) = 0
             iUpperBc_D = nIJK_D(1:nDim); iUpperBc_D(Dim2_) = 0
             jLowerBc_D = 1;              jLowerBc_D(Dim2_) = nJ
             jUpperBc_D = nIJK_D(1:nDim)

             call HYPRE_SStructGridSetNeighborPart( i8Grid, &
                  iPart, iLowerBc_D, iUpperBc_D, &
                  jPart, jLowerBc_D, jUpperBc_D, &
                  iIndexMap_D, iIndexDir_D, iError)
          end if

          ! +J neighbor
          if(DiLevelNei_IIIB(0,+1,0,iBlock) == 0)then
             jPart = iPart_A(iNodeNei_IIIB(1,3,1,iBlock))
             if(DoTestMe) write(*,*)'+J iPart, jPart=', iPart, jPart

             iLowerBc_D(Dim2_) = nJ+1
             iUpperBc_D(Dim2_) = nJ+1
             jLowerBc_D(Dim2_) = 1
             jUpperBc_D(Dim2_) = 1
             
             call HYPRE_SStructGridSetNeighborPart( i8Grid, &
                  iPart, iLowerBc_D, iUpperBc_D, &
                  jPart, jLowerBc_D, jUpperBc_D, &
                  iIndexMap_D, iIndexDir_D, iError)
          end if

       end if

       if(nK > 1)then
          ! -K neighbor
          if(DiLevelNei_IIIB(0,0,-1,iBlock) == 0)then
             jPart = iPart_A(iNodeNei_IIIB(1,1,0,iBlock))
             if(DoTestMe) write(*,*)'-K iPart, jPart=', iPart, jPart

             iLowerBc_D = 1;              iLowerBc_D(Dim3_) = 0
             iUpperBc_D = nIJK_D(1:nDim); iUpperBc_D(Dim3_) = 0
             jLowerBc_D = 1;              jLowerBc_D(Dim3_) = nK
             jUpperBc_D = nIJK_D(1:nDim)

             call HYPRE_SStructGridSetNeighborPart( i8Grid, &
                  iPart, iLowerBc_D, iUpperBc_D, &
                  jPart, jLowerBc_D, jUpperBc_D, &
                  iIndexMap_D, iIndexDir_D, iError)
          end if

          ! +K neighbor
          if(DiLevelNei_IIIB(0,0,+1,iBlock) == 0)then
             jPart = iPart_A(iNodeNei_IIIB(1,1,3,iBlock))
             if(DoTestMe) write(*,*)'+K iPart, jPart=', iPart, jPart

             iLowerBc_D(Dim3_) = nK+1
             iUpperBc_D(Dim3_) = nK+1
             jLowerBc_D(Dim3_) = 1
             jUpperBc_D(Dim3_) = 1
             
             call HYPRE_SStructGridSetNeighborPart( i8Grid, &
                  iPart, iLowerBc_D, iUpperBc_D, &
                  jPart, jLowerBc_D, jUpperBc_D, &
                  iIndexMap_D, iIndexDir_D, iError)
          end if

       end if

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

    ! Tell the graph which stencil to use for each variable on each part 
    do jPart = 0, nPart - 1
       call HYPRE_SStructGraphSetStencil(i8Graph, jPart, iVar, i8Stencil, &
            iError)
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
    call HYPRE_BoomerAMGSetMaxIter(i8Precond, 1, iError)
    call HYPRE_BoomerAMGSetTol(i8Precond, 0.0, iError)

    ! Print AMG solution info
    iPrintLevel = 0
    if(DoTest)iPrintLevel =2
    call HYPRE_BoomerAMGSetPrintLevel(i8Precond, iPrintLevel, iError)
    call HYPRE_BoomerAMGSetCoarsenType(i8Precond, 6, iError)

    ! Sym G.S./Jacobi hybrid
    call HYPRE_BoomerAMGSetRelaxType(i8Precond, 6, iError)
    call HYPRE_BoomerAMGSetNumSweeps(i8Precond, 1, iError)
    if(DoTestMe)write(*,*) NameSub,' HYPRE_BoomerAMGSetNumSweeps done'

    if(DoTestMe)write(*,*) NameSub,' finished'

  end subroutine hypre_initialize

  !===========================================================================

  subroutine hypre_set_matrix_block(iImplBlock, Jacobian_CI)

    ! Set the maxtrix elements corresponding to iImplBlock

    use ModMain, ONLY: TypeBc_I

    integer, intent(in):: iImplBlock
    real,    intent(inout):: Jacobian_CI(nI,nJ,nK,nStencil)

    integer:: iValue, i, j, k, iPart, iBlock, iError

    !------------------------------------------------------------------------
    iPart = iPart_I(iImplBlock)
    iBlock = impl2iblk(iImplBlock)

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

    call HYPRE_SStructMatrixSetBoxValues(i8A, iPart, iLower_D, iUpper_D, &
         iVar, nStencil, iStencil_I, Value_I, iError)

  end subroutine hypre_set_matrix_block

  !============================================================================
  subroutine hypre_set_matrix

    integer:: iError

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'hypre_setup_amg'
    !-------------------------------------------------------------------------
    ! Assemble matrix
    call HYPRE_SStructMatrixAssemble(i8A, iError)

    if(.not.DoInitHypreAmg) RETURN

    ! Pass matrix to the solvers
    call HYPRE_SStructMatrixGetObject(i8A, i8ParA, iError)
    if(DoTestMe)write(*,*) NameSub,' HYPRE_SStructMatrixGetObject done'

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

    integer:: i, iImplBlock, iPart, iError
    real, allocatable:: Value_I(:)

    logical, parameter:: DoDebug = .false.

    character(len=*), parameter:: NameSub = 'hypre_preconditioner'
    !-------------------------------------------------------------------------

    if(DoDebug)write(*,*) NameSub,' starting n, maxval, minval, y_I(1)=', &
         n, maxval(y_I), minval(y_I), y_I(1)

    ! Preconditioning: y'= AMG.y

    ! Set initial guess value to zero
    allocate(Value_I(nIJK))
    Value_I = 0.0

    ! Set y_I as the RHS
    i = 1
    do iImplBlock = 1, nImplBlock
       iPart = iPart_I(iImplBlock)

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
       iPart = iPart_I(iImplBlock)

       call HYPRE_SStructVectorGetBoxValues(i8X, iPart, &
            iLower_D, iUpper_D, iVar, y_I(i), iError)

       i = i + nIJK
    end do

    if(DoDebug)write(*,*) NameSub,' finished'

  end subroutine hypre_preconditioner

end module ModImplHypre

