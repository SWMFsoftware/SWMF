module ModHdf5

  use ModIO, ONLY: NamePlotDir, MaxFile, nplotvarmax
  use ModMpi
  use hdf5
  implicit none

  private ! except
  public hdf5_setup
  public write_var_hdf5
  public write_plot_hdf5

  SAVE

  character (len=80) :: filename
  integer :: plotFileNumber, checkpointFileNumber

  logical:: hdf5_writeRestartFile
!  integer:: hdf5_plotArrayBegin = 1, hdf5_plotArrayEnd
  integer:: hdf5_currentBlockIdx = 1, localNumBlocks

  ! Constants used by the FLASH hdf5 IO package

  ! Variables for dimensional indexing
  integer :: IAXIS = 1, JAXIS = 2, KAXIS = 3
  integer :: LOW = 1, HIGH = 2

  ! Refinement tree information arrays
  integer, allocatable:: which_child(:), lrefine(:), nodetype(:), bflags(:,:)
  real, allocatable:: bnd_box(:,:,:), bsize(:,:), coord(:,:)

  ! FLASH-formatted unkowns array
  real, allocatable :: unk(:,:,:,:,:)

  ! FLASH-formatted grid variables
  integer, allocatable :: gr_gid(:,:)
  integer :: gr_globalOffset, gr_globalNumBlocks

  ! FLASH IO variables
  integer:: io_outputSplitNum
  integer:: io_splitNumBlks

  character (len=20) :: io_flashRelease
  character (len=80) :: io_buildDate, io_buildDir, io_buildMachine, &
       io_fileCreationTime, io_setupTimeStamp,     &
       io_buildTimeStamp
  character (len=400) :: io_setupCall, io_cflags, io_fflags

  integer:: io_numRealParms, io_numIntParms, io_numStrParms, io_numLogParms
  integer:: io_numRealScalars, io_numIntScalars
  integer:: io_numStrScalars, io_numLogScalars

  integer:: io_intParmValues(1), io_intScalarValues(11)
  integer:: io_logToIntParmValues(1), io_logToIntScalarValues(1)
  logical:: io_logParmValues(1), io_logScalarValues(1)
  real   :: io_realParmValues(1), io_realScalarValues(1)

  character (len=80):: io_intParmNames(1)
  character (len=80):: io_intScalarNames(11)    
  character (len=80):: io_realParmNames(1)
  character (len=80):: io_realScalarNames(1)
  character (len=80):: io_strParmNames(1)
  character (len=80):: io_strScalarNames(1)
  character (len=80):: io_strParmValues(1)
  character (len=80):: io_strScalarValues(1)
  character (len=80):: io_logParmNames(1)
  character (len=80):: io_logScalarNames(1)
  character (len=80):: io_logParmNamesPrev(1)

  character (len=5) :: io_plotVarStr(nplotvarmax)
  character (len=5) :: io_unklabels(nplotvarmax)


  character(len=80):: io_geometry

  integer:: io_nPlotVars
  integer:: io_comm

  integer:: NUNK_VARS

contains

  !=====================================================================
  ! write_plot_hdf5
  !
  !   This subroutine is called by write_plot_common and handles setting
  !   up all FLASH adapter modules and calling the appropriate
  !   subroutines to write an hdf5 plotfile in the FLASH format.
  !=====================================================================
  subroutine write_plot_hdf5(filename, plotVarNames, nPlotVar, ifile)

    use ModProcMH, ONLY : iProc, nProc
    use ModIO, ONLY : nplotvarmax
    use ModMain, ONLY : nBlockMax
    use ModSize, ONLY: nI, nJ, nK
    use hdf5
    character (len=80), intent(in)  :: filename
    character (len=10), intent(in)  :: plotVarNames(nplotvarmax)
    integer, intent(in)  :: nPlotVar, ifile

    integer  :: fileID, error

    !--------------------------------------------------------------------------

    if (iProc == 0) write (*,*) '===================================='
    if (iProc == 0) write (*,*) 'Writing hdf5 parallel plotfile'
    if (iProc == 0) write (*,*) '------------------------------------'

    ! Setup this module to write a plotfile
    hdf5_writeRestartFile = .false.
    call hdf5_setupGrid(iProc)
    call hdf5_setupIOVars(plotVarNames, nPlotVar)
    call hdf5_setupTree()

!    hdf5_plotArrayEnd = nPlotVar

    if (iProc == 0) write (*,*) '  opening hdf5 file'
    !initialize the Fortran Hdf5 library
    call h5open_f(error)
    ! Open the hdf5 file and write data to it

!---------------------------------------------------------------------    
    fileID = -1
    call H5InitFile(fileID, filename)!, io_comm)
       if(fileID == -1) then
       if (iProc == 0) write (*,*)  "Error: unable to initialize file"
       call stop_mpi("unable to initialize hdf5 file")
    end if

   if (iProc == 0) write (*,*) '  writing data'
    call hdf5_writeData(iProc, nProc, fileID)

    if (iProc == 0) write (*,*) '  closing file'
    call h5garbage_collect_f(error)
    call h5fclose_f(fileID,error)
    !closing the hdf5 interface
    call h5close_f(error)
     if (error == -1) write (*,*) 'h5fclose_f failed'
    !Free MPI info too?

   
    if (iProc == 0) write (*,*) '===================================='

    ! Deallocate the variable array
    deallocate(unk)

  end subroutine write_plot_hdf5

  !=====================================================================

  subroutine write_var_hdf5(PlotVar, nPlotVar, iBLK)

    !  Handles adding one block of variable data that for plotting
    !  in the next plotfile.

    use ModMain, ONLY : nI, nJ, nK, nBlockMax
    use ModIO, ONLY: nPlotVarMax
    use hdf5
    use BATL_lib, ONLY : MaxBlock
    use ModProcMH, ONLY : iProc
    integer, intent(in) :: nPlotVar, iBLK
    real, intent(in) :: PlotVar(-1:nI+2,-1:nJ+2,-1:nK+2,nPlotVarMax)

    integer:: ii, jj, kk, iPlot
    !----------------------------------------------------------------------

    ! Allocate the storage array if it has not been allocated for
    ! this plotfile yet
    if (hdf5_currentBlockIdx > MaxBlock) &
        hdf5_currentBlockIdx = 1

    if (.not. allocated(unk)) then
       allocate(unk(nPlotVar, nI, nJ, nK, MaxBlock))
    end if
               
    do iPlot = 1, nPlotVar
       do ii = 1, nI; do jj = 1, nJ; do kk = 1, nK
          unk(iPlot, ii, jj, kk, hdf5_currentBlockIdx) = &
               PlotVar(ii, jj, kk, iPlot)


       end do; end do; end do
    end do

    hdf5_currentBlockIdx = hdf5_currentBlockIdx + 1

  end subroutine write_var_hdf5

  !=====================================================================

  subroutine hdf5_setup

    !   This subroutine is called by BATS_setup to allocate the dynamic
    !   arrays for various modules used by the hdf5 plotting module.

    use BATL_lib, ONLY : nDim, MaxBlock
    !---------------------------------------------------------------------  

    ! Allocate tree data arrays
    allocate(bnd_box(2, 3, MaxBlock))
    allocate(bsize(3, MaxBlock))
    allocate(coord(3, MaxBlock))
    allocate(which_child(MaxBlock))
    allocate(lrefine(MaxBlock))
    allocate(nodetype(MaxBlock))
    allocate(bflags(1, MaxBlock))

    ! Allocate grid data array
    allocate(gr_gid(2*nDim+1+2**nDim, MaxBlock))

  end subroutine hdf5_setup

  !=====================================================================

  subroutine hdf5_writeData(myPE, numProcs, fileID) 

    ! Handles gathering any additional information needed for plotting
    ! and then calling the appropriate C functions to produce an hdf5
    ! output file.
    use ModSize, only : MaxBlock, nI, nJ, nK
    use ModAdvance, only : State_VGB
    use hdf5
    !---------------------------------------------------------------------  

    integer, intent(in) :: myPE, fileID
    integer, intent(in) :: numProcs

    character (len=80) :: buff
    character(len=16) :: numToStr  

    integer :: blockID, sizeGid(2)
    integer :: i, lb, j, AllocateStatus

    integer, allocatable :: procnumber(:) 

    ! allocate storage to hold the 3,coordinate information and bounding box
    ! information
    real, allocatable :: globalVarMin(:), globalVarMax(:)
    real, allocatable :: globalFaceXMin(:), globalFaceXMax(:)
    real, allocatable :: globalFaceYMin(:), globalFaceYMax(:)
    real, allocatable :: globalFaceZMin(:), globalFaceZMax(:)

    real, allocatable :: unkBuf(:,:,:,:)


    logical :: writeGuardCells = .true.

    integer :: blkLimits(2,3), blkLimitsGC(2, 3)

    !Adjust our offset so that we can eliminate redundant data in a group
    integer :: localOffset, localRank, splitOffset
    integer :: error

    call  h5open_f(error)

    ! If we are writing a restart file, fill in the unkown variable array
    !  here
    if (hdf5_writeRestartFile) then
       unk(:,:,:,:,:) = State_VGB(1:NUNK_VARS, 1:nI, 1:nJ, 1:nK, :) !Replaced NUNK_VARS with nVar because
        !State_VGB is only allocated to nVar
    end if


    if(localNumBlocks == 0) then
       call stop_mpi("io_writeData: localNumBlocks == 0")
    end if

    !Find our local offset for split-file IO
    if (io_outputSplitNum > 1) then

       call MPI_ALLREDUCE(gr_globalOffset, splitOffset, 1, MPI_INTEGER, &
            MPI_MIN, io_comm, error)
       localOffset = gr_globalOffset - splitOffset
       !find number of blocks for a file
       call MPI_ALLREDUCE(localNumBlocks, io_splitNumBlks, 1, MPI_INTEGER,&
            MPI_SUM, io_comm, error)
    else
       localOffset = gr_globalOffset
       io_splitNumBlks = gr_globalNumBlocks
    end if

! 
!  Creates a compound datatype to hold things like file format version, Crash version, File creation
!  time, build date, etc...
!     !If we are writing a checkpoint file write all of the unk vars.  
!     !If writing a plotfile then only certain variables are written



       ! write the header info
   call writeHdf5Header(MyPE, io_nPlotVars, fileID, &  !removed io_geometry
            io_plotVarStr, io_setupCall, io_fileCreationTime, io_flashRelease,&
            io_buildDate, io_buildDir, io_buildMachine, io_cflags, io_fflags, &
            io_setupTimeStamp, io_buildTimeStamp)

 !Writes the int, scalar, string, runtime parameters.
    ! write the runtime parameters
    call writeHdf5RuntimeParams(mype, &
         fileID, &
         io_numRealParms, &
         io_realParmNames, &
         io_realParmValues, &
         io_numIntParms, &
         io_intParmNames, &
         io_intParmValues, &
         io_numStrParms, &
         io_strParmNames, &
         io_strParmValues, &
         io_numLogParms, &
         io_logParmNames, &
         io_logToIntParmValues)
         
!
    ! write the scalars
    call writeHdf5Scalars(mype, &
         fileID, &
         io_numRealScalars, &
         io_realScalarNames, &
         io_realScalarValues, &
         io_numIntScalars, &
         io_intScalarNames, &
         io_intScalarValues, &
         io_numStrScalars, &
         io_strScalarNames, &
         io_strScalarValues, &
         io_numLogScalars, &
         io_logScalarNames, &
         io_logToIntScalarValues)

    ! write the nodetype
    call writeHdf5Rank1Int(myPE, &
         fileID, &
         nodetype, &
         localNumBlocks, &
         io_splitNumBlks, &
         localOffset, "node type")

    
    ! write lrefine
    call writeHdf5Rank1Int(myPE, &
         fileID, &
         lrefine, &
         localNumBlocks, &
         io_splitNumBlks, &
         localOffset, "refine level")


    !data struct new to PM3.  which_child definitely
    !needed for restart. Adding bflags for completeness.
    ! write which child
    call writeHdf5Rank1Int(myPE, &
         fileID, &
         which_child, &
         localNumBlocks, &
         io_splitNumBlks, &
         localOffset, "which child")

    ! Write bflags
    call writeHdf5Rank1Int(myPE, &
         fileID, &
         bflags(1,:), &
         localNumBlocks, &
         io_splitNumBlks, &
         localOffset, "bflags")
    
    sizeGid = shape(gr_gid)
    ! Write the global ID
    call writeHdf5Rank2Int(myPE, &
         fileID, &
         gr_gid, &
         localNumBlocks, &
         io_splitNumBlks, &
         localOffset, "gid", sizeGid(1))


    !write the processor number. 
    allocate(procnumber(localnumblocks))
    do lb = 1, localnumblocks
       procnumber(lb) = mype
    end do

    call writehdf5rank1int(mype, &
         fileid,  & 
         procnumber,  & 
         localnumblocks,  & 
         io_splitnumblks,  & 
         localoffset, "processor number")

    deallocate(procnumber)

    ! write the bnd_box
    call writeHdf5Rank3Real(myPE, &
        fileID, bnd_box, &   !See that io_comm is also taken care of
        localNumBlocks, &
        io_splitNumBlks, &
        localOffset, "bounding box", 2, 3)
    
    !write the block center coordinates
    !(not to be confused with the cell coordinates)
   call writeHdf5Rank2Real(myPE, &
        fileID,  & 
        coord,  & 
        localNumBlocks,  &   
        io_splitNumBlks,  & 
        localOffset, "coordinates", 3)

    !write the physical size of each block
    call writeHdf5Rank2Real(myPE, &
        fileID,  & 
        bsize,  & 
        localNumBlocks,  & 
        io_splitNumBlks,  & 
        localOffset, "block size", 3)
     
     
     allocate(globalVarMin(NUNK_VARS)) !Replaced NUNK_VARS w/nVar
     allocate(globalVarMax(NUNK_VARS)) !Replaced NUNK_VARS w/nVar

    !get the max and minimum variables
     call hdf5_getVarExtrema(NUNK_VARS, globalVarMin, globalVarMax, 2) !Replaced NUNK_VARS w/nVar

    !store the unknowns
        allocate(unkBuf(nI, nJ, nK, MaxBlock))
     do i = 1, NUNK_VARS 


       unkBuf(1:nI,1:nJ,1:nK,1:MaxBlock) = &
            unk(i,1:nI, 1:nJ, 1:nK,1:MaxBlock) 

      call writeHdf5Unknowns(myPE, &
           fileID, & 
           nI, & 
           nJ, & 
           nK, & 
           unkBuf, & 
           globalVarMin(i), &
           globalVarMax(i), &
           io_unklabels(i), &
           localNumBlocks, &
           io_splitNumBlks,  & 
           localOffset)

    end do

       deallocate(unkBuf)
    deallocate(globalVarMin)
    deallocate(globalVarMax)
    

    ! Removed segment for plotting scratch variables

    ! Removed segment for plotting face variables


    return

  end subroutine hdf5_writeData

  !=====================================================================
  ! hdf5_initFile, hdf5_closeFile
  !
  !    These subroutines handle calling the appropriate Fortran functions
  !    for opening/closing hdf5 files.
  !=====================================================================
!     use hdf5
! 
!     character (len=80), intent(in) :: filename
! !---------------------------------------------------------------------    
!     fileID = -1
!     call ioH5InitFile(fileID, filename, io_comm, io_outputSplitNum)
!        if(fileID == -1) then
!        if (iProc == 0) write (*,*)  "Error: unable to initialize file"
!        call stop_mpi("unable to initialize hdf5 file")
!     end if
! 
!   end subroutine hdf5_initFile
! 
  !=====================================================================

subroutine hdf5_setupTree

    ! Initializes FLASH tree arrays for plotting

    use BATL_tree, ONLY : nChild, nNodeUsed, iTree_IA, iNode_B, Level_, &
         Parent_, Child0_
    use BATL_lib,    ONLY : CoordMin_DB, CoordMax_DB, CellSize_DB
    use ModMain,     ONLY : nBlockMax, nBlock, unusedBLK
    use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK
    use ModSize  
    integer  :: iBLK, iNode, iParent, iChild, idx
    integer :: halfI = 1, halfJ = 1, halfK = 1
    !---------------------------------------------------------------------
    idx = 1

    if (nI > 1) halfI = nI/2
    if (nJ > 1) halfJ = nJ/2
    if (nK > 1) halfK = nK/2

    do iBLK = 1, nBlockMax
       if (unusedBLK(iBLK)) CYCLE

       ! The coord array stores cell-centered coordinates        
       coord(IAXIS, idx) = x_BLK(halfI, halfJ, halfK, iBLK)
       coord(JAXIS, idx) = y_BLK(halfI, halfJ, halfK, iBLK)
       coord(KAXIS, idx) = z_BLK(halfI, halfJ, halfK, iBLK)


       bnd_box(LOW,:,idx) = CoordMin_DB(:,iBLK)
       bnd_box(HIGH,:,idx) = CoordMax_DB(:,iBLK)
       bsize(:,idx) = CellSize_DB(:,iBLK)

       ! nodetype holds no info for now
       nodetype(idx) = 2

       iNode = iNode_B(idx)

       ! In CRASH, refinement indices start at 0 as opposed to 1 with FLASH. 
       ! Add 1 to the value to comply with the FLASH standard.
       lrefine(idx) = iTree_IA(Level_, iNode) + 1

       ! Find which child this node is
       iParent = iTree_IA(Parent_, iNode)
       which_child(:) = -1
       if (iParent >= 0) then
          do iChild = 1, nChild
             if (iTree_IA(Child0_+iChild, iParent) == iNode) &
                  which_child(idx) = iChild
          end do
       end if

       bflags(1, idx) = -1

       idx = idx+1

    end do

  end subroutine hdf5_setupTree

  !=====================================================================
  ! hdf5_setIOVars
  !
  !    Initializes some of the IO variables imported from FLASH. Some are 
  !    initialized as constants here rather than in their declaration
  !    because their implementation may change.
  !=====================================================================
  subroutine hdf5_setupIOVars(plotVarNames, nPlotVar)

    use BATL_lib, ONLY : nDim, nDimAmr, TypeGeometry, nI, nJ, nK
    use BATL_tree, ONLY : nInfo => ChildLast_, nNode
    use ModMain, ONLY : nBlockAll, Time_Simulation
    use ModVarIndexes, ONLY : NameVar_V, nVar, nFluid

    ! Arguments
    character (len=10), intent(in), optional :: plotVarNames(nplotvarmax)
    integer, intent(in), optional :: nPlotVar

    integer  :: iVar, iLen, labelLeng
    !--------------------------------------------------------------------------
    io_flashRelease = 'flash release'

    io_geometry = 'io geometry'

    ! Set compiler and setup string info as empty
    io_buildDate = 'build date'
    io_buildDir = 'build dir'
    io_buildMachine = 'build machine'
    io_fileCreationTime = 'file creation time'
    io_setupTimeStamp = 'setup time stamp'
    io_buildTimeStamp = 'build time stamp'
    io_setupCall = 'setup call'
    io_cflags ='ccflags'
    io_fflags = 'fflags'

    io_outputSplitNum = 1 ! file splitting not supported

    ! FLASH makes a distinction between plot var labels and unk var labels, but
    ! CRASH will treat both the same.  The Visit plugin needs 4 character variale names
    ! with null termination so we format the variable names here.  Underscores are used 
    ! of spaces because the Visit plugin does not recognise the variable if its name includes
    ! spaces.

    
    if (.not. hdf5_writeRestartFile) then
       do iVar = 1, nPlotVar

          io_plotVarStr(iVar) =  trim(plotVarNames(iVar)(1:4))


            labelLeng = len_trim(io_plotVarStr(iVar))
            if (labelLeng .ne. 4) then
                !        label(labelLeng:4) = 'P'
                do iLen = labelLeng + 1, 4 
                    io_plotVarStr(iVar)(iLen:) = "_"
                end do
            end if
            io_plotVarStr(iVar)(5:) = CHAR(0)




          io_unklabels(iVar) = io_plotVarStr(iVar)
       end do
       io_nPlotVars = nPlotVar
    else
       io_nPlotVars = nVar + nFluid
       do iVar = 1, io_nPlotVars
          io_plotVarStr(iVar) = trim(NameVar_V(iVar))
 
            labelLeng = len_trim(io_plotVarStr(iVar))
            if (labelLeng .ne. 4) then
                !        label(labelLeng:4) = 'P'
                do iLen = labelLeng + 1, 4 
                    io_plotVarStr(iVar)(iLen:) = "_"
                end do
            end if
            io_plotVarStr(iVar)(5:) = CHAR(0)

         
          io_unklabels(iVar) = io_plotVarStr(iVar)
       end do
    end if

    ! Output of runtime parameters, scalars, and logfiles 
    ! is not implemented yet, so clear the corresponding variables
    io_numRealParms = 0
    io_numIntParms = 0
    io_numStrParms = 0
    io_numLogParms = 0
    io_numRealScalars = 1
    io_numIntScalars = 11
    io_numStrScalars = 0
    io_numLogScalars = 0
    ! probably don't need to alloc string and value arrays...

    ! The FLASH Visit plugin requires some scalar metadata. 
    ! This should be moved to the empty io_prepareLists subroutine eventually.
    io_intScalarValues(1) = nI
    io_intScalarNames(1) = 'nxb' 
    io_intScalarValues(2) = nJ
    io_intScalarNames(2) = 'nyb'
    io_intScalarValues(3) = nK
    io_intScalarNames(3) = 'nzb'
    io_intScalarValues(4) = nDim
    io_intScalarNames(4) = 'dimensionality'
    io_intScalarValues(5) = nBlockAll
    io_intScalarNames(5) = 'globalnumblocks'
    io_intScalarValues(6) = io_outputSplitNum-1
    io_intScalarNames(6) = 'splitnumblocks'
    io_intScalarValues(7) = 0
    io_intScalarNames(7) = 'splitnumparticles'

    ! Restart files require node and tree info
    io_intScalarValues(8) = nNode
    io_intScalarNames(8) = 'nNode'
    io_intScalarValues(9) = nInfo
    io_intScalarNames(9) = 'nInfo'
    io_intScalarValues(10) = nDim
    io_intScalarNames(10) = 'nDim'
    io_intScalarValues(11) = nDimAmr
    io_intScalarNames(11) = 'nDimAmr'


    io_realScalarValues(1) = Time_Simulation
    io_realScalarNames(1) = 'time'

    ! io_chkGuardCells = .false.

    ! we need more MPI initialization for split file IO
    io_comm = MPI_COMM_WORLD 

    NUNK_VARS = nVar !+ nFluid

  end subroutine hdf5_setupIOVars

  !=====================================================================

  subroutine hdf5_setupGrid(myPE)

    ! Initializes the grid data arrays derived from the FLASH code.

    use BATL_tree,  ONLY : nChild, nNodeUsed, iTree_IA, iNodeNei_IIIB, &
         Child0_, Block_, ChildLast_, Child1_,       &
         Parent_, Proc_, Coord1_, CoordLast_
    use ModMain,    ONLY : nBlockMax, nBlockAll, unusedBLK, nBlockAll, &
         global_block_number
    use BATL_lib,   ONLY : iNode_B
    use ModProcMH,  ONLY : iProc, nProc, iComm
    use ModSize,    ONLY : nI, nJ, nK
    use BATL_size,  ONLY : nDim
    use ModParallel,ONLY : iProc_A

    integer, intent(IN)  :: myPE

    integer :: i, j, k, iNode, iBLK, iGid, error, idx, nUnused
    integer :: bn, pe
    integer :: localIdx(nBlockMax)
    integer, allocatable :: blocksPerProc(:), offsetPerProc(:)
    integer, allocatable :: localIdx_pe(:,:)
    !---------------------------------------------------------------------  

    gr_globalNumBlocks = nBlockAll

    allocate(blocksPerProc(nProc))
    allocate(offsetPerProc(0:nProc-1))
    allocate(localIdx_pe(0:nProc-1, nBlockMax))

    ! Get the number of blocks local to this processor and the mapping
    ! from block index to output index (where unused blocks are skipped)
    call hdf5_getLocalBlks(localNumBlocks, localIdx)

    ! Collect results across all processors
    call MPI_Allgather(localNumBlocks, 1, MPI_INTEGER, blocksPerProc, &
         1, MPI_INTEGER, iComm, error)
    call MPI_Allgather(localIdx, nBlockMax, MPI_INTEGER, localIdx_pe, &
         nBlockMax, MPI_INTEGER, iComm, error)

    offsetPerProc(0) = 0
    if (nProc > 1)&
        offsetPerProc(1) = blocksPerProc(1)
    if (nProc > 2) then
        do i = 2, nProc-1
           offsetPerProc(i) = offsetPerProc(i-2) + blocksPerProc(i-2)
        end do
    end if
    gr_globalOffset = offsetPerProc(iProc)

    deallocate(blocksPerProc)
    gr_gid(:,:) = -32

    ! keep track of the number of unused blocks we've passed
    nUnused = 0

    idx = 1
    do iBLK = 1, nBlockMax
       if (unusedBLK(iBLK)) then
          nUnused = nUnused + 1
          CYCLE
       end if

       iNode = iNode_B(iBLK)

       ! Begin by filling in the neighbor cells for this cell

       if (iNodeNei_IIIB(0, 1, 1, iBLK) > 0) then
          bn = iTree_IA(Block_, iNodeNei_IIIB(0, 1, 1, iBLK))
          pe = iTree_IA(Proc_, iNodeNei_IIIB(0, 1, 1, iBLK))
          gr_gid(1, idx) = localIdx_pe(pe, bn) + offsetPerProc(pe)
       end if

       if (iNodeNei_IIIB(3, 1, 1, iBLK) > 0) then
          bn = iTree_IA(Block_, iNodeNei_IIIB(3, 1, 1, iBLK))
          pe = iTree_IA(Proc_, iNodeNei_IIIB(3, 1, 1, iBLK))
          gr_gid(2, idx) = localIdx_pe(pe, bn) + offsetPerProc(pe)
       end if

       iGid = 3

       if (nJ > 1) then
          if (iNodeNei_IIIB(1, 0, 1, iBLK) > 0) then
             bn = iTree_IA(Block_, iNodeNei_IIIB(1, 0, 1, iBLK))
             pe = iTree_IA(Proc_, iNodeNei_IIIB(1, 0, 1, iBLK))
             gr_gid(3, idx) = localIdx_pe(pe, bn) + offsetPerProc(pe)
          end if

          if (iNodeNei_IIIB(1, 3, 1, iBLK) > 0) then
             bn = iTree_IA(Block_, iNodeNei_IIIB(1, 3, 1, iBLK))
             pe = iTree_IA(Proc_, iNodeNei_IIIB(1, 3, 1, iBLK))
             gr_gid(4, idx) = localIdx_pe(pe, bn) + offsetPerProc(pe)
          end if

          iGid = 5
       end if

       if (nK > 1) then
          if (iNodeNei_IIIB(1, 1, 0, iBLK) > 0) then
             bn = iTree_IA(Block_, iNodeNei_IIIB(1, 1, 0, iBLK))
             pe = iTree_IA(Proc_, iNodeNei_IIIB(1, 1, 0, iBLK))
             gr_gid(5, idx) = localIdx_pe(pe, bn) + offsetPerProc(pe)
          end if

          if (iNodeNei_IIIB(1, 1, 3, iBLK) > 0) then
             bn = iTree_IA(Block_, iNodeNei_IIIB(1, 1, 3, iBLK))
             pe = iTree_IA(Proc_, iNodeNei_IIIB(1, 1, 3, iBLK))
             gr_gid(6, idx) = localIdx_pe(pe, bn) + offsetPerProc(pe)
          end if

          iGid = 7
       end if

       ! Store parent's node number and children's node numbers
       do i = 0, 2**nDim
          ! No AMR tree info stored yet...
          !if (iTree_IA(Child0_+i, iNode) > 0) then
          !  gr_gid(iGid, idx) = iTree_IA(Block_, iTree_IA(Child0_+i, iNode)) &
          !    + offsetPerProc(iTree_IA(Proc_, iTree_IA(Child0_+i, iNode)))
          !else
          gr_gid(iGid, idx) = -1
          !end if

          iGid = iGid+1
       end do

       idx = idx + 1
    end do

    deallocate(offsetPerProc)
    deallocate(localIdx_pe)

  end subroutine hdf5_setupGrid

  !=====================================================================

  subroutine hdf5_getLocalBlks(numBlocks, localIdx)

    ! Returns the number of blocks used on the local processor as well
    ! as an array mapping the index of each used block to its output
    ! array index (where unused blocks are skipped). 

    use ModMain, ONLY : nBlockMax, unusedBLK, nBlockALL
    use ModProcMH,  ONLY : iProc

    integer, intent(out) :: numBlocks
    integer, dimension(nBlockMax), intent(out), optional :: localIdx
    integer  :: i, nBlockUsed
    !---------------------------------------------------------------------

    nBlockUsed = 0
    localIdx(:) = -1
    do i = 1, nBlockMax
       if (unusedBLK(i)) CYCLE

       nBlockUsed = nBlockUsed+1
       localIdx(i) = nBlockUsed
    end do

    numBlocks = nBlockUsed
    return
  end subroutine hdf5_getLocalBlks

  !=====================================================================

  subroutine hdf5_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, &
       gridDataStruct)

    ! Returns the index limits for a block, both with and without guard cells.

    use ModSize

    integer,intent(IN) :: blockId
    integer, dimension(2,3), intent(OUT) :: blkLimits,blkLimitsGC
    integer, optional, intent(IN) :: gridDataStruct

    !---------------------------------------------------------------------

    blkLimits(LOW,:) = 1
    blkLimitsGC(LOW, :) = 1-gcn

    blkLimits(HIGH, IAXIS) = nI
    blkLimits(HIGH, JAXIS) = nJ
    blkLimits(HIGH, KAXIS) = nK

    blkLimitsGC(HIGH, IAXIS) = nI+gcn
    blkLimitsGC(HIGH, JAXIS) = nJ+gcn
    blkLimitsGC(HIGH, KAXIS) = nK+gcn

  end subroutine hdf5_getBlkIndexLimits

  !=====================================================================

  subroutine hdf5_getVarExtrema( &
       nvars, globalVarMin, globalVarMax, gridDataStruct)
    
    ! Returns the global variable extrema for each unkown.

    use ModMain
    use ModAdvance, ONLY : State_VGB
    use ModSize

    integer, intent(in) :: nvars, gridDataStruct
    real, DIMENSION(nvars), INTENT(out) :: globalVarMin, globalVarMax

    real, allocatable  :: varMin(:), varMax(:)
    integer  :: iBLK, i, j, k, iVar, error
    !--------------------------------------------------------------------------

    if (nvars > 0) then

       allocate(varMin(nvars))
       allocate(varMax(nvars))

       varMin(:) = huge(varMin)
       varMax(:) = -huge(varMax)

       do iBLK = 1, nBlockMax
          if (unusedBLK(iBLK)) cycle

          do i=1,nI; do j=1,nJ; do k=1,nK

             do iVar = 1, nvars
                varMin(iVar) = min(varMin(iVar), State_VGB(iVar,i,j,k,iBLK))
                varMax(iVar) = max(varMax(iVar), State_VGB(iVar,i,j,k,iBLK))
             end do

          end do; end do; end do;

       end do

       call MPI_AllReduce(varMin, globalVarMin, nvars, MPI_DOUBLE_PRECISION, &
            MPI_MIN, MPI_COMM_WORLD, error)

       call MPI_AllReduce(varMax, globalVarMax, nvars, MPI_DOUBLE_PRECISION, &
            MPI_MAX, MPI_COMM_WORLD, error)

       deallocate(varMin)
       deallocate(varMax)

    end if
    return

  end subroutine hdf5_getVarExtrema


!subroutines that call hdf5 libraries

subroutine H5InitFile(fileID, filename)

    use hdf5
    use ModProcMH, Only : iComm
    character (len=80), intent(in) :: filename
    integer :: error, accTemplate
    integer, intent(inout) :: fileID
    !integer, intent(in) :: io_comm
    INTEGER :: mpierror       ! MPI error flag
    INTEGER :: info,dxpList
    INTEGER :: mpi_size, mpi_rank

    !Initalize currentBlock andplotArrayBegin for this file

!    call h5open_f(error)                    

    info = MPI_INFO_NULL

    !create MPI info Object


    !Create file access propertty list
    call h5pcreate_f(H5P_FILE_ACCESS_F, accTemplate, error)
    CALL h5pset_fapl_mpio_f(accTemplate, iComm, info, error)
    ! Create the file collectively.
    
    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, fileID, error, access_prp = accTemplate)
    CALL h5pclose_f(accTemplate, error)
    if (error == -1) &
        call stop_mpi("error in subroutine H5InitFile. Error marker 1")
 
    if (error == -1) fileID = -1
end subroutine

 !==================================================================================
!Write Runtime Parameters

subroutine writeHdf5RuntimeParams(mype, &
         fileID, &
         io_numRealParams, &
         io_realParamNames, &
         io_realParamValues, &
         io_numIntParams, &
         io_intParamNames, &
         io_intParamValues, &
         io_numStrParams, &
         io_strParamNames, &
         io_strParamValues, &
         io_numLogParams, &
         io_logParamNames, &
         io_logToIntParamValues)

    use hdf5
    integer, intent(in) :: MyPE, fileID, io_numRealParams, io_numIntParams, io_numStrParams
    integer, intent(in) :: io_numLogParams
    character (len=80), intent(in) :: io_realParamNames(io_numRealParams)
    character (len=80), intent(in) :: io_intParamNames(io_numIntParams)
    character (len=80), intent(in) :: io_strParamNames(io_numStrParams)
    character (len=80), intent(in) :: io_logParamNames(io_numLogParams)
    character (len=80), intent(in) :: io_strParamValues(io_numStrParams)
    integer, intent(in) :: io_intParamValues(io_numIntParams), io_logToIntParamValues(io_numLogParams)
    real, intent(in) :: io_realParamValues(io_numRealParams)
    integer :: error, i 
    integer(HSIZE_T) :: nameSize, valueSize, compoundTypeSize, offset
    integer(HSIZE_T) :: dimens1D(1)
    integer(HID_T) :: dataspace, realListType, strListType, intListType, nameListType, valueListType
    integer(HID_T) :: logListType, dataset, nameType, plist_id
    character ::  datasetNames(4)
    datasetNames(1) = "real runtime parameters"
    datasetNames(2) = "integer runtime parameters"
    datasetNames(3) = "string runtime parameters"
    datasetNames(4) = "logical runtime parameters"

    call h5open_f(error)                    
    nameSize = 80 !byte size for names since names are given 80 characters
    if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 1")
     
    !datatype for the names
   
    !======Real runtime params===
       
    if (io_numRealParams > 0) then
        dimens1D(1) = io_numRealParams
        
        call h5tcopy_f(H5T_NATIVE_CHARACTER, nameType, error)
        nameSize = 80
        call h5tset_size_f(nameType, nameSize, error)

        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 2")
        
        
 
        !create the write buffer
     
      !get the byte size for the compound datatype
        call h5tget_size_f(H5T_NATIVE_DOUBLE, valueSize, error)
        call h5tget_size_f(nameType, nameSize, error)
        compoundTypeSize = (valueSize + nameSize) * dimens1D(1)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 3")
        
     
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        CALL h5pset_preserve_f(plist_id, .TRUE., error)
        
        !create a simple 1d dataspace
        call h5screate_simple_f(1, dimens1D, dataspace, error)

        call h5tcreate_f(H5T_COMPOUND_F, compoundTypeSize, realListType, error) 
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 4")
        
            
        offset = 0
        CALL h5tinsert_f(realListType, "name", offset, nameType, error)
        offset = nameSize
        CALL h5tinsert_f(realListType, "value", offset, H5T_NATIVE_DOUBLE, error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 5")
        
     

        !The dataspace for the compound type
        
        offset = 0
        call h5tcreate_f(H5T_COMPOUND_F, nameSize, nameListType, error)
        call h5tinsert_f(nameListType, "name",  offset, nameType, error)
        offset = 0
        call h5tcreate_f(H5T_COMPOUND_F, valueSize, valueListType, error)
        call h5tinsert_f(valueListType, "value",  offset, H5T_NATIVE_DOUBLE, error)

        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 6")
     
        call h5dcreate_f(fileID, datasetNames(1), realListType, dataspace, dataset, &
            error)
        
        if (myPE ==0) then
            call h5dwrite_f(dataset, nameListType, io_realParamNames, dimens1D, H5S_ALL_F,&
                     H5S_ALL_F, xfer_prp = plist_id)
            call h5dwrite_f(dataset, valueListType, io_realParamValues, dimens1D, H5S_ALL_F,&
                     H5S_ALL_F, xfer_prp = plist_id)
        if (error == -1) & 
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 7")
            
        end if
        call h5tclose_f(nameType, error)
        call h5tclose_f(nameListType, error)
        call h5tclose_f(valueListType, error)
        call h5tclose_f(realListType, error)
        call h5dclose_f(dataset, error)
        call h5sclose_f(dataspace, error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 8")
        
    end if

    !=====Integer runtime params======

    if (io_numIntParams > 0) then
        dimens1D(1) = io_numIntParams
        
        call h5tcopy_f(H5T_NATIVE_CHARACTER, nameType, error)
        nameSize = 80
        call h5tset_size_f(nameType, nameSize, error)

        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 9")

        
 
        !create the write buffer
     
      !get the byte size for the compound datatype
        call h5tget_size_f(H5T_NATIVE_INTEGER, valueSize, error)
        call h5tget_size_f(nameType, nameSize, error)
        compoundTypeSize = (valueSize + nameSize) * dimens1D(1)

        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        CALL h5pset_preserve_f(plist_id, .TRUE., error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 10")
            
        !create a simple 1d dataspace
        call h5screate_simple_f(1, dimens1D, dataspace, error)

        call h5tcreate_f(H5T_COMPOUND_F, compoundTypeSize, intListType, error) 
        
        offset = 0
        CALL h5tinsert_f(intListType, "name", offset, H5T_NATIVE_CHARACTER, error)
        offset = nameSize
        CALL h5tinsert_f(intListType, "value", offset, H5T_NATIVE_INTEGER, error)
        if (error == -1) & 
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 11")
     

        !The dataspace for the compound type
        
        offset = 0
        call h5tcreate_f(H5T_COMPOUND_F, nameSize, nameListType, error)
        call h5tinsert_f(nameListType, "name",  offset, nameType, error)
        offset = 0
        call h5tcreate_f(H5T_COMPOUND_F, valueSize, valueListType, error)
        call h5tinsert_f(valueListType, "value",  offset, H5T_NATIVE_INTEGER, error)


        call h5dcreate_f(fileID, datasetNames(1), intListType, dataspace, dataset, &
            error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 12")
            
        if (myPE ==0) then
            call h5dwrite_f(dataset, nameListType, io_intParamNames, dimens1D, H5S_ALL_F,&
                     H5S_ALL_F, xfer_prp = plist_id)
            call h5dwrite_f(dataset, valueListType, io_intParamValues, dimens1D, H5S_ALL_F,&
                     H5S_ALL_F, xfer_prp = plist_id)
            if (error == -1) &
                call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 13")
 
        end if
        call h5tclose_f(nameType, error)
        call h5tclose_f(nameListType, error)
        call h5tclose_f(valueListType, error)
        call h5tclose_f(intListType, error)
        call h5dclose_f(dataset, error)
        call h5sclose_f(dataspace, error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 14")
     
    end if

    !=========Logical Runtime Params============

    if (io_numLogParams > 0) then
        dimens1D(1) = io_numLogParams
        
        call h5tcopy_f(H5T_NATIVE_CHARACTER, nameType, error)
        nameSize = 80
        call h5tset_size_f(nameType, nameSize, error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 15")
     
    
 
        !create the write buffer
     
      !get the byte size for the compound datatype
        call h5tget_size_f(H5T_NATIVE_INTEGER, valueSize, error)
        call h5tget_size_f(nameType, nameSize, error)
        compoundTypeSize = (valueSize + nameSize) * dimens1D(1)

        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        CALL h5pset_preserve_f(plist_id, .TRUE., error)
        if (error == -1) & 
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 16")
            
        !create a simple 1d dataspace
        call h5screate_simple_f(1, dimens1D, dataspace, error)

        call h5tcreate_f(H5T_COMPOUND_F, compoundTypeSize, logListType, error) 
        
        offset = 0
        CALL h5tinsert_f(logListType, "name", offset, nameType, error)
        offset = nameSize
        CALL h5tinsert_f(logListType, "value", offset, H5T_NATIVE_INTEGER, error)
        if (error == -1) & 
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 17")
        
     

        !The dataspace for the compound type
        
        offset = 0
        call h5tcreate_f(H5T_COMPOUND_F, nameSize, nameListType, error)
        call h5tinsert_f(nameListType, "name",  offset, nameType, error)
        offset = 0
        call h5tcreate_f(H5T_COMPOUND_F, valueSize, valueListType, error)
        call h5tinsert_f(valueListType, "value",  offset, H5T_NATIVE_INTEGER, error)


        call h5dcreate_f(fileID, datasetNames(1), logListType, dataspace, dataset, &
            error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 18")
            
        if (myPE ==0) then
            call h5dwrite_f(dataset, nameListType, io_logParamNames, dimens1D, H5S_ALL_F,&
                     H5S_ALL_F, xfer_prp = plist_id)
            call h5dwrite_f(dataset, valueListType, io_logToIntParamValues, dimens1D, H5S_ALL_F,&
                     H5S_ALL_F, xfer_prp = plist_id)
            if (error == -1) &
                call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 19")
         
        end if
        call h5tclose_f(nameType, error)
        call h5tclose_f(nameListType, error)
        call h5tclose_f(valueListType, error)
        call h5tclose_f(logListType, error)
        call h5dclose_f(dataset, error)
        call h5sclose_f(dataspace, error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 20")
     
   end if

    !=========String Runtime Params===========

    if (io_numStrParams > 0) then
        dimens1D(1) = io_numStrParams
        
        call h5tcopy_f(H5T_NATIVE_CHARACTER, nameType, error)
        nameSize = 80
        call h5tset_size_f(nameType, nameSize, error)

        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 21")
        
 
        !create the write buffer
     
      !get the byte size for the compound datatype
        call h5tget_size_f(nameType, nameSize, error)
        valueSize = nameSize
        compoundTypeSize = (valueSize + nameSize) * dimens1D(1) !values of the same type as names

        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        CALL h5pset_preserve_f(plist_id, .TRUE., error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 22")
            
        !create a simple 1d dataspace
        call h5screate_simple_f(1, dimens1D, dataspace, error)

        call h5tcreate_f(H5T_COMPOUND_F, compoundTypeSize, strListType, error) 
        
        offset = 0
        CALL h5tinsert_f(strListType, "name", offset, nameType, error)
        offset = nameSize
        CALL h5tinsert_f(strListType, "value", offset, nameType, error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 23")
     

        !The dataspace for the compound type
        
        offset = 0
        call h5tcreate_f(H5T_COMPOUND_F, nameSize, nameListType, error)
        call h5tinsert_f(nameListType, "name",  offset, nameType, error)
        offset = 0
        call h5tcreate_f(H5T_COMPOUND_F, valueSize, valueListType, error)
        call h5tinsert_f(valueListType, "value",  offset, nameType, error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 24")
     

        call h5dcreate_f(fileID, datasetNames(1), strListType, dataspace, dataset, &
            error)
        
        if (myPE ==0) then
            call h5dwrite_f(dataset, nameListType, io_strParamNames, dimens1D, H5S_ALL_F,&
                     H5S_ALL_F, xfer_prp = plist_id)
            call h5dwrite_f(dataset, valueListType, io_strParamValues, dimens1D, H5S_ALL_F,&
                     H5S_ALL_F, xfer_prp = plist_id)
            if (error == -1) &
                call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 25")
     
        end if


        call h5tclose_f(nameType, error)
        call h5tclose_f(nameListType, error)
        call h5tclose_f(valueListType, error)
        call h5tclose_f(strListType, error)
        call h5dclose_f(dataset, error)
        call h5sclose_f(dataspace, error)
     if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5RuntimeParams. Error marker 26")
 
    end if

end subroutine

!=========================================================================================

        
subroutine writeHdf5Header(MyPE, io_nPlotVars, fileID, &  !removed io_geometry
            io_plotVarStr, io_setupCall, io_fileCreationTime, io_flashRelease,&
            io_buildDate, io_buildDir, io_buildMachine, io_cflags, io_fflags, &
            io_setupTimeStamp, io_buildTimeStamp)

    use hdf5
    integer :: FileFormatVersion, ivar, iLen, labelLeng
    integer :: error, MyPE, io_nPlotVars
    character (len=80), intent(in) :: io_fileCreationTime
    character (len=20), intent(in) :: io_flashRelease
    character (len=80), intent(in) :: io_buildDate
    character (len=80), intent(in) :: io_buildDir
    character (len=80), intent(in) :: io_buildMachine
    character (len=400), intent(in) :: io_cflags
    character (len=400), intent(in) :: io_fflags 
    character (len=400), intent(in) :: io_setupCall
    character (len=80), intent(in) :: io_setupTimeStamp
    character (len=80), intent(in) :: io_buildTimeStamp
    character (len=5), intent(in) :: io_plotVarStr(io_nPlotVars)
    integer(HID_T) ::  stringTypeFH, stringTypeTwenty, stringTypeEighty, stringTypeFour, FFVType
    integer(HID_T) :: SCType, FCTType, FRType,  BDType, BDirType, BMType, CFType, FFType 
    integer(HID_T) :: STSType, BTSType, FVType, stringTypeOne
    integer(HID_T) :: memspace,  dataspace, dataset, simInfoType, rank, plist_id
    integer(HID_T), intent(in) :: fileID
    integer(HSIZE_T) :: simInfoSize, dimens1D(1), dimens2D(2), intSize, sizeOne
    integer(HSIZE_T) :: sizeFH, sizeEighty, sizeFour, sizeTwenty, strOffset, intOffset, offset
    character (len=5) :: UnknownNameArray(1,io_nPlotVars)
   FileFormatVersion = 9

    call h5open_f(error)

    sizeFH = 400
    call h5tcopy_f(H5T_NATIVE_CHARACTER, stringTypeFH, error)
    call h5tset_size_f(stringTypeFH, sizeFH, error)
 !   call h5tset_pad_f(stringTypeFH, H5T_PAD_ZERO_F, H5T_PAD_ZERO_F, error)
    sizeTwenty = 20
    call h5tcopy_f(H5T_NATIVE_CHARACTER, stringTypeTwenty, error)
    call h5tset_size_f(stringTypeTwenty, sizeTwenty, error)
 !   call h5tset_pad_f(stringTypeTwenty, H5T_PAD_ZERO_F, H5T_PAD_ZERO_F, error)
    sizeEighty = 80
    call h5tcopy_f(H5T_NATIVE_CHARACTER, stringTypeEighty, error)
    call h5tset_size_f(stringTypeEighty, sizeEighty, error)
!    call h5tset_pad_f(stringTypeEighty, H5T_PAD_ZERO_F, H5T_PAD_ZERO_F, error)
    sizeOne = 1
    call h5tcopy_f(H5T_NATIVE_CHARACTER, stringTypeOne, error)
    call h5tset_size_f(stringTypeOne, sizeOne, error)
 
    if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5Header. Error marker 1")
 

    rank = 1
    dimens1D(1) = 1

    
    ! Create write buffer
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    CALL h5pset_preserve_f(plist_id, .TRUE., error)
    if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5Header. Error marker 2")
    
    call h5tget_size_f(H5T_NATIVE_INTEGER, intSize, error)
    call h5tget_size_f(stringTypeFH, sizeFH, error)
    call h5tget_size_f(stringTypeTwenty, sizeTwenty, error)
    call h5tget_size_f(stringTypeEighty, sizeEighty, error)
    call h5tget_size_f(stringTypeOne, sizeOne, error)

   ! simInfoSize = sizeFH*2 + sizeEighty*8 + sizeTwenty + intSize
    simInfoSize = 400*3 +80*8 + 20 + 4 
    if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5Header. Error marker 3")
 

    call h5tget_offset_f(stringTypeFH, strOffset, error)
    call h5tget_offset_f(H5T_NATIVE_INTEGER, intOffset, error)

    call h5screate_simple_f(rank, dimens1D, dataspace, error)
    call h5tcreate_f(H5T_COMPOUND_F, 2*simInfoSize, simInfoType, error)
    if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5Header. Error marker 4")
    
    offset = 0
    call h5tinsert_f(simInfoType, "file format version", offset , H5T_NATIVE_INTEGER,&
        error)
    offset = offset + intSize
    
    call h5tinsert_f(simInfoType, "setup call",  offset, stringTypeFH, error)
    offset = offset + sizeFH
    call h5tinsert_f(simInfoType, "file creation time", offset, stringTypeEighty, error)
    offset = offset + sizeEighty
    call h5tinsert_f(simInfoType, "flash release",  offset, stringTypeTwenty, error)
    offset = offset + sizeTwenty
    call h5tinsert_f(simInfoType, "build date",  offset,  stringTypeEighty, error)
    offset = offset + sizeEighty
    call h5tinsert_f(simInfoType, "build dir",  offset, stringTypeEighty, error)
    offset = offset + sizeEighty
    call h5tinsert_f(simInfoType, "build machine",  offset, stringTypeEighty, error)
    offset = offset + sizeEighty
    call h5tinsert_f(simInfoType, "cflags",  offset, stringTypeFH, error)
    offset = offset + sizeFH
    call h5tinsert_f(simInfoType, "fflags",  offset, stringTypeFH, error)
    offset = offset + sizeFH
    call h5tinsert_f(simInfoType, "setup time stamp", offset, stringTypeEighty, error)
    offset = offset + sizeEighty
    call h5tinsert_f(simInfoType, "build time stamp", offset, stringTypeEighty, error)
    if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5Header. Error marker 5")
 

    call h5dcreate_f(fileID, "sim info", simInfoType, dataspace, dataset, error, H5P_DEFAULT_F) 
    offset = 0
    if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5Header. Error marker 6")
 
!In Fortran each compound datatype 
    call h5tcreate_f(H5T_COMPOUND_F, intSize, FFVType, error)
    call h5tinsert_f(FFVType, "file format version", offset , H5T_NATIVE_INTEGER,&
        error)

    offset = 0
    call h5tcreate_f(H5T_COMPOUND_F, sizeFH, SCType, error)
    call h5tinsert_f(SCType, "setup call",  offset, stringTypeFH, error)

    offset = 0
    call h5tcreate_f(H5T_COMPOUND_F, sizeEighty, FCTType, error)
    call h5tinsert_f(FCTType, "file creation time", offset, stringTypeEighty, error)

    offset = 0
    call h5tcreate_f(H5T_COMPOUND_F, sizeTwenty, FRType, error)
    call h5tinsert_f(FRType, "flash release",  offset, stringTypeTwenty, error)

    offset = 0
    call h5tcreate_f(H5T_COMPOUND_F, sizeEighty, BDType, error)
    call h5tinsert_f(BDType, "build date",  offset,  stringTypeEighty, error)

    offset = 0
    call h5tcreate_f(H5T_COMPOUND_F, sizeEighty, BDirType, error)
    call h5tinsert_f(BDirType, "build dir",  offset, stringTypeEighty, error)

    offset = 0
    call h5tcreate_f(H5T_COMPOUND_F, sizeEighty, BMType, error)
    call h5tinsert_f(BMType, "build machine",  offset, stringTypeEighty, error)

    offset = 0
    call h5tcreate_f(H5T_COMPOUND_F, sizeFH, CFType, error)
    call h5tinsert_f(CFType, "cflags",  offset, stringTypeFH, error)

    offset = 0
    call h5tcreate_f(H5T_COMPOUND_F, sizeFH, FFType, error)
    call h5tinsert_f(FFType, "fflags",  offset, stringTypeFH, error)

    offset = 0
    call h5tcreate_f(H5T_COMPOUND_F, sizeEighty, STSType, error)
    call h5tinsert_f(STSType, "setup time stamp", offset, stringTypeEighty, error)

    offset = 0
    call h5tcreate_f(H5T_COMPOUND_F, sizeEighty, BTSType, error)
    call h5tinsert_f(BTSType, "build time stamp", offset, stringTypeEighty, error)
    if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5Header. Error marker 7")
 

    if (myPE == 0 ) then
        call h5dwrite_f(dataset, FFVType, FileFormatVersion, dimens1D, error, H5S_ALL_F,&
             H5S_ALL_F, xfer_prp = plist_id)
        call h5dwrite_f(dataset, SCType, io_setupCall, dimens1D, error, H5S_ALL_F,&
             H5S_ALL_F, xfer_prp = plist_id)
        call h5dwrite_f(dataset, FCTType, io_fileCreationTime, dimens1D, error, H5S_ALL_F,&
             H5S_ALL_F, xfer_prp = plist_id)
        call h5dwrite_f(dataset, FRType, io_flashRelease, dimens1D, error, H5S_ALL_F,&
             H5S_ALL_F, xfer_prp = plist_id)
        call h5dwrite_f(dataset, BDType, io_buildDate, dimens1D, error, H5S_ALL_F,&
             H5S_ALL_F, xfer_prp = plist_id)
        call h5dwrite_f(dataset, BDirType, io_buildDir, dimens1D, error, H5S_ALL_F,&
             H5S_ALL_F, xfer_prp = plist_id)
        call h5dwrite_f(dataset, BMType, io_buildMachine, dimens1D, error, H5S_ALL_F,&
             H5S_ALL_F, xfer_prp = plist_id)
        call h5dwrite_f(dataset, CFType, io_cflags, dimens1D, error, H5S_ALL_F,&
             H5S_ALL_F, xfer_prp = plist_id)
        call h5dwrite_f(dataset, FFType, io_fflags, dimens1D, error, H5S_ALL_F,&
             H5S_ALL_F, xfer_prp = plist_id)
        call h5dwrite_f(dataset, STSType, io_setupTimeStamp, dimens1D, error, H5S_ALL_F,&
             H5S_ALL_F, xfer_prp = plist_id)
        call h5dwrite_f(dataset, BTSType, io_buildTimeStamp, dimens1D, error, H5S_ALL_F,&
             H5S_ALL_F, xfer_prp = plist_id)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Header. Error marker 8")
     
    end if



    call h5tclose_f(FFVType ,error)
    call h5tclose_f(SCType ,error)
    call h5tclose_f(FRType ,error)
    call h5tclose_f(BDType ,error)
    call h5tclose_f(BDirType ,error)
    call h5tclose_f(BMType ,error)
    call h5tclose_f(CFType ,error)
    call h5tclose_f(FFType ,error)
    call h5tclose_f(STSType ,error)
    call h5tclose_f(BTSType ,error)
    call h5tclose_f(stringTypeFH, error)
    call h5tclose_f(stringTypeTwenty, error)
    call h5tclose_f(stringTypeEighty, error)
    call h5dclose_f(dataset, error)
    call h5tclose_f(simInfoType, error)
    call h5sclose_f(dataspace, error)
    if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5Header. Error marker 9")
 
   if (io_nPlotVars > 0) then
    rank = 2  
    dimens2D(1) = 1 
    dimens2d(2) = io_nPlotVars
    
    UnknownNameArray(1,1:io_nPlotVars) = io_plotVarStr(1:io_nPlotVars)

   sizeFour = 5
    call h5tcopy_f(H5T_NATIVE_CHARACTER, stringTypeFour, error)
    call h5tset_size_f(stringTypeFour, sizeFour, error)
    if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5Header. Error marker 10")
 
    call h5screate_simple_f(rank, dimens2D, dataspace, error)
    call h5dcreate_f(fileID, "unknown names", stringTypeFour, dataspace, dataset, error)
     if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5Header. Error marker 11")
 
    if (MyPE == 0) &
        call h5dwrite_f(dataset, stringTypeFour, UnknownNameArray, dimens2D, error, H5S_ALL_F,&
        H5S_ALL_F, H5P_DEFAULT_F)
     if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5Header. Error marker 11")


 
    call h5tclose_f(stringTypeFour, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
     if (error == -1)& 
        call stop_mpi("error in subroutine writeHdf5Header. Error marker 13")
 
  end if

end subroutine
!==================================================================================================
!Writes scalars

subroutine writeHdf5Scalars(mype, &
         fileID, &
         io_numRealScalars, &
         io_realScalarNames, &
         io_realScalarValues, &
         io_numIntScalars, &
         io_intScalarNames, &
         io_intScalarValues, &
         io_numStrScalars, &
         io_strScalarNames, &
         io_strScalarValues, &
         io_numLogScalars, &
         io_logScalarNames, &
         io_logToIntScalarValues)

    use hdf5
    integer, intent(in) :: MyPE, fileID, io_numRealScalars, io_numIntScalars, io_numStrScalars
    integer, intent(in) :: io_numLogScalars
    character (len=80), intent(in) :: io_realScalarNames(io_numRealScalars)
    character (len=80), intent(in) :: io_intScalarNames(io_numIntScalars)
    character (len=80), intent(in) :: io_strScalarNames(io_numStrScalars)
    character (len=80), intent(in) :: io_logScalarNames(io_numLogScalars)
    character (len=80), intent(in) :: io_strScalarValues(io_numStrScalars)
    integer, intent(in) :: io_intScalarValues(io_numIntScalars), io_logToIntScalarValues(io_numLogScalars)
    real, intent(in) :: io_realScalarValues(io_numRealScalars)
    integer :: error, i 
    integer(HSIZE_T) :: nameSize, valueSize, compoundTypeSize, offset
    integer(HSIZE_T) :: dimens1D(1)
    integer(HID_T) :: dataspace, realListType, strListType, intListType, nameListType, valueListType
    integer(HID_T) :: logListType, dataset, nameType, plist_id
    character (len=15)::  datasetNames(4)
    datasetNames(1) = "real scalars"
    datasetNames(2) = "integer scalars"
    datasetNames(3) = "string scalars"
    datasetNames(4) = "logical scalars"

    call h5open_f(error)                    
    nameSize = 80 !byte size for names since names are given 80 characters
     
    !datatype for the names
    if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 1")
   
    !======Real runtime params===
       
    if (io_numRealScalars > 0) then
        dimens1D(1) = io_numRealScalars
        
        call h5tcopy_f(H5T_NATIVE_CHARACTER, nameType, error)
        nameSize = 80
        call h5tset_size_f(nameType, nameSize, error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 2")
 
           
 
        !create the write buffer
     
      !get the byte size for the compound datatype
        call h5tget_size_f(H5T_NATIVE_DOUBLE, valueSize, error)
        call h5tget_size_f(nameType, nameSize, error)
        compoundTypeSize = valueSize + nameSize
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 3")
     
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        CALL h5pset_preserve_f(plist_id, .TRUE., error)
        
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 4")
     

               !create a simple 1d dataspace
        call h5screate_simple_f(1, dimens1D, dataspace, error)

        call h5tcreate_f(H5T_COMPOUND_F, compoundTypeSize, realListType, error) 
        
         if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 5")
     
        offset = 0
        CALL h5tinsert_f(realListType, "name", offset, nameType, error)
        offset = nameSize
        CALL h5tinsert_f(realListType, "value", offset, H5T_NATIVE_DOUBLE, error)


        !The dataspace for the compound type
        
        offset = 0
        call h5tcreate_f(H5T_COMPOUND_F, nameSize, nameListType, error)
        call h5tinsert_f(nameListType, "name",  offset, nameType, error)
        offset = 0
        call h5tcreate_f(H5T_COMPOUND_F, valueSize, valueListType, error)
        call h5tinsert_f(valueListType, "value",  offset, H5T_NATIVE_DOUBLE, error)

        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 6")
     
        call h5dcreate_f(fileID, datasetNames(1), realListType, dataspace, dataset, &
            error)
        
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 7")
     
        if (myPE ==0) then
            call h5dwrite_f(dataset, nameListType, io_realScalarNames, dimens1D, H5S_ALL_F,&
                     H5S_ALL_F, xfer_prp = plist_id)
            call h5dwrite_f(dataset, valueListType, io_realScalarValues, dimens1D, H5S_ALL_F,&
                H5S_ALL_F, xfer_prp = plist_id)
            if (error == -1) &
                call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 8")
         
        end if
        call h5tclose_f(nameType, error)
        call h5tclose_f(nameListType, error)
        call h5tclose_f(valueListType, error)
        call h5tclose_f(realListType, error)
        call h5dclose_f(dataset, error)
        call h5sclose_f(dataspace, error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 9")
     
    end if

    !=====Integer runtime params======

    if (io_numIntScalars > 0) then
        dimens1D(1) = io_numIntScalars
        
        call h5tcopy_f(H5T_NATIVE_CHARACTER, nameType, error)
        nameSize = 80
        call h5tset_size_f(nameType, nameSize, error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 2")
 
           
 
        !create the write buffer
     
      !get the byte size for the compound datatype
        call h5tget_size_f(H5T_NATIVE_INTEGER, valueSize, error)
        call h5tget_size_f(nameType, nameSize, error)
        compoundTypeSize = valueSize + nameSize
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 3")
     
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        CALL h5pset_preserve_f(plist_id, .TRUE., error)
        
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 4")
     

               !create a simple 1d dataspace
        call h5screate_simple_f(1, dimens1D, dataspace, error)

        call h5tcreate_f(H5T_COMPOUND_F, compoundTypeSize, intListType, error) 
        
         if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 5")
     
        offset = 0
        CALL h5tinsert_f(intListType, "name", offset, nameType, error)
        offset = nameSize
        CALL h5tinsert_f(intListType, "value", offset, H5T_NATIVE_INTEGER, error)


        !The dataspace for the compound type
        
        offset = 0
        call h5tcreate_f(H5T_COMPOUND_F, nameSize, nameListType, error)
        call h5tinsert_f(nameListType, "name",  offset, nameType, error)
        offset = 0
        call h5tcreate_f(H5T_COMPOUND_F, valueSize, valueListType, error)
        call h5tinsert_f(valueListType, "value",  offset, H5T_NATIVE_INTEGER, error)

        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 6")
     
        call h5dcreate_f(fileID, datasetNames(2), intListType, dataspace, dataset, &
            error)
        
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 7")
     
        if (myPE ==0) then
            call h5dwrite_f(dataset, nameListType, io_intScalarNames, dimens1D, H5S_ALL_F,&
                     H5S_ALL_F, xfer_prp = plist_id)
            call h5dwrite_f(dataset, valueListType, io_intScalarValues, dimens1D, H5S_ALL_F,&
                H5S_ALL_F, xfer_prp = plist_id)
            if (error == -1) &
                call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 8")
         
        end if
        call h5tclose_f(nameType, error)
        call h5tclose_f(nameListType, error)
        call h5tclose_f(valueListType, error)
        call h5tclose_f(intListType, error)
        call h5dclose_f(dataset, error)
        call h5sclose_f(dataspace, error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 9")
     
    end if

    !=========Logical Runtime Scalars============

    if (io_numLogScalars > 0) then
        dimens1D(1) = io_numLogScalars
        
        call h5tcopy_f(H5T_NATIVE_CHARACTER, nameType, error)
        nameSize = 80
        call h5tset_size_f(nameType, nameSize, error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 2")
 
           
 
        !create the write buffer
     
      !get the byte size for the compound datatype
        call h5tget_size_f(H5T_NATIVE_INTEGER, valueSize, error)
        call h5tget_size_f(nameType, nameSize, error)
        compoundTypeSize = valueSize + nameSize
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 3")
     
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        CALL h5pset_preserve_f(plist_id, .TRUE., error)
        
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 4")
     

               !create a simple 1d dataspace
        call h5screate_simple_f(1, dimens1D, dataspace, error)

        call h5tcreate_f(H5T_COMPOUND_F, compoundTypeSize, logListType, error) 
        
         if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 5")
     
        offset = 0
        CALL h5tinsert_f(logListType, "name", offset, nameType, error)
        offset = nameSize
        CALL h5tinsert_f(logListType, "value", offset, H5T_NATIVE_INTEGER, error)


        !The dataspace for the compound type
        
        offset = 0
        call h5tcreate_f(H5T_COMPOUND_F, nameSize, nameListType, error)
        call h5tinsert_f(nameListType, "name",  offset, nameType, error)
        offset = 0
        call h5tcreate_f(H5T_COMPOUND_F, valueSize, valueListType, error)
        call h5tinsert_f(valueListType, "value",  offset, H5T_NATIVE_INTEGER, error)

        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 6")
     
        call h5dcreate_f(fileID, datasetNames(3), logListType, dataspace, dataset, &
            error)
        
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 7")
     
        if (myPE ==0) then
            call h5dwrite_f(dataset, nameListType, io_logScalarNames, dimens1D, H5S_ALL_F,&
                     H5S_ALL_F, xfer_prp = plist_id)
            call h5dwrite_f(dataset, valueListType, io_logToIntScalarValues, dimens1D, H5S_ALL_F,&
                H5S_ALL_F, xfer_prp = plist_id)
            if (error == -1) &
                call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 8")
         
        end if
        call h5tclose_f(nameType, error)
        call h5tclose_f(nameListType, error)
        call h5tclose_f(valueListType, error)
        call h5tclose_f(logListType, error)
        call h5dclose_f(dataset, error)
        call h5sclose_f(dataspace, error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 9")
        end if

    !=========String Runtime Scalars===========

    if (io_numStrScalars > 0) then
        dimens1D(1) = io_numStrScalars
        
        call h5tcopy_f(H5T_NATIVE_CHARACTER, nameType, error)
        nameSize = 80
        call h5tset_size_f(nameType, nameSize, error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 2")
 
           
 
        !create the write buffer
     
      !get the byte size for the compound datatype
        call h5tget_size_f(nameType, valueSize, error)
        call h5tget_size_f(nameType, nameSize, error)
        compoundTypeSize = valueSize + nameSize
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 3")
     
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        CALL h5pset_preserve_f(plist_id, .TRUE., error)
        
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 4")
     

               !create a simple 1d dataspace
        call h5screate_simple_f(1, dimens1D, dataspace, error)

        call h5tcreate_f(H5T_COMPOUND_F, compoundTypeSize, strListType, error) 
        
         if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 5")
     
        offset = 0
        CALL h5tinsert_f(strListType, "name", offset, nameType, error)
        offset = nameSize
        CALL h5tinsert_f(strListType, "value", offset, nameType, error)


        !The dataspace for the compound type
        
        offset = 0
        call h5tcreate_f(H5T_COMPOUND_F, nameSize, nameListType, error)
        call h5tinsert_f(nameListType, "name",  offset, nameType, error)
        offset = 0
        call h5tcreate_f(H5T_COMPOUND_F, valueSize, valueListType, error)
        call h5tinsert_f(valueListType, "value",  offset, nameType, error)

        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 6")
     
        call h5dcreate_f(fileID, datasetNames(3), strListType, dataspace, dataset, &
            error)
        
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 7")
     
        if (myPE ==0) then
            call h5dwrite_f(dataset, nameListType, io_strScalarNames, dimens1D, H5S_ALL_F,&
                     H5S_ALL_F, xfer_prp = plist_id)
            call h5dwrite_f(dataset, valueListType, io_strScalarValues, dimens1D, H5S_ALL_F,&
                H5S_ALL_F, xfer_prp = plist_id)
            if (error == -1) &
                call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 8")
         
        end if
        call h5tclose_f(nameType, error)
        call h5tclose_f(nameListType, error)
        call h5tclose_f(valueListType, error)
        call h5tclose_f(strListType, error)
        call h5dclose_f(dataset, error)
        call h5sclose_f(dataspace, error)
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Scalars. Error marker 9")    
  
    end if

end subroutine




!======================================y===========================================================
!Writes rank 1 integer data
subroutine writeHdf5Rank1Int(myPE, fileID, dataBuff, localNumBlocks, io_splitNumBlks,&
    localOffset, description)

    use hdf5
    implicit none

    integer :: rank
    integer :: error
    integer(HID_T) :: dataset, plist_id
    integer(HSIZE_T) :: dimens1D(1)
    integer(HID_T) :: dataspace
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: start1d(1)
    integer(HSIZE_T) :: stride1d(1)
    integer(HSIZE_T) :: count1d(1)
    integer(HID_T), intent(in) :: fileID
    integer, intent(in) :: myPE 
    integer, intent(in) :: localNumBlocks
    integer, intent(in) :: io_splitNumBlks
    integer, intent(in) :: dataBuff(:)
    integer, intent(in) :: localOffset
    character (len=*), intent(in) :: description
    !Set the dimensions of the dataset
    rank = 1
    dimens1D(1) = io_splitNumBlks
    call h5open_f(error)
    if (error == -1) &
    
        call stop_mpi("error in subroutine writeHdf5Rank1Int. Error marker 1")
 

    call h5screate_simple_f(rank, dimens1D(1), dataspace, error) 
  !create the dataset
        call h5dcreate_f(fileID, description, H5T_NATIVE_INTEGER, dataspace, dataset, error)

         call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
         call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5Rank1Int. Error marker 2")

 
    start1D(1) = localOffset 
    stride1D(1) = 1
    count1D(1) = localNumBlocks
    
    if (localNumBlocks > 0 ) then    
        !create the hyperslab.  This will differ on the different processors
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start1D, count1D, error, &
            stride = stride1D)
!         !create the memory space
        if (error == -1)& 
            call stop_mpi("error in subroutine writeHdf5Rank1Int. Error marker 3")
     
        dimens1D(1) = localNumBlocks
        call h5screate_simple_f(rank, dimens1D, memspace, error)
         !Write the data
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Rank1Int. Error marker 4")
     
        call h5dwrite_f(dataset, H5T_NATIVE_INTEGER, dataBuff, dimens1D, error, &
            mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
         end if
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Rank1Int. Error marker 5")
     
    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
    if (error == -1)& 
        call stop_mpi("error in subroutine writeHdf5Rank1Int. Error marker 6")

end subroutine
!===========================================================================================
! Rank 2 Real

subroutine writeHdf5Rank2Real(myPE, fileID, dataBuff, localNumBlocks, io_splitNumBlks,&
    localOffset, description, dimens)

    use hdf5
    implicit none

    integer :: rank, dimens
    integer :: error
    integer(HID_T) :: dataset, plist_id
    integer(HSIZE_T) :: dimens2D(2)
    integer(HID_T) :: dataspace
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: start2D(2)
    integer(HSIZE_T) :: stride2D(2)
    integer(HSIZE_T) :: count2D(2)
    integer(HID_T), intent(in) :: fileID
    integer, intent(in) :: myPE 
    integer, intent(in) :: localNumBlocks
    integer, intent(in) :: io_splitNumBlks
    real, intent(in) :: dataBuff(:,:)
    integer, intent(in) :: localOffset
    character (len=*), intent(in) :: description
    !Set the dimensions of the dataset
    rank = 2
    dimens2D(1) = dimens
    dimens2D(2) = io_splitNumBlks
    call h5open_f(error)
    if (error == -1) &
    
        call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 1")
 

    call h5screate_simple_f(rank, dimens2D, dataspace, error) 
  !create the dataset
        call h5dcreate_f(fileID, description, H5T_NATIVE_DOUBLE, dataspace, dataset, error)

         call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
         call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 2")

    start2d(1) = 0
    start2D(2) = localOffset 
    stride2D(:) = 1
    count2D(1) = dimens2D(1)
    count2D(2) = localNumBlocks
    
    if (localNumBlocks > 0 ) then    
        !create the hyperslab.  This will differ on the different processors
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start2D, count2D, error, &
            stride = stride2D)
!         !create the memory space
        if (error == -1)& 
            call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 3")
     
        dimens2D(2) = localNumBlocks
        call h5screate_simple_f(rank, dimens2D, memspace, error)
         !Write the data
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 4")
     
        call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dataBuff, dimens2D, error, &
            mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
         end if
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 5")
     
    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
    if (error == -1)& 
        call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 6")

end subroutine

!============================================================================================
!
subroutine writeHdf5Rank2Int(myPE, fileID, dataBuff, localNumBlocks, io_splitNumBlks,&
    localOffset, description, dimens)

    use hdf5
    implicit none

    integer :: rank, dimens
    integer :: error
    integer(HID_T) :: dataset, plist_id
    integer(HSIZE_T) :: dimens2D(2)
    integer(HID_T) :: dataspace
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: start2D(2)
    integer(HSIZE_T) :: stride2D(2)
    integer(HSIZE_T) :: count2D(2)
    integer(HID_T), intent(in) :: fileID
    integer, intent(in) :: myPE 
    integer, intent(in) :: localNumBlocks
    integer, intent(in) :: io_splitNumBlks
    integer, intent(in) :: dataBuff(:,:)
    integer, intent(in) :: localOffset
    character (len=*), intent(in) :: description
    !Set the dimensions of the dataset
    rank = 2
    dimens2D(1) = dimens
    dimens2D(2) = io_splitNumBlks
    call h5open_f(error)
    if (error == -1) &
    
        call stop_mpi("error in subroutine writeHdf5Rank2Int. Error marker 1")
 

    call h5screate_simple_f(rank, dimens2D, dataspace, error) 
  !create the dataset
        call h5dcreate_f(fileID, description, H5T_NATIVE_INTEGER, dataspace, dataset, error)

         call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
         call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5Rank2Int. Error marker 2")

    start2d(1) = 0
    start2D(2) = localOffset 
    stride2D(:) = 1
    count2D(1) = dimens2D(1)
    count2D(2) = localNumBlocks
    
    if (localNumBlocks > 0 ) then    
        !create the hyperslab.  This will differ on the different processors
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start2D, count2D, error, &
            stride = stride2D)
!         !create the memory space
        if (error == -1)& 
            call stop_mpi("error in subroutine writeHdf5Rank2Int. Error marker 3")
     
        dimens2D(2) = localNumBlocks
        call h5screate_simple_f(rank, dimens2D, memspace, error)
         !Write the data
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Rank2Int. Error marker 4")
     
        call h5dwrite_f(dataset, H5T_NATIVE_INTEGER, dataBuff, dimens2D, error, &
            mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
         end if
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Rank2Int. Error marker 5")
     
    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
    if (error == -1)& 
        call stop_mpi("error in subroutine writeHdf5Rank2Int. Error marker 6")

end subroutine

 

!==================================================================================
! Writes rank 3 real

subroutine writeHdf5Rank3Real(myPE, fileID, dataBuff, localNumBlocks, io_splitNumBlks,&
    localOffset, description, dimen1, dimen2)

    use hdf5
    implicit none

    integer :: rank, dimen1, dimen2
    integer :: error
    integer(HID_T) :: dataset, plist_id
    integer(HSIZE_T) :: dimens3D(3)
    integer(HID_T) :: dataspace
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: start3D(3)
    integer(HSIZE_T) :: stride3D(3)
    integer(HSIZE_T) :: count3D(3)
    integer(HID_T), intent(in) :: fileID
    integer, intent(in) :: myPE 
    integer, intent(in) :: localNumBlocks
    integer, intent(in) :: io_splitNumBlks
    real, intent(in) :: dataBuff(:,:,:)
    integer, intent(in) :: localOffset
    character (len=*), intent(in) :: description
    !Set the dimensions of the dataset
    rank = 3
    dimens3D(1) = dimen1
    dimens3D(2) = dimen2
    dimens3D(3) = io_splitNumBlks
    call h5open_f(error)
    if (error == -1) &
    
        call stop_mpi("error in subroutine writeHdf5Rank3Real. Error marker 1")
 

    call h5screate_simple_f(rank, dimens3D, dataspace, error) 
  !create the dataset
        call h5dcreate_f(fileID, description, H5T_NATIVE_DOUBLE, dataspace, dataset, error)

         call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
         call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    if (error == -1) &
        call stop_mpi("error in subroutine writeHdf5Rank3Real. Error marker 3")

    start3D(1) = 0
    start3D(2) = 0
    start3D(3) = localOffset 
    
    stride3D(1) = 1
    stride3D(2) = 1
    stride3D(3) = 1
    
    count3D(1) = dimens3D(1)
    count3D(2) = dimens3D(2)
    count3D(3) = localNumBlocks
    
    if (localNumBlocks > 0 ) then    
        !create the hyperslab.  This will differ on the different processors
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start3D, count3D, error, &
            stride = stride3D)
!         !create the memory space
        if (error == -1)& 
            call stop_mpi("error in subroutine writeHdf5Rank3Real. Error marker 3")
     
        dimens3D(3) = localNumBlocks
        call h5screate_simple_f(rank, dimens3D, memspace, error)
         !Write the data
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Rank3Real. Error marker 4")
     
        call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dataBuff, dimens3D, error, &
            mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
         end if
        if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Rank3Real. Error marker 5")
     
    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
    if (error == -1)& 
        call stop_mpi("error in subroutine writeHdf5Rank3Real. Error marker 6")

end subroutine



! !===========================================================================================
! !Writes Unknowns
! 
subroutine writeHdf5Unknowns(myPE, fileID, nI, nJ, nK, unkBuf, globalVarMin, globalVarMax, &
           io_unklabel, localNumBlocks, io_splitNumBlks, localOffset)

    use hdf5
    implicit none
    
    integer(HID_T), intent(in) :: fileID, myPE
    integer, intent(in) :: localOffset
    integer, intent(in) :: nI, nJ, nK, localNumBlocks, io_splitNumBlks
    character (len=5), intent(in) :: io_unkLabel
    real, intent(in) :: globalVarMin, globalVarMax, unkBuf(nI,nJ,nK,localNumBlocks)
    integer(HID_T) :: dataset, dataspace, attribute, attributeSpace
    integer(HID_T) :: plist_id, memspace, datasetPlist
    integer :: rank, error
    integer(HSIZE_T) :: start4D(4), stride4D(4),  count4D(4)
    integer(HSIZE_T) :: dimens4D(4),  dimens1D(1), dimens5d(5)

    rank = 4
   
    dimens4D(1) = nI
    dimens4D(2) = nJ
    dimens4D(3) = nK
    dimens4D(4) = io_splitNumBlks
  
    call h5screate_simple_f(rank, dimens4D, dataspace, error)
    call h5pcreate_f(H5P_DATASET_CREATE_F, datasetPlist, error)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call h5dcreate_f(fileID, io_unkLabel, H5T_NATIVE_DOUBLE, dataspace, dataset, error, &
        dcpl_id = datasetPlist)
   
    if (error == -1)&
        call stop_mpi("error in subroutine writeHdf5Unknowns. Error marker 1")

 
    
    if (error == -1)&
        call stop_mpi("error in subroutine writeHdf5Unknowns. Error marker 2")
    !Write the attributes
    !Minimum
    dimens1D=(1)
    call h5screate_simple_f(1, dimens1D, attributeSpace, error)
   
    if (myPE == 0 .and. localOffset .ne.  0 ) then
        !do nothing
    else
        call h5acreate_f(dataset, "minimum", H5T_NATIVE_DOUBLE, attributeSpace, attribute,&
            error, H5P_DEFAULT_F)
        call h5awrite_f(attribute, H5T_NATIVE_DOUBLE, globalVarMin, dimens1D, error)
        call h5sclose_f(attributeSpace, error)
        call h5aclose_f(attribute, error)
 
    end if
    if (error == -1)&
        call stop_mpi("error in subroutine writeHdf5Unknowns. Error marker 3")


    if (error == -1)&
        call stop_mpi("error in subroutine writeHdf5Unknowns. Error marker 4")
    
   if (error == -1)&
        call stop_mpi("error in subroutine writeHdf5Unknowns. Error marker 5")
    
    !maximum
     dimens1D=(1)
    call h5screate_simple_f(1, dimens1D, attributeSpace, error)   
    if (myPE == 0 .and. localOffset .ne.  0 ) then
        !do nothing
    else
        call h5acreate_f(dataset, "maximum", H5T_NATIVE_DOUBLE, attributeSpace, attribute,&
            error, H5P_DEFAULT_F)
        call h5awrite_f(attribute, H5T_NATIVE_DOUBLE, globalVarMax, dimens1D, error)
        call h5sclose_f(attributeSpace, error)
        call h5aclose_f(attribute, error)
 
    end if
    if (error == -1)&
        call stop_mpi("error in subroutine writeHdf5Unknowns. Error marker 6")

    if (error == -1)&
        call stop_mpi("error in subroutine writeHdf5Unknowns. Error marker 8")
    
    rank = 4

    dimens4D(4) = localNumBlocks
    
    start4D(1) = 0
    start4D(2) = 0
    start4D(3) = 0
    start4D(4) = localOffset

    stride4D(1) = 1
    stride4D(2) = 1
    stride4D(3) = 1
    stride4D(4) = 1

    count4D(1) = dimens4D(1)
    count4D(2) = dimens4D(2)
    count4D(3) = dimens4D(3)
    count4D(4) = localNumBlocks
    
    dimens5D(1) = 1
    dimens5D(2) = nI
    dimens5D(3) = nJ
    dimens5D(4) = nK
    dimens5D(5) = localNumBlocks

    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start4D, count4D, error, &
        stride = stride4D)
    if (error == -1)&
        call stop_mpi("error in subroutine writeHdf5Unknowns. Error marker 9")
    rank = 4
    call h5screate_simple_f(rank, dimens4D, memspace, error)

    if (error == -1)&
        call stop_mpi("error in subroutine writeHdf5Unknowns. Error marker 10")

    call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, unkBuf, dimens4D, error,&
        mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)


    if (error == -1)&
        call stop_mpi("error in subroutine writeHdf5Unknowns. Error marker 11")
    
    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
    if (error == -1)&
        call stop_mpi("error in subroutine writeHdf5Unknowns. Error marker 12")
end subroutine
 end module ModHdf5 
