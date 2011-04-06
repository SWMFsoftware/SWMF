module ModHdf5
  use ModIO, ONLY: NamePlotDir, MaxFile, nplotvarmax

  include "mpif.h"

  character (len=80) :: filename
  integer :: plotFileNumber, checkpointFileNumber

  logical, save :: hdf5_writeRestartFile
  integer, save :: hdf5_plotArrayBegin = 1, hdf5_plotArrayEnd
  integer, save :: hdf5_currentBlockIdx = 1, localNumBlocks

  !! Constants used by the FLASH HDF5 IO package

  !! Variables for dimensional indexing
  integer, save  :: IAXIS = 1, JAXIS = 2, KAXIS = 3
  integer, save  :: LOW = 1, HIGH = 2

  !! Refinement tree information arrays
  integer, allocatable, save :: which_child(:), lrefine(:), nodetype(:), bflags(:,:)
  real, allocatable, save :: bnd_box(:,:,:), bsize(:,:), coord(:,:)

  !! FLASH-formatted unkowns array
  real, allocatable, save  :: unk(:,:,:,:,:)

  !! FLASH-formatted grid variables
  integer, allocatable, save  :: gr_gid(:,:)
  integer, save  :: gr_globalOffset, gr_globalNumBlocks

  !! FLASH IO variables
  integer, save :: io_outputSplitNum
  integer, save :: io_splitNumBlks

  character (len=20), save  :: io_flashRelease
  character (len=80), save  :: io_buildDate, io_buildDir, io_buildMachine, &
       io_fileCreationTime, io_setupTimeStamp,     &
       io_buildTimeStamp
  character (len=400), save  :: io_setupCall, io_cflags, io_fflags

  integer, save :: io_numRealParms, io_numIntParms, io_numStrParms, io_numLogParms
  integer :: io_numRealScalars, io_numIntScalars, io_numStrScalars, io_numLogScalars

  integer, save :: io_intParmValues(1), io_intScalarValues(11)
  integer, save :: io_logToIntParmValues(1), io_logToIntScalarValues(1)
  logical, save :: io_logParmValues(1), io_logScalarValues(1)
  real, save    :: io_realParmValues(1), io_realScalarValues(1)

  character (len=80), save :: io_intParmNames(1)
  character (len=80), save :: io_intScalarNames(11)    
  character (len=80), save :: io_realParmNames(1)
  character (len=80), save :: io_realScalarNames(1)
  character (len=80), save :: io_strParmNames(1)
  character (len=80), save :: io_strScalarNames(1)
  character (len=80), save :: io_strParmValues(1)
  character (len=80), save :: io_strScalarValues(1)
  character (len=80), save :: io_logParmNames(1)
  character (len=80), save :: io_logScalarNames(1)
  character (len=80), save :: io_logParmNamesPrev(1)

  character (len=4), save :: io_plotVarStr(nplotvarmax)
  character (len=4), save :: io_plotGridVarStr(nplotvarmax)
  character (len=4), save :: io_unklabels(nplotvarmax)


  character(len=80), save :: io_geometry

  integer, save :: io_nPlotVars
  integer, save :: io_comm

  integer, save :: NUNK_VARS

contains

  !=====================================================================
  ! write_plot_hdf5
  !
  !   This subroutine is called by write_plot_common and handles setting
  !   up all FLASH adapter modules and calling the appropriate
  !   subroutines to write an HDF5 plotfile in the FLASH format.
  !=====================================================================
  subroutine write_plot_hdf5(filename, plotVarNames, nPlotVar, ifile)

    use ModProcMH, ONLY : iProc, nProc
    use ModIO, ONLY : nplotvarmax
    use ModMain, ONLY : nBlockMax
    use ModSize, ONLY: nI, nJ, nK

    implicit none

    character (len=80), intent(in)  :: filename
    character (len=10), intent(in)  :: plotVarNames(nplotvarmax)
    integer, intent(in)  :: nPlotVar, ifile

    integer  :: fileID
    !-------------------------------------------------------------------------------------------

    if (iProc == 0) write (*,*) '===================================='
    if (iProc == 0) write (*,*) 'Writing HDF5 parallel plotfile'
    if (iProc == 0) write (*,*) '------------------------------------'

    !! Setup this module to write a plotfile
    hdf5_writeRestartFile = .false.
    call hdf5_setupGrid(iProc)
    call hdf5_setupIOVars(plotVarNames, nPlotVar)
    call hdf5_setupTree()

    hdf5_plotArrayEnd = nPlotVar

    if (iProc == 0) write (*,*) '  opening HDF5 file'

    !! Open the HDF5 file and write data to it
    call hdf5_initFile(fileID, filename)
    if (fileID .eq. -1) call stop_mpi('could not open HDF5 file for plotting')

    if (iProc == 0) write (*,*) '  writing data'
    call hdf5_writeData(iProc, nProc, fileID)

    if (iProc == 0) write (*,*) '  closing file'
    call hdf5_closeFile(fileID)
    if (iProc == 0) write (*,*) '===================================='

    !! Deallocate the variable array
    deallocate(unk)

  end subroutine write_plot_hdf5

  !=====================================================================
  ! write_var_hdf5
  !
  !  Handles adding one block of variable data that for plotting
  !  in the next plotfile.
  !=====================================================================
  subroutine write_var_hdf5(PlotVar, nPlotVar, iBLK)

    use ModMain, ONLY : nI, nJ, nK, nBlockMax
    use ModIO, ONLY: nPlotVarMax

    integer, intent(in) :: nPlotVar, iBLK
    real, intent(in) :: PlotVar(-1:nI+2,-1:nJ+2,-1:nK+2,nPlotVarMax)
    !----------------------------------------------------------------------

    !! Allocate the storage array if it has not been allocated for
    !! this plotfile yet
    if (.not. allocated(unk)) then
       allocate(unk(nPlotVar, nI, nJ, nK, nBlockMax))
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
  ! hdf5_setup
  !
  !   This subroutine is called by BATS_setup to allocate the dynamic
  !   arrays for various modules used by the HDF5 plotting module.
  !=====================================================================
  subroutine hdf5_setup()

    use BATL_tree, ONLY : nChild
    use BATL_size, ONLY : nDim
    use ModMain, ONLY : nBlock, nBlockMax
    use ModIO, ONLY : nPlotVarMax
    use ModSize  

    implicit none  
    !---------------------------------------------------------------------  

    !! Allocate tree data arrays
    allocate(bnd_box(2, 3, nBLK))
    allocate(bsize(3, nBLK))
    allocate(coord(3, nBLK))
    allocate(which_child(nBLK))
    allocate(lrefine(nBLK))
    allocate(nodetype(nBLK))
    allocate(bflags(1, nBLK))

    !! Allocate grid data array
    allocate(gr_gid(2*nDim+1+2**nDim, nBlockMax))

  end subroutine hdf5_setup

  !=====================================================================
  ! hdf5_writeData
  !
  !    Handles gathering any additional information needed for plotting
  !    and then calling the appropriate C functions to produce an HDF5
  !    output file.
  !=====================================================================
  subroutine hdf5_writeData(myPE, numProcs, fileID) 

    use ModSize, only : MaxBlock, nI, nJ, nK
    use ModAdvance, only : State_VGB

    implicit none
    !---------------------------------------------------------------------  

    integer, intent(in) :: myPE, fileID
    integer, intent(in) :: numProcs

    character (len=80) :: buff
    character(len=16) :: numToStr  

    integer :: blockID
    integer :: i, lb, j, AllocateStatus

    integer, allocatable :: procnumber(:) 

    ! allocate storage to hold a single variable information
    ! this should only be a small memory overhead
    integer, parameter :: single = SELECTED_REAL_KIND(p=6)

    ! allocate storage to hold the coordinate information and bounding box
    ! information
    real (kind=single) :: tmpSingle(3,MaxBlock)
    real (kind=single) :: bndSingle(2,3,MaxBlock)
    real (kind=single) :: spMax, spMin

    real, allocatable :: globalVarMin(:), globalVarMax(:)
    real, allocatable :: globalFaceXMin(:), globalFaceXMax(:)
    real, allocatable :: globalFaceYMin(:), globalFaceYMax(:)
    real, allocatable :: globalFaceZMin(:), globalFaceZMax(:)

    real, allocatable :: unkBuf(:,:,:,:,:)

    real (kind=single),allocatable :: unkt(:,:,:,:,:)

    logical :: writeGuardCells = .true.

    integer :: blkLimits(2,3), blkLimitsGC(2, 3)

    !Adjust our offset so that we can eliminate redundant data in a group
    integer :: localOffset, localRank, splitOffset
    integer :: ierr


    !! If we are writing a restart file, fill in the unkown variable array
    !!  here
    if (hdf5_writeRestartFile) then
       unk(:,:,:,:,:) = State_VGB(1:NUNK_VARS, 1:nI, 1:nJ, 1:nK, :)
    end if


    if(localNumBlocks == 0) then
       call stop_mpi("io_writeData: localNumBlocks == 0")
    end if

    !Find our local offset for split-file IO
    if (io_outputSplitNum > 1) then

       call MPI_ALLREDUCE(gr_globalOffset, splitOffset, 1, MPI_INTEGER, &
            MPI_MIN, io_comm, ierr)
       localOffset = gr_globalOffset - splitOffset
       !find number of blocks for a file
       call MPI_ALLREDUCE(localNumBlocks, io_splitNumBlks, 1, MPI_INTEGER,&
            MPI_SUM, io_comm, ierr)
    else
       localOffset = gr_globalOffset
       io_splitNumBlks = gr_globalNumBlocks
    end if

    !If we are writing a checkpoint file write all of the unk vars.  
    !If writing a plotfile then only certain variables are written
    if(hdf5_writeRestartFile) then

       !! write the header info
       call io_h5write_header(MyPE, NUNK_VARS, fileID, io_geometry, &
            io_unklabels, io_setupCall, io_fileCreationTime, io_flashRelease, &
            io_buildDate, io_buildDir, io_buildMachine, io_cflags, io_fflags, &
            io_setupTimeStamp, io_buildTimeStamp, io_outputSplitNum)

    else

       call io_h5write_header(MyPE, io_nPlotVars, fileID, io_geometry, &
            io_plotVarStr, io_setupCall, io_fileCreationTime, io_flashRelease, &
            io_buildDate, io_buildDir, io_buildMachine, io_cflags, io_fflags, &
            io_setupTimeStamp, io_buildTimeStamp, io_outputSplitNum)

    end if

    !! write the runtime parameters
    call io_h5write_runtime_parameters(mype, &
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
         io_logToIntParmValues, &
         io_outputSplitNum)

    !! write the scalars
    call io_h5write_scalars(mype, &
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
         io_logToIntScalarValues, &
         io_outputSplitNum)

    !! write the nodetype
    call io_h5write_nodetype(myPE, &
         fileID, &
         nodetype, &
         localNumBlocks, &
         io_splitNumBlks, &
         localOffset)


    !! write lrefine
    call io_h5write_lrefine(myPE, &
         fileID, &
         lrefine, &
         localNumBlocks, &
         io_splitNumBlks, &
         localOffset)


    !data struct new to PM3.  which_child definitely
    !needed for restart. Adding bflags for completeness.
    !! write which child
    call io_h5write_which_child(myPE, &
         fileID, &
         which_child, &
         localNumBlocks, &
         io_splitNumBlks, &
         localOffset)

!!!Write bflags
    call io_h5write_bflags(myPE, &
         fileID, &
         bflags, &
         localNumBlocks, &
         io_splitNumBlks, &
         localOffset)

!!!Write the global ID
    call io_h5write_gid(myPE, &
         fileID, &
         gr_gid, &
         localNumBlocks, &
         io_splitNumBlks, &
         localOffset)


    !!write the processor number. 
    allocate(procnumber(localNumBlocks))
    do lb = 1, localNumBlocks
       procnumber(lb) = MyPE
    end do

    call io_h5write_procnumber(myPE, &
         fileID,  & 
         procnumber,  & 
         localNumBlocks,  & 
         io_splitNumBlks,  & 
         localOffset)

    deallocate(procnumber)

    !! write the bnd_box
    if(hdf5_writeRestartFile) then
       call io_h5write_bndbox(myPE, &
            fileID, &
            bnd_box, &
            localNumBlocks, &
            io_splitNumBlks, &
            localOffset)
    else
       do blockID = 1, localNumBlocks
          bndSingle(:,:,blockID) = real(bnd_box(:,:,blockID), kind = single)
       enddo
       call io_h5write_bndbox_sp(myPE, fileID, & 
            bndSingle, & 
            localNumBlocks,  & 
            io_splitNumBlks,  & 
            localOffset)
    end if

    !!write the block center coordinates
    !!(not to be confused with the cell coordinates)
    if(hdf5_writeRestartFile) then
       call io_h5write_coords(myPE, &
            fileID,  & 
            coord,  & 
            localNumBlocks,  &   
            io_splitNumBlks,  & 
            localOffset)
    else
       do blockID = 1, localNumBlocks
          tmpSingle(:,blockID) = real(coord(:,blockID), kind = single)
       enddo
       call io_h5write_coords_sp(myPE, fileID,  & 
            tmpSingle,  & 
            localNumBlocks,  & 
            io_splitNumBlks,  & 
            localOffset)
    end if


    !!write the physical size of each block
    if(hdf5_writeRestartFile) then
       call io_h5write_blksize(myPE, &
            fileID,  & 
            bsize,  & 
            localNumBlocks,  & 
            io_splitNumBlks,  & 
            localOffset)
    else
       do blockID = 1, localNumBlocks
          tmpSingle(:,blockID) = real(bsize(:,blockID), kind = single)
       enddo
       ! store the block size
       call io_h5write_blksize_sp(myPE, fileID,  & 
            tmpSingle,  & 
            localNumBlocks,  & 
            io_splitNumBlks,  & 
            localOffset)
    end if

    allocate(globalVarMin(NUNK_VARS))
    allocate(globalVarMax(NUNK_VARS))

    !get the max and minimum variables
    call hdf5_getVarExtrema(NUNK_VARS, globalVarMin, globalVarMax, 2)

    !store the unknowns
    if (hdf5_writeRestartFile) then
       allocate(unkBuf(1, nI, nJ, nK, MaxBlock))
    else
       allocate(unkt(1, nI, nJ, nK, MaxBlock))
    end if

    do i = hdf5_plotArrayBegin, hdf5_plotArrayEnd 

       if(hdf5_writeRestartFile) then

          unkBuf(1,1:nI,1:nJ,1:nK,1:MaxBlock) = &
               unk(i,1:nI, 1:nJ, 1:nK,1:MaxBlock) 

          call io_h5write_unknowns(myPE, &
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
       else

          unkt(1,:,:,:,1:localNumBlocks) = real(unk(i, & 
               1:nI, & 
               1:nJ, & 
               1:nK, &
               1:localNumBlocks), & 
               kind = single)

          spMax = real(globalVarMax(i), kind = single)
          spMin = real(globalVarMin(i), kind = single)

          call io_h5write_unknowns_sp(myPE, &
               fileID, & 
               nI,   & 
               nJ,   & 
               nK, & 
               spMin, &
               spMax, &
               unkt,          & 
               io_unklabels(i),  & 
               localNumBlocks,  & 
               io_splitNumBlks,  & 
               localOffset)

       endif

    end do

    if (hdf5_writeRestartFile) then
       deallocate(unkBuf)
    else
       deallocate(unkt)
    end if
    deallocate(globalVarMin)
    deallocate(globalVarMax)

    !! Removed segment for plotting scratch variables

    !! Removed segment for plotting face variables


    return

  end subroutine hdf5_writeData

  !=====================================================================
  ! hdf5_initFile, hdf5_closeFile
  !
  !    These subroutines handle calling the appropriate C functions
  !    for opening/closing HDF5 files.
  !=====================================================================
  subroutine hdf5_initFile(fileID, filename)
    use ModProcMH, ONLY : iProc

    implicit none

    integer, intent(inout) :: fileID
    character (len=80), intent(in) :: filename
    !---------------------------------------------------------------------    
    fileID = -1
    call io_h5init_file(fileID, filename, io_comm, io_outputSplitNum)

    if(fileID == -1) then
       if (iProc == 0) write (*,*)  "Error: unable to initialize file"
       call stop_mpi("unable to initialize hdf5 file")
    end if

  end subroutine hdf5_initFile

  subroutine hdf5_closeFile(fileID)

    implicit none

    integer, intent(in) :: fileID
    call io_h5close_file(fileID)

  end subroutine hdf5_closeFile

  !=====================================================================
  ! hdf5_setupTree
  !
  !  Initializes FLASH tree arrays for plotting
  !=====================================================================
  subroutine hdf5_setupTree

    use BATL_tree,     ONLY : nChild, nNodeUsed, iTree_IA, iNode_B, Level_, &
         Parent_, Child0_
    use BATL_grid,     ONLY : CoordMin_DB, CoordMax_DB, CellSize_DB
    use ModMain,       ONLY : nBlockMax, nBlock, unusedBLK
    use ModGeometry,  ONLY : x_BLK, y_BLK, z_BLK
    use ModSize  

    implicit none

    integer  :: iBLK, iNode, iParent, iChild, idx
    integer :: halfI = 1, halfJ = 1, halfK = 1
    !---------------------------------------------------------------------
    idx = 1

    if (nI > 1) halfI = nI/2
    if (nJ > 1) halfJ = nJ/2
    if (nK > 1) halfK = nK/2

    do iBLK = 1, nBlockMax
       if (unusedBLK(iBLK)) CYCLE

       !! The coord array stores cell-centered coordinates        
       coord(IAXIS, idx) = x_BLK(halfI, halfJ, halfK, iBLK)
       coord(JAXIS, idx) = y_BLK(halfI, halfJ, halfK, iBLK)
       coord(KAXIS, idx) = z_BLK(halfI, halfJ, halfK, iBLK)


       bnd_box(LOW,:,idx) = CoordMin_DB(:,iBLK)
       bnd_box(HIGH,:,idx) = CoordMax_DB(:,iBLK)
       bsize(:,idx) = CellSize_DB(:,iBLK)

       !! nodetype holds no info for now
       nodetype(idx) = 2

       iNode = iNode_B(idx)

       !! In CRASH, refinement indices start at 0 as opposed to 1 with FLASH. Add
       !! 1 to the value to comply with the FLASH standard.
       lrefine(idx) = iTree_IA(Level_, iNode) + 1

       !! Find which child this node is
       iParent = iTree_IA(Parent_, iNode)
       which_child(:) = -1
       if (iParent >= 0) then
          do iChild = 1, nChild
             if (iTree_IA(Child0_+iChild, iParent) .eq. iNode) which_child(idx) = iChild
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

    use BATL_geometry, ONLY : TypeGeometry
    use BATL_size, ONLY : nDim, nDimAmr
    use BATL_tree, ONLY : nInfo => ChildLast_, nNode
    use ModMain, ONLY : nBlockAll, Time_Simulation
    use ModVarIndexes, ONLY : NameVar_V, nVar, nFluid

    include "mpif.h"

    !! Arguments
    character (len=10), intent(in), optional :: plotVarNames(nplotvarmax)
    integer, intent(in), optional :: nPlotVar

    integer  :: iVar
    !-------------------------------------------------------------------------------------------
    io_flashRelease = 'CRASH'

    io_geometry = TypeGeometry

    !! Set compiler and setup string info as empty
    io_buildDate = ''
    io_buildDir = ''
    io_buildMachine = ''
    io_fileCreationTime = ''
    io_setupTimeStamp = ''
    io_buildTimeStamp = ''
    io_setupCall = ''
    io_cflags =''
    io_fflags = ''

    io_outputSplitNum = 1 ! file splitting not supported

    !! FLASH makes a distinction between plot var labels and unk var labels, but
    !! CRASH will treat both the same
    if (.not. hdf5_writeRestartFile) then
       do iVar = 1, nPlotVar
          io_plotVarStr(iVar) = trim(plotVarNames(iVar)(1:4))
          io_unklabels(iVar) = trim(plotVarNames(iVar)(1:4))
       end do
       io_nPlotVars = nPlotVar
    else
       io_nPlotVars = nVar + nFluid
       do iVar = 1, io_nPlotVars
          io_plotVarStr(iVar) = trim(NameVar_V(iVar))
          io_unklabels(iVar) = trim(NameVar_V(iVar))
       end do
    end if

    !! Output of runtime parameters, scalars, and logfiles is not implemented yet,
    !! so clear the corresponding variables
    io_numRealParms = 0
    io_numIntParms = 0
    io_numStrParms = 0
    io_numLogParms = 0
    io_numRealScalars = 1
    io_numIntScalars = 11
    io_numStrScalars = 0
    io_numLogScalars = 0
    ! probably don't need to alloc string and value arrays...

    !! The FLASH Visit plugin requires some scalar metadata. This should be moved
    !! to the empty io_prepareLists subroutine eventually.
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

    !! Restart files require node and tree info
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

    io_chkGuardCells = .false.

    io_comm = MPI_COMM_WORLD ! we need more MPI initialization for split file IO

    NUNK_VARS = nVar + nFluid

  end subroutine hdf5_setupIOVars

  !=====================================================================
  ! hdf5_setupGrid
  !
  !   Initializes the grid data arrays derived from the FLASH code.
  !=====================================================================
  subroutine hdf5_setupGrid(myPE)

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

    implicit none

    integer, intent(IN)  :: myPE

    integer :: i, j, k, iNode, iBLK, iGid, ierr, idx, nUnused
    integer :: bn, pe
    integer :: localIdx(nBlockMax)
    integer, allocatable :: blocksPerProc(:), offsetPerProc(:)
    integer, allocatable :: localIdx_pe(:,:)
    !---------------------------------------------------------------------  

    gr_globalNumBlocks = nBlockAll

    allocate(blocksPerProc(nProc))
    allocate(offsetPerProc(0:nProc-1))
    allocate(localIdx_pe(0:nProc-1, nBlockMax))

    !! Get the number of blocks local to this processor and the mapping
    !! from block index to output index (where unused blocks are skipped)
    call hdf5_getLocalBlks(localNumBlocks, localIdx)

    !! Collect results across all processors
    call MPI_Allgather(localNumBlocks, 1, MPI_INTEGER, blocksPerProc, 1, MPI_INTEGER, iComm, ierr)
    call MPI_Allgather(localIdx, nBlockMax, MPI_INTEGER, localIdx_pe, nBlockMax, MPI_INTEGER, iComm, ierr)

    offsetPerProc(:) = 0

    do i = 0, nProc-2
       offsetPerProc(i+1) = offsetPerProc(i) + blocksPerProc(i+1)
    end do

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

       !! Begin by filling in the neighbor cells for this cell

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

       !! Store parent's node number and children's node numbers
       do i = 0, 2**nDim
          !! No AMR tree info stored yet...
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
  ! hdf5_getLocalBlks
  !
  !   Returns the number of blocks used on the local processor as well
  !   as an array mapping the index of each used block to its output
  !   array index (where unused blocks are skipped). 
  !=====================================================================
  subroutine hdf5_getLocalBlks(numBlocks, localIdx)

    use ModMain, ONLY : nBlockMax, unusedBLK, nBlockALL
    use ModProcMH,  ONLY : iProc

    implicit none

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
  ! hdf5_getBlkIndexLimits
  !
  !   Returns the index limits for a block, both with and without guard
  !    cells.
  !=====================================================================
  subroutine hdf5_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, gridDataStruct)

    use ModSize

    implicit none

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
  ! hdf5_getVarExtrema
  !
  !   Returns the global variable extrema for each unkown.
  !=====================================================================
  subroutine hdf5_getVarExtrema(nvars, globalVarMin, globalVarMax, gridDataStruct)

    use ModMain
    use ModAdvance, ONLY : State_VGB
    use ModSize

    implicit none

    integer, intent(in) :: nvars, gridDataStruct
    real, DIMENSION(nvars), INTENT(out) :: globalVarMin, globalVarMax

    real, allocatable  :: varMin(:), varMax(:)
    integer  :: iBLK, i, j, k, iVar, ierr
    !-----------------------------------------------------------------------------

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
            MPI_MIN, MPI_COMM_WORLD, ierr)

       call MPI_AllReduce(varMax, globalVarMax, nvars, MPI_DOUBLE_PRECISION, &
            MPI_MAX, MPI_COMM_WORLD, ierr)

       deallocate(varMin)
       deallocate(varMax)

    end if
    return

  end subroutine hdf5_getVarExtrema

end module ModHdf5
