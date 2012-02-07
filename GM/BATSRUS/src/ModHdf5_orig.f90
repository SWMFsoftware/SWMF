module ModHdf5

  use ModProcMH, ONLY: iProc, nProc, iComm 
  use hdf5

  implicit none
  SAVE


  private ! except
  public write_var_hdf5
  public write_plot_hdf5

  real,  allocatable :: PlotVarIdx(:,:,:,:,:)
  integer, allocatable :: UsedBlocks(:)
  integer, allocatable :: usednodes(:)
  integer, allocatable :: UsedBlocksTemp(:)

  integer, parameter:: lNameVar = 10
  integer, parameter:: lNameH5 = lNameVar + 1

  character (len=lNameH5), allocatable :: UnknownNameArray(:)

  integer :: nBlocksUsed, nBlkUsedGlobal, iBlk

  real, allocatable :: Position_DA(:,:)

contains

  subroutine write_var_hdf5(PlotVar, nPlotVar, H5Index, iBlk)
    use ModMain, only : nBlockMax, nI,nJ,nK


    !  Handles adding one block of variable data that for plotting
    !  in the next plotfile.

    integer, intent(in) :: nPlotVar, H5Index, iBlk
    real, intent(in) :: PlotVar(-1:nI+2,-1:nJ+2,-1:nK+2,nPlotVar)
    integer:: ii, jj, kk, iPlot, error
    !----------------------------------------------------------------------
    ! Allocate the storage array if it has not been allocated for
    ! this plotfile yet

!!! Use nBlock if possible

    if (.not. allocated(PlotVarIdx)) then
       allocate(PlotVarIdx(nI, nJ, nK, nBlockMax, nPlotVar))
       allocate(UsedBlocksTemp(nBlockMax))

    end if


    !translates iBlk to index for plotfile where unused blocks are skipped  
    UsedBlocksTemp(H5Index) = iBLK
    nBlocksUsed = H5Index

!!! Efficient loop ordering: do k = ...; do j = ...; do i = ...

    do iPlot = 1, nPlotVar
       do ii = 1, nI; do jj = 1, nJ; do kk = 1, nK
          PlotVarIdx(ii, jj, kk, H5Index, iPlot) = &
               PlotVar(ii, jj, kk, iPlot)


       end do; end do; end do
    end do

  end subroutine write_var_hdf5
  !=============================================================
  !=============================================================
!   integer function TestFunction2D(a,b)
!     integer, intent(in) :: a(2),b(2)
! 
!     TestFunction2D = 0
!     if(a(2) < b(2)) then
!        TestFunction2D = -1
!     elseif(b(2) < a(2)) then
!        TestFunction2D = 1
!     else
!        if(a(1) < b(1)) then
!           TestFunction2D = -1
!        elseif(b(1) < a(1)) then
!           TestFunction2D = 1
!        end if
!     end if
!   end function TestFunction2D
! 
!   !=============================================================
!   !=============================================================
!   integer function TestFunction3D(a,b)
!     integer, intent(in) :: a(2),b(2)
! 
!     TestFunction3D = 0
!     if(a(3) < b(3)) then
!        TestFunction3D = -1
!     elseif(b(3) < a(3)) then
!        TestFunction3D = 1
!     else
!        if(a(2) < b(2)) then
!           TestFunction3D = -1
!        elseif(b(2) < a(2)) then
!           TestFunction3D = 1
!        else
!           if(a(1) < b(1)) then
!              TestFunction3D = -1
!           elseif(b(1) < a(1)) then
!              TestFunction3D = 1
!           end if
!        end if
!     end if
!   end function TestFunction3D
!   !=============================================================
!   !=============================================================
!   integer function TestFunction1D(a,b)
!     integer, intent(in) :: a(2),b(2)
! 
!     TestFunction1D = 0
!     if(a(1) < b(1)) then
!        TestFunction1D = -1
!     elseif(b(1) < a(1)) then
!        TestFunction1D = 1
!     end if
!   end function TestFunction1D
    
! 
!     logical function sort_test(i,j)
!         use BATL_tree, only : iTree_IA, Level_, Coord1_,Coord2_,Coord3_,&
!             nDim, MaxLevel 
!         integer, intent(in) :: i, j
!         integer :: valueI(3), valueJ(3), iGlobal, jGlobal
!         
!         iGlobal = GlobalUsedNodes(i)
!         jGlobal = GlobalUsedNodes(j)
!         valueI = 
!         if (nDim == 3) then
! 
            
       
  logical function sort_test(i, j)
    use BATL_lib, ONLY : nDim 

    integer, intent(in) :: i,j
    if(ndim == 3) then
        sort_test = 0
        if(Position_DA(3,i) < Position_DA(3,j)) then
           sort_test = .false.
        elseif(Position_DA(3,j) < Position_DA(3,i)) then
           sort_test = .true.
        else
           if(Position_DA(2,i) < Position_DA(2,j)) then
              sort_test = .false.
           elseif(Position_DA(2,j) < Position_DA(2,i)) then
              sort_test = .true.
           else
              if(Position_DA(1,i) < Position_DA(1,j)) then
                 sort_test = .false.
              elseif(Position_DA(1,j) < Position_DA(1,i)) then
                 sort_test = .true.
              end if
           end if
        end if
    elseif(nDim == 2) then
       if(Position_DA(2,i) < Position_DA(2,j)) then
          sort_test = .false.
       elseif(Position_DA(2,j) < Position_DA(2,i)) then
          sort_test = .true.
       else
          if(Position_DA(1,i) < Position_DA(1,j)) then
             sort_test = .false.
          elseif(Position_DA(1,j) < Position_DA(1,i)) then
             sort_test = .true.
          end if
       end if
    elseif(nDim == 1) then
        if(Position_DA(1,i) < Position_DA(1,j)) then
            sort_test = .false.
        elseif(Position_DA(1,j) < Position_DA(1,i)) then
             sort_test = .true.
        end if

    else
        sort_test = .false.
    end if
    end function sort_test
  !=================================================================
  !================================================================

  subroutine write_plot_hdf5(filename, plotVarNames, plotVarUnits,&
    nPlotVar,xmin, xmax, ymin, ymax, zmin, zmax, nBLKcells)

    use BATL_tree, only: iNode_B, iTree_IA, Coord1_, Coord2_, Coord3_,&
        Level_,iMortonNode_A, MaxNode, get_tree_position, nNodeUsed,&
        Block_, Proc_,iNodeMorton_I, Unused_BP
    use ModMpi
    use ModMain, only : nI,nJ,nK, nBlockMax
    use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK
    use BATL_lib, ONLY : nDim, MaxDim
    use ModSort, only : sort_quick_func
    use BATL_lib, ONLY : CoordMin_DB

    integer,                 intent(in):: nPlotVar
    character(len=80),       intent(in):: filename
    character(len=lNameVar), intent(in):: plotVarNames(nPlotVar)
    character(len=lNameVar),  intent(in):: plotVarUnits(nPlotVar)
    real,                    intent(in):: xMin, xMax, yMin, yMax, zMin, zMax
    integer,                 intent(in):: nBLKcells

    integer :: offsetPerProc(0:nProc-1), blocksPerProc(0:nProc-1)
    integer :: fileID, error, procIdx, lastCoord, iNode, iLen, iVar 
    integer :: labelLeng, i, PosMinType, iMorton
    real, allocatable :: coordinates(:,:)
    integer, allocatable :: minLogicalExtents(:,:)
    integer, allocatable :: procNum(:)
    integer, allocatable :: VisItIndx(:)
    integer, allocatable :: NumUnsedBlkTo(:,:)
    real :: PositionMin_D(MaxDim), PositionMax_D(MaxDim)
    !--------------------------------------------------------------------------

    if (iProc == 0) write (*,*) '===================================='
    if (iProc == 0) write (*,*) 'Writing BATL hdf5 parallel plotfile'
    if (iProc == 0) write (*,*) '------------------------------------'
    if (iProc == 0) write (*,*) '  opening hdf5 file'


    !initialize the Fortran Hdf5 library
    call h5open_f(error)
    ! Open the hdf5 file and write data to it

    !---------------------------------------------------------------------    
    fileID = -1
    call hdf5_init_file(fileID, filename)
    if(fileID == -1) then
       if (iProc == 0) write (*,*)  "Error: unable to initialize file"
       call stop_mpi("unable to initialize hdf5 file")
    end if
    !get processor block offsets for writes

    call MPI_Allgather(nBlocksUsed, 1, MPI_INTEGER, blocksPerProc, &
         1, MPI_INTEGER, iComm, error)

    offsetPerProc(0) = 0
    nBlkUsedGlobal = blocksPerProc(0) 
    if (nProc > 1) then
       do procIdx = 1, nProc-1
          offsetPerProc(procIdx) = offsetPerProc(procIdx-1) &
               + blocksPerProc(procIdx-1)
          nBlkUsedGlobal = nBlkUsedGlobal + blocksPerProc(procIdx)
       end do
    end if

    !The UsedBlocksTemp array will not serve our purpose because it will
    !have trailing negative numbers.  Move it to an array of the correct
    !size and use it to create an array of used Nodes.
    allocate(UsedBlocks(nBlocksUsed))
    UsedBlocks(1:nBlocksUsed) = UsedBlocksTemp(1:nBlocksUsed)
    deallocate(UsedBlocksTemp)

    allocate(UsedNodes(nBlocksUsed))
    UsedNodes(1:nBlocksUsed) = iNode_B(UsedBlocks(1:nBlocksUsed))

    allocate(NumUnsedBlkTo(nBlockMax,0:nProc-1))
    NumUnsedBlkTo = 0
    do procIdx = 0, nProc-1
        if(Unused_BP(1, procIdx))&
            NumUnsedBlkTo(1, procIdx) = 1
        do iBlk = 2, nBlockMax
            if (Unused_BP(iBlk,procIdx)) then
                NumUnsedBlkTo(iBlk,procIdx)=NumUnsedBlkTo(iBlk-1,ProcIdx)&
                    +1
            else
                NumUnsedBlkTo(iBlk,procIdx)=NumUnsedBlkTo(iBlk-1,procIdx)
            end if
        end do
    end do

    allocate(Position_DA(nDim,nNodeUsed))
    do iMorton = 1, nNodeUsed
        iNode = iNodeMorton_I(iMorton)
        iBlk = iTree_IA(Block_,iNode)
        procIdx = iTree_IA(Proc_,iNode) 
        i = offsetPerProc(procIdx) + iBlk - NumUnsedBlkTo(iBlk,procIdx)
        call get_tree_position(iNode, PositionMin_D, PositionMax_D)
        Position_DA(:,i) = PositionMin_D(1:nDim)
    end do
    deallocate(NumUnsedBlkTo)

    allocate(VisItIndx(nBlkUsedGlobal))        

    call sort_quick_func(nBlkUsedGlobal, sort_test, VisItIndx)

    deallocate(Position_DA)
    call writeHdf5Rank1Integer(fileID, VisItIndx, nBlkUsedGlobal,&
       nBlkUsedGlobal, 0,"VisIt Index")




    deallocate(VisItIndx)
    

    allocate(unknownNameArray(nPlotVar))
    do iVar = 1, nPlotVar
       UnknownNameArray(iVar) = plotVarNames(iVar)
       !The VisIt plugin needs null padded names.
       labelLeng = len_trim(UnknownNameArray(iVar))
       do iLen = labelLeng + 1,lNameH5 
          UnknownNameArray(iVar)(iLen:iLen) = CHAR(0)
       end do
    end do

    call write_plot_string(nPlotVar,UnknownNameArray,"plotVarNames",&
         fileID)

    call write_plot_vars(nPlotVar, offsetPerProc(iProc), fileID)

    do iVar = 1, nPlotVar
       UnknownNameArray(iVar) = plotVarUnits(iVar)
       !The VisIt plugin needs null padded names.
       labelLeng = len_trim(UnknownNameArray(iVar))
       do iLen = labelLeng + 1,lNameH5 
          UnknownNameArray(iVar)(iLen:iLen) = CHAR(0)
       end do
    end do
    call write_plot_string(nPlotVar, UnknownNameArray, "plotVarUnits",&
        fileID)


    deallocate(PlotVarIdx)
    deallocate(unknownNameArray)
    call write_bounding_box(offsetPerProc(iProc), fileID)
    !    call write_block_center_coords(offsetPerProc(iProc),fileID)

    !write_min_logical_extents also write the VisIt index
    allocate(coordinates(ndim,nBlocksUsed))
    do iBlk = 1,nBlocksUsed
       coordinates(1, iBlk) = 0.5*(x_BLK(nI,nJ,nK,UsedBlocks(iBLK))+&
            x_BLK(1,1,1,UsedBlocks(iBLK)))
       if (nDim >=2)&
            coordinates(2, iBlk) = 0.5*(y_BLK(nI,nJ,nK,&
            UsedBlocks(iBLK))+ y_BLK(1,1,1,UsedBlocks(iBLK)))
       if (ndim ==3)&
            coordinates(3, iBlk) = 0.5*&
            (z_BLK(nI,nJ,nK,UsedBlocks(iBLK))+&
            z_BLK(1,1,1,UsedBlocks(iBLK)))
    end do
    !
    call writeHdf5Rank2Real(fileID, coordinates, nBlocksUsed,&
         nBlkUsedGlobal, offsetPerProc(iProc), "coordinates", nDim)

    deallocate(coordinates)
    allocate(minLogicalExtents(ndim, nBlocksUsed))
    lastCoord = Coord1_ + nDim-1
    do iBLK = 1, nBlocksUsed
       iNode = UsedNodes(iBLK)
       minLogicalExtents(:,iBlk) = iTree_IA(Coord1_:LastCoord,iNode)&
            -1.0
    end do
    !
    call writeHdf5Rank2Integer(fileID, minLogicalExtents, nBlocksUsed,&
         nBlkUsedGlobal, offsetPerProc(iProc), "MinLogicalExtents", nDim)
    !    call write_iMortonNode_A(fileID,offsetPerProc(iProc))

    deallocate(minLogicalExtents)
    call writeHdf5Rank1Integer(fileID, iMortonNode_A(UsedNodes), &
         nBlocksUsed, nBlkUsedGlobal, offsetPerProc(iProc),"iMortonNode_A")

    deallocate(usedBlocks)
    call writeHdf5Rank1Integer(fileID, iTree_IA(Level_,UsedNodes), &
         nBlocksUsed, nBlkUsedGlobal, offsetPerProc(iProc),"refine level")
    deallocate(usedNodes)

    allocate(procNum(nBlocksUsed))
    procNum = iProc
    call writeHdf5Rank1Integer(fileID, procNum,nBlocksUsed, &
         nBlkUsedGlobal, offsetPerProc(iProc),"Processor Number")
    deallocate(procNum)

    call write_integer_metadata(fileID, nPlotVar)
    call write_real_metadata(FileID,xmin,xmax,ymin,ymax, zmin,zmax)
    if (iProc == 0) write (*,*) '  closing file'
    !Close the file
    call h5garbage_collect_f(error)
    call h5fclose_f(fileID,error)
    !closing the hdf5 interface
    call h5close_f(error)
    if (error == -1) write (*,*) 'h5fclose_f failed'

    if (iProc == 0) write (*,*) '===================================='

    ! Deallocate the variable array

  end subroutine write_plot_hdf5

  !=====================================================================
  !=====================================================================

  subroutine writeHdf5Rank1Integer(fileID, dataBuff, localNumBlocks,&
       globalNumBlocks, localOffset, description)

    implicit none

    integer :: rank, dimens
    integer :: error
    integer(HID_T) :: dataset, plist_id
    integer(HSIZE_T) :: dimens1D(2)
    integer(HID_T) :: dataspace
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: start1D(2)
    integer(HSIZE_T) :: stride1D(2)
    integer(HSIZE_T) :: count1D(2)
    integer(HID_T), intent(in) :: fileID
    integer, intent(in) :: localNumBlocks, globalNumBlocks
    integer, intent(in) :: dataBuff(:)
    integer, intent(in) :: localOffset
    character (len=*), intent(in) :: description
    !Set the dimensions of the dataset
    rank = 1
    dimens1D(1) = globalNumBlocks
    call h5open_f(error)
    if (error == -1) &
         
         call stop_mpi("error in subroutine writeHdf5Rank1Real. Error marker 1")


    call h5screate_simple_f(rank, dimens1D, dataspace, error) 
    !create the dataset
    call h5dcreate_f(fileID, description, H5T_NATIVE_INTEGER, dataspace, dataset, error)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Rank1Real. Error marker 2")

    start1D(1) = localOffset 
    stride1D(:) = 1
    count1D(1) = localNumBlocks

    if (localNumBlocks > 0 ) then    
       !create the hyperslab.  This will differ on the different processors
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start1D, count1D, error, &
            stride = stride1D)
       !         !create the memory space
       if (error == -1)& 
            call stop_mpi("error in subroutine writeHdf5Rank1Real. Error marker 3")

       dimens1D(1) = localNumBlocks
       call h5screate_simple_f(rank, dimens1D, memspace, error)
       !Write the data
       if (error == -1) &
            call stop_mpi("error in subroutine writeHdf5Rank1Real. Error marker 4")

       call h5dwrite_f(dataset, H5T_NATIVE_INTEGER, dataBuff, dimens1D, error, &
            mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
    end if
    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Rank1Real. Error marker 5")

    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)

    if (error == -1)& 
         call stop_mpi("error in subroutine writeHdf5Rank1Real. Error marker 6")

  end subroutine writeHdf5Rank1Integer

  !=====================================================================
  !=====================================================================

  subroutine writeHdf5Rank2Integer(fileID, dataBuff, localNumBlocks,&
       globalNumBlocks, localOffset, description, dimens)

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
    integer, intent(in) :: localNumBlocks, globalNumBlocks
    integer, intent(in) :: dataBuff(:,:)
    integer, intent(in) :: localOffset
    character (len=*), intent(in) :: description
    !Set the dimensions of the dataset
    rank = 2
    dimens2D(1) = dimens
    dimens2D(2) = globalNumBlocks
    call h5open_f(error)
    if (error == -1) &
         
         call stop_mpi("error in subroutine writeHdf5Rank2Real. Error marker 1")


    call h5screate_simple_f(rank, dimens2D, dataspace, error) 
    !create the dataset
    call h5dcreate_f(fileID, description, H5T_NATIVE_INTEGER, dataspace, dataset, error)

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

       call h5dwrite_f(dataset, H5T_NATIVE_INTEGER, dataBuff, dimens2D, error, &
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

  end subroutine writeHdf5Rank2Integer

  !======================================================================
  !=====================================================================

  subroutine writeHdf5Rank2Real(fileID, dataBuff, localNumBlocks,&
       globalNumBlocks, localOffset, description, dimens)

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
    integer, intent(in) :: localNumBlocks, globalNumBlocks
    real, intent(in) :: dataBuff(:,:)
    integer, intent(in) :: localOffset
    character (len=*), intent(in) :: description
    !Set the dimensions of the dataset
    rank = 2
    dimens2D(1) = dimens
    dimens2D(2) = globalNumBlocks
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

  end subroutine writeHdf5Rank2Real

  !======================================================================
  !=====================================================================
  subroutine hdf5_init_file(fileID, filename)

    use ModMpi

    character (len=80), intent(in) :: filename
    integer :: error, accTemplate
    integer, intent(inout) :: fileID
    INTEGER :: mpierror       ! MPI error flag
    INTEGER :: dxpList
!!!    INTEGER :: mpi_size, mpi_rank, comm

    !Initalize currentBlock andplotArrayBegin for this file

    !    call h5open_f(error)                    

    !create MPI info Object

    !Create file access propertty list
    call h5pcreate_f(H5P_FILE_ACCESS_F, accTemplate, error)
    CALL h5pset_fapl_mpio_f(accTemplate, iComm, MPI_INFO_NULL, error)
    ! Create the file collectively.

    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, fileID, error,&
         access_prp = accTemplate)
    CALL h5pclose_f(accTemplate, error)
    if (error == -1) &
         call stop_mpi(&
         "error in subroutine hdf5_init_file. Error marker 1")

    if (error == -1) fileID = -1
  end subroutine hdf5_init_file

  !=====================================================================
  !=====================================================================

  subroutine write_plot_string(nPlotVar,dataBuff, description, fileID)

    integer, intent(in) :: nPlotVar
    character (len=11), intent(in) :: dataBuff(nPlotVar)
    integer (HID_T), intent(in) :: fileID
    character (len=*), intent(in) :: description
    integer (HID_T) :: dataType, dataset, dataspace
    integer (HSIZE_T) :: sizeString, dimens1D(1) 

    integer :: error, rank, iLen, iVar, labelleng

    rank = 1  
    dimens1d(1) = nPlotVar

    sizeString = lNameH5


    call h5tcopy_f(H5T_NATIVE_CHARACTER, datatype, error)
    call h5tset_size_f(datatype, sizeString, error)
    if (error == -1) &
         call stop_mpi(&
         "error in subroutine writeHdf5Header. Error marker 10")
    call h5screate_simple_f(rank, dimens1D, dataspace, error)
    call h5dcreate_f(&
         fileID, description, datatype, dataspace, dataset, error)
    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Header. Error marker 11")

    if (iProc == 0) &
         call h5dwrite_f(dataset, datatype, UnknownNameArray, dimens1D, error, H5S_ALL_F,&
         H5S_ALL_F, H5P_DEFAULT_F)
    if (error == -1) &
         call stop_mpi("error in subroutine writeHdf5Header. Error marker 11")



    call h5tclose_f(datatype, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)

    if (error == -1)& 
         call stop_mpi(&
         "error in subroutine writeHdf5Header. Error marker 13")
  end subroutine write_plot_string
  !=====================================================================
  !=====================================================================

  subroutine write_plot_vars(nPlotVar, iProcOffset,fileID)

    use ModMain, only : nI,nJ,nK

    integer, intent(in) :: nPlotVar, iProcOffset, fileID
    real :: globalVarMin, globalVarMax
    integer :: iVar,rank, error
    integer(HSIZE_T) :: dimens4D(4), start4D(4), stride4D(4)
    integer(HSIZE_T) :: dimens1D(1)
    integer(HID_T) :: dataset, dataspace, attribute, attributeSpace
    integer(HID_T) :: memspace, plist_id



    do iVar=1, nPlotVar
       call get_var_extrema_hdf5(iVar, globalVarMin, globalVarMax)
       rank = 4
       dimens4D(1) = nI
       dimens4D(2) = nJ
       dimens4D(3) = nK
       dimens4D(4) = nBlkUsedGlobal

       call h5screate_simple_f(rank, dimens4D, dataspace, error) 
       call h5dcreate_f(fileID, UnknownNameArray(iVar),&
            &H5T_NATIVE_DOUBLE,dataspace, dataset, error)
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
       call h5pset_dxpl_mpio_f(plist_id,&
            H5FD_MPIO_COLLECTIVE_F, error)
       !Write the attributes
       !Minimum
       dimens1D=(1)
       call h5screate_simple_f(1, dimens1D, attributeSpace, error)

       call h5acreate_f(dataset, "minimum", H5T_NATIVE_DOUBLE, attributeSpace, attribute,&
            error, H5P_DEFAULT_F)
       call h5awrite_f(attribute, H5T_NATIVE_DOUBLE, globalVarMin, dimens1D, error)
       call h5sclose_f(attributeSpace, error)
       call h5aclose_f(attribute, error)
       !Maximum
       dimens1D=(1)
       call h5screate_simple_f(1, dimens1D, attributeSpace, error)

       call h5acreate_f(dataset, "maximum", H5T_NATIVE_DOUBLE, attributeSpace, attribute,&
            error, H5P_DEFAULT_F)
       call h5awrite_f(attribute, H5T_NATIVE_DOUBLE, globalVarMax, dimens1D, error)
       call h5sclose_f(attributeSpace, error)
       call h5aclose_f(attribute, error)

       dimens4D(4) = nBlocksUsed

       start4D(1) = 0
       start4D(2) = 0
       start4D(3) = 0
       start4D(4) = iProcOffset

       stride4D = 1
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F,&
            start4D, dimens4D, error, stride = stride4D)
       if (error == -1)&
            call stop_mpi(&
            "error in subroutine init_grid. Error marker 1")
       call h5screate_simple_f(rank, dimens4D, memspace, error)

       if (error == -1)&
            call stop_mpi("error in subroutine init_grid. Error marker 2")

       call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, &
            PlotVarIdx(:,:,:,:,iVar),&
            dimens4D, error, mem_space_id = memspace,&
            file_space_id = dataspace)
       if (error == -1)&
            call stop_mpi(&
            "error in subroutine init_grid. Error marker 3")

       call h5sclose_f(memspace, error)
       if (error == -1)&
            call stop_mpi(&
            "error in subroutine init_grid. Error marker 4")

       call h5pclose_f(plist_id, error)
       call h5sclose_f(dataspace, error)
       call h5dclose_f(dataset, error)

    end do

  end subroutine write_plot_vars


  !=================================================================
  !================================================================

  subroutine get_var_extrema_hdf5( iVar, globalVarMin, globalVarMax)

    use ModMpi

    ! Returns the global variable extrema for each unkown.

    integer, intent(in) :: iVar
    real, intent(out) :: globalVarMin, globalVarMax

    real :: varMin, varMax
    integer :: error
    !------------------------------------------------------------------


    varMin = huge(varMin)
    varMax = -huge(varMax)

    varMin = minval(PlotVarIdx(:,:,:,:,iVar))
    varMax = maxval(PlotVarIdx(:,:,:,:,iVar))

    call MPI_AllReduce(varMin, globalVarMin, 1,&
         MPI_DOUBLE_PRECISION, MPI_MIN, iComm, error)

    call MPI_AllReduce(varMax, globalVarMax, 1,&
         MPI_DOUBLE_PRECISION, MPI_MAX, iComm, error)

  end subroutine get_var_extrema_hdf5

  !=================================================================
  !================================================================

  subroutine write_bounding_box(iProcOffset, fileID)

    use BATL_lib, ONLY : CoordMin_DB, CoordMax_DB, nDim
    use ModProcMH, ONLY : iProc

    integer, intent(in) :: iProcOffset
    integer (HID_T), intent(in) :: fileID

    integer(HID_T) :: dataset, dataspace, attribute, attributeSpace
    integer(HID_T) :: plist_id, memspace

    integer(HSIZE_T) :: start3D(3), stride3D(3), dimens3D(3)

    integer :: error, rank


    rank = 3
    dimens3D(1) = 2
    dimens3D(2) = nDim
    dimens3D(3) = nBlkUsedGlobal

    call h5screate_simple_f(rank, dimens3D, dataspace, error) 
    call h5dcreate_f(fileID, "bounding box",&
         &H5T_NATIVE_DOUBLE,dataspace, dataset, error)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id,&
         H5FD_MPIO_COLLECTIVE_F, error)

    dimens3D(1) = 1        
    dimens3D(3) = nBlocksUsed

    start3D(1) = 0
    start3D(2) = 0
    start3D(3) = iProcOffset

    stride3D(1) = 1
    stride3D(2) = 1
    stride3D(3) = 1

    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F,&
         start3D, dimens3D, error, stride = stride3D)
    if (error == -1)&
         call stop_mpi(&
         "error in subroutine init_grid. Error marker 1")
    call h5screate_simple_f(rank, dimens3D, memspace, error)

    if (error == -1)&
         call stop_mpi("error in subroutine init_grid. Error marker 2")

    call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, &
         CoordMin_DB(1:nDim, UsedBlocks),&
         dimens3D, error, mem_space_id = memspace,&
         file_space_id = dataspace)
    if (error == -1)&
         call stop_mpi(&
         "error in subroutine init_grid. Error marker 3")

    call h5sclose_f(memspace, error)
    if (error == -1)&
         call stop_mpi(&
         "error in subroutine init_grid. Error marker 4")

    start3D(1) = 1
    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F,&
         start3D, dimens3D, error, stride = stride3D)
    if (error == -1)&
         call stop_mpi(&
         "error in subroutine init_grid. Error marker 1.1")
    call h5screate_simple_f(rank, dimens3D, memspace, error)

    if (error == -1)&
         call stop_mpi("error in subroutine init_grid. Error marker 2.1")


    call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, &
         CoordMax_DB(1:nDim, UsedBlocks),&
         dimens3D, error, mem_space_id = memspace,&
         file_space_id = dataspace)
    if (error == -1)&
         call stop_mpi(&
         "error in subroutine init_grid. Error marker 3.1")

    call h5sclose_f(memspace, error)
    if (error == -1)&
         call stop_mpi(&
         "error in subroutine init_grid. Error marker 4.1")


    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)


  end subroutine write_bounding_box
  !=======================================================================
  !=======================================================================


  subroutine write_real_metadata(FileID,xmin,xmax,ymin,ymax, zmin,zmax)

    use ModMain, only : Time_Simulation, CodeVersion
    use ModProcMH, only : iProc
    use BATL_lib, only : nLevel, nDim

    integer (HID_T), intent(in) :: fileID
    integer(HID_T) :: dataset, dataspace, attribute, attributeSpace
    integer :: levelFactorArray(0:nLevel, nBlocksUsed)
    integer :: minLogicalExtents(ndim,nBlocksUsed)
    integer(HSIZE_T) :: dimens1D(1)
    integer :: error, rank, iData

    real, intent(in) :: xmin,xmax,ymin,ymax, zmin,zmax

    real :: RealMetaData(7)
    !    allocate(attName(nAtts))
    iData = 1
    !    attName(1) = 'Simulation Time'
    RealMetaData(iData) = Time_Simulation
    iData = iData + 1   
    !    attName(2) = 'xmin'
    RealMetaData(iData) = xmin
    iData = iData + 1
    !    attName(3) = 'xmax'
    RealMetaData(iData) = xmax
    iData = iData + 1
    !    attName(4) = 'ymin'
    RealMetaData(iData) = ymin
    iData = iData + 1
    !    attName(5) = 'ymax'
    RealMetaData(iData) = ymax
    iData = iData + 1
    !    attName(6) = 'zmin'
    RealMetaData(iData) = zmin

    iData = iData + 1
    !    attName(7) = 'zmax'
    RealMetaData(iData) = zmax
    !-------------------------------------------------------------------
    !write the real Metadata
    rank = 1
    dimens1D = iData

    call h5screate_simple_f(rank, dimens1D, dataspace, error) 
    call h5dcreate_f(fileID, "Real Metadata",&
         H5T_NATIVE_DOUBLE,dataspace, dataset, error)


    if (iProc == 0)&
         call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, &
         RealMetaData,dimens1D, error, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
    if (error == -1)&
         call stop_mpi(&
         "error in subroutine init_grid. Error marker 3")


    dimens1D = 1    
    call h5screate_simple_f(1, dimens1D, attributeSpace, error)

    !   Include a code version attribute in case anyone wants to know what
    !   made the file
    call h5acreate_f(dataset, "Code Version", H5T_NATIVE_DOUBLE,&
         attributeSpace, attribute,error, H5P_DEFAULT_F)
    call h5awrite_f(attribute, H5T_NATIVE_DOUBLE, CodeVersion,&
         dimens1D, error)
    call h5sclose_f(attributeSpace, error)
    call h5aclose_f(attribute, error)

    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
  end subroutine write_real_metadata

  !=================================================================
  !===============================================================``=
  subroutine write_integer_metadata(fileID, nPlotVar)

    use ModProcMH, only : nProc, iProc
    use BATL_lib, only : nDim, nDimAmr, nLevel,IsCartesian, IsRzGeometry,&
         IsSpherical,IsCylindrical,IsPeriodic_D  
    use ModMain, only : n_step, nI, nJ, nK
    integer (HID_T), intent(in) :: fileID, nPlotVar
    integer(HID_T) :: dataset, dataspace
    integer(HSIZE_T) :: dimens1D(1)

    integer :: error, rank, i, iData
    integer :: IntegerMetaData(14)
    integer :: FFV = 1 ! file format version
    iData = 1

    IntegerMetaData(iData) = FFV
    !    attName(1) = 'File Format Version'
    iData = iData + 1
    IntegerMetaData(iData) = n_step 
    !    attName(2) = "Time Step"
    iData = iData + 1
    IntegerMetaData(iData) = nDim
    !    attName(3) = 'nDim'
    iData = iData + 1
    IntegerMetaData(iData) = nBlkUsedGlobal
    !    attName(4) = 'globalNumBlocks'
    iData = iData + 1    
    IntegerMetaData(iData) = nProc
    !    attName(5) = 'numProcessors'
    iData = iData + 1

    IntegerMetaData(iData) = nI*nJ*nK
    !    attName(6) = 'nBLKcells'
    iData = iData + 1

    IntegerMetaData(iData) = nLevel
    !    attName(7) = 'nLevel'
    iData = iData + 1
    IntegerMetaData(iData) = nI
    !    attName(8) = 'nxb' 
    iData = iData + 1
    IntegerMetaData(iData) = nJ
    !    attName(9) = 'nyb'
    iData = iData + 1
    IntegerMetaData(iData) = nK
    !    attName(10) = 'nzb'

    !HDF5 has no boolean datatype so any logical parameters should go here.

    iData = iData + 1
    if(IsCartesian) then
       IntegerMetaData(iData) = 0
    elseif(IsRzGeometry) then
       IntegerMetaData(iData) = 1
    elseif(IsSpherical) then
       IntegerMetaData(iData) = 2
    elseif(IsCylindrical) then
       IntegerMetaData(iData) = 3
    endif

    !as of 2/3/2012 this is not implimented in the plugin but it probably
    !should be in the future
    do i = 1, 3
       if(IsPeriodic_D(i)) then
          IntegerMetaData(iData+i) = 1
       else
          IntegerMetaData(iData+i) = 0
       endif
    end do
    iData = iData + 3



    !-------------------------------------------------------------------
    !write the integer Metadata
    rank = 1
    dimens1D = iData

    call h5screate_simple_f(rank, dimens1D, dataspace, error) 
    call h5dcreate_f(fileID, "Integer Metadata",&
         H5T_NATIVE_INTEGER,dataspace, dataset, error)

    if (iProc == 0)&   
         call h5dwrite_f(dataset,H5T_NATIVE_INTEGER,IntegerMetaData,dimens1D,&
         error, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
    if (error == -1)&
         call stop_mpi(&
         "error in subroutine init_grid. Error marker 3")

    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)
  end subroutine write_integer_metadata

end module ModHdf5
