module ModHdf5Utils

  use hdf5
  use ModMpiOrig

  implicit none
  private!except
  save

  public:: losHdf5PlotWrite
  public:: writeHdf5Rank2Real
  public:: writeHdf5Rank4Real
  public:: writeHdf5Rank1Integer
  public:: writeHdf5Rank2Integer
  public:: writeHdf5Rank3Real
  public:: hdf5_init_file
  public:: write_plot_string
  public:: write_hdf5_attribute
contains 
  subroutine losHdf5PlotWrite(fileName, TypePosition, TypeStatus, StringHeader,&
       nStep, nBlocksUsed,nBlkUsedGlobal, offset,Time, nDimOut, nParam, nVar, n_D,  NameVar, PlotVarNameList,&
       unitList,coordMin, coordMax,iCoordMin, bBox, plotVarBlk, codeVersion, mpiComm, iProc, nProc)


    character (len=*), intent(in) :: fileName
    character (len=*), intent(in) :: TypePosition
    character (len=*), intent(in) :: PlotVarNameList(:)
    character (len=*), intent(in) :: unitList(:)
    character (len=10), intent(in) :: TypeStatus
    character (len=500), intent(in) :: StringHeader
    character (len=500), intent(in) :: NameVar
    character (len=501) :: HeaderString

    integer, intent(in) :: nDimOut, nParam, nVar, n_D(3), mpiComm, nStep, nBlocksUsed,nBlkUsedGlobal, offset
    integer, intent(in) :: iCoordMin(:,:), iProc, nProc
    real, intent(in) :: Time, coordMin(:), coordMax(:), codeVersion
    real, intent(in) ::plotVarBlk(:,:,:,:,:), bbox(:,:,:)

    integer (HID_T) :: fileID
    integer(HID_T) :: dataset, dataspace
    integer(HSIZE_T) :: dimens1D(1)
    integer :: error, rank, iData, IntegerMetaData(16), iLen, iVar
    integer :: labelLeng, numCells,i,j,k,n, iBlk
    integer, parameter :: FFV = 1

    real :: RealMetaData(7),varMin, varMax
    real, allocatable :: coordinates(:,:)
    integer, allocatable :: procNumAndLevel(:)
    integer, parameter :: lnameh5 = 11
    character (len=lnameh5) :: UnknownNameArray(nVar)
    logical, parameter :: collectiveWrite = .true.
    ! 
    !     write(*,*)  filename, TypePosition, TypeStatus
    !     write(*,*)  StringHeader
    !     write(*,*)  NameVar
    !     write(*,*)  "nStep=",nStep, "Time=", Time, "nDimOut=", nDimOut, "nParam=",nParam, "nVar=",nVar, "n_D=", n_D,&
    !                 "PlotVarNameList="
    !     do iVar = 1, nVar
    !         write(*,*) (PlotVarNameList(iVar))
    !     end do
    !   
    !     
    call h5open_f(error)
    fileID = -1
    call hdf5_init_file(fileID, fileName, mpiComm)
    if(fileID == -1) then
       write (*,*)  "Error: unable to initialize file"
       call CON_stop("unable to initialize hdf5 file")
    end if
    do iVar = 1, nVar
       UnknownNameArray(iVar) = trim(PlotVarNameList(iVar))
       !The VisIt plugin needs null padded names.
       labelLeng = len_trim(UnknownNameArray(iVar))
       do iLen = labelLeng + 1,lNameH5 
          UnknownNameArray(iVar)(iLen:iLen) = CHAR(0)
       end do
    end do
    call write_plot_string(nVar,lNameH5, UnknownNameArray,"plotVarNames",&
         fileID)
    do iVar = 1, nVar
       varMin = minVal(plotVarBlk(:,:,:,:,iVar))
       varMax = maxVal(plotVarBlk(:,:,:,:,iVar))
       call writeHdf5Rank4Real(fileID, plotVarBlk(:,:,:,:,iVar), nBlocksUsed,&
            nBlkUsedGlobal, offset, UnknownNameArray(iVar), &
            n_D(1), n_D(2), n_D(3), .true., varMin, varMax,&
            collectiveWrite)
    end do
    call writeHdf5Rank3Real(fileID, bBox, nBlocksUsed,nBlkUsedGlobal,&
         offset, "bounding box", 2, nDimOut,collectiveWrite)
    allocate(coordinates(nDimOut, nBlocksUsed))
    do iBlk=1,nBlocksUsed
       coordinates(1:nDimOut, iBlk) = .5*(bbox(1,1:nDimOut, iBlk) + bbox(2,1:nDimOut, iBlk))
    end do
    call writeHdf5Rank2Real(fileID, coordinates, nBlocksUsed,&
         nBlkUsedGlobal, offset, "coordinates", nDimOut, collectiveWrite)   

    deallocate(coordinates)

    ! 
    !     call writeHdf5Rank2Real(fileID, Coord_ID, nCells,&
    !          nCells, 0, "cell coordinates", nDimOut, collectiveWrite)
    !     
    do iVar = 1, nVar
       UnknownNameArray(iVar) = trim(unitList(iVar))
       !The VisIt plugin needs null padded names.
       labelLeng = len_trim(UnknownNameArray(iVar))
       do iLen = labelLeng + 1,lNameH5 
          UnknownNameArray(iVar)(iLen:iLen) = CHAR(0)
       end do
    end do
    call write_plot_string(nVar,lNameH5, UnknownNameArray,"plotVarUnits",&
         fileID)
    UnknownNameArray(1) = "X-Axis"
    UnknownNameArray(2) = "Y-Axis"
    UnknownNameArray(3) = "Z-Axis"
    do iVar = 1, nDimOut
       labelLeng = len_trim(UnknownNameArray(iVar))
       do iLen = labelLeng + 1,lNameH5 
          UnknownNameArray(iVar)(iLen:iLen) = CHAR(0)
       end do
    end do
    call write_plot_string(nDimOut, lNameH5, UnknownNameArray(1:nDimOut), "Axis Labels",&
         fileID)
    HeaderString = trim(StringHeader)
    labelLeng = len_trim(StringHeader)
    iLen = labelLeng + 1
    HeaderString(iLen:iLen) = CHAR(0)
    call write_plot_string(1, iLen, HeaderString, "Header",fileID)

    call writeHdf5Rank2Integer(fileID, iCoordMin, nBlocksUsed,&
         nBlkUsedGlobal, offset, "MinLogicalExtents", nDimOut, collectiveWrite)


    allocate(procnumAndLevel(nBlocksUsed))
    procnumAndLevel = 1
    call writeHdf5Rank1Integer(fileID, procnumAndLevel, &
         nBlocksUsed, nBlkUsedGlobal, offset,"refine level", collectiveWrite)
    procNumAndLevel = iProc
    call writeHdf5Rank1Integer(fileID, procNumAndLevel,nBlocksUsed, &
         nBlkUsedGlobal, offset,"Processor Number", collectiveWrite)
    deallocate(procnumAndLevel)
    !    allocate(attName(nAtts))
    iData = 1
    !    attName(1) = 'Simulation Time'
    RealMetaData(iData) = Time
    iData = iData + 1
    RealMetaData(iData) = coordMin(1)
    iData = iData + 1
    RealMetaData(iData) = coordMax(1)
    iData = iData + 1
    RealMetaData(iData) = coordMin(2)
    iData = iData + 1
    RealMetaData(iData) = coordMax(2)
    iData = iData + 1
    if(nDimOut < 3) then
       RealMetaData(iData) = 0!coordMin(3)
       iData = iData + 1
       RealMetaData(iData) = 0!coordMax(3)   
    else
       RealMetaData(iData) = coordMin(3)
       iData = iData + 1
       RealMetaData(iData) = coordMax(3)
    end if

    !-------------------------------------------------------------------
    !write the real Metadata
    rank = 1
    dimens1D = iData

    call h5screate_simple_f(rank, dimens1D, dataspace, error) 
    call h5dcreate_f(fileID, "Real Plot Metadata",&
         H5T_NATIVE_DOUBLE,dataspace, dataset, error)


    call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, &
         RealMetaData,dimens1D, error, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
    if (error == -1)&
         call CON_stop(&
         "error in subroutine init_grid. Error marker 3")
    call write_hdf5_attribute("Code Version", dataset, H5T_NATIVE_DOUBLE,&
         realAtt=CodeVersion)

    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)




    iData = 1

    IntegerMetaData(iData) = FFV
    !    attName(1) = 'File Format Version'
    iData = iData + 1
    IntegerMetaData(iData) = nStep 
    !    attName(2) = "Time Step"
    iData = iData + 1
    IntegerMetaData(iData) = nDimOut
    !    attName(3) = 'nDim'
    iData = iData + 1
    IntegerMetaData(iData) = 0
    !    attName(3) = 'nDimAMR'
    iData = iData + 1
    IntegerMetaData(iData) = nBlkUsedGlobal
    !    attName(4) = 'globalNumBlocks'
    iData = iData + 1    
    IntegerMetaData(iData) = nProc
    !    attName(5) = 'numProcessors'
    iData = iData + 1

    IntegerMetaData(iData) = 1 ! no AMR levels
    !    attName(7) = 'nLevel'
    iData = iData + 1
    IntegerMetaData(iData) = n_D(1)
    iData = iData + 1
    IntegerMetaData(iData) = n_D(2)
    iData = iData + 1
    IntegerMetaData(iData) = n_D(3)
    iData = iData + 1


    IntegerMetaData(iData) = 0 !I think these are always Cartesian. Refer to geometry types in
    !GM/BATSRUS/src/ModHdf5_orig for geometry types.
    iData = iData + 1
    !as of 2/3/2012 this is not implimented in the plugin but it probably
    !should be in the future
    do i = 1, 3
       IntegerMetaData(iData) = 0 ! none of the axis are periodic
       iData = iData + 1
    end do

    integerMetaData(iData) = 1 !Tells VisIt that this is a cut file, which to VisIt means that
    ! there is no morton curve index to be read.
    iData = iData + 1
    integerMetaData(iData) = nVar

    !-------------------------------------------------------------------
    !write the integer Metadata
    rank = 1
    dimens1D = iData

    call h5screate_simple_f(rank, dimens1D, dataspace, error) 
    call h5dcreate_f(fileID, "Integer Plot Metadata",&
         H5T_NATIVE_INTEGER,dataspace, dataset, error)

    !    if (iProc == 0)&   
    call h5dwrite_f(dataset,H5T_NATIVE_INTEGER,IntegerMetaData,dimens1D,&
         error, H5S_ALL_F, H5S_ALL_F, H5P_DEFAULT_F)
    if (error == -1)&
         call CON_stop(&
         "error in subroutine init_grid. Error marker 3")
    write (*,*) "closing file"
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)


    call h5garbage_collect_f(error)
    call h5fclose_f(fileID,error)
    !closing the hdf5 interface
    call h5close_f(error)
    if (error == -1) write (*,*) 'h5fclose_f failed'



  end subroutine losHdf5PlotWrite
  !These are the same write routines that used to live in GM/BATSRUS/src/ModHdf5.  I put them here so write_plot_file
  !could access them.
  subroutine writeHdf5Rank3Real(fileID, dataBuff, localNumBlocks,&
       globalNumBlocks, localOffset, description, nIplot, nJplot, collectiveWrite)

    implicit none

    integer :: rank
    integer, intent(in) :: nIplot,nJplot
    integer :: error
    integer(HID_T) :: dataset, plist_id
    integer(HSIZE_T) :: dimens3D(3)
    integer(HID_T) :: dataspace
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: start3D(3),stride3D(3),count3D(3)
    integer(HID_T), intent(in) :: fileID
    integer, intent(in) :: localNumBlocks, globalNumBlocks
    real, intent(in) :: dataBuff(:,:,:)
    integer, intent(in) :: localOffset
    character (len=*), intent(in) :: description
    logical, intent(in) :: collectiveWrite
    integer(HSIZE_T) :: one

    !Set the dimensions of the dataset

    rank = 3
    dimens3D(1) = nIplot
    dimens3D(2) = nJplot
    dimens3D(3) = globalNumBlocks
    call h5open_f(error)
    if (error == -1) &
         call CON_stop("error in subroutine writeHdf5Rank3Real. Error marker 1")

    !     if (localNumBlocks == 0 .and. (.not. collectiveWrite)) return
    call h5screate_simple_f(rank, dimens3D, dataspace, error) 
    !create the dataset
    call h5dcreate_f(fileID, description, H5T_NATIVE_DOUBLE, dataspace, dataset, error)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (collectiveWrite) then
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    else
       !          call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    end if

    if (error == -1) &
         call CON_stop("error in subroutine writeHdf5Rank3Real. Error marker 2")

    if (localNumBlocks == 0) then
       if (collectiveWrite) then
          one = 0
          call h5screate_simple_f(1, (/one/), memspace, error)
          call h5sselect_none_f(memspace,error)
          call h5sselect_none_f(dataspace,error)
          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dataBuff, dimens3D, error, &
               mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
       else
          one = 1
          call h5screate_simple_f(1, (/one/), memspace, error)
       end if
    else

       start3D(1:2) = 0
       start3D(3) = localOffset 
       stride3D = 1
       count3D(1:2) = dimens3D(1:2)
       count3D(3) = localNumBlocks


       !create the hyperslab.  This will differ on the different processors
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start3D, count3D, error, &
            stride = stride3D)
       !         !create the memory space
       if (error == -1)& 
            call CON_stop("error in subroutine writeHdf5Rank3Real. Error marker 3")

       dimens3D(3) = localNumBlocks
       call h5screate_simple_f(rank, dimens3D, memspace, error)
       !Write the data
       if (error == -1) &
            call CON_stop("error in subroutine writeHdf5Rank3Real. Error marker 3")

       call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dataBuff, dimens3D, error, &
            mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
       if (error == -1) &
            call CON_stop("error in subroutine writeHdf5Rank3Real. Error marker 5")
    end if
    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)

    if (error == -1)& 
         call CON_stop("error in subroutine writeHdf5Rank3Real. Error marker 6")

  end subroutine writeHdf5Rank3Real



  !======================================================================
  !=====================================================================
  subroutine writeHdf5Rank1Integer(fileID, dataBuff, localNumBlocks,&
       globalNumBlocks, localOffset, description, collectiveWrite)

    implicit none

    logical, intent(in) :: collectiveWrite
    integer :: rank
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
    integer(HSIZE_T) ::  one
    !Set the dimensions of the dataset
    rank = 1
    dimens1D(1) = globalNumBlocks
    call h5open_f(error)
    if (error == -1) &
         call CON_stop("error in subroutine writeHdf5Rank1Integer. Error marker 1")

    call h5screate_simple_f(rank, dimens1D, dataspace, error) 
    !create the dataset
    call h5dcreate_f(fileID, description, H5T_NATIVE_INTEGER, dataspace, dataset, error)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (collectiveWrite) then
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    else
       !nothing = INDEPENDENT  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    end if


    if (localNumBlocks == 0) then
       if (collectiveWrite) then
          one = 0
          call h5screate_simple_f(1, (/one/), memspace, error)
          call h5sselect_none_f(memspace,error)
          call h5sselect_none_f(dataspace,error)
          call h5dwrite_f(dataset, H5T_NATIVE_INTEGER, dataBuff, dimens1D, error, &
               mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
       else
          one = 1
          call h5screate_simple_f(1, (/one/), memspace, error)
       end if
    else
       if (error == -1) &
            call CON_stop("error in subroutine writeHdf5Rank1Integer. Error marker 2")

       start1D(1) = localOffset 
       stride1D(:) = 1
       count1D(1) = localNumBlocks

       !create the hyperslab.  This will differ on the different processors
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start1D, count1D, error, &
            stride = stride1D)
       !         !create the memory space
       if (error == -1)& 
            call CON_stop("error in subroutine writeHdf5Rank1Integer. Error marker 3")

       dimens1D(1) = localNumBlocks
       call h5screate_simple_f(rank, dimens1D, memspace, error)
       !Write the data
       if (error == -1) &
            call CON_stop("error in subroutine writeHdf5Rank1Integer. Error marker 4")

       call h5dwrite_f(dataset, H5T_NATIVE_INTEGER, dataBuff, dimens1D, error, &
            mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)


       if (error == -1) &
            call CON_stop("error in subroutine writeHdf5Rank1Integer. Error marker 5")
    end if
    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)

    if (error == -1)& 
         call CON_stop("error in subroutine writeHdf5Rank1Integer. Error marker 6")

  end subroutine writeHdf5Rank1Integer

  !======================================================================
  !=====================================================================


  subroutine writeHdf5Rank1Real(fileID, dataBuff, localNumBlocks,&
       globalNumBlocks, localOffset, description, hasMinMax, varMin, varMax, collectiveWrite)

    implicit none

    logical, intent(in) :: collectiveWrite
    integer :: rank
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
    real, intent(in) :: dataBuff(:), varMin, varMax
    logical, intent(in) :: hasMinMax
    integer, intent(in) :: localOffset
    character (len=*), intent(in) :: description
    integer(HSIZE_T) ::  one
    !Set the dimensions of the dataset
    rank = 1
    dimens1D(1) = globalNumBlocks
    call h5open_f(error)
    if (error == -1) &
         call CON_stop("error in subroutine writeHdf5Rank1Real. Error marker 1")

    call h5screate_simple_f(rank, dimens1D, dataspace, error) 
    !create the dataset
    call h5dcreate_f(fileID, description, H5T_NATIVE_DOUBLE, dataspace, dataset, error)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (collectiveWrite) then
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    else
       !nothing = INDEPENDENT  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    end if

    if (hasMinMax) then
       !Minimum
       call write_hdf5_attribute("minimum", dataset, H5T_NATIVE_DOUBLE, realAtt = varMin)
       call write_hdf5_attribute("maximum", dataset, H5T_NATIVE_DOUBLE, realAtt = varMax)
    end if



    if (localNumBlocks == 0) then
       if (collectiveWrite) then
          one = 0
          call h5screate_simple_f(1, (/one/), memspace, error)
          call h5sselect_none_f(memspace,error)
          call h5sselect_none_f(dataspace,error)
          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dataBuff, dimens1D, error, &
               mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
       else
          one = 1
          call h5screate_simple_f(1, (/one/), memspace, error)
       end if
    else
       if (error == -1) &
            call CON_stop("error in subroutine writeHdf5Rank1Real. Error marker 2")

       start1D(1) = localOffset 
       stride1D(:) = 1
       count1D(1) = localNumBlocks

       !create the hyperslab.  This will differ on the different processors
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start1D, count1D, error, &
            stride = stride1D)
       !         !create the memory space
       if (error == -1)& 
            call CON_stop("error in subroutine writeHdf5Rank1Real. Error marker 3")

       dimens1D(1) = localNumBlocks
       call h5screate_simple_f(rank, dimens1D, memspace, error)
       !Write the data
       if (error == -1) &
            call CON_stop("error in subroutine writeHdf5Rank1Real. Error marker 4")

       call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dataBuff, dimens1D, error, &
            mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)


       if (error == -1) &
            call CON_stop("error in subroutine writeHdf5Rank1Real. Error marker 5")
    end if
    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)

    if (error == -1)& 
         call CON_stop("error in subroutine writeHdf5Rank1Real. Error marker 6")

  end subroutine writeHdf5Rank1Real

  !=====================================================================
  !=====================================================================

  subroutine writeHdf5Rank2Integer(fileID, dataBuff, localNumBlocks,&
       globalNumBlocks, localOffset, description, dimens, collectiveWrite)

    implicit none

    logical, intent(in) :: collectiveWrite
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
    integer(HSIZE_T) :: one

    !Set the dimensions of the dataset
    rank = 2
    dimens2D(1) = dimens
    dimens2D(2) = globalNumBlocks
    call h5open_f(error)
    if (error == -1) &
         
         call CON_stop("error in subroutine writeHdf5Rank2Integer. Error marker 1")

    !     if (localNumBlocks == 0 .and. (.not. collectiveWrite)) return
    call h5screate_simple_f(rank, dimens2D, dataspace, error) 
    !create the dataset
    call h5dcreate_f(fileID, description, H5T_NATIVE_INTEGER, dataspace, dataset, error)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (collectiveWrite) then
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    else
       !         call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    end if


    if (localNumBlocks == 0) then
       if (collectiveWrite) then
          one = 0
          call h5screate_simple_f(1, (/one/), memspace, error)
          call h5sselect_none_f(memspace,error)
          call h5sselect_none_f(dataspace,error)
          call h5dwrite_f(dataset, H5T_NATIVE_INTEGER, dataBuff, dimens2D, error, &
               mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
       else
          one = 1
          call h5screate_simple_f(1, (/one/), memspace, error)
       end if
    else
       if (error == -1) &
            call CON_stop("error in subroutine writeHdf5Rank2Integer. Error marker 2")

       start2d(1) = 0
       start2D(2) = localOffset 
       stride2D(:) = 1
       count2D(1) = dimens2D(1)
       count2D(2) = localNumBlocks

       !create the hyperslab.  This will differ on the different processors
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start2D, count2D, error, &
            stride = stride2D)
       !         !create the memory space
       if (error == -1)& 
            call CON_stop("error in subroutine writeHdf5Rank2Integer. Error marker 3")

       dimens2D(2) = localNumBlocks
       call h5screate_simple_f(rank, dimens2D, memspace, error)
       !Write the data
       if (error == -1) &
            call CON_stop("error in subroutine writeHdf5Rank2Integer. Error marker 4")

       call h5dwrite_f(dataset, H5T_NATIVE_INTEGER, dataBuff, dimens2D, error, &
            mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
       if (error == -1) &
            call CON_stop("error in subroutine writeHdf5Rank2Integer. Error marker 5")
    end if

    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)

    if (error == -1)& 
         call CON_stop("error in subroutine writeHdf5Rank2Integer. Error marker 6")

  end subroutine writeHdf5Rank2Integer

  !======================================================================
  !=====================================================================

  subroutine writeHdf5Rank2Real(fileID, dataBuff, localNumBlocks,&
       globalNumBlocks, localOffset, description, dimens, collectiveWrite)

    implicit none

    logical, intent(in) :: collectiveWrite
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
    integer(HSIZE_T) :: one, nullSize(1) = 1

    !Set the dimensions of the dataset
    rank = 2
    dimens2D(1) = dimens
    dimens2D(2) = globalNumBlocks
    call h5open_f(error)
    if (error == -1) &
         
         call CON_stop("error in subroutine writeHdf5Rank2Real. Error marker 1")
    !     if (localNumBlocks == 0 .and. (.not. collectiveWrite)) return

    call h5screate_simple_f(rank, dimens2D, dataspace, error) 
    !create the dataset
    call h5dcreate_f(fileID, description, H5T_NATIVE_DOUBLE, dataspace, dataset, error)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (collectiveWrite) then
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    else
       !         call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    end if


    if (localNumBlocks == 0) then
       if (collectiveWrite) then
          one = 0
          call h5screate_simple_f(1, (/one/), memspace, error)
          call h5sselect_none_f(memspace,error)
          call h5sselect_none_f(dataspace,error)
          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dataBuff, dimens2D, error, &
               mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
       else
          one = 1
          call h5screate_simple_f(1, (/one/), memspace, error)
       end if
    else
       if (error == -1) &
            call CON_stop("error in subroutine writeHdf5Rank2Real. Error marker 2")

       start2d(1) = 0
       start2D(2) = localOffset 
       stride2D(:) = 1
       count2D(1) = dimens2D(1)
       count2D(2) = localNumBlocks

       !create the hyperslab.  This will differ on the different processors
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start2D, count2D, error, &
            stride = stride2D)
       !         !create the memory space
       if (error == -1)& 
            call CON_stop("error in subroutine writeHdf5Rank2Real. Error marker 3")

       dimens2D(2) = localNumBlocks
       call h5screate_simple_f(rank, dimens2D, memspace, error)
       !Write the data
       if (error == -1) &
            call CON_stop("error in subroutine writeHdf5Rank2Real. Error marker 4")

       call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dataBuff, dimens2D, error, &
            mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
       if (error == -1) &
            call CON_stop("error in subroutine writeHdf5Rank2Real. Error marker 5")
    end if

    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)

    if (error == -1)& 
         call CON_stop("error in subroutine writeHdf5Rank2Real. Error marker 6")

  end subroutine writeHdf5Rank2Real

  !======================================================================
  !=====================================================================
  subroutine hdf5_init_file(fileID, filename, mpiComm)

    character (len=80), intent(in) :: filename
    integer :: error, accTemplate
    integer, optional, intent(in) :: mpiComm
    integer, intent(inout) :: fileID
    INTEGER :: mpierror       ! MPI error flag
    INTEGER :: dxpList
!!!    INTEGER :: mpi_size, mpi_rank, comm

    !Initalize currentBlock andplotArrayBegin for this file

    !    call h5open_f(error)                    

    !create MPI info Object

    !Create file access propertty list
    call h5pcreate_f(H5P_FILE_ACCESS_F, accTemplate, error)
    if(present(mpiComm))&
         CALL h5pset_fapl_mpio_f(accTemplate, mpiComm, MPI_INFO_NULL, error)
    ! Create the file collectively.

    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, fileID, error,&
         access_prp = accTemplate)
    CALL h5pclose_f(accTemplate, error)
    if (error == -1) &
         call CON_stop(&
         "error in subroutine hdf5_init_file. Error marker 1")

    if (error == -1) fileID = -1
  end subroutine hdf5_init_file

  !=====================================================================
  !=====================================================================

  subroutine write_plot_string(nPlotVar, stringLen, dataBuff, description, fileID)

    integer, intent(in) :: nPlotVar, stringLen
    character (len=*), intent(in) :: dataBuff(nPlotVar)
    integer (HID_T), intent(in) :: fileID
    character (len=*), intent(in) :: description
    integer (HID_T) :: dataType, dataset, dataspace
    integer (HSIZE_T) :: sizeString, dimens1D(1) 

    integer :: error, rank

    rank = 1  
    dimens1d(1) = nPlotVar

    sizeString = stringLen


    call h5tcopy_f(H5T_NATIVE_CHARACTER, datatype, error)
    call h5tset_size_f(datatype, sizeString, error)
    if (error == -1) &
         call CON_stop(&
         "error in subroutine writeHdf5Header. Error marker 10")
    call h5screate_simple_f(rank, dimens1D, dataspace, error)
    call h5dcreate_f(&
         fileID, description, datatype, dataspace, dataset, error)
    if (error == -1) &
         call CON_stop("error in subroutine writeHdf5Header. Error marker 11")

    call h5dwrite_f(dataset, datatype, databuff, dimens1D, error, H5S_ALL_F,&
         H5S_ALL_F, H5P_DEFAULT_F)
    if (error == -1) &
         call CON_stop("error in subroutine writeHdf5Header. Error marker 11")



    call h5tclose_f(datatype, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)

    if (error == -1)& 
         call CON_stop(&
         "error in subroutine writeHdf5Header. Error marker 13")
  end subroutine write_plot_string

  !=====================================================================
  !=====================================================================

  subroutine writeHdf5Rank4Real(fileID, dataBuff, localNumBlocks,&
       globalNumBlocks, localOffset, description, nIplot, nJplot, nKplot,&
       minMax ,vMin, vMax, collectiveWrite)

    implicit none

    logical, intent(in) :: collectiveWrite
    integer :: rank
    integer, intent(in) :: nIplot,nJplot,nKplot
    integer :: error
    integer(HID_T) :: dataset, plist_id
    integer(HSIZE_T) :: dimens4D(4), dimens1D(1)
    integer(HID_T) :: dataspace
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: start4D(4)
    integer(HSIZE_T) :: stride4D(4)
    integer(HSIZE_T) :: count4D(4)
    integer(HID_T), intent(in) :: fileID
    integer, intent(in) :: localNumBlocks, globalNumBlocks
    real, intent(in) :: dataBuff(:,:,:,:)
    integer, intent(in) :: localOffset
    character (len=*), intent(in) :: description
    logical, intent(in) :: minMax
    real, intent(in) :: vMin, vMax
    integer(HSIZE_T) :: one

    !Set the dimensions of the dataset

    rank = 4
    dimens4D(1) = nIplot
    dimens4D(2) = nJplot
    dimens4D(3) = nKplot
    dimens4D(4) = globalNumBlocks
    call h5open_f(error)
    if (error == -1) &
         call CON_stop("error in subroutine writeHdf5Rank4Real. Error marker 1")


    call h5screate_simple_f(rank, dimens4D, dataspace, error) 
    !create the dataset
    call h5dcreate_f(fileID, description, H5T_NATIVE_DOUBLE, dataspace, dataset, error)

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (collectiveWrite) then
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    else
       !          call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    end if

    if (error == -1) &
         call CON_stop("error in subroutine writeHdf5Rank4Real. Error marker 2")

    if (minMax) then
       call write_hdf5_attribute("minimum", dataset, H5T_NATIVE_DOUBLE, realAtt = vMin)
       call write_hdf5_attribute("maximum", dataset, H5T_NATIVE_DOUBLE, realAtt = vMax)
    end if

    if (localNumBlocks == 0) then
       if (collectiveWrite) then
          one = 0
          call h5screate_simple_f(1, (/one/), memspace, error)
          call h5sselect_none_f(memspace,error)
          call h5sselect_none_f(dataspace,error)
          call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dataBuff, dimens4D, error, &
               mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
       else
          one = 1
          call h5screate_simple_f(1, (/one/), memspace, error)
       end if
    else

       start4D(1:3) = 0
       start4D(4) = localOffset 
       stride4D(:) = 1
       count4D(1:3) = dimens4D(1:3)
       count4D(4) = localNumBlocks


       !create the hyperslab.  This will differ on the different processors
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start4D, count4D, error, &
            stride = stride4D)
       !         !create the memory space
       if (error == -1)& 
            call CON_stop("error in subroutine writeHdf5Rank4Real. Error marker 3")

       dimens4D(4) = localNumBlocks
       call h5screate_simple_f(rank, dimens4D, memspace, error)
       !Write the data
       if (error == -1) &
            call CON_stop("error in subroutine writeHdf5Rank4Real. Error marker 4")
       call h5dwrite_f(dataset, H5T_NATIVE_DOUBLE, dataBuff, dimens4D, error, &
            mem_space_id = memspace, file_space_id = dataspace, xfer_prp = plist_id)
       if (error == -1) &
            call CON_stop("error in subroutine writeHdf5Rank4Real. Error marker 5")
    end if

    call h5sclose_f(memspace,error)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dataset, error)

    if (error == -1)& 
         call CON_stop("error in subroutine writeHdf5Rank4Real. Error marker 6")

  end subroutine writeHdf5Rank4Real

  subroutine write_hdf5_attribute(attName, dataset, datatype, realAtt,&
       intAtt)
    integer(HID_T), intent(in) :: dataset, datatype
    character (len=*), intent(in) :: attName
    real, optional, intent(in) :: realAtt
    real, optional, intent(in) :: intAtt
    integer(HID_T) :: attributeSpace, attribute
    integer(HSIZE_T) :: dimens1D(1)
    integer :: error

    dimens1D=(1)
    call h5screate_simple_f(1, dimens1D, attributeSpace, error)

    call h5acreate_f(dataset, attname, datatype,&
         attributeSpace, attribute,&
         error, H5P_DEFAULT_F)
    if (Present(realAtt)) then
       call h5awrite_f(attribute, datatype, realAtt, dimens1D, error)
    elseif (present(intAtt)) then
       call h5awrite_f(attribute, datatype, realAtt, dimens1D, error)
    endif
    call h5sclose_f(attributeSpace, error)
    call h5aclose_f(attribute, error)

  end subroutine write_hdf5_attribute
end module ModHdf5Utils
