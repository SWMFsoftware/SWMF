!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModWriteTecplot

  ! Save cell centered data into Tecplot files
  !
  ! Save cell centers (same as 3D IDL plots, but no cell size).
  ! For single data file: 
  !     record size is nPlotVar*nIJK
  !     record index is iRec = iMortonNode_A(iNode)
  ! For single connectivity file:
  !     record size is 8 integers (global cell indexes of brick)
  !     record index is 
  !     iRec = iRec0_A(iNode) + i + (j-1)*nIPoint + (k-1)*nIJPoint
  !  where i goes from i0 to i1, j from j0 to j1, k from k0 to k1.
  !     and nIPoint = i1-i0+1, nIJPoint = nIPoints*(j1-j0+1)
  !     and total number of connectivities for a block is 
  !        nIJKPoint = nIJPoint*(k1-k0+1)
  !  Connectivity at block boundaries with no resolution change
  !     is written by the side towards the negative direction
  !  Connectivity at resolution changes is written by the finer side.
  !  So
  !      i0 =  0   if block on the left is coarser and 1 otherwise
  !      i1 = nI-1 if block on right is finer or no block, nI otherwise
  !      same for j0, j1, k0 and k1
  ! The connectivity brick contains i:i+1,j:j+1,k:k+1 cell centers.
  !      Indexes belonging to ghost cells are taken from the neighbor.
  !      These could be obtained by message passing the cell indexes (simple!)
  !      Or can be calculated. 
  !      Fine ghost cell is replaced with the index of the coarse cell contain

  ! Count number of used blocks for each processor and save number of 
  !    preceeding blocks: nBlockBefore_BP
  ! Count number of connectivity bricks to be written before this node: iRec0_A
  !  inside a cut region?
  !
  ! Possibly limit to a cut region (?) Things get complicated.

  implicit none

  SAVE

  private ! except

  public:: write_tecplot_data    ! write plot variables for 1 block
  public:: write_tecplot_head    ! write header information
  public:: write_tecplot_connect ! write connectivity file(s)
  public:: write_tecplot_setinfo ! set some variables to be written as AUX
  public:: write_tecplot_auxdata ! write auxiliary information

  character(len=23), public :: textDateTime
  character(len=22), public :: textNandT
  character, public, parameter:: CharNewLine = char(10)
  integer,   public, parameter:: lRecConnect = 8*11+1 ! for (8i11,\) format
  
  ! Local variables
  character (len=23) :: textDateTime0

  ! Total number of connectivity bricks
  integer:: nBrickAll

  integer:: iRecData = -1
  character(len=80) :: StringFormat

contains
  !===========================================================================
  subroutine write_tecplot_data(iBlock, nPlotVar, PlotVar_GV)

    use BATL_lib,  ONLY: nI, nJ, nK, nIJK, &
         MinI, MaxI, MinJ, MaxJ, MinK, MaxK, Xyz_DGB, iProc
    use ModIO,     ONLY: nPlotVarMax, DoSaveOneTecFile
    use ModIoUnit, ONLY: UnitTmp_
    use ModAdvance,ONLY: iTypeAdvance_BP, SkippedBlock_
    use ModMain,   ONLY: nBlockMax, BlkTest

    integer, intent(in):: iBlock, nPlotVar
    real,    intent(in):: PlotVar_GV(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nPlotvarMax)

    integer:: i, j, k

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'write_tecplot_data'
    !--------------------------------------------------------------------------
    if(iBlock==BlkTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if
    if(DoTestMe) write(*,*) NameSub,' starting with nPlotVar=', nPlotVar

    if(DoSaveOneTecFile)then
       if(iRecData < 0)then
          ! Initialize the record index based on number of used blocks 
          iRecData = 0
          if(iProc > 0) iRecData = nIJK* &
               count(iTypeAdvance_BP(1:nBlockMax,0:iProc-1) /= SkippedBlock_)
          write(StringFormat, '(a,i2,a)') "(", nPlotVar+3, "(ES14.6), a)"
       end if
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          iRecData = iRecData + 1
          write(UnitTmp_, StringFormat, REC=iRecData) &
               Xyz_DGB(1:3,i,j,k,iBlock),             &
               PlotVar_GV(i,j,k,1:nPlotVar),          &
               CharNewLine
       end do; end do; end do
    else
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          write(UnitTmp_,"(50(ES14.6))")            &
               Xyz_DGB(1:3,i,j,k,iBlock),           &
               PlotVar_GV(i,j,k,1:nPlotVar)
       end do; end do; end do
    end if

    if(DoTestMe) write(*,*) NameSub,' finished'

  end subroutine write_tecplot_data

  !========================================================================
  subroutine write_tecplot_connect(iFile, NameFile)
    
    use ModProcMH,    ONLY: iProc, nProc, iComm
    use ModAdvance,   ONLY: iTypeAdvance_BP, SkippedBlock_
    use ModIO,        ONLY: DoSaveOneTecFile
    use ModIoUnit,    ONLY: UnitTmp_
    use ModUtilities, ONLY: open_file, close_file
    use ModMain,      ONLY: nBlockMax
    use ModMpi,       ONLY: MPI_SUM, mpi_reduce_integer_scalar, &
         MPI_allgather, MPI_INTEGER
    use BATL_lib,     ONLY: nI, nJ, nK, nIJK, MaxBlock, nBlock, Unused_B, &
         DiLevelNei_IIIB, message_pass_cell, set_tree_periodic

    ! Write out connectivity file

    integer,          intent(in):: iFile     ! index of plot file
    character(len=*), intent(in):: NameFile  ! name of connectivity file

    integer:: nBlockBefore, iCell, iBlock, iError

    integer:: i0, i1, j0, j1, k0, k1, i, j, k

    ! This is a real array for sake of message passing only
    real, allocatable:: CellIndex_GB(:,:,:,:)
    integer, allocatable:: iCell_G(:,:,:)

    ! Number of bricks for this and all processors
    integer:: nBrick
    integer, allocatable:: nBrick_P(:)

    ! Multiple stages may be needed
    integer:: iStage, nStage

    ! Record index
    integer:: iRec

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'write_tecplot_connect'
    !------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    if(DoTestMe)write(*,*) NameSub,' starting with NameFile=', NameFile

    allocate( &
         CellIndex_GB(0:nI+1,0:nJ+1,0:nK+1,MaxBlock), &
         iCell_G(0:nI+1,0:nJ+1,0:nK+1))

    ! count number of cells written by processors before this one
    nBlockBefore = 0
    if(iProc > 0) nBlockBefore = &
         count(iTypeAdvance_BP(1:nBlockMax,0:iProc-1) /= SkippedBlock_)

    if(DoTestMe)write(*,*) NameSub,' nBlockBefore=', nBlockBefore

    ! initial cell index
    iCell = nIJK*nBlockBefore

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       CellIndex_GB(:,:,:,iBlock) = 0
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          iCell = iCell + 1
          CellIndex_GB(i,j,k,iBlock) = iCell
       end do; end do; end do
    end do

    ! Set the global cell indexes for the ghost cells. First order prolongation
    ! is used, so fine ghost cells are set to the index of the coarse cell.
    ! Note that the coarse ghost cells are not going to be used.
    call message_pass_cell(1, CellIndex_GB, nProlongOrderIn=1)

    ! switch off "fake" periodicity so there are no connections
    call set_tree_periodic(.false.)

    if(DoSaveOneTecFile)then
       nStage = 2
       nBrick = 0
       allocate(nBrick_P(0:nProc-1))
       call open_file(File=NameFile, ACCESS='direct', RECL=lRecConnect)
    else
       nStage = 1
       nBrickAll = 0
       ! Open connectivity file
       call open_file(File=NameFile)
    end if

    do iStage = 1, nStage
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          ! Get the index limits for the lower left corner
          ! of the connectivity bricks.

          ! Start connectivity brick from index 1 unless coarser left neighbor
          i0 = 1
          if(DiLevelNei_IIIB(-1,0,0,iBlock) == 1) i0 = 0
          j0 = 1
          if(DiLevelNei_IIIB(0,-1,0,iBlock) == 1) j0 = 0
          k0 = 1
          if(DiLevelNei_IIIB(0,0,-1,iBlock) == 1) k0 = 0

          ! Finish at nI unless there is no block on the right or it is finer
          i1 = nI
          if(DiLevelNei_IIIB(1,0,0,iBlock) < 0) i1 = nI-1
          j1 = nJ
          if(DiLevelNei_IIIB(0,1,0,iBlock) < 0) j1 = nJ-1
          k1 = nK
          if(DiLevelNei_IIIB(0,0,1,iBlock) < 0) k1 = nK-1

          if(iStage < nStage)then
             ! Simply count the number of bricks for this processor
             nBrick = nBrick + (i1-i0+1)*(j1-j0+1)*(k1-k0+1)
             CYCLE
          end if

          iCell_G = nint(CellIndex_GB(:,:,:,iBlock))

          if(DoSaveOneTecFile)then

             do i = i0, i1; do j = j0, j1; do k = k0, k1
                iRec      = iRec + 1
                write(UnitTmp_,'(8i11,a)', REC=iRec) &
                     iCell_G(i  ,j  ,k  ), &
                     iCell_G(i+1,j  ,k  ), &
                     iCell_G(i+1,j+1,k  ), &
                     iCell_G(i  ,j+1,k  ), &
                     iCell_G(i  ,j  ,k+1), &
                     iCell_G(i+1,j  ,k+1), &
                     iCell_G(i+1,j+1,k+1), &
                     iCell_G(i  ,j+1,k+1), &
                     CharNewLine
             end do; end do; end do
          else
             do i = i0, i1; do j = j0, j1; do k = k0, k1
                nBrickAll = nBrickAll + 1
                write(UnitTmp_,'(8i11)') &
                     iCell_G(i  ,j  ,k  ), &
                     iCell_G(i+1,j  ,k  ), &
                     iCell_G(i+1,j+1,k  ), &
                     iCell_G(i  ,j+1,k  ), &
                     iCell_G(i  ,j  ,k+1), &
                     iCell_G(i+1,j  ,k+1), &
                     iCell_G(i+1,j+1,k+1), &
                     iCell_G(i  ,j+1,k+1)
             end do; end do; end do
          end if 
       end do ! iBlock
       if(iStage < nStage)then
          ! Collect number of bricks from all processors
          call MPI_allgather(nBrick, 1, MPI_INTEGER, &
               nBrick_P, 1, MPI_INTEGER, iComm, iError)
          ! Add up number of bricks on previous prorecssors
          iRec = 0
          if(iProc > 0) iRec = sum(nBrick_P(0:iProc-1))
       end if
    end do ! iStage
    call close_file

    ! Reset periodicity as it was
    call set_tree_periodic(.true.)

    ! Reset iRec for next plot
    iRecData = -1

    ! Calculate total number of bricks
    if(DoSaveOneTecFile)then
       if(iProc == 0) nBrickAll = sum(nBrick_P)
       deallocate(nBrick_P)
    elseif(nProc > 1)then
       call mpi_reduce_integer_scalar(nBrickAll, MPI_SUM, 0, iComm, iError)
    end if

    deallocate(CellIndex_GB, iCell_G)

    if(DoTestMe)write(*,*) NameSub,' done with nBrickAll=', nBrickAll

  end subroutine write_tecplot_connect

  !========================================================================
  subroutine write_tecplot_head(NameFile, StringUnit)

    use ModProcMH,    ONLY: iProc
    use ModIoUnit,    ONLY: UnitTmp_
    use ModUtilities, ONLY: open_file, close_file
    use BATL_lib,     ONLY: nNodeUsed, nIJK

    ! Write out tecplot header file

    character(len=*), intent(in):: NameFile
    character(len=*), intent(in):: StringUnit

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'write_tecplot_connect'
    !---------------------------------------------------------------------
    call write_tecplot_setinfo
    if(iProc /= 0) RETURN

    call open_file(File=NameFile)
    write(UnitTmp_,'(a)')'TITLE="BATSRUS: 3D Data, '//textDateTime//'"'
    write(UnitTmp_,'(a)') trim(StringUnit)
    write(UnitTmp_,'(a,a,i12,a,i12,a)') &
         'ZONE T="3D   '//textNandT//'"', &
         ', N=',nNodeUsed*nIJK, &
         ', E=',nBrickAll,      &
         ', F=FEPOINT, ET=BRICK'
    call write_tecplot_auxdata
    call close_file
    
  end subroutine write_tecplot_head

  !=========================================================================

  subroutine write_tecplot_auxdata(iUnitIn)

    use ModMain, ONLY : nI,nJ,nK, &
         nBlockALL, time_accurate,n_step, &
         nOrder, UseRotatingBc,           &
         TypeCoordSystem, CodeVersion, nTrueCellsALL
    use ModFaceValue, ONLY: TypeLimiter, BetaLimiter
    use ModMain, ONLY: boris_correction
    use ModPhysics, ONLY : &
         ThetaTilt, Rbody, boris_cLIGHT_factor, BodyNDim_I, Gamma_I
    use ModAdvance, ONLY : FluxType
    use ModMultiFluid, ONLY: IonFirst_
    use ModIoUnit, ONLY: UnitTmp_
    use ModIO, ONLY: StringDateOrTime
    use ModNumConst, ONLY : cRadToDeg
    use BATL_lib, ONLY: nProc, nIJK

    integer, intent(in), optional :: iUnitIn

    character(len=8)  :: real_date
    character(len=10) :: real_time
    integer :: iUnitHere
    character(len=500) :: stmp

    !-------------------------------------------------------------------------
    if (present(iUnitIn)) then
       iUnitHere = iUnitIn
    else
       iUnitHere = UnitTmp_
    end if

    !BLOCKS
    write(stmp,'(i12,3(a,i2))')nBlockALL,'  ',nI,' x',nJ,' x',nK
    write(iUnitHere,'(a,a,a)') 'AUXDATA BLOCKS="',trim(adjustl(stmp)),'"'

    !BODYDENSITY
    write(stmp,'(f12.2)')BodyNDim_I(IonFirst_)
    write(iUnitHere,'(a,a,a)') &
         'AUXDATA BODYNUMDENSITY="',trim(adjustl(stmp)),'"'

    !BORIS
    if(boris_correction)then
       write(stmp,'(a,f8.4)')'T ',boris_cLIGHT_factor
    else
       write(stmp,'(a)')'F'
    end if
    write(iUnitHere,'(a,a,a)') 'AUXDATA BORIS="',trim(adjustl(stmp)),'"'

    !BTHETATILT
    write(stmp,'(f12.4)')ThetaTilt*cRadToDeg
    write(iUnitHere,'(a,a,a)') 'AUXDATA BTHETATILT="',trim(adjustl(stmp)),'"'

    !CELLS
    write(stmp,'(i12)')nBlockALL*nIJK
    write(iUnitHere,'(a,a,a)') 'AUXDATA CELLS="',trim(adjustl(stmp)),'"'

    !CELLSUSED
    write(stmp,'(i12)')nTrueCellsALL
    write(iUnitHere,'(a,a,a)') 'AUXDATA CELLSUSED="',trim(adjustl(stmp)),'"'

    !CODEVERSION
    write(stmp,'(a,f5.2)')'BATSRUS',CodeVersion
    write(iUnitHere,'(a,a,a)') 'AUXDATA CODEVERSION="',trim(adjustl(stmp)),'"'

    !COORDSYSTEM
    write(stmp,'(a)')TypeCoordSystem
    write(iUnitHere,'(a,a,a)') 'AUXDATA COORDSYSTEM="',trim(adjustl(stmp)),'"'

    !COROTATION
    if(UseRotatingBc)then
       write(stmp,'(a)')'T'
    else
       write(stmp,'(a)')'F'
    end if
    write(iUnitHere,'(a,a,a)') 'AUXDATA COROTATION="',trim(adjustl(stmp)),'"'

    !FLUXTYPE
    write(stmp,'(a)')FluxType
    write(iUnitHere,'(a,a,a)') 'AUXDATA FLUXTYPE="',trim(adjustl(stmp)),'"'

    !GAMMA
    write(stmp,'(100(f14.6))')Gamma_I(1)
    write(iUnitHere,'(a,a,a)') 'AUXDATA GAMMA="',trim(adjustl(stmp)),'"'

    !ITER
    write(stmp,'(i12)')n_step
    write(iUnitHere,'(a,a,a)') 'AUXDATA ITER="',trim(adjustl(stmp)),'"'

    !NPROC
    write(stmp,'(i12)')nProc
    write(iUnitHere,'(a,a,a)') 'AUXDATA NPROC="',trim(adjustl(stmp)),'"'

    !ORDER
    if(nOrder > 1)then
       write(stmp,'(i12,a,f8.5)') &
            nOrder,' '//trim(TypeLimiter)//', beta=',BetaLimiter
    else
       write(stmp,'(i12)') nOrder
    end if
    write(iUnitHere,'(a,a,a)') 'AUXDATA ORDER="',trim(adjustl(stmp)),'"'

    !RBODY
    write(stmp,'(f12.2)')rBody
    write(iUnitHere,'(a,a,a)') 'AUXDATA RBODY="',trim(adjustl(stmp)),'"'

    !SAVEDATE
    call Date_and_time (real_date, real_time)
    write(stmp,'(a11,a4,a1,a2,a1,a2, a4,a2,a1,a2,a1,a2)') &
         'Save Date: ', real_date(1:4),'/',real_date(5:6),'/',real_date(7:8), &
         ' at ',  real_time(1:2),':',real_time(3:4),':',real_time(5:6)
    write(iUnitHere,'(a,a,a)') 'AUXDATA SAVEDATE="',trim(adjustl(stmp)),'"'

    !TIMEEVENT
    write(stmp,'(a)')textDateTime
    write(iUnitHere,'(a,a,a)') 'AUXDATA TIMEEVENT="',trim(adjustl(stmp)),'"'

    !TIMEEVENTSTART
    write(stmp,'(a)')textDateTime0
    write(iUnitHere,'(a,a,a)') 'AUXDATA TIMEEVENTSTART="',trim(adjustl(stmp)),'"'

    !TIMESIM
    if(time_accurate)then
       write(stmp,'(a)')'T='// &
            StringDateOrTime(1:4)//":"// &
            StringDateOrTime(5:6)//":"// &
            StringDateOrTime(7:8)
    else
       write(stmp,'(a)')'T= N/A'
    end if
    write(iUnitHere,'(a,a,a)') 'AUXDATA TIMESIM="',trim(adjustl(stmp)),'"'

    !TIMESIMSHORT
    if(time_accurate)then
       write(stmp,'(a)')'T='// &
            StringDateOrTime(1:4)//":"// &
            StringDateOrTime(5:6)
    else
       write(stmp,'(a)')'T= SS'
    end if
    write(iUnitHere,'(a,a,a)') 'AUXDATA TIMESIMSHORT="',trim(adjustl(stmp)),'"'

  end subroutine write_tecplot_auxdata
  !===========================================================================
  subroutine write_tecplot_setinfo

    use ModProcMH, ONLY: iProc
    use ModMain, ONLY: n_step, time_accurate
    use ModIO, ONLY: StringDateOrTime

    integer :: iTime0_I(7),iTime_I(7)
    character (len=80) :: format
    !-----------------------------------------------------------------------
    call count_true_cells

    if(iProc /= 0) RETURN

    ! Create text string for zone name like 'N=0002000 T=0000:05:00'
    if(time_accurate)then
       call get_time_string
       write(textNandT,'(a,i7.7,a)') "N=",n_step," T="// &
            StringDateOrTime(1:4)//":"// &
            StringDateOrTime(5:6)//":"// &
            StringDateOrTime(7:8)
    else
       write(textNandT,'(a,i7.7)') &
            "N=",n_step
    end if

    format='(i4.4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2,":",i2.2,".",i3.3)'
    call get_date_time_start(iTime0_I)
    call get_date_time(iTime_I)
    write(textDateTime0,format) iTime0_I
    write(textDateTime ,format) iTime_I

  end subroutine write_tecplot_setinfo

end module ModWriteTecplot
