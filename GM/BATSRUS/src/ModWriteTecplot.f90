!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModWriteTecplot

  use BATL_lib, ONLY: &
       test_start, test_stop, lVerbose

  ! Save cell centered data into Tecplot files
  !
  ! Save cell centers (same as 3D IDL plots, but no cell size).
  ! For single data file:
  !     record size is nDim+nPlotVar reals + 1 newline
  !     record index is going from 1 to nPointAll ordered per processor
  ! For single connectivity file:
  !     record size is 2**nDim integers (global cell indexes) + 1 newline
  !     record index is going from 1 to nBrickAll ordered per processor
  ! Connectivity at block boundaries with no resolution change
  !     is written by the side towards the negative direction
  ! Connectivity at resolution changes is written by the finer side.
  ! The connectivity brick contains i:i+iPlotDim,j:j+jPlotDim,k:k+kPlotDim
  !     cell centers.
  !     Indexes belonging to ghost cells are taken from the neighbor.
  !     These are obtained by message passing the cell indexes.

  use BATL_size, ONLY: nDim

  implicit none

  SAVE

  private ! except

  public:: write_tecplot_data      ! write plot variables for 1 block
  public:: write_tecplot_head      ! write header information
  public:: write_tecplot_connect   ! write connectivity file(s)
  public:: write_tecplot_setinfo   ! set some variables to be written as AUX
  public:: write_tecplot_auxdata   ! write auxiliary information
  public:: write_tecplot_node_data ! write node averaged tecplot output
  public:: assign_node_numbers     ! assign node numbers for node averaged plot
  public:: set_tecplot_var_string  ! set variable and unit names

  character(len=23), public :: textDateTime
  character(len=22), public :: textNandT
  character, public, parameter :: CharNewLine = char(10)
  integer,   public :: lRecConnect = 2**nDim*11+1

  ! Dimensionality of plot
  integer, public :: nPlotDim = nDim

  ! Local variables
  character (len=23) :: textDateTime0

  ! Index limits inside the block based on the cut
  integer:: IjkMin_D(3), IjkMax_D(3)

  ! Cuts
  logical:: DoCut

  ! Generalized coordinate limits of the cut box/slice
  real:: CutMin_D(3), CutMax_D(3)

  ! The dimension normal to the cut plane (=0 for no cut or 3D cut box)
  integer:: iCutDim

  ! iPlot_D is 1 if the plot has non-zero width in that dimension
  integer:: iPlot_D(3)

  ! Global cell index. It is a real array for sake of message passing only.
  real, allocatable:: CellIndex_GB(:,:,:,:)

  ! Total number of cells written into file
  integer:: nPointAll

  ! Total number of connectivity bricks
  integer:: nBrickAll

  character(len=80) :: StringFormat

contains
  !============================================================================
  subroutine write_tecplot_data(iBlock, nPlotVar, PlotVar_GV)

    use BATL_lib,  ONLY: MaxDim, nDim, nJ, nK, nIjk_D, &
         MinI, MaxI, MinJ, MaxJ, MinK, MaxK, Xyz_DGB, &
         r_, Phi_, Theta_, Lat_, CoordMin_DB, CoordMax_DB, &
         IsAnyAxis, IsLatitudeAxis, IsSphericalAxis, IsCylindricalAxis
    use ModNumConst, ONLY: cPi, cHalfPi
    use ModIO,     ONLY: nPlotVarMax, DoSaveOneTecFile
    use ModIoUnit, ONLY: UnitTmp_

    integer, intent(in):: iBlock, nPlotVar
    real,    intent(in):: PlotVar_GV(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nPlotvarMax)

    integer:: i, j, k, iMin, iMax, jMin, jMax, kMin, kMax
    integer:: iRecData

    real:: Xyz_D(MaxDim)
    real, allocatable:: PlotVar_V(:)

    ! Interpolation
    integer:: Di, Dj, Dk
    real:: CoefL, CoefR

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_tecplot_data'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(DoTest) write(*,*) NameSub,' starting with nPlotVar=', nPlotVar

    IjkMin_D = 1
    IjkMax_D = nIjk_D

    if(DoCut)then
       ! Check if block is inside cut and set interpolation info if needed
       if(do_skip_block(iBlock, Di, Dj, Dk, CoefL, CoefR)) RETURN
    endif

    iMin = IjkMin_D(1); iMax = IjkMax_D(1)
    jMin = IjkMin_D(2); jMax = IjkMax_D(2)
    kMin = IjkMin_D(3); kMax = IjkMax_D(3)

    allocate(PlotVar_V(nPlotVar))

    if(DoSaveOneTecFile)then
       write(StringFormat, '(a,i2,a)') "(", nDim+nPlotVar, "(ES14.6), a)"
       do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
          ! Record index is the global cell index
          iRecData = nint(CellIndex_GB(i,j,k,iBlock))
          ! Skip points outside the cut
          if(iRecData == 0) CYCLE
          call set_xyz_state
          write(UnitTmp_, StringFormat, REC=iRecData) &
               Xyz_D(1:nDim), PlotVar_V, CharNewLine
       end do; end do; end do
    else
       do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
          ! Skip points outside the cut
          if(CellIndex_GB(i,j,k,iBlock) == 0.0) CYCLE
          call set_xyz_state
          write(UnitTmp_,"(50(ES14.6))") Xyz_D(1:nDim), PlotVar_V
       end do; end do; end do
    end if

    deallocate(PlotVar_V)

    if(DoTest) write(*,*) NameSub,' finished'

    call test_stop(NameSub, DoTest, iBlock)
  contains
    !==========================================================================
    subroutine set_xyz_state

      ! Interpolate variables and coordinates to cut planes
      ! Set the coordinates and fix them to fill the hole at the pole
      !------------------------------------------------------------------------
      if(iCutDim == 0)then
         Xyz_D = Xyz_DGB(:,i,j,k,iBlock)
         PlotVar_V = PlotVar_GV(i,j,k,1:nPlotVar)
      else
         ! Interpolate using precalculated coefficients
         Xyz_D     = CoefL*Xyz_DGB(:,i  ,j,k,iBlock) &
              +      CoefR*Xyz_DGB(:,i+Di,j+Dj,k+Dk,iBlock)
         PlotVar_V = CoefL*PlotVar_GV(i  ,j,k,1:nPlotVar) &
              +      CoefR*PlotVar_GV(i+Di,j+Dj,k+Dk,1:nPlotVar)
      end if

      ! No need to push points to the pole if there is no axis
      ! or the 2D cut is along the Phi direction
      if(.not.IsAnyAxis .or. iCutDim == Phi_) RETURN

      if(IsLatitudeAxis)then
         if(  k==1  .and. CoordMin_DB(Lat_,iBlock) <= -cHalfPi .or. &
              k==nK .and. CoordMax_DB(Lat_,iBlock) >= +cHalfPi) &
              Xyz_D(1:2) = 0.0
      elseif(IsSphericalAxis)then
         if(  j==1  .and. CoordMin_DB(Theta_,iBlock) <= 0.0 .or. &
              j==nJ .and. CoordMax_DB(Theta_,iBlock) >= cPi) &
              Xyz_D(1:2) = 0.0
      elseif(IsCylindricalAxis)then
         if(  i==1 .and. CoordMin_DB(r_,iBlock) <= 0.0) &
              Xyz_D(1:2) = 0.0
      end if
    end subroutine set_xyz_state
    !==========================================================================

  end subroutine write_tecplot_data
  !============================================================================

  subroutine write_tecplot_connect(iFile, NameFile)

    use ModProcMH,    ONLY: iProc, nProc, iComm
    use ModAdvance,   ONLY: iTypeAdvance_BP, SkippedBlock_
    use ModIO,        ONLY: DoSaveOneTecFile, plot_type1, plot_range
    use ModIoUnit,    ONLY: UnitTmp_
    use ModUtilities, ONLY: open_file, close_file
    use ModMain,      ONLY: nBlockMax
    use ModMpi,       ONLY: MPI_SUM, mpi_reduce_integer_scalar, &
         MPI_allgather, MPI_INTEGER
    use ModNumConst,  ONLY: cHalfPi, cPi
    use BATL_lib,     ONLY: nI, nJ, nK, nIJK, k0_, nKp1_, &
         MaxBlock, nBlock, Unused_B, nNodeUsed, &
         CoordMin_DB, CoordMax_DB, Xyz_DGB, &
         IsAnyAxis, IsLatitudeAxis, IsSphericalAxis, IsCylindricalAxis, &
         r_, Phi_, Theta_, Lat_, &
         DiLevelNei_IIIB, message_pass_cell, set_tree_periodic

    ! Write out connectivity file

    integer,          intent(in):: iFile     ! index of plot file
    character(len=*), intent(in):: NameFile  ! name of connectivity file

    integer:: nBlockBefore, iCell, iBlock, iError

    integer:: i0, i1, j0, j1, k0, k1, i, j, k

    ! Cell index to be written into the connectivity file
    integer, allocatable:: iCell_G(:,:,:)

    ! Number of bricks for this and all processors
    integer:: nBrick
    integer, allocatable:: nBrick_P(:)

    ! Multiple stages may be needed
    integer:: iStage, nStage

    ! Cut related variables
    integer:: iPlotDim, jPlotDim, kPlotDim
    logical:: IsPlotDim1, IsPlotDim2, IsPlotDim3

    integer, allocatable:: nCell_P(:)

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_tecplot_connect'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    DoCut = nDim == 3 .and. plot_type1(1:2) /= '3d'

    ! Set iPlotDim_D to 0 in ignored dimensions
    iPlot_D = 0
    iPlot_D(1:nDim) = 1
    nPlotDim = nDim

    ! Set the normal direction of a 2D cut plane to 0 as default (no cut)
    iCutDim = 0

    if(DoTest)write(*,*) NameSub,' starting with NameFile, DoCut=', &
         NameFile, DoCut

    allocate( &
         CellIndex_GB(0:nI+1,0:nJ+1,k0_:nKp1_,MaxBlock), &
         iCell_G(0:nI+1,0:nJ+1,k0_:nKp1_))

    if(DoCut)then
       CutMin_D = plot_range(1:5:2,iFile)
       CutMax_D = plot_range(2:6:2,iFile)

       ! Set iPlot_D to 0 in 0 width cut direction
       where(CutMin_D == CutMax_D) iPlot_D = 0
       nPlotDim = sum(iPlot_D)

       ! Index of the dimension normal to the slice
       if(nPlotDim < nDim)then
          do iCutDim = 1, nDim
             if(iPlot_D(iCutDim) == 0) EXIT
          end do
       end if

       if(DoTest .or. nPlotDim == 1)then
          write(*,*) NameSub,' CutMin_D=',  CutMin_D
          write(*,*) NameSub,' CutMax_D=',  CutMax_D
          write(*,*) NameSub,' iPlot_D =', iPlot_D
          write(*,*) NameSub,' iCutDim =', iCutDim

          if(nPlotDim == 1) &
               call stop_mpi(NameSub//': only 2D cuts are possible here')
       end if

       ! For more than 1 processors:
       ! in stage 1 count number of cells written by this processor
       ! In stage 2 set global cell index
       nStage = min(2, nProc)

       iCell = 0
       do iStage = 1, nStage
          do iBlock = 1, nBlock
             if(iStage == nStage) CellIndex_GB(:,:,:,iBlock) = 0.0

             ! Check if block is inside cut and set IjkMin_D, IjkMax_D...
             if(do_skip_block(iBlock)) CYCLE

             if(iStage < nStage)then
                ! Just count the number of cells inside the cut
                iCell = iCell + product(max(0, IjkMax_D - IjkMin_D + 1))
                CYCLE
             end if

             ! Set global cell index for each cell inside the cut
             do k = IjkMin_D(3), IjkMax_D(3)
                do j = IjkMin_D(2), IjkMax_D(2)
                   do i = IjkMin_D(1), IjkMax_D(1)
                      iCell = iCell + 1
                      CellIndex_GB(i,j,k,iBlock) = iCell
                   end do
                end do
             end do
          end do
          if(iStage < nStage)then
             ! Collect the number of cells written by other processors
             allocate(nCell_P(0:nProc-1))
             call MPI_allgather(iCell, 1, MPI_INTEGER, &
                  nCell_P, 1, MPI_INTEGER, iComm, iError)
             if(iProc == 0)then
                iCell = 0
                ! Store total number of points saved
                nPointAll = sum(nCell_P)
             else
                ! Add up number of cells written by the previous processors
                iCell = sum(nCell_P(0:iProc-1))
             end if
             deallocate(nCell_P)
          end if
       end do
       ! Store total number of points saved when running on 1 processor
       if(nStage == 1) nPointAll = iCell
    else
       ! Full 2D or 3D plot
       ! count number of cells written by processors before this one
       nBlockBefore = 0
       if(iProc > 0) nBlockBefore = &
            count(iTypeAdvance_BP(1:nBlockMax,0:iProc-1) /= SkippedBlock_)

       if(DoTest)write(*,*) NameSub,' nBlockBefore=', nBlockBefore

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
    end if

    ! Set the global cell indexes for the ghost cells. First order prolongation
    ! is used, so fine ghost cells are set to the index of the coarse cell.
    ! Note that the coarse ghost cells are not going to be used.
    ! Just to be safe, we use the minimum index which is better than the
    ! average index.
    call message_pass_cell(1, CellIndex_GB, nProlongOrderIn=1, &
         NameOperatorIn='min')

    ! switch off "fake" periodicity so there are no connections
    call set_tree_periodic(.false.)

    ! Some useful scalars
    iPlotDim = iPlot_D(1); jPlotDim = iPlot_D(2); kPlotDim = iPlot_D(3)
    IsPlotDim1 = iPlotDim>0; IsPlotDim2 = jPlotDim>0; IsPlotDim3 = kPlotDim>0

    if(DoSaveOneTecFile)then
       ! Two stages are needed to figure out the global brick indexes
       nStage = 2
       allocate(nBrick_P(0:nProc-1))

       ! Open connectivity file as direct access
       ! Record length for connectivity file
       lRecConnect = 2**nPlotDim*11+1
       call open_file(File=NameFile, ACCESS='direct', RECL=lRecConnect, &
            iComm=iComm, NameCaller=NameSub//'_direct_connect')
    else
       nStage = 1
       ! Open connectivity file
       call open_file(File=NameFile,  NameCaller=NameSub//'_connect')
    end if

    nBrick = 0
    do iStage = 1, nStage
       do iBlock = 1, nBlock
          ! Check if block is used and inside cut
          if(do_skip_block(iBlock)) CYCLE

          ! Get the index limits for the lower left corner
          ! of the connectivity bricks.

          ! Start connectivity brick from index 1 unless coarser left neighbor
          i0 = 1
          if(IsPlotDim1 .and. DiLevelNei_IIIB(-1,0,0,iBlock) == 1) i0 = 0
          j0 = 1
          if(IsPlotDim2 .and. DiLevelNei_IIIB(0,-1,0,iBlock) == 1) j0 = 0
          k0 = 1
          if(IsPlotDim3 .and. DiLevelNei_IIIB(0,0,-1,iBlock) == 1) k0 = 0

          ! Finish at nI unless there is no block on the right or it is finer
          i1 = nI
          if(IsPlotDim1 .and. DiLevelNei_IIIB(1,0,0,iBlock) < 0) i1 = nI-1
          j1 = nJ
          if(IsPlotDim2 .and. DiLevelNei_IIIB(0,1,0,iBlock) < 0) j1 = nJ-1
          k1 = nK
          if(IsPlotDim3 .and. DiLevelNei_IIIB(0,0,1,iBlock) < 0) k1 = nK-1

          ! For cuts in the Phi direction we want connectivity across poles
          ! Both sides are at either minimum or maximum so we use the
          ! Y coordinate to select which side takes care of the connection
          ! We assume no resolution change across the poles for simplicity
          if(iCutDim == Phi_ .and. IsAnyAxis &
               .and. Xyz_DGB(2,1,1,1,iBlock) > 0.0)then
             if(IsLatitudeAxis)then
                if(CoordMin_DB(Lat_,iBlock) <= -cHalfPi) k0 = 0
                if(CoordMax_DB(Lat_,iBlock) >= +cHalfPi) k1 = nK - 1
             elseif(IsSphericalAxis)then
                if(CoordMin_DB(Theta_,iBlock) <= 0.0) j0 = 0
                if(CoordMax_DB(Theta_,iBlock) >= cPi) j1 = nJ - 1
             elseif(IsCylindricalAxis)then
                if(CoordMin_DB(r_,iBlock) <= 0.0) i0 = 0
             end if
          end if

          iCell_G = nint(CellIndex_GB(:,:,:,iBlock))

          ! Coarse ghost cells should not be used.
          ! The face ghost cells are surely not used due to above settings
          ! The "left" edge and corner ghost cells are only used if the
          ! left face neighbor is coarser, so these cannot be prolonged cells.
          ! So only "right" edge and corner ghost cells have to be zeroed out.

          if(IsPlotDim1 .and. IsPlotDim2 .and. &
               DiLevelNei_IIIB(1,1,0,iBlock) < 0) iCell_G(nI+1,nJ+1,:) = 0
          if(IsPlotDim1 .and. IsPlotDim3 .and. &
               DiLevelNei_IIIB(1,0,1,iBlock) < 0) iCell_G(nI+1,:,nK+1) = 0
          if(IsPlotDim2 .and. IsPlotDim3 .and. &
               DiLevelNei_IIIB(0,1,1,iBlock) < 0) iCell_G(:,nJ+1,nK+1) = 0
          if(nPlotDim == 3 .and. &
               DiLevelNei_IIIB(1,1,1,iBlock) < 0) iCell_G(nI+1,nJ+1,nK+1) = 0

          ! Loop over the "lower-left" corner of the bricks
          do k = k0, k1; do j = j0, j1; do i = i0, i1
             ! Skip bricks that are not fully inside/usable
             if(any(iCell_G(i:i+iPlotDim,j:j+jPlotDim,k:k+kPlotDim)==0))&
                  CYCLE
             nBrick = nBrick + 1

             ! In stage 1 only count bricks
             if(iStage < nStage) CYCLE

             if(nPlotDim == 3)then
                if(DoSaveOneTecFile)then
                   write(UnitTmp_,'(8i11,a)', REC=nBrick) &
                        iCell_G(i  ,j  ,k  ), &
                        iCell_G(i+1,j  ,k  ), &
                        iCell_G(i+1,j+1,k  ), &
                        iCell_G(i  ,j+1,k  ), &
                        iCell_G(i  ,j  ,k+1), &
                        iCell_G(i+1,j  ,k+1), &
                        iCell_G(i+1,j+1,k+1), &
                        iCell_G(i  ,j+1,k+1), &
                        CharNewLine
                else
                   write(UnitTmp_,'(8i11)') &
                        iCell_G(i  ,j  ,k  ), &
                        iCell_G(i+1,j  ,k  ), &
                        iCell_G(i+1,j+1,k  ), &
                        iCell_G(i  ,j+1,k  ), &
                        iCell_G(i  ,j  ,k+1), &
                        iCell_G(i+1,j  ,k+1), &
                        iCell_G(i+1,j+1,k+1), &
                        iCell_G(i  ,j+1,k+1)
                end if
             elseif(.not.IsPlotDim3)then
                if(DoSaveOneTecFile)then
                   write(UnitTmp_,'(4i11,a)', REC=nBrick) &
                        iCell_G(i  ,j  ,k), &
                        iCell_G(i+1,j  ,k), &
                        iCell_G(i+1,j+1,k), &
                        iCell_G(i  ,j+1,k), &
                        CharNewLine
                else
                   write(UnitTmp_,'(4i11)') &
                        iCell_G(i  ,j  ,k), &
                        iCell_G(i+1,j  ,k), &
                        iCell_G(i+1,j+1,k), &
                        iCell_G(i  ,j+1,k)
                end if
             elseif(.not.IsPlotDim2)then
                if(DoSaveOneTecFile)then
                   write(UnitTmp_,'(4i11,a)', REC=nBrick) &
                        iCell_G(i  ,j,k  ), &
                        iCell_G(i+1,j,k  ), &
                        iCell_G(i+1,j,k+1), &
                        iCell_G(i  ,j,k+1), &
                        CharNewLine
                else
                   write(UnitTmp_,'(4i11)') &
                        iCell_G(i  ,j,k  ), &
                        iCell_G(i+1,j,k  ), &
                        iCell_G(i+1,j,k+1), &
                        iCell_G(i  ,j,k+1)
                end if
             elseif(.not.IsPlotDim1)then
                if(DoSaveOneTecFile)then
                   write(UnitTmp_,'(4i11,a)', REC=nBrick) &
                        iCell_G(i,j  ,k  ), &
                        iCell_G(i,j+1,k  ), &
                        iCell_G(i,j+1,k+1), &
                        iCell_G(i,j  ,k+1), &
                        CharNewLine
                else
                   write(UnitTmp_,'(4i11)') &
                        iCell_G(i,j  ,k  ), &
                        iCell_G(i,j+1,k  ), &
                        iCell_G(i,j+1,k+1), &
                        iCell_G(i,j  ,k+1)
                end if
             end if
          end do; end do; end do
       end do ! iBlock
       if(iStage < nStage)then
          ! Collect number of bricks from all processors
          call MPI_allgather(nBrick, 1, MPI_INTEGER, &
               nBrick_P, 1, MPI_INTEGER, iComm, iError)
          ! Add up number of bricks on previous prorecssors
          nBrick = 0
          if(iProc > 0) nBrick = sum(nBrick_P(0:iProc-1))
       end if
    end do ! iStage
    call close_file

    ! Reset periodicity as it was
    call set_tree_periodic(.true.)

    ! Calculate and store total number of bricks for header file
    if(DoSaveOneTecFile)then
       if(iProc == 0) nBrickAll = sum(nBrick_P)
       deallocate(nBrick_P)
    else
       nBrickAll = nBrick
       if(nProc > 1) call mpi_reduce_integer_scalar(nBrickAll, &
            MPI_SUM, 0, iComm, iError)
    end if

    ! Calculate and store total number of cells for header file
    if(.not. DoCut .and. iProc == 0) nPointAll = nNodeUsed*nIJK

    deallocate(iCell_G)

    if(DoTest)write(*,*) NameSub,' done with nPointAll, nBrickAll=', &
         nPointAll, nBrickAll

    call test_stop(NameSub, DoTest)
  end subroutine write_tecplot_connect
  !============================================================================

  subroutine write_tecplot_head(NameFile, StringUnit)

    use ModProcMH,    ONLY: iProc
    use ModIoUnit,    ONLY: UnitTmp_
    use ModUtilities, ONLY: open_file, close_file

    ! Write out tecplot header file

    character(len=*), intent(in):: NameFile
    character(len=*), intent(in):: StringUnit

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_tecplot_head'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(allocated(CellIndex_GB)) deallocate(CellIndex_GB)

    call write_tecplot_setinfo
    if(iProc /= 0) RETURN

    call open_file(File=NameFile)
    if(DoCut)then
       write(UnitTmp_,'(a)')'TITLE="BATSRUS: cut Data, '//textDateTime//'"'
    else
       write(UnitTmp_,'(a,i1,a)') &
            'TITLE="BATSRUS: ', nDim,'D Data,'//textDateTime//'"'
    end if
    write(UnitTmp_,'(a)') trim(StringUnit)
    select case(nPlotDim)
    case(2)
       write(UnitTmp_,'(a,a,i12,a,i12,a)') &
            'ZONE T="2D   '//textNandT//'"', &
            ', N=', nPointAll, &
            ', E=', nBrickAll, &
            ', F=FEPOINT, ET=QUADRILATERAL'
    case(3)
       write(UnitTmp_,'(a,a,i12,a,i12,a)') &
            'ZONE T="3D   '//textNandT//'"', &
            ', N=', nPointAll, &
            ', E=', nBrickAll, &
            ', F=FEPOINT, ET=BRICK'
    end select
    call write_tecplot_auxdata
    call close_file

    call test_stop(NameSub, DoTest)
  end subroutine write_tecplot_head
  !============================================================================

  subroutine write_tecplot_auxdata(iUnitIn)

    use ModMain, ONLY : nI,nJ,nK, &
         nBlockALL, time_accurate,n_step, &
         nOrder, UseRotatingBc,           &
         TypeCoordSystem, CodeVersion
    use ModGeometry, ONLY: nTrueCells
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

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_tecplot_auxdata'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if (present(iUnitIn)) then
       iUnitHere = iUnitIn
    else
       iUnitHere = UnitTmp_
    end if

    ! BLOCKS
    write(stmp,'(i12,3(a,i2))')nBlockALL,'  ',nI,' x',nJ,' x',nK
    write(iUnitHere,'(a,a,a)') 'AUXDATA BLOCKS="',trim(adjustl(stmp)),'"'

    ! BODYDENSITY
    write(stmp,'(f12.2)')BodyNDim_I(IonFirst_)
    write(iUnitHere,'(a,a,a)') &
         'AUXDATA BODYNUMDENSITY="',trim(adjustl(stmp)),'"'

    ! BORIS
    if(boris_correction)then
       write(stmp,'(a,f8.4)')'T ',boris_cLIGHT_factor
    else
       write(stmp,'(a)')'F'
    end if
    write(iUnitHere,'(a,a,a)') 'AUXDATA BORIS="',trim(adjustl(stmp)),'"'

    ! BTHETATILT
    write(stmp,'(f12.4)')ThetaTilt*cRadToDeg
    write(iUnitHere,'(a,a,a)') 'AUXDATA BTHETATILT="',trim(adjustl(stmp)),'"'

    ! CELLS
    write(stmp,'(i12)')nBlockALL*nIJK
    write(iUnitHere,'(a,a,a)') 'AUXDATA CELLS="',trim(adjustl(stmp)),'"'

    ! CELLSUSED
    write(stmp,'(i12)')nTrueCells
    write(iUnitHere,'(a,a,a)') 'AUXDATA CELLSUSED="',trim(adjustl(stmp)),'"'

    ! CODEVERSION
    write(stmp,'(a,f5.2)')'BATSRUS',CodeVersion
    write(iUnitHere,'(a,a,a)') 'AUXDATA CODEVERSION="',trim(adjustl(stmp)),'"'

    ! COORDSYSTEM
    write(stmp,'(a)')TypeCoordSystem
    write(iUnitHere,'(a,a,a)') 'AUXDATA COORDSYSTEM="',trim(adjustl(stmp)),'"'

    ! COROTATION
    if(UseRotatingBc)then
       write(stmp,'(a)')'T'
    else
       write(stmp,'(a)')'F'
    end if
    write(iUnitHere,'(a,a,a)') 'AUXDATA COROTATION="',trim(adjustl(stmp)),'"'

    ! FLUXTYPE
    write(stmp,'(a)')FluxType
    write(iUnitHere,'(a,a,a)') 'AUXDATA FLUXTYPE="',trim(adjustl(stmp)),'"'

    ! GAMMA
    write(stmp,'(100(f14.6))')Gamma_I(1)
    write(iUnitHere,'(a,a,a)') 'AUXDATA GAMMA="',trim(adjustl(stmp)),'"'

    ! ITER
    write(stmp,'(i12)')n_step
    write(iUnitHere,'(a,a,a)') 'AUXDATA ITER="',trim(adjustl(stmp)),'"'

    ! NPROC
    write(stmp,'(i12)')nProc
    write(iUnitHere,'(a,a,a)') 'AUXDATA NPROC="',trim(adjustl(stmp)),'"'

    ! ORDER
    if(nOrder > 1)then
       write(stmp,'(i12,a,f8.5)') &
            nOrder,' '//trim(TypeLimiter)//', beta=',BetaLimiter
    else
       write(stmp,'(i12)') nOrder
    end if
    write(iUnitHere,'(a,a,a)') 'AUXDATA ORDER="',trim(adjustl(stmp)),'"'

    ! RBODY
    write(stmp,'(f12.2)')rBody
    write(iUnitHere,'(a,a,a)') 'AUXDATA RBODY="',trim(adjustl(stmp)),'"'

    ! SAVEDATE
    call Date_and_time (real_date, real_time)
    write(stmp,'(a11,a4,a1,a2,a1,a2, a4,a2,a1,a2,a1,a2)') &
         'Save Date: ', real_date(1:4),'/',real_date(5:6),'/',real_date(7:8), &
         ' at ',  real_time(1:2),':',real_time(3:4),':',real_time(5:6)
    write(iUnitHere,'(a,a,a)') 'AUXDATA SAVEDATE="',trim(adjustl(stmp)),'"'

    ! TIMEEVENT
    write(stmp,'(a)')textDateTime
    write(iUnitHere,'(a,a,a)') 'AUXDATA TIMEEVENT="',trim(adjustl(stmp)),'"'

    ! TIMEEVENTSTART
    write(stmp,'(a)')textDateTime0
    write(iUnitHere,'(a,a,a)') 'AUXDATA TIMEEVENTSTART="',trim(adjustl(stmp)),'"'

    ! TIMESIM
    if(time_accurate)then
       write(stmp,'(a)')'T='// &
            StringDateOrTime(1:4)//":"// &
            StringDateOrTime(5:6)//":"// &
            StringDateOrTime(7:8)
    else
       write(stmp,'(a)')'T= N/A'
    end if
    write(iUnitHere,'(a,a,a)') 'AUXDATA TIMESIM="',trim(adjustl(stmp)),'"'

    ! TIMESIMSHORT
    if(time_accurate)then
       write(stmp,'(a)')'T='// &
            StringDateOrTime(1:4)//":"// &
            StringDateOrTime(5:6)
    else
       write(stmp,'(a)')'T= SS'
    end if
    write(iUnitHere,'(a,a,a)') 'AUXDATA TIMESIMSHORT="',trim(adjustl(stmp)),'"'

    call test_stop(NameSub, DoTest)
  end subroutine write_tecplot_auxdata
  !============================================================================
  subroutine write_tecplot_setinfo

    use ModProcMH, ONLY: iProc
    use ModMain, ONLY: n_step, time_accurate, iStartTime_I
    use ModIO, ONLY: StringDateOrTime
    use ModGeometry, ONLY: count_true_cells

    integer :: iTime_I(7)
    character (len=80) :: format
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_tecplot_setinfo'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
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
    call get_date_time(iTime_I)
    write(textDateTime0,format) iStartTime_I
    write(textDateTime ,format) iTime_I

    call test_stop(NameSub, DoTest)
  end subroutine write_tecplot_setinfo
  !============================================================================
  logical function do_skip_block(iBlock, Di, Dj, Dk, CoefL, CoefR)

    ! Set the logical value to false if the block is unused
    ! or if it is outside the cut.
    ! If the block intersects the cut, then set IjkMin_D and IjkMax_D
    ! index limits.
    ! If the index shifts Di, Dj, Dk and the interpolation coefficients
    ! CoefL and CoefR are present
    ! and iCutDim is not zero, set the index shifts and
    ! interpolation coefficients to interpolate onto the cut plane.

    use BATL_lib, ONLY: Unused_B, CoordMin_DB, CoordMax_DB, CellSize_DB, &
         nIjk_D, Phi_
    use ModNumConst, ONLY: i_DD, cPi, cHalfPi

    integer, intent(in):: iBlock

    integer, intent(out), optional:: Di, Dj, Dk
    real,    intent(out), optional:: CoefL, CoefR

    ! Arrays describing block geometry
    real:: BlockMin_D(3), BlockMax_D(3), CellSize_D(3)

    real:: PhiCut, PhiBlock, CutCoordNorm

    !--------------------------------------------------------------------------
    do_skip_block = UnUsed_B(iBlock)
    if(do_skip_block) RETURN

    ! There is nothing else to check if there is no cut
    if(.not.DoCut) RETURN

    ! Block coordinates
    BlockMin_D = CoordMin_DB(:,iBlock)
    BlockMax_D = CoordMax_DB(:,iBlock)

    ! Shift Phi coordinate for 2D cuts along the Phi coordinate
    ! e.g. x=0 and y=0 cuts of a spherical grid.
    if(iCutDim == Phi_)then
       ! Blocks "far" from the phi cut are shifted 180 degrees towards cut
       ! This will capture both sides of the sphere/cylinder and is the
       ! expected behavior for x=0 and y=0 cuts. But also works for other.
       PhiCut = CutMin_D(Phi_)
       PhiBlock = 0.5*(BlockMin_D(Phi_)+BlockMax_D(Phi_))

       if(PhiBlock - PhiCut > cHalfPi)then
          BlockMin_D(Phi_) = BlockMin_D(Phi_) - cPi
          BlockMax_D(Phi_) = BlockMax_D(Phi_) - cPi
       elseif(PhiCut - PhiBlock > cHalfPi)then
          BlockMin_D(Phi_) = BlockMin_D(Phi_) + cPi
          BlockMax_D(Phi_) = BlockMax_D(Phi_) + cPi
       end if

    end if

    ! Check if block is inside the cut
    do_skip_block = any(BlockMin_D > CutMax_D) .or. any(BlockMax_D < CutMin_D)
    if(do_skip_block) RETURN

    ! Calculate start and ending indexes of cut
    ! The +-0.01 is used to avoid rounding errors and make sure
    ! that there are at least 1 cells in zero-width directions,
    ! and 2 cells in non-zero width directions.
    CellSize_D = CellSize_DB(:,iBlock)
    IjkMin_D = max(1, &
         nint(-0.01 +          (CutMin_D - BlockMin_D)/CellSize_D))
    IjkMax_D = min(nIJK_D, &
         nint(0.01 + iPlot_D + (CutMax_D - BlockMin_D)/CellSize_D))

    ! Nothing else to do if there is no 2D cut or Di is not present
    if(iCutDim == 0 .or. .not. present(Di)) RETURN

    ! Cell index shift in the normal direction for interpolation
    Di = i_DD(1,iCutDim); Dj = i_DD(2,iCutDim); Dk = i_DD(3,iCutDim)

    ! Normalized (to block index) coordinate of the cut
    CutCoordNorm = 0.5 + &
         (CutMin_D(iCutDim) - BlockMin_D(iCutDim))/CellSize_D(iCutDim)
    CoefR = CutCoordNorm - IjkMin_D(iCutDim)
    CoefL = 1.0 - CoefR

  end function do_skip_block
  !============================================================================

  subroutine write_tecplot_node_data( &
       iFile, nPlotVar, PlotVarBlk, PlotVarNodes_VNB, PlotXYZNodes_DNB, &
       unitstr_TEC, xmin, xmax, ymin, ymax, zmin, zmax, iUnit)

    ! NOTE:
    ! This routine assumes that the blocks are sorted on PEs by their global
    ! block number, ie blocks 1 to n on PE 0, blocks n+1 to n+m on PE 1,
    ! etc.

    use ModProcMH
    use ModMain, ONLY : nI,nJ,nK, nBlock, nBlockALL
    use ModPhysics, ONLY : No2Io_V, UnitX_, &
         ThetaTilt
    use ModAdvance, ONLY : iTypeAdvance_B, SkippedBlock_
    use ModIO
    use ModIoUnit, ONLY: UnitTmp_, UnitTmp2_
    use ModNodes, ONLY: nNodeAll, NodeNumberGlobal_NB, NodeUniqueGlobal_NB
    use BATL_lib, ONLY: IsCartesianGrid, IsRLonLat,                    &
         nNodeUsed, iNodeMorton_I, iTree_IA, Block_, Proc_,            &
         Xyz_DGB, MinI, MaxI, MinJ, MaxJ, MinK, MaxK,                  &
         find_grid_block, iMortonNode_A, iNode_B
    use ModMpi

    ! Arguments
    integer, intent(in) :: ifile, nPlotVar
    character (LEN=1000), intent(in) :: unitstr_TEC
    real, intent(in) :: PlotVarBLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nPlotVarMax)
    real, intent(in) :: PlotVarNodes_VNB(nPlotVarMax,nI+1,nJ+1,nK+1,nBlock)
    real, intent(in) :: PlotXYZNodes_DNB(3,nI+1,nJ+1,nK+1,nBlock)
    real, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax
    integer, intent(in) :: iUnit

    ! Local Variables
    integer :: i,j,k, cut1,cut2, iPE,iBlock, iBlockAll, iNode, nBlockCuts, iError
    real :: CutValue, factor1,factor2
    integer, allocatable, dimension(:) :: BlockCut
    character(len=500) :: stmp

    ! parameters for saving 3d tecplot in a single file
    character (len=80) :: formatData
    integer :: iRec

    integer::ic1,ic2,jc1,jc2,kc1,kc2, nCuts, nCutsTotal
    real :: XarbP,YarbP,ZarbP, XarbNormal,YarbNormal,ZarbNormal, Xp,Yp,Zp
    real, dimension(3,1:nI+1,1:nJ+1,1:nK+1) :: NodeXYZ_DN
    logical :: okdebug

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_tecplot_node_data'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(DoTest)write(*,*) plot_type1,plot_type1(1:3)

    call write_tecplot_setinfo

    select case(plot_type1(1:3))
    case('blk')
       call find_grid_block(plot_point(:,iFile), iPE, iBlock)
       if(iPE /= iProc) RETURN

       write(UnitTmp_,'(a)')'TITLE="BATSRUS: BLK Only, '//textDateTime//'"'
       write(UnitTmp_,'(a)')trim(unitstr_TEC)
       write(UnitTmp_,'(a,i8,a,i8,a,i8,a)') &
            'ZONE T="BLK Only '//textNandT//'", I=',nI+4,&
            ', J=',nJ+4,', K=',nK+4,', F=POINT'
       call write_tecplot_auxdata
       ! DEBUGBLK
       write(stmp,'(i12)')iBlock
       write(UnitTmp_,'(a,a,a)') 'AUXDATA DEBUGBLK="',trim(adjustl(stmp)),'"'
       ! DEBUGPROC
       write(stmp,'(i12)')iProc
       write(UnitTmp_,'(a,a,a)') 'AUXDATA DEBUGPROC="',trim(adjustl(stmp)),'"'
       ! Write cell values
       do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
          if (plot_dimensional(ifile)) then
             write(UnitTmp_,fmt="(50(ES14.6))") &
                  Xyz_DGB(:,i,j,k,iBlock)*No2Io_V(UnitX_), &
                  PlotVarBlk(i,j,k,1:nPlotVar)
          else
             write(UnitTmp_,fmt="(50(ES14.6))") &
                  Xyz_DGB(:,i,j,k,iBlock), &
                  PlotVarBlk(i,j,k,1:nPlotVar)
          end if
       end do; end do; end do
    case('1d_')
       if(iProc==0)then
          ! Write file header
          write(UnitTmp_,'(a)') &
               'TITLE="BATSRUS: 1D Block Data, '//textDateTime//'"'
          write(UnitTmp_,'(a)')trim(unitstr_TEC)
          write(UnitTmp_,'(a,a,i12,a,a)') &
               'ZONE T="1D   '//textNandT//'"', &
               ', I=',nBlockALL,', J=1, K=1,', &
               ', ZONETYPE=ORDERED, DATAPACKING=POINT'
          call write_tecplot_auxdata
       end if
       !================================= 1d ============================
       do iBlockAll = 1, nNodeUsed
          iNode = iNodeMorton_I(iBlockAll)
          iBlock  = iTree_IA(Block_,iNode)
          iPE   = iTree_IA(Proc_,iNode)
          if(iProc==iPE)then
             ! Write point values
             call fill_NodeXYZ
             write(UnitTmp_,fmt="(50(ES14.6))") &
                  (NodeXYZ_DN(1,1,1,1)+NodeXYZ_DN(1,nI+1,1,1))/2., &
                  (NodeXYZ_DN(2,1,1,1)+NodeXYZ_DN(2,1,nJ+1,1))/2., &
                  (NodeXYZ_DN(3,1,1,1)+NodeXYZ_DN(3,1,1,nK+1))/2., &
                  PlotVarNodes_VNB(1:nPlotVar,2,2,2,iBlock)
          end if
       end do
    case('3d_')
       if(iProc==0)then
          ! For DoSaveOneTecFile = T, iUnit is assoicated with the header
          ! file. iUnit = UnitTMp_ for DoSaveOneTecFile = F
          ! Write file header
          write(iUnit,'(a)')'TITLE="BATSRUS: 3D Data, '//textDateTime//'"'
          write(iUnit,'(a)')trim(unitstr_TEC)
          write(iUnit,'(a,a,i12,a,i12,a)') &
               'ZONE T="3D   '//textNandT//'"', &
               ', N=',nNodeALL, &
               ', E=',nBlockALL*nIJK, &
               ', F=FEPOINT, ET=BRICK'
          call write_tecplot_auxdata(iUnit)
       end if
       !================================= 3d ============================
       if (DoSaveOneTecFile) then
          ! output format for tecplot data: coordinate info + nPlotVar
          write(formatData, '(a,i2.2,a)') "(", nPlotVar+3, "(ES14.6), a)"
          if(DoTest)write(*,*) 'formatData =', formatData
          do iBlock = 1, nBlock
             if(iTypeAdvance_B(iBlock) == SkippedBlock_) CYCLE
             call fill_NodeXYZ

             ! Write point values
             do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
                if(NodeUniqueGlobal_NB(i,j,k,iBlock))then
                   iRec = NodeNumberGlobal_NB(i,j,k,iBlock)
                   write(UnitTmp_, FMT=formatData, REC=iRec) &
                        NodeXYZ_DN(1:3,i,j,k),                  &
                        PlotVarNodes_VNB(1:nPlotVar,i,j,k,iBlock),&
                        CharNewLine
                end if
             end do; end do; end do

             ! Write point connectivity.
             ! Initalize record index based on Morton ordering
             iRec = nIJK*(iMortonNode_A(iNode_B(iBlock)) - 1)
             do k=1,nK; do j=1,nJ; do i=1,nI
                iRec = iRec + 1
                write(UnitTmp2_,'(8i11,a)', REC=iRec) &
                     NodeNumberGlobal_NB(i  ,j  ,k  ,iBlock), &
                     NodeNumberGlobal_NB(i+1,j  ,k  ,iBlock), &
                     NodeNumberGlobal_NB(i+1,j+1,k  ,iBlock), &
                     NodeNumberGlobal_NB(i  ,j+1,k  ,iBlock), &
                     NodeNumberGlobal_NB(i  ,j  ,k+1,iBlock), &
                     NodeNumberGlobal_NB(i+1,j  ,k+1,iBlock), &
                     NodeNumberGlobal_NB(i+1,j+1,k+1,iBlock), &
                     NodeNumberGlobal_NB(i  ,j+1,k+1,iBlock), CharNewLine
             end do; end do; end do
          end do
       else
          do iBlock = 1, nBlock
             if(iTypeAdvance_B(iBlock) == SkippedBlock_) CYCLE
             ! Write point values
             call fill_NodeXYZ
             do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
                if(NodeUniqueGlobal_NB(i,j,k,iBlock))then
                   write(UnitTmp_,fmt="(50(ES14.6))") &
                        NodeXYZ_DN(1:3,i,j,k),       &
                        PlotVarNodes_VNB(1:nPlotVar,i,j,k,iBlock)
                end if
             end do; end do; end do
             ! Write point connectivity
             do k=1,nK; do j=1,nJ; do i=1,nI
                write(UnitTmp2_,'(8i11)') &
                     NodeNumberGlobal_NB(i  ,j  ,k  ,iBlock), &
                     NodeNumberGlobal_NB(i+1,j  ,k  ,iBlock), &
                     NodeNumberGlobal_NB(i+1,j+1,k  ,iBlock), &
                     NodeNumberGlobal_NB(i  ,j+1,k  ,iBlock), &
                     NodeNumberGlobal_NB(i  ,j  ,k+1,iBlock), &
                     NodeNumberGlobal_NB(i+1,j  ,k+1,iBlock), &
                     NodeNumberGlobal_NB(i+1,j+1,k+1,iBlock), &
                     NodeNumberGlobal_NB(i  ,j+1,k+1,iBlock)
             end do; end do; end do
          end do
       end if
    case('cut','x=0','y=0','z=0')
       !================================ cut ============================
       ! Allocate memory for storing the blocks that are cut
       allocate(BlockCut(nBlockALL))
       BlockCut=0
       nBlockCuts=0

       if((xmax-xmin)<(ymax-ymin) .and. (xmax-xmin)<(zmax-zmin))then
          ! X Slice
          CutValue = 0.5*(xmin+xmax)
          if(plot_type1(1:3) == 'x=0') CutValue = 0.
          if(IsCartesianGrid)then
             ! First loop to count nodes and cells
             do iBlockAll = 1, nNodeUsed
                iNode = iNodeMorton_I(iBlockAll)
                iBlock  = iTree_IA(Block_,iNode)
                iPE   = iTree_IA(Proc_,iNode)
                if(iProc==iPE)then
                   if ( CutValue> PlotXYZNodes_DNB(1,1   ,1,1,iBlock) .and. &
                        CutValue<=PlotXYZNodes_DNB(1,1+nI,1,1,iBlock)  )then
                      nBlockCuts=nBlockCuts+1
                      BlockCut(iBlockALL)=nBlockCuts
                   end if
                end if
                call MPI_Bcast(nBlockCuts,1,MPI_Integer,iPE,iComm,iError)
             end do
             if(iProc==0)then
                ! Write file header
                write(UnitTmp_,'(a)')'TITLE="BATSRUS: Cut X Data, ' &
                     //textDateTime//'"'
                write(UnitTmp_,'(a)')trim(unitstr_TEC)
                write(UnitTmp_,'(a,a,i12,a,i12,a)') &
                     'ZONE T="2D X '//textNandT//'"', &
                     ', N=',nBlockCuts*((nJ+1)*(nK+1)), &
                     ', E=',nBlockCuts*((nJ  )*(nK  )), &
                     ', F=FEPOINT, ET=QUADRILATERAL'
                call write_tecplot_auxdata
             end if
             ! Now loop to write values
             do iBlockAll = 1, nNodeUsed
                iNode = iNodeMorton_I(iBlockAll)
                iBlock  = iTree_IA(Block_,iNode)
                iPE   = iTree_IA(Proc_,iNode)
                if(iProc==iPE)then
                   if ( CutValue> PlotXYZNodes_DNB(1,1   ,1,1,iBlock) .and. &
                        CutValue<=PlotXYZNodes_DNB(1,1+nI,1,1,iBlock) )then
                      ! Find cut interpolation factors
                      do i=1,nI
                         if ( CutValue> PlotXYZNodes_DNB(1,i  ,1,1,iBlock) .and.&
                              CutValue<=PlotXYZNodes_DNB(1,i+1,1,1,iBlock)  )then
                            cut1=i
                            cut2=i+1
                            factor2 = &
                                 (CutValue - PlotXYZNodes_DNB(1,i,1,1,iBlock))/ &
                                 ( PlotXYZNodes_DNB(1,i+1,1,1,iBlock) &
                                 - PlotXYZNodes_DNB(1,i,1,1,iBlock))
                            factor1=1.-factor2
                            EXIT
                         end if
                      end do
                      ! Write point values
                      call fill_NodeXYZ
                      do k=1,1+nK; do j=1,1+nJ
                         write(UnitTmp_,fmt="(50(ES14.6))") &
                              (factor1*NodeXYZ_DN(1:3,cut1,j,k)+ &
                              factor2*NodeXYZ_DN(1:3,cut2,j,k)), &
                              (factor1*PlotVarNodes_VNB(1:nPlotVar,cut1,j,k,iBlock)+ &
                              factor2*PlotVarNodes_VNB(1:nPlotVar,cut2,j,k,iBlock))
                      end do; end do
                      ! Write point connectivity
                      do k=1,nK; do j=1,nJ
                         write(UnitTmp2_,'(4(i8,1x))') &
                              ((BlockCut(iBlockALL)-1)*(nJ+1)*(nK+1)) &
                              + (k-1)*(nJ+1)+j, &
                              ((BlockCut(iBlockALL)-1)*(nJ+1)*(nK+1)) &
                              + (k-1)*(nJ+1)+j+1, &
                              ((BlockCut(iBlockALL)-1)*(nJ+1)*(nK+1)) &
                              + (k  )*(nJ+1)+j+1, &
                              ((BlockCut(iBlockALL)-1)*(nJ+1)*(nK+1)) &
                              + (k  )*(nJ+1)+j
                      end do; end do
                   end if
                end if
             end do
          elseif(IsRLonLat)then
             ! First loop to count nodes and cells
             do iBlockAll = 1, nNodeUsed
                iNode = iNodeMorton_I(iBlockAll)
                iBlock  = iTree_IA(Block_,iNode)
                iPE   = iTree_IA(Proc_,iNode)
                if(iProc==iPE)then
                   if(  (CutValue - PlotXYZNodes_DNB(1,1,1   ,2,iBlock))* &
                        (CutValue - PlotXYZNodes_DNB(1,1,1+nJ,2,iBlock)) <=0)then
                      nBlockCuts=nBlockCuts+1
                      BlockCut(iBlockALL)=nBlockCuts
                   end if
                end if
                call MPI_Bcast(nBlockCuts,1,MPI_Integer,iPE,iComm,iError)
             end do
             if(iProc==0)then
                ! Write file header
                write(UnitTmp_,'(a)')'TITLE="BATSRUS: Cut X Data, ' &
                     //textDateTime//'"'
                write(UnitTmp_,'(a)')trim(unitstr_TEC)
                write(UnitTmp_,'(a,a,i12,a,i12,a)')       &
                     'ZONE T="2D X '//textNandT//'"',   &
                     ', N=',nBlockCuts*((nI+1)*(nK+1)), &
                     ', E=',nBlockCuts*((nI  )*(nK  )), &
                     ', F=FEPOINT, ET=QUADRILATERAL'
                call write_tecplot_auxdata
             end if
             ! Now loop to write values
             do iBlockAll = 1, nNodeUsed
                iNode = iNodeMorton_I(iBlockAll)
                iBlock  = iTree_IA(Block_,iNode)
                iPE   = iTree_IA(Proc_,iNode)
                if(iProc==iPE)then
                   if ( (CutValue - PlotXYZNodes_DNB(1,1,1   ,2,iBlock))* &
                        (CutValue - PlotXYZNodes_DNB(1,1,1+nJ,2,iBlock)) <=0)then
                      ! Find cut interpolation factors
                      do j=1, nJ
                         if ( (CutValue - PlotXYZNodes_DNB(1,1,j  ,2,iBlock))* &
                              (CutValue - PlotXYZNodes_DNB(1,1,j+1,2,iBlock)) &
                              <= 0) then
                            cut1=j
                            cut2=j+1
                            factor2 = &
                                 (CutValue - PlotXYZNodes_DNB(1,1,j,2,iBlock))/ &
                                 ( PlotXYZNodes_DNB(1,1,j+1,2,iBlock) &
                                 - PlotXYZNodes_DNB(1,1,j,2,iBlock))
                            factor1=1.-factor2
                            EXIT
                         end if
                      end do
                      ! Write point values
                      call fill_NodeXYZ
                      do k=1,1+nK; do i=1,1+nI
                         write(UnitTmp_,fmt="(50(ES14.6))") &
                              (factor1*NodeXYZ_DN(1:3,i,cut1,k)+ &
                              factor2*NodeXYZ_DN(1:3,i,cut2,k)), &
                              (factor1*PlotVarNodes_VNB(1:nPlotVar,i,cut1,k,iBlock)+ &
                              factor2*PlotVarNodes_VNB(1:nPlotVar,i,cut2,k,iBlock))
                      end do; end do
                      ! Write point connectivity
                      do k=1,nK; do i=1,nI
                         write(UnitTmp2_,'(4(i8,1x))') &
                              ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) &
                              + (k-1)*(nI+1)+i, &
                              ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) &
                              + (k-1)*(nI+1)+i+1, &
                              ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) &
                              + (k  )*(nI+1)+i+1, &
                              ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) &
                              + (k  )*(nI+1)+i
                      end do; end do
                   end if
                end if
             end do
          end if

       elseif((ymax-ymin)<(zmax-zmin))then
          ! Y Slice
          CutValue = 0.5*(ymin+ymax)
          if(plot_type1(1:3) == 'y=0') CutValue = 0.
          if(IsCartesianGrid)then
             ! First loop to count nodes and cells
             do iBlockAll = 1, nNodeUsed
                iNode = iNodeMorton_I(iBlockAll)
                iBlock  = iTree_IA(Block_,iNode)
                iPE   = iTree_IA(Proc_,iNode)
                if(iProc==iPE)then
                   if ( CutValue> PlotXYZNodes_DNB(2,1,1   ,1,iBlock) .and. &
                        CutValue<=PlotXYZNodes_DNB(2,1,1+nJ,1,iBlock)  )then
                      nBlockCuts=nBlockCuts+1
                      BlockCut(iBlockALL)=nBlockCuts
                   end if
                end if
                call MPI_Bcast(nBlockCuts,1,MPI_Integer,iPE,iComm,iError)
             end do
             if(iProc==0)then
                ! Write file header
                write(UnitTmp_,'(a)')'TITLE="BATSRUS: Cut Y Data, ' &
                     //textDateTime//'"'
                write(UnitTmp_,'(a)')unitstr_TEC(1:len_trim(unitstr_TEC))
                write(UnitTmp_,'(a,a,i12,a,i12,a)') &
                     'ZONE T="2D Y '//textNandT//'"', &
                     ', N=',nBlockCuts*((nI+1)*(nK+1)), &
                     ', E=',nBlockCuts*((nI  )*(nK  )), &
                     ', F=FEPOINT, ET=QUADRILATERAL'
                call write_tecplot_auxdata
             end if
             ! Now loop to write values
             do iBlockAll = 1, nNodeUsed
                iNode = iNodeMorton_I(iBlockAll)
                iBlock  = iTree_IA(Block_,iNode)
                iPE   = iTree_IA(Proc_,iNode)
                if(iProc==iPE)then
                   if ( CutValue> PlotXYZNodes_DNB(2,1,1   ,1,iBlock) .and. &
                        CutValue<=PlotXYZNodes_DNB(2,1,1+nJ,1,iBlock)  )then
                      ! Find cut interpolation factors
                      do j=1,nJ
                         if ( CutValue> PlotXYZNodes_DNB(2,1,j  ,1,iBlock) .and.&
                              CutValue<=PlotXYZNodes_DNB(2,1,j+1,1,iBlock)  )then
                            cut1=j
                            cut2=j+1
                            factor2=(CutValue-PlotXYZNodes_DNB(2,1,j,1,iBlock))/&
                                 ( PlotXYZNodes_DNB(2,1,j+1,1,iBlock) &
                                 - PlotXYZNodes_DNB(2,1,j,1,iBlock))
                            factor1=1.-factor2
                            EXIT
                         end if
                      end do
                      ! Write point values
                      call fill_NodeXYZ
                      do k=1,1+nK; do i=1,1+nI
                         write(UnitTmp_,fmt="(50(ES14.6))") &
                              (factor1*NodeXYZ_DN(1:3,i,cut1,k)+ &
                              factor2*NodeXYZ_DN(1:3,i,cut2,k)), &
                              (factor1*PlotVarNodes_VNB(1:nPlotVar,i,cut1,k,iBlock)+ &
                              factor2*PlotVarNodes_VNB(1:nPlotVar,i,cut2,k,iBlock))
                      end do; end do
                      ! Write point connectivity
                      do k=1,nK; do i=1,nI
                         write(UnitTmp2_,'(4(i8,1x))') &
                              ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) &
                              + (k-1)*(nI+1)+i, &
                              ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) &
                              + (k-1)*(nI+1)+i+1, &
                              ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) &
                              + (k  )*(nI+1)+i+1, &
                              ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) &
                              + (k  )*(nI+1)+i
                      end do; end do
                   end if
                end if
             end do
          else if(IsRLonLat) then
             ! First loop to count nodes and cells
             do iBlockAll = 1, nNodeUsed
                iNode = iNodeMorton_I(iBlockAll)
                iBlock  = iTree_IA(Block_,iNode)
                iPE   = iTree_IA(Proc_,iNode)
                if(iProc==iPE)then
                   if ( (CutValue - PlotXYZNodes_DNB(2,1,1   ,2,iBlock))* &
                        (CutValue - PlotXYZNodes_DNB(2,1,1+nJ,2,iBlock)) <=0)then
                      nBlockCuts=nBlockCuts+1
                      BlockCut(iBlockALL)=nBlockCuts
                   end if
                end if
                call MPI_Bcast(nBlockCuts,1,MPI_Integer,iPE,iComm,iError)
             end do
             if(iProc==0)then
                ! Write file header
                write(UnitTmp_,'(a)')'TITLE="BATSRUS: Cut Y Data, ' &
                     //textDateTime//'"'
                write(UnitTmp_,'(a)')unitstr_TEC(1:len_trim(unitstr_TEC))
                write(UnitTmp_,'(a,a,i12,a,i12,a)') &
                     'ZONE T="2D Y '//textNandT//'"', &
                     ', N=',nBlockCuts*((nI+1)*(nK+1)), &
                     ', E=',nBlockCuts*((nI  )*(nK  )), &
                     ', F=FEPOINT, ET=QUADRILATERAL'
                call write_tecplot_auxdata
             end if
             ! Now loop to write values
             do iBlockAll = 1, nNodeUsed
                iNode = iNodeMorton_I(iBlockAll)
                iBlock  = iTree_IA(Block_,iNode)
                iPE   = iTree_IA(Proc_,iNode)
                if(iProc==iPE)then
                   if ( (CutValue - PlotXYZNodes_DNB(2,1,1   ,2,iBlock))* &
                        (CutValue - PlotXYZNodes_DNB(2,1,1+nJ,2,iBlock)) <=0)then
                      ! Find cut interpolation factors
                      do j=1,nJ
                         if ( (CutValue - PlotXYZNodes_DNB(2,1,j  ,2,iBlock))* &
                              (CutValue - PlotXYZNodes_DNB(2,1,j+1,2,iBlock)) &
                              <= 0)then
                            cut1=j
                            cut2=j+1
                            factor2 = &
                                 (CutValue - PlotXYZNodes_DNB(2,1,j,2,iBlock))/ &
                                 ( PlotXYZNodes_DNB(2,1,j+1,2,iBlock) &
                                 - PlotXYZNodes_DNB(2,1,j,2,iBlock))
                            factor1=1.-factor2
                            EXIT
                         end if
                      end do
                      ! Write point values
                      call fill_NodeXYZ
                      do k=1,1+nK; do i=1,1+nI
                         write(UnitTmp_,fmt="(50(ES14.6))") &
                              (factor1*NodeXYZ_DN(1:3,i,cut1,k)+ &
                              factor2*NodeXYZ_DN(1:3,i,cut2,k)), &
                              (factor1 &
                              *PlotVarNodes_VNB(1:nPlotVar,i,cut1,k,iBlock)+ &
                              factor2 &
                              *PlotVarNodes_VNB(1:nPlotVar,i,cut2,k,iBlock))
                      end do; end do
                      ! Write point connectivity
                      do k=1,nK; do i=1,nI
                         write(UnitTmp2_,'(4(i8,1x))') &
                              ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) &
                              + (k-1)*(nI+1)+i, &
                              ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) &
                              + (k-1)*(nI+1)+i+1, &
                              ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) &
                              + (k  )*(nI+1)+i+1, &
                              ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) &
                              + (k  )*(nI+1)+i
                      end do; end do
                   end if
                end if
             end do
          end if

       else
          ! Z Slice
          CutValue = 0.5*(zmin+zmax)
          if(plot_type1(1:3) == 'z=0') CutValue = 0.
          ! First loop to count nodes and cells
          do iBlockAll = 1, nNodeUsed
             iNode = iNodeMorton_I(iBlockAll)
             iBlock  = iTree_IA(Block_,iNode)
             iPE   = iTree_IA(Proc_,iNode)
             if(iProc==iPE)then
                if ( CutValue> PlotXYZNodes_DNB(3,1,1,1   ,iBlock) .and. &
                     CutValue<=PlotXYZNodes_DNB(3,1,1,1+nK,iBlock)  )then
                   nBlockCuts=nBlockCuts+1
                   BlockCut(iBlockALL)=nBlockCuts
                end if
             end if
             call MPI_Bcast(nBlockCuts,1,MPI_Integer,iPE,iComm,iError)
          end do
          if(iProc==0)then
             ! Write file header
             write(UnitTmp_,'(a)')'TITLE="BATSRUS: Cut Z Data, ' &
                  //textDateTime//'"'
             write(UnitTmp_,'(a)')unitstr_TEC(1:len_trim(unitstr_TEC))
             write(UnitTmp_,'(a,a,i12,a,i12,a)') &
                  'ZONE T="2D Z '//textNandT//'"', &
                  ', N=',nBlockCuts*((nI+1)*(nJ+1)), &
                  ', E=',nBlockCuts*((nI  )*(nJ  )), &
                  ', F=FEPOINT, ET=QUADRILATERAL'
             call write_tecplot_auxdata
          end if
          ! Now loop to write values
          do iBlockAll = 1, nNodeUsed
             iNode = iNodeMorton_I(iBlockAll)
             iBlock  = iTree_IA(Block_,iNode)
             iPE   = iTree_IA(Proc_,iNode)
             if(iProc==iPE)then
                if ( CutValue> PlotXYZNodes_DNB(3,1,1,1   ,iBlock) .and. &
                     CutValue<=PlotXYZNodes_DNB(3,1,1,1+nK,iBlock)  )then
                   ! Find cut interpolation factors
                   do k=1,nK
                      if ( CutValue> PlotXYZNodes_DNB(3,1,1,k  ,iBlock) .and. &
                           CutValue<=PlotXYZNodes_DNB(3,1,1,k+1,iBlock)  )then
                         cut1=k
                         cut2=k+1
                         factor2=(CutValue-PlotXYZNodes_DNB(3,1,1,k,iBlock))/ &
                              ( PlotXYZNodes_DNB(3,1,1,k+1,iBlock) &
                              - PlotXYZNodes_DNB(3,1,1,k,iBlock))
                         factor1=1.-factor2
                         EXIT
                      end if
                   end do
                   ! Write point values
                   call fill_NodeXYZ
                   do j=1,1+nJ; do i=1,1+nI
                      write(UnitTmp_,fmt="(50(ES14.6))") &
                           (factor1*NodeXYZ_DN(1:3,i,j,cut1)+ &
                           factor2*NodeXYZ_DN(1:3,i,j,cut2)), &
                           (factor1 &
                           *PlotVarNodes_VNB(1:nPlotVar,i,j,cut1,iBlock)+ &
                           factor2 &
                           *PlotVarNodes_VNB(1:nPlotVar,i,j,cut2,iBlock))
                   end do; end do
                   ! Write point connectivity
                   do j=1,nJ; do i=1,nI
                      write(UnitTmp2_,'(4(i8,1x))') &
                           ((BlockCut(iBlockALL)-1)*(nI+1)*(nJ+1)) &
                           + (j-1)*(nI+1)+i, &
                           ((BlockCut(iBlockALL)-1)*(nI+1)*(nJ+1)) &
                           + (j-1)*(nI+1)+i+1, &
                           ((BlockCut(iBlockALL)-1)*(nI+1)*(nJ+1)) &
                           + (j  )*(nI+1)+i+1, &
                           ((BlockCut(iBlockALL)-1)*(nI+1)*(nJ+1)) &
                           + (j  )*(nI+1)+i
                   end do; end do
                end if
             end if
          end do
       end if
       deallocate(BlockCut)
    case('slc','dpl')
       !================================ arbitrary slices ===============
       okdebug=.false.

       ! XarbP,YarbP,ZarbP                    point on plane
       ! XarbNormal,YarbNormal,ZarbNormal     normal for cut
       ! ic1,jc1,kc1,ic2,jc2,kc2              two opposite corner indices

       if (plot_type1(1:3)=='slc')then
          ! Point-Normal cut plot
          XarbP=plot_point(1,ifile); XarbNormal=plot_normal(1,ifile)
          YarbP=plot_point(2,ifile); YarbNormal=plot_normal(2,ifile)
          ZarbP=plot_point(3,ifile); ZarbNormal=plot_normal(3,ifile)
       else
          ! Dipole cut plot
          XarbP=0.; XarbNormal=-sin(ThetaTilt)
          YarbP=0.; YarbNormal=0.
          ZarbP=0.; ZarbNormal= cos(ThetaTilt)
       end if

       ! First loop to count cuts
       nBlockCuts=0
       do iBlock = 1, nBlock
          if(iTypeAdvance_B(iBlock) == SkippedBlock_) CYCLE
          ic1=1; ic2=1+nI
          jc1=1; jc2=1+nJ
          kc1=1; kc2=1+nK
          call find_cuts(-1)
          if ( nCuts>0 )then
             ! count up number of cuts
             do i=1,nI; do j=1,nJ; do k=1,nK
                ic1=i; ic2=i+1
                jc1=j; jc2=j+1
                kc1=k; kc2=k+1
                call find_cuts(0)
                nBlockCuts=nBlockCuts+nCuts
             end do; end do; end do
          end if
       end do
       call MPI_reduce(nBlockCuts, nCutsTotal, 1, MPI_INTEGER, MPI_SUM, 0, &
            iComm, iError)

       ! Write file header
       if(iProc==0)then
          if (plot_type1(1:3)=='slc')then
             write(UnitTmp_,'(a)')'TITLE="BATSRUS: Slice, '//textDateTime//'"'
             write(UnitTmp_,'(a)')trim(unitstr_TEC)
             write(UnitTmp_,'(a,i8,a)') &
                  'ZONE T="Slice '//textNandT//'", I=', nCutsTotal,&
                  ', J=1, K=1, F=POINT'
          else
             write(UnitTmp_,'(a)')'TITLE="BATSRUS: Dipole Cut, '// &
                  textDateTime//'"'
             write(UnitTmp_,'(a)')trim(unitstr_TEC)
             write(UnitTmp_,'(a,i8,a)') &
                  'ZONE T="Dipole Cut '//textNandT//'", I=', nCutsTotal,&
                  ', J=1, K=1, F=POINT'
          end if
          call write_tecplot_auxdata
       end if

       ! Now loop to write values
       do iBlock = 1, nBlock
          if(iTypeAdvance_B(iBlock) == SkippedBlock_) CYCLE
          ic1=1; ic2=1+nI
          jc1=1; jc2=1+nJ
          kc1=1; kc2=1+nK
          call find_cuts(-1)
          if ( nCuts>0 )then
             ! write the cuts
             call fill_NodeXYZ
             do i=1,nI; do j=1,nJ; do k=1,nK
                ic1=i; ic2=i+1
                jc1=j; jc2=j+1
                kc1=k; kc2=k+1
                call find_cuts(1)
             end do; end do; end do
          end if
       end do
    case default
       write(*,*) NameSub,': ERROR: Unknown plot_type='//plot_type1
    end select

    call test_stop(NameSub, DoTest)
  contains
    !==========================================================================

    ! Assumes iBlock value is correct
    subroutine fill_nodeXYZ
      ! Fill array with position (optionally dimensioned)
      !------------------------------------------------------------------------
      if (plot_dimensional(ifile)) then
         NodeXYZ_DN(1:3,:,:,:)=PlotXYZNodes_DNB(1:3,:,:,:,iBlock)*No2Io_V(UnitX_)
      else
         NodeXYZ_DN(1:3,:,:,:)=PlotXYZNodes_DNB(1:3,:,:,:,iBlock)
      end if
    end subroutine fill_nodeXYZ
    !==========================================================================

    ! iopt =-1 check all edges to see if cut
    !      = 0 count cuts only
    !      = 1 find cuts and write to disk
    subroutine find_cuts(iopt)
      integer, intent(in) :: iopt
      integer :: ic,jc,kc

      !------------------------------------------------------------------------
      nCuts=0

      ! Check edges.
      ! X edges
      if (XarbNormal>0.01) then
         ic=ic1; jc=jc1; kc=kc1
         do jc=jc1,jc2,jc2-jc1; do kc=kc1,kc2,kc2-kc1
            if(iopt>-1 .and. (jc==jc1 .or. kc==kc1) .and. (jc/=0 .and. kc/=0))&
                 CYCLE
            Yp = PlotXYZNodes_DNB(2,ic,jc,kc,iBlock)
            Zp = PlotXYZNodes_DNB(3,ic,jc,kc,iBlock)
            Xp = XarbP &
                 - (YarbNormal*(Yp-YarbP) + ZarbNormal*(Zp-ZarbP))/XarbNormal
            if ( Xp> PlotXYZNodes_DNB(1,ic1,jc,kc,iBlock) .and. &
                 Xp<=PlotXYZNodes_DNB(1,ic2,jc,kc,iBlock) )then
               if(okdebug)write(*,*)'x-cut:',iopt,Xp,Yp,Zp
               if(iopt==-1)then
                  nCuts=1; RETURN
               end if
               ! Cycle if outside of clipping box
               if ( Xp<xmin .or. Yp<ymin .or. Zp<zmin .or. &
                    Xp>xmax .or. Yp>ymax .or. Zp>zmax) CYCLE
               nCuts=nCuts+1
               if (iopt>0) then
                  ! Write point values
                  factor2 = (Xp-PlotXYZNodes_DNB(1,ic1,jc,kc,iBlock))/ &
                       ( PlotXYZNodes_DNB(1,ic2,jc,kc,iBlock) &
                       - PlotXYZNodes_DNB(1,ic1,jc,kc,iBlock))
                  factor1=1.-factor2
                  write(UnitTmp_,fmt="(50(ES14.6))") &
                       (factor1*NodeXYZ_DN(:, ic1,jc,kc)+ &
                       factor2*NodeXYZ_DN(:, ic2,jc,kc)), &
                       (factor1*PlotVarNodes_VNB(1:nPlotVar,ic1,jc,kc,iBlock)+ &
                       factor2*PlotVarNodes_VNB(1:nPlotVar,ic2,jc,kc,iBlock))
                  if(okdebug)write(*,*)'  i=',ic1,'-',ic2,' j=',jc,' k=',kc
               end if
            end if
         end do; end do
      end if
      ! Y edges
      if (YarbNormal>0.01) then
         ic=ic1; jc=jc1; kc=kc1
         do ic=ic1,ic2,ic2-ic1; do kc=kc1,kc2,kc2-kc1
            if(iopt>-1 .and. (ic==ic1 .or. kc==kc1) .and. (ic/=0 .and. kc/=0))&
                 CYCLE
            Xp = PlotXYZNodes_DNB(1,ic,jc,kc,iBlock)
            Zp = PlotXYZNodes_DNB(3,ic,jc,kc,iBlock)
            Yp = YarbP &
                 - (XarbNormal*(Xp-XarbP) + ZarbNormal*(Zp-ZarbP))/YarbNormal
            if ( Yp> PlotXYZNodes_DNB(2,ic,jc1,kc,iBlock) .and. &
                 Yp<=PlotXYZNodes_DNB(2,ic,jc2,kc,iBlock) )then
               if(okdebug)write(*,*)'y-cut:',iopt,Xp,Yp,Zp
               if(iopt==-1)then
                  nCuts=1; RETURN
               end if
               ! Cycle if outside of clipping box
               if ( Xp<xmin .or. Yp<ymin .or. Zp<zmin .or. &
                    Xp>xmax .or. Yp>ymax .or. Zp>zmax) CYCLE
               nCuts=nCuts+1
               if (iopt>0) then
                  ! Write point values
                  factor2=(Yp-PlotXYZNodes_DNB(2,ic,jc1,kc,iBlock))/ &
                       ( PlotXYZNodes_DNB(2,ic,jc2,kc,iBlock) &
                       - PlotXYZNodes_DNB(2,ic,jc1,kc,iBlock))
                  factor1 = 1 - factor2
                  write(UnitTmp_,fmt="(50(ES14.6))") &
                       (factor1*NodeXYZ_DN(:, ic,jc1,kc)+ &
                       factor2*NodeXYZ_DN(:, ic,jc2,kc)), &
                       (factor1*PlotVarNodes_VNB(1:nPlotVar,ic,jc1,kc,iBlock)+ &
                       factor2*PlotVarNodes_VNB(1:nPlotVar,ic,jc2,kc,iBlock))
                  if(okdebug)write(*,*)'  i=',ic,' j=',jc1,'-',jc2,' k=',kc
               end if
            end if
         end do; end do
      end if
      ! Z edges
      if (ZarbNormal>0.01) then
         ic=ic1; jc=jc1; kc=kc1
         do ic=ic1,ic2,ic2-ic1; do jc=jc1,jc2,jc2-jc1
            if(iopt>-1 .and. (ic==ic1 .or. jc==jc1) .and. (ic/=0 .and. jc/=0))&
                 CYCLE
            Xp = PlotXYZNodes_DNB(1,ic,jc,kc,iBlock)
            Yp = PlotXYZNodes_DNB(2,ic,jc,kc,iBlock)
            Zp = ZarbP &
                 - (XarbNormal*(Xp-XarbP) + YarbNormal*(Yp-YarbP))/ZarbNormal
            if ( Zp> PlotXYZNodes_DNB(3,ic,jc,kc1,iBlock) .and. &
                 Zp<=PlotXYZNodes_DNB(3,ic,jc,kc2,iBlock) )then
               if(okdebug)write(*,*)'z-cut:',iopt,Xp,Yp,Zp
               if(iopt==-1)then
                  nCuts=1; RETURN
               end if
               ! Cycle if outside of clipping box
               if ( Xp<xmin .or. Yp<ymin .or. Zp<zmin .or. &
                    Xp>xmax .or. Yp>ymax .or. Zp>zmax) CYCLE
               nCuts=nCuts+1
               if (iopt>0) then
                  ! Write point values
                  factor2 = (Zp - PlotXYZNodes_DNB(3,ic,jc,kc1,iBlock))/ &
                       ( PlotXYZNodes_DNB(3,ic,jc,kc2,iBlock) &
                       - PlotXYZNodes_DNB(3,ic,jc,kc1,iBlock))
                  factor1=1.-factor2
                  write(UnitTmp_,fmt="(50(ES14.6))") &
                       (factor1*NodeXYZ_DN(:, ic,jc,kc1)+ &
                       factor2*NodeXYZ_DN(:, ic,jc,kc2)), &
                       (factor1*PlotVarNodes_VNB(1:nPlotVar,ic,jc,kc1,iBlock)+ &
                       factor2*PlotVarNodes_VNB(1:nPlotVar,ic,jc,kc2,iBlock))
                  if(okdebug)write(*,*)'  i=',ic,' j=',jc,' k=',kc1,'-',kc2
               end if
            end if
         end do; end do
      end if

    end subroutine find_cuts
    !==========================================================================

  end subroutine write_tecplot_node_data
  !============================================================================

  subroutine assign_node_numbers

    use ModProcMH
    use ModIO, ONLY: write_prefix, iUnitOut
    use ModMain, ONLY: nBlock, nBlockMax, nBlockALL
    use ModAdvance,  ONLY: iTypeAdvance_B, iTypeAdvance_BP, SkippedBlock_
    use ModNodes
    use ModMpi
    use BATL_lib, ONLY: message_pass_node

    integer, parameter :: NodesPerBlock=(nI+1)*(nJ+1)*(nK+1)
    integer :: iBlockStart
    integer :: i, j, k, iNode, iBlock, iError, iPE, iTag
    integer :: nOffset, nOffsetPrevious
    integer, allocatable, dimension(:) :: NodeOffset, NodeOffsetMax, nOffset_P
    real, allocatable, dimension(:,:,:,:,:) :: IndexNode_VNB
    logical :: DoAllReduce=.true.
    integer :: iStatus(MPI_STATUS_SIZE)

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'assign_node_numbers'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    ! Write information to the screen
    if(iProc==0.and.lVerbose>0)then
       call write_prefix; write(iUnitOut,*)'Starting assign_node_numbers ...'
    end if

    ! Initialize all node numbers to zero
    NodeNumberLocal_NB=0

    ! Number of nodes on each block (maximum)
    nNodeALL=nBlockALL*NodesPerBlock

    ! Count number of used blocks on all processors with rank below this one
    iBlockStart = 0
    if(iProc > 0) iBlockStart = &
         count(iTypeAdvance_BP(1:nBlockMax,0:iProc-1) /= SkippedBlock_)

    iNode = iBlockStart*NodesPerBlock

    ! Loop to assign local and global node numbers
    TREE1: do iBlock  = 1, nBlock
       if(iTypeAdvance_B(iBlock) == SkippedBlock_) CYCLE
       do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
          iNode = iNode+1
          NodeNumberLocal_NB(i,j,k,iBlock)= iNode
       end do; end do; end do
    end do TREE1
    NodeNumberGlobal_NB = NodeNumberLocal_NB

    ! Set logical array
    NodeUniqueGlobal_NB = NodeNumberGlobal_NB>0

    ! Assign value to internal passing variable and do message pass
    !  NOTE: convert integer to real for message pass first

    ! Done a evel one, with allocate and dealocate. NEED to be fixed
    allocate(IndexNode_VNB(1,nI+1,nJ+1,nK+1,nBlock))
    IndexNode_VNB(1,:,:,:,:) = real(NodeNumberGlobal_NB(:,:,:,1:nBlock))
    call message_pass_node(1,IndexNode_VNB, &
         NameOperatorIn='Min', UsePeriodicCoordIn = .true.)
    NodeNumberGlobal_NB(:,:,:,1:nBlock) = nint(IndexNode_VNB(1,:,:,:,:))
    deallocate(IndexNode_VNB)

    ! Allocate memory for storing the node offsets
    allocate( NodeOffset   (nBlockALL*NodesPerBlock))
    allocate( NodeOffsetMax(nBlockALL*NodesPerBlock))
    NodeOffset=0

    ! Loop to compute node offsets
    nOffset=0
    TREE2: do iBlock  = 1, nBlock
       if(iTypeAdvance_B(iBlock) == SkippedBlock_) CYCLE
       do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
          if(  NodeNumberLocal_NB(i,j,k,iBlock) > &
               NodeNumberGlobal_NB(i,j,k,iBlock))then
             nOffset = nOffset+1
             NodeUniqueGlobal_NB(i,j,k,iBlock) = .false.
          end if
          NodeOffset(NodeNumberLocal_NB(i,j,k,iBlock)) = nOffset
       end do; end do; end do
    end do TREE2

    ! Collect offsets from all the PEs
    allocate(nOffset_P(0:nProc-1))
    call MPI_allgather(nOffset, 1, MPI_INTEGER, nOffset_P, 1, MPI_INTEGER, &
         iComm, iError)

    ! Add up the offsets on processors with lower rank
    nOffsetPrevious = 0
    if(iProc > 0) nOffsetPrevious = sum(nOffset_P(0:iProc-1))

    ! Increase the offset on this processor by nOffsetPrevious
    do iBlock  = 1, nBlock
       if(iTypeAdvance_B(iBlock) == SkippedBlock_) CYCLE
       do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
          iNode = NodeNumberLocal_NB(i,j,k,iBlock)
          NodeOffset(iNode) = NodeOffset(iNode) + nOffsetPrevious
       end do; end do; end do
    end do

    ! Gather offsets from all PE-s.
    ! NodeOffset was initialized to 0 so MPI_MAX works.
    if(DoAllReduce)then
       call MPI_allreduce(NodeOffset,NodeOffsetMax,nBlockALL*NodesPerBlock, &
            MPI_INTEGER,MPI_MAX,iComm,iError)
       NodeOffset = NodeOffsetMax
       nNodeALL   = nNodeALL - sum(nOffset_P)
    else
       if(iProc == 0) then
          do iPE=1,nProc-1
             iTag = iPE
             call MPI_recv(NodeOffsetMax,nBlockALL*NodesPerBlock, &
                  MPI_INTEGER,iPE,itag,iComm,iStatus,iError)
             NodeOffset = max(NodeOffset,NodeOffsetMax)
          end do
       else
          itag = iProc
          call MPI_send(NodeOffset,nBlockALL*NodesPerBlock, &
               MPI_INTEGER,0,itag,iComm,iError)
       end if
       call MPI_Bcast(NodeOffset, nBlockALL*NodesPerBlock, MPI_INTEGER, 0, &
            iComm, iError)
    end if

    ! Loop to fix NodeNumberGlobal_NB for offset
    TREE3: do iBlock  = 1, nBlock
       if(iTypeAdvance_B(iBlock) == SkippedBlock_) CYCLE
       do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
          NodeNumberGlobal_NB(i,j,k,iBlock) = NodeNumberGlobal_NB(i,j,k,iBlock) &
               - NodeOffset(NodeNumberGlobal_NB(i,j,k,iBlock))
          if(NodeNumberGlobal_NB(i,j,k,iBlock)>nNodeALL &
               .or. NodeNumberGlobal_NB(i,j,k,iBlock)<1)then
             ! Error in numbering, report values and stop.
             write(*,*)'ERROR: Global node numbering problem.', &
                  ' PE=',iProc,' BLK=',iBlock,' ijk=',i,j,k
             write(*,*)'  NodeNumberGlobal_NB=',&
                  NodeNumberGlobal_NB(i,j,k,iBlock)
             write(*,*)'  NodeOffset           =',&
                  NodeOffset(NodeNumberGlobal_NB(i,j,k,iBlock))
             write(*,*)'  nBlockALL=',nBlockALL,&
                  ' NodesPerBlock=',NodesPerBlock,&
                  ' unreduced total=',nBlockALL*NodesPerBlock,&
                  ' nNodeALL=',nNodeALL
             call stop_mpi('message_pass_nodes: error in numbering')
          end if
       end do; end do; end do
    end do TREE3

    ! Deallocate memory when done with it
    deallocate(NodeOffset, NodeOffsetMax, nOffset_P)

    ! Write information to the screen
    if(iProc==0)then
       call write_prefix; write(iUnitOUt,*) &
            ' nBlockALL=',nBlockALL,' NodesPerBlock=',NodesPerBlock, &
            ' unreduced total=',nBlockALL*NodesPerBlock,' nNodeALL=',nNodeALL
    end if

    call test_stop(NameSub, DoTest)
  end subroutine assign_node_numbers
  !============================================================================

  subroutine set_tecplot_var_string( &
       iFile, nPlotVar, NamePlotVar_V, StringVarTec)

    ! Set StringVarTec with variable names and units
    ! based on NamePlotVar_V array

    use ModPhysics
    use ModUtilities,  ONLY: lower_case
    use ModIO,         ONLY: plot_dimensional, plot_type1
    use ModVarIndexes, ONLY: IsMhd
    use ModIO,         ONLY: NameVarUserTec_I, NameUnitUserTec_I
    use ModMultiFluid, ONLY: extract_fluid_name, iFluid, NameFluid
    use BATL_lib,      ONLY: nDim

    ! Arguments

    integer, intent(in)              :: nPlotVar, iFile
    character (len=*), intent(in)    :: NamePlotVar_V(nPlotVar)
    character (len=1500), intent(out) :: StringVarTec

    character (len=20) :: NameTecFluid
    character (len=20) :: String, NamePlotVar, NameTecVar, NameUnit
    integer            :: iPlotVar, iVar

    ! Coordinate names and units
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'set_tecplot_var_string'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(plot_type1(1:3) == 'box')then
       if (plot_dimensional(iFile)) then
          StringVarTec = 'VARIABLES ="R '// trim(NameTecUnit_V(UnitX_)) &
               // '", "", "'
       else
          StringVarTec = 'VARIABLES ="x", "y", "z'
       end if

    elseif(plot_type1(1:3) == 'shl')then
       if (plot_dimensional(iFile)) then
          StringVarTec = 'VARIABLES ="R '// trim(NameTecUnit_V(UnitX_)) &
               // '", "Lon [deg]", "Lat [deg]'
       else
          StringVarTec = 'VARIABLES ="R", "Lon", "Lat'
       end if

    else
       if (plot_dimensional(iFile)) then
          StringVarTec = 'VARIABLES ="X ' // trim(NameTecUnit_V(UnitX_)) &
               // '", "Y ' // trim(NameTecUnit_V(UnitX_))
          if(nDim==3) StringVarTec = trim(StringVarTec) &
               // '", "Z ' // trim(NameTecUnit_V(UnitX_))
       else
          if(nDim==2) StringVarTec = 'VARIABLES = "X", "Y'
          if(nDim==3) StringVarTec = 'VARIABLES = "X", "Y", "Z'
       end if

    end if

    do iPlotVar = 1, nPlotVar

       NamePlotVar = NamePlotVar_V(iPlotVar)
       call lower_case(NamePlotVar)
       String = NamePlotVar
       call extract_fluid_name(String)
       NameTecFluid = ''
       if(iFluid > 1 .or. .not. IsMhd) NameTecFluid = '^'//NameFluid

       ! Default value for NameUnit is empty string
       NameUnit = ''

       select case(String)
       case('rho')
          NameTecVar = 'Rho'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitRho_)
       case('rhoux','mx')
          NameTecVar = 'Rho U_x'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitRhoU_)
       case('rhouy','my')
          NameTecVar = 'Rho U_y'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitRhoU_)
       case('rhouz','mz')
          NameTecVar = 'Rho U_z'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitRhoU_)
       case('bx')
          NameTecVar = 'B_x'
          NameUnit   = NameTecUnit_V(UnitB_)
       case('by')
          NameTecVar = 'B_y'
          NameUnit   = NameTecUnit_V(UnitB_)
       case('bz')
          NameTecVar = 'B_z'
          NameUnit   = NameTecUnit_V(UnitB_)
          ! face centered magnetic field
       case('bxl') ! east
          NameTecVar = 'B_e'
          NameUnit   = NameTecUnit_V(UnitB_)
       case('bxr') ! west
          NameTecVar = 'B_w'
          NameUnit   = NameTecUnit_V(UnitB_)
       case('byl') ! south
          NameTecVar = 'B_s'
          NameUnit   = NameTecUnit_V(UnitB_)
       case('byr') ! north
          NameTecVar = 'B_n'
          NameUnit   = NameTecUnit_V(UnitB_)
       case('bzl') ! bottom
          NameTecVar = 'B_b'
          NameUnit   = NameTecUnit_V(UnitB_)
       case('bzr') ! top
          NameTecVar = 'B_t'
          NameUnit   = NameTecUnit_V(UnitB_)
          !
       case('e')
          NameTecVar = 'E'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitEnergydens_)
       case('p','pth')
          NameTecVar = 'P'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitP_)
       case('ppar')
          NameTecVar = 'P_par'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitP_)
       case('pperp')
          NameTecVar = 'P_perp'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitP_)
       case('n')
          NameTecVar = 'n'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitN_)
       case('t','temp')
          NameTecVar = 'T'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitTemperature_)
       case('ux')
          NameTecVar = 'U_x'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitU_)
       case('uy')
          NameTecVar = 'U_y'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitU_)
       case('uz', 'uzrot')
          NameTecVar = 'U_z'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitU_)
       case('uxrot')
          NameTecVar = 'U_x_rot'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitU_)
       case('uyrot')
          NameTecVar = 'U_y_rot'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitU_)
       case('ur')
          NameTecVar = 'U_r'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitU_)
       case('rhour','mr')
          NameTecVar = 'Rho U_r'//NameTecFluid
          NameUnit   = NameTecUnit_V(UnitRhoU_)
       case('br')
          NameTecVar = 'B_r'
          NameUnit   = NameTecUnit_V(UnitB_)
       case('b1x')
          NameTecVar = 'B1_x'
          NameUnit   = NameTecUnit_V(UnitB_)
       case('b1y')
          NameTecVar = 'B1_y'
          NameUnit   = NameTecUnit_V(UnitB_)
       case('b1z')
          NameTecVar = 'B1_z'
          NameUnit   = NameTecUnit_V(UnitB_)
       case('b1r')
          NameTecVar = 'B1_r'
          NameUnit   = NameTecUnit_V(UnitB_)
       case('jx')
          NameTecVar = 'J_x'
          NameUnit   = NameTecUnit_V(UnitJ_)
       case('jy')
          NameTecVar = 'J_y'
          NameUnit   = NameTecUnit_V(UnitJ_)
       case('jz')
          NameTecVar = 'J_z'
          NameUnit   = NameTecUnit_V(UnitJ_)
       case('jr')
          NameTecVar = 'J_r'
          NameUnit   = NameTecUnit_V(UnitJ_)
       case('ex')
          NameTecVar = 'E_x'
          NameUnit   = NameTecUnit_V(UnitElectric_)
       case('ey')
          NameTecVar = 'E_y'
          NameUnit   = NameTecUnit_V(UnitElectric_)
       case('ez')
          NameTecVar = 'E_z'
          NameUnit   = NameTecUnit_V(UnitElectric_)
       case('expot')
          NameTecVar = 'Epot_x'
          NameUnit   = NameTecUnit_V(UnitElectric_)
       case('eypot')
          NameTecVar = 'Epot_y'
          NameUnit   = NameTecUnit_V(UnitElectric_)
       case('ezpot')
          NameTecVar = 'Epot_z'
          NameUnit   = NameTecUnit_V(UnitElectric_)
       case('exind')
          NameTecVar = 'Eind_x'
          NameUnit   = NameTecUnit_V(UnitElectric_)
       case('eyind')
          NameTecVar = 'Eind_y'
          NameUnit   = NameTecUnit_V(UnitElectric_)
       case('ezind')
          NameTecVar = 'Eind_z'
          NameUnit   = NameTecUnit_V(UnitElectric_)
       case('gradpex')
          NameTecVar = 'GradPe_x'
          NameUnit   = 'nPa/m'
       case('gradpey')
          NameTecVar = 'GradPe_y'
          NameUnit   = 'nPa/m'
       case('gradpez')
          NameTecVar = 'GradPe_z'
          NameUnit   = 'nPa/m'
       case('gradper')
          NameTecVar = 'GradPe_r'
          NameUnit   = 'nPa/m'
       case('pote')
          NameTecVar = 'PotE'
          NameUnit   = 'V'
       case('dive')
          NameTecVar = 'div(E)'
          NameUnit   = 'V/m^2'
       case('er')
          NameTecVar = 'E_r'
          NameUnit   = NameTecUnit_V(UnitElectric_)
       case('pvecx')
          NameTecVar = 'S_x'
          NameUnit   = NameTecUnit_V(UnitPoynting_)
       case('pvecy')
          NameTecVar = 'S_y'
          NameUnit   = NameTecUnit_V(UnitPoynting_)
       case('pvecz')
          NameTecVar = 'S_z'
          NameUnit   = NameTecUnit_V(UnitPoynting_)
       case('pvecr')
          NameTecVar = 'S_r'
          NameUnit   = NameTecUnit_V(UnitPoynting_)
       case('b2ur')
          NameTecVar = 'B^2/mu_0 U_r'
          NameUnit   = NameTecUnit_V(UnitPoynting_)
       case('divb', 'divb_cd', 'divb_ct', 'absdivb')
          NameTecVar = 'div B'
          NameUnit   = NameTecUnit_V(UnitDivB_)
       case('theta1')
          NameTecVar = 'theta_1'
          NameUnit   = NameTecUnit_V(UnitAngle_)
       case('phi1')
          NameTecVar = 'phi_1'
          NameUnit   = NameTecUnit_V(UnitAngle_)
       case('theta2')
          NameTecVar = 'theta_2'
          NameUnit   = NameTecUnit_V(UnitAngle_)
       case('phi2')
          NameTecVar = 'phi_2'
          NameUnit   = NameTecUnit_V(UnitAngle_)
       case('status')
          NameTecVar = 'Status'
       case('f1x','f1y','f1z','f2x','f2y','f2z')
          NameTecVar = NamePlotVar
       case('x','y','z','r','dx','dy','dz','req1','req2')
          NameTecVar = String
          NameUnit   = NameTecUnit_V(UnitX_)
       case('dvol')
          NameTecVar = 'dvol'
          NameUnit   = trim(NameTecUnit_V(UnitX_))//'^3'
       case('dt')
          NameTecVar = 'dt'
          NameUnit   = NameTecUnit_V(UnitT_)
       case('dtblk')
          NameTecVar = 'dtblk'
          NameUnit   = NameTecUnit_V(UnitT_)
       case('impl')
          NameTecVar = 'impl'
       case('proc')
          NameTecVar = 'PE #'
       case('blk')
          NameTecVar = 'Block #'
       case('node')
          NameTecVar = 'Node #'
       case('ew','erad')
          NameTecVar = String
          NameUnit   = NameTecUnit_V(UnitEnergydens_)
       case default
          ! Set the default or user defined values
          NameTecVar = NameVarUserTec_I(iPlotVar)
          NameUnit   = NameUnitUserTec_I(iPlotVar)

          ! Try to find the plot variable among the basic variables
          do iVar = 1, nVar
             if(NamePlotVar /= NameVarLower_V(iVar)) CYCLE
             NameUnit = NameUnitUserTec_V(iVar)
             EXIT
          end do
       end select

       StringVarTec = trim(StringVarTec) // '", "' // NameTecVar

       if (plot_dimensional(iFile)) &
            StringVarTec = trim(StringVarTec) // ' ' //NameUnit

    end do

    ! Append a closing double quote
    StringVarTec = trim(StringVarTec) // '"'

    call test_stop(NameSub, DoTest)
  end subroutine set_tecplot_var_string
  !============================================================================

end module ModWriteTecplot
!==============================================================================
