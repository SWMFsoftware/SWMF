!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine write_plot_tec(iFile, nPlotVar, PlotVarBlk, PlotVarNodes_VNB, &
     PlotXYZNodes_DNB, unitstr_TEC, xmin, xmax, ymin, ymax, zmin, zmax,  &
     iUnit)

  !NOTE: This routine assumes that the blocks are sorted on PEs by their global
  !       block number, ie blocks 1 to n on PE 0, blocks n+1 to n+m on PE 1,
  !       etc.

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
  use ModWriteTecplot, ONLY: textDateTime, textNandT, CharNewLine, &
       write_tecplot_setinfo, write_tecplot_auxdata

  implicit none

  ! Arguments  
  integer, intent(in) :: ifile, nPlotVar
  character (LEN=1000), intent(in) :: unitstr_TEC
  real, intent(in) :: PlotVarBLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nPlotVarMax)
  real, intent(in) :: PlotVarNodes_VNB(nPlotVarMax,nI+1,nJ+1,nK+1,nBlock)
  real, intent(in) :: PlotXYZNodes_DNB(3,nI+1,nJ+1,nK+1,nBlock)
  real, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax
  integer, intent(in) :: iUnit

  ! Local Variables
  integer :: i,j,k, cut1,cut2, iPE,iBLK, iBlockAll, iNode, nBlockCuts, iError
  real :: CutValue, factor1,factor2
  logical :: oktest,oktest_me
  integer, allocatable, dimension(:) :: BlockCut
  character(len=500) :: stmp

  ! parameters for saving 3d tecplot in a single file
  character (len=80) :: formatData
  integer :: iRec

  integer::ic1,ic2,jc1,jc2,kc1,kc2, nCuts, nCutsTotal
  real :: XarbP,YarbP,ZarbP, XarbNormal,YarbNormal,ZarbNormal, Xp,Yp,Zp
  real, dimension(3,1:nI+1,1:nJ+1,1:nK+1) :: NodeXYZ_DN
  logical :: okdebug
  !----------------------------------------------------------------------------

  call set_oktest('write_plot_tec',oktest,oktest_me)
  if(oktest_me)write(*,*) plot_type1,plot_type1(1:3)  

  call write_tecplot_setinfo

  select case(plot_type1(1:3))
  case('blk')
     call find_grid_block(plot_point(:,iFile), iPE, iBlk)
     if(iPE /= iProc) RETURN
           
     write(UnitTmp_,'(a)')'TITLE="BATSRUS: BLK Only, '//textDateTime//'"'
     write(UnitTmp_,'(a)')trim(unitstr_TEC)
     write(UnitTmp_,'(a,i8,a,i8,a,i8,a)') &
          'ZONE T="BLK Only '//textNandT//'", I=',nI+4,&
          ', J=',nJ+4,', K=',nK+4,', F=POINT'
     call write_tecplot_auxdata
     !DEBUGBLK
     write(stmp,'(i12)')iBLK
     write(UnitTmp_,'(a,a,a)') 'AUXDATA DEBUGBLK="',trim(adjustl(stmp)),'"'
     !DEBUGPROC
     write(stmp,'(i12)')iProc
     write(UnitTmp_,'(a,a,a)') 'AUXDATA DEBUGPROC="',trim(adjustl(stmp)),'"'
     ! Write cell values
     do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
        if (plot_dimensional(ifile)) then
           write(UnitTmp_,fmt="(50(ES14.6))") &
                Xyz_DGB(:,i,j,k,iBLK)*No2Io_V(UnitX_), &
                PlotVarBlk(i,j,k,1:nPlotVar)
        else
           write(UnitTmp_,fmt="(50(ES14.6))") &
                Xyz_DGB(:,i,j,k,iBLK), &
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
        iBlk  = iTree_IA(Block_,iNode)
        iPE   = iTree_IA(Proc_,iNode)
        if(iProc==iPE)then
           ! Write point values
           call fill_NodeXYZ
           write(UnitTmp_,fmt="(50(ES14.6))") &
                (NodeXYZ_DN(1,1,1,1)+NodeXYZ_DN(1,nI+1,1,1))/2., &
                (NodeXYZ_DN(2,1,1,1)+NodeXYZ_DN(2,1,nJ+1,1))/2., &
                (NodeXYZ_DN(3,1,1,1)+NodeXYZ_DN(3,1,1,nK+1))/2., &
                PlotVarNodes_VNB(1:nPlotVar,2,2,2,iBLK)
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
        ! output format for tecplot data, including coordinate info + nPlotVar
        write(formatData, '(a,i2.2,a)') "(", nPlotVar+3, "(ES14.6), a)"
        if(oktest_me)write(*,*) 'formatData =', formatData
        do iBLK = 1, nBlock
           if(iTypeAdvance_B(iBlk) == SkippedBlock_) CYCLE
           call fill_NodeXYZ

           ! Write point values
           do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
              if(NodeUniqueGlobal_NB(i,j,k,iBLK))then
                 iRec = NodeNumberGlobal_NB(i,j,k,iBLK)
                 write(UnitTmp_, FMT=formatData, REC=iRec) &
                      NodeXYZ_DN(1:3,i,j,k),                  &
                      PlotVarNodes_VNB(1:nPlotVar,i,j,k,iBLK),&
                      CharNewLine
              end if
           end do; end do; end do

           ! Write point connectivity.
           ! Initalize record index based on Morton ordering
           iRec = nIJK*(iMortonNode_A(iNode_B(iBLK)) - 1)
           do k=1,nK; do j=1,nJ; do i=1,nI
              iRec = iRec + 1
              write(UnitTmp2_,'(8i11,a)', REC=iRec) &
                   NodeNumberGlobal_NB(i  ,j  ,k  ,iBLK), &
                   NodeNumberGlobal_NB(i+1,j  ,k  ,iBLK), &
                   NodeNumberGlobal_NB(i+1,j+1,k  ,iBLK), &
                   NodeNumberGlobal_NB(i  ,j+1,k  ,iBLK), &
                   NodeNumberGlobal_NB(i  ,j  ,k+1,iBLK), &
                   NodeNumberGlobal_NB(i+1,j  ,k+1,iBLK), &
                   NodeNumberGlobal_NB(i+1,j+1,k+1,iBLK), &
                   NodeNumberGlobal_NB(i  ,j+1,k+1,iBLK), CharNewLine
           end do; end do; end do
        end do
     else
        do iBLK = 1, nBlock
           if(iTypeAdvance_B(iBlk) == SkippedBlock_) CYCLE
           ! Write point values
           call fill_NodeXYZ
           do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
              if(NodeUniqueGlobal_NB(i,j,k,iBLK))then
                 write(UnitTmp_,fmt="(50(ES14.6))") &
                      NodeXYZ_DN(1:3,i,j,k),       &
                      PlotVarNodes_VNB(1:nPlotVar,i,j,k,iBLK)
              end if
           end do; end do; end do
           ! Write point connectivity
           do k=1,nK; do j=1,nJ; do i=1,nI
              write(UnitTmp2_,'(8i11)') &
                   NodeNumberGlobal_NB(i  ,j  ,k  ,iBLK), &
                   NodeNumberGlobal_NB(i+1,j  ,k  ,iBLK), &
                   NodeNumberGlobal_NB(i+1,j+1,k  ,iBLK), &
                   NodeNumberGlobal_NB(i  ,j+1,k  ,iBLK), &
                   NodeNumberGlobal_NB(i  ,j  ,k+1,iBLK), &
                   NodeNumberGlobal_NB(i+1,j  ,k+1,iBLK), &
                   NodeNumberGlobal_NB(i+1,j+1,k+1,iBLK), &
                   NodeNumberGlobal_NB(i  ,j+1,k+1,iBLK)
           end do; end do; end do
        end do
     end if
  case('cut','x=0','y=0','z=0')
     !================================ cut ============================
     ! Allocate memory for storing the blocks that are cut
     allocate( BlockCut(nBlockALL), stat=iError ); call alloc_check(iError,"BlockCut")
     BlockCut=0
     nBlockCuts=0

     if((xmax-xmin)<(ymax-ymin) .and. (xmax-xmin)<(zmax-zmin))then
        !X Slice
        CutValue = 0.5*(xmin+xmax)
        if(plot_type1(1:3) == 'x=0') CutValue = 0.
        if(IsCartesianGrid)then             
           ! First loop to count nodes and cells        
           do iBlockAll = 1, nNodeUsed
              iNode = iNodeMorton_I(iBlockAll)
              iBlk  = iTree_IA(Block_,iNode)
              iPE   = iTree_IA(Proc_,iNode)
              if(iProc==iPE)then
                 if ( CutValue> PlotXYZNodes_DNB(1,1   ,1,1,iBLK) .and. &
                      CutValue<=PlotXYZNodes_DNB(1,1+nI,1,1,iBLK)  )then
                    nBlockCuts=nBlockCuts+1
                    BlockCut(iBlockALL)=nBlockCuts
                 end if
              end if
              call MPI_Bcast(nBlockCuts,1,MPI_Integer,iPE,iComm,iError)
           end do
           if(iProc==0)then
              ! Write file header
              write(UnitTmp_,'(a)')'TITLE="BATSRUS: Cut X Data, '//textDateTime//'"'
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
              iBlk  = iTree_IA(Block_,iNode)
              iPE   = iTree_IA(Proc_,iNode)
              if(iProc==iPE)then
                 if ( CutValue> PlotXYZNodes_DNB(1,1   ,1,1,iBLK) .and. &
                      CutValue<=PlotXYZNodes_DNB(1,1+nI,1,1,iBLK) )then
                    ! Find cut interpolation factors
                    do i=1,nI
                       if ( CutValue> PlotXYZNodes_DNB(1,i  ,1,1,iBLK) .and. &
                            CutValue<=PlotXYZNodes_DNB(1,i+1,1,1,iBLK)  )then
                          cut1=i
                          cut2=i+1
                          factor2=(CutValue-PlotXYZNodes_DNB(1,i,1,1,iBLK))/ &
                               (PlotXYZNodes_DNB(1,i+1,1,1,iBLK)-PlotXYZNodes_DNB(1,i,1,1,iBLK))
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
                            (factor1*PlotVarNodes_VNB(1:nPlotVar,cut1,j,k,iBLK)+ &
                            factor2*PlotVarNodes_VNB(1:nPlotVar,cut2,j,k,iBLK))
                    end do; end do
                    ! Write point connectivity
                    do k=1,nK; do j=1,nJ
                       write(UnitTmp2_,'(4(i8,1x))') &
                            ((BlockCut(iBlockALL)-1)*(nJ+1)*(nK+1)) + (k-1)*(nJ+1)+j, &
                            ((BlockCut(iBlockALL)-1)*(nJ+1)*(nK+1)) + (k-1)*(nJ+1)+j+1, &
                            ((BlockCut(iBlockALL)-1)*(nJ+1)*(nK+1)) + (k  )*(nJ+1)+j+1, &
                            ((BlockCut(iBlockALL)-1)*(nJ+1)*(nK+1)) + (k  )*(nJ+1)+j
                    end do; end do
                 end if
              end if
           end do
        elseif(IsRLonLat)then
           ! First loop to count nodes and cells        
           do iBlockAll = 1, nNodeUsed
              iNode = iNodeMorton_I(iBlockAll)
              iBlk  = iTree_IA(Block_,iNode)
              iPE   = iTree_IA(Proc_,iNode)
              if(iProc==iPE)then
                 if ( (CutValue - PlotXYZNodes_DNB(1,1,1   ,2,iBLK))* &
                      (CutValue - PlotXYZNodes_DNB(1,1,1+nJ,2,iBLK)) <= 0.0  )then
                    nBlockCuts=nBlockCuts+1
                    BlockCut(iBlockALL)=nBlockCuts
                 end if
              end if
              call MPI_Bcast(nBlockCuts,1,MPI_Integer,iPE,iComm,iError)
           end do
           if(iProc==0)then
              ! Write file header
              write(UnitTmp_,'(a)')'TITLE="BATSRUS: Cut X Data, '//textDateTime//'"'
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
              iBlk  = iTree_IA(Block_,iNode)
              iPE   = iTree_IA(Proc_,iNode)
              if(iProc==iPE)then
                 if ( (CutValue - PlotXYZNodes_DNB(1,1,1   ,2,iBLK))* &
                      (CutValue - PlotXYZNodes_DNB(1,1,1+nJ,2,iBLK)) <= 0.0  )then
                    ! Find cut interpolation factors
                    do j=1, nJ
                       if ( (CutValue - PlotXYZNodes_DNB(1,1,j  ,2,iBLK))* &
                            (CutValue - PlotXYZNodes_DNB(1,1,j+1,2,iBLK)) <= 0.0  )then
                          cut1=j
                          cut2=j+1
                          factor2=(CutValue-PlotXYZNodes_DNB(1,1,j,2,iBLK))/ &
                               (PlotXYZNodes_DNB(1,1,j+1,2,iBLK)-PlotXYZNodes_DNB(1,1,j,2,iBLK))
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
                            (factor1*PlotVarNodes_VNB(1:nPlotVar,i,cut1,k,iBLK)+ &
                            factor2*PlotVarNodes_VNB(1:nPlotVar,i,cut2,k,iBLK))
                    end do; end do
                    ! Write point connectivity
                    do k=1,nK; do i=1,nI
                       write(UnitTmp2_,'(4(i8,1x))') &
                            ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) + (k-1)*(nI+1)+i, &
                            ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) + (k-1)*(nI+1)+i+1, &
                            ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) + (k  )*(nI+1)+i+1, &
                            ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) + (k  )*(nI+1)+i
                    end do; end do
                 end if
              end if
           end do
        end if

     elseif((ymax-ymin)<(zmax-zmin))then
        !Y Slice
        CutValue = 0.5*(ymin+ymax)
        if(plot_type1(1:3) == 'y=0') CutValue = 0.
        if(IsCartesianGrid)then                   
           ! First loop to count nodes and cells
           do iBlockAll = 1, nNodeUsed
              iNode = iNodeMorton_I(iBlockAll)
              iBlk  = iTree_IA(Block_,iNode)
              iPE   = iTree_IA(Proc_,iNode)
              if(iProc==iPE)then
                 if ( CutValue> PlotXYZNodes_DNB(2,1,1   ,1,iBLK) .and. &
                      CutValue<=PlotXYZNodes_DNB(2,1,1+nJ,1,iBLK)  )then
                    nBlockCuts=nBlockCuts+1
                    BlockCut(iBlockALL)=nBlockCuts
                 end if
              end if
              call MPI_Bcast(nBlockCuts,1,MPI_Integer,iPE,iComm,iError)
           end do
           if(iProc==0)then
              ! Write file header
              write(UnitTmp_,'(a)')'TITLE="BATSRUS: Cut Y Data, '//textDateTime//'"'
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
              iBlk  = iTree_IA(Block_,iNode)
              iPE   = iTree_IA(Proc_,iNode)
              if(iProc==iPE)then
                 if ( CutValue> PlotXYZNodes_DNB(2,1,1   ,1,iBLK) .and. &
                      CutValue<=PlotXYZNodes_DNB(2,1,1+nJ,1,iBLK)  )then
                    ! Find cut interpolation factors
                    do j=1,nJ
                       if ( CutValue> PlotXYZNodes_DNB(2,1,j  ,1,iBLK) .and. &
                            CutValue<=PlotXYZNodes_DNB(2,1,j+1,1,iBLK)  )then
                          cut1=j
                          cut2=j+1
                          factor2=(CutValue-PlotXYZNodes_DNB(2,1,j,1,iBLK))/ &
                               (PlotXYZNodes_DNB(2,1,j+1,1,iBLK)-PlotXYZNodes_DNB(2,1,j,1,iBLK))
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
                            (factor1*PlotVarNodes_VNB(1:nPlotVar,i,cut1,k,iBLK)+ &
                            factor2*PlotVarNodes_VNB(1:nPlotVar,i,cut2,k,iBLK))
                    end do; end do
                    ! Write point connectivity
                    do k=1,nK; do i=1,nI
                       write(UnitTmp2_,'(4(i8,1x))') &
                            ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) + (k-1)*(nI+1)+i, &
                            ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) + (k-1)*(nI+1)+i+1, &
                            ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) + (k  )*(nI+1)+i+1, &
                            ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) + (k  )*(nI+1)+i
                    end do; end do
                 end if
              end if
           end do
        else if(IsRLonLat) then 
           ! First loop to count nodes and cells
           do iBlockAll = 1, nNodeUsed
              iNode = iNodeMorton_I(iBlockAll)
              iBlk  = iTree_IA(Block_,iNode)
              iPE   = iTree_IA(Proc_,iNode)
              if(iProc==iPE)then
                 if ( (CutValue - PlotXYZNodes_DNB(2,1,1   ,2,iBLK))* &
                      (CutValue - PlotXYZNodes_DNB(2,1,1+nJ,2,iBLK)) <= 0.0  )then
                    nBlockCuts=nBlockCuts+1
                    BlockCut(iBlockALL)=nBlockCuts
                 end if
              end if
              call MPI_Bcast(nBlockCuts,1,MPI_Integer,iPE,iComm,iError)
           end do
           if(iProc==0)then
              ! Write file header
              write(UnitTmp_,'(a)')'TITLE="BATSRUS: Cut Y Data, '//textDateTime//'"'
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
              iBlk  = iTree_IA(Block_,iNode)
              iPE   = iTree_IA(Proc_,iNode)
              if(iProc==iPE)then
                 if ( (CutValue - PlotXYZNodes_DNB(2,1,1   ,2,iBLK))* &
                      (CutValue - PlotXYZNodes_DNB(2,1,1+nJ,2,iBLK)) <= 0.0  )then
                    ! Find cut interpolation factors
                    do j=1,nJ
                       if ( (CutValue - PlotXYZNodes_DNB(2,1,j  ,2,iBLK))* &
                            (CutValue - PlotXYZNodes_DNB(2,1,j+1,2,iBLK)) <= 0.0  )then
                          cut1=j
                          cut2=j+1
                          factor2=(CutValue-PlotXYZNodes_DNB(2,1,j,2,iBLK))/ &
                               (PlotXYZNodes_DNB(2,1,j+1,2,iBLK)-PlotXYZNodes_DNB(2,1,j,2,iBLK))
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
                            (factor1*PlotVarNodes_VNB(1:nPlotVar,i,cut1,k,iBLK)+ &
                            factor2*PlotVarNodes_VNB(1:nPlotVar,i,cut2,k,iBLK))
                    end do; end do
                    ! Write point connectivity
                    do k=1,nK; do i=1,nI
                       write(UnitTmp2_,'(4(i8,1x))') &
                            ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) + (k-1)*(nI+1)+i, &
                            ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) + (k-1)*(nI+1)+i+1, &
                            ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) + (k  )*(nI+1)+i+1, &
                            ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) + (k  )*(nI+1)+i
                    end do; end do
                 end if
              end if
           end do
        end if

     else
        !Z Slice
        CutValue = 0.5*(zmin+zmax)
        if(plot_type1(1:3) == 'z=0') CutValue = 0.
        ! First loop to count nodes and cells
        do iBlockAll = 1, nNodeUsed
           iNode = iNodeMorton_I(iBlockAll)
           iBlk  = iTree_IA(Block_,iNode)
           iPE   = iTree_IA(Proc_,iNode)
           if(iProc==iPE)then
              if ( CutValue> PlotXYZNodes_DNB(3,1,1,1   ,iBLK) .and. &
                   CutValue<=PlotXYZNodes_DNB(3,1,1,1+nK,iBLK)  )then
                 nBlockCuts=nBlockCuts+1
                 BlockCut(iBlockALL)=nBlockCuts
              end if
           end if
           call MPI_Bcast(nBlockCuts,1,MPI_Integer,iPE,iComm,iError)
        end do
        if(iProc==0)then
           ! Write file header
           write(UnitTmp_,'(a)')'TITLE="BATSRUS: Cut Z Data, '//textDateTime//'"'
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
           iBlk  = iTree_IA(Block_,iNode)
           iPE   = iTree_IA(Proc_,iNode)
           if(iProc==iPE)then
              if ( CutValue> PlotXYZNodes_DNB(3,1,1,1   ,iBLK) .and. &
                   CutValue<=PlotXYZNodes_DNB(3,1,1,1+nK,iBLK)  )then
                 ! Find cut interpolation factors
                 do k=1,nK
                    if ( CutValue> PlotXYZNodes_DNB(3,1,1,k  ,iBLK) .and. &
                         CutValue<=PlotXYZNodes_DNB(3,1,1,k+1,iBLK)  )then
                       cut1=k
                       cut2=k+1
                       factor2=(CutValue-PlotXYZNodes_DNB(3,1,1,k,iBLK))/ &
                            (PlotXYZNodes_DNB(3,1,1,k+1,iBLK)-PlotXYZNodes_DNB(3,1,1,k,iBLK))
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
                         (factor1*PlotVarNodes_VNB(1:nPlotVar,i,j,cut1,iBLK)+ &
                         factor2*PlotVarNodes_VNB(1:nPlotVar,i,j,cut2,iBLK))
                 end do; end do
                 ! Write point connectivity
                 do j=1,nJ; do i=1,nI
                    write(UnitTmp2_,'(4(i8,1x))') &
                         ((BlockCut(iBlockALL)-1)*(nI+1)*(nJ+1)) + (j-1)*(nI+1)+i, &
                         ((BlockCut(iBlockALL)-1)*(nI+1)*(nJ+1)) + (j-1)*(nI+1)+i+1, &
                         ((BlockCut(iBlockALL)-1)*(nI+1)*(nJ+1)) + (j  )*(nI+1)+i+1, &
                         ((BlockCut(iBlockALL)-1)*(nI+1)*(nJ+1)) + (j  )*(nI+1)+i
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
        !Point-Normal cut plot
        XarbP=plot_point(1,ifile); XarbNormal=plot_normal(1,ifile)
        YarbP=plot_point(2,ifile); YarbNormal=plot_normal(2,ifile)
        ZarbP=plot_point(3,ifile); ZarbNormal=plot_normal(3,ifile)
     else
        !Dipole cut plot
        XarbP=0.; XarbNormal=-sin(ThetaTilt)
        YarbP=0.; YarbNormal=0.
        ZarbP=0.; ZarbNormal= cos(ThetaTilt)
     end if

     ! First loop to count cuts
     nBlockCuts=0
     do iBLK = 1, nBlock
        if(iTypeAdvance_B(iBlk) == SkippedBlock_) CYCLE
        ic1=1; ic2=1+nI 
        jc1=1; jc2=1+nJ
        kc1=1; kc2=1+nK
        call find_cuts(-1)
        if ( nCuts>0 )then
           !count up number of cuts
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
     do iBLK = 1, nBlock
        if(iTypeAdvance_B(iBlk) == SkippedBlock_) CYCLE
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
     write(*,*)'Error in write_plot_tec: Unknown plot_type='//plot_type1
  end select

contains

  ! Assumes iBLK value is correct
  subroutine fill_nodeXYZ
    ! Fill array with position (optionally dimensioned)
    if (plot_dimensional(ifile)) then
       NodeXYZ_DN(1:3,:,:,:)=PlotXYZNodes_DNB(1:3,:,:,:,iBLK)*No2Io_V(UnitX_)
    else
       NodeXYZ_DN(1:3,:,:,:)=PlotXYZNodes_DNB(1:3,:,:,:,iBLK)
    end if
  end subroutine fill_nodeXYZ

  ! iopt =-1 check all edges to see if cut
  !      = 0 count cuts only
  !      = 1 find cuts and write to disk
  subroutine find_cuts(iopt)
    integer, intent(in) :: iopt
    integer :: ic,jc,kc

    nCuts=0

    !Check edges.
    ! X edges
    if (XarbNormal>0.01) then
       ic=ic1; jc=jc1; kc=kc1
       do jc=jc1,jc2,jc2-jc1; do kc=kc1,kc2,kc2-kc1
          if(iopt>-1 .and. (jc==jc1 .or. kc==kc1) .and. (jc/=0 .and. kc/=0)) CYCLE
          Yp=PlotXYZNodes_DNB(2,ic,jc,kc,iBLK)
          Zp=PlotXYZNodes_DNB(3,ic,jc,kc,iBLK)
          Xp=XarbP-( YarbNormal*(Yp-YarbP) + ZarbNormal*(Zp-ZarbP) )/XarbNormal
          if ( Xp> PlotXYZNodes_DNB(1,ic1,jc,kc,iBLK) .and. &
               Xp<=PlotXYZNodes_DNB(1,ic2,jc,kc,iBLK) )then
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
                factor2=(Xp-PlotXYZNodes_DNB(1,ic1,jc,kc,iBLK))/ &
                     (PlotXYZNodes_DNB(1,ic2,jc,kc,iBLK)-PlotXYZNodes_DNB(1,ic1,jc,kc,iBLK))
                factor1=1.-factor2
                write(UnitTmp_,fmt="(50(ES14.6))") &
                     (factor1*NodeXYZ_DN(:, ic1,jc,kc)+ &
                     factor2*NodeXYZ_DN(:, ic2,jc,kc)), &
                     (factor1*PlotVarNodes_VNB(1:nPlotVar,ic1,jc,kc,iBLK)+ &
                     factor2*PlotVarNodes_VNB(1:nPlotVar,ic2,jc,kc,iBLK))
                if(okdebug)write(*,*)'  i=',ic1,'-',ic2,' j=',jc,' k=',kc
             end if
          end if
       end do; end do
    end if
    ! Y edges
    if (YarbNormal>0.01) then
       ic=ic1; jc=jc1; kc=kc1
       do ic=ic1,ic2,ic2-ic1; do kc=kc1,kc2,kc2-kc1
          if(iopt>-1 .and. (ic==ic1 .or. kc==kc1) .and. (ic/=0 .and. kc/=0)) CYCLE
          Xp=PlotXYZNodes_DNB(1,ic,jc,kc,iBLK)
          Zp=PlotXYZNodes_DNB(3,ic,jc,kc,iBLK)
          Yp=YarbP-( XarbNormal*(Xp-XarbP) + ZarbNormal*(Zp-ZarbP) )/YarbNormal
          if ( Yp> PlotXYZNodes_DNB(2,ic,jc1,kc,iBLK) .and. &
               Yp<=PlotXYZNodes_DNB(2,ic,jc2,kc,iBLK) )then
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
                factor2=(Yp-PlotXYZNodes_DNB(2,ic,jc1,kc,iBLK))/ &
                     (PlotXYZNodes_DNB(2,ic,jc2,kc,iBLK)-PlotXYZNodes_DNB(2,ic,jc1,kc,iBLK))
                factor1=1.-factor2
                write(UnitTmp_,fmt="(50(ES14.6))") &
                     (factor1*NodeXYZ_DN(:, ic,jc1,kc)+ &
                     factor2*NodeXYZ_DN(:, ic,jc2,kc)), &
                     (factor1*PlotVarNodes_VNB(1:nPlotVar,ic,jc1,kc,iBLK)+ &
                     factor2*PlotVarNodes_VNB(1:nPlotVar,ic,jc2,kc,iBLK))
                if(okdebug)write(*,*)'  i=',ic,' j=',jc1,'-',jc2,' k=',kc
             end if
          end if
       end do; end do
    end if
    ! Z edges
    if (ZarbNormal>0.01) then
       ic=ic1; jc=jc1; kc=kc1
       do ic=ic1,ic2,ic2-ic1; do jc=jc1,jc2,jc2-jc1
          if(iopt>-1 .and. (ic==ic1 .or. jc==jc1) .and. (ic/=0 .and. jc/=0)) CYCLE
          Xp=PlotXYZNodes_DNB(1,ic,jc,kc,iBLK)
          Yp=PlotXYZNodes_DNB(2,ic,jc,kc,iBLK)
          Zp=ZarbP-( XarbNormal*(Xp-XarbP) + YarbNormal*(Yp-YarbP) )/ZarbNormal
          if ( Zp> PlotXYZNodes_DNB(3,ic,jc,kc1,iBLK) .and. &
               Zp<=PlotXYZNodes_DNB(3,ic,jc,kc2,iBLK) )then
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
                factor2=(Zp-PlotXYZNodes_DNB(3,ic,jc,kc1,iBLK))/ &
                     (PlotXYZNodes_DNB(3,ic,jc,kc2,iBLK)-PlotXYZNodes_DNB(3,ic,jc,kc1,iBLK))
                factor1=1.-factor2
                write(UnitTmp_,fmt="(50(ES14.6))") &
                     (factor1*NodeXYZ_DN(:, ic,jc,kc1)+ &
                     factor2*NodeXYZ_DN(:, ic,jc,kc2)), &
                     (factor1*PlotVarNodes_VNB(1:nPlotVar,ic,jc,kc1,iBLK)+ &
                     factor2*PlotVarNodes_VNB(1:nPlotVar,ic,jc,kc2,iBLK))
                if(okdebug)write(*,*)'  i=',ic,' j=',jc,' k=',kc1,'-',kc2
             end if
          end if
       end do; end do
    end if

  end subroutine find_cuts

end subroutine write_plot_tec

!==========================================================================
subroutine assign_node_numbers
  use ModProcMH
  use ModIO, ONLY: write_prefix, iUnitOut
  use ModMain, ONLY : lVerbose, nBlock, nBlockMax, nBlockALL
  use ModAdvance,  ONLY: iTypeAdvance_B, iTypeAdvance_BP, SkippedBlock_
  use ModNodes
  use ModMpi
  use BATL_lib, ONLY: message_pass_node
  implicit none

  integer, parameter :: NodesPerBlock=(nI+1)*(nJ+1)*(nK+1)
  integer :: iBlockStart
  integer :: i, j, k, iNode, iBLK, iError, iPE, iTag
  integer :: nOffset, nOffsetPrevious
  integer, allocatable, dimension(:) :: NodeOffset, NodeOffsetMax, nOffset_P
  real, allocatable, dimension(:,:,:,:,:) :: IndexNode_VNB
  logical :: DoAllReduce=.true.
  integer :: iStatus(MPI_STATUS_SIZE)

  !-------------------------------------------------------------------------

  ! Write information to the screen
  if(iProc==0.and.lVerbose>0)then
     call write_prefix; write(iUnitOut,*)'Starting assign_node_numbers ...'
  end if

  ! Initialize all node numbers to zero
  NodeNumberLocal_NB=0

  ! Number of nodes on each block (maximum)
  nNodeALL=nBlockALL*NodesPerBlock

  ! Count number of used blocks on all processors with rank less than this one
  iBlockStart = 0
  if(iProc > 0) iBlockStart = &
       count(iTypeAdvance_BP(1:nBlockMax,0:iProc-1) /= SkippedBlock_)

  iNode = iBlockStart*NodesPerBlock

  ! Loop to assign local and global node numbers
  TREE1: do iBlk  = 1, nBlock
     if(iTypeAdvance_B(iBlk) == SkippedBlock_) CYCLE
     do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
        iNode = iNode+1
        NodeNumberLocal_NB(i,j,k,iBlk)= iNode
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

  !Allocate memory for storing the node offsets
  allocate( NodeOffset   (nBlockALL*NodesPerBlock), stat=iError)
  call alloc_check(iError,"NodeOffset")
  allocate( NodeOffsetMax(nBlockALL*NodesPerBlock), stat=iError)
  call alloc_check(iError,"NodeOffsetMax")
  NodeOffset=0

  ! Loop to compute node offsets
  nOffset=0
  TREE2: do iBLK  = 1, nBlock
     if(iTypeAdvance_B(iBLK) == SkippedBlock_) CYCLE
     do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
        if(NodeNumberLocal_NB(i,j,k,iBLK) > NodeNumberGlobal_NB(i,j,k,iBLK))then
           nOffset = nOffset+1
           NodeUniqueGlobal_NB(i,j,k,iBLK) = .false.
        end if
        NodeOffset(NodeNumberLocal_NB(i,j,k,iBLK)) = nOffset
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
  do iBLK  = 1, nBlock
     if(iTypeAdvance_B(iBLK) == SkippedBlock_) CYCLE
     do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
        iNode = NodeNumberLocal_NB(i,j,k,iBLK)
        NodeOffset(iNode) = NodeOffset(iNode) + nOffsetPrevious
     end do; end do; end do
  end do

  ! Gather offsets from all PE-s. NodeOffset was initialized to 0 so MPI_MAX works.
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
     call MPI_Bcast(NodeOffset,nBlockALL*NodesPerBlock,MPI_Integer,0,iComm,iError)
  end if

  ! Loop to fix NodeNumberGlobal_NB for offset
  TREE3: do iBlk  = 1, nBlock
     if(iTypeAdvance_B(iBLK) == SkippedBlock_) CYCLE
     do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
        NodeNumberGlobal_NB(i,j,k,iBLK) = NodeNumberGlobal_NB(i,j,k,iBLK) &
             - NodeOffset(NodeNumberGlobal_NB(i,j,k,iBLK))
        if(NodeNumberGlobal_NB(i,j,k,iBLK)>nNodeALL &
             .or. NodeNumberGlobal_NB(i,j,k,iBLK)<1)then
           ! Error in numbering, report values and stop.
           write(*,*)'ERROR: Global node numbering problem.', &
                ' PE=',iProc,' BLK=',iBLK,' ijk=',i,j,k
           write(*,*)'  NodeNumberGlobal_NB=',&
                NodeNumberGlobal_NB(i,j,k,iBLK)
           write(*,*)'  NodeOffset           =',&
                NodeOffset(NodeNumberGlobal_NB(i,j,k,iBLK))
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

end subroutine assign_node_numbers
