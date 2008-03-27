!^CFG COPYRIGHT UM
subroutine write_plot_tec(ifile,nPlotVar,PlotVarBlk,PlotVarNodes_NBI,unitstr_TEC,&
     xmin,xmax,ymin,ymax,zmin,zmax)
  !
  !NOTE: This routine assumes that the blocks are sorted on PEs by their global
  !       block number, ie blocks 1 to n on PE 0, blocks n+1 to n+m on PE 1,
  !       etc.
  !
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,globalBLK,global_block_number, nBlock, &
       nBlockALL,nBlockMax, time_accurate,n_step,&
       nOrder, UseRotatingBc, BlkTest, ProcTest, &
       TypeCoordSystem, CodeVersion
  use ModFaceValue, ONLY: TypeLimiter, BetaLimiter
  use ModMain, ONLY: boris_correction                     !^CFG IF BORISCORR
  use ModCovariant, ONLY: UseCovariant, TypeGeometry      
  use ModParallel, ONLY : iBlock_A, iProc_A
  use ModPhysics, ONLY : No2Io_V, UnitX_, &
       ThetaTilt, Rbody, boris_cLIGHT_factor, BodyNDim_I, g
  use ModAdvance, ONLY : FluxType, iTypeAdvance_B, SkippedBlock_
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK
  use ModIO
  use ModNodes
  use ModNumConst, ONLY : cRadToDeg
  use ModMpi
  implicit none

  ! Arguments  
  integer, intent(in) :: ifile, nPlotVar
  character (LEN=500), intent(in) :: unitstr_TEC
  real, intent(in) :: PlotVarBLK(-1:nI+2,-1:nJ+2,-1:nK+2,nPlotVarMax)
  real, intent(in) :: PlotVarNodes_NBI(1:1+nI,1:1+nJ,1:1+nK,nBLK,nPlotVarMax)
  real, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax

  ! Local Variables
  integer :: i,j,k, cut1,cut2, iPE,iBLK, iBlockALL, nBlockCuts, iError
  real :: CutValue, factor1,factor2
  logical :: oktest,oktest_me
  integer, allocatable, dimension(:) :: BlockCut
  character (len=22) :: textNandT
  character (len=23) :: textDateTime0,textDateTime
  character (len=80) :: format
  character(len=80) :: stmp

  integer :: iTime0_I(7),iTime_I(7)

  integer ic,ic1,ic2, jc,jc1,jc2, kc,kc1,kc2, nCuts, nCutsTotal
  real :: XarbP,YarbP,ZarbP, XarbNormal,YarbNormal,ZarbNormal, Xp,Yp,Zp
  real, dimension(1:nI+1,1:nJ+1,1:nK+1,3) :: NodeXYZ_N
  logical :: okdebug
  !----------------------------------------------------------------------------

  call set_oktest('write_plot_tec',oktest,oktest_me)
  if(oktest_me)write(*,*) plot_type1,plot_type1(1:3)  

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

  write(format,*)'(i4.4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2,":",i2.2,".",i3.3)'
  call get_date_time_start(iTime0_I)
  call get_date_time(iTime_I)
  write(textDateTime0,format) iTime0_I
  write(textDateTime ,format) iTime_I

  select case(plot_type1(1:3))
  case('blk')
     do iBLK = 1, nBlock
        if(iTypeAdvance_B(iBlk) == SkippedBlock_) CYCLE
        if ( plot_point(1,ifile)> NodeX_NB(1   ,1   ,1   ,iBLK) .and. &
             plot_point(1,ifile)<=NodeX_NB(1+nI,1+nJ,1+nK,iBLK) .and. &
             plot_point(2,ifile)> NodeY_NB(1   ,1   ,1   ,iBLK) .and. &
             plot_point(2,ifile)<=NodeY_NB(1+nI,1+nJ,1+nK,iBLK) .and. &
             plot_point(3,ifile)> NodeZ_NB(1   ,1   ,1   ,iBLK) .and. &
             plot_point(3,ifile)<=NodeZ_NB(1+nI,1+nJ,1+nK,iBLK) )then
           write(unit_tmp,'(a)')'TITLE="BATSRUS: BLK Only, '//textDateTime//'"'
           write(unit_tmp,'(a)')trim(unitstr_TEC)
           write(unit_tmp,'(a,i8,a,i8,a,i8,a)') &
                'ZONE T="BLK Only '//textNandT//'", I=',nI+4,&
                ', J=',nJ+4,', K=',nK+4,', F=POINT'
           call write_auxdata
           !DEBUGBLK
           write(stmp,'(i12)')iBLK
           write(unit_tmp,'(a,a,a)') 'AUXDATA DEBUGBLK="',trim(adjustl(stmp)),'"'
           !DEBUGPROC
           write(stmp,'(i12)')iProc
           write(unit_tmp,'(a,a,a)') 'AUXDATA DEBUGPROC="',trim(adjustl(stmp)),'"'
           ! Write cell values
           do k=-1,nK+2; do j=-1,nJ+2; do i=-1,nI+2
              if (plot_dimensional(ifile)) then
                 write(unit_tmp,fmt="(30(E14.6))") &
                      x_BLK(i,j,k,iBLK)*No2Io_V(UnitX_), &
                      y_BLK(i,j,k,iBLK)*No2Io_V(UnitX_), &
                      z_BLK(i,j,k,iBLK)*No2Io_V(UnitX_), &
                      PlotVarBlk(i,j,k,1:nPlotVar)
              else
                 write(unit_tmp,fmt="(30(E14.6))") &
                      x_BLK(i,j,k,iBLK), &
                      y_BLK(i,j,k,iBLK), &
                      z_BLK(i,j,k,iBLK), &
                      PlotVarBlk(i,j,k,1:nPlotVar)
              end if
           end do; end do; end do
        end if
     end do
  case('1d_')
     if(iProc==0)then
        ! Write file header
        write(unit_tmp,'(a)')'TITLE="BATSRUS: 1D Block Data, '//textDateTime//'"'
        write(unit_tmp,'(a)')trim(unitstr_TEC)
        write(unit_tmp,'(a,a,i8,a,a)') &
             'ZONE T="1D   '//textNandT//'"', &
             ', I=',nBlockALL,', J=1, K=1,', &
             ', ZONETYPE=ORDERED, DATAPACKING=POINT'
        call write_auxdata
     end if
     !================================= 1d ============================
     do iBlockALL  = 1, nBlockALL
        iBLK = iBlock_A(iBlockALL)
        iPE  = iProc_A(iBlockALL)
        if(iProc==iPE)then
           ! Write point values
           call fill_NodeXYZ
           write(unit_tmp,fmt="(30(E14.6))") &
                (NodeXYZ_N(1,1,1,1)+NodeXYZ_N(nI+1,1,1,1))/2., &
                (NodeXYZ_N(1,1,1,2)+NodeXYZ_N(1,nJ+1,1,2))/2., &
                (NodeXYZ_N(1,1,1,3)+NodeXYZ_N(1,1,nK+1,3))/2., &
                PlotVarNodes_NBI(2,2,2,iBLK,1:nPlotVar)
        end if
     end do
  case('3d_')
     if(iProc==0)then
        ! Write file header
        write(unit_tmp,'(a)')'TITLE="BATSRUS: 3D Data, '//textDateTime//'"'
        write(unit_tmp,'(a)')trim(unitstr_TEC)
        write(unit_tmp,'(a,a,i8,a,i8,a)') &
             'ZONE T="3D   '//textNandT//'"', &
             ', N=',nNodeALL, &
             ', E=',nBlockALL*((nI  )*(nJ  )*(nK  )), &
             ', F=FEPOINT, ET=BRICK'
        call write_auxdata
     end if
     !================================= 3d ============================
     do iBLK = 1, nBlock
        if(iTypeAdvance_B(iBlk) == SkippedBlock_) CYCLE
        ! Write point values
        call fill_NodeXYZ
        do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
           if(NodeUniqueGlobal_NB(i,j,k,iBLK))then
              write(unit_tmp,fmt="(30(E14.6))") &
                   NodeXYZ_N(i,j,k,1:3),PlotVarNodes_NBI(i,j,k,iBLK,1:nPlotVar)
           end if
        end do; end do; end do
        ! Write point connectivity
        do k=1,nK; do j=1,nJ; do i=1,nI
           write(unit_tmp2,'(8(i8,1x))') &
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
  case('cut','x=0','y=0','z=0')
     !================================ cut ============================
     ! Allocate memory for storing the blocks that are cut
     allocate( BlockCut(nBlockALL), stat=iError ); call alloc_check(iError,"BlockCut")
     BlockCut=0
     nBlockCuts=0

     if((xmax-xmin)<(ymax-ymin) .and. (xmax-xmin)<(zmax-zmin))then
        !X Slice
        CutValue = 0.5*(xmin+xmax)
        if(.not.UseCovariant)then             
           ! First loop to count nodes and cells        
           do iBlockALL  = 1, nBlockALL
              iBLK = iBlock_A(iBlockALL)
              iPE  = iProc_A(iBlockALL)
              if(iProc==iPE)then
                 if ( CutValue> NodeX_NB(1   ,1,1,iBLK) .and. &
                      CutValue<=NodeX_NB(1+nI,1,1,iBLK)  )then
                    nBlockCuts=nBlockCuts+1
                    BlockCut(iBlockALL)=nBlockCuts
                 end if
              end if
              call MPI_Bcast(nBlockCuts,1,MPI_Integer,iPE,iComm,iError)
           end do
           if(iProc==0)then
              ! Write file header
              write(unit_tmp,'(a)')'TITLE="BATSRUS: Cut X Data, '//textDateTime//'"'
              write(unit_tmp,'(a)')trim(unitstr_TEC)
              write(unit_tmp,'(a,a,i8,a,i8,a)') &
                   'ZONE T="2D X '//textNandT//'"', &
                   ', N=',nBlockCuts*((nJ+1)*(nK+1)), &
                   ', E=',nBlockCuts*((nJ  )*(nK  )), &
                   ', F=FEPOINT, ET=QUADRILATERAL'
              call write_auxdata
           end if
           ! Now loop to write values
           do iBlockALL  = 1, nBlockALL
              iBLK = iBlock_A(iBlockALL)
              iPE  = iProc_A(iBlockALL)
              if(iProc==iPE)then
                 if ( CutValue> NodeX_NB(1   ,1,1,iBLK) .and. &
                      CutValue<=NodeX_NB(1+nI,1,1,iBLK) )then
                    ! Find cut interpolation factors
                    do i=1,nI
                       if ( CutValue> NodeX_NB(i  ,1,1,iBLK) .and. &
                            CutValue<=NodeX_NB(i+1,1,1,iBLK)  )then
                          cut1=i
                          cut2=i+1
                          factor2=(CutValue-NodeX_NB(i,1,1,iBLK))/ &
                               (NodeX_NB(i+1,1,1,iBLK)-NodeX_NB(i,1,1,iBLK))
                          factor1=1.-factor2
                          EXIT
                       end if
                    end do
                    ! Write point values
                    call fill_NodeXYZ
                    do k=1,1+nK; do j=1,1+nJ
                       write(unit_tmp,fmt="(30(E14.6))") &
                            (factor1*NodeXYZ_N(cut1,j,k,1:3)+ &
                            factor2*NodeXYZ_N(cut2,j,k,1:3)), &
                            (factor1*PlotVarNodes_NBI(cut1,j,k,iBLK,1:nPlotVar)+ &
                            factor2*PlotVarNodes_NBI(cut2,j,k,iBLK,1:nPlotVar))
                    end do; end do
                    ! Write point connectivity
                    do k=1,nK; do j=1,nJ
                       write(unit_tmp2,'(4(i8,1x))') &
                            ((BlockCut(iBlockALL)-1)*(nJ+1)*(nK+1)) + (k-1)*(nJ+1)+j, &
                            ((BlockCut(iBlockALL)-1)*(nJ+1)*(nK+1)) + (k-1)*(nJ+1)+j+1, &
                            ((BlockCut(iBlockALL)-1)*(nJ+1)*(nK+1)) + (k  )*(nJ+1)+j+1, &
                            ((BlockCut(iBlockALL)-1)*(nJ+1)*(nK+1)) + (k  )*(nJ+1)+j
                    end do; end do
                 end if
              end if
           end do
        else if(TypeGeometry == 'spherical_lnr' .or. & 
             TypeGeometry == 'spherical') then 
           ! First loop to count nodes and cells        
           do iBlockALL  = 1, nBlockALL
              iBLK = iBlock_A(iBlockALL)
              iPE  = iProc_A(iBlockALL)
              if(iProc==iPE)then
                 if ( (CutValue - NodeX_NB(1,1   ,2,iBLK))* &
                      (CutValue - NodeX_NB(1,1+nJ,2,iBLK)) <= 0.0  )then
                    nBlockCuts=nBlockCuts+1
                    BlockCut(iBlockALL)=nBlockCuts
                 end if
              end if
              call MPI_Bcast(nBlockCuts,1,MPI_Integer,iPE,iComm,iError)
           end do
           if(iProc==0)then
              ! Write file header
              write(unit_tmp,'(a)')'TITLE="BATSRUS: Cut X Data, '//textDateTime//'"'
              write(unit_tmp,'(a)')trim(unitstr_TEC)
              write(unit_tmp,'(a,a,i8,a,i8,a)')       &
                   'ZONE T="2D X '//textNandT//'"',   &
                   ', N=',nBlockCuts*((nI+1)*(nK+1)), &
                   ', E=',nBlockCuts*((nI  )*(nK  )), &
                   ', F=FEPOINT, ET=QUADRILATERAL'
              call write_auxdata
           end if
           ! Now loop to write values
           do iBlockALL  = 1, nBlockALL
              iBLK = iBlock_A(iBlockALL)
              iPE  = iProc_A(iBlockALL)
              if(iProc==iPE)then
                 if ( (CutValue - NodeX_NB(1,1   ,2,iBLK))* &
                      (CutValue - NodeX_NB(1,1+nJ,2,iBLK)) <= 0.0  )then
                    ! Find cut interpolation factors
                    do j=1, nJ
                       if ( (CutValue - NodeX_NB(1,j  ,2,iBLK))* &
                            (CutValue - NodeX_NB(1,j+1,2,iBLK)) <= 0.0  )then
                          cut1=j
                          cut2=j+1
                          factor2=(CutValue-NodeX_NB(1,j,2,iBLK))/ &
                               (NodeX_NB(1,j+1,2,iBLK)-NodeX_NB(1,j,2,iBLK))
                          factor1=1.-factor2
                          EXIT
                       end if
                    end do
                    ! Write point values
                    call fill_NodeXYZ                    
                    do k=1,1+nK; do i=1,1+nI
                       write(unit_tmp,fmt="(30(E14.6))") &
                            (factor1*NodeXYZ_N(i,cut1,k,1:3)+ &
                            factor2*NodeXYZ_N(i,cut2,k,1:3)), &
                            (factor1*PlotVarNodes_NBI(i,cut1,k,iBLK,1:nPlotVar)+ &
                            factor2*PlotVarNodes_NBI(i,cut2,k,iBLK,1:nPlotVar))
                    end do; end do
                    ! Write point connectivity
                    do k=1,nK; do i=1,nI
                       write(unit_tmp2,'(4(i8,1x))') &
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
        if(.not.UseCovariant)then                   
           ! First loop to count nodes and cells
           do iBlockALL  = 1, nBlockALL
              iBLK = iBlock_A(iBlockALL)
              iPE  = iProc_A(iBlockALL)
              if(iProc==iPE)then
                 if ( CutValue> NodeY_NB(1,1   ,1,iBLK) .and. &
                      CutValue<=NodeY_NB(1,1+nJ,1,iBLK)  )then
                    nBlockCuts=nBlockCuts+1
                    BlockCut(iBlockALL)=nBlockCuts
                 end if
              end if
              call MPI_Bcast(nBlockCuts,1,MPI_Integer,iPE,iComm,iError)
           end do
           if(iProc==0)then
              ! Write file header
              write(unit_tmp,'(a)')'TITLE="BATSRUS: Cut Y Data, '//textDateTime//'"'
              write(unit_tmp,'(a)')unitstr_TEC(1:len_trim(unitstr_TEC))
              write(unit_tmp,'(a,a,i8,a,i8,a)') &
                   'ZONE T="2D Y '//textNandT//'"', &
                   ', N=',nBlockCuts*((nI+1)*(nK+1)), &
                   ', E=',nBlockCuts*((nI  )*(nK  )), &
                   ', F=FEPOINT, ET=QUADRILATERAL'
              call write_auxdata
           end if
           ! Now loop to write values
           do iBlockALL  = 1, nBlockALL
              iBLK = iBlock_A(iBlockALL)
              iPE  = iProc_A(iBlockALL)
              if(iProc==iPE)then
                 if ( CutValue> NodeY_NB(1,1   ,1,iBLK) .and. &
                      CutValue<=NodeY_NB(1,1+nJ,1,iBLK)  )then
                    ! Find cut interpolation factors
                    do j=1,nJ
                       if ( CutValue> NodeY_NB(1,j  ,1,iBLK) .and. &
                            CutValue<=NodeY_NB(1,j+1,1,iBLK)  )then
                          cut1=j
                          cut2=j+1
                          factor2=(CutValue-NodeY_NB(1,j,1,iBLK))/ &
                               (NodeY_NB(1,j+1,1,iBLK)-NodeY_NB(1,j,1,iBLK))
                          factor1=1.-factor2
                          EXIT
                       end if
                    end do
                    ! Write point values
                    call fill_NodeXYZ
                    do k=1,1+nK; do i=1,1+nI
                       write(unit_tmp,fmt="(30(E14.6))") &
                            (factor1*NodeXYZ_N(i,cut1,k,1:3)+ &
                            factor2*NodeXYZ_N(i,cut2,k,1:3)), &
                            (factor1*PlotVarNodes_NBI(i,cut1,k,iBLK,1:nPlotVar)+ &
                            factor2*PlotVarNodes_NBI(i,cut2,k,iBLK,1:nPlotVar))
                    end do; end do
                    ! Write point connectivity
                    do k=1,nK; do i=1,nI
                       write(unit_tmp2,'(4(i8,1x))') &
                            ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) + (k-1)*(nI+1)+i, &
                            ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) + (k-1)*(nI+1)+i+1, &
                            ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) + (k  )*(nI+1)+i+1, &
                            ((BlockCut(iBlockALL)-1)*(nI+1)*(nK+1)) + (k  )*(nI+1)+i
                    end do; end do
                 end if
              end if
           end do
        else if(TypeGeometry == 'spherical_lnr' .or. & 
             TypeGeometry == 'spherical') then 
           ! First loop to count nodes and cells
           do iBlockALL  = 1, nBlockALL
              iBLK = iBlock_A(iBlockALL)
              iPE  = iProc_A(iBlockALL)
              if(iProc==iPE)then
                 if ( (CutValue - NodeY_NB(1,1   ,2,iBLK))* &
                      (CutValue - NodeY_NB(1,1+nJ,2,iBLK)) <= 0.0  )then
                    nBlockCuts=nBlockCuts+1
                    BlockCut(iBlockALL)=nBlockCuts
                 end if              
              end if
              call MPI_Bcast(nBlockCuts,1,MPI_Integer,iPE,iComm,iError)
           end do
           if(iProc==0)then
              ! Write file header
              write(unit_tmp,'(a)')'TITLE="BATSRUS: Cut Y Data, '//textDateTime//'"'
              write(unit_tmp,'(a)')unitstr_TEC(1:len_trim(unitstr_TEC))
              write(unit_tmp,'(a,a,i8,a,i8,a)') &
                   'ZONE T="2D Y '//textNandT//'"', &
                   ', N=',nBlockCuts*((nI+1)*(nK+1)), &
                   ', E=',nBlockCuts*((nI  )*(nK  )), &
                   ', F=FEPOINT, ET=QUADRILATERAL'
              call write_auxdata
           end if
           ! Now loop to write values
           do iBlockALL  = 1, nBlockALL
              iBLK = iBlock_A(iBlockALL)
              iPE  = iProc_A(iBlockALL)
              if(iProc==iPE)then
                 if ( (CutValue - NodeY_NB(1,1   ,2,iBLK))* &
                      (CutValue - NodeY_NB(1,1+nJ,2,iBLK)) <= 0.0  )then
                    ! Find cut interpolation factors
                    do j=1,nJ
                       if ( (CutValue - NodeY_NB(1,j  ,2,iBLK))* &
                            (CutValue - NodeY_NB(1,j+1,2,iBLK)) <= 0.0  )then
                          cut1=j
                          cut2=j+1
                          factor2=(CutValue-NodeY_NB(1,j,2,iBLK))/ &
                               (NodeY_NB(1,j+1,2,iBLK)-NodeY_NB(1,j,2,iBLK))
                          factor1=1.-factor2
                          EXIT
                       end if
                    end do
                    ! Write point values
                    call fill_NodeXYZ
                    do k=1,1+nK; do i=1,1+nI
                       write(unit_tmp,fmt="(30(E14.6))") &
                            (factor1*NodeXYZ_N(i,cut1,k,1:3)+ &
                            factor2*NodeXYZ_N(i,cut2,k,1:3)), &
                            (factor1*PlotVarNodes_NBI(i,cut1,k,iBLK,1:nPlotVar)+ &
                            factor2*PlotVarNodes_NBI(i,cut2,k,iBLK,1:nPlotVar))
                    end do; end do
                    ! Write point connectivity
                    do k=1,nK; do i=1,nI
                       write(unit_tmp2,'(4(i8,1x))') &
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
        ! First loop to count nodes and cells
        do iBlockALL  = 1, nBlockALL
           iBLK = iBlock_A(iBlockALL)
           iPE  = iProc_A(iBlockALL)
           if(iProc==iPE)then
              if ( CutValue> NodeZ_NB(1,1,1   ,iBLK) .and. &
                   CutValue<=NodeZ_NB(1,1,1+nK,iBLK)  )then
                 nBlockCuts=nBlockCuts+1
                 BlockCut(iBlockALL)=nBlockCuts
              end if
           end if
           call MPI_Bcast(nBlockCuts,1,MPI_Integer,iPE,iComm,iError)
        end do
        if(iProc==0)then
           ! Write file header
           write(unit_tmp,'(a)')'TITLE="BATSRUS: Cut Z Data, '//textDateTime//'"'
           write(unit_tmp,'(a)')unitstr_TEC(1:len_trim(unitstr_TEC))
           write(unit_tmp,'(a,a,i8,a,i8,a)') &
                'ZONE T="2D Z '//textNandT//'"', &
                ', N=',nBlockCuts*((nI+1)*(nJ+1)), &
                ', E=',nBlockCuts*((nI  )*(nJ  )), &
                ', F=FEPOINT, ET=QUADRILATERAL'
           call write_auxdata
        end if
        ! Now loop to write values
        do iBlockALL  = 1, nBlockALL
           iBLK = iBlock_A(iBlockALL)
           iPE  = iProc_A(iBlockALL)
           if(iProc==iPE)then
              if ( CutValue> NodeZ_NB(1,1,1   ,iBLK) .and. &
                   CutValue<=NodeZ_NB(1,1,1+nK,iBLK)  )then
                 ! Find cut interpolation factors
                 do k=1,nK
                    if ( CutValue> NodeZ_NB(1,1,k  ,iBLK) .and. &
                         CutValue<=NodeZ_NB(1,1,k+1,iBLK)  )then
                       cut1=k
                       cut2=k+1
                       factor2=(CutValue-NodeZ_NB(1,1,k,iBLK))/ &
                            (NodeZ_NB(1,1,k+1,iBLK)-NodeZ_NB(1,1,k,iBLK))
                       factor1=1.-factor2
                       EXIT
                    end if
                 end do
                 ! Write point values
                 call fill_NodeXYZ
                 do j=1,1+nJ; do i=1,1+nI
                    write(unit_tmp,fmt="(30(E14.6))") &
                         (factor1*NodeXYZ_N(i,j,cut1,1:3)+ &
                          factor2*NodeXYZ_N(i,j,cut2,1:3)), &
                         (factor1*PlotVarNodes_NBI(i,j,cut1,iBLK,1:nPlotVar)+ &
                          factor2*PlotVarNodes_NBI(i,j,cut2,iBLK,1:nPlotVar))
                 end do; end do
                 ! Write point connectivity
                 do j=1,nJ; do i=1,nI
                    write(unit_tmp2,'(4(i8,1x))') &
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
           write(unit_tmp,'(a)')'TITLE="BATSRUS: Slice, '//textDateTime//'"'
           write(unit_tmp,'(a)')trim(unitstr_TEC)
           write(unit_tmp,'(a,i8,a)') &
                'ZONE T="Slice '//textNandT//'", I=', nCutsTotal,&
                ', J=1, K=1, F=POINT'
        else
           write(unit_tmp,'(a)')'TITLE="BATSRUS: Dipole Cut, '// &
                textDateTime//'"'
           write(unit_tmp,'(a)')trim(unitstr_TEC)
           write(unit_tmp,'(a,i8,a)') &
                'ZONE T="Dipole Cut '//textNandT//'", I=', nCutsTotal,&
                ', J=1, K=1, F=POINT'
        end if
        call write_auxdata
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
       NodeXYZ_N(:,:,:,1)=NodeX_NB(:,:,:,iBLK)*No2Io_V(UnitX_)
       NodeXYZ_N(:,:,:,2)=NodeY_NB(:,:,:,iBLK)*No2Io_V(UnitX_)
       NodeXYZ_N(:,:,:,3)=NodeZ_NB(:,:,:,iBLK)*No2Io_V(UnitX_)
    else
       NodeXYZ_N(:,:,:,1)=NodeX_NB(:,:,:,iBLK)
       NodeXYZ_N(:,:,:,2)=NodeY_NB(:,:,:,iBLK)
       NodeXYZ_N(:,:,:,3)=NodeZ_NB(:,:,:,iBLK)
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
          Yp=NodeY_NB(ic,jc,kc,iBLK)
          Zp=NodeZ_NB(ic,jc,kc,iBLK)
          Xp=XarbP-( YarbNormal*(Yp-YarbP) + ZarbNormal*(Zp-ZarbP) )/XarbNormal
          if ( Xp> NodeX_NB(ic1,jc,kc,iBLK) .and. &
               Xp<=NodeX_NB(ic2,jc,kc,iBLK) )then
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
                factor2=(Xp-NodeX_NB(ic1,jc,kc,iBLK))/ &
                     (NodeX_NB(ic2,jc,kc,iBLK)-NodeX_NB(ic1,jc,kc,iBLK))
                factor1=1.-factor2
                write(unit_tmp,fmt="(30(E14.6))") &
                     (factor1*NodeXYZ_N( ic1,jc,kc,:)+ &
                      factor2*NodeXYZ_N( ic2,jc,kc,:)), &
                     (factor1*PlotVarNodes_NBI(ic1,jc,kc,iBLK,1:nPlotVar)+ &
                      factor2*PlotVarNodes_NBI(ic2,jc,kc,iBLK,1:nPlotVar))
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
          Xp=NodeX_NB(ic,jc,kc,iBLK)
          Zp=NodeZ_NB(ic,jc,kc,iBLK)
          Yp=YarbP-( XarbNormal*(Xp-XarbP) + ZarbNormal*(Zp-ZarbP) )/YarbNormal
          if ( Yp> NodeY_NB(ic,jc1,kc,iBLK) .and. &
               Yp<=NodeY_NB(ic,jc2,kc,iBLK) )then
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
                factor2=(Yp-NodeY_NB(ic,jc1,kc,iBLK))/ &
                     (NodeY_NB(ic,jc2,kc,iBLK)-NodeY_NB(ic,jc1,kc,iBLK))
                factor1=1.-factor2
                write(unit_tmp,fmt="(30(E14.6))") &
                     (factor1*NodeXYZ_N( ic,jc1,kc,:)+ &
                      factor2*NodeXYZ_N( ic,jc2,kc,:)), &
                     (factor1*PlotVarNodes_NBI(ic,jc1,kc,iBLK,1:nPlotVar)+ &
                      factor2*PlotVarNodes_NBI(ic,jc2,kc,iBLK,1:nPlotVar))
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
          Xp=NodeX_NB(ic,jc,kc,iBLK)
          Yp=NodeY_NB(ic,jc,kc,iBLK)
          Zp=ZarbP-( XarbNormal*(Xp-XarbP) + YarbNormal*(Yp-YarbP) )/ZarbNormal
          if ( Zp> NodeZ_NB(ic,jc,kc1,iBLK) .and. &
               Zp<=NodeZ_NB(ic,jc,kc2,iBLK) )then
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
                factor2=(Zp-NodeZ_NB(ic,jc,kc1,iBLK))/ &
                     (NodeZ_NB(ic,jc,kc2,iBLK)-NodeZ_NB(ic,jc,kc1,iBLK))
                factor1=1.-factor2
                write(unit_tmp,fmt="(30(E14.6))") &
                     (factor1*NodeXYZ_N( ic,jc,kc1,:)+ &
                      factor2*NodeXYZ_N( ic,jc,kc2,:)), &
                     (factor1*PlotVarNodes_NBI(ic,jc,kc1,iBLK,1:nPlotVar)+ &
                      factor2*PlotVarNodes_NBI(ic,jc,kc2,iBLK,1:nPlotVar))
                if(okdebug)write(*,*)'  i=',ic,' j=',jc,' k=',kc1,'-',kc2
             end if
          end if
       end do; end do
    end if

  end subroutine find_cuts

  subroutine write_auxdata
    use ModMultiFluid, ONLY: IonFirst_
    character(len=8)  :: real_date
    character(len=10) :: real_time

    !BLOCKS
    write(stmp,'(i12,3(a,i2))')nBlockALL,'  ',nI,' x',nJ,' x',nK
    write(unit_tmp,'(a,a,a)') 'AUXDATA BLOCKS="',trim(adjustl(stmp)),'"'

    !BODYDENSITY
    write(stmp,'(f12.2)')BodyNDim_I(IonFirst_)
    write(unit_tmp,'(a,a,a)') &
         'AUXDATA BODYNUMDENSITY="',trim(adjustl(stmp)),'"'

    !^CFG IF BORISCORR BEGIN
    !BORIS
    if(boris_correction)then
       write(stmp,'(a,f8.4)')'T ',boris_cLIGHT_factor
    else
       write(stmp,'(a)')'F'
    end if
    write(unit_tmp,'(a,a,a)') 'AUXDATA BORIS="',trim(adjustl(stmp)),'"'
    !^CFG END BORISCORR

    !BTHETATILT
    write(stmp,'(f12.4)')ThetaTilt*cRadToDeg
    write(unit_tmp,'(a,a,a)') 'AUXDATA BTHETATILT="',trim(adjustl(stmp)),'"'

    !CELLS
    write(stmp,'(i12)')nBlockALL*nI*nJ*nK
    write(unit_tmp,'(a,a,a)') 'AUXDATA CELLS="',trim(adjustl(stmp)),'"'

    !CODEVERSION
    write(stmp,'(a,f5.2)')'BATSRUS',CodeVersion
    write(unit_tmp,'(a,a,a)') 'AUXDATA CODEVERSION="',trim(adjustl(stmp)),'"'

    !COORDSYSTEM
    write(stmp,'(a)')TypeCoordSystem
    write(unit_tmp,'(a,a,a)') 'AUXDATA COORDSYSTEM="',trim(adjustl(stmp)),'"'

    !COROTATION
    if(UseRotatingBc)then
       write(stmp,'(a)')'T'
    else
       write(stmp,'(a)')'F'
    end if
    write(unit_tmp,'(a,a,a)') 'AUXDATA COROTATION="',trim(adjustl(stmp)),'"'

    !FLUXTYPE
    write(stmp,'(a)')FluxType
    write(unit_tmp,'(a,a,a)') 'AUXDATA FLUXTYPE="',trim(adjustl(stmp)),'"'

    !GAMMA
    write(stmp,'(f14.6)')g
    write(unit_tmp,'(a,a,a)') 'AUXDATA GAMMA="',trim(adjustl(stmp)),'"'

    !ITER
    write(stmp,'(i12)')n_step
    write(unit_tmp,'(a,a,a)') 'AUXDATA ITER="',trim(adjustl(stmp)),'"'

    !NPROC
    write(stmp,'(i12)')nProc
    write(unit_tmp,'(a,a,a)') 'AUXDATA NPROC="',trim(adjustl(stmp)),'"'

    !ORDER
    if(nORDER==2)then
       write(stmp,'(i12,a,f8.5)') &
            nOrder,' '//trim(TypeLimiter)//', beta=',BetaLimiter
    else
       write(stmp,'(i12)')nORDER
    end if
    write(unit_tmp,'(a,a,a)') 'AUXDATA ORDER="',trim(adjustl(stmp)),'"'

    !RBODY
    write(stmp,'(f12.2)')rBody
    write(unit_tmp,'(a,a,a)') 'AUXDATA RBODY="',trim(adjustl(stmp)),'"'

    !SAVEDATE
    call Date_and_time (real_date, real_time)
    write(stmp,'(a11,a4,a1,a2,a1,a2, a4,a2,a1,a2,a1,a2)') &
         'Save Date: ', real_date(1:4),'/',real_date(5:6),'/',real_date(7:8), &
         ' at ',  real_time(1:2),':',real_time(3:4),':',real_time(5:6)
    write(unit_tmp,'(a,a,a)') 'AUXDATA SAVEDATE="',trim(adjustl(stmp)),'"'

    !TIMEEVENT
    write(stmp,'(a)')textDateTime
    write(unit_tmp,'(a,a,a)') 'AUXDATA TIMEEVENT="',trim(adjustl(stmp)),'"'

    !TIMEEVENTSTART
    write(stmp,'(a)')textDateTime0
    write(unit_tmp,'(a,a,a)') 'AUXDATA TIMEEVENTSTART="',trim(adjustl(stmp)),'"'

    !TIMESIM
    if(time_accurate)then
       write(stmp,'(a)')'T='// &
            StringDateOrTime(1:4)//":"// &
            StringDateOrTime(5:6)//":"// &
            StringDateOrTime(7:8)
    else
       write(stmp,'(a)')'T= N/A'
    end if
    write(unit_tmp,'(a,a,a)') 'AUXDATA TIMESIM="',trim(adjustl(stmp)),'"'

  end subroutine write_auxdata

end subroutine write_plot_tec
