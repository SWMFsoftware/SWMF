!^CFG COPYRIGHT UM
subroutine write_plot_tec(ifile,nplotvar,plotvarnodes,nplotvarmax,unitstr_TEC,&
     xmin,xmax,ymin,ymax,zmin,zmax)
  !
  !NOTE: This routine assumes that the blocks are sorted on PEs by their global
  !       block number, ie blocks 1 to n on PE 0, blocks n+1 to n+m on PE 1,
  !       etc.
  !
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,globalBLK,global_block_number, &
       nBlockALL,nBlockMax, StringTimeH4M2S2,time_accurate,n_step
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,true_cell
  use ModParallel, ONLY : UseCorners, iBlock_A, iProc_A
  use ModPhysics, ONLY : unitUSER_x
  use ModIO
  use ModNodes
  use ModMpi
  implicit none

  ! Arguments  
  integer, intent(in) :: ifile, nplotvar, nplotvarmax
  character (LEN=500), intent(in) :: unitstr_TEC
  real, intent(in) :: PlotVarNodes(0:nI,0:nJ,0:nK,nBLK,nplotvarmax)
  real, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax

  ! Local Variables
  integer :: i,j,k, cut1,cut2, iPE,iBLK, iBlockALL, nBlockCuts, iError
  real :: CutValue, factor1,factor2
  logical :: oktest,oktest_me
  integer, allocatable, dimension(:) :: BlockCut
  character (len=22) :: textNandT
  character (len=23) :: textDateTime

  integer :: iTime_I(7)
  !----------------------------------------------------------------------------

  call set_oktest('write_plot_tec',oktest,oktest_me)
  if(oktest_me)write(*,*) plot_type1,plot_type1(1:3)  

  ! Create text string for zone name like 'N=0002000 T=0000:05:00'
  if(time_accurate)then
     call get_time_string
     write(textNandT,'(a,i7.7,a)') "N=",n_step," T="// &
          StringTimeH4M2S2(1:4)//":"// &
          StringTimeH4M2S2(5:6)//":"// &
          StringTimeH4M2S2(7:8)
  else
     write(textNandT,'(a,i7.7)') &
          "N=",n_step
  end if

  call get_date_time(iTime_I)
  write(textDateTime,'(i4.4,"-",5(i2.2,"-"),i3.3)') iTime_I

  select case(plot_type1(1:3))
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
     end if
     do iBlockALL  = 1, nBlockALL
        iBLK = iBlock_A(iBlockALL)
        iPE  = iProc_A(iBlockALL)
        if(iProc==iPE)then
           !================================= 3d ============================
           ! Write point values
           do k=0,nK; do j=0,nJ; do i=0,nI
              if(NodeUniqueGlobal_IIIB(i,j,k,iBLK))then
                 if (plot_dimensional(ifile)) then
                    write(unit_tmp,fmt="(30(E14.6))") &
                         NodeX_IIIB(i,j,k,iBLK)*unitUSER_x, &
                         NodeY_IIIB(i,j,k,iBLK)*unitUSER_x, &
                         NodeZ_IIIB(i,j,k,iBLK)*unitUSER_x, &
                         PlotVarNodes(i,j,k,iBLK,1:nplotvar)
                 else
                    write(unit_tmp,fmt="(30(E14.6))") &
                         NodeX_IIIB(i,j,k,iBLK), &
                         NodeY_IIIB(i,j,k,iBLK), &
                         NodeZ_IIIB(i,j,k,iBLK), &
                         PlotVarNodes(i,j,k,iBLK,1:nplotvar)
                 end if
              end if
           end do; end do; end do
           ! Write point connectivity
           do k=1,nK; do j=1,nJ; do i=1,nI
              write(unit_tmp2,'(8(i8,1x))') &
                   NodeNumberGlobal_IIIB(i-1,j-1,k-1,iBLK), &
                   NodeNumberGlobal_IIIB(i  ,j-1,k-1,iBLK), &
                   NodeNumberGlobal_IIIB(i  ,j  ,k-1,iBLK), &
                   NodeNumberGlobal_IIIB(i-1,j  ,k-1,iBLK), &
                   NodeNumberGlobal_IIIB(i-1,j-1,k  ,iBLK), &
                   NodeNumberGlobal_IIIB(i  ,j-1,k  ,iBLK), &
                   NodeNumberGlobal_IIIB(i  ,j  ,k  ,iBLK), &
                   NodeNumberGlobal_IIIB(i-1,j  ,k  ,iBLK)
           end do; end do; end do
        end if
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
        ! First loop to count nodes and cells
        do iBlockALL  = 1, nBlockALL
           iBLK = iBlock_A(iBlockALL)
           iPE  = iProc_A(iBlockALL)
           if(iProc==iPE)then
              if ( CutValue> NodeX_IIIB(0 ,0,0,iBLK) .and. &
                   CutValue<=NodeX_IIIB(nI,0,0,iBLK) )then
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
        end if
        ! Now loop to write values
        do iBlockALL  = 1, nBlockALL
           iBLK = iBlock_A(iBlockALL)
           iPE  = iProc_A(iBlockALL)
           if(iProc==iPE)then
              if ( CutValue> NodeX_IIIB(0 ,0,0,iBLK) .and. &
                   CutValue<=NodeX_IIIB(nI,0,0,iBLK) )then
                 ! Find cut interpolation factors
                 do i=1,nI
                    if ( CutValue> NodeX_IIIB(i-1,0,0,iBLK) .and. &
                         CutValue<=NodeX_IIIB(i  ,0,0,iBLK)  )then
                       cut1=i-1
                       cut2=i
                       factor2=(CutValue-NodeX_IIIB(i-1,0,0,iBLK))/ &
                            (NodeX_IIIB(i,0,0,iBLK)-NodeX_IIIB(i-1,0,0,iBLK))
                       factor1=1.-factor2
                       EXIT
                    end if
                 end do
                 ! Write point values
                 do k=0,nK; do j=0,nJ
                    if (plot_dimensional(ifile)) then
                       write(unit_tmp,fmt="(30(E14.6))") &
                            (factor1*NodeX_IIIB(cut1,j,k,iBLK)+factor2*NodeX_IIIB(cut2,j,k,iBLK))*unitUSER_x, &
                            (factor1*NodeY_IIIB(cut1,j,k,iBLK)+factor2*NodeY_IIIB(cut2,j,k,iBLK))*unitUSER_x, &
                            (factor1*NodeZ_IIIB(cut1,j,k,iBLK)+factor2*NodeZ_IIIB(cut2,j,k,iBLK))*unitUSER_x, &
                            (factor1*PlotVarNodes(cut1,j,k,iBLK,1:nplotvar)+ &
                            factor2*PlotVarNodes(cut2,j,k,iBLK,1:nplotvar))
                    else
                       write(unit_tmp,fmt="(30(E14.6))") &
                            (factor1*NodeX_IIIB(cut1,j,k,iBLK)+factor2*NodeX_IIIB(cut2,j,k,iBLK)), &
                            (factor1*NodeY_IIIB(cut1,j,k,iBLK)+factor2*NodeY_IIIB(cut2,j,k,iBLK)), &
                            (factor1*NodeZ_IIIB(cut1,j,k,iBLK)+factor2*NodeZ_IIIB(cut2,j,k,iBLK)), &
                            (factor1*PlotVarNodes(cut1,j,k,iBLK,1:nplotvar)+ &
                            factor2*PlotVarNodes(cut2,j,k,iBLK,1:nplotvar))
                    end if
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
     elseif((ymax-ymin)<(zmax-zmin))then
        !Y Slice
        CutValue = 0.5*(ymin+ymax)
        ! First loop to count nodes and cells
        do iBlockALL  = 1, nBlockALL
           iBLK = iBlock_A(iBlockALL)
           iPE  = iProc_A(iBlockALL)
           if(iProc==iPE)then
              if ( CutValue> NodeY_IIIB(0,0 ,0,iBLK) .and. &
                   CutValue<=NodeY_IIIB(0,nJ,0,iBLK) )then
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
        end if
        ! Now loop to write values
        do iBlockALL  = 1, nBlockALL
           iBLK = iBlock_A(iBlockALL)
           iPE  = iProc_A(iBlockALL)
           if(iProc==iPE)then
              if ( CutValue> NodeY_IIIB(0,0 ,0,iBLK) .and. &
                   CutValue<=NodeY_IIIB(0,nJ,0,iBLK) )then
                 ! Find cut interpolation factors
                 do j=1,nJ
                    if ( CutValue> NodeY_IIIB(0,j-1,0,iBLK) .and. &
                         CutValue<=NodeY_IIIB(0,j  ,0,iBLK)  )then
                       cut1=j-1
                       cut2=j
                       factor2=(CutValue-NodeY_IIIB(0,j-1,0,iBLK))/ &
                            (NodeY_IIIB(0,j,0,iBLK)-NodeY_IIIB(0,j-1,0,iBLK))
                       factor1=1.-factor2
                       EXIT
                    end if
                 end do
                 ! Write point values
                 do k=0,nK; do i=0,nI
                    if (plot_dimensional(ifile)) then
                       write(unit_tmp,fmt="(30(E14.6))") &
                            (factor1*NodeX_IIIB(i,cut1,k,iBLK)+factor2*NodeX_IIIB(i,cut2,k,iBLK))*unitUSER_x, &
                            (factor1*NodeY_IIIB(i,cut1,k,iBLK)+factor2*NodeY_IIIB(i,cut2,k,iBLK))*unitUSER_x, &
                            (factor1*NodeZ_IIIB(i,cut1,k,iBLK)+factor2*NodeZ_IIIB(i,cut2,k,iBLK))*unitUSER_x, &
                            (factor1*PlotVarNodes(i,cut1,k,iBLK,1:nplotvar)+ &
                            factor2*PlotVarNodes(i,cut2,k,iBLK,1:nplotvar))
                    else
                       write(unit_tmp,fmt="(30(E14.6))") &
                            (factor1*NodeX_IIIB(i,cut1,k,iBLK)+factor2*NodeX_IIIB(i,cut2,k,iBLK)), &
                            (factor1*NodeY_IIIB(i,cut1,k,iBLK)+factor2*NodeY_IIIB(i,cut2,k,iBLK)), &
                            (factor1*NodeZ_IIIB(i,cut1,k,iBLK)+factor2*NodeZ_IIIB(i,cut2,k,iBLK)), &
                            (factor1*PlotVarNodes(i,cut1,k,iBLK,1:nplotvar)+ &
                            factor2*PlotVarNodes(i,cut2,k,iBLK,1:nplotvar))
                    end if
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
     else
        !Z Slice
        CutValue = 0.5*(zmin+zmax)
        ! First loop to count nodes and cells
        do iBlockALL  = 1, nBlockALL
           iBLK = iBlock_A(iBlockALL)
           iPE  = iProc_A(iBlockALL)
           if(iProc==iPE)then
              if ( CutValue> NodeZ_IIIB(0,0, 0,iBLK) .and. &
                   CutValue<=NodeZ_IIIB(0,0,nK,iBLK) )then
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
        end if
        ! Now loop to write values
        do iBlockALL  = 1, nBlockALL
           iBLK = iBlock_A(iBlockALL)
           iPE  = iProc_A(iBlockALL)
           if(iProc==iPE)then
              if ( CutValue> NodeZ_IIIB(0,0, 0,iBLK) .and. &
                   CutValue<=NodeZ_IIIB(0,0,nK,iBLK) )then
                 ! Find cut interpolation factors
                 do k=1,nK
                    if ( CutValue> NodeZ_IIIB(0,0,k-1,iBLK) .and. &
                         CutValue<=NodeZ_IIIB(0,0,k  ,iBLK)  )then
                       cut1=k-1
                       cut2=k
                       factor2=(CutValue-NodeZ_IIIB(0,0,k-1,iBLK))/ &
                            (NodeZ_IIIB(0,0,k,iBLK)-NodeZ_IIIB(0,0,k-1,iBLK))
                       factor1=1.-factor2
                       EXIT
                    end if
                 end do
                 ! Write point values
                 do j=0,nJ; do i=0,nI
                    if (plot_dimensional(ifile)) then
                       write(unit_tmp,fmt="(30(E14.6))") &
                            (factor1*NodeX_IIIB(i,j,cut1,iBLK)+factor2*NodeX_IIIB(i,j,cut2,iBLK))*unitUSER_x, &
                            (factor1*NodeY_IIIB(i,j,cut1,iBLK)+factor2*NodeY_IIIB(i,j,cut2,iBLK))*unitUSER_x, &
                            (factor1*NodeZ_IIIB(i,j,cut1,iBLK)+factor2*NodeZ_IIIB(i,j,cut2,iBLK))*unitUSER_x, &
                            (factor1*PlotVarNodes(i,j,cut1,iBLK,1:nplotvar)+ &
                            factor2*PlotVarNodes(i,j,cut2,iBLK,1:nplotvar))
                    else
                       write(unit_tmp,fmt="(30(E14.6))") &
                            (factor1*NodeX_IIIB(i,j,cut1,iBLK)+factor2*NodeX_IIIB(i,j,cut2,iBLK)), &
                            (factor1*NodeY_IIIB(i,j,cut1,iBLK)+factor2*NodeY_IIIB(i,j,cut2,iBLK)), &
                            (factor1*NodeZ_IIIB(i,j,cut1,iBLK)+factor2*NodeZ_IIIB(i,j,cut2,iBLK)), &
                            (factor1*PlotVarNodes(i,j,cut1,iBLK,1:nplotvar)+ &
                            factor2*PlotVarNodes(i,j,cut2,iBLK,1:nplotvar))
                    end if
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
  case default
     write(*,*)'Error in write_plot_tec: Unknown plot_type='//plot_type1
  end select

end subroutine write_plot_tec
