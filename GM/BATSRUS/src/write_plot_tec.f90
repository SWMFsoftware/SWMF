!^CFG COPYRIGHT UM
subroutine write_plot_tec(ifile,nplotvar,plotvarnodes,unitstr_TEC,&
     xmin,xmax,ymin,ymax,zmin,zmax)
  !
  !NOTE: This routine assumes that the blocks are sorted on PEs by their global
  !       block number, ie blocks 1 to n on PE 0, blocks n+1 to n+m on PE 1,
  !       etc.
  !
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,globalBLK,global_block_number, &
       nBlockALL,nBlockMax, StringTimeH4M2S2,time_accurate,n_step,&
       nOrder, limiter_type,betalimiter, UseRotatingBc, &
       TypeCoordSystem, problem_type, StringProblemType_I, CodeVersion
  use ModMain, ONLY: boris_correction                     !^CFG IF BORISCORR
  use ModParallel, ONLY : iBlock_A, iProc_A
  use ModPhysics, ONLY : unitUSER_x, thetaTilt, Rbody, boris_cLIGHT_factor, &
       Body_rho_dim, g
  use ModAdvance, ONLY : FluxType
  use ModIO
  use ModNodes
  use ModNumConst, ONLY : cRadToDeg
  use ModMpi
  implicit none

  ! Arguments  
  integer, intent(in) :: ifile, nplotvar
  character (LEN=500), intent(in) :: unitstr_TEC
  real, intent(in) :: PlotVarNodes(0:nI,0:nJ,0:nK,nBLK,nplotvarmax)
  real, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax

  ! Local Variables
  integer :: i,j,k, cut1,cut2, iPE,iBLK, iBlockALL, nBlockCuts, iError
  real :: CutValue, factor1,factor2
  logical :: oktest,oktest_me
  integer, allocatable, dimension(:) :: BlockCut
  character (len=22) :: textNandT
  character (len=23) :: textDateTime0,textDateTime
  character (len=80) :: format

  integer :: iTime0_I(7),iTime_I(7)
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

  write(format,*)'(i4.4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2,":",i2.2,".",i3.3)'
  call get_date_time_start(iTime0_I)
  call get_date_time(iTime_I)
  write(textDateTime0,format) iTime0_I
  write(textDateTime ,format) iTime_I

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
        call write_auxdata
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
           call write_auxdata
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
           call write_auxdata
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
           call write_auxdata
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

contains

  subroutine write_auxdata
    character(len=8)  :: real_date
    character(len=10) :: real_time
    character(len=80) :: stmp

    !BLOCKS
    write(stmp,'(i12,3(a,i2))')nBlockALL,'  ',nI,' x',nJ,' x',nK
    write(unit_tmp,'(a,a,a)') 'AUXDATA BLOCKS="',trim(adjustl(stmp)),'"'

    !BODYDENSITY
    write(stmp,'(f12.2)')Body_rho_dim
    write(unit_tmp,'(a,a,a)') 'AUXDATA RBODY="',trim(adjustl(stmp)),'"'

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
       if(limiter_type=='beta')then
          write(stmp,'(i12,a,e13.5)')nOrder,', beta=',BetaLimiter
       else
          write(stmp,'(i12,a)')nORDER,' '//trim(limiter_type)
       end if
    else
       write(stmp,'(i12)')nORDER
    end if
    write(unit_tmp,'(a,a,a)') 'AUXDATA ORDER="',trim(adjustl(stmp)),'"'

    !PROBLEMTYPE
    write(stmp,'(i12,a)')problem_type,' '//trim(StringProblemType_I(problem_type))
    write(unit_tmp,'(a,a,a)') 'AUXDATA PROBLEMTYPE="',trim(adjustl(stmp)),'"'

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
            StringTimeH4M2S2(1:4)//":"// &
            StringTimeH4M2S2(5:6)//":"// &
            StringTimeH4M2S2(7:8)
    else
       write(stmp,'(a)')'T= N/A'
    end if
    write(unit_tmp,'(a,a,a)') 'AUXDATA TIMESIM="',trim(adjustl(stmp)),'"'

  end subroutine write_auxdata

end subroutine write_plot_tec
