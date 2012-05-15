!^CFG COPYRIGHT UM

module ModMessagePass

  implicit none

  logical:: DoOneCoarserLayer = .true.

contains
  ! moved form file exchange_messages.f90 
  subroutine exchange_messages(DoResChangeOnlyIn, UseOrder2In)
    use ModProcMH
    use ModMain, ONLY : nI, nJ, nK, gcn, nBlockMax, nBlock, unusedBLK, &
         TypeBc_I, time_loop, UseB, &
         UseConstrainB,&              !^CFG IF CONSTRAINB 
         UseProjection,&              !^CFG IF PROJECTION
         UseDivbDiffusion,&           !^CFG IF DIVBDIFFUSE
         time_simulation,nOrder,prolong_order,optimize_message_pass
    use ModVarIndexes
    use ModAdvance, ONLY : State_VGB,divB1_GB
    use ModGeometry, ONLY : far_field_BCs_BLK        
    use ModPhysics, ONLY : ShockSlope
    use ModFaceValue,ONLY: UseAccurateResChange
    use ModEnergy,   ONLY: calc_energy_ghost, correctP

    use BATL_lib, ONLY: message_pass_cell
    use ModParallel, ONLY: BLKneighborLEV ! we should use BATL:DiNeiLevel_IIIB
    use ModAMR, ONLY: DoProfileAmr
    use ModMpi

    logical, optional, intent(in) :: DoResChangeOnlyIn, UseOrder2In

    integer :: iBlock
    logical :: DoRestrictFace, DoOneLayer, DoTwoCoarseLayers
    logical :: DoCorners, DoFaces
    logical :: UseOrder2=.false.
    integer :: nWidth, nCoarseLayer
    logical :: DoResChangeOnly

    logical:: DoTest, DoTestMe, DoTime, DoTimeMe
    character (len=*), parameter :: NameSub = 'exchange_messages'
    !--------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)
    call set_oktest('time_exchange', DoTime, DoTimeMe)

    !!^CFG IF DEBUGGING BEGIN
    ! call testmessage_pass_nodes
    ! call time_message_passing
    !!^CFG END DEBUGGING

    DoResChangeOnly = .false.
    if(present(DoResChangeOnlyIn)) DoResChangeOnly = DoResChangeOnlyIn

    UseOrder2=.false.
    if(present(UseOrder2In)) UseOrder2 = UseOrder2In

    DoRestrictFace = prolong_order==1
    if(UseConstrainB) DoRestrictFace = .false.   !^CFG IF CONSTRAINB

    DoTwoCoarseLayers = &
         nOrder==2 .and. prolong_order==1 .and. .not. DoOneCoarserLayer


    if(DoTestMe)write(*,*) NameSub, &
         ': DoResChangeOnly, UseOrder2, DoRestrictFace, DoTwoCoarseLayers=',&
         DoResChangeOnly, UseOrder2, DoRestrictFace, DoTwoCoarseLayers

    call timing_start('exch_msgs')
    ! Ensure that energy and pressure are consistent and positive in real cells
    !if(prolong_order==2)then     !^CFG IF NOT PROJECTION
    if(.not.DoResChangeOnly) then
       do iBlock = 1, nBlock
          if (unusedBLK(iBlock)) CYCLE
          if (far_field_BCs_BLK(iBlock) .and. prolong_order==2)&
               call set_outer_BCs(iBlock,time_simulation,.false.)        
          if(UseConstrainB)call correctP(iBlock)   !^CFG IF CONSTRAINB
          if(UseProjection)call correctP(iBlock)   !^CFG IF PROJECTION
       end do
       !end if                       !^CFG IF NOT PROJECTION
    end if

    if (UseOrder2) then
       call message_pass_cell(nVar, State_VGB,&
            DoResChangeOnlyIn=DoResChangeOnlyIn)
       if(.not.DoResChangeOnly) &
            call fix_boundary_ghost_cells(DoRestrictFace)
    elseif (optimize_message_pass=='all') then
       ! If ShockSlope is not zero then even the first order scheme needs 
       ! two ghost cell layers to fill in the corner cells at the sheared BCs.
       DoOneLayer = nOrder == 1 .and. ShockSlope == 0.0

       nWidth = 2;       if(DoOneLayer)        nWidth = 1
       nCoarseLayer = 1; if(DoTwoCoarseLayers) nCoarseLayer = 2
       call message_pass_cell(nVar, State_VGB, &
            nWidthIn=nWidth, nProlongOrderIn=1, &
            nCoarseLayerIn=nCoarseLayer, DoRestrictFaceIn = DoRestrictFace,&
            DoResChangeOnlyIn=DoResChangeOnlyIn)
       if(.not.DoResChangeOnly) &
            call fix_boundary_ghost_cells(DoRestrictFace)
    else
       ! Do not pass corners if not necessary
       DoFaces = .not.(nOrder == 2 .and. UseAccurateResChange)
       ! Pass one layer if possible
       DoOneLayer = nOrder == 1
       nWidth = 2;       if(DoOneLayer)        nWidth = 1
       nCoarseLayer = 1; if(DoTwoCoarseLayers) nCoarseLayer = 2
       call message_pass_cell(nVar, State_VGB, &
            nWidthIn=nWidth, &
            nProlongOrderIn=1, &
            nCoarseLayerIn=nCoarseLayer,&
            DoSendCornerIn=.not.DoFaces, &
            DoRestrictFaceIn=DoRestrictFace,&
            DoResChangeOnlyIn=DoResChangeOnlyIn)
       if(.not.DoResChangeOnly) &
            call fix_boundary_ghost_cells(DoRestrictFace)
    end if

    if(DoProfileAmr) call timing_start('E and P')

    do iBlock = 1, nBlock
       if (unusedBLK(iBlock)) CYCLE

       ! The corner ghost cells outside the domain are updated
       ! from the ghost cells inside the domain, so the outer 
       ! boundary condition have to be reapplied.
       if(.not.DoResChangeOnly &
            .or. any(abs(BLKneighborLEV(:,:,:,iBlock)) == 1) )then
          if (far_field_BCs_BLK(iBlock)) &
               call set_outer_BCs(iBlock, time_simulation, .false.) 
          if(time_loop.and. any(TypeBc_I=='buffergrid'))&
               call fill_in_from_buffer(iBlock)
       end if

       call calc_energy_ghost(iBlock, DoResChangeOnlyIn=DoResChangeOnlyIn)
    end do

    if(DoProfileAmr) call timing_stop('E and P')

    call timing_stop('exch_msgs')
    if(DoTime)call timing_show('exch_msgs',1)

    if(DoTestMe)write(*,*) NameSub,' finished'

  end subroutine exchange_messages

  !============================================================================
  ! moved form file exchange_messages.f90 
  subroutine fill_in_from_buffer(iBlock)
    use ModGeometry,ONLY:R_BLK,x_BLK,y_BLK,z_BLK
    use ModMain,ONLY:nI,nJ,nK,gcn,rBuffMin,rBuffMax,nDim,x_,y_,z_
    use ModAdvance,ONLY:nVar,State_VGB,rho_,rhoUx_,rhoUz_,Ux_,Uz_
    use ModProcMH,ONLY:iProc
    implicit none
    integer,intent(in)::iBlock
    integer::i,j,k
    real::X_D(nDim)
    logical::DoWrite=.true.
    if(DoWrite)then
       DoWrite=.false.
       if(iProc==0)then
          write(*,*)'Fill in the cells near the inner boundary from the buffer'
       end if
    end if

    do k=1-gcn,nK+gcn
       do j=1-gcn,nJ+gcn
          do i=1-gcn,nI+gcn
             if(R_BLK(i,j,k,iBlock)>rBuffMax.or.R_BLK(i,j,k,iBlock)<rBuffMin)&
                  CYCLE
             X_D(x_)=x_BLK(i,j,k,iBlock)
             X_D(y_)=y_BLK(i,j,k,iBlock)
             X_D(z_)=z_BLK(i,j,k,iBlock)
             !Get interpolated values from buffer grid:
             call get_from_spher_buffer_grid(&
                  X_D,nVar,State_VGB(:,i,j,k,iBlock))
             !Transform primitive variables to conservative ones:
             State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)=&
                  State_VGB(Ux_:Uz_,i,j,k,iBlock)*&
                  State_VGB(rho_   ,i,j,k,iBlock)
          end do
       end do
    end do
  end subroutine fill_in_from_buffer

end module ModMessagePass
