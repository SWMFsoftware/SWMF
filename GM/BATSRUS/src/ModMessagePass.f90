!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!This code is a copyright protected software (c) 2002- University of Michigan

module ModMessagePass

  implicit none

  logical:: DoOneCoarserLayer = .true.

contains
  ! moved form file exchange_messages.f90 
  subroutine exchange_messages(DoResChangeOnlyIn, UseOrder2In)

    use ModCellBoundary, ONLY: set_cell_boundary
    use ModProcMH
    use ModMain, ONLY : nBlock, Unused_B, &
         TypeBc_I, time_loop, &
         UseConstrainB, UseProjection, &
         nOrder, nOrderProlong, optimize_message_pass, &
         UseHighResChange
    use ModVarIndexes
    use ModAdvance,  ONLY: State_VGB
    use ModGeometry, ONLY: far_field_BCs_BLK        
    use ModPhysics,  ONLY: ShockSlope
    use ModFaceValue,ONLY: UseAccurateResChange
    use ModEnergy,   ONLY: calc_energy_ghost, correctP

    use BATL_lib, ONLY: message_pass_cell, DiLevelNei_IIIB, nG
    use ModMpi

    ! Fill ghost cells at res. change only
    logical, optional, intent(in) :: DoResChangeOnlyIn 

    ! Use 2nd order prolongation to fill
    logical, optional, intent(in) :: UseOrder2In

    logical :: UseOrder2
    logical :: DoResChangeOnly

    logical :: DoRestrictFace, DoTwoCoarseLayers, DoSendCorner
    integer :: nWidth, nCoarseLayer

    integer :: iBlock

    logical:: DoTest, DoTestMe, DoTime, DoTimeMe
    character (len=*), parameter :: NameSub = 'exchange_messages'
    !--------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)
    call set_oktest('time_exchange', DoTime, DoTimeMe)

    DoResChangeOnly = .false.
    if(present(DoResChangeOnlyIn)) DoResChangeOnly = DoResChangeOnlyIn

    UseOrder2=.false.
    if(present(UseOrder2In)) UseOrder2 = UseOrder2In

    DoRestrictFace = nOrderProlong==1
    if(UseConstrainB) DoRestrictFace = .false.

    DoTwoCoarseLayers = &
         nOrder>1 .and. nOrderProlong==1 .and. .not. DoOneCoarserLayer

    if(DoTestMe)write(*,*) NameSub, &
         ': DoResChangeOnly, UseOrder2, DoRestrictFace, DoTwoCoarseLayers=',&
         DoResChangeOnly, UseOrder2, DoRestrictFace, DoTwoCoarseLayers

    call timing_start('exch_msgs')

    ! Ensure that energy and pressure are consistent and positive in real cells
    if(.not.DoResChangeOnly) then
       do iBlock = 1, nBlock
          if (Unused_B(iBlock)) CYCLE
          if (far_field_BCs_BLK(iBlock) .and. nOrderProlong==2) call &
               set_cell_boundary(nG, iBlock, nVar, State_VGB(:,:,:,:,iBlock))
          if(UseConstrainB)call correctP(iBlock)
          if(UseProjection)call correctP(iBlock)
       end do
    end if

    if (UseOrder2 .or. nOrderProlong > 1) then
       call message_pass_cell(nVar, State_VGB,&
            DoResChangeOnlyIn=DoResChangeOnlyIn)
       if(.not.DoResChangeOnly) call fix_boundary_ghost_cells
    elseif (optimize_message_pass=='all') then
       ! If ShockSlope is not zero then even the first order scheme needs 
       ! all ghost cell layers to fill in the corner cells at the sheared BCs.
       nWidth = nG; if(nOrder == 1 .and. ShockSlope == 0.0)  nWidth = 1
       nCoarseLayer = 1; if(DoTwoCoarseLayers) nCoarseLayer = 2
       call message_pass_cell(nVar, State_VGB, &
            nWidthIn=nWidth, nProlongOrderIn=1, &
            nCoarseLayerIn=nCoarseLayer, DoRestrictFaceIn = DoRestrictFace,&
            DoResChangeOnlyIn=DoResChangeOnlyIn, &
            UseHighResChangeIn=UseHighResChange)
       if(.not.DoResChangeOnly) call fix_boundary_ghost_cells
    else
       ! Pass corners if necessary
       DoSendCorner = nOrder > 1 .and. UseAccurateResChange
       ! Pass one layer if possible
       nWidth = nG;      if(nOrder == 1)       nWidth = 1
       nCoarseLayer = 1; if(DoTwoCoarseLayers) nCoarseLayer = 2
       call message_pass_cell(nVar, State_VGB, &
            nWidthIn=nWidth, &
            nProlongOrderIn=1, &
            nCoarseLayerIn=nCoarseLayer,&
            DoSendCornerIn=DoSendCorner, &
            DoRestrictFaceIn=DoRestrictFace,&
            DoResChangeOnlyIn=DoResChangeOnlyIn,&
            UseHighResChangeIn=UseHighResChange)
       if(.not.DoResChangeOnly) call fix_boundary_ghost_cells
    end if

    do iBlock = 1, nBlock
       if (Unused_B(iBlock)) CYCLE

       ! The corner ghost cells outside the domain are updated
       ! from the ghost cells inside the domain, so the outer 
       ! boundary condition have to be reapplied.
       if(.not.DoResChangeOnly &
            .or. any(abs(DiLevelNei_IIIB(:,:,:,iBlock)) == 1) )then
          if (far_field_BCs_BLK(iBlock)) call set_cell_boundary( &
               nG, iBlock, nVar, State_VGB(:,:,:,:,iBlock))
          if(time_loop.and. any(TypeBc_I=='buffergrid'))&
               call fill_in_from_buffer(iBlock)
       end if

       call calc_energy_ghost(iBlock, DoResChangeOnlyIn=DoResChangeOnlyIn)
    end do

    call timing_stop('exch_msgs')
    if(DoTime)call timing_show('exch_msgs',1)

    if(DoTestMe)write(*,*) NameSub,' finished'

  end subroutine exchange_messages

  !============================================================================
  subroutine fill_in_from_buffer(iBlock)
  
    use ModGeometry,ONLY: R_BLK
    use ModMain,    ONLY: BufferMin_D, BufferMax_D
    use ModAdvance, ONLY: nVar, State_VGB, Rho_, RhoUx_, RhoUz_, Ux_, Uz_
    use ModProcMH,  ONLY: iProc
    use BATL_lib,   ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, Xyz_DGB
    implicit none
    integer,intent(in)::iBlock

    integer:: i, j, k
    logical:: DoWrite=.true.
    !------------------------------------------------------------------------

    if(DoWrite)then
       DoWrite=.false.
       if(iProc==0)then
          write(*,*)'Fill in the cells near the inner boundary from the buffer'
       end if
    end if

    do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
       if(  R_BLK(i,j,k,iBlock) > BufferMax_D(1) .or. &
            R_BLK(i,j,k,iBlock) < BufferMin_D(1))&
            CYCLE

       ! Get interpolated values from buffer grid:
       call get_from_spher_buffer_grid(&
            Xyz_DGB(:,i,j,k,iBlock), nVar, State_VGB(:,i,j,k,iBlock))

       ! Transform primitive variables to conservative ones:
       State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = &
            State_VGB(Rho_,i,j,k,iBlock)*State_VGB(Ux_:Uz_,i,j,k,iBlock)

    end do; end do; end do

  end subroutine fill_in_from_buffer

end module ModMessagePass
