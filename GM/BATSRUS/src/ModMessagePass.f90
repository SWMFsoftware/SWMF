!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModMessagePass

  implicit none

  logical:: DoOneCoarserLayer = .true.

contains
  ! moved form file exchange_messages.f90 
  subroutine exchange_messages(DoResChangeOnlyIn, UseOrder2In)

    use ModCellBoundary, ONLY: set_cell_boundary, set_edge_corner_ghost
    use ModProcMH
    use ModMain, ONLY : nBlock, Unused_B, &
         TypeBc_I, time_loop, &
         UseConstrainB, UseProjection, &
         nOrder, nOrderProlong, optimize_message_pass, &
         UseHighResChange
    use ModVarIndexes
    use ModAdvance,  ONLY: State_VGB
    use ModGeometry, ONLY: far_field_BCs_BLK        
    use ModPhysics,  ONLY: ShockSlope, nVectorVar, iVectorVar_I
    use ModFaceValue,ONLY: UseAccurateResChange
    use ModEnergy,   ONLY: calc_energy_ghost, correctP
    use ModCoordTransform, ONLY: rot_xyz_sph

    use BATL_lib, ONLY: message_pass_cell, DiLevelNei_IIIB, nG, &
         MinI, MaxI, MinJ, MaxJ, MinK, MaxK, Xyz_DGB, &
         IsSpherical, IsRLonLat, IsPeriodic_D, IsPeriodicCoord_D
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

    logical :: UseHighResChangeNow

    !!! TO BE GENERALIZED
    logical:: IsPeriodicWedge = .false.
    integer:: iVector, iVar, i, j, k
    real   :: XyzSph_DD(3,3)

    logical:: DoTest, DoTestMe, DoTime, DoTimeMe
    character (len=*), parameter :: NameSub = 'exchange_messages'
    !--------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)
    call set_oktest('time_exchange', DoTime, DoTimeMe)

    ! This way of doing periodic BC for wedge is not perfect.
    ! It does not work for AMR or semi-implicit scheme with vectors.
    ! But it works for a number of simple but useful applications.
    ! A periodic wedge BC is needed if there is a periodic boundary condition
    ! for a non-periodic angular coordinate.
    ! In this case the vector variables are convered to spherical components
    ! during the message passing.
    IsPeriodicWedge = (IsSpherical .or. IsRLonLat) .and. &
         (any(IsPeriodic_D(2:3) .and. .not. IsPeriodicCoord_D(2:3)))

    DoResChangeOnly = .false.
    if(present(DoResChangeOnlyIn)) DoResChangeOnly = DoResChangeOnlyIn

    UseOrder2=.false.
    if(present(UseOrder2In)) UseOrder2 = UseOrder2In

    DoRestrictFace = nOrderProlong==1
    if(UseConstrainB) DoRestrictFace = .false.

    DoTwoCoarseLayers = &
         nOrder>1 .and. nOrderProlong==1 .and. .not. DoOneCoarserLayer

    UseHighResChangeNow = nOrder==5 .and. UseHighResChange

    if(DoTestMe)write(*,*) NameSub, &
         ': DoResChangeOnly, UseOrder2, DoRestrictFace, DoTwoCoarseLayers=',&
         DoResChangeOnly, UseOrder2, DoRestrictFace, DoTwoCoarseLayers

    call timing_start('exch_msgs')

    ! Ensure that energy and pressure are consistent and positive in real cells
    if(.not.DoResChangeOnly) then
       do iBlock = 1, nBlock
          if (Unused_B(iBlock)) CYCLE
          if (far_field_BCs_BLK(iBlock) .and. &
               (nOrderProlong==2 .or. UseHighResChangeNow)) then
             call set_cell_boundary&
                  (nG, iBlock, nVar, State_VGB(:,:,:,:,iBlock))
             if(UseHighResChangeNow) &
                  call set_edge_corner_ghost&
                  (nG,iBlock,nVar,State_VGB(:,:,:,:,iBlock))
          endif
          if(UseConstrainB)call correctP(iBlock)
          if(UseProjection)call correctP(iBlock)
       end do
    end if

    if(IsPeriodicWedge)then
       do iBlock = 1, nBlock
          if (Unused_B(iBlock)) CYCLE
          ! Skip blocks not at the boundary !!!
          do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
             XyzSph_DD = rot_xyz_sph(Xyz_DGB(:,i,j,k,iBlock))
             do iVector = 1, nVectorVar
                iVar = iVectorVar_I(iVector)
                State_VGB(iVar:iVar+2,i,j,k,iBlock) = &
                     matmul(State_VGB(iVar:iVar+2,i,j,k,iBlock), XyzSph_DD)
             end do
          end do; end do; enddo
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
            UseHighResChangeIn=UseHighResChangeNow,&
            DefaultState_V=DefaultState_V)
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
            UseHighResChangeIn=UseHighResChangeNow,&
            DefaultState_V=DefaultState_V)
       if(.not.DoResChangeOnly) call fix_boundary_ghost_cells
    end if

    if(IsPeriodicWedge)then
       do iBlock = 1, nBlock
          if (Unused_B(iBlock)) CYCLE
          ! Skip blocks not at the boundary !!!
          do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
             XyzSph_DD = rot_xyz_sph(Xyz_DGB(:,i,j,k,iBlock))
             do iVector = 1, nVectorVar
                iVar = iVectorVar_I(iVector)
                State_VGB(iVar:iVar+2,i,j,k,iBlock) = &
                     matmul(XyzSph_DD, State_VGB(iVar:iVar+2,i,j,k,iBlock))
             end do
          end do; end do; enddo
       end do
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
