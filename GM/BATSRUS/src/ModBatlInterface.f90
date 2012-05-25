module ModBatlInterface

  implicit none

  logical, public :: UseBatlTest = .true.

contains
  !===========================================================================
  subroutine set_batsrus_grid

    use BATL_lib, ONLY: nNodeUsed, nBlock, Unused_B, Unused_BP, &
         iProc, iComm,&
         IsNewDecomposition, IsNewTree

    use ModMain, ONLY: nBlockAll, nBlockBats => nBlock, nBlockMax, &
         iNewGrid, iNewDecomposition

    use ModPartSteady, ONLY: UsePartSteady

    use ModGeometry, ONLY: dx_BLK, MinDxValue, MaxDxValue

    use ModAdvance, ONLY: iTypeAdvance_B, iTypeAdvance_BP, &
         SkippedBlock_, ExplBlock_
    use ModMpi
    use ModIO, ONLY: restart

    integer:: iBlock, iError
    real   :: DxMin, DxMax
    !-------------------------------------------------------------------------

    ! Tell if the grid and/or the tree has changed
    if(IsNewDecomposition) iNewDecomposition = mod(iNewDecomposition+1,10000)
    if(IsNewTree) iNewGrid = mod( iNewGrid+1, 10000)

    if( IsNewDecomposition .or. IsNewTree .or. restart) then

       ! If regrid_batl is not called in each time step, we say that the
       ! grid/tree is not changed from the view of BATL
       IsNewDecomposition = .false.
       IsNewTree          = .false.

       nBlockAll  = nNodeUsed
       nBlockBats = nBlock
       call MPI_ALLREDUCE(nBlock, nBlockMax, 1, MPI_INTEGER, MPI_MAX, &
            iComm, iError)

       if(.not.UsePartSteady)then
          where(Unused_BP(1:nBlockMax,:))
             iTypeAdvance_BP(1:nBlockMax,:) = SkippedBlock_
          elsewhere
             iTypeAdvance_BP(1:nBlockMax,:) = ExplBlock_
          end where
          iTypeAdvance_B(1:nBlock) = iTypeAdvance_BP(1:nBlock,iProc)
       end if

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          call set_batsrus_block(iBlock)
       end do

       DxMin = minval(dx_BLK(1:nBlock), MASK=(.not.Unused_B(1:nBlock)))
       DxMax = maxval(dx_BLK(1:nBlock), MASK=(.not.Unused_B(1:nBlock)))
       call MPI_allreduce(DxMin, minDXvalue,1,MPI_REAL,MPI_MIN,iComm,iError)
       call MPI_allreduce(DxMax, maxDXvalue,1,MPI_REAL,MPI_MAX,iComm,iError)

    end if

  end subroutine set_batsrus_grid
  !===========================================================================
  subroutine set_batsrus_block(iBlock)

    use BATL_lib, ONLY: nDim, x_, y_, z_, &
         nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         Xyz_DGB, CellVolume_B, CellVolume_GB, &
         CellSize_DB, CoordMin_DB, &
         CellFace_DB, CellFace_DFB, FaceNormal_DDFB, &
         iNode_B, iNodeNei_IIIB, DiLevelNei_IIIB, &
         iTree_IA, Block_, Proc_, Unset_, &
         IsCartesian, IsRzGeometry
    use ModGeometry, ONLY: &
         XyzStart_BLK, &
         x_BLK, y_BLK, z_BLK, r_BLK, rMin_BLK, Cv_BLK, vInv_CB, &
         dx_BLK, dy_BLK, dz_BLK, fax_BLK, fay_BLK, faz_BLK, &
         FaceAreaI_DFB, FaceAreaJ_DFB, FaceAreaK_DFB, &
         FaceArea2MinI_B, FaceArea2MinJ_B, FaceArea2MinK_B

    use ModParallel, ONLY: BLKneighborLEV,  neiLEV, neiBLK, neiPE, &
         neiLeast, neiLwest, neiLsouth, neiLnorth, neiLbot, neiLtop, &
         neiBeast, neiBwest, neiBsouth, neiBnorth, neiBbot, neiBtop, &
         neiPeast, neiPwest, neiPsouth, neiPnorth, neiPbot, neiPtop

    integer, intent(in):: iBlock

    ! Convert from BATL to BATSRUS ordering of subfaces. 

    integer, parameter:: iOrder_I(4) = (/1,3,2,4/)
    integer:: iNodeNei, iNodeNei_I(4)
    integer:: i, j, k
    !-------------------------------------------------------------------------
    BLKneighborLEV(:,:,:,iBlock) = DiLevelNei_IIIB(:,:,:,iBlock)

    neiLeast(iBlock)  = BLKneighborLEV(-1,0,0,iBlock)
    neiLwest(iBlock)  = BLKneighborLEV(+1,0,0,iBlock)
    neiLsouth(iBlock) = BLKneighborLEV(0,-1,0,iBlock)
    neiLnorth(iBlock) = BLKneighborLEV(0,+1,0,iBlock)
    neiLbot(iBlock)   = BLKneighborLEV(0,0,-1,iBlock)
    neiLtop(iBlock)   = BLKneighborLEV(0,0,+1,iBlock)

    neiLEV(1,iBlock)  = neiLeast(iBlock)
    neiLEV(2,iBlock)  = neiLwest(iBlock)
    neiLEV(3,iBlock)  = neiLsouth(iBlock)
    neiLEV(4,iBlock)  = neiLnorth(iBlock)
    neiLEV(5,iBlock)  = neiLbot(iBlock)
    neiLEV(6,iBlock)  = neiLtop(iBlock)

    ! neiBeast ... neiPbot are used in 
    ! ModFaceValue::correct_monotone_restrict and
    ! ModConserveFlux::apply_cons_flux

    select case(DiLevelNei_IIIB(-1,0,0,iBlock))
    case(Unset_)
       neiBeast(:,iBlock)  = Unset_
       neiPeast(:,iBlock)  = Unset_
    case(-1)
       iNodeNei_I = pack(iNodeNei_IIIB(0,1:2,1:2,iBlock),.true.)
       iNodeNei_I = iNodeNei_I(iOrder_I)
       if(nDim < 3) where(iNodeNei_I == Unset_) iNodeNei_I = iNode_B(iBlock)
       neiBeast(:,iBlock)  = iTree_IA(Block_,iNodeNei_I)
       neiPeast(:,iBlock)  = iTree_IA(Proc_,iNodeNei_I)
    case default
       iNodeNei = iNodeNei_IIIB(0,1,1,iBlock)
       neiBeast(:,iBlock)  = iTree_IA(Block_,iNodeNei)
       neiPeast(:,iBlock)  = iTree_IA(Proc_,iNodeNei)
    end select

    select case(DiLevelNei_IIIB(+1,0,0,iBlock))
    case(Unset_)
       neiBwest(:,iBlock)  = Unset_
       neiPwest(:,iBlock)  = Unset_
    case(-1)
       iNodeNei_I = pack(iNodeNei_IIIB(3,1:2,1:2,iBlock),.true.)
       iNodeNei_I = iNodeNei_I(iOrder_I)
       if(nDim < 3) where(iNodeNei_I == Unset_) iNodeNei_I = iNode_B(iBlock)
       neiBwest(:,iBlock)  = iTree_IA(Block_,iNodeNei_I)
       neiPwest(:,iBlock)  = iTree_IA(Proc_,iNodeNei_I)
    case default
       iNodeNei = iNodeNei_IIIB(3,1,1,iBlock)
       neiBwest(:,iBlock)  = iTree_IA(Block_,iNodeNei)
       neiPwest(:,iBlock)  = iTree_IA(Proc_,iNodeNei)
    end select

    select case(DiLevelNei_IIIB(0,-1,0,iBlock))
    case(Unset_)
       neiBsouth(:,iBlock)  = Unset_
       neiPsouth(:,iBlock)  = Unset_
    case(-1)
       iNodeNei_I = pack(iNodeNei_IIIB(1:2,0,1:2,iBlock),.true.)
       iNodeNei_I = iNodeNei_I(iOrder_I)
       if(nDim < 3) where(iNodeNei_I == Unset_) iNodeNei_I = iNode_B(iBlock)
       neiBsouth(:,iBlock)  = iTree_IA(Block_,iNodeNei_I)
       neiPsouth(:,iBlock)  = iTree_IA(Proc_,iNodeNei_I)
    case default
       iNodeNei = iNodeNei_IIIB(1,0,1,iBlock)
       neiBsouth(:,iBlock)  = iTree_IA(Block_,iNodeNei)
       neiPsouth(:,iBlock)  = iTree_IA(Proc_,iNodeNei)
    end select

    select case(DiLevelNei_IIIB(0,+1,0,iBlock))
    case(Unset_)
       neiBnorth(:,iBlock)  = Unset_
       neiPnorth(:,iBlock)  = Unset_
    case(-1)
       iNodeNei_I = pack(iNodeNei_IIIB(1:2,3,1:2,iBlock),.true.)
       iNodeNei_I = iNodeNei_I(iOrder_I)
       if(nDim < 3) where(iNodeNei_I == Unset_) iNodeNei_I = iNode_B(iBlock)
       neiBnorth(:,iBlock)  = iTree_IA(Block_,iNodeNei_I)
       neiPnorth(:,iBlock)  = iTree_IA(Proc_,iNodeNei_I)
    case default
       iNodeNei = iNodeNei_IIIB(1,3,1,iBlock)
       neiBnorth(:,iBlock)  = iTree_IA(Block_,iNodeNei)
       neiPnorth(:,iBlock)  = iTree_IA(Proc_,iNodeNei)
    end select

    select case(DiLevelNei_IIIB(0,0,-1,iBlock))
    case(Unset_ )
       neiBbot(:,iBlock)  = Unset_
       neiPbot(:,iBlock)  = Unset_
    case(-1)
       iNodeNei_I = pack(iNodeNei_IIIB(1:2,1:2,0,iBlock),.true.)
       neiBbot(:,iBlock)  = iTree_IA(Block_,iNodeNei_I)
       neiPbot(:,iBlock)  = iTree_IA(Proc_,iNodeNei_I)
    case default
       iNodeNei = iNodeNei_IIIB(1,1,0,iBlock)
       neiBbot(:,iBlock)  = iTree_IA(Block_,iNodeNei)
       neiPbot(:,iBlock)  = iTree_IA(Proc_,iNodeNei)
    end select

    select case(DiLevelNei_IIIB(0,0,+1,iBlock))
    case(Unset_)
       neiBtop(:,iBlock)  = Unset_
       neiPtop(:,iBlock)  = Unset_
    case(-1)
       iNodeNei_I = pack(iNodeNei_IIIB(1:2,1:2,3,iBlock),.true.)
       neiBtop(:,iBlock)  = iTree_IA(Block_,iNodeNei_I)
       neiPtop(:,iBlock)  = iTree_IA(Proc_,iNodeNei_I)
    case default
       iNodeNei = iNodeNei_IIIB(1,1,3,iBlock)
       neiBtop(:,iBlock)  = iTree_IA(Block_,iNodeNei)
       neiPtop(:,iBlock)  = iTree_IA(Proc_,iNodeNei)
    end select

    ! neiBLK and neiPE are used in ray_pass, constrain_B, ModPartSteady
    neiBLK(:,1,iBlock) = neiBeast(:,iBlock)
    neiBLK(:,2,iBlock) = neiBwest(:,iBlock)
    neiBLK(:,3,iBlock) = neiBsouth(:,iBlock)
    neiBLK(:,4,iBlock) = neiBnorth(:,iBlock)
    neiBLK(:,5,iBlock) = neiBbot(:,iBlock)
    neiBLK(:,6,iBlock) = neiBtop(:,iBlock)

    neiPE(:,1,iBlock) = neiPeast(:,iBlock)
    neiPE(:,2,iBlock) = neiPwest(:,iBlock)
    neiPE(:,3,iBlock) = neiPsouth(:,iBlock)
    neiPE(:,4,iBlock) = neiPnorth(:,iBlock)
    neiPE(:,5,iBlock) = neiPbot(:,iBlock)
    neiPE(:,6,iBlock) = neiPtop(:,iBlock)

    XyzStart_BLK(:,iBlock) = CoordMin_DB(:,iBlock) + 0.5*CellSize_DB(:,iBlock)

    dx_BLK(iBlock) = CellSize_DB(1,iBlock)
    dy_BLK(iBlock) = CellSize_DB(2,iBlock)
    dz_BLK(iBlock) = CellSize_DB(3,iBlock)

    if(UseBatlTest)then
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          x_BLK(i,j,k,iBlock) = Xyz_DGB(1,i,j,k,iBlock)
          y_BLK(i,j,k,iBlock) = Xyz_DGB(2,i,j,k,iBlock)
          z_BLK(i,j,k,iBlock) = Xyz_DGB(3,i,j,k,iBlock)
          r_BLK(i,j,k,iBlock) = sqrt(sum(Xyz_DGB(1:nDim,i,j,k,iBlock)**2))
       end do; end do; end do

       Rmin_BLK(iBlock) = minval(r_BLK(:,:,:,iBlock))

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          vInv_CB(i,j,k,iBlock) = 1/CellVolume_GB(i,j,k,iBlock)
       end do; end do; end do

       if(IsCartesian)then
          fAx_BLK(iBlock) = CellFace_DB(1,iBlock)
          fAy_BLK(iBlock) = CellFace_DB(2,iBlock)
          fAz_BLK(iBlock) = CellFace_DB(3,iBlock)
          cV_BLK(iBlock)  = CellVolume_B(iBlock)
       end if

    end if

    call fix_block_geometry(iBlock)

    if(.not.IsCartesian .and. UseBatlTest .or. IsRzGeometry)then
       FaceArea2MinI_B(iBlock) = 1e-30
       FaceArea2MinJ_B(iBlock) = 1e-30
       FaceArea2MinK_B(iBlock) = 1e-30

       ! Initialize for ignored dimensions
       FaceAreaI_DFB(:,:,:,:,iBlock) = 0.0
       FaceAreaJ_DFB(:,:,:,:,iBlock) = 0.0
       FaceAreaK_DFB(:,:,:,:,iBlock) = 0.0

       if(IsRzGeometry)then
          ! This is like Cartesian except for the areas in R (=x) and Z (=y)
          FaceAreaI_DFB(1,:,:,:,iBlock) = CellFace_DFB(1,:,1:nJ,1:nK,iBlock)
          FaceAreaJ_DFB(2,:,:,:,iBlock) = CellFace_DFB(2,1:nI,:,1:nK,iBlock)
       else
          do k=1, nK; do j=1,nJ; do i=1,nI+1
             FaceAreaI_DFB(1:nDim,i,j,k,iBlock) = &
                  FaceNormal_DDFB(:,x_,i,j,k,iBlock)
          end do; end do; end do
          if(nDim > 1)then
             do k=1, nK; do j=1,nJ+1; do i=1,nI
                FaceAreaJ_DFB(1:nDim,i,j,k,iBlock) = &
                     FaceNormal_DDFB(:,y_,i,j,k,iBlock)
             end do; end do; end do
          end if
          if(nDim > 2)then
             do k=1, nK+1; do j=1,nJ; do i=1,nI
                FaceAreaK_DFB(1:nDim,i,j,k,iBlock) = &
                     FaceNormal_DDFB(:,z_,i,j,k,iBlock)
             end do; end do; end do
          end if

       endif

    end if

  end subroutine set_batsrus_block
  !============================================================================
  subroutine set_batsrus_state

    ! Here we should fix B0 and other things

    use BATL_lib, ONLY: nBlock, iAmrChange_B, AmrMoved_, Unused_B,&
         set_amr_geometry
    use ModEnergy, ONLY: calc_energy_ghost
    use ModResistivity,   ONLY: UseResistivity
    use ModUser,    ONLY : user_specify_refinement

    integer:: iBlock
    !-------------------------------------------------------------------------

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       ! If nothing happened to the block, no need to do anything
       if(iAmrChange_B(iBlock) < AmrMoved_) CYCLE

       ! Update all kinds of extra block variables
       call calc_other_vars(iBlock)
       call calc_energy_ghost(iBlock)
       call set_amr_geometry(iBlock,&
            user_amr_geometry=user_specify_refinement)
       if(UseResistivity) call set_resistivity(iBlock)
    end do

  end subroutine set_batsrus_state
  !============================================================================
  subroutine calc_other_vars(iBlock)

    use ModAdvance,  ONLY: State_VGB, nVar, &
         fbody_x_BLK, fbody_y_BLK, fbody_z_BLK, &
         B0_DGB, B0ResChange_DXSB, B0ResChange_DYSB, B0ResChange_DZSB
    use ModPhysics,  ONLY : CellState_VI
    use ModGeometry, ONLY: body_BLK, true_cell
    use ModMain,     ONLY: west_, TypeBC_I, body1_, UseB0, UseGravity, &
         UseRotatingFrame, globalBLK
    use ModParallel, ONLY: neiLwest, NOBLK
    use ModConserveFlux, ONLY: init_cons_flux
    use BATL_size, ONLY: nI, MinI, MaxI, MinJ, MaxJ, MinK, MaxK

    integer, intent(in) :: iBlock

    integer:: i, j, k
    !--------------------------------------------------------------------------
    ! Initialize variables for flux conservation
    call init_cons_flux(iBlock)

    ! Set B0
    if(UseB0)then
       B0_DGB(:,:,:,:,iBlock) = 0.0
       B0ResChange_DXSB(:,:,:,:,iBlock) = 0.0
       B0ResChange_DYSB(:,:,:,:,iBlock) = 0.0
       B0ResChange_DZSB(:,:,:,:,iBlock) = 0.0

       call set_b0(iBlock)
    end if

    ! Initialize and set body forces
    if(allocated(fbody_x_BLK))then
       fbody_x_BLK(:,:,:,iBlock) = 0.0
       fbody_y_BLK(:,:,:,iBlock) = 0.0
       fbody_z_BLK(:,:,:,iBlock) = 0.0
    end if

    globalBLK = iBlock
    if(UseGravity.or.UseRotatingFrame) call body_force_averages

    ! Set values inside body
    if(body_BLK(iBlock)) then
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          if(true_cell(i,j,k,iBlock)) CYCLE
          State_VGB(1:nVar,i,j,k,iBlock) = CellState_VI(1:nVar,body1_)
       end do;end do; end do     
    end if

    ! For coupled (IH->GM) boundary condition fill in ghost cells
    ! with the first physical cell, because IH may not couple after AMR
    if(TypeBc_I(west_)=='coupled' .and. neiLwest(iBlock)==NOBLK)then
       State_VGB(:,nI+1,:,:,iBlock) = State_VGB(:,nI,:,:,iBlock)
       State_VGB(:,nI+2,:,:,iBlock) = State_VGB(:,nI,:,:,iBlock)
    endif

  end subroutine calc_other_vars

end module ModBatlInterface
