!^CFG COPYRIGHT UM
subroutine get_shifts(iCube,iShift,jShift,kShift)
  use ModCube, ONLY: iShiftChild_DI,nCell2_D
  implicit none
  integer,intent(in) :: iCube
  integer,intent(out):: iShift, jShift, kShift

  iShift = iShiftChild_DI(1,iCube)*nCell2_D(1)
  jShift = iShiftChild_DI(2,iCube)*nCell2_D(2)
  kShift = iShiftChild_DI(3,iCube)*nCell2_D(3)

end subroutine get_shifts


!**************************REFINEMENT**********************************
!\
! For each of the PEs in this group, create the necessary
! refined solution blocks.
!/
subroutine create_refined_soln_blocks(nPEsRefBlk, PEsRefBlk, refined_PE, &
     refined_BLK)

  use ModProcMH
  use ModMain
  use ModVarIndexes,ONLY:nVar
  use ModCT                                             !^CFG IF CONSTRAINB
  use ModAdvance, ONLY :  State_VGB, iTypeAdvance_B
  use ModAMR, ONLY : &
       local_cube,local_cubeBLK
  use ModParallel, ONLY: &
       neiLEV
  use ModGeometry, ONLY : &
       dx_BLK,dy_BLK,dz_BLK,dxyz,xyzStart,xyzStart_BLK,true_cell
  use ModMpi
  use ModEnergy, ONLY: calc_energy_ghost

  implicit none

  integer, intent(in) :: nPEsRefBlk, PEsRefBlk(9), refined_PE, refined_BLK

  integer :: icube, iBLK, isize, iError

  real    :: BxFaceFine_XQS(1:nJ,1:nK,4,2)  !^CFG IF CONSTRAINB BEGIN
  real    :: ByFaceFine_YQS(1:nI,1:nK,4,2)
  real    :: BzFaceFine_ZQS(1:nI,1:nJ,4,2)
  logical :: IsFinerNei_E(east_:top_)       !^CFG END CONSTRAINB

  ! Buffer for broadcasting some scalars to the refined blocks
  integer, parameter :: nSizeBcastBuffer = 8
  real               :: BcastBuffer_I(nSizeBcastBuffer)

  real, dimension(1:nVar,&
      1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn):: sol
  real, dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn):: sol1,sol2,sol3
  logical, dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn) :: INtrue_cell

  !---------------------------------------------------------------------------

  ! Return if processor is not needed here
  if(.not.any(PEsRefBlk(1:nPEsRefBlk)==iProc)) return

  !\
  ! Send coarse block geometry and time stepping info to all PEs owning
  ! the eight newly created refined blocks.
  !/
  if (iProc == refined_PE) then
     BcastBuffer_I(1)   = 0.5*Dt_BLK(refined_BLK)
     BcastBuffer_I(2)   = iTypeAdvance_B(refined_BLK)
     BcastBuffer_I(3)   = Dx_BLK(refined_BLK)
     BcastBuffer_I(4)   = Dy_BLK(refined_BLK)
     BcastBuffer_I(5)   = Dz_BLK(refined_BLK)
     BcastBuffer_I(6:8) = XyzStart_BLK(:,refined_BLK)
  end if
  if (nPEsRefBlk > 1) &
       call cube_bcast_r(nPEsRefBlk,PEsRefBlk,BcastBuffer_I,nSizeBcastBuffer)
  do iCube=1,8
     if(iProc==local_cube(iCube)) then
        iBLK = local_cubeBLK(iCube)
        Dt_BLK(iBLK)         =      BcastBuffer_I(1)
        iTypeAdvance_B(iBLK) = nint(BcastBuffer_I(2))
        Dxyz                 =      BcastBuffer_I(3:5)
        XyzStart             =      BcastBuffer_I(6:8)
     end if
  end do
  !\
  ! If inside the time loop, must also assign solution to
  ! newly created refined blocks.
  !/
  if (time_loop) then
     !\
     ! Send original coarse block solution to all PEs owning
     ! eight newly created refined blocks and prolong this
     ! coarse grid solution to refined solution blocks.
     !/
     isize =(nI+2*gcn)*(nJ+2*gcn)*(nK+2*gcn)

     if (iProc == refined_PE) INtrue_cell = true_cell(:,:,:,refined_BLK)

     !\
     ! Set the geometry parameters for the eight newly 
     ! created refined blocks. This changes true_cell !!!
     !/
     call set_refined_block_geometry

     ! Pass INtrue_cell
     if (nPEsRefBlk > 1) &
          call cube_bcast_l(nPEsRefBlk,PEsRefBlk,INtrue_cell,isize)

     if(nVar>=1)then
        isize=nVar*isize
        ! Pass and prolong rho-p
        if (iProc == refined_PE) &
             sol = State_VGB(&
             1:nVar,:,:,:,refined_BLK)
        if (nPEsRefBlk > 1) &
             call cube_bcast_r(nPEsRefBlk,PEsRefBlk,sol,isize)
        do icube = 1, 8
           if (iProc == local_cube(icube))then
              call prolong_block(2, &
                   INtrue_cell, sol, icube, &
                   State_VGB(&
                   1:nVar,:,:,:,local_cubeBLK(icube)))
           end if
        end do
     end if


     if(UseConstrainB)then                   !^CFG IF CONSTRAINB BEGIN
        isize =(nI+2*gcn)*(nJ+2*gcn)*(nK+2*gcn)
        ! Pass and prolong Bface
        if (iProc == refined_PE)then
           sol1           = Bxface_BLK(:,:,:,refined_BLK)
           sol2           = Byface_BLK(:,:,:,refined_BLK)
           sol3           = Bzface_BLK(:,:,:,refined_BLK)
           IsFinerNei_E   = (neiLEV(:,refined_BLK) == -1)
           if(any(IsFinerNei_E(east_:west_)))&
                BxFaceFine_XQS = BxFaceFine_XQSB(:,:,:,:,refined_BLK)
           if(any(IsFinerNei_E(south_:north_)))&
                ByFaceFine_YQS = ByFaceFine_YQSB(:,:,:,:,refined_BLK)
           if(any(IsFinerNei_E(bot_:top_)))&
                BzFaceFine_ZQS = BzFaceFine_ZQSB(:,:,:,:,refined_BLK)
        end if
        if (nPEsRefBlk > 1)then
           call cube_bcast_r(nPEsRefBlk,PEsRefBlk,sol1,isize)
           call cube_bcast_r(nPEsRefBlk,PEsRefBlk,sol2,isize)
           call cube_bcast_r(nPEsRefBlk,PEsRefBlk,sol3,isize)
           call cube_bcast_l(nPEsRefBlk,PEsRefBlk,IsFinerNei_E,6)
           if(any(IsFinerNei_E(east_:west_)))&
                call cube_bcast_r(nPEsRefBlk,PEsRefBlk,BxFaceFine_XQS,nJ*nK*8)
           if(any(IsFinerNei_E(south_:north_)))&
                call cube_bcast_r(nPEsRefBlk,PEsRefBlk,ByFaceFine_YQS,nI*nK*8)
           if(any(IsFinerNei_E(bot_:top_)))&
                call cube_bcast_r(nPEsRefBlk,PEsRefBlk,BzFaceFine_ZQS,nI*nJ*8)
        end if
        do icube = 1, 8
           if (iProc /= local_cube(icube))CYCLE
           iBLK=local_cubeBLK(icube)
           call prolong_b_face(sol1,sol2,sol3,&
                BxFaceFine_XQS,ByFaceFine_YQS,BzFaceFine_ZQS,&
                IsFinerNei_E,icube,iBLK,&
                BxFace_BLK(:,:,:,iBLK),&
                ByFace_BLK(:,:,:,iBLK),&
                BzFace_BLK(:,:,:,iBLK))
        end do
     end if                                    !^CFG END CONSTRAINB
     !\
     ! For each of these PEs, prolong coarse grid solution
     ! to refined solution blocks.
     !/
     do icube = 1, 8
        if (iProc == local_cube(icube)) then

           iBLK=local_cubeBLK(icube)

           if(UseConstrainB)call Bface2Bcenter(iBLK) !^CFG IF CONSTRAINB

           call fix_soln_block(iBLK)

           call calc_energy_ghost(iBLK)

           call calc_other_soln_vars(iBLK)
        end if
     end do
  else
     !\
     ! If not in time loop, set the geometry parameters for the eight newly 
     ! created refined blocks only.
     !/
     call set_refined_block_geometry
  end if

end subroutine create_refined_soln_blocks

!==============================================================================

subroutine  cube_bcast_r(nPar,Par,var,isize)
  use ModProcMH
  use ModMpi
  implicit none

  integer, intent(in) :: nPar,isize
  integer, intent(in), dimension (9) :: Par
  real, intent(inout), dimension (isize) :: var

  integer :: itag, request, iBLK, i, iError
  integer :: number_send_requests, send_requests(8), &
       status(MPI_STATUS_SIZE, 8), RECstatus(MPI_STATUS_SIZE)  

  number_send_requests = nPar-1
  itag = 10

  if (iProc == Par(1)) then 
     do i=2,nPar
        call MPI_isend(var,isize, MPI_REAL, Par(i), &
             itag, iComm, send_requests(i-1), iError)
     end do

     if (number_send_requests > 0) then
        call MPI_waitall(number_send_requests, &
             send_requests, Status, iError)
     end if
  else
     var = -888888.
     call MPI_recv(var,isize, MPI_REAL, Par(1), &
          itag, iComm, RECstatus, iError)
  end if

end subroutine cube_bcast_r

!==============================================================================

subroutine  cube_bcast_l(nPar,Par,var,isize)
  use ModProcMH
  use ModMpi
  implicit none

  integer, intent(in) :: nPar,isize
  integer, intent(in), dimension (9) :: Par
  logical, intent(inout), dimension (isize) :: var

  integer :: itag, request, iBLK, i, iError
  integer :: number_send_requests, send_requests(8), &
       status(MPI_STATUS_SIZE, 8), RECstatus(MPI_STATUS_SIZE)  

  number_send_requests=nPar-1
  itag = 10

  if (iProc == Par(1)) then 
     do i=2,nPar
        call MPI_isend(var,isize, MPI_LOGICAL, Par(i), &
             itag, iComm, send_requests(i-1), iError)
     end do

     if (number_send_requests > 0) then
        call MPI_waitall(number_send_requests, &
             send_requests, status, iError)
     end if
  else
     var = .false.
     call MPI_recv(var,isize, MPI_LOGICAL, Par(1), &
          itag, iComm, RECstatus, iError)
  end if

end subroutine cube_bcast_l

!==============================================================================

subroutine prolong_block(nOrderIn, IsTrueCellC_G, VarC_VG, iCube, VarF_VG)
  use ModMain
  use ModNumConst
  use ModAdvance,ONLY:nVar
  implicit none

  ! prolong coarse VarC_VG solution into fine VarF_VG

  integer, intent(in) :: nOrderIn, iCube
  logical, intent(in), dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn) :: &
       IsTrueCellC_G
  real, intent(in) :: VarC_VG(nVar,&
       1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn) 
  real, intent(out):: VarF_VG(nVar,&
       1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn) 

  integer :: i,j,k,iF,jF,kF,iShift,jShift,kShift,idir,nTrueNeighbors, iVar
  logical :: UseLimiter_D(3)
  real    :: VarC_DI(nVar,3,-1:1)
  real    :: SlopeRight,SlopeLeft,SlopeSign,SlopeLim_D(nVar,3)
  real,dimension(nVar)    :: GradX_V,GradY_V,GradZ_V

  !---------------------------------------------------------------------------

  call get_shifts(iCube,iShift,jShift,kShift) 
  VarF_VG = cOne ! To fill Ghostcells

  !\
  ! Prolong coarse grid solution to finer block.
  !/
  do k=1+kShift,nK/2+kShift
  do j=1+jShift,nJ/2+jShift
  do i=1+iShift,nI/2+iShift
     iF = 2*(i-iShift)-1
     jF = 2*(j-jShift)-1
     kF = 2*(k-kShift)-1

     ! Averaging for cells whose parent cell was inside the body
     if(.not.IsTrueCellC_G(i,j,k))then
        ! Count true cells in the vicinity of this false cell
        nTrueNeighbors=count(IsTrueCellC_G(i-1:i+1,j,k))+ &
                       count(IsTrueCellC_G(i,j-1:j+1,k))+ &
                       count(IsTrueCellC_G(i,j,k-1:k+1))
        
        ! Take average of surrounding true cells if any
        if(nTrueNeighbors>0)then
           do iVar=1,nVar
              VarF_VG(iVar,iF:iF+1,jF:jF+1,kF:kF+1) = &
                (sum(VarC_VG(iVar,i-1:i+1,j,k),&
                MASK=IsTrueCellC_G(i-1:i+1,j,k))&
                +sum(VarC_VG(iVar,i,j-1:j+1,k),&
                MASK=IsTrueCellC_G(i,j-1:j+1,k))&
                +sum(VarC_VG(iVar,i,j,k-1:k+1),&
                MASK=IsTrueCellC_G(i,j,k-1:k+1))&
                )/nTrueNeighbors
           end do
        end if
        CYCLE
     end if
     
     ! Do first order prolongation for cells 
     ! whose parents were outside the body
     do iVar=1,nVar
        VarF_VG(iVar,iF:iF+1,jF:jF+1,kF:kF+1) = VarC_VG(iVar,i,j,k)
     end do
     if(nOrder==1 .or. nOrderIn==1) CYCLE
     
     !Do second order prolongation if possible
     
     ! First check to see if any of the points used to
     ! calculate the gradients are inside the body.  
     
     UseLimiter_D(1)= all(IsTrueCellC_G(i-1:i+1,j,k))
     UseLimiter_D(2)= all(IsTrueCellC_G(i,j-1:j+1,k))
     UseLimiter_D(3)= all(IsTrueCellC_G(i,j,k-1:k+1))
     
     ! Load the state variables into an array.
     do iVar=1,nVar
        VarC_DI(iVar,1,-1:1) =   VarC_VG(iVar,i-1:i+1,j,k)
        VarC_DI(iVar,2,-1:1) =   VarC_VG(iVar,i,j-1:j+1,k)
        VarC_DI(iVar,3,-1:1) =   VarC_VG(iVar,i,j,k-1:k+1)
     end do
     ! Compute the limiter and get the gradients
     do idir=1,3
        if (UseLimiter_D(idir))then
           ! MINMOD limiter
           do iVar=1,nVar
              SlopeRight = -VarC_DI(iVar,idir,0)+VarC_DI(iVar,idir,1)
              SlopeLeft  = +VarC_DI(iVar,idir,0)-VarC_DI(iVar,idir,-1)
              SlopeSign  = sign(cOne,SlopeRight)
              SlopeLim_D(iVar,idir)=SlopeSign*max(cZero,min(abs(SlopeRight),&
                   SlopeSign*SlopeLeft))
           end do
        else
           SlopeLim_D(:,idir)=cZero
        end if
     end do
     
     GradX_V = 0.25*SlopeLim_D(:,1)
     GradY_V = 0.25*SlopeLim_D(:,2)
     GradZ_V = 0.25*SlopeLim_D(:,3)
     
     VarF_VG(:,iF  ,jF  ,kF  ) = VarC_VG(:,i,j,k) - GradX_V - GradY_V - GradZ_V
     VarF_VG(:,iF+1,jF  ,kF  ) = VarC_VG(:,i,j,k) + GradX_V - GradY_V - GradZ_V
     VarF_VG(:,iF  ,jF+1,kF  ) = VarC_VG(:,i,j,k) - GradX_V + GradY_V - GradZ_V
     VarF_VG(:,iF+1,jF+1,kF  ) = VarC_VG(:,i,j,k) + GradX_V + GradY_V - GradZ_V
     VarF_VG(:,iF  ,jF  ,kF+1) = VarC_VG(:,i,j,k) - GradX_V - GradY_V + GradZ_V
     VarF_VG(:,iF+1,jF  ,kF+1) = VarC_VG(:,i,j,k) + GradX_V - GradY_V + GradZ_V
     VarF_VG(:,iF  ,jF+1,kF+1) = VarC_VG(:,i,j,k) - GradX_V + GradY_V + GradZ_V
     VarF_VG(:,iF+1,jF+1,kF+1) = VarC_VG(:,i,j,k) + GradX_V + GradY_V + GradZ_V
     
  end do
  end do
  end do

end subroutine prolong_block

!==============================================================================

subroutine set_refined_block_geometry
  use ModProcMH
  use ModGeometry, ONLY : dx_BLK,dy_BLK,dz_BLK,dxyz,xyzStart,xyzStart_BLK
  use ModNumConst
  use ModAMR, ONLY: local_cube,local_cubeBLK,iShiftChild_DI,nCell2_D 
  implicit none

  integer :: iPE,iBLK

  !\
  ! Loop over each of eight processors of the newly refined block
  ! and assign the geometry for the eight refined blocks.
  !/
  do iPE =  1, 8 

     if (iProc /= local_cube(iPE)) CYCLE
     !\
     ! Get the multiblock level for the current block.
     !/
     iBLK = local_cubeBLK(iPE)
     !\
     ! Set refined block corner offsets.
     !/
     XyzStart_BLK(:,iBLK) = XyzStart(:) &
          + (iShiftChild_DI(:,iPE)*nCell2_D - 0.25)*DXyz(:)
     !\
     ! Assign refined block geometry parameters.
     !/
     dx_BLK(iBLK) = dxyz(1)*cHalf
     dy_BLK(iBLK) = dxyz(2)*cHalf
     dz_BLK(iBLK) = dxyz(3)*cHalf

     call fix_block_geometry(iBLK)
  end do

end subroutine set_refined_block_geometry

!******************************************************************************



!*****************************COARSENING***************************************
!\
! Using each of the PEs in this group, create the necessary
! coarse solution block on the appropriate PE.
!/
subroutine create_coarse_soln_block(nPEsCrseBlk, PEsCrseBlk)

  use ModProcMH
  use ModMain
  use ModVarIndexes,ONLY:nVar
  use ModAdvance, ONLY : State_VGB, tmp1_BLK, &
       iTypeAdvance_B, SkippedBlock_
  use ModCT, ONLY : Bxface_BLK,Byface_BLK,Bzface_BLK        !^CFG IF CONSTRAINB
  use ModAMR, ONLY : local_cube,local_cubeBLK
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,R_BLK, &
       dx_BLK,dy_BLK,dz_BLK,dxyz,xyzStart_BLK,true_cell,true_BLK
  use ModGeometry, ONLY :  R2_BLK                           !^CFG IF SECONDBODY
  use ModNumConst
  use ModMpi
  use ModEnergy, ONLY: calc_energy_ghost

  implicit none

  integer, intent(in) :: nPEsCrseBlk, PEsCrseBlk(8)

  integer :: remaining_PE, remaining_BLK, icube, i, iBLK
  integer :: iError,iTag
  integer :: send_requests(7), receive_requests(7), &
       number_receive_requests, &
       Status_II(MPI_STATUS_SIZE, 7), Status_I(MPI_STATUS_SIZE)

  real    :: DtToReduce, DtBuff(8)
  integer :: MaxTypeAdvance, iTypeAdvance_I(8)

  logical :: IsTrueCellF_C(nI,nJ,nK)
  !--------------------------------------------------------------------------

  ! Return if processor is not needed here
  if(.not.any(PEsCrseBlk(1:nPEsCrseBlk) == iProc)) return

  !\
  ! Set block and processor numbers for remaining
  ! solution block (1 of 8).
  !/
  remaining_PE = local_cube(1)
  remaining_BLK = local_cubeBLK(1)

  ! set coarse block geometry
  if (iProc == remaining_PE) then 
     dxyz(1) = dx_BLK(remaining_BLK)
     dxyz(2) = dy_BLK(remaining_BLK)
     dxyz(3) = dz_BLK(remaining_BLK)
  end if

  !\
  ! Re-assign the geometry parameters for the eight 
  ! blocks involved in forming the coarse block.
  !/
  do icube = 1, 8
     iBLK = local_cubeBLK(icube)
     if (icube == 1 .and. &
          iProc == remaining_PE) then
        !\
        ! Set the geometry parameters for the newly 
        ! created coarse block.
        !/

        ! Save true_cell for this fine block !!!
        IsTrueCellF_C = true_cell(1:nI,1:nJ,1:nK,iBLK)

        call set_coarse_block_geometry(remaining_BLK)

        unusedBLK(remaining_BLK) = .false.

     else if (iProc == local_cube(icube)) then
        !\
        ! Reset the geometry parameters to the default values
        ! for the other seven blocks that are no longer used.
        !/
        Dx_BLK(iBLK)         = -777777.
        Dy_BLK(iBLK)         = -777777.
        Dz_BLK(iBLK)         = -777777.
        XyzStart_BLK(:,iBLK) = -777777.

        x_BLK(:,:,:,iBLK)    = -777777.
        y_BLK(:,:,:,iBLK)    = -777777.
        z_BLK(:,:,:,iBLK)    = -777777.
        R_BLK(:,:,:,iBLK)    = -777777.
        R2_BLK(:,:,:,iBLK)   = -777777.   !^CFG IF SECONDBODY

        UnusedBlk(iBLK)      = .true.
        iTypeAdvance_B(iBLK) = SkippedBlock_

     end if
  end do

  !\
  ! If inside the time loop, must also assign solution to
  ! newly created coarse block.
  !/
  if (time_loop) then
     if(nVar>=1)call assign_coarse_blk_soln
     if(UseConstrainB)then                        !^CFG IF CONSTRAINB BEGIN
        call assign_coarse_face_soln(Bxface_BLK,5)
        call assign_coarse_face_soln(Byface_BLK,6)
        call assign_coarse_face_soln(Bzface_BLK,7)
     end if                                       !^CFG END CONSTRAINB
 
     !\ 
     ! Calculating dt_BLK(Coarse)=2*min(dt_BLK(Ref))
     !/
     DtToReduce = minval(dt_BLK(local_cubeBLK(1:8)),&
          MASK=local_cube(1:8)==iProc)
     iTag=10
     if(iProc/=remaining_PE) then 
        call MPI_isend(DtToReduce,1,MPI_REAL,remaining_PE,&
             iTag,iComm,send_requests(1),iError)
        call MPI_waitall(1,send_requests,Status_II,iError)
     else
        DtBuff(1) = DtToReduce
        do i=2,nPEsCrseBlk
           call MPI_recv(DtBuff(i),1,MPI_REAL,PEsCrseBlk(i),&
                iTag,iComm,Status_I,iError)
        end do
        dt_BLK(remaining_BLK) = minval(DtBuff(1:nPEsCrseBlk))*2
     end if

     !\ 
     ! Calculating iTypeAdvance_B(Coarse)=maxval(iTypeAdvance_I(Refined))
     !/
     MaxTypeAdvance = maxval(iTypeAdvance_B(local_cubeBLK(1:8)),&
          MASK=local_cube(1:8)==iProc)

     if(iProc/=remaining_PE) then 
        call MPI_isend(MaxTypeAdvance,1,MPI_INTEGER,remaining_PE,&
             iTag,iComm,send_requests(1),iError)
        call MPI_waitall(1,send_requests,Status_II,iError)
     else
        iTypeAdvance_I(1) = MaxTypeAdvance
        do i=2,nPEsCrseBlk
           call MPI_recv(iTypeAdvance_I(i),1,MPI_INTEGER,PEsCrseBlk(i),&
                iTag,iComm,Status_I,iError)
        end do
        iTypeAdvance_B(remaining_BLK) = maxval(iTypeAdvance_I(1:nPEsCrseBlk))
     end if

     if (iProc == remaining_PE) then

        if(UseConstrainB)call Bface2Bcenter(remaining_BLK)  !^CFG IF CONSTRAINB

        call fix_soln_block(remaining_BLK)

        call calc_energy_ghost(remaining_BLK)

        call calc_other_soln_vars(remaining_BLK)
     end if

  end if

contains

  subroutine assign_coarse_blk_soln

    real, dimension(nVar,1:nI/2, 1:nJ/2, 1:nK/2, 8) ::&
         restricted_soln_blks

    integer,parameter :: isize=(nI/2)*(nJ/2)*(nK/2)*nVar
    integer :: number_send_requests, send_requests(7) 
    integer :: iBLK, iShift,jShift,kShift
    integer :: i, j, k, iF, jF, kF, iVar
    logical :: IsTrue_III(2,2,2)
    !-----------------------------------------------------------------------

    number_send_requests = 0

    do icube = 1, 8
       if (iProc /= local_cube(icube)) CYCLE
       iBLK = local_cubeBLK(icube)
       if(true_BLK(iBLK))then
          restricted_soln_blks(&
               1:nVar,1:nI/2,1:nJ/2,1:nK/2,icube) = 0.125*(&
               State_VGB(1:nVar,1:nI:2, 1:nJ:2, 1:nK:2,iBLK)+  &
               State_VGB(1:nVar,2:nI:2, 1:nJ:2, 1:nK:2,iBLK)+  &
               State_VGB(1:nVar,1:nI:2, 2:nJ:2, 1:nK:2,iBLK)+  &
               State_VGB(1:nVar,1:nI:2, 1:nJ:2, 2:nK:2,iBLK)+  &
               State_VGB(1:nVar,2:nI:2, 2:nJ:2, 1:nK:2,iBLK)+  &
               State_VGB(1:nVar,2:nI:2, 1:nJ:2, 2:nK:2,iBLK)+  &
               State_VGB(1:nVar,1:nI:2, 2:nJ:2, 2:nK:2,iBLK)+  &
               State_VGB(1:nVar,2:nI:2, 2:nJ:2, 2:nK:2,iBLK))
       else
          ! Exclude non-true (body) cells from the restriction 
          ! if there are any true cells
          do       k=1,nK/2; kF = 2*k-1 
             do    j=1,nJ/2; jF = 2*j-1
                do i=1,nI/2; iF = 2*i-1

                   ! Set true cell info for the fine cells to be restricted
                   if(iCube == 1)then
                      ! For the remaining fine block use the saved array
                      IsTrue_III = IsTrueCellF_C(iF:iF+1,jF:jF+1,kF:kF+1)
                   else
                      ! For the other fine blocks use the original array
                      IsTrue_III = true_cell(iF:iF+1,jF:jF+1,kF:kF+1,iBLK)
                   end if

                   if(all(IsTrue_III) .or. .not. any(IsTrue_III)) then
                      do iVar = 1, nVar
                         restricted_soln_blks(iVar,i,j,k,iCube) = 0.125 * &
                              sum(State_VGB(iVar,iF:iF+1,jF:jF+1,kF:kF+1,iBLK))
                      end do
                   else
                      do iVar = 1, nVar
                         restricted_soln_blks(iVar,i,j,k,iCube) = &
                              sum(State_VGB(iVar,iF:iF+1,jF:jF+1,kF:kF+1,iBLK)&
                              , MASK = IsTrue_III)/count(IsTrue_III)
                      end do
                   end if
                end do
             end do
          end do
       end if
       if (icube > 1 .and. iProc /= remaining_PE) then
          itag = iBLK
          number_send_requests = number_send_requests + 1
          call MPI_isend(restricted_soln_blks(1,1,1,1,icube), &
               isize, MPI_REAL, remaining_PE,itag, iComm,&
               send_requests(number_send_requests), iError)
       end if
    end do

    if (number_send_requests > 0) then
       call MPI_waitall(number_send_requests, &
            send_requests, Status_II, iError)
    end if

    number_receive_requests = 0

    if (iProc == remaining_PE) then ! remaining coarse block
       do icube = 2, 8
          if (local_cube(icube) /= remaining_PE) then
             itag = local_cubeBLK(icube)
             number_receive_requests = number_receive_requests + 1
             call MPI_irecv(restricted_soln_blks(1,1,1,1,icube), &
                  isize, MPI_REAL, local_cube(icube), itag, iComm,&
                  receive_requests(number_receive_requests), iError)
          end if
       end do

       if (number_receive_requests > 0) then
          call MPI_waitall(number_receive_requests, &
               receive_requests, Status_II, iError)
       end if
       !\
       ! Assign default value (to get corners).
       !/
       State_VGB(1:nVar,:,:,:,remaining_BLK) = cOne

       !\
       ! Assign coarse solution from the eight restricted solution quadrants of
       ! original fine blocks.
       !/
       do iCube = 1,8
          call get_shifts(iCube,iShift,jShift,kShift)
          State_VGB(1:nVar,          &
               1+iShift:nI/2+iShift, &
               1+jShift:nJ/2+jShift, &
               1+kShift:nK/2+kShift, &
               remaining_BLK)        &
               = restricted_soln_blks(:,1:nI/2,1:nJ/2,1:nK/2,icube)
       end do
       
    end if ! remaining coarse block

  end subroutine assign_coarse_blk_soln
end subroutine create_coarse_soln_block

!==============================================================================

subroutine set_coarse_block_geometry(iBLK)
  use ModNumConst
  use ModGeometry, ONLY : dx_BLK,dy_BLK,dz_BLK,dxyz,xyzStart_BLK
  use ModAMR, ONLY: iShiftChild_DI,nCell2_D
  implicit none

  integer, intent(in) :: iBLK

  dxyz = 2*dxyz

  ! Former the first child becomes the coarse block. Now his coordinates
  ! are in xyzStart_BLK(:,iBLK)
  XyzStart_BLK(:,iBLK)=XyzStart_BLK(:,iBLK)-&
       (iShiftChild_DI(:,1)*nCell2_D - 0.25)*DXyz(:)

  !\
  ! Set dx, dy, and dz for the new coarsened block.

  dx_BLK(iBLK) =dxyz(1)
  dy_BLK(iBLK) =dxyz(2)
  dz_BLK(iBLK) =dxyz(3)

  call fix_block_geometry(iBLK)

end subroutine set_coarse_block_geometry


!***************************************************************************
!*****************************COMMON ROUTINES*******************************
subroutine fix_soln_block(iBLK)
  use ModAdvance, ONLY : State_VGB, nVar
  use ModPhysics, ONLY : CellState_VI
  use ModGeometry, ONLY: body_BLK, true_cell
  use ModMain, ONLY: west_, TypeBC_I, nI, nJ, nK, gcn, body1_
  use ModParallel, ONLY: neiLwest, NOBLK
  implicit none

  integer, intent(in) :: iBLK
  integer::iVar,i,j,k

  if(body_BLK(iBLK)) then
     do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
        if(true_cell(i,j,k,iBLK))CYCLE
        State_VGB(1:nVar,i,j,k,iBLK) = CellState_VI(1:nVar,body1_)
     end do;end do; end do     
  end if

  ! For coupled (IH->GM) boundary condition fill in ghost cells
  ! with the first physical cell, because IH may not couple after AMR
  if(TypeBc_I(west_)=='coupled' .and. neiLwest(iBLK)==NOBLK)then
     State_VGB(:,nI+1,:,:,iBLK)   = State_VGB(:,nI,:,:,iBLK)
     State_VGB(:,nI+2,:,:,iBLK)   = State_VGB(:,nI,:,:,iBLK)
  endif

end subroutine fix_soln_block

!==============================================================================

subroutine calc_other_soln_vars(iBLK)
  use ModMain
  use ModAdvance, ONLY : fbody_x_BLK,fbody_y_BLK,fbody_z_BLK, &
       B0_DGB, B0ResChange_DXSB, B0ResChange_DYSB, B0ResChange_DZSB
  use ModConserveFlux, ONLY: init_cons_flux
  use ModImplicit, ONLY: UsePartImplicit             !^CFG IF IMPLICIT
  implicit none

  integer, intent(in) :: iBLK

  if(UseB0)then
     B0_DGB(:,:,:,:,iBLK) = 0.00
     B0ResChange_DXSB(:,:,:,:,iBLK) = 0.00
     B0ResChange_DYSB(:,:,:,:,iBLK) = 0.00
     B0ResChange_DZSB(:,:,:,:,iBLK) = 0.00

     call set_b0(iBLK)
  end if

  fbody_x_BLK(:,:,:,iBLK) = 0.00
  fbody_y_BLK(:,:,:,iBLK) = 0.00
  fbody_z_BLK(:,:,:,iBLK) = 0.00

  call init_cons_flux(iBLK)

  globalBLK = iBLK

  if(UseGravity.or.UseRotatingFrame) call body_force_averages

end subroutine calc_other_soln_vars
