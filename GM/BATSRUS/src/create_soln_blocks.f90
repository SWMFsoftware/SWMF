!^CFG COPYRIGHT UM
subroutine get_shifts(iCube,iShift,jShift,kShift)
  use ModAMR, ONLY: iShiftChild_ID
  implicit none
  integer,intent(in) :: iCube
  integer,intent(out):: iShift, jShift, kShift

  iShift = iShiftChild_ID(iCube,1)
  jShift = iShiftChild_ID(iCube,2)
  kShift = iShiftChild_ID(iCube,3)

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
  use ModAdvance, ONLY :  State_VGB
  use ModAMR, ONLY : &
       local_cube,local_cubeBLK
  use ModParallel, ONLY: &
       neiLEV
  use ModGeometry, ONLY : &
       dx_BLK,dy_BLK,dz_BLK,dxyz,xyzStart,xyzStart_BLK,true_cell
  use ModMpi
  implicit none

  integer, intent(in) :: nPEsRefBlk, PEsRefBlk(9), refined_PE, refined_BLK

  integer :: icube, iBLK, isize, iError

  real    :: BxFaceFine_XQS(1:nJ,1:nK,4,2)  !^CFG IF CONSTRAINB BEGIN
  real    :: ByFaceFine_YQS(1:nI,1:nK,4,2)
  real    :: BzFaceFine_ZQS(1:nI,1:nJ,4,2)
  logical :: IsFinerNei_E(east_:top_)       !^CFG END CONSTRAINB

  real,dimension(1) :: DtBroadcast
  real, dimension(1:nVar,&
      1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn):: sol
  real, dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn):: sol1,sol2,sol3
  logical, dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn) :: INtrue_cell

  !---------------------------------------------------------------------------

  ! Return if processor is not needed here
  if(.not.any(PEsRefBlk(1:nPEsRefBlk)==iProc)) return

  !\
  ! Send coarse block geometry to all PEs owning
  ! eight newly created refined blocks.
  !/
  if (iProc == refined_PE) then
     dxyz(1) = dx_BLK(refined_BLK)
     dxyz(2) = dy_BLK(refined_BLK)
     dxyz(3) = dz_BLK(refined_BLK)
     xyzStart(:)=xyzStart_BLK(:,refined_BLK)
     DtBroadcast(1) = dt_BLK(refined_BLK)/2
  end if
  if (nPEsRefBlk > 1) then
     call cube_bcast_r(nPEsRefBlk,PEsRefBlk,xyzStart,3)
     call cube_bcast_r(nPEsRefBlk,PEsRefBlk,dxyz,3)
     call cube_bcast_r(nPEsRefBlk,PEsRefBlk,DtBroadcast,1)
  end if
  do iCube=1,8
     if(iProc==local_cube(iCube))&
          dt_BLK(local_cubeBLK(iCube))= DtBroadcast(1)
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

           globalBLK=local_cubeBLK(icube)

           if(UseConstrainB)call Bface2Bcenter !^CFG IF CONSTRAINB

           call fix_soln_block(globalBLK)

           call correctE

           call calc_other_soln_vars(globalBLK)
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
             send_requests(1), status(1,1), iError)
     end if
  else
     var = -888888.
     call MPI_recv(var,isize, MPI_REAL, Par(1), &
          itag, iComm,RECstatus,iError)
  end if

end subroutine cube_bcast_r


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
             send_requests(1), status(1,1), iError)
     end if
  else
     var = .false.
     call MPI_recv(var,isize, MPI_LOGICAL, Par(1), &
          itag, iComm, RECstatus, iError)
  end if

end subroutine cube_bcast_l


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

subroutine set_refined_block_geometry
  use ModProcMH
  use ModGeometry, ONLY : dx_BLK,dy_BLK,dz_BLK,dxyz,xyzStart,xyzStart_BLK
  use ModNumConst
  use ModAMR, ONLY: local_cube,local_cubeBLK,iShiftChild_ID 
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
     xyzStart_BLK(:,iBLK)=xyzStart(:) +(iShiftChild_ID(iPE,:) - cQuarter)*dxyz(:)
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
  use ModAdvance, ONLY : &
       State_VGB,tmp1_BLK
  use ModCT, ONLY : Bxface_BLK,Byface_BLK,Bzface_BLK          !^CFG IF CONSTRAINB
  use ModAMR, ONLY : local_cube,local_cubeBLK
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,R_BLK, &
       dx_BLK,dy_BLK,dz_BLK,dxyz,xyzStart_BLK
  use ModGeometry, ONLY :  R2_BLK                             !^CFG IF SECONDBODY
  use ModNumConst
  use ModMpi
  implicit none

  integer, intent(in) :: nPEsCrseBlk, PEsCrseBlk(8)

  integer :: remaining_PE, remaining_BLK, icube, i
  integer :: iError,iTag
  integer :: send_request, receive_requests(7), number_receive_requests, &
       status(MPI_STATUS_SIZE, 7)

  real:: DtToReduce, DtBuff(8)

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
     if (icube == 1 .and. &
          iProc == remaining_PE) then
        !\
        ! Set the geometry parameters for the newly 
        ! created coarse block.
        !/
        call set_coarse_block_geometry(remaining_BLK)

        unusedBLK(remaining_BLK) = .false.

     else if (iProc == local_cube(icube)) then
        !\
        ! Reset the geometry parameters for the other
        ! seven blocks that are nolonger used to defaults
        ! values.
        !/
        dx_BLK(local_cubeBLK(icube)) = -777777.
        dy_BLK(local_cubeBLK(icube)) = -777777.
        dz_BLK(local_cubeBLK(icube)) = -777777.
        xyzStart_BLK(:,local_cubeBLK(icube)) = -777777.

        x_BLK(:,:,:,local_cubeBLK(icube)) = -777777.
        y_BLK(:,:,:,local_cubeBLK(icube)) = -777777.
        z_BLK(:,:,:,local_cubeBLK(icube)) = -777777.
        R_BLK(:,:,:,local_cubeBLK(icube)) = -777777.
        R2_BLK(:,:,:,local_cubeBLK(icube)) = -777777.   !^CFG IF SECONDBODY

        unusedBLK(local_cubeBLK(icube)) = .true.

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
     !Calculating dt_BLK(Coarse)=2*min(dt_BLK(Ref))
     !/

     DtToReduce = minval(dt_BLK(local_cubeBLK(1:8)),&
          MASK=local_cube(1:8)==iProc)
     iTag=10
     if(iProc/=remaining_PE) then 
        call MPI_isend(DtToReduce,1,MPI_REAL,remaining_PE,&
             iTag,iComm,send_request,iError)
        call MPI_waitall(1,send_request,status(1,1),iError)
     else
        DtBuff(1) = DtToReduce
        do i=2,nPEsCrseBlk
           call MPI_recv(DtBuff(i),1,MPI_REAL,PEsCrseBlk(i),&
                iTag,iComm,status(1,1),iError)
        end do
        dt_BLK(remaining_BLK) = minval(DtBuff(1:nPEsCrseBlk))*2
     end if

     if (iProc == remaining_PE) then

        globalBLK=remaining_BLK

        if(UseConstrainB)call Bface2Bcenter  !^CFG IF CONSTRAINB

        call fix_soln_block(remaining_BLK)

        call correctE

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
    number_send_requests = 0

    do icube = 1, 8
       if (iProc == local_cube(icube)) then
          iBLK = local_cubeBLK(icube)
          restricted_soln_blks(&
               1:nVar,1:nI/2,1:nJ/2,1:nK/2,icube)=cEighth*(&
               State_VGB(1:nVar,1:nI:2, 1:nJ:2, 1:nK:2,iBLK)+  &
               State_VGB(1:nVar,2:nI:2, 1:nJ:2, 1:nK:2,iBLK)+  &
               State_VGB(1:nVar,1:nI:2, 2:nJ:2, 1:nK:2,iBLK)+  &
               State_VGB(1:nVar,1:nI:2, 1:nJ:2, 2:nK:2,iBLK)+  &
               State_VGB(1:nVar,2:nI:2, 2:nJ:2, 1:nK:2,iBLK)+  &
               State_VGB(1:nVar,2:nI:2, 1:nJ:2, 2:nK:2,iBLK)+  &
               State_VGB(1:nVar,1:nI:2, 2:nJ:2, 2:nK:2,iBLK)+  &
               State_VGB(1:nVar,2:nI:2, 2:nJ:2, 2:nK:2,iBLK))
          if (icube > 1 .and. iProc .ne. remaining_PE) then
             itag = iBLK
             number_send_requests = number_send_requests + 1
             call MPI_isend(restricted_soln_blks(1,1,1,1,icube), &
                  isize, MPI_REAL, remaining_PE,itag, iComm,&
                  send_requests(number_send_requests), iError)
          end if
       end if
    end do

    if (number_send_requests > 0) then
       call MPI_waitall(number_send_requests, &
            send_requests(1), &
            status(1,1), iError)
    end if

    number_receive_requests = 0

    if (iProc == remaining_PE) then ! remaining coarse block
       do icube = 2, 8
          if (local_cube(icube) .ne. remaining_PE) then
             itag = local_cubeBLK(icube)
             number_receive_requests = number_receive_requests + 1
             call MPI_irecv(restricted_soln_blks(1,1,1,1,icube), &
                  isize, MPI_REAL, local_cube(icube), itag, iComm,&
                  receive_requests(number_receive_requests), iError)
          end if
       end do

       if (number_receive_requests > 0) then
          call MPI_waitall(number_receive_requests, &
               receive_requests(1), &
               status(1,1), iError)
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
     State_VGB(1:nVar,&
          1+iShift:nI/2+iShift ,1+jShift:nJ/2+jShift,1+kShift:nK/2+kShift,&
           remaining_BLK) =  restricted_soln_blks(:,1:nI/2,1:nJ/2,1:nK/2,icube)
  end do
       
    end if ! remaining coarse block

  end subroutine assign_coarse_blk_soln
end subroutine create_coarse_soln_block

subroutine set_coarse_block_geometry(iBLK)
  use ModNumConst
  use ModGeometry, ONLY : dx_BLK,dy_BLK,dz_BLK,dxyz,xyzStart_BLK
  use ModAMR, ONLY: iShiftChild_ID
  implicit none

  integer, intent(in) :: iBLK

  dxyz=cTwo*dxyz
  ! Former the first child becomes the coarse block. Now his coordinates
  ! are in xyzStart_BLK(:,iBLK)
  xyzStart_BLK(:,iBLK)=xyzStart_BLK(:,iBLK)-&
       (iShiftChild_ID(1,:) - cQuarter)*dxyz(:)

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

subroutine calc_other_soln_vars(iBLK)
  use ModMain
  use ModAdvance, ONLY : fbody_x_BLK,fbody_y_BLK,fbody_z_BLK,qheat_BLK, &
       B0xCell_BLK,B0yCell_BLK,B0zCell_BLK, &
       B0xFace_x_BLK,B0yFace_x_BLK,B0zFace_x_BLK, &
       B0xFace_y_BLK,B0yFace_y_BLK,B0zFace_y_BLK, &
       B0xFace_z_BLK,B0yFace_z_BLK,B0zFace_z_BLK
  use ModImplicit,ONLY: UsePartImplicit             !^CFG IF IMPLICIT
  implicit none

  integer, intent(in) :: iBLK

  B0xCell_BLK(:,:,:,iBLK) = 0.00
  B0yCell_BLK(:,:,:,iBLK) = 0.00
  B0zCell_BLK(:,:,:,iBLK) = 0.00

  B0xFace_x_BLK(:,:,:,iBLK) = 0.00
  B0yFace_x_BLK(:,:,:,iBLK) = 0.00
  B0zFace_x_BLK(:,:,:,iBLK) = 0.00
  B0xFace_y_BLK(:,:,:,iBLK) = 0.00
  B0yFace_y_BLK(:,:,:,iBLK) = 0.00
  B0zFace_y_BLK(:,:,:,iBLK) = 0.00
  B0xFace_z_BLK(:,:,:,iBLK) = 0.00
  B0yFace_z_BLK(:,:,:,iBLK) = 0.00
  B0zFace_z_BLK(:,:,:,iBLK) = 0.00


  call set_b0(iBLK)

  fbody_x_BLK(:,:,:,iBLK) = 0.00
  fbody_y_BLK(:,:,:,iBLK) = 0.00
  fbody_z_BLK(:,:,:,iBLK) = 0.00

  if(UsePartImplicit)&                              !^CFG IF IMPLICIT
       call init_conservative_facefluxes(iBLK)      !^CFG IF IMPLICIT  

  globalBLK = iBLK

  if(UseGravity.or.UseRotatingFrame) call body_force_averages

  qheat_BLK(:,:,:,iBLK) = 0.00

  if(UseUserHeating)     call user_heat_source      !^CFG IF USERFILES


end subroutine calc_other_soln_vars
