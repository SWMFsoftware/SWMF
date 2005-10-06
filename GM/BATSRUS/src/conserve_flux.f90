!^CFG COPYRIGHT UM

subroutine init_conservative_facefluxes(iBLK)
  use ModAdvance, ONLY: CorrectedFlux_VXB, &
       CorrectedFlux_VYB, CorrectedFlux_VZB
  use ModNumConst,ONLY:cZero
  implicit none
  integer,intent(in)::iBLK
  !--------------------------------------------------------------------------
  CorrectedFlux_VXB(:,:,:,:,iBLK)=cZero
  CorrectedFlux_VYB(:,:,:,:,iBLK)=cZero
  CorrectedFlux_VZB(:,:,:,:,iBLK)=cZero
end subroutine init_conservative_facefluxes

!============================================================================

! Save conservative flux.
subroutine save_conservative_facefluxes
  use ModProcMH
  use ModMain
  use ModNumConst
  use ModAdvance
  use ModGeometry,ONLY: UseCovariant                   !^CFG IF COVARIANT
  use ModParallel, ONLY : &
       neiLtop, neiLbot, neiLeast, neiLwest, neiLnorth, neiLsouth, &
       neiPtop, neiPbot, neiPeast, neiPwest, neiPnorth, neiPsouth, &
       neiBtop, neiBbot, neiBeast, neiBwest, neiBnorth, neiBsouth
  implicit none

  logical :: oktest, oktest_me, oktest_row
  integer :: lFaceFrom,lFaceTo,i,j,k
  !--------------------------------------------------------------------------

  if(iProc==PROCtest .and. globalBLK==BLKtest)then
     call set_oktest('save_conserv_flux',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if

  if (neiLeast(globalBLK)==+1) then
     lFaceFrom=1
     lFaceTo=1
     call save_corrected_flux_x
  end if

  if (neiLwest(globalBLK)==+1) then
     lFaceFrom=1+nI
     lFaceTo=2
     call save_corrected_flux_x
  end if

  if (neiLsouth(globalBLK)==+1) then
     lFaceFrom=1
     lFaceTo=1
     call save_corrected_flux_y
  end if

  if (neiLnorth(globalBLK)==+1) then
     lFaceFrom=1+nJ
     lFaceTo=2
     call save_corrected_flux_y
  end if

  if (neiLbot(globalBLK)==+1) then
     lFaceFrom=1
     lFaceTo=1
     call save_corrected_flux_z
  end if

  if (neiLtop(globalBLK)==+1) then
     lFaceFrom=nK+1
     lFaceTo=2
     call save_corrected_flux_z
  end if

Contains

  !===========================================================================

  subroutine save_corrected_flux_x
    do k=1,nK;do j=1,nJ
       CorrectedFlux_VXB(1:nVar-1,j,k,lFaceTo,globalBLK)&
            =Flux_VX(1:nVar-1, lFaceFrom,j,k)
       CorrectedFlux_VXB(nVar,j,k,lFaceTo,globalBLK)&
            =Flux_VX(Energy_, lFaceFrom,j,k)

       CorrectedFlux_VXB(Vdt_,   j,k,lFaceTo,globalBLK)= &
            VdtFace_x(             lFaceFrom,j,k)
       !^CFG IF NOT COVARIANT BEGIN
       CorrectedFlux_VXB(BnL_,   j,k,lFaceTo,globalBLK)= &
            LeftState_VX(Bx_,     lFaceFrom,j,k)*cQuarter
       CorrectedFlux_VXB(BnR_,   j,k,lFaceTo,globalBLK)= &
            RightState_VX(Bx_,     lFaceFrom,j,k)*cQuarter
       !^CFG END COVARIANT
    end do;end do
    if(UseCovariant)call save_bn_faceI(&   !^CFG IF COVARIANT
         lFaceTo,lFaceFrom,globalBLK)      !^CFG IF COVARIANT
  end subroutine save_corrected_flux_x

  !===========================================================================

  subroutine save_corrected_flux_y
    do k=1,nK;do i=1,nI
       CorrectedFlux_VYB(1:nVar-1,  i,k,lFaceTo,globalBLK)= &
            Flux_VY(1:nVar-1,i,lFaceFrom,k)
       CorrectedFlux_VYB(nVar,  i,k,lFaceTo,globalBLK)= &
            Flux_VY(Energy_,i,lFaceFrom,k)
       CorrectedFlux_VYB(Vdt_,     i,k,lFaceTo,globalBLK)= &
            VdtFace_y(i,lFaceFrom,k)
       !^CFG IF NOT COVARIANT BEGIN
       CorrectedFlux_VYB(BnL_,  i,k,lFaceTo,globalBLK)= &
            LeftState_VY(By_,i,lFaceFrom,k)*cQuarter
       CorrectedFlux_VYB(BnR_,  i,k,lFaceTo,globalBLK)= &
            RightState_VY(By_,i,lFaceFrom,k)*cQuarter
       !^CFG END COVARIANT
    end do; end do
    if(UseCovariant)call save_bn_faceJ(&  !^CFG IF COVARIANT
         lFaceTo,lFaceFrom,globalBLK)     !^CFG IF COVARIANT
  end subroutine save_corrected_flux_y

  !===========================================================================

  subroutine save_corrected_flux_z
    do j=1,nJ;do i=1,nI
       CorrectedFlux_VZB(1:nVar-1,  i,j,lFaceTo,globalBLK)= &
            Flux_VZ(1:nVar-1,i,j,lFaceFrom)
       CorrectedFlux_VZB(nVar,  i,j,lFaceTo,globalBLK)= &
            Flux_VZ(Energy_,i,j,lFaceFrom)

       CorrectedFlux_VZB(Vdt_,  i,j,lFaceTo,globalBLK)= &
            VdtFace_z(i,j,lFaceFrom)
       !^CFG IF NOT COVARIANT BEGIN
       CorrectedFlux_VZB(BnL_,i,j,lFaceTo,globalBLK)= &
            LeftState_VZ(Bz_,i,j,lFaceFrom)*cQuarter
       CorrectedFlux_VZB(BnR_,i,j,lFaceTo,globalBLK)= &
            RightState_VZ(Bz_,i,j,lFaceFrom)*cQuarter
       !^CFG END COVARIANT
    end do; end do
    if(UseCovariant)call save_bn_faceK(& !^CFG IF COVARIANT
         lFaceTo,lFaceFrom,globalBLK)    !^CFG IF COVARIANT
  end subroutine save_corrected_flux_z
end subroutine save_conservative_facefluxes

!==============================================================================

subroutine apply_cons_face_flux
  use ModProcMH
  use ModMain
  use ModAdvance
  use ModVarIndexes
  use ModGeometry, ONLY: body_BLK, true_cell
  use ModGeometry, ONLY: UseCovariant    !^CFG IF COVARIANT
  use ModAMR
  use ModParallel, ONLY : &
       neiLtop, neiLbot, neiLeast, neiLwest, neiLnorth, neiLsouth, &
       neiPtop, neiPbot, neiPeast, neiPwest, neiPnorth, neiPsouth, &
       neiBtop, neiBbot, neiBeast, neiBwest, neiBnorth, neiBsouth
  use ModPhysics
  use ModNumConst
  implicit none

  integer :: i,j,k
  logical :: oktest, oktest_me, oktest_row
  integer :: lFaceFrom, lFaceTo
  !--------------------------------------------------------------------------

  if(iProc==PROCtest .and. globalBLK==BLKtest)then
     call set_oktest('apply_cons_face_flux',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if

  if (neiLeast(globalBLK)==-1)then
     !^CFG IF IMPLICIT BEGIN
     if ( .not.unusedBlock_BP(neiBeast(1,globalBLK),neiPeast(1,globalBLK)).and.&
          .not.unusedBlock_BP(neiBeast(2,globalBLK),neiPeast(2,globalBLK)).and.&
          .not.unusedBlock_BP(neiBeast(3,globalBLK),neiPeast(3,globalBLK)).and.&
          .not.unusedBlock_BP(neiBeast(4,globalBLK),neiPeast(4,globalBLK)))then
        !^CFG END IMPLICIT
        lFaceTo=1
        lFaceFrom=1
        call apply_corrected_flux_x
     end if                                 !^CFG IF IMPLICIT
  end if

  if (neiLwest(globalBLK)==-1) then
     !^CFG IF IMPLICIT BEGIN
     if ( .not.unusedBlock_BP(neiBwest(1,globalBLK),neiPwest(1,globalBLK)).and.&
          .not.unusedBlock_BP(neiBwest(2,globalBLK),neiPwest(2,globalBLK)).and.&
          .not.unusedBlock_BP(neiBwest(3,globalBLK),neiPwest(3,globalBLK)).and.&
          .not.unusedBlock_BP(neiBwest(4,globalBLK),neiPwest(4,globalBLK)))then
        !^CFG END IMPLICIT
        lFaceTo=nI+1
        lFaceFrom=2
        call apply_corrected_flux_x
     end if                                                     !^CFG IF IMPLICIT
  end if

  if (neiLsouth(globalBLK)==-1) then
     !^CFG IF IMPLICIT BEGIN
     if ( .not.unusedBlock_BP(neiBsouth(1,globalBLK),neiPsouth(1,globalBLK)).and.&
          .not.unusedBlock_BP(neiBsouth(2,globalBLK),neiPsouth(2,globalBLK)).and.&
          .not.unusedBlock_BP(neiBsouth(3,globalBLK),neiPsouth(3,globalBLK)).and.&
          .not.unusedBlock_BP(neiBsouth(4,globalBLK),neiPsouth(4,globalBLK)))then
        !^CFG END IMPLICIT
        lFaceTo=1
        lFaceFrom=1
        call apply_corrected_flux_y
     end if                                         !^CFG IF IMPLICIT
  end if

  if (neiLnorth(globalBLK)==-1) then
     !^CFG IF IMPLICIT BEGIN
     if ( .not.unusedBlock_BP(neiBnorth(1,globalBLK),neiPnorth(1,globalBLK)).and.&
          .not.unusedBlock_BP(neiBnorth(2,globalBLK),neiPnorth(2,globalBLK)).and.&
          .not.unusedBlock_BP(neiBnorth(3,globalBLK),neiPnorth(3,globalBLK)).and.&
          .not.unusedBlock_BP(neiBnorth(4,globalBLK),neiPnorth(4,globalBLK)))then
        !^CFG END IMPLICIT
        lFaceTo=nJ+1
        lFaceFrom=2
        call apply_corrected_flux_y
     end if                                          !^CFG IF IMPLICIT
  end if

  if (neiLbot(globalBLK)==-1) then
     !^CFG IF IMPLICIT BEGIN
     if ( .not.unusedBlock_BP(neiBbot(1,globalBLK),neiPbot(1,globalBLK)).and. &
          .not.unusedBlock_BP(neiBbot(2,globalBLK),neiPbot(2,globalBLK)).and. &
          .not.unusedBlock_BP(neiBbot(3,globalBLK),neiPbot(3,globalBLK)).and. &
          .not.unusedBlock_BP(neiBbot(4,globalBLK),neiPbot(4,globalBLK)))then
        !^CFG END IMPLICIT
        lFaceTo=1
        lFaceFrom=1
        call apply_corrected_flux_z
     end if                                            !^CFG IF IMPLICIT
  end if

  if (neiLtop(globalBLK)==-1) then
     !^CFG IF IMPLICIT BEGIN
     if ( .not.unusedBlock_BP(neiBtop(1,globalBLK),neiPtop(1,globalBLK)).and. &
          .not.unusedBlock_BP(neiBtop(2,globalBLK),neiPtop(2,globalBLK)).and. &
          .not.unusedBlock_BP(neiBtop(3,globalBLK),neiPtop(3,globalBLK)).and. &
          .not.unusedBlock_BP(neiBtop(4,globalBLK),neiPtop(4,globalBLK))) then
        !^CFG END IMPLICIT
        lFaceTo=nK+1
        lFaceFrom=2
        call apply_corrected_flux_z
     end if                                   !^CFG IF IMPLICIT
  end if

Contains

  !============================================================================

  subroutine apply_corrected_flux_x
    if (.not. body_BLK(globalBLK)) then
       do k = 1, nK; do j = 1, nJ
          Flux_VX(1:nVar-1,lFaceTo,j,k) = &
               CorrectedFlux_VXB(1:nVar-1,j,k,lFaceFrom,globalBLK)
          Flux_VX(Energy_,lFaceTo,j,k) = &
               CorrectedFlux_VXB(nVar,j,k,lFaceFrom,globalBLK)
          VdtFace_x(lFaceTo,j,k) = &
               CorrectedFlux_VXB(Vdt_,j,k,lFaceFrom,globalBLK)
          if(UseCovariant)CYCLE    !^CFG IF COVARIANT 
          !^CFG IF NOT COVARIANT BEGIN
          LeftState_VX(Bx_,lFaceTo,j,k) = &
               CorrectedFlux_VXB(BnL_,j,k,lFaceFrom,globalBLK)
          RightState_VX(Bx_,lFaceTo,j,k) = &
               CorrectedFlux_VXB(BnR_,j,k,lFaceFrom,globalBLK)
          !^CFG END COVARIANT
       end do;end do
    else
       do k = 1, nK; do j = 1, nJ
          if (all(true_cell(lFaceTo-1:lFaceTo,j,k,globalBLK))) then
             Flux_VX(1:nVar-1,lFaceTo,j,k) = &
                  CorrectedFlux_VXB(1:nVar-1,j,k,lFaceFrom,globalBLK)
             Flux_VX(Energy_,lFaceTo,j,k) = &
                  CorrectedFlux_VXB(nVar,j,k,lFaceFrom,globalBLK)
             VdtFace_x(lFaceTo,j,k) = &
                  CorrectedFlux_VXB(Vdt_,j,k,lFaceFrom,globalBLK)
             if(UseCovariant)CYCLE       !^CFG IF COVARIANT 
             !^CFG IF NOT COVARIANT BEGIN
             LeftState_VX(Bx_,lFaceTo,j,k) = &
                  CorrectedFlux_VXB(BnL_,j,k,lFaceFrom,globalBLK)
             RightState_VX(Bx_,lFaceTo,j,k) = &
                  CorrectedFlux_VXB(BnR_,j,k,lFaceFrom,globalBLK)
             !^CFG END COVARIANT
          end if
       end do; end do 
    end if
    if(UseCovariant)call apply_bn_faceI(&      !^CFG IF COVARIANT
         lFaceFrom,lFaceTo,globalBLK)          !^CFG IF COVARIANT
  end subroutine apply_corrected_flux_x

  !============================================================================

  subroutine apply_corrected_flux_y
    if (.not. body_BLK(globalBLK)) then
       do k = 1, nK; do i = 1, nI
          Flux_VY(1:nVar-1,i,lFaceTo,k) = &
               CorrectedFlux_VYB(1:nVar-1,i,k,lFaceFrom,globalBLK)
          Flux_VY(Energy_,i,lFaceTo,k) = &
               CorrectedFlux_VYB(nVar,i,k,lFaceFrom,globalBLK)
          VdtFace_y(i,lFaceTo,k)= &
               CorrectedFlux_VYB(Vdt_,i,k,lFaceFrom,globalBLK)
          if(UseCovariant)CYCLE          !^CFG IF COVARIANT 
          !^CFG IF NOT COVARIANT BEGIN
          LeftState_VY(By_,i,lFaceTo,k) = &
               CorrectedFlux_VYB(BnL_,i,k,lFaceFrom,globalBLK)
          RightState_VY(By_,i,lFaceTo,k) = &
               CorrectedFlux_VYB(BnR_,i,k,lFaceFrom,globalBLK)
          !^CFG END COVARIANT
       end do; end do 
    else
       do k = 1, nK; do i = 1, nI
          if (all(true_cell(i,lFaceTo-1:lFaceTo,k,globalBLK))) then
             Flux_VY(1:nVar-1,i,lFaceTo,k) = &
                  CorrectedFlux_VYB(1:nVar-1,i,k,lFaceFrom,globalBLK)
             Flux_VY(Energy_,i,lFaceTo,k) = &
                  CorrectedFlux_VYB(nVar,i,k,lFaceFrom,globalBLK)
             VdtFace_y(i,lFaceTo,k)= &
                  CorrectedFlux_VYB(Vdt_,i,k,lFaceFrom,globalBLK)
             if(UseCovariant)CYCLE     !^CFG IF COVARIANT 
             !^CFG IF NOT COVARIANT BEGIN
             LeftState_VY(By_,i,lFaceTo,k) = &
                  CorrectedFlux_VYB(BnL_,i,k,lFaceFrom,globalBLK)
             RightState_VY(By_,i,lFaceTo,k) = &
                  CorrectedFlux_VYB(BnR_,i,k,lFaceFrom,globalBLK)
             !^CFG END COVARIANT
          end if
       end do; end do 
    end if
    if(UseCovariant)call apply_bn_faceJ(&      !^CFG IF COVARIANT
         lFaceFrom,lFaceTo,globalBLK)          !^CFG IF COVARIANT
  end subroutine apply_corrected_flux_y

  !============================================================================

  subroutine apply_corrected_flux_z
    if (.not. body_BLK(globalBLK)) then
       do j = 1, nJ; do i = 1, nI
          Flux_VZ(1:nVar-1,i,j,lFaceTo) = &
               CorrectedFlux_VZB(1:nVar-1,i,j,lFaceFrom,globalBLK)
          Flux_VZ(Energy_,i,j,lFaceTo) = &
               CorrectedFlux_VZB(nVar,i,j,lFaceFrom,globalBLK)
          VdtFace_z(i,j,lFaceTo) = &
               CorrectedFlux_VZB(Vdt_,i,j,lFaceFrom,globalBLK)
          if(UseCovariant)CYCLE      !^CFG IF COVARIANT 
          !^CFG IF NOT COVARIANT BEGIN
          LeftState_VZ(Bz_,i,j,lFaceTo) = &
               CorrectedFlux_VZB(BnL_,i,j,lFaceFrom,globalBLK)
          RightState_VZ(Bz_,i,j,lFaceTo) = &
               CorrectedFlux_VZB(BnR_,i,j,lFaceFrom,globalBLK)
          !^CFG END COVARIANT
       end do; end do 
    else
       do j = 1, nJ; do i = 1, nI
          if (all(true_cell(i,j,lFaceTo-1:lFaceTo,globalBLK))) then
             Flux_VZ(1:nVar-1,i,j,lFaceTo) = &
                  CorrectedFlux_VZB(1:nVar-1,i,j,lFaceFrom,globalBLK)
             Flux_VZ(Energy_,i,j,lFaceTo) = &
                  CorrectedFlux_VZB(nVar,i,j,lFaceFrom,globalBLK)
             VdtFace_z(i,j,lFaceTo) = &
                  CorrectedFlux_VZB(Vdt_,i,j,lFaceFrom,globalBLK)
             if(UseCovariant)CYCLE    !^CFG IF COVARIANT 
             !^CFG IF NOT COVARIANT BEGIN
             LeftState_VZ(Bz_,i,j,lFaceTo) = &
                  CorrectedFlux_VZB(BnL_,i,j,lFaceFrom,globalBLK)
             RightState_VZ(Bz_,i,j,lFaceTo) = &
                  CorrectedFlux_VZB(BnR_,i,j,lFaceFrom,globalBLK)
             !^CFG END COVARIANT
          end if
       end do; end do 
    end if
    if(UseCovariant)call apply_bn_faceK(&       !^CFG IF COVARIANT
         lFaceFrom,lFaceTo,globalBLK)           !^CFG IF COVARIANT
  end subroutine apply_corrected_flux_z

end subroutine apply_cons_face_flux
