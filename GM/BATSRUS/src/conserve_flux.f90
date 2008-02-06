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
  use ModGeometry,ONLY: UseCovariant                   
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
      
       CorrectedFlux_VXB(BnL_,   j,k,lFaceTo,globalBLK)= &
            LeftState_VX(Bx_,     lFaceFrom,j,k)*cQuarter
       CorrectedFlux_VXB(BnR_,   j,k,lFaceTo,globalBLK)= &
            RightState_VX(Bx_,     lFaceFrom,j,k)*cQuarter
    end do;end do
    if(UseCovariant)call save_bn_faceI(&   
         lFaceTo,lFaceFrom,globalBLK)      
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
       CorrectedFlux_VYB(BnL_,  i,k,lFaceTo,globalBLK)= &
            LeftState_VY(By_,i,lFaceFrom,k)*cQuarter
       CorrectedFlux_VYB(BnR_,  i,k,lFaceTo,globalBLK)= &
            RightState_VY(By_,i,lFaceFrom,k)*cQuarter
    end do; end do
    if(UseCovariant)call save_bn_faceJ(&  
         lFaceTo,lFaceFrom,globalBLK)     
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
       CorrectedFlux_VZB(BnL_,i,j,lFaceTo,globalBLK)= &
            LeftState_VZ(Bz_,i,j,lFaceFrom)*cQuarter
       CorrectedFlux_VZB(BnR_,i,j,lFaceTo,globalBLK)= &
            RightState_VZ(Bz_,i,j,lFaceFrom)*cQuarter
    end do; end do
    if(UseCovariant)call save_bn_faceK(& 
         lFaceTo,lFaceFrom,globalBLK)    
  end subroutine save_corrected_flux_z
end subroutine save_conservative_facefluxes

!==============================================================================

subroutine apply_cons_face_flux
  use ModProcMH
  use ModMain
  use ModAdvance
  use ModVarIndexes
  use ModGeometry, ONLY: body_BLK, true_cell
  use ModGeometry, ONLY: UseCovariant    
  use ModAMR
  use ModParallel, ONLY : &
       neiLtop, neiLbot, neiLeast, neiLwest, neiLnorth, neiLsouth, &
       neiPtop, neiPbot, neiPeast, neiPwest, neiPnorth, neiPsouth, &
       neiBtop, neiBbot, neiBeast, neiBwest, neiBnorth, neiBsouth
  use ModPhysics
  use ModNumConst
  implicit none

  integer :: i,j,k,iBlock
  logical :: oktest, oktest_me, oktest_row
  integer :: lFaceFrom, lFaceTo
  !--------------------------------------------------------------------------

  if(iProc==PROCtest .and. globalBLK==BLKtest)then
     call set_oktest('apply_cons_face_flux',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if

  iBlock = globalBLK

  if (neiLeast(iBlock)==-1)then
     if ( .not.unusedBlock_BP(neiBeast(1,iBlock),neiPeast(1,iBlock)).and.&
          .not.unusedBlock_BP(neiBeast(2,iBlock),neiPeast(2,iBlock)).and.&
          .not.unusedBlock_BP(neiBeast(3,iBlock),neiPeast(3,iBlock)).and.&
          .not.unusedBlock_BP(neiBeast(4,iBlock),neiPeast(4,iBlock)))then
        lFaceTo=1
        lFaceFrom=1
        call apply_corrected_flux_x
     end if
  end if

  if (neiLwest(iBlock)==-1) then
     if ( .not.unusedBlock_BP(neiBwest(1,iBlock),neiPwest(1,iBlock)).and.&
          .not.unusedBlock_BP(neiBwest(2,iBlock),neiPwest(2,iBlock)).and.&
          .not.unusedBlock_BP(neiBwest(3,iBlock),neiPwest(3,iBlock)).and.&
          .not.unusedBlock_BP(neiBwest(4,iBlock),neiPwest(4,iBlock)))then
        lFaceTo=nI+1
        lFaceFrom=2
        call apply_corrected_flux_x
     end if
  end if

  if (neiLsouth(iBlock)==-1) then
     if(.not.unusedBlock_BP(neiBsouth(1,iBlock),neiPsouth(1,iBlock)).and.&
          .not.unusedBlock_BP(neiBsouth(2,iBlock),neiPsouth(2,iBlock)).and.&
          .not.unusedBlock_BP(neiBsouth(3,iBlock),neiPsouth(3,iBlock)).and.&
          .not.unusedBlock_BP(neiBsouth(4,iBlock),neiPsouth(4,iBlock)))then
        lFaceTo=1
        lFaceFrom=1
        call apply_corrected_flux_y
     end if
  end if

  if (neiLnorth(iBlock)==-1) then
     if ( .not.unusedBlock_BP(neiBnorth(1,iBlock),neiPnorth(1,iBlock)).and.&
          .not.unusedBlock_BP(neiBnorth(2,iBlock),neiPnorth(2,iBlock)).and.&
          .not.unusedBlock_BP(neiBnorth(3,iBlock),neiPnorth(3,iBlock)).and.&
          .not.unusedBlock_BP(neiBnorth(4,iBlock),neiPnorth(4,iBlock)))then
        lFaceTo=nJ+1
        lFaceFrom=2
        call apply_corrected_flux_y
     end if
  end if

  if (neiLbot(iBlock)==-1) then
     if ( .not.unusedBlock_BP(neiBbot(1,iBlock),neiPbot(1,iBlock)).and. &
          .not.unusedBlock_BP(neiBbot(2,iBlock),neiPbot(2,iBlock)).and. &
          .not.unusedBlock_BP(neiBbot(3,iBlock),neiPbot(3,iBlock)).and. &
          .not.unusedBlock_BP(neiBbot(4,iBlock),neiPbot(4,iBlock)))then
        lFaceTo=1
        lFaceFrom=1
        call apply_corrected_flux_z
     end if
  end if

  if (neiLtop(iBlock)==-1) then
     if ( .not.unusedBlock_BP(neiBtop(1,iBlock),neiPtop(1,iBlock)).and. &
          .not.unusedBlock_BP(neiBtop(2,iBlock),neiPtop(2,iBlock)).and. &
          .not.unusedBlock_BP(neiBtop(3,iBlock),neiPtop(3,iBlock)).and. &
          .not.unusedBlock_BP(neiBtop(4,iBlock),neiPtop(4,iBlock))) then
        lFaceTo=nK+1
        lFaceFrom=2
        call apply_corrected_flux_z
     end if
  end if

contains

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
          if(UseCovariant)CYCLE    
          LeftState_VX(Bx_,lFaceTo,j,k) = &
               CorrectedFlux_VXB(BnL_,j,k,lFaceFrom,globalBLK)
          RightState_VX(Bx_,lFaceTo,j,k) = &
               CorrectedFlux_VXB(BnR_,j,k,lFaceFrom,globalBLK)
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
             if(UseCovariant)CYCLE       
             LeftState_VX(Bx_,lFaceTo,j,k) = &
                  CorrectedFlux_VXB(BnL_,j,k,lFaceFrom,globalBLK)
             RightState_VX(Bx_,lFaceTo,j,k) = &
                  CorrectedFlux_VXB(BnR_,j,k,lFaceFrom,globalBLK)
          end if
       end do; end do 
    end if
    if(UseCovariant)call apply_bn_faceI(&     
         lFaceFrom,lFaceTo,globalBLK)          
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
          if(UseCovariant)CYCLE          
          LeftState_VY(By_,i,lFaceTo,k) = &
               CorrectedFlux_VYB(BnL_,i,k,lFaceFrom,globalBLK)
          RightState_VY(By_,i,lFaceTo,k) = &
               CorrectedFlux_VYB(BnR_,i,k,lFaceFrom,globalBLK)
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
             if(UseCovariant)CYCLE    
             LeftState_VY(By_,i,lFaceTo,k) = &
                  CorrectedFlux_VYB(BnL_,i,k,lFaceFrom,globalBLK)
             RightState_VY(By_,i,lFaceTo,k) = &
                  CorrectedFlux_VYB(BnR_,i,k,lFaceFrom,globalBLK)
          end if
       end do; end do 
    end if
    if(UseCovariant)call apply_bn_faceJ(&      
         lFaceFrom,lFaceTo,globalBLK)          
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
          if(UseCovariant)CYCLE       
          LeftState_VZ(Bz_,i,j,lFaceTo) = &
               CorrectedFlux_VZB(BnL_,i,j,lFaceFrom,globalBLK)
          RightState_VZ(Bz_,i,j,lFaceTo) = &
               CorrectedFlux_VZB(BnR_,i,j,lFaceFrom,globalBLK)
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
             if(UseCovariant)CYCLE   
             LeftState_VZ(Bz_,i,j,lFaceTo) = &
                  CorrectedFlux_VZB(BnL_,i,j,lFaceFrom,globalBLK)
             RightState_VZ(Bz_,i,j,lFaceTo) = &
                  CorrectedFlux_VZB(BnR_,i,j,lFaceFrom,globalBLK)
          end if
       end do; end do 
    end if
    if(UseCovariant)call apply_bn_faceK(&      
         lFaceFrom,lFaceTo,globalBLK)           
  end subroutine apply_corrected_flux_z

end subroutine apply_cons_face_flux
