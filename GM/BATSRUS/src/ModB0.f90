module ModB0

  ! The magnetic field can be split into an analytic and numeric part.
  ! The analytic part is B0. It should satisfy div(B0) = 0, 
  ! and usually curl(B0) = 0 except if UseCurlB0 = .true.
  ! Even if curl(B0) = 0 analytically, the numerical representation
  ! may have a finite curl.
  ! It is often constant in time dB0/dt = 0, but not necessarily.
  ! This module provides data and methods for B0.
  
  use ModSize

  use ModMain, ONLY: UseB0, UseB0Source, UseCurlB0 !!! should be here

  implicit none
  SAVE

  !!! private ! except

  ! Cell-centered B0 field vector
  real, public, allocatable:: B0_DGB(:,:,:,:,:)

  ! Face-centered B0 field arrays for one block
  real, public, allocatable:: B0_DX(:,:,:,:), B0_DY(:,:,:,:), B0_DZ(:,:,:,:)

  ! The numerical curl and divergence of B0 for one block
  real, public, allocatable :: CurlB0_DC(:,:,:,:)
  real, public, allocatable :: DivB0_C(:,:,:)


  ! Local variables

  ! B0 field at the resolution change interfaces
  real, allocatable :: B0ResChange_DXSB(:,:,:,:,:)
  real, allocatable :: B0ResChange_DYSB(:,:,:,:,:)
  real, allocatable :: B0ResChange_DZSB(:,:,:,:,:)

contains
  !===========================================================================
  subroutine init_mod_b0

    use ModMain, ONLY: UseConstrainB
    !------------------------------------------------------------------------

    if(.not.allocated(B0_DGB))then
       allocate( &
            B0_DGB(nDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), &
            B0ResChange_DXSB(nDim,nJ,nK,1:2,MaxBlock),           &
            B0ResChange_DYSB(nDim,nI,nK,3:4,MaxBlock),           &
            B0ResChange_DZSB(nDim,nI,nJ,5:6,MaxBlock))
       B0_DGB           = 0.0
       B0ResChange_DXSB = 0.0
       B0ResChange_DYSB = 0.0
       B0ResChange_DZSB = 0.0

       if(UseConstrainB)then
          ! The current implementation of CT requires fluxes
          ! between ghost cells. Not a good solution.
          allocate( &
               B0_DX(MaxDim,nI+1,0:nJ+1,0:nK+1), &
               B0_DY(MaxDim,0:nI+1,nJ+1,0:nK+1), &
               B0_DZ(MaxDim,0:nI+1,0:nJ+1,nK+1))
       else
          allocate( &
               B0_DX(MaxDim,nI+1,nJ,nK), &
               B0_DY(MaxDim,nI,nJ+1,nK), &
               B0_DZ(MaxDim,nI,nJ,nK+1))
       end if
    end if

    if(UseB0Source .and. .not.allocated(DivB0_C)) &
         allocate(DivB0_C(nI,nJ,nK))

    if((UseCurlB0 .or. UseB0Source) .and. .not.allocated(CurlB0_DC)) &
         allocate(CurlB0_DC(3,nI,nJ,nK))

  end subroutine init_mod_b0
  !===========================================================================
  subroutine clean_mod_b0

    if(allocated(B0_DGB))    deallocate(B0_DGB, &
         B0ResChange_DXSB, B0ResChange_DYSB, B0ResChange_DZSB)

    if(allocated(DivB0_C))   deallocate(DivB0_C)
    if(allocated(CurlB0_DC)) deallocate(CurlB0_DC)

  end subroutine clean_mod_b0
  !===========================================================================
  subroutine set_b0_face(iBlock)

    ! Calculate the face centered B0 for block iBlock

    use ModParallel, ONLY : NeiLev
    integer,intent(in)::iBlock
    !-------------------------------------------------------------------------
    if(.not.UseB0) RETURN

    ! Face averages (X)
    B0_DX = 0.5*(B0_DGB(:,0:nI  ,1:nJ,1:nK,iBlock) &
         +       B0_DGB(:,1:nI+1,1:nJ,1:nK,iBlock))

    ! Face averages (Y)
    B0_DY = 0.5*(B0_DGB(:,1:nI,0:nJ  ,1:nK,iBlock) &
         +       B0_DGB(:,1:nI,1:nJ+1,1:nK,iBlock))

    ! Face averages (Z)
    B0_DZ = 0.5*(B0_DGB(:,1:nI,1:nJ,0:nK  ,iBlock) &
         +       B0_DGB(:,1:nI,1:nJ,1:nK+1,iBlock))

    ! Correct B0 at resolution changes
    if(NeiLev(1,iBlock) == -1) &
         B0_DX(:,1   ,1:nJ,1:nK) = B0ResChange_DXSB(:,:,:,1,iBlock)
    if(NeiLev(2,iBlock) == -1) &
         B0_DX(:,nI+1,1:nJ,1:nK) = B0ResChange_DXSB(:,:,:,2,iBlock)
    if(NeiLev(3,iBlock) == -1) &
         B0_DY(:,1:nI,1   ,1:nK) = B0ResChange_DYSB(:,:,:,3,iBlock)
    if(NeiLev(4,iBlock) == -1) &
         B0_DY(:,1:nI,1+nJ,1:nK) = B0ResChange_DYSB(:,:,:,4,iBlock)
    if(NeiLev(5,iBlock) == -1) &
         B0_DZ(:,1:nI,1:nJ,1   ) = B0ResChange_DZSB(:,:,:,5,iBlock)
    if(NeiLev(6,iBlock) == -1) &
         B0_DZ(:,1:nI,1:nJ,1+nK) = B0ResChange_DZSB(:,:,:,6,iBlock)

  end subroutine set_b0_face
  !===========================================================================
  subroutine set_b0_reschange

    ! Set face centered B0 at resolution changes. Use the face area weighted
    ! average of the fine side B0 for the coarce face. This works for 
    ! non-Cartesian grids too, because all the fine and coarse face normal
    ! vectors are parallel with each other (see algorithm in BATL_grid).

    use BATL_lib, ONLY: nDim, nBlock, Unused_B, DiLevelNei_IIIB, &
         IsCartesian, CellFace_DB, CellFace_DFB, message_pass_face

    integer:: i, j, k, iBlock
    real:: Coef

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'set_b0_reschange'
    !------------------------------------------------------------------------
    ! There are no resolution changes in 1D
    if(nDim == 1) RETURN

    if(.not.UseB0) RETURN
    
    call set_oktest(NameSub, DoTest, DoTestMe)
    if(DoTest)write(*,*)NameSub,' starting'


    ! For Cartesian grid take 1/8-th of the contributing fine B0 values.
    ! For non-Cartesian grid the averaged cell values are weighted by face area
    if(IsCartesian) Coef = 0.125

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       if(DiLevelNei_IIIB(-1,0,0,iBlock) == +1) then
          do k = 1, nK; do j = 1, nJ
             if(.not.IsCartesian) Coef = 0.5*CellFace_DFB(1,1,j,k,iBlock)
             B0ResChange_DXSB(:,j,k,1,iBlock) =     &
                  Coef*( B0_DGB(:,0,j,k,iBlock) &
                  +      B0_DGB(:,1,j,k,iBlock) )
          end do; end do
       end if

       if(DiLevelNei_IIIB(+1,0,0,iBlock) == +1) then
          do k = 1, nK; do j = 1, nJ
             if(.not.IsCartesian) Coef = 0.5*CellFace_DFB(1,nI+1,j,k,iBlock)
             B0ResChange_DXSB(:,j,k,2,iBlock) =    &
                  Coef*( B0_DGB(:,nI  ,j,k,iBlock) &
                  +      B0_DGB(:,nI+1,j,k,iBlock) )
          end do; end do
       end if

       if(DiLevelNei_IIIB(0,-1,0,iBlock) == +1) then
          do k = 1, nK; do i = 1, nI
             if(.not.IsCartesian) Coef = 0.5*CellFace_DFB(2,i,1,k,iBlock)
             B0ResChange_DYSB(:,i,k,3,iBlock) = &
                  Coef*( B0_DGB(:,i,0,k,iBlock) &
                  +      B0_DGB(:,i,1,k,iBlock) )
          end do; end do
       end if

       if(DiLevelNei_IIIB(0,+1,0,iBlock) == +1) then
          do k = 1, nK; do i = 1, nI
             if(.not.IsCartesian) Coef = 0.5*CellFace_DFB(2,i,nJ+1,k,iBlock)
             B0ResChange_DYSB(:,i,k,4,iBlock) =    &
                  Coef*( B0_DGB(:,i,nJ  ,k,iBlock) &
                  +      B0_DGB(:,i,nJ+1,k,iBlock) )
          end do; end do
       end if

       if(nDim < 3) CYCLE

       if(DiLevelNei_IIIB(0,0,-1,iBlock) == +1) then
          do j = 1, nJ; do i = 1, nI
             if(.not.IsCartesian) Coef = 0.5*CellFace_DFB(3,i,j,1,iBlock)
             B0ResChange_DZSB(:,i,j,5,iBlock) = &
                  Coef*( B0_DGB(:,i,j,0,iBlock) &
                  +      B0_DGB(:,i,j,1,iBlock) )
          end do; end do
       end if

       if(DiLevelNei_IIIB(0,0,+1,iBlock) == +1)  then
          do j = 1, nJ; do i = 1, nI
             if(.not.IsCartesian) Coef = 0.5*CellFace_DFB(3,i,j,nK+1,iBlock)
             B0ResChange_DZSB(:,i,j,6,iBlock) =    &
                  Coef*( B0_DGB(:,i,j,nK  ,iBlock) &
                  +      B0_DGB(:,i,j,nK+1,iBlock) )
          end do; end do
       end if
    end do

    call message_pass_face(                                       &
         3, B0ResChange_DXSB, B0ResChange_DYSB, B0ResChange_DZSB, &
         DoSubtractIn=.false.)

    if(IsCartesian) RETURN

    ! For non-Cartesian grid divide B0ResChange by the coarse face area
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       if(DiLevelNei_IIIB(-1,0,0,iBlock) == -1) then
          do k = 1, nK; do j = 1, nJ
             B0ResChange_DXSB(:,j,k,1,iBlock) =    &
                  B0ResChange_DXSB(:,j,k,1,iBlock) &
                  /max(1e-30, CellFace_DFB(1,1,j,k,iBlock))
          end do; end do
       end if

       if(DiLevelNei_IIIB(+1,0,0,iBlock) == -1) then
          do k = 1, nK; do j = 1, nJ
             B0ResChange_DXSB(:,j,k,2,iBlock) =    &
                  B0ResChange_DXSB(:,j,k,2,iBlock) &
                  /max(1e-30, CellFace_DFB(1,nI+1,j,k,iBlock))
          end do; end do
       end if

       if(DiLevelNei_IIIB(0,-1,0,iBlock) == -1) then
          do k = 1, nK; do i = 1, nI
             B0ResChange_DYSB(:,i,k,3,iBlock) =    &
                  B0ResChange_DYSB(:,i,k,3,iBlock) &
                  /max(1e-30, CellFace_DFB(2,i,1,k,iBlock))
          end do; end do
       end if

       if(DiLevelNei_IIIB(0,+1,0,iBlock) == -1) then
          do k = 1, nK; do i = 1, nI
             B0ResChange_DYSB(:,i,k,4,iBlock) =    &
                  B0ResChange_DYSB(:,i,k,4,iBlock) &
                  /max(1e-30, CellFace_DFB(2,i,nJ+1,k,iBlock))
          end do; end do
       end if

       if(nDim < 3) CYCLE

       if(DiLevelNei_IIIB(0,0,-1,iBlock) == -1) then
          do j = 1, nJ; do i = 1, nI
             B0ResChange_DZSB(:,i,j,5,iBlock) =    &
                  B0ResChange_DZSB(:,i,j,5,iBlock) &
                  /max(1e-30, CellFace_DFB(3,i,j,1,iBlock))
          end do; end do
       end if


       if(DiLevelNei_IIIB(0,0,+1,iBlock) == -1) then
          do j = 1, nJ; do i = 1, nI
             B0ResChange_DZSB(:,i,j,6,iBlock) =    &
                  B0ResChange_DZSB(:,i,j,6,iBlock) &
                  /max(1e-30, CellFace_DFB(3,i,j,nK+1,iBlock))
          end do; end do
       end if
    end do

    if(DoTest)write(*,*)NameSub,' finished'

  end subroutine set_b0_reschange

  !===========================================================================
  subroutine set_b0_source(iBlock)

    ! Calculate div(B0) and curl(B0) for block iBlock

    use BATL_lib, ONLY: IsCartesian, &
         CellSize_DB, FaceNormal_DDFB, CellVolume_GB
    use ModCoordTransform, ONLY: cross_product

    integer, intent(in):: iBlock

    integer:: i, j, k
    real:: DxInv, DyInv, DzInv

    character(len=*), parameter:: NameSub = 'set_b0_source'
    !-----------------------------------------------------------------------
    if(.not.(UseB0 .and. (UseB0Source .or. UseCurlB0))) RETURN

    ! set B0_DX, B0_DY, B0_DZ for this block
    call set_b0_face(iBlock)

    if(IsCartesian)then
       ! Cartesian case
       DxInv = 1/CellSize_DB(x_,iBlock)
       DyInv = 1/CellSize_DB(y_,iBlock)
       DzInv = 1/CellSize_DB(z_,iBlock)

       if(UseB0Source)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             DivB0_C(i,j,k)= &
                  DxInv*(B0_DX(x_,i+1,j,k) - B0_DX(x_,i,j,k)) + &
                  DyInv*(B0_DY(y_,i,j+1,k) - B0_DY(y_,i,j,k)) + &
                  DzInv*(B0_DZ(z_,i,j,k+1) - B0_DZ(z_,i,j,k))
          end do; end do; end do
       end if

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          CurlB0_DC(z_,i,j,k) = &
               DxInv*(B0_DX(y_,i+1,j,k) - B0_DX(y_,i,j,k)) - &
               DyInv*(B0_DY(x_,i,j+1,k) - B0_DY(x_,i,j,k))
          
          CurlB0_DC(y_,i,j,k) = &
               DzInv*(B0_DZ(x_,i,j,k+1) - B0_DZ(x_,i,j,k)) - &
               DxInv*(B0_DX(z_,i+1,j,k) - B0_DX(z_,i,j,k))
          
          CurlB0_DC(x_,i,j,k) = & 
               DyInv*(B0_DY(z_,i,j+1,k) - B0_DY(z_,i,j,k)) - &
               DzInv*(B0_DZ(y_,i,j,k+1) - B0_DZ(y_,i,j,k))
       end do; end do; end do
    else
       if(UseB0Source)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             DivB0_C(i,j,k) = &
                  + sum(FaceNormal_DDFB(:,1,i+1,j,k,iBlock)*B0_DX(:,i+1,j,k)) &
                  - sum(FaceNormal_DDFB(:,1,i  ,j,k,iBlock)*B0_DX(:,i  ,j,k)) &
                  + sum(FaceNormal_DDFB(:,2,i,j+1,k,iBlock)*B0_DY(:,i,j+1,k)) &
                  - sum(FaceNormal_DDFB(:,2,i,j  ,k,iBlock)*B0_DY(:,i,j  ,k))

             if(nDim == 3) DivB0_C(i,j,k) = DivB0_C(i,j,k) &
                  + sum(FaceNormal_DDFB(:,3,i,j,k+1,iBlock)*B0_DZ(:,i,j,k+1)) &
                  - sum(FaceNormal_DDFB(:,3,i,j,k  ,iBlock)*B0_DZ(:,i,j,k  ))

             DivB0_C(i,j,k) = DivB0_C(i,j,k)/CellVolume_GB(i,j,k,iBlock)
          end do; end do; end do
       end if

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          CurlB0_DC(:,i,j,k) =                                          &
               + cross_product(                                         &
               FaceNormal_DDFB(:,1,i+1,j,k,iBlock), B0_DX(:,i+1,j,k))   &
               - cross_product(                                         &
               FaceNormal_DDFB(:,1,i  ,j,k,iBlock), B0_DX(:,i  ,j,k))   &
               + cross_product(                                         &
               FaceNormal_DDFB(:,2,i,j+1,k,iBlock), B0_DY(:,i,j+1,k))   &
               - cross_product(                                         &
               FaceNormal_DDFB(:,2,i,j  ,k,iBlock), B0_DY(:,i,j  ,k))

          if(nDim == 3) CurlB0_DC(:,i,j,k) = CurlB0_DC(:,i,j,k)         &
               + cross_product(                                         &
               FaceNormal_DDFB(:,3,i,j,k+1,iBlock), B0_DZ(:,i,j,k+1))   &
               - cross_product(                                         &
               FaceNormal_DDFB(:,3,i,j,k  ,iBlock), B0_DZ(:,i,j,k  ))

          CurlB0_DC(:,i,j,k) = CurlB0_DC(:,i,j,k)/CellVolume_GB(i,j,k,iBlock)
       end do; end do; end do

    endif

  end subroutine set_b0_source

end module ModB0
