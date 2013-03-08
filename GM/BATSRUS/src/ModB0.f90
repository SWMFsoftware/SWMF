module ModB0

  ! The magnetic field can be split into an analytic and numeric part.
  ! The analytic part is B0. It should satisfy div(B0) = 0, 
  ! and usually curl(B0) = 0 except if UseCurlB0 = .true.
  ! Even if curl(B0) = 0 analytically, the numerical representation
  ! may have a finite curl.
  ! It is often constant in time dB0/dt = 0, but not necessarily.
  ! This module provides data and methods for B0.

  use BATL_size, ONLY: MaxDim, nDim, MaxBlock, nI, nJ, nK, &
       MinI, MaxI, MinJ, MaxJ, MinK, MaxK
  use ModMain, ONLY: UseB, UseB0, UseConstrainB

  implicit none
  SAVE

  private ! except

  public:: init_mod_b0      ! initialize B0 module
  public:: clean_mod_b0     ! clean B0 module
  public:: read_b0_param    ! read UseB0, UseB0Source, UseCurlB0
  public:: set_b0_cell      ! set cell centered B0_DGB
  public:: set_b0_reschange ! make face centered B0 consistent at res.change
  public:: set_b0_face      ! set face centered B0_DX,Y,Z
  public:: set_b0_source    ! set DivB0 and CurlB0
  public:: get_b0           ! get B0 at an arbitrary point

  ! If B0 varies with time it should be update at some frequency
  logical, public:: DoUpdateB0  = UseB
  real,    public:: DtUpdateB0  = 0.0001

  ! Use source terms related to finite div(B0) and curl(B0)
  logical, public:: UseB0Source = UseB

  ! Use source terms related to finite curl(B0) for non-force-free B0 field
  ! If UseCurlB0 is true, then UseB0Source must be true!
  logical, public:: UseCurlB0 = .false.

  ! Radius within which the B0 field is curl free (analytically)
  real, public:: rCurrentFreeB0 = -1.0

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
  subroutine read_b0_param(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'read_b0_param'
    !------------------------------------------------------------------------

    select case(NameCommand)
    case("#USEB0")
       call read_var('UseB0', UseB0)

    case("#DIVBSOURCE")
       call read_var('UseB0Source', UseB0Source)
       
    case("#USECURLB0")
       call read_var('UseCurlB0', UseCurlB0)
       if(UseCurlB0)call read_var('rCurrentFreeB0', rCurrentFreeB0)

    case default
       call stop_mpi(NameSub//': unknown command='//NameCommand)
    end select

    if(UseCurlB0) UseB0Source = .true.

  end subroutine read_b0_param
  !===========================================================================
  subroutine init_mod_b0

    if(.not.allocated(B0_DGB))then
       allocate( &
            B0_DGB(MaxDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock), &
            B0ResChange_DXSB(MaxDim,nJ,nK,1:2,MaxBlock),           &
            B0ResChange_DYSB(MaxDim,nI,nK,3:4,MaxBlock),           &
            B0ResChange_DZSB(MaxDim,nI,nJ,5:6,MaxBlock))
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

    if(allocated(B0_DGB)) deallocate(B0_DGB, &
         B0ResChange_DXSB, B0ResChange_DYSB, B0ResChange_DZSB)

    if(allocated(DivB0_C))   deallocate(DivB0_C)
    if(allocated(CurlB0_DC)) deallocate(CurlB0_DC)

  end subroutine clean_mod_b0
  !===========================================================================
  subroutine set_b0_cell(iBlock)

    ! Calculate the cell centered B0 for block iBlock

    use ModProcMH, ONLY: iProc
    use ModMain,   ONLY: ProcTest, BlkTest, iTest, jTest, kTest
    use BATL_lib,  ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, Xyz_DGB

    integer, intent(in) :: iBlock

    integer :: i, j, k

    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------
    if(iProc==PROCtest.and.iBlock==BLKtest)then
       call set_oktest('set_b0_cell',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    endif

    do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
       call get_b0(Xyz_DGB(:,i,j,k,iBlock), B0_DGB(:,i,j,k,iBlock))
    end do; end do; end do

    if(DoTestMe)write(*,*)'B0*Cell_BLK=',&
         B0_DGB(:,Itest,Jtest,Ktest,BLKtest)

  end subroutine set_b0_cell

  !===========================================================================

  subroutine set_b0_face(iBlock)

    ! Calculate the face centered B0 for block iBlock

    use ModParallel, ONLY: NeiLev
    use BATL_lib,    ONLY: nDim

    integer,intent(in)::iBlock
    !-------------------------------------------------------------------------
    if(.not.UseB0) RETURN

    ! Average cell centered B0_DGB to the face centered B0_DX,Y,Z arrays
    if(UseConstrainB)then
       B0_DX = 0.5*(B0_DGB(:,0:nI  ,0:nJ+1,0:nK+1,iBlock) &
            +       B0_DGB(:,1:nI+1,0:nJ+1,0:nK+1,iBlock))

       B0_DY = 0.5*(B0_DGB(:,0:nI+1,0:nJ  ,0:nK+1,iBlock) &
            +       B0_DGB(:,0:nI+1,1:nJ+1,0:nK+1,iBlock))

       B0_DZ = 0.5*(B0_DGB(:,0:nI+1,0:nJ+1,0:nK  ,iBlock) &
            +       B0_DGB(:,0:nI+1,0:nJ+1,1:nK+1,iBlock))

    else
       B0_DX = 0.5*(B0_DGB(:,0:nI  ,1:nJ,1:nK,iBlock) &
            +       B0_DGB(:,1:nI+1,1:nJ,1:nK,iBlock))

       if(nDim >= 2) &
            B0_DY = 0.5*(B0_DGB(:,1:nI,0:nJ  ,1:nK,iBlock) &
            +            B0_DGB(:,1:nI,1:nJ+1,1:nK,iBlock))

       if(nDim >= 3) &
            B0_DZ = 0.5*(B0_DGB(:,1:nI,1:nJ,0:nK  ,iBlock) &
               +         B0_DGB(:,1:nI,1:nJ,1:nK+1,iBlock))
    end if
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
         IsCartesian, CellFace_DFB, message_pass_face

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
    if(IsCartesian .and. nDim==2) Coef = 0.25
    if(IsCartesian .and. nDim==3) Coef = 0.125

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

    use ModSize,  ONLY: x_, y_, z_
    use BATL_lib, ONLY: IsCartesian, IsRzGeometry, &
         CellSize_DB, CellFace_DFB, FaceNormal_DDFB, CellVolume_GB, Xyz_DGB
    use ModCoordTransform, ONLY: cross_product

    integer, intent(in):: iBlock

    integer:: i, j, k
    real:: DxInv, DyInv, DzInv

    character(len=*), parameter:: NameSub = 'set_b0_source'
    !-----------------------------------------------------------------------
    if(.not.(UseB0 .and. UseB0Source)) RETURN

    ! set B0_DX, B0_DY, B0_DZ for this block
    call set_b0_face(iBlock)

    if(IsCartesian)then
       if(nDim == 2)then
          DxInv = 1/CellSize_DB(x_,iBlock)
          DyInv = 1/CellSize_DB(y_,iBlock)

          k = 1
          do j = 1, nJ; do i = 1, nI
             DivB0_C(i,j,k)= &
                  DxInv*(B0_DX(x_,i+1,j,k) - B0_DX(x_,i,j,k)) + &
                  DyInv*(B0_DY(y_,i,j+1,k) - B0_DY(y_,i,j,k))

             CurlB0_DC(x_,i,j,k) = & 
                  +DyInv*(B0_DY(z_,i,j+1,k) - B0_DY(z_,i,j,k))

             CurlB0_DC(y_,i,j,k) = &
                  -DxInv*(B0_DX(z_,i+1,j,k) - B0_DX(z_,i,j,k))

             CurlB0_DC(z_,i,j,k) = &
                  DxInv*(B0_DX(y_,i+1,j,k) - B0_DX(y_,i,j,k)) - &
                  DyInv*(B0_DY(x_,i,j+1,k) - B0_DY(x_,i,j,k))

          end do; end do
       else
          DxInv = 1/CellSize_DB(x_,iBlock)
          DyInv = 1/CellSize_DB(y_,iBlock)
          DzInv = 1/CellSize_DB(z_,iBlock)

          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             DivB0_C(i,j,k)= &
                  DxInv*(B0_DX(x_,i+1,j,k) - B0_DX(x_,i,j,k)) + &
                  DyInv*(B0_DY(y_,i,j+1,k) - B0_DY(y_,i,j,k)) + &
                  DzInv*(B0_DZ(z_,i,j,k+1) - B0_DZ(z_,i,j,k))

             CurlB0_DC(x_,i,j,k) = & 
                  DyInv*(B0_DY(z_,i,j+1,k) - B0_DY(z_,i,j,k)) - &
                  DzInv*(B0_DZ(y_,i,j,k+1) - B0_DZ(y_,i,j,k))

             CurlB0_DC(y_,i,j,k) = &
                  DzInv*(B0_DZ(x_,i,j,k+1) - B0_DZ(x_,i,j,k)) - &
                  DxInv*(B0_DX(z_,i+1,j,k) - B0_DX(z_,i,j,k))

             CurlB0_DC(z_,i,j,k) = &
                  DxInv*(B0_DX(y_,i+1,j,k) - B0_DX(y_,i,j,k)) - &
                  DyInv*(B0_DY(x_,i,j+1,k) - B0_DY(x_,i,j,k))

          end do; end do; end do
       end if
    elseif(IsRzGeometry)then
       k = 1
       do j = 1, nJ; do i = 1, nI
          DivB0_C(i,j,k)= ( &
               + CellFace_DFB(x_,i+1,j,k,iBlock)*B0_DX(x_,i+1,j,k)   &
               - CellFace_DFB(x_,i  ,j,k,iBlock)*B0_DX(x_,i  ,j,k)   &
               + CellFace_DFB(y_,i,j+1,k,iBlock)*B0_DY(y_,i,j+1,k)   &
               - CellFace_DFB(y_,i,j  ,k,iBlock)*B0_DY(y_,i,j  ,k) ) &
               /CellVolume_GB(i,j,k,iBlock)

          CurlB0_DC(x_,i,j,k) = ( &
               + CellFace_DFB(y_,i,j+1,k,iBlock)*B0_DY(z_,i,j+1,k)   &
               - CellFace_DFB(y_,i,j  ,k,iBlock)*B0_DY(z_,i,j  ,k) ) &
               /CellVolume_GB(i,j,k,iBlock)

          CurlB0_DC(y_,i,j,k) = &
               -DxInv*(B0_DX(z_,i+1,j,k) - B0_DX(z_,i,j,k))

          CurlB0_DC(z_,i,j,k) = &
               DxInv*(B0_DX(y_,i+1,j,k) - B0_DX(y_,i,j,k)) - ( &
               + CellFace_DFB(y_,i,j+1,k,iBlock)*B0_DY(x_,i,j+1,k)   &
               - CellFace_DFB(y_,i,j  ,k,iBlock)*B0_DY(x_,i,j  ,k) ) &
               /CellVolume_GB(i,j,k,iBlock) &
               + B0_DGB(x_,i,j,k,iBlock)/Xyz_DGB(y_,i,j,k,iBlock)
       end do; end do
    else
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
  !============================================================================
  subroutine get_b0(Xyz_D, B0_D)

    ! Get B0 at location Xyz_D

    use ModMain,          ONLY : Time_Simulation, NameThisComp, &
         TypeCoordSystem, IsStandAlone
    use ModPhysics,       ONLY: Si2No_V, UnitB_, DipoleStrengthSi
    use CON_planet_field, ONLY: get_planet_field
    use ModMain,          ONLY: UseBody2
    use ModMain,          ONLY: UseUserB0, UseMagnetogram
    use ModMagnetogram,   ONLY: get_magnetogram_field

    real, intent(in) :: Xyz_D(3)
    real, intent(out):: B0_D(3)
    !-------------------------------------------------------------------------

    if(UseMagnetogram)then
       call get_magnetogram_field(Xyz_D(1), Xyz_D(2), Xyz_D(3), B0_D)
       B0_D = B0_D*Si2No_V(UnitB_)
    elseif(UseUserB0)then
       call user_get_b0(Xyz_D(1), Xyz_D(2), Xyz_D(3), B0_D)
    elseif(IsStandAlone .and. DipoleStrengthSi==0.0)then
       B0_D = 0.0
       RETURN
    elseif(NameThisComp=='GM' .and. nDim == 3)then
       call get_planet_field(Time_Simulation, Xyz_D, TypeCoordSystem//' NORM',&
            B0_D)
       B0_D = B0_D*Si2No_V(UnitB_)
    else 
       call get_b0_multipole(Xyz_D, B0_D) 
    end if
    if(UseBody2)call add_b0_body2(Xyz_D, B0_D)

  end subroutine get_b0

  !============================================================================

  subroutine get_b0_multipole(Xyz_D, B0_D)

    ! Calculate the multipole based B0_D field at location Xyz_D.

    ! Note that ThetaTilt is positive when the north magnetic pole
    ! points AWAY from the sun.

    use ModPhysics, ONLY: Bdp, CosThetaTilt, SinThetaTilt, cTiny, &
         rBody, Qqp, Oop
    
    real, intent(in) :: Xyz_D(3)
    real, intent(out):: B0_D(3)
    
    integer :: i, j, k, l
    real :: x, y, r2, rInv, r2Inv, r3Inv, r5Inv, r7Inv
    real :: XyzTilt_D(3), b_D(3), Bx, By
    real :: Dp, Tmp1, Tmp2, Dipole_D(3)

    logical :: DoQuadrupole, DoOctupole
    !-------------------------------------------------------------------------
    ! Determine radial distance and powers of it
    r2 = sum(Xyz_D(1:nDim)**2)

    ! Avoid calculating B0 inside a critical radius = 1.E-6*Rbody
    if(r2 <= (cTiny*rBody)**2)then
       B0_D = 0.0
       RETURN
    end if

    if(nDim == 2) then
       !
       ! 2D magnetic dipole is implemented from F.R. Cardoso et.al.'s paper
       ! "2D MHD simulation of the magnetic dipole tilt and IMF influence 
       ! on the magnetosphere"
       ! except that they define the dipole pointing down, while this one 
       ! points upward.

       x = Xyz_D(1)
       y = Xyz_D(2)
       Dp = Bdp/r2**2

       ! Dipole aligned with the +Y axis
       Bx = Dp*2*x*y
       By = Dp*(y**2 - x**2)

       ! Rotate the field with -ThetaTilt around the Z axis
       B0_D(1) =  CosThetaTilt*Bx + SinThetaTilt*By
       B0_D(2) = -SinThetaTilt*Bx + CosThetaTilt*By
       B0_D(3) = 0.0
       RETURN
    end if
    
    rInv  = 1.0/sqrt(r2)
    r2Inv = rInv**2
    r3Inv = rInv*r2Inv
    
    ! Compute dipole moment of the intrinsic magnetic field B0.

    Dipole_D = (/ -SinThetaTilt*Bdp, 0.0, CosThetaTilt*Bdp /) 

    Dp = 3*sum(Dipole_D*Xyz_D)*r2Inv

    B0_D = (Dp*Xyz_D - Dipole_D)*r3Inv

    DoQuadrupole=any(abs(Qqp)>cTiny)
    DoOctupole  =any(abs(Oop)>cTiny)

    if(DoQuadrupole .or. DoOctupole)then
       ! Compute the xx's in the tilted reference frame aligned with
       ! the magnetic field.
       XyzTilt_D(1) = CosThetaTilt*Xyz_D(1) + SinThetaTilt*Xyz_D(3)
       XyzTilt_D(2) = Xyz_D(2)
       XyzTilt_D(3)= -SinThetaTilt*Xyz_D(1) + CosThetaTilt*Xyz_D(3)

       r5Inv = r3Inv*r2Inv
       r7Inv = r5Inv*r2Inv
    end if

    if(DoQuadrupole)then
       ! Compute quadrupole moment of the intrinsic 
       ! magnetic field B0.
       do k=1,3
          Tmp1 = 0.0
          Tmp2 = 0.0
          do i=1,3
             Tmp1 = Tmp1 + Qqp(k,i)*XyzTilt_D(i)
             do j=1,3
                Tmp2 = Tmp2 + Qqp(i,j)*XyzTilt_D(i)*XyzTilt_D(j)*XyzTilt_D(k)
             end do
          end do
          b_D(k) = 2.5*Tmp2*r7Inv - Tmp1*r5Inv
       end do

       B0_D(1) = B0_D(1) + CosThetaTilt*b_D(1) - SinThetaTilt*b_D(3) 
       B0_D(2) = B0_D(2) + b_D(2)
       B0_D(3) = B0_D(3) + SinThetaTilt*b_D(1) + CosThetaTilt*b_D(3)
    end if

    if(DoOctupole)then
       ! Compute octupole moment of the intrinsic 
       ! magnetic field B0.
       do k = 1, 3
          Tmp1 = 0.0
          Tmp2 = 0.0
          do i = 1, 3
             do j = 1, 3
                Tmp1 = Tmp1 + Oop(i,j,k)*XyzTilt_D(i)*XyzTilt_D(j)
                do l = 1, 3
                   Tmp2 = Tmp2 + Oop(i,j,l) &
                        *XyzTilt_D(i)*XyzTilt_D(j)*XyzTilt_D(l)*XyzTilt_D(k)
                end do
             end do
          end do
          b_D(k) = 7.0*Tmp2*r7Inv - 3.0*Tmp1*r5Inv
       end do

       B0_D(1) = B0_D(1) + CosThetaTilt*b_D(1) - SinThetaTilt*b_D(3) 
       B0_D(2) = B0_D(2) + b_D(2)
       B0_D(3) = B0_D(3) + SinThetaTilt*b_D(1) + CosThetaTilt*b_D(3)
    end if

  end subroutine get_b0_multipole


  !============================================================================
  subroutine add_b0_body2(XyzIn_D, B0_D)

    !\
    ! If there is a second body that has a magnetic field the contribution
    ! to the field from the second body should be computed here (inside the
    ! if block.
    !/
    use ModPhysics
    use ModNumConst, ONLY: cTiny

    real, intent(in)   :: XyzIn_D(3)
    real, intent(inout):: B0_D(3)

    real :: Xyz_D(3),R0,rInv,r2Inv,r3Inv,Dp
    !--------------------------------------------------------------------------
    !\
    ! Determine normalized relative coordinates and radial distance from body 2
    !/
    Xyz_D = (XyzIn_D - (/xBody2, yBody2, zBody2/))/rBody2

    R0 = sqrt(sum(Xyz_D**2))

    ! Avoid calculating B0 inside a critical normalized radius = cTiny
    if(R0 <= cTiny) RETURN

    rInv  = 1/R0
    r2Inv = rInv**2
    r3Inv = rInv*r2Inv

    ! Add dipole field of the second body
    Dp = sum(BdpBody2_D*Xyz_D)*3*r2Inv

    B0_D = B0_D + (Dp*Xyz_D - BdpBody2_D)*r3Inv

  end subroutine add_b0_body2


end module ModB0
