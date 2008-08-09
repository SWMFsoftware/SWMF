module ModB0
  use ModSize
  use ModUtilities,ONLY: check_allocate
  implicit none
  SAVE
  !\
  !Cell-centered B0 field vector
  !/
  real,allocatable :: B0_DGB(:,:,:,:,:)

  !\
  ! Face-centered intrinsic magnetic field arrays.
  !/
  real :: B0_DX(nDim,2-gcn:nI+gcn, 0:nJ+1, 0:nK+1)=0.00 

  real :: B0_DY(nDim,0:nI+1, 2-gcn:nJ+gcn, 0:nK+1)=0.00

  real :: B0_DZ(nDim,0:nI+1, 0:nJ+1, 2-gcn:nK+gcn)=0.00

  !\
  ! B0 field at the resolution interfaces
  !/
  real,allocatable :: B0ResChange_DXSB(:,:,:,:,:)
  real,allocatable :: B0ResChange_DYSB(:,:,:,:,:)
  real,allocatable :: B0ResChange_DZSB(:,:,:,:,:)

  real,allocatable :: CurlB0_DCB(:,:,:,:,:)
  real,allocatable :: DivB0_CB(  :,:,:,:)
  real,allocatable :: NormB0_CB( :,:,:,:)
contains
  !=========================================================!
  subroutine allocate_b0_arrays
    integer::iError
    !--------------!
    if(allocated(B0_DGB))return
    allocate(B0_DGB(nDim,1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn, MaxBlock),&
         B0ResChange_DXSB(nDim,nJ,nK,East_ :West_ ,MaxBlock),            &
         B0ResChange_DYSB(nDim,nI,nK,South_:North_,MaxBlock),            &
         B0ResChange_DZSB(nDim,nI,nJ,Bot_  :Top_  ,MaxBlock),            &
         stat=iError)
    call check_allocate(iError,'B0 arrays are not allocated')
    B0_DGB=0.00 
    B0ResChange_DXSB=0.00;B0ResChange_DYSB=0.00;B0ResChange_DZSB=0.00
  end subroutine allocate_b0_arrays
  !=========================================================!
  subroutine allocate_b0_source_arrays
    integer::iError
    !--------------!
    if(allocated(DivB0_CB))return
    allocate(DivB0_CB(nI, nJ, nK, MaxBlock),&
         CurlB0_DCB(3,nI, nJ, nK, MaxBlock),stat=iError)
    call check_allocate(iError,'B0 source arrays are not allocated')
    DivB0_CB=0.00;CurlB0_DCB=0.00
  end subroutine allocate_b0_source_arrays
  !=========================================================!
  subroutine allocate_b0_norm
    integer::iError
    !--------------!
    if(allocated(NormB0_CB))return
    allocate(NormB0_CB(nI, nJ, nK, MaxBlock),stat=iError)
    call check_allocate(iError,'B0 norm array is not allocated')
    NormB0_CB=0.00
  end subroutine allocate_b0_norm
  !=========================================================!
  subroutine set_b0_face(iBlock)
    use ModParallel, ONLY : neiLev,&
         neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
    use ModMain, ONLY: UseB0
    integer,intent(in)::iBlock
    if(.not.UseB0)then
       B0_DX=0.00; B0_DY=0.00; B0_DZ=0.00; return
    end if
    !\
    ! Face averages (X)
    !/
    B0_DX = 0.50 * (B0_DGB(:,2-gcn:nI+gcn  ,0:nJ+1,0:nK+1,iBlock)+&
         B0_DGB(:,1-gcn:nI+gcn-1,0:nJ+1,0:nK+1,iBlock))
    !\
    ! Face averages (Y)
    !/
    B0_DY = 0.50 * (B0_DGB(:,0:nI+1,2-gcn:nJ+gcn  ,0:nK+1,iBlock)+&
         B0_DGB(:,0:nI+1,1-gcn:nJ+gcn-1,0:nK+1,iBlock))

    !\
    ! Face averages (Z)
    !/
    B0_DZ = 0.50 * (B0_DGB(:,0:nI+1,0:nJ+1,2-gcn:nK+gcn  ,iBlock)+&
         B0_DGB(:,0:nI+1,0:nJ+1,1-gcn:nK+gcn-1,iBlock))

    if(all(neiLev(:,iBlock)/=-1))return
    if (neiLeast(iBlock)==-1) &
         B0_DX(:,1   ,1:nJ,1:nK)=B0ResChange_DXSB(:,:,:,East_,iBlock)
    if (neiLwest(iBlock)==-1) &
         B0_DX(:,nI+1,1:nJ,1:nK)=B0ResChange_DXSB(:,:,:,West_,iBlock)
    if (neiLsouth(iBlock)==-1) &
         B0_DY(:,1:nI,1   ,1:nK)=B0ResChange_DYSB(:,:,:,South_,iBlock)
    if (neiLnorth(iBlock)==-1) &
         B0_DY(:,1:nI,1+nJ,1:nK)=B0ResChange_DYSB(:,:,:,North_,iBlock)
    if (neiLbot(iBlock)==-1) &
         B0_DZ(:,1:nI,1:nJ,1   )=B0ResChange_DZSB(:,:,:,Bot_,iBlock)
    if (neiLtop(iBlock)==-1) &
         B0_DZ(:,1:nI,1:nJ,1+nK)=B0ResChange_DZSB(:,:,:,Top_,iBlock)
  end subroutine set_b0_face
end module ModB0
