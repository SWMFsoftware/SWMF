!^CFG COPYRIGHT UM
!==============================================================================
module ModHeatConduction

  implicit none
  save

  private ! except

  ! Public methods
  public :: read_heatconduction_param
  public :: init_heat_conduction
  public :: get_heat_flux
  public :: impl_jac_heat_conduction

  ! Logical for adding parallel heat conduction
  Logical, public :: UseParallelConduction = .false.
  logical, public :: IsNewBlockHeatConduction = .true.

  ! electron temperature used for calculating field parallel heat flux
  real, allocatable :: Te_G(:,:,:)

contains

  !============================================================================

  subroutine read_heatconduction_param(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand
    character(len=*), parameter :: NameSub = 'read_heatconduction_param'
    !--------------------------------------------------------------------------

    select case(NameCommand)
    case("#PARALLELCONDUCTION")
       call read_var('UseParallelConduction', UseParallelConduction)
    case default
       call stop_mpi(NameSub//' invalid NameCommand='//NameCommand)
    end select

  end subroutine read_heatconduction_param

  !============================================================================

  subroutine init_heat_conduction

    use ModSize, ONLY: nI, nJ, nK

    character(len=*), parameter :: NameSub = 'init_heat_conduction'
    !--------------------------------------------------------------------------

    if(allocated(Te_G)) RETURN

    allocate(Te_G(-1:nI+2,-1:nJ+2,-1:nK+2))

  end subroutine init_heat_conduction

  !============================================================================

  subroutine get_heat_flux(iDir, i, j, k, iBlock, State_V, HeatFlux_D)

    use ModB0,         ONLY: B0_DX, B0_DY, B0_DZ
    use ModMain,       ONLY: UseB0
    use ModNumConst,   ONLY: cTolerance
    use ModVarIndexes, ONLY: nVar, Bx_, Bz_

    integer, intent(in) :: iDir, i, j, k, iBlock
    real,    intent(in) :: State_V(nVar)
    real,    intent(out):: HeatFlux_D(3)

    real :: B_D(3), Bunit_D(3), Bnorm
    real :: FaceGrad_D(3)
    !--------------------------------------------------------------------------

    if(UseB0)then
       select case(iDir)
       case(1)
          B_D = State_V(Bx_:Bz_) + B0_DX(:,i,j,k)
       case(2)
          B_D = State_V(Bx_:Bz_) + B0_DY(:,i,j,k)
       case(3)
          B_D = State_V(Bx_:Bz_) + B0_DZ(:,i,j,k)
       end select
    else
       B_D = State_V(Bx_:Bz_)
    end if

    ! The magnetic field should nowhere be zero. The following fix will
    ! turn the magnitude of the field direction to zero.
    Bnorm = max(sqrt(sum(B_D**2)),cTolerance)
    Bunit_D = B_D/Bnorm

    call calc_face_gradient(iDir, i, j, k, iBlock, FaceGrad_D)

    HeatFlux_D = -Bunit_D*dot_product(Bunit_D,FaceGrad_D)

  end subroutine get_heat_flux

  !============================================================================

  subroutine set_block_scalar(Scalar_G, iBlock)

    ! correct the ghostcells of the given scalar field on iBlock

    use ModParallel, ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot, BlkNeighborLev, NOBLK
    use ModSize,     ONLY: nI, nJ, nK

    integer, intent(in) :: iBlock
    real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2), intent(inout) :: Scalar_G

    real, parameter :: c0 = 0.5, p0 = 1./6., F1 = 1./3.

    integer :: i1, j1, k1, i2, j2, k2, iC, jC, kC, iSide, jSide, kSide
    integer :: iL, iR, jL, jR, kL, kR
    integer :: ip, jp, kp

    real    :: Scalar1_G(0:nI+1,0:nJ+1,0:nK+1) ! temporary array
    logical :: IsEqualLevel_G(0:nI+1,0:nJ+1,0:nK+1)
    !--------------------------------------------------------------------------
    IsNewBlockHeatConduction = .false.

    Scalar1_G = Scalar_G(0:nI+1,0:nJ+1,0:nK+1)

    do kSide = -1,1; do jSide = -1,1; do iSide = -1,1
       if(iSide==0)then
          iL = 1; iR = nI
       elseif(iSide==1)then
          iL = nI+1; iR = iL
       else
          iL = 0; iR = iL
       end if
       if(jSide==0)then
          jL = 1; jR = nJ
       elseif(jSide==1)then
          jL = nJ+1; jR = jL
       else
          jL = 0; jR = jL
       end if
       if(kSide==0)then
          kL = 1; kR = nK
       elseif(kSide==1)then
          kL = nK+1; kR = kL
       else
          kL = 0; kR = kL
       end if
       if( BlkNeighborLev(iSide, jSide, kSide, iBlock) == 0 )then
          IsEqualLevel_G(iL:iR,jL:jR,kL:kR) = .true.
       else
          IsEqualLevel_G(iL:iR,jL:jR,kL:kR) = .false.
       end if
    end do; end do; end do

    ! Do six faces
    if(NeiLeast(iBlock) == 1)then
       do k1=1, nK, 2; do j1=1, nJ, 2; do k2 = k1,k1+1; do j2 = j1,j1+1
          jp = 3*j2 - 2*j1 -1
          kp = 3*k2 - 2*k1 -1
          if(IsEqualLevel_G(0,jp,kp))then
             Scalar_G(0,j2,k2) = c0*Scalar1_G(0,j2,k2) &
                  + 0.25*Scalar1_G(0,jp,kp) + 0.25*Scalar_G(1,j2,k2)
          else
             Scalar_G(0,j2,k2) = c0*Scalar1_G(0,j2,k2) &
                  + p0*Scalar1_G(0,jp,kp) + F1*Scalar_G(1,j2,k2)
          end if
       end do; end do; end do; end do
    end if

    if(NeiLwest(iBlock) == 1)then
       do k1=1, nK, 2; do j1=1, nJ, 2; do k2 = k1,k1+1; do j2 = j1,j1+1
          jp = 3*j2 - 2*j1 -1
          kp = 3*k2 - 2*k1 -1
          if(IsEqualLevel_G(nI+1,jp,kp))then
             Scalar_G(nI+1,j2,k2) = c0*Scalar1_G(nI+1,j2,k2) &
                  + 0.25*Scalar1_G(nI+1,jp,kp) + 0.25*Scalar_G(nI,j2,k2)
          else
             Scalar_G(nI+1,j2,k2) = c0*Scalar1_G(nI+1,j2,k2) &
                  + p0*Scalar1_G(nI+1,jp,kp) + F1*Scalar_G(nI,j2,k2)
          end if
       end do; end do; end do; end do
    end if

    if(NeiLsouth(iBlock) == 1)then
       do k1=1, nK, 2; do i1=1, nI, 2; do k2 = k1,k1+1; do i2 = i1,i1+1
          ip = 3*i2 - 2*i1 -1
          kp = 3*k2 - 2*k1 -1
          if(IsEqualLevel_G(ip,0,kp))then
             Scalar_G(i2,0,k2) = c0*Scalar1_G(i2,0,k2) &
                  + 0.25*Scalar1_G(ip,0,kp) + 0.25*Scalar_G(i2,1,k2)
          else
             Scalar_G(i2,0,k2) = c0*Scalar1_G(i2,0,k2) &
                  + p0*Scalar1_G(ip,0,kp) + F1*Scalar_G(i2,1,k2)
          end if
       end do; end do; end do; end do
    end if

    if(NeiLnorth(iBlock) == 1)then
       do k1=1, nK, 2; do i1=1, nI, 2; do k2 = k1,k1+1; do i2 = i1,i1+1
          ip = 3*i2 - 2*i1 -1
          kp = 3*k2 - 2*k1 -1
          if(IsEqualLevel_G(ip,nJ+1,kp))then
             Scalar_G(i2,nJ+1,k2) = c0*Scalar1_G(i2,nJ+1,k2) &
                  + 0.25*Scalar1_G(ip,nJ+1,kp) + 0.25*Scalar_G(i2,nJ,k2)
          else
             Scalar_G(i2,nJ+1,k2) = c0*Scalar1_G(i2,nJ+1,k2) &
                  + p0*Scalar1_G(ip,nJ+1,kp) + F1*Scalar_G(i2,nJ,k2)
          end if
       end do; end do; end do; end do
    end if

    if(NeiLbot(iBlock) == 1)then
       do j1=1, nJ, 2; do i1=1, nI, 2; do j2 = j1,j1+1; do i2 = i1,i1+1
          ip = 3*i2 - 2*i1 -1
          jp = 3*j2 - 2*j1 -1
          if(IsEqualLevel_G(ip,jp,0))then
             Scalar_G(i2,j2,0) = c0*Scalar1_G(i2,j2,0) &
                  + 0.25*Scalar1_G(ip,jp,0) + 0.25*Scalar_G(i2,j2,1)
          else
             Scalar_G(i2,j2,0) = c0*Scalar1_G(i2,j2,0) &
                  + p0*Scalar1_G(ip,jp,0) + F1*Scalar_G(i2,j2,1)
          end if
       end do; end do; end do; end do
    end if

    if(NeiLtop(iBlock) == 1)then
       do j1=1, nJ, 2; do i1=1, nI, 2; do j2 = j1,j1+1; do i2 = i1,i1+1
          ip = 3*i2 - 2*i1 -1
          jp = 3*j2 - 2*j1 -1
          if(IsEqualLevel_G(ip,jp,nK+1))then
             Scalar_G(i2,j2,nK+1) = c0*Scalar1_G(i2,j2,nK+1) &
                  + 0.25*Scalar1_G(ip,jp,nK+1) + 0.25*Scalar_G(i2,j2,nK)
          else
             Scalar_G(i2,j2,nK+1) = c0*Scalar1_G(i2,j2,nK+1) &
                  + p0*Scalar1_G(ip,jp,nK+1) + F1*Scalar_G(i2,j2,nK)
          end if
       end do; end do; end do; end do
    end if

    ! Do 12 edges
    ! 4 X edges
    do kSide = -1,1,2; do jSide = -1,1,2
       if(  BlkNeighborLev(0, jSide, kSide, iBlock) /= 1 .and. .not. ( &
            BlkNeighborLev(0, jSide, kSide, iBlock) == NOBLK .and. ( &
            BlkNeighborLev(0, jSide, 0, iBlock) == 1 .or. &
            BlkNeighborLev(0, 0, kSide, iBlock) == 1))) CYCLE

       j1=1; if(jSide==1) j1=nJ; jC = j1+jSide
       k1=1; if(kSide==1) k1=nK; kC = k1+kSide
       do i1 = 1,nI,2; do i2 = i1, i1+1
          ip = 3*i2 - 2*i1 -1
          if(IsEqualLevel_G(ip,jC,kC))then
             Scalar_G(i2,jC,kC) = c0*Scalar1_G(i2,jC,kC) &
                  + 0.25*Scalar1_G(ip,jC,kC) + 0.25*Scalar_G(i2,j1,k1)
          else
             Scalar_G(i2,jC,kC) = c0*Scalar1_G(i2,jC,kC) &
                  + p0*Scalar1_G(ip,jC,kC) + F1*Scalar_G(i2,j1,k1)
          end if
       end do; end do
    end do; end do
    ! 4 Y edges
    do kSide = -1,1,2; do iSide = -1,1,2
       if(  BlkNeighborLev(iSide, 0, kSide, iBlock) /= 1 .and. .not. ( &
            BlkNeighborLev(iSide, 0, kSide, iBlock) == NOBLK .and. ( &
            BlkNeighborLev(iSide, 0, 0, iBlock) == 1 .or. &
            BlkNeighborLev(0, 0, kSide, iBlock) == 1))) CYCLE

       i1=1; if(iSide==1) i1=nI; iC = i1+iSide
       k1=1; if(kSide==1) k1=nK; kC = k1+kSide
       do j1 = 1, nJ, 2; do j2 = j1, j1+1
          jp = 3*j2 - 2*j1 -1
          if(IsEqualLevel_G(iC,jp,kC))then
             Scalar_G(iC,j2,kC) = c0*Scalar1_G(iC,j2,kC) &
                  + 0.25*Scalar1_G(iC,jp,kC) + 0.25*Scalar_G(i1,j2,k1)
          else
             Scalar_G(iC,j2,kC) = c0*Scalar1_G(iC,j2,kC) &
                  + p0*Scalar1_G(iC,jp,kC) + F1*Scalar_G(i1,j2,k1)
          end if
       end do; end do
    end do; end do
    ! 4 Z edges
    do jSide = -1,1,2; do iSide = -1,1,2
       if(  BlkNeighborLev(iSide, jSide, 0, iBlock) /= 1.and. .not. ( &
            BlkNeighborLev(iSide, jSide, 0, iBlock) == NOBLK .and. ( &
            BlkNeighborLev(iSide, 0, 0, iBlock) == 1 .or. &
            BlkNeighborLev(0, jSide, 0, iBlock) == 1))) CYCLE

       i1=1; if(iSide==1) i1=nI; iC = i1+iSide
       j1=1; if(jSide==1) j1=nJ; jC = j1+jSide
       do k1 = 1, nK, 2 ; do k2 = k1, k1 + 1
          kp = 3*k2 - 2*k1 -1
          if(IsEqualLevel_G(iC,jC,kp))then
             Scalar_G(iC,jC,k2) = c0*Scalar1_G(iC,jC,k2) &
                  + 0.25*Scalar1_G(iC,jC,kp) + 0.25*Scalar_G(i1,j1,k2)
          else
             Scalar_G(iC,jC,k2) = c0*Scalar1_G(iC,jC,k2) &
                  + p0*Scalar1_G(iC,jC,kp) + F1*Scalar_G(i1,j1,k2)
          end if
       end do; end do         
    end do; end do

  end subroutine set_block_scalar

  !============================================================================

  subroutine calc_face_gradient(iDir, i, j, k, iBlock, FaceGrad_D) 

    use ModAdvance,    ONLY: State_VGB
    use ModCovariant,  ONLY: UseCovariant
    use ModGeometry,   ONLY: Dx_Blk, Dy_Blk, Dz_Blk
     use ModMain,      ONLY: x_, y_, z_
    use ModParallel,   ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot, BlkNeighborLev
    use ModSize,       ONLY: nI, nJ, nK
    use ModVarIndexes, ONLY: p_, Rho_

    integer, intent(in) :: iDir, i, j, k, iBlock
    real, intent(out) :: FaceGrad_D(3)

    integer :: iL, iR, jL, jR, kL, kR
    real :: Ax, Bx, Cx, Ay, By, Cy, Az, Bz, Cz

    Real :: InvDx, InvDy, InvDz
    !--------------------------------------------------------------------------
    InvDx = 1.0/Dx_Blk(iBlock)
    InvDy = 1.0/Dy_Blk(iBlock)
    InvDz = 1.0/Dz_Blk(iBlock)

    if(IsNewBlockHeatConduction)then
       Te_G = State_VGB(p_,:,:,:,iBlock)/State_VGB(Rho_,:,:,:,iBlock)
       call set_block_scalar(Te_G, iBlock)
       if(UseCovariant)then
          ! call ...
       end if
    end if

    ! Central difference with averaging in orthogonal direction
    iR = i+1; iL = i-1; 
    jR = j+1; jL = j-1; 
    kR = k+1; kL = k-1; 

    Ax = -0.25*InvDx; Bx = 0.0; Cx = +0.25*InvDx
    Ay = -0.25*InvDy; By = 0.0; Cy = +0.25*InvDy
    Az = -0.25*InvDz; Bz = 0.0; Cz = +0.25*InvDz

    if(i==1)then
       if(NeiLeast(iBlock)==-1 &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev(-1,-1, 0,iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev(-1, 1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev(-1, 0,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev(-1, 0, 1,iBlock)==-1)) &
            )then
          iL = i+1; iR = i+2; Ax=InvDx; Bx=-0.75*InvDx; Cx=-0.25*InvDx
       end if
    elseif(i==nI)then
       if(NeiLwest(iBlock)==-1 &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev( 1,-1, 0,iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev( 1, 1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev( 1, 0,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev( 1, 0, 1,iBlock)==-1)) &
            )then
          iL = i-1; iR = i-2; Ax=-InvDx; Bx=0.75*InvDx; Cx=0.25*InvDx
       end if
    end if

    if(j==1)then
       if(NeiLsouth(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1,-1, 0,iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1,-1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev( 0,-1,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev( 0,-1, 1,iBlock)==-1)) &
            )then
          jL = j+1; jR = j+2; Ay=InvDy; By=-0.75*InvDy; Cy=-0.25*InvDy
       end if
    elseif(j==nJ)then
       if(NeiLnorth(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1, 1, 0,iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1, 1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev( 0, 1,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev( 0, 1, 1,iBlock)==-1))&
            )then
          jL = j-1; jR = j-2; Ay=-InvDy; By=0.75*InvDy; Cy=0.25*InvDy
       end if
    end if

    if(k==1)then
       if(NeiLbot(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1, 0,-1,iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1, 0,-1,iBlock)==-1)) &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev( 0,-1,-1,iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev( 0, 1,-1,iBlock)==-1)) &
            )then
          kL = k+1; kR = k+2; Az=InvDz; Bz=-0.75*InvDz; Cz=-0.25*InvDz
       end if
    elseif(k==nK)then
       if(NeiLtop(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1, 0, 1,iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1, 0, 1,iBlock)==-1)) &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev( 0,-1, 1,iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev( 0, 1, 1,iBlock)==-1)) &
            )then
          kL = k-1; kR = k-2; Az=-InvDz; Bz=0.75*InvDz; Cz=0.25*InvDz
       end if
    end if

    ! Use central difference to get gradient at face
    if(UseCovariant)then
       call calc_covariant_gradient
    else
       call calc_cartesian_gradient
    end if

  contains

    !==========================================================================

    subroutine calc_cartesian_gradient

      select case(iDir)
      case(x_)
         FaceGrad_D(x_) = InvDx*(Te_G(i,j,k) - Te_G(i-1,j,k))
         FaceGrad_D(y_) = &
              + Ay*(Te_G(i-1,jL,k) + Te_G(i,jL,k)) &
              + By*(Te_G(i-1,j ,k) + Te_G(i,j ,k)) &
              + Cy*(Te_G(i-1,jR,k) + Te_G(i,jR,k))
         FaceGrad_D(z_) = &
              + Az*(Te_G(i-1,j,kL) + Te_G(i,j,kL)) &
              + Bz*(Te_G(i-1,j,k ) + Te_G(i,j,k )) &
              + Cz*(Te_G(i-1,j,kR) + Te_G(i,j,kR))
      case(y_)
         FaceGrad_D(x_) = &
              + Ax*(Te_G(iL,j-1,k) + Te_G(iL,j,k)) &
              + Bx*(Te_G(i ,j-1,k) + Te_G(i ,j,k)) &
              + Cx*(Te_G(iR,j-1,k) + Te_G(iR,j,k))
         FaceGrad_D(y_) = InvDy*(Te_G(i,j,k) - Te_G(i,j-1,k))
         FaceGrad_D(z_) = &
              + Az*(Te_G(i,j-1,kL) + Te_G(i,j,kL)) &
              + Bz*(Te_G(i,j-1,k ) + Te_G(i,j,k )) &
              + Cz*(Te_G(i,j-1,kR) + Te_G(i,j,kR))
      case(z_)
         FaceGrad_D(x_) = &
              + Ax*(Te_G(iL,j,k-1) + Te_G(iL,j,k)) &
              + Bx*(Te_G(i ,j,k-1) + Te_G(i ,j,k)) &
              + Cx*(Te_G(iR,j,k-1) + Te_G(iR,j,k))
         FaceGrad_D(y_) = &
              + Ay*(Te_G(i,jL,k-1) + Te_G(i,jL,k))  &
              + By*(Te_G(i,j ,k-1) + Te_G(i,j ,k))  &
              + Cy*(Te_G(i,jR,k-1) + Te_G(i,jR,k))
         FaceGrad_D(z_) = InvDz*(Te_G(i,j,k) - Te_G(i,j,k-1))
      case default
         write(*,*)'Error in calc_cartesian_gradient: iDir=',iDir
         call stop_mpi('DEBUG')
      end select

    end subroutine calc_cartesian_gradient

    !==========================================================================

    subroutine calc_covariant_gradient

      call stop_mpi('calc_covariant_gradient not yet implemented')

    end subroutine calc_covariant_gradient

  end subroutine calc_face_gradient

    !==========================================================================

  subroutine impl_jac_heat_conduction(iBlock, nVar, Jacobian_VVCI)

    use ModMain, ONLY: nI, nJ, nK, nDim

    integer, parameter :: nStencil = 2*nDim + 1

    integer, intent(in) :: iBlock, nVar
    real, intent(inout) :: Jacobian_VVCI(nVar,nVar,nI,nJ,nK,nStencil)
    !--------------------------------------------------------------------------


  end subroutine impl_jac_heat_conduction

end module ModHeatConduction
