!^CFG COPYRIGHT UM
!============================================================================
module ModGrayDiffusion

  use ModSize,       ONLY: nI, nJ, nK
  use ModVarIndexes, ONLY: p_

  implicit none
  save

  ! This module is needed for the order(u/c) gray diffusion approximation for
  ! radiation-hydrodynamics.

  private !except

  ! Public methods
  public :: get_radiation_energy_flux
  public :: calc_source_gray_diffusion

  ! Index for radiation energy
  integer, parameter, public :: Eradiation_ = p_ - 2

  ! Logical for adding Gray Diffusion
  logical, public :: UseGrayDiffusion        = .false.
  logical, public :: IsNewBlockGrayDiffusion = .true.

  ! Parameters for radiation flux limiter
  logical,           public :: UseRadFluxLimiter  = .false.
  character(len=20), public :: TypeRadFluxLimiter = 'larsen'

  ! Local variables
  real :: Erad_G(-1:nI+2,-1:nJ+2,-1:nK+2)

  real, parameter :: GammaRel = 4.0/3.0

contains

  !==========================================================================

  subroutine get_radiation_energy_flux( &
       iDir, i, j, k, iBlock, State_V, DiffRad, EradFlux_D)

    !\
    ! Calculate the diffusion part of the radiation energy flux.
    !/
    use ModAdvance,    ONLY: State_VGB
    use ModConst,      ONLY: cLightSpeed
    use ModNumConst,   ONLY: cTolerance
    use ModPhysics,    ONLY: Si2No_V, UnitX_, UnitU_
    use ModUser,       ONLY: user_material_properties
    use ModVarIndexes, ONLY: nVar

    integer, intent(in) :: iDir, i, j, k, iBlock
    real,    intent(in) :: State_V(nVar)
    real,    intent(out):: DiffRad
    real,    intent(out):: EradFlux_D(3)

    real :: RosselandMeanOpacitySi
    real :: FaceGrad_D(3), Erad, Grad2ByErad2
    !------------------------------------------------------------------------

    call user_material_properties(State_V, &
         RosselandMeanOpacitySi = RosselandMeanOpacitySi)

    DiffRad = cLightSpeed/(3.0*RosselandMeanOpacitySi) &
         *Si2No_V(UnitU_)*Si2No_V(UnitX_)


    EradFlux_D = 0.0
    call calc_face_gradient(iDir, i, j, k, iBlock, FaceGrad_D)

    if(UseRadFluxLimiter)then

       Grad2ByErad2 = sum(FaceGrad_D**2) &
            /(State_V(Eradiation_)+cTolerance)**2

       select case(TypeRadFluxLimiter)
       case("sum")
          EradFlux_D(iDir) = -FaceGrad_D(iDir) &
               /(1.0/DiffRad+sqrt(Grad2ByErad2))
       case("max")
          EradFlux_D(iDir) = -FaceGrad_D(iDir) &
               /max(1.0/DiffRad,sqrt(Grad2ByErad2))
       case("larsen")
          EradFlux_D(iDir) = -FaceGrad_D(iDir) &
               /sqrt(1.0/DiffRad**2+Grad2ByErad2)
       end select

    else
       EradFlux_D(iDir) = -DiffRad*FaceGrad_D(iDir)
    end if

  end subroutine get_radiation_energy_flux

  !==========================================================================

  subroutine set_block_scalar(Scalar_G, iBlock)

    use ModParallel, ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot, BlkNeighborLev, NOBLK

    integer, intent(in) :: iBlock
    real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2), intent(inout) :: Scalar_G

    real, parameter :: C1 = 8./15., F1 = 2./3., F2 = -1./5.
    real, parameter :: p0 = 5./32., m0 = -3./32., c0 = 15./16.

    integer :: i1, j1, k1, i2, j2, k2, iC, jC, kC, iSide, jSide, kSide
    integer :: iL, iR, jL, jR, kL, kR
    integer :: ip, im, jp, jm, kp, km

    real :: Scalar1_G(-1:nI+2,-1:nJ+2,-1:nK+2) ! temporary array
    real :: Scalarc ! interpolated coarse cell Scalar
    !------------------------------------------------------------------------
    IsNewBlockGrayDiffusion = .false.

    Scalar1_G = Scalar_G

    ! Fix ghost edge and corner ghost cells
    do kSide = -1,1; do jSide = -1,1; do iSide = -1,1
       ! If the corner or edge is not coarser but the face is
       ! then average out the 8 fine cells so that the
       ! general interpolation formulas remain 2nd order
       ! accurate
       if(  BlkNeighborLev(iSide,jSide,kSide,iBlock) /= 1 .and. &
            BlkNeighborLev(iSide,jSide,kSide,iBlock) /= NOBLK .and. ( &
            BlkNeighborLev(iSide,0,0,iBlock) == 1 .or. &
            BlkNeighborLev(0,jSide,0,iBlock) == 1 .or. &
            BlkNeighborLev(0,0,kSide,iBlock) == 1)) then

          if(iSide==0)then
             iL = 1; iR = nI
          elseif(iSide==1)then
             iL = nI+1; iR=nI+2
          else
             iL = -1; iR = 0
          end if
          if(jSide==0)then
             jL = 1; jR = nJ
          elseif(jSide==1)then
             jL = nJ+1; jR=nJ+2
          else
             jL = -1; jR = 0
          end if
          if(kSide==0)then
             kL = 1; kR = nK
          elseif(kSide==1)then
             kL = nK+1; kR=nK+2
          else
             kL = -1; kR = 0
          end if

          do k1=kL,kR,2; do j1=jL,jR,2; do i1=iL,iR,2
             Scalar1_G(i1:i1+1,j1:j1+1,k1:k1+1)= &
                  0.125*sum(Scalar_G(i1:i1+1,j1:j1+1,k1:k1+1))
          end do; end do; end do

       end if
    end do; end do; end do

    ! Do six faces
    if(NeiLeast(iBlock) == 1)then
       do k1=1, nK, 2; do j1=1, nJ, 2; do k2 = k1,k1+1; do j2 = j1,j1+1
          jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
          kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
          Scalarc = c0*Scalar1_G(0,j2,k2) &
               +    p0*Scalar1_G(0,jp,kp) &
               +    m0*Scalar1_G(0,jm,km)
          Scalar_G(0,j2,k2) = &
               C1*Scalarc + F1*Scalar_G(1,j2,k2) + F2*Scalar_G(2,j2,k2)
       end do; end do; end do; end do
    end if

    if(NeiLwest(iBlock) == 1)then
       do k1=1, nK, 2; do j1=1, nJ, 2; do k2 = k1,k1+1; do j2 = j1,j1+1
          jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
          kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
          Scalarc = c0*Scalar1_G(nI+1,j2,k2) &
               +    p0*Scalar1_G(nI+1,jp,kp) &
               +    m0*Scalar1_G(nI+1,jm,km)
          Scalar_G(nI+1,j2,k2)= &
               C1*Scalarc + F1*Scalar_G(nI,j2,k2) + F2*Scalar_G(nI-1,j2,k2)
       end do; end do; end do; end do
    end if

    if(NeiLsouth(iBlock) == 1)then
       do k1=1, nK, 2; do i1=1, nI, 2; do k2 = k1,k1+1; do i2 = i1,i1+1
          ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
          kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
          Scalarc = c0*Scalar1_G(i2,0,k2) &
               +    p0*Scalar1_G(ip,0,kp) &
               +    m0*Scalar1_G(im,0,km)
          Scalar_G(i2,0,k2) = &
               C1*Scalarc + F1*Scalar_G(i2,1,k2) + F2*Scalar_G(i2,2,k2)
       end do; end do; end do; end do
    end if

    if(NeiLnorth(iBlock) == 1)then
       do k1=1, nK, 2; do i1=1, nI, 2; do k2 = k1,k1+1; do i2 = i1,i1+1
          ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
          kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
          Scalarc = c0*Scalar1_G(i2,nJ+1,k2) &
               +    p0*Scalar1_G(ip,nJ+1,kp) &
               +    m0*Scalar1_G(im,nJ+1,km)
          Scalar_G(i2,nJ+1,k2) = &
               C1*Scalarc + F1*Scalar_G(i2,nJ,k2) + F2*Scalar_G(i2,nJ-1,k2)
       end do; end do; end do; end do
    end if

    if(NeiLbot(iBlock) == 1)then
       do j1=1, nJ, 2; do i1=1, nI, 2; do j2 = j1,j1+1; do i2 = i1,i1+1
          ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
          jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
          Scalarc = c0*Scalar1_G(i2,j2,0) &
               +    p0*Scalar1_G(ip,jp,0) &
               +    m0*Scalar1_G(im,jm,0)
          Scalar_G(i2,j2,0) = &
               C1*Scalarc + F1*Scalar_G(i2,j2,1) + F2*Scalar_G(i2,j2,2)
       end do; end do; end do; end do
    end if

    if(NeiLtop(iBlock) == 1)then
       do j1=1, nJ, 2; do i1=1, nI, 2; do j2 = j1,j1+1; do i2 = i1,i1+1
          ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
          jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
          Scalarc = c0*Scalar1_G(i2,j2,nK+1) &
               +    p0*Scalar1_G(ip,jp,nK+1) &
               +    m0*Scalar1_G(im,jm,nK+1)
          Scalar_G(i2,j2,nK+1) = &
               C1*Scalarc + F1*Scalar_G(i2,j2,nK) + F2*Scalar_G(i2,j2,nK-1)
       end do; end do; end do; end do
    end if

    ! Do 12 edges
    ! 4 X edges
    do kSide = -1,1,2; do jSide = -1,1,2
       if(  BlkNeighborLev(0, jSide, kSide, iBlock) /= 1 .and. .not. ( &
            BlkNeighborLev(0, jSide, kSide, iBlock) == NOBLK .and. ( &
            BlkNeighborLev(0, jSide, 0, iBlock) == 1 .or. &
            BlkNeighborLev(0, 0, kSide, iBlock) == 1))) CYCLE

       j1=1; if(jSide==1) j1=nJ; j2 = j1-jSide; jC = j1+jSide
       k1=1; if(kSide==1) k1=nK; k2 = k1-kSide; kC = k1+kSide
       do i1 = 1,nI,2; do i2 = i1, i1+1
          ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
          Scalarc = c0*Scalar1_G(i2,jC,kC) &
               +    p0*Scalar1_G(ip,jC,kC) &
               +    m0*Scalar1_G(im,jC,kC)
          Scalar_G(i2,jC,kC) = &
               C1*Scalarc + F1*Scalar_G(i2,j1,k1) + F2*Scalar_G(i2,j2,k2)
       end do; end do
    end do; end do
    ! 4 Y edges
    do kSide = -1,1,2; do iSide = -1,1,2
       if(  BlkNeighborLev(iSide, 0, kSide, iBlock) /= 1 .and. .not. ( &
            BlkNeighborLev(iSide, 0, kSide, iBlock) == NOBLK .and. ( &
            BlkNeighborLev(iSide, 0, 0, iBlock) == 1 .or. &
            BlkNeighborLev(0, 0, kSide, iBlock) == 1))) CYCLE

       i1=1; if(iSide==1) i1=nI; i2 = i1-iSide; iC = i1+iSide
       k1=1; if(kSide==1) k1=nK; k2 = k1-kSide; kC = k1+kSide
       do j1 = 1, nJ, 2; do j2 = j1, j1+1
          jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
          Scalarc = c0*Scalar1_G(iC,j2,kC) &
               +    p0*Scalar1_G(iC,jp,kC) &
               +    m0*Scalar1_G(iC,jm,kC)
          Scalar_G(iC,j2,kC) = &
               C1*Scalarc + F1*Scalar_G(i1,j2,k1) + F2*Scalar_G(i2,j2,k2)
       end do; end do
    end do; end do
    ! 4 Z edges
    do jSide = -1,1,2; do iSide = -1,1,2
       if(  BlkNeighborLev(iSide, jSide, 0, iBlock) /= 1.and. .not. ( &
            BlkNeighborLev(iSide, jSide, 0, iBlock) == NOBLK .and. ( &
            BlkNeighborLev(iSide, 0, 0, iBlock) == 1 .or. &
            BlkNeighborLev(0, jSide, 0, iBlock) == 1))) CYCLE

       i1=1; if(iSide==1) i1=nI; i2 = i1-iSide; iC = i1+iSide
       j1=1; if(jSide==1) j1=nJ; j2 = j1-jSide; jC = j1+jSide
       do k1 = 1, nK, 2 ; do k2 = k1, k1 + 1
          kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
          Scalarc = c0*Scalar1_G(iC,jC,k2) &
               +    p0*Scalar1_G(iC,jC,kp) &
               +    m0*Scalar1_G(iC,jC,km)
          Scalar_G(iC,jC,k2) = &
               C1*Scalarc + F1*Scalar_G(i1,j1,k2) + F2*Scalar_G(i2,j2,k2)
       end do; end do         
    end do; end do

  end subroutine set_block_scalar

  !==========================================================================

  subroutine calc_face_gradient(iDir, i, j, k, iBlock, FaceGrad_D) 

    use ModAdvance,   ONLY: State_VGB
    use ModCovariant, ONLY: UseCovariant
    use ModGeometry,  ONLY: Dx_Blk, Dy_Blk, Dz_Blk
    use ModMain,      ONLY: x_, y_, z_
    use ModParallel,  ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot, BlkNeighborLev

    integer, intent(in) :: iDir, i, j, k, iBlock
    real, intent(out) :: FaceGrad_D(3)

    integer :: iL, iR, jL, jR, kL, kR
    real :: Ax, Bx, Cx, Ay, By, Cy, Az, Bz, Cz

    Real :: InvDx, InvDy, InvDz
    !------------------------------------------------------------------------
    InvDx = 1.0/Dx_Blk(iBlock)
    InvDy = 1.0/Dy_Blk(iBlock)
    InvDz = 1.0/Dz_Blk(iBlock)

    if(IsNewBlockGrayDiffusion)then
       Erad_G = State_VGB(Eradiation_,:,:,:,iBlock)
       call set_block_scalar(Erad_G, iBlock)
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
    !========================================================================

    subroutine calc_cartesian_gradient

      select case(iDir)
      case(x_)
         FaceGrad_D(x_) = InvDx*(Erad_G(i,j,k) - Erad_G(i-1,j,k))
         FaceGrad_D(y_) = &
              + Ay*(Erad_G(i-1,jL,k) + Erad_G(i,jL,k)) &
              + By*(Erad_G(i-1,j ,k) + Erad_G(i,j ,k)) &
              + Cy*(Erad_G(i-1,jR,k) + Erad_G(i,jR,k))
         FaceGrad_D(z_) = &
              + Az*(Erad_G(i-1,j,kL) + Erad_G(i,j,kL)) &
              + Bz*(Erad_G(i-1,j,k ) + Erad_G(i,j,k )) &
              + Cz*(Erad_G(i-1,j,kR) + Erad_G(i,j,kR))
      case(y_)
         FaceGrad_D(x_) = &
              + Ax*(Erad_G(iL,j-1,k) + Erad_G(iL,j,k)) &
              + Bx*(Erad_G(i ,j-1,k) + Erad_G(i ,j,k)) &
              + Cx*(Erad_G(iR,j-1,k) + Erad_G(iR,j,k))
         FaceGrad_D(y_) = InvDy*(Erad_G(i,j,k) - Erad_G(i,j-1,k))
         FaceGrad_D(z_) = &
              + Az*(Erad_G(i,j-1,kL) + Erad_G(i,j,kL)) &
              + Bz*(Erad_G(i,j-1,k ) + Erad_G(i,j,k )) &
              + Cz*(Erad_G(i,j-1,kR) + Erad_G(i,j,kR))
      case(z_)
         FaceGrad_D(x_) = &
              + Ax*(Erad_G(iL,j,k-1) + Erad_G(iL,j,k)) &
              + Bx*(Erad_G(i ,j,k-1) + Erad_G(i ,j,k)) &
              + Cx*(Erad_G(iR,j,k-1) + Erad_G(iR,j,k))
         FaceGrad_D(y_) = &
              + Ay*(Erad_G(i,jL,k-1) + Erad_G(i,jL,k))  &
              + By*(Erad_G(i,j ,k-1) + Erad_G(i,j ,k))  &
              + Cy*(Erad_G(i,jR,k-1) + Erad_G(i,jR,k))
         FaceGrad_D(z_) = InvDz*(Erad_G(i,j,k) - Erad_G(i,j,k-1))
      case default
         write(*,*)'Error in calc_cartesian_gradient: iDir=',iDir
         call stop_mpi('DEBUG')
      end select

    end subroutine calc_cartesian_gradient

    !========================================================================

    subroutine calc_covariant_gradient

      call stop_mpi('calc_covariant_gradient not yet implemented')

    end subroutine calc_covariant_gradient

  end subroutine calc_face_gradient

  !==========================================================================

  subroutine calc_source_gray_diffusion(iBlock)

    use ModAdvance,    ONLY: State_VGB, Source_VC, &
         uDotArea_XI, uDotArea_YI, uDotArea_ZI
    use ModConst,      ONLY: cRadiation, cLightSpeed
    use ModGeometry,   ONLY: vInv_CB
    use ModPhysics,    ONLY: Si2No_V, UnitEnergyDens_, UnitT_
    use ModUser,       ONLY: user_material_properties
    use ModVarIndexes, ONLY: Energy_, NameVar_V

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: TeSi, AbsorptionOpacitySi, AbsorptionOpacity
    real :: RadCompression, AbsorptionEmission
    character(len=*), parameter:: NameSub = "calc_source_gray_diffusion"
    !------------------------------------------------------------------------

    ! Make sure that Eradiation_ is correct
    if(NameVar_V(Eradiation_) /= "Erad") call stop_mpi(NameSub// &
         ": incorrect index for Erad variable in ModEquation")

    do k=1,nK; do j=1,nJ; do i=1,nI

       call user_material_properties(State_VGB(1,i,j,k,iBlock), &
            TeSi = TeSi, &
            AbsorptionOpacitySi = AbsorptionOpacitySi)

       AbsorptionOpacity = AbsorptionOpacitySi*cLightSpeed &
            /Si2No_V(UnitT_)

       ! Adiabatic compression of radiation energy by fluid velocity (fluid 1)
       ! (GammaRel-1)*Erad*Div(U)
       RadCompression = (GammaRel-1.0)*State_VGB(Eradiation_,i,j,k,iBlock) &
            *vInv_CB(i,j,k,iBlock)&
            *(uDotArea_XI(i+1,j,k,1) - uDotArea_XI(i,j,k,1) &
             +uDotArea_YI(i,j+1,k,1) - uDotArea_YI(i,j,k,1) &
             +uDotArea_ZI(i,j,k+1,1) - uDotArea_ZI(i,j,k,1))

       ! Source term due to absorption and emission
       ! Sigma_a*(cRadiation*Te**4-Erad)
       AbsorptionEmission = AbsorptionOpacity &
            *(cRadiation*TeSi**4*Si2No_V(UnitEnergyDens_) &
            - State_VGB(Eradiation_,i,j,k,iBlock))

       ! dErad/dt = - adiabatic compression + AbsorptionEmission
       Source_VC(Eradiation_,i,j,k) = Source_VC(Eradiation_,i,j,k) &
            - RadCompression + AbsorptionEmission

       ! dE/dt    = + adiabatic compression - AbsorptionEmission
       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
            + RadCompression - AbsorptionEmission

    end do; end do; end do

  end subroutine calc_source_gray_diffusion

end module ModGrayDiffusion
