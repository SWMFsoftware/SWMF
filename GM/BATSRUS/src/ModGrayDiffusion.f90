!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
!============================================================================
module ModGrayDiffusion

  use ModVarIndexes, ONLY: p_
  use ModProcMH, ONLY: iProc
  !use ModMain, ONLY: iTest, jTest, kTest, BlkTest


  implicit none
  save

  ! This module is needed for the order(u/c) gray diffusion approximation for
  ! radiation-hydrodynamics.

  private !except

  ! Public methods
  public :: init_gray_diffusion
  public :: get_radiation_energy_flux
  public :: calc_source_gray_diffusion
  public :: get_impl_gray_diff_state
  public :: get_gray_diffusion_bc
  public :: get_gray_diffusion_rhs
  public :: get_gray_diff_jacobian
  public :: update_impl_gray_diff

  ! Logical for adding Gray Diffusion
  logical, public :: IsNewBlockGrayDiffusion = .true.
  logical, public :: IsNewTimestepGrayDiffusion = .true.

  ! Coefficients for two-temperature electron-radiation model
  real, allocatable, public :: DiffCoef_VFDB(:,:,:,:,:,:)
  real, allocatable         :: RelaxationCoef_CB(:,:,:,:)
  real, allocatable         :: RelaxSemiCoef_VCB(:,:,:,:,:)
  real, allocatable         :: DiffSemiCoef_VGB(:,:,:,:,:)
  real, allocatable         :: PointSemiCoef_VCB(:,:,:,:,:)

  ! Index which vars involve diffusion
  integer, allocatable :: iDiff_I(:)
  integer :: nDiff

  ! Number of relaxation coefficients
  integer, allocatable :: iRelax_I(:)
  integer :: nRelax

  ! Number of point implicit coefficients
  integer :: nPoint
  ! Index which variable involve point implicit calculation
  integer :: iPoint
  ! Named indices for point implicit coefficients
  integer, parameter :: Relax_ = 1, Planck_ = 2

  ! radiation energy used for calculating radiative energy flux
  real, allocatable :: Erad_G(:,:,:)

  real, parameter :: GammaRel = 4.0/3.0

  ! local index for extra internal energy to keep the compiler happy
  integer, parameter :: EintExtra_ = p_ - 1

  ! Named indices for semi-implicit variables
  integer :: TeImpl_, EradImpl_

  real :: EradMin

contains

  !============================================================================

  subroutine init_gray_diffusion

    use ModAdvance,     ONLY: Eradiation_
    use ModMain,        ONLY: UseGrayDiffusion
    use ModSize,        ONLY: nI, nJ, nK, MaxBlock, nDim
    use ModVarIndexes,  ONLY: NameVar_V
    use ModImplicit,    ONLY: UseSemiImplicit, UseFullImplicit, &
         TypeSemiImplicit
    use ModPhysics,     ONLY: Si2No_V, UnitTemperature_, cRadiationNo
    use ModTemperature, ONLY: TradMinSi

    character(len=*), parameter :: NameSub = "init_gray_diffusion"
    !------------------------------------------------------------------------

    if(allocated(Erad_G)) RETURN

    ! Make sure that Eradiation_ is correct
    if(UseGrayDiffusion)then
       if(NameVar_V(Eradiation_) /= "Erad") call stop_mpi(NameSub// &
            ": incorrect index for Erad variable in ModEquation")

       EradMin = cRadiationNo*(TradMinSi*Si2No_V(UnitTemperature_))**4
    end if

    allocate(Erad_G(-1:nI+2,-1:nJ+2,-1:nK+2))

    if(UseFullImplicit)then
       nDiff = 1
       allocate(RelaxationCoef_CB(1:nI,1:nJ,1:nK,MaxBlock))
    end if
       
    if(UseSemiImplicit)then

       select case(TypeSemiImplicit)
       case('radiation')
          EradImpl_ = 1
          nDiff = 1
          allocate(iDiff_I(nDiff))
          iDiff_I(1) = EradImpl_
          nRelax = 0
          nPoint = 2
          iPoint = EradImpl_
       case('radcond')
          TeImpl_ = 1; EradImpl_ = 2
          nDiff = 2
          allocate(iDiff_I(nDiff))
          iDiff_I = (/ TeImpl_, EradImpl_ /)
          nRelax = 1
          allocate(iRelax_I(nRelax))
          iRelax_I(1) = EradImpl_
          nPoint = 0
       case('cond')
          TeImpl_ = 1
          nDiff = 1
          allocate(iDiff_I(nDiff))
          iDiff_I(1) = TeImpl_
          nRelax = 0
          nPoint = 0
       end select

       if(nRelax>0)then
          allocate(RelaxSemiCoef_VCB(nRelax,nI,nJ,nK,MaxBlock))
          RelaxSemiCoef_VCB = 0.0
       end if

       if(nPoint>0)then
          allocate(PointSemiCoef_VCB(nPoint,nI,nJ,nK,MaxBlock))
          PointSemiCoef_VCB = 0.0
       end if

       allocate(DiffSemiCoef_VGB(nDiff,-1:nI+2,-1:nJ+2,-1:nK+2,MaxBlock))
       DiffSemiCoef_VGB = 0.0
    end if

    allocate(DiffCoef_VFDB(nDiff,1:nI+1,1:nJ+1,1:nK+1,nDim,MaxBlock))
    DiffCoef_VFDB = 0.0 ! make sure all elements are initialized

  end subroutine init_gray_diffusion

  !==========================================================================

  subroutine get_radiation_energy_flux( &
       iDir, i, j, k, iBlock, State_V, EradFlux_D)

    !\
    ! Calculate the diffusion part of the radiation energy flux.
    !/
    use ModAdvance,    ONLY: Eradiation_
    use ModPhysics,    ONLY: Si2No_V, UnitX_, Clight
    use ModTemperature,ONLY: UseRadFluxLimiter, TypeRadFluxLimiter
    use ModUser,       ONLY: user_material_properties
    use ModVarIndexes, ONLY: nVar

    integer, intent(in) :: iDir, i, j, k, iBlock
    real,    intent(in) :: State_V(nVar)
    real,    intent(out):: EradFlux_D(3)

    real :: RosselandMeanOpacitySi, DiffRad
    real :: RosselandMeanOpacity
    real :: FaceGrad_D(3), Grad2ByErad2
    !------------------------------------------------------------------------

    call calc_face_gradient(iDir, i, j, k, iBlock, FaceGrad_D)

    if(IsNewTimestepGrayDiffusion)then

       call user_material_properties(State_V, &
            RosselandMeanOpacitySiOut = RosselandMeanOpacitySi)

       RosselandMeanOpacity = RosselandMeanOpacitySi/Si2No_V(UnitX_)

       if(UseRadFluxLimiter)then
          Grad2ByErad2 = sum(FaceGrad_D**2)/State_V(Eradiation_)**2

          select case(TypeRadFluxLimiter)
          case("sum")
             DiffRad = Clight/(3*RosselandMeanOpacity + sqrt(Grad2ByErad2))
          case("max")
             DiffRad = Clight/max(3*RosselandMeanOpacity,sqrt(Grad2ByErad2))
          case("larsen")
             DiffRad = Clight/sqrt(9*RosselandMeanOpacity**2 + Grad2ByErad2)
          end select
       else
          DiffRad = Clight/(3*RosselandMeanOpacity)
       end if

       DiffCoef_VFDB(1,i,j,k,iDir,iBlock) = DiffRad
    else
       DiffRad = DiffCoef_VFDB(1,i,j,k,iDir,iBlock)
    end if

    EradFlux_D = 0.0
    EradFlux_D(iDir) = -DiffRad*FaceGrad_D(iDir)

  end subroutine get_radiation_energy_flux

  !==========================================================================

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
    !------------------------------------------------------------------------
    IsNewBlockGrayDiffusion = .false.

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

  !==========================================================================

  subroutine calc_face_gradient(iDir, i, j, k, iBlock, FaceGrad_D) 

    use ModAdvance,   ONLY: State_VGB, Eradiation_
    use ModCovariant, ONLY: UseCovariant
    use ModGeometry,  ONLY: Dx_Blk, Dy_Blk, Dz_Blk
    use ModMain,      ONLY: x_, y_, z_
    use ModParallel,  ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot, BlkNeighborLev
    use ModSize,      ONLY: nI, nJ, nK

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
         uDotArea_XI, uDotArea_YI, uDotArea_ZI, Eradiation_
    use ModConst,      ONLY: cLightSpeed
    use ModGeometry,   ONLY: vInv_CB, y_BLK, TypeGeometry
    use ModImplicit,   ONLY: UseFullImplicit
    use ModPhysics,    ONLY: cRadiationNo, Si2No_V, UnitTemperature_, UnitT_
    use ModMain,       ONLY: nI, nJ, nK, UseGrayDiffusion
    use ModUser,       ONLY: user_material_properties
    use ModVarIndexes, ONLY: Energy_, RhoUy_

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: TeSi, Te, DivU
    real :: RadCompression, AbsorptionEmission, PlanckOpacitySi
    character(len=*), parameter:: NameSub = "calc_source_gray_diffusion"
    !------------------------------------------------------------------------

    do k=1,nK; do j=1,nJ; do i=1,nI

       DivU = vInv_CB(i,j,k,iBlock)* &
            ( uDotArea_XI(i+1,j,k,1) - uDotArea_XI(i,j,k,1) &
            + uDotArea_YI(i,j+1,k,1) - uDotArea_YI(i,j,k,1) &
            + uDotArea_ZI(i,j,k+1,1) - uDotArea_ZI(i,j,k,1) )

       ! Adiabatic compression of radiation energy by fluid velocity (fluid 1)
       ! (GammaRel-1)*Erad*Div(U)
       RadCompression = (GammaRel-1.0)*State_VGB(Eradiation_,i,j,k,iBlock)*DivU

       ! dErad/dt = - adiabatic compression
       Source_VC(Eradiation_,i,j,k) = Source_VC(Eradiation_,i,j,k) &
            - RadCompression

       ! dE/dt    = + adiabatic compression
       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + RadCompression

       if(.not.UseFullImplicit) CYCLE

       if(IsNewTimestepGrayDiffusion)then
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               AbsorptionOpacitySiOut = PlanckOpacitySi)

          RelaxationCoef_CB(i,j,k,iBlock) = &
               PlanckOpacitySi*cLightSpeed/Si2No_V(UnitT_)
       end if

       call user_material_properties(State_VGB(:,i,j,k,iBlock), &
            TeSiOut = TeSi)

       Te = TeSi*Si2No_V(UnitTemperature_)

       ! Source term due to absorption and emission
       ! Sigma_a*(cRadiation*Te**4-Erad)
       AbsorptionEmission =  RelaxationCoef_CB(i,j,k,iBlock) &
            *(cRadiationNo*Te**4 - State_VGB(Eradiation_,i,j,k,iBlock))

       ! dErad/dt = + AbsorptionEmission
       Source_VC(Eradiation_,i,j,k) = Source_VC(Eradiation_,i,j,k) &
            + AbsorptionEmission

       ! dE/dt    = - AbsorptionEmission
       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
            - AbsorptionEmission

    end do; end do; end do

    if(TypeGeometry=='rz' .and. UseGrayDiffusion)then
       ! Add "geometrical source term" p/r to the radial momentum equation
       ! The "radial" direction is along the Y axis
       ! NOTE: here we have to use signed radial distance!
       do k=1,nK; do j=1, nJ; do i=1, nI
          Source_VC(RhoUy_,i,j,k) = Source_VC(RhoUy_,i,j,k) &
               + (1./3.)*State_VGB(Eradiation_,i,j,k,iBlock) &
               / y_BLK(i,j,k,iBlock)
       end do; end do; end do
    end if

  end subroutine calc_source_gray_diffusion

  !============================================================================
  ! Semi-implicit interface for Erad and Eint
  !============================================================================

  subroutine get_impl_gray_diff_state(StateImpl_VGB,DconsDsemi_VCB)

    use ModAdvance,  ONLY: Eradiation_, State_VGB
    use ModImplicit, ONLY: nw, nImplBlk, impl2iBlk, kr, TypeSemiImplicit, &
         ImplCoeff
    use ModMain,     ONLY: nDim, x_, y_, nI, nJ, nK, MaxImplBlk, Dt
    use ModPhysics,  ONLY: inv_gm1, Clight, cRadiationNo, &
         Si2No_V, UnitTemperature_, UnitEnergyDens_, UnitX_, UnitU_
    use ModUser,     ONLY: user_material_properties
    use ModTemperature, ONLY: UseRadFluxLimiter, TypeRadFluxLimiter
    use ModGeometry, ONLY: dx_BLK, dy_BLK, dz_BLK, vInv_CB, &
         UseCovariant, TypeGeometry, FaceAreaI_DFB, FaceAreaJ_DFB
    use ModParallel, ONLY: NOBLK, NeiLev
    use ModMessagePass, ONLY: message_pass_dir

    real, intent(out) :: StateImpl_VGB(nw,0:nI+1,0:nJ+1,0:nK+1,MaxImplBlk)
    real, intent(inout) :: DconsDsemi_VCB(nw,nI,nJ,nK,MaxImplBlk)

    integer :: iImplBlock, iBlock, i, j, k
    real :: PlanckOpacitySi, PlanckOpacity, CvSi, Cv, TeSi, Te
    real :: RosselandMeanOpacitySi, RosselandMeanOpacity
    real :: HeatConductionCoefSi, HeatConductionCoef
    real :: Grad2ByErad2, DiffRad, InvDx2, InvDy2, InvDz2
    real :: InvDx, InvDy

    integer :: iDim, Di, Dj, Dk, iDiff, nDimInUse
    real :: Coeff, Dxyz_D(3)

    character(len=*), parameter:: NameSub='get_impl_gray_diff_state'
    !--------------------------------------------------------------------------

    do iImplBlock = 1, nImplBLK

       iBlock = impl2iBLK(iImplBlock)
       IsNewBlockGrayDiffusion = .true.

       if(TypeSemiImplicit=='radiation' .or. TypeSemiImplicit=='radcond')then
          do k = 0, nK+1; do j = 0, nJ+1; do i = 0, nI+1
             StateImpl_VGB(EradImpl_,i,j,k,iImplBlock) = &
                  State_VGB(Eradiation_,i,j,k,iBlock)
          end do; end do; end do
       end if
       if(TypeSemiImplicit=='radcond')then
          do k = 0, nK+1; do j = 0, nJ+1; do i = 0, nI+1
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  TeSiOut = TeSi)
             Te = TeSi*Si2No_V(UnitTemperature_)
             StateImpl_VGB(TeImpl_,i,j,k,iImplBlock) = cRadiationNo*Te**4
          end do; end do; end do
       elseif(TypeSemiImplicit=='cond')then
          do k = 0, nK+1; do j = 0, nJ+1; do i = 0, nI+1
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  TeSiOut = TeSi)
             StateImpl_VGB(TeImpl_,i,j,k,iImplBlock) = &
                  TeSi*Si2No_V(UnitTemperature_)
          end do; end do; end do
       end if

       InvDx2 = 0.5/dx_BLK(iBlock)
       InvDy2 = 0.5/dy_BLK(iBlock)
       InvDz2 = 0.5/dz_BLK(iBlock)

       ! calculate coefficients for linearized energy exchange and diffusion
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               AbsorptionOpacitySiOut = PlanckOpacitySi, &
               RosselandMeanOpacitySiOut = RosselandMeanOpacitySi, &
               CvSiOut = CvSi, TeSiOut = TeSi, &
               HeatConductionCoefSiOut = HeatConductionCoefSi)

          PlanckOpacity = PlanckOpacitySi/Si2No_V(UnitX_)
          Cv = CvSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
          Te = TeSi*Si2No_V(UnitTemperature_)

          select case(TypeSemiImplicit)
          case('radiation')
             ! This coefficient is cR'' = cR/(1+dt*cR*dPlanck/dEint)
             PointSemiCoef_VCB(Relax_,i,j,k,iBlock) = Clight*PlanckOpacity  &
                  /(1 + ImplCoeff*Dt*Clight*PlanckOpacity &
                  *4.0*cRadiationNo*Te**3 / Cv)

             ! This is just the Planck function at time level * saved
             PointSemiCoef_VCB(Planck_,i,j,k,iBlock) = cRadiationNo*Te**4
          case('radcond')
             RelaxSemiCoef_VCB(1,i,j,k,iBlock) = Clight*PlanckOpacity

             DconsDsemi_VCB(TeImpl_,i,j,k,iImplBlock) = &
                  Cv/(4.0*cRadiationNo*Te**3)
          case('cond')
             DconsDsemi_VCB(TeImpl_,i,j,k,iImplBlock) = Cv
          end select

          call get_diffusion_coef

       end do; end do; end do

       if(NeiLev(1,iBlock) == NOBLK)then
          i = 0
          do k = 1, nK; do j = 1, nJ
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  RosselandMeanOpacitySiOut = RosselandMeanOpacitySi, &
                  HeatConductionCoefSiOut = HeatConductionCoefSi, &
                  TeSiOut = TeSi)
             call get_diffusion_coef
          end do; end do
       end if
       if(NeiLev(2,iBlock) == NOBLK)then
          i = nI + 1
          do k = 1, nK; do j = 1, nJ
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  RosselandMeanOpacitySiOut = RosselandMeanOpacitySi, &
                  HeatConductionCoefSiOut = HeatConductionCoefSi, &
                  TeSiOut = TeSi)
             call get_diffusion_coef
          end do; end do
       end if
       if(NeiLev(3,iBlock) == NOBLK)then
          j = 0
          do k = 1, nK; do i = 1, nI
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  RosselandMeanOpacitySiOut = RosselandMeanOpacitySi, &
                  HeatConductionCoefSiOut = HeatConductionCoefSi, &
                  TeSiOut = TeSi)
             call get_diffusion_coef
          end do; end do
       end if
       if(NeiLev(4,iBlock) == NOBLK)then
          j = nJ + 1
          do k = 1, nK; do i = 1, nI
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  RosselandMeanOpacitySiOut = RosselandMeanOpacitySi, &
                  HeatConductionCoefSiOut = HeatConductionCoefSi, &
                  TeSiOut = TeSi)
             call get_diffusion_coef
          end do; end do
       end if
       if(NeiLev(5,iBlock) == NOBLK)then
          k = 0
          do j = 1, nJ; do i = 1, nI
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  RosselandMeanOpacitySiOut = RosselandMeanOpacitySi, &
                  HeatConductionCoefSiOut = HeatConductionCoefSi, &
                  TeSiOut = TeSi)
             call get_diffusion_coef
          end do; end do
       end if
       if(NeiLev(6,iBlock) == NOBLK)then
          k = nK + 1
          do j = 1, nJ; do i = 1, nI
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  RosselandMeanOpacitySiOut = RosselandMeanOpacitySi, &
                  HeatConductionCoefSiOut = HeatConductionCoefSi, &
                  TeSiOut = TeSi)
             call get_diffusion_coef
          end do; end do
       end if

    end do

    nDimInUse = 3; if(TypeGeometry == 'rz') nDimInUse = 2

    ! Message pass to fill in ghost cells 
    call message_pass_dir(iDirMin=1,iDirMax=3,Width=1, &
         SendCorners=.false.,ProlongOrder=1,nVar=nDiff, &
         Sol_VGB=DiffSemiCoef_VGB,restrictface=.true.)

    do iImplBlock = 1, nImplBLK
       iBlock = impl2iBLK(iImplBlock)
       ! Calculate face averaged values. Include geometric factors.

       call face_equal(1,2,nI,1,nJ,1,nK)
       if(NeiLev(1,iBlock)==0.or.NeiLev(1,iBlock)==NOBLK)then
          call face_equal(1,1,1,1,nJ,1,nK)
       else if(NeiLev(1,iBlock)==-1)then
          call face_left_coarse2fine(1,1,1,1,nJ,1,nK)
       else if(NeiLev(1,iBlock)==1)then
          call face_left_fine2coarse(1,1,1,1,nJ,1,nK)
       end if
       if(NeiLev(2,iBlock)==0.or.NeiLev(2,iBlock)==NOBLK)then
          call face_equal(1,nI+1,nI+1,1,nJ,1,nK)
       else if(NeiLev(2,iBlock)==-1)then
          call face_right_coarse2fine(1,nI+1,nI+1,1,nJ,1,nK)
       else if(NeiLev(2,iBlock)==1)then
          call face_right_fine2coarse(1,nI+1,nI+1,1,nJ,1,nK)
       end if

       call face_equal(2,1,nI,2,nJ,1,nK)
       if(NeiLev(3,iBlock)==0.or.NeiLev(3,iBlock)==NOBLK)then
          call face_equal(2,1,nI,1,1,1,nK)
       else if(NeiLev(3,iBlock)==-1)then
          call face_left_coarse2fine(2,1,nI,1,1,1,nK)
       else if(NeiLev(3,iBlock)==1)then
          call face_left_fine2coarse(2,1,nI,1,1,1,nK)
       end if
       if(NeiLev(4,iBlock)==0.or.NeiLev(4,iBlock)==NOBLK)then
          call face_equal(2,1,nI,nJ+1,nJ+1,1,nK)
       else if(NeiLev(4,iBlock)==-1)then
          call face_right_coarse2fine(2,1,nI,nJ+1,nJ+1,1,nK)
       else if(NeiLev(4,iBlock)==1)then
          call face_right_fine2coarse(2,1,nI,nJ+1,nJ+1,1,nK)
       end if

       if(nDimInUse==3)then
          call face_equal(3,1,nI,1,nJ,2,nK)
          if(NeiLev(5,iBlock)==0.or.NeiLev(5,iBlock)==NOBLK)then
             call face_equal(3,1,nI,1,nJ,1,1)
          else if(NeiLev(5,iBlock)==-1)then
             call face_left_coarse2fine(3,1,nI,1,nJ,1,1)
          else if(NeiLev(5,iBlock)==1)then
             call face_left_fine2coarse(3,1,nI,1,nJ,1,1)
          end if
          if(NeiLev(6,iBlock)==0.or.NeiLev(6,iBlock)==NOBLK)then
             call face_equal(3,1,nI,1,nJ,nK+1,nK+1)
          else if(NeiLev(6,iBlock)==-1)then
             call face_right_coarse2fine(3,1,nI,1,nJ,nK+1,nK+1)
          else if(NeiLev(6,iBlock)==1)then
             call face_right_fine2coarse(3,1,nI,1,nJ,nK+1,nK+1)
          end if
       end if

       if(.not.UseCovariant)then
          Dxyz_D = (/dx_BLK(iBlock), dy_BLK(iBlock), dz_Blk(iBlock)/)
          do iDim = 1, nDim
             ! FaceYZ/dx = Volume/dx^2
             Coeff = 1.0 / (Dxyz_D(iDim)**2 * vInv_CB(1,1,1,iBlock))
             Di = kr(iDim,1); Dj = kr(iDim,2); Dk = kr(iDim,3)
             do k=1,nK+Dk; do j=1,nJ+Dj; do i=1,nI+Di
                do iDiff = 1, nDiff
                   DiffCoef_VFDB(iDiff,i,j,k,iDim,iBlock) = &
                        Coeff*DiffCoef_VFDB(iDiff,i,j,k,iDim,iBlock)
                end do
             enddo; enddo; enddo
          end do

       elseif(TypeGeometry == 'rz')then

          InvDx = 1.0/Dx_Blk(iBlock)
          do k=1,nK; do j=1,nJ; do i=1,nI+1
             do iDiff = 1, nDiff
                DiffCoef_VFDB(iDiff,i,j,k,x_,iBlock) = &
                     InvDx*FaceAreaI_DFB(x_,i,j,k,iBlock) &
                     *DiffCoef_VFDB(iDiff,i,j,k,x_,iBlock)
             end do
          end do; end do; end do

          InvDy = 1.0/Dy_Blk(iBlock)
          do k=1,nK; do j=1,nJ+1; do i=1,nI
             do iDiff = 1, nDiff
                DiffCoef_VFDB(iDiff,i,j,k,y_,iBlock) = &
                     InvDy*FaceAreaJ_DFB(y_,i,j,k,iBlock) &
                     *DiffCoef_VFDB(iDiff,i,j,k,y_,iBlock)
             end do
          end do; end do; end do
       else
          call stop_mpi(NameSub//': unimplemented TypeGeometry=//TypeGeometry')
       end if
    end do

  contains

    subroutine face_equal(iDim,iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iDim, iMin, iMax, jMin, jMax, kMin, kMax

      !------------------------------------------------------------------------

      Di = kr(iDim,1); Dj = kr(iDim,2); Dk = kr(iDim,3)
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         DiffCoef_VFDB(:,i,j,k,iDim,iBlock) = 0.5*( &
              DiffSemiCoef_VGB(:,i-Di,j-Dj,k-Dk,iBlock) &
              + DiffSemiCoef_VGB(:,i,j,k,iBlock))
      enddo; enddo; enddo

    end subroutine face_equal

    !==========================================================================

    subroutine face_left_coarse2fine(iDim,iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iDim, iMin, iMax, jMin, jMax, kMin, kMax

      !------------------------------------------------------------------------

      Di = kr(iDim,1); Dj = kr(iDim,2); Dk = kr(iDim,3)
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         DiffCoef_VFDB(:,i,j,k,iDim,iBlock) =  &
              (DiffSemiCoef_VGB(:,i,j,k,iBlock) &
              + 2.0*DiffSemiCoef_VGB(:,i-Di,j-Dj,k-Dk,iBlock))*4.0/9.0
      enddo; enddo; enddo

    end subroutine face_left_coarse2fine

    !==========================================================================

    subroutine face_right_coarse2fine(iDim,iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iDim, iMin, iMax, jMin, jMax, kMin, kMax

      !------------------------------------------------------------------------

      Di = kr(iDim,1); Dj = kr(iDim,2); Dk = kr(iDim,3)
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         DiffCoef_VFDB(:,i,j,k,iDim,iBlock) =  &
              (DiffSemiCoef_VGB(:,i-Di,j-Dj,k-Dk,iBlock) &
              + 2.0*DiffSemiCoef_VGB(:,i,j,k,iBlock))*4.0/9.0
      enddo; enddo; enddo

    end subroutine face_right_coarse2fine

    !==========================================================================

    subroutine face_left_fine2coarse(iDim,iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iDim, iMin, iMax, jMin, jMax, kMin, kMax

      integer :: iShift, jShift, kShift
      !------------------------------------------------------------------------

      Di = kr(iDim,1); Dj = kr(iDim,2); Dk = kr(iDim,3)
      iShift = 1-Di; jShift = 1-Dj; kShift = 1-Dk
      do k=kMin,kMax,2-Dk; do j=jMin,jMax,2-Dj; do i=iMin,iMax,2-Di
         do iDiff = 1, nDiff
            DiffCoef_VFDB(iDiff,i:i+iShift,j:j+jShift,k:k+kShift,iDim,iBlock) &
                 = (DiffSemiCoef_VGB(iDiff,i-Di,j-Dj,k-Dk,iBlock) &
                 + 0.5*sum(DiffSemiCoef_VGB( &
                 iDiff,i:i+iShift,j:j+jShift,k:k+kShift,iBlock)))*2.0/9.0
         end do
      enddo; enddo; enddo

    end subroutine face_left_fine2coarse

    !==========================================================================

    subroutine face_right_fine2coarse(iDim,iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iDim, iMin, iMax, jMin, jMax, kMin, kMax

      integer :: iShift, jShift, kShift, i1, j1, k1
      !------------------------------------------------------------------------

      Di = kr(iDim,1); Dj = kr(iDim,2); Dk = kr(iDim,3)
      iShift = 1-Di; jShift = 1-Dj; kShift = 1-Dk
      do k=kMin,kMax,2-Dk; do j=jMin,jMax,2-Dj; do i=iMin,iMax,2-Di
         i1=i-Di; j1=j-Dj; k1=k-Dk
         do iDiff = 1, nDiff
            DiffCoef_VFDB(iDiff,i:i+iShift,j:j+jShift,k:k+kShift,iDim,iBlock) &
                 = (DiffSemiCoef_VGB(iDiff,i,j,k,iBlock) &
                 + 0.5*sum(DiffSemiCoef_VGB( &
                 iDiff,i1:i1+iShift,j1:j1+jShift,k1:k1+kShift,iBlock)))*2.0/9.0
         end do
      enddo; enddo; enddo

    end subroutine face_right_fine2coarse

    !==========================================================================

    subroutine get_diffusion_coef

      !------------------------------------------------------------------------

      if(TypeSemiImplicit=='radiation' .or. TypeSemiImplicit=='radcond')then
         RosselandMeanOpacity = RosselandMeanOpacitySi/Si2No_V(UnitX_)

         ! Calculate the cell centered diffusion coefficients
         if(UseRadFluxLimiter)then

            if(IsNewBlockGrayDiffusion)then
               Erad_G = State_VGB(Eradiation_,:,:,:,iBlock)
               call set_block_scalar(Erad_G, iBlock)
            end if
        
            Grad2ByErad2 = &
                 (((Erad_G(i+1,j,k) - Erad_G(i-1,j,k))*InvDx2)**2 &
                 +((Erad_G(i,j+1,k) - Erad_G(i,j-1,k))*InvDy2)**2 &
                 +((Erad_G(i,j,k+1) - Erad_G(i,j,k-1))*InvDz2)**2 &
                 )/ Erad_G(i,j,k)**2

            select case(TypeRadFluxLimiter)
            case("sum")
               DiffRad = Clight/(3*RosselandMeanOpacity + sqrt(Grad2ByErad2))
            case("max")
               DiffRad = Clight/max(3*RosselandMeanOpacity,sqrt(Grad2ByErad2))
            case("larsen")
               DiffRad = Clight/sqrt(9*RosselandMeanOpacity**2 + Grad2ByErad2)
            end select
         else
            DiffRad = Clight/(3*RosselandMeanOpacity)
         end if

         ! Store it for message passing
         DiffSemiCoef_VGB(EradImpl_,i,j,k,iBlock) = DiffRad
      end if

      if(TypeSemiImplicit=='cond')then
         HeatConductionCoef = HeatConductionCoefSi &
              *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_) &
              *Si2No_V(UnitU_)*Si2No_V(UnitX_)

         DiffSemiCoef_VGB(TeImpl_,i,j,k,iBlock) = HeatConductionCoef
      elseif(TypeSemiImplicit=='radcond')then
         HeatConductionCoef = HeatConductionCoefSi &
              *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_) &
              *Si2No_V(UnitU_)*Si2No_V(UnitX_)
         Te = TeSi*Si2No_V(UnitTemperature_)

         DiffSemiCoef_VGB(TeImpl_,i,j,k,iBlock) = &
              HeatConductionCoef/(4.0*cRadiationNo*Te**3)
      end if

    end subroutine get_diffusion_coef

  end subroutine get_impl_gray_diff_state

  !============================================================================

  subroutine get_gray_diffusion_bc(iBlock, IsLinear)

    use ModGeometry, ONLY: TypeGeometry
    use ModImplicit, ONLY: StateSemi_VGB
    use ModMain,     ONLY: nI, nJ, nK, TypeBc_I
    use ModParallel, ONLY: NOBLK, NeiLev
    use ModUser,     ONLY: user_set_outerbcs

    integer, intent(in) :: iBlock
    logical, intent(in) :: IsLinear

    logical :: IsFound
    character(len=20), parameter :: TypeUserBc='usersemi'
    character(len=*), parameter :: NameSub='get_gray_diffusion_bc'
    !--------------------------------------------------------------------------

    if(NeiLev(1,iBlock) == NOBLK)then
       if(IsLinear .and. TypeBc_I(1) /= 'reflect')then
          StateSemi_VGB(:,0,:,:,iBlock) = 0.0
       elseif(TypeBc_I(1) == 'user')then
          IsFound = .false.
          call user_set_outerbcs(iBlock,1,TypeUserBc,IsFound)
          if(.not. IsFound) call stop_mpi( &
               NameSub//': unknown TypeBc_I='//TypeBc_I(1))
       else
          StateSemi_VGB(:,0,:,:,iBlock) = StateSemi_VGB(:,1,:,:,iBlock)
       end if
    end if
    if(NeiLev(2,iBlock) == NOBLK)then
       if(IsLinear .and. TypeBc_I(2) /= 'reflect')then
          StateSemi_VGB(:,nI+1,:,:,iBlock) = 0.0
       elseif(TypeBc_I(2) == 'user')then
          IsFound = .false.
          call user_set_outerbcs(iBlock,2,TypeUserBc,IsFound)
          if(.not. IsFound) call stop_mpi( &
               NameSub//': unknown TypeBc_I='//TypeBc_I(2))
       else
          StateSemi_VGB(:,nI+1,:,:,iBlock) = StateSemi_VGB(:,nI,:,:,iBlock)
       end if
    end if
    if(NeiLev(3,iBlock) == NOBLK)then
       if(TypeBc_I(3) == 'shear')then
          StateSemi_VGB(:,:,0,:,iBlock) = StateSemi_VGB(:,:,1,:,iBlock)
          call semi_bc_shear(3)
       elseif(IsLinear .and. TypeBc_I(3) /= 'reflect')then
          StateSemi_VGB(:,:,0,:,iBlock) = 0.0
       elseif(TypeBc_I(3) == 'user')then
          IsFound = .false.
          call user_set_outerbcs(iBlock,3,TypeUserBc,IsFound)
          if(.not. IsFound) call stop_mpi( &
               NameSub//': unknown TypeBc_I='//TypeBc_I(3))
       else
          StateSemi_VGB(:,:,0,:,iBlock) = StateSemi_VGB(:,:,1,:,iBlock)
       end if
    end if
    if(NeiLev(4,iBlock) == NOBLK) then
       if(TypeBc_I(4) == 'shear')then
          StateSemi_VGB(:,:,nJ+1,:,iBlock) = StateSemi_VGB(:,:,nJ,:,iBlock)
          call semi_bc_shear(4)
       elseif(IsLinear .and. TypeBc_I(4) /= 'reflect')then
          StateSemi_VGB(:,:,nJ+1,:,iBlock) = 0.0
       elseif(TypeBc_I(4) == 'user')then
          IsFound = .false.
          call user_set_outerbcs(iBlock,4,TypeUserBc,IsFound)
          if(.not. IsFound) call stop_mpi( &
               NameSub//': unknown TypeBc_I='//TypeBc_I(4))
       else
          StateSemi_VGB(:,:,nJ+1,:,iBlock) = StateSemi_VGB(:,:,nJ,:,iBlock)
       end if
    end if
    if(NeiLev(5,iBlock) == NOBLK) then
       if(IsLinear .and. TypeBc_I(5) /= 'reflect')then
          StateSemi_VGB(:,:,:,0,iBlock) = 0.0
       elseif(TypeBc_I(5) == 'user')then
          IsFound = .false.
          call user_set_outerbcs(iBlock,5,TypeUserBc,IsFound)
          if(.not. IsFound) call stop_mpi( &
               NameSub//': unknown TypeBc_I='//TypeBc_I(5))
       else
          StateSemi_VGB(:,:,:,0,iBlock) = StateSemi_VGB(:,:,:,1,iBlock)
       end if
    end if
    if(NeiLev(6,iBlock) == NOBLK)then
       if(IsLinear .and. TypeBc_I(6) /= 'reflect')then
          StateSemi_VGB(:,:,:,nK+1,iBlock) = 0.0
       elseif(TypeBc_I(6) == 'user')then
          IsFound = .false.
          call user_set_outerbcs(iBlock,6,TypeUserBc,IsFound)
          if(.not. IsFound) call stop_mpi( &
               NameSub//': unknown TypeBc_I='//TypeBc_I(6))
       else
          StateSemi_VGB(:,:,:,nK+1,iBlock) = StateSemi_VGB(:,:,:,nK,iBlock)
       end if
    end if

    if(NeiLev(1,iBlock)==1) call correct_left_ghostcell(1,0,0,1,nJ,1,nK)
    if(NeiLev(2,iBlock)==1) call correct_right_ghostcell(1,nI+1,nI+1,1,nJ,1,nK)
    if(NeiLev(3,iBlock)==1) call correct_left_ghostcell(2,1,nI,0,0,1,nK)
    if(NeiLev(4,iBlock)==1) call correct_right_ghostcell(2,1,nI,nJ+1,nJ+1,1,nK)
    if(TypeGeometry /= 'rz')then
       if(NeiLev(5,iBlock)==1) call correct_left_ghostcell(3,1,nI,1,nJ,0,0)
       if(NeiLev(6,iBlock)==1) &
            call correct_right_ghostcell(3,1,nI,1,nJ,nK+1,nK+1)
    end if

  contains

    subroutine semi_bc_shear(iSide)

      use ModNumConst, ONLY: cTiny
      use ModPhysics, ONLY: ShockSlope
      use ModSize, ONLY: south_, north_

      integer, intent(in) :: iSide

      integer :: Dn, iDiff, iVar
      !------------------------------------------------------------------------

      ! If the shock is not tilted, there is nothing to do
      if(abs(ShockSlope)<cTiny) RETURN

      do iDiff = 1, nDiff
         iVar = iDiff_I(iDiff)

         ! Shear according to ShockSlope
         if(ShockSlope < -cTiny)then
            call stop_mpi('ShockSlope must be positive!')
         elseif(ShockSlope > 1.0)then
            call stop_mpi('ShockSlope > 1 not allowed!')
         else
            ! ShockSlope <= 1
            Dn = nint(1.0/ShockSlope)
            if(abs(Dn-1.0/ShockSlope)>cTiny)call stop_mpi( &
                 'ShockSlope <= 1 should be the inverse of a round number!')
            select case(iSide)
               ! Shift parallel to X by 1, but copy from distance Dn in Y
            case(south_)
               StateSemi_VGB(iVar,1:nI,0,:,iBlock) = &
                    StateSemi_VGB(iVar,0:nI-1,Dn,:,iBlock)
            case(north_)
               StateSemi_VGB(iVar,1:nI,nJ+1,:,iBlock) = &
                    StateSemi_VGB(iVar,2:nI+1,nJ+1-Dn,:,iBlock)
            end select
         end if
      end do

    end subroutine semi_bc_shear

    !==========================================================================

    subroutine correct_left_ghostcell(iDim,iMin,iMax,jMin,jMax,kMin,kMax)

      use ModImplicit, ONLY: kr

      integer, intent(in) :: iDim, iMin, iMax, jMin, jMax, kMin, kMax

      integer :: i, j, k, iShift, jShift, kShift, Di, Dj, Dk, i1, j1, k1
      integer :: iDiff, iVar
      !------------------------------------------------------------------------

      Di = kr(iDim,1); Dj = kr(iDim,2); Dk = kr(iDim,3)
      iShift = 1-Di; jShift = 1-Dj; kShift = 1-Dk
      do k=kMin,kMax,2-Dk; do j=jMin,jMax,2-Dj; do i=iMin,iMax,2-Di
         i1=i+Di; j1=j+Dj; k1=k+Dk
         do iDiff = 1, nDiff
            iVar = iDiff_I(iDiff)
            StateSemi_VGB(iVar,i:i+iShift,j:j+jShift,k:k+kShift,iBlock) = &
                 StateSemi_VGB(iVar,i:i+iShift,j:j+jShift,k:k+kShift,iBlock) &
                 + StateSemi_VGB( &
                 iVar,i1:i1+iShift,j1:j1+jShift,k1:k1+kShift,iBlock) &
                 -0.25*sum(StateSemi_VGB( &
                 iVar,i1:i1+iShift,j1:j1+jShift,k1:k1+kShift,iBlock))
         end do
      enddo; enddo; enddo

    end subroutine correct_left_ghostcell

    !==========================================================================

    subroutine correct_right_ghostcell(iDim,iMin,iMax,jMin,jMax,kMin,kMax)

      use ModImplicit, ONLY: kr

      integer, intent(in) :: iDim, iMin, iMax, jMin, jMax, kMin, kMax

      integer :: i, j, k, iShift, jShift, kShift, Di, Dj, Dk, i1, j1, k1
      integer :: iDiff, iVar
      !------------------------------------------------------------------------

      Di = kr(iDim,1); Dj = kr(iDim,2); Dk = kr(iDim,3)
      iShift = 1-Di; jShift = 1-Dj; kShift = 1-Dk
      do k=kMin,kMax,2-Dk; do j=jMin,jMax,2-Dj; do i=iMin,iMax,2-Di
         i1=i-Di; j1=j-Dj; k1=k-Dk
         do iDiff = 1, nDiff
            iVar = iDiff_I(iDiff)
            StateSemi_VGB(iVar,i:i+iShift,j:j+jShift,k:k+kShift,iBlock) = &
                 StateSemi_VGB(iVar,i:i+iShift,j:j+jShift,k:k+kShift,iBlock) &
                 + StateSemi_VGB( &
                 iVar,i1:i1+iShift,j1:j1+jShift,k1:k1+kShift,iBlock) &
                 -0.25*sum(StateSemi_VGB( &
                 iVar,i1:i1+iShift,j1:j1+jShift,k1:k1+kShift,iBlock))
         end do
      enddo; enddo; enddo

    end subroutine correct_right_ghostcell

  end subroutine get_gray_diffusion_bc

  !============================================================================

  subroutine get_gray_diffusion_rhs(iBlock, StateImpl_VG, Rhs_VC, IsLinear)

    use ModGeometry, ONLY: TypeGeometry, vInv_CB
    use ModImplicit, ONLY: nw
    use ModMain,     ONLY: nI, nJ, nK

    integer, intent(in) :: iBlock
    real, intent(inout) :: StateImpl_VG(nw,-1:nI+2,-1:nJ+2,-1:nK+2)
    real, intent(out)   :: Rhs_VC(nw,nI,nJ,nK)
    logical, intent(in) :: IsLinear

    real :: Te, EnergyExchange
    integer :: i, j, k, iDiff, iRelax, iVar
    character(len=*), parameter :: NameSub='get_gray_diffusion_rhs'
    !--------------------------------------------------------------------------

    !!! Rhs_VC = 0.0

    if(TypeGeometry == 'rz')then
       ! No flux from Z direction
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          do iDiff = 1, nDiff
             iVar = iDiff_I(iDiff)
             Rhs_VC(iVar,i,j,k) = &
                  vInv_CB(i,j,k,iBlock) * ( &
                  + DiffCoef_VFDB(iDiff,i+1,j,k,1,iBlock)* &
                  (   StateImpl_VG(iVar,i+1,j,k)   &
                  -   StateImpl_VG(iVar,i  ,j,k))  &
                  - DiffCoef_VFDB(iDiff,i  ,j,k,1,iBlock)* &
                  (   StateImpl_VG(iVar,i  ,j,k)   &
                  -   StateImpl_VG(iVar,i-1,j,k))  &
                  + DiffCoef_VFDB(iDiff,i,j+1,k,2,iBlock)* &
                  (   StateImpl_VG(iVar,i,j+1,k)   &
                  -   StateImpl_VG(iVar,i,j  ,k))  &
                  - DiffCoef_VFDB(iDiff,i,j  ,k,2,iBlock)* &
                  (   StateImpl_VG(iVar,i,j  ,k)   &
                  -   StateImpl_VG(iVar,i,j-1,k)) )
          end do
       end do; end do; end do
    else
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          do iDiff = 1, nDiff
             iVar = iDiff_I(iDiff)
             Rhs_VC(iVar,i,j,k) = &
                  vInv_CB(i,j,k,iBlock) * ( &
                  DiffCoef_VFDB(iDiff,i+1,j,k,1,iBlock)*   &
                  (   StateImpl_VG(iVar,i+1,j,k)   &
                  -   StateImpl_VG(iVar,i  ,j,k))  &
                  - DiffCoef_VFDB(iDiff,i  ,j,k,1,iBlock)* &
                  (   StateImpl_VG(iVar,i  ,j,k)   &
                  -   StateImpl_VG(iVar,i-1,j,k))  &
                  + DiffCoef_VFDB(iDiff,i,j+1,k,2,iBlock)* &
                  (   StateImpl_VG(iVar,i,j+1,k)   &
                  -   StateImpl_VG(iVar,i,j  ,k))  &
                  - DiffCoef_VFDB(iDiff,i,j  ,k,2,iBlock)* &
                  (   StateImpl_VG(iVar,i,j  ,k)   &
                  -   StateImpl_VG(iVar,i,j-1,k))  &
                  + DiffCoef_VFDB(iDiff,i,j,k+1,3,iBlock)* &
                  (   StateImpl_VG(iVar,i,j,k+1)   &
                  -   StateImpl_VG(iVar,i,j,k  ))  &
                  - DiffCoef_VFDB(iDiff,i,j,k  ,3,iBlock)* &
                  (   StateImpl_VG(iVar,i,j,k  )   &
                  -   StateImpl_VG(iVar,i,j,k-1)) )
          end do
       end do; end do; end do

    end if

    ! Source terms due to energy exchange
    if(nRelax>0)then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          do iRelax = 1, nRelax
             iVar = iRelax_I(iRelax)

             EnergyExchange = RelaxSemiCoef_VCB(iRelax,i,j,k,iBlock) &
                  *(StateImpl_VG(TeImpl_,i,j,k) - StateImpl_VG(iVar,i,j,k))

             ! dEvar/dt = + EnergyExchange
             Rhs_VC(iVar,i,j,k)    = Rhs_VC(iVar,i,j,k)    + EnergyExchange

             ! dEe/dt   = - EnergyExchange
             Rhs_VC(TeImpl_,i,j,k) = Rhs_VC(TeImpl_,i,j,k) - EnergyExchange
          end do
       end do; end do; end do
    end if
       
    ! Point implicit source terms due to energy exchange
    if(nPoint>0)then
       if(IsLinear)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             Rhs_VC(iPoint,i,j,k) = Rhs_VC(iPoint,i,j,k) &
                  - PointSemiCoef_VCB(Relax_,i,j,k,iBlock) &
                  *StateImpl_VG(iPoint,i,j,k)
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             EnergyExchange = PointSemiCoef_VCB(Relax_,i,j,k,iBlock) &
                  * (PointSemiCoef_VCB(Planck_,i,j,k,iBlock) &
                  -  StateImpl_VG(iPoint,i,j,k))

             Rhs_VC(iPoint,i,j,k) = Rhs_VC(iPoint,i,j,k) + EnergyExchange
          end do; end do; end do
       end if
    end if

  end subroutine get_gray_diffusion_rhs

  !============================================================================

  subroutine get_gray_diff_jacobian(iBlock, nVar, Jacobian_VVCI)

    use ModGeometry, ONLY: vInv_CB
    use ModImplicit, ONLY: kr, TypeSemiImplicit
    use ModMain,     ONLY: nI, nJ, nK, nDim

    integer, parameter:: nStencil = 2*nDim + 1

    integer, intent(in) :: iBlock, nVar
    real, intent(out) :: Jacobian_VVCI(nVar,nVar,nI,nJ,nK,nStencil)

    integer :: iVar, i, j, k, iDim, Di, Dj, Dk, iDiff, iRelax
    real :: DiffLeft, DiffRight
    real :: RelaxCoef
    !--------------------------------------------------------------------------

    ! All elements have to be set
    Jacobian_VVCI(:,:,:,:,:,:) = 0.0

    if(nPoint>0)then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          ! dSvar/dVar (diagonal)
          Jacobian_VVCI(iPoint,iPoint,i,j,k,1) = &
               -PointSemiCoef_VCB(Relax_,i,j,k,iBlock)
       end do; end do; end do
    end if

    if(nRelax>0)then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          do iRelax = 1, nRelax
             iVar = iRelax_I(iRelax)

             RelaxCoef = RelaxSemiCoef_VCB(iRelax,i,j,k,iBlock)

             ! dSvar/dVar (diagonal)
             Jacobian_VVCI(iVar,iVar,i,j,k,1) = -RelaxCoef

             ! dSe/dVar (off diagonal)
             Jacobian_VVCI(TeImpl_,iVar,i,j,k,1) = +RelaxCoef

             ! dSe/daTe^4 (diagonal)
             Jacobian_VVCI(TeImpl_,TeImpl_,i,j,k,1) = -RelaxCoef

             ! dSvar/daTe^4 (off diagonal)
             Jacobian_VVCI(iVar,TeImpl_,i,j,k,1) = +RelaxCoef
          end do
       end do; end do; end do
    end if

    do iDim = 1, nDim
       Di = kr(iDim,1); Dj = kr(iDim,2); Dk = kr(iDim,3)
       do k=1,nK; do j=1,nJ; do i=1,nI
          do iDiff = 1, nDiff
             iVar = iDiff_I(iDiff)
             DiffLeft = vInv_CB(i,j,k,iBlock) &
                  *DiffCoef_VFDB(iDiff,i,j,k,iDim,iBlock)
             DiffRight = vInv_CB(i,j,k,iBlock) &
                  *DiffCoef_VFDB(iDiff,i+Di,j+Dj,k+Dk,iDim,iBlock)
             Jacobian_VVCI(iVar,iVar,i,j,k,1) = &
                  Jacobian_VVCI(iVar,iVar,i,j,k,1) - (DiffLeft + DiffRight)

             if(iDim==1.and.i==1 .or. iDim==2.and.j==1 .or. iDim==3.and.k==1)&
                  DiffLeft = 0.0
             if(iDim==1.and.i==nI .or. iDim==2.and.j==nJ &
                  .or. iDim==3.and.k==nK) DiffRight = 0.0

             Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim)   = DiffLeft
             Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim+1) = DiffRight
          end do
       end do; end do; end do
    end do

  end subroutine get_gray_diff_jacobian

  !============================================================================

  subroutine update_impl_gray_diff(iBlock, iImplBlock, StateImpl_VG)

    use ModAdvance,  ONLY: State_VGB, Rho_, p_, Eradiation_
    use ModEnergy,   ONLY: calc_energy_cell
    use ModImplicit, ONLY: nw, TypeSemiImplicit, DconsDsemi_VCB, ImplOld_VCB, &
         ImplCoeff
    use ModMain,     ONLY: nI, nJ, nK, dt
    use ModPhysics,  ONLY: inv_gm1, No2Si_V, Si2No_V, UnitEnergyDens_, UnitP_,&
         UnitRho_, UnitTemperature_
    use ModUser,     ONLY: user_material_properties

    integer, intent(in) :: iBlock, iImplBlock
    real, intent(in) :: StateImpl_VG(nw,nI,nJ,nK)

    integer :: i, j, k
    real :: Einternal, EinternalSi, PressureSi, AbsorptionEmission

    character(len=*), parameter :: NameSub = 'update_impl_gray_diff'
    !--------------------------------------------------------------------------

    if(TypeSemiImplicit=='radiation' .or. TypeSemiImplicit=='radcond')then
       do k = 1,nK; do j = 1,nJ; do i = 1,nI
          State_VGB(Eradiation_,i,j,k,iBlock) = &
               max(EradMin, StateImpl_VG(EradImpl_,i,j,k))

          if(State_VGB(Eradiation_,i,j,k,iBlock) < 0.0)then
             write(*,*)NameSub,': ERROR EradMin, EradOrig=', &
                  EradMin, StateImpl_VG(EradImpl_,i,j,k)

             write(*,*)NameSub,': ERROR negative Erad =', &
                  State_VGB(Eradiation_,i,j,k,iBlock)
             write(*,*)NameSub,': ERROR at i,j,k,iBlock=', i, j, k, iBlock
             call stop_mpi(NameSub//' negative Erad')
          end if
       end do; end do; end do
    end if
    
    do k = 1,nK; do j = 1,nJ; do i = 1,nI
       if(TypeSemiImplicit=='radiation')then
          AbsorptionEmission = ImplCoeff &
               *PointSemiCoef_VCB(Relax_,i,j,k,iBlock) &
               * (PointSemiCoef_VCB(Planck_,i,j,k,iBlock) &
               -  State_VGB(Eradiation_,i,j,k,iBlock)) &
               + (1.0-ImplCoeff)*PointSemiCoef_VCB(Relax_,i,j,k,iBlock) &
               *(PointSemiCoef_VCB(Planck_,i,j,k,iBlock) &
               - ImplOld_VCB(EradImpl_,i,j,k,iBlock))
               
          Einternal = inv_gm1*State_VGB(p_,i,j,k,iBlock) &
               + State_VGB(EintExtra_,i,j,k,iBlock) &
               - dt*AbsorptionEmission
       else
          Einternal = inv_gm1*State_VGB(p_,i,j,k,iBlock) &
               + State_VGB(EintExtra_,i,j,k,iBlock) &
               + DconsDsemi_VCB(TeImpl_,i,j,k,iImplBlock) &
               *( StateImpl_VG(TeImpl_,i,j,k) &
               -  ImplOld_VCB(TeImpl_,i,j,k,iBlock) )
       end if
       EinternalSi = Einternal*No2Si_V(UnitEnergyDens_)

       if(Einternal < 0.0)then
          write(*,*)NameSub,': ERROR Rho, p, TOrigSi=', &
               State_VGB(Rho_,i,j,k,iBlock)*No2Si_V(UnitRho_), &
               State_VGB(p_,i,j,k,iBlock)*No2Si_V(UnitP_), &
               State_VGB(p_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock) &
               *No2Si_V(UnitTemperature_)

          write(*,*)NameSub,': ERROR negative Eint=', Einternal
          write(*,*)NameSub,': ERROR at i,j,k,iBlock=', i, j, k, iBlock
          call stop_mpi(NameSub//' negative Eint')
       end if

       call user_material_properties(State_VGB(:,i,j,k,iBlock), &
            EinternalSiIn = EinternalSi, PressureSiOut = PressureSi)

       State_VGB(p_,i,j,k,iBlock) = PressureSi*Si2No_V(UnitP_)

       State_VGB(EintExtra_,i,j,k,iBlock) = &
            Einternal - inv_gm1*State_VGB(p_,i,j,k,iBlock)

    end do; end do; end do

    call calc_energy_cell(iBlock)

  end subroutine update_impl_gray_diff

end module ModGrayDiffusion
