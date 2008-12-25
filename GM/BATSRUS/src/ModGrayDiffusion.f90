!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
!============================================================================
module ModGrayDiffusion

  use ModSize,       ONLY: nI, nJ, nK, nDim
  use ModVarIndexes, ONLY: p_

  implicit none
  save

  ! This module is needed for the order(u/c) gray diffusion approximation for
  ! radiation-hydrodynamics.

  private !except

  ! Public methods
  public :: init_gray_diffusion
  public :: set_frozen_coefficients
  public :: get_radiation_energy_flux
  public :: calc_source_gray_diffusion
  public :: get_impl_gray_diff_state
  public :: get_gray_diffusion_rhs
  public :: update_impl_gray_diff

  ! Logical for adding Gray Diffusion
  logical, public :: IsNewBlockGrayDiffusion = .true.
  logical, public :: IsNewTimestepGrayDiffusion = .true.
  logical, public :: DoUpdateFrozenCoefficients = .false.

  ! Parameters for radiation flux limiter
  logical,           public :: UseRadFluxLimiter  = .false.
  character(len=20), public :: TypeRadFluxLimiter = 'larsen'

  ! Coefficients for two-temperature electron-radiation model
  real, allocatable, public :: DiffusionRad_FDB(:,:,:,:,:)
  real, allocatable         :: RelaxationCoef_CB(:,:,:,:)
  real, allocatable         :: SpecificHeat_VCB(:,:,:,:,:)

  ! Internal energy needed for energy update
  real, allocatable :: Eint_VCB(:,:,:,:,:)

  ! Local variables
  real :: Erad_G(0:nI+1,0:nJ+1,0:nK+1)

  real, parameter :: GammaRel = 4.0/3.0

  ! local index for extra internal energy to keep the compiler happy
  integer, parameter :: EintExtra_ = p_ - 1

  ! named indexes for semi-implicit scheme
  integer, parameter :: TeImpl_ = 1, TradImpl_ = 2

contains

  !==========================================================================
  subroutine init_gray_diffusion

    use ModImplicit, ONLY: nw, UseSemiImplicit
    use ModSize,     ONLY: nBlk

    !------------------------------------------------------------------------

    if(.not.allocated(DiffusionRad_FDB)) allocate( &
         DiffusionRad_FDB(1:nI+1,1:nJ+1,1:nK+1,nDim,nBlk), &
         RelaxationCoef_CB(1:nI,1:nJ,1:nK,nBlk) )

    if(UseSemiImplicit)then
       if(.not.allocated(SpecificHeat_VCB)) allocate( &
            SpecificHeat_VCB(1:nw,1:nI,1:nJ,1:nK,nBlk), &
            Eint_VCB(1:nw,1:nI,1:nJ,1:nK,nBlk) )
    end if

  end subroutine init_gray_diffusion

  !==========================================================================
  subroutine set_frozen_coefficients(iBlock)

    use ModAdvance,    ONLY: State_VGB, Eradiation_
    use ModConst,      ONLY: cLightSpeed
    use ModImplicit,   ONLY: UseSemiImplicit
    use ModPhysics,    ONLY: Si2No_V, UnitT_, UnitEnergyDens_, &
         UnitTemperature_, cRadiationNo
    use ModUser,       ONLY: user_material_properties

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: PlanckOpacitySi
    real :: CvSi, TeSi, Te, Trad
    !------------------------------------------------------------------------

    do k=1,nK; do j=1,nJ; do i=1,nI
       call user_material_properties(State_VGB(:,i,j,k,iBlock), &
            AbsorptionOpacitySiOut = PlanckOpacitySi, &
            TeSiOut = TeSi, CvSiOut = CvSi)

       RelaxationCoef_CB(i,j,k,iBlock) = &
            PlanckOpacitySi*cLightSpeed/Si2No_V(UnitT_)

       if(.not.UseSemiImplicit) CYCLE

       Te = TeSi*Si2No_V(UnitTemperature_)
       Trad = sqrt(sqrt(State_VGB(Eradiation_,i,j,k,iBlock)/cRadiationNo))

       RelaxationCoef_CB(i,j,k,iBlock) = RelaxationCoef_CB(i,j,k,iBlock) &
            *cRadiationNo*(Te+Trad)*(Te**2+Trad**2)

       SpecificHeat_VCB(TeImpl_,i,j,k,iBlock) = CvSi &
            *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)

       SpecificHeat_VCB(TradImpl_,i,j,k,iBlock) = &
            4.0*cRadiationNo*Trad**3

    end do; end do; end do

  end subroutine set_frozen_coefficients
  !==========================================================================

  subroutine get_radiation_energy_flux( &
       iDir, i, j, k, iBlock, State_V, DiffRad, EradFlux_D)

    !\
    ! Calculate the diffusion part of the radiation energy flux.
    !/
    use ModAdvance,    ONLY: State_VGB, Eradiation_
    use ModConst,      ONLY: cLightSpeed
    use ModImplicit,   ONLY: UseSemiImplicit
    use ModPhysics,    ONLY: Si2No_V, UnitX_, UnitU_, cRadiationNo
    use ModUser,       ONLY: user_material_properties
    use ModVarIndexes, ONLY: nVar

    integer, intent(in) :: iDir, i, j, k, iBlock
    real,    intent(in) :: State_V(nVar)
    real,    intent(out):: DiffRad
    real,    intent(out):: EradFlux_D(3)

    real :: RosselandMeanOpacitySi
    real :: FaceGrad_D(3), Erad, Grad2ByErad2, Trad
    !------------------------------------------------------------------------

    call calc_face_gradient(iDir, i, j, k, iBlock, FaceGrad_D)


    if(DoUpdateFrozenCoefficients)then

       call user_material_properties(State_V, &
            RosselandMeanOpacitySiOut = RosselandMeanOpacitySi)

       DiffRad = cLightSpeed/(3.0*RosselandMeanOpacitySi) &
            *Si2No_V(UnitU_)*Si2No_V(UnitX_)

       if(UseRadFluxLimiter)then

          Grad2ByErad2 = sum(FaceGrad_D**2)/State_V(Eradiation_)**2

          select case(TypeRadFluxLimiter)
          case("sum")
             DiffRad = 1.0/(1.0/DiffRad+sqrt(Grad2ByErad2))
          case("max")
             DiffRad = 1.0/max(1.0/DiffRad,sqrt(Grad2ByErad2))
          case("larsen")
             DiffRad = 1.0/sqrt(1.0/DiffRad**2+Grad2ByErad2)
          end select

       end if

       if(UseSemiImplicit)then
          Trad = sqrt(sqrt(State_V(Eradiation_)/cRadiationNo))
          DiffRad = DiffRad*4.0*cRadiationNo*Trad**3
       end if
       DiffusionRad_FDB(i,j,k,iDir,iBlock) = DiffRad
    else
       DiffRad = DiffusionRad_FDB(i,j,k,iDir,iBlock)
    end if

    EradFlux_D = 0.0
    EradFlux_D(iDir) = -DiffRad*FaceGrad_D(iDir)

  end subroutine get_radiation_energy_flux

  !==========================================================================

  subroutine set_block_scalar(Scalar_G, iBlock)

    ! correct the ghostcells of the given scalar field on iBlock

    use ModParallel, ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot, BlkNeighborLev, NOBLK

    integer, intent(in) :: iBlock
    real, dimension(0:nI+1,0:nJ+1,0:nK+1), intent(inout) :: Scalar_G

    real, parameter :: c0 = 0.5, p0 = 1./6., F1 = 1./3.

    integer :: i1, j1, k1, i2, j2, k2, iC, jC, kC, iSide, jSide, kSide
    integer :: iL, iR, jL, jR, kL, kR
    integer :: ip, jp, kp

    real    :: Scalar1_G(0:nI+1,0:nJ+1,0:nK+1) ! temporary array
    logical :: IsEqualLevel_G(0:nI+1,0:nJ+1,0:nK+1)
    !------------------------------------------------------------------------
    IsNewBlockGrayDiffusion = .false.

    Scalar1_G = Scalar_G

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
       Erad_G = State_VGB(Eradiation_,0:nI+1,0:nJ+1,0:nK+1,iBlock)
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
    use ModConst,      ONLY: cRadiation
    use ModGeometry,   ONLY: vInv_CB
    use ModPhysics,    ONLY: Si2No_V, UnitEnergyDens_
    use ModUser,       ONLY: user_material_properties
    use ModVarIndexes, ONLY: Energy_, NameVar_V
    use ModImplicit,   ONLY: UseSemiImplicit

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: TeSi
    real :: RadCompression, AbsorptionEmission
    character(len=*), parameter:: NameSub = "calc_source_gray_diffusion"
    !------------------------------------------------------------------------

    ! Make sure that Eradiation_ is correct
    if(NameVar_V(Eradiation_) /= "Erad") call stop_mpi(NameSub// &
         ": incorrect index for Erad variable in ModEquation")

    do k=1,nK; do j=1,nJ; do i=1,nI

       ! Adiabatic compression of radiation energy by fluid velocity (fluid 1)
       ! (GammaRel-1)*Erad*Div(U)
       RadCompression = (GammaRel-1.0)*State_VGB(Eradiation_,i,j,k,iBlock) &
            *vInv_CB(i,j,k,iBlock)&
            *(uDotArea_XI(i+1,j,k,1) - uDotArea_XI(i,j,k,1) &
             +uDotArea_YI(i,j+1,k,1) - uDotArea_YI(i,j,k,1) &
             +uDotArea_ZI(i,j,k+1,1) - uDotArea_ZI(i,j,k,1))

       ! dErad/dt = - adiabatic compression + AbsorptionEmission
       Source_VC(Eradiation_,i,j,k) = Source_VC(Eradiation_,i,j,k) &
            - RadCompression

       ! dE/dt    = + adiabatic compression - AbsorptionEmission
       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + RadCompression

       ! The absorption-emission term will be added by the implicit solver
       if(UseSemiImplicit) CYCLE

       call user_material_properties(State_VGB(:,i,j,k,iBlock), &
            TeSiOut = TeSi)

       ! Source term due to absorption and emission
       ! Sigma_a*(cRadiation*Te**4-Erad)
       AbsorptionEmission =  RelaxationCoef_CB(i,j,k,iBlock)&
            *(cRadiation*TeSi**4*Si2No_V(UnitEnergyDens_) &
            - State_VGB(Eradiation_,i,j,k,iBlock))

       ! dErad/dt = - adiabatic compression + AbsorptionEmission
       Source_VC(Eradiation_,i,j,k) = Source_VC(Eradiation_,i,j,k) &
            + AbsorptionEmission

       ! dE/dt    = + adiabatic compression - AbsorptionEmission
       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
            - AbsorptionEmission

    end do; end do; end do

  end subroutine calc_source_gray_diffusion

  !==========================================================================

  subroutine get_impl_gray_diff_state(ImplState_GVB)

    use ModAdvance,  ONLY: Eradiation_, State_VGB
    use ModImplicit, ONLY: nw, nImplBlk, impl2iBlk
    use ModPhysics,  ONLY: inv_gm1, cRadiationNo, Si2No_V, UnitTemperature_
    use ModSize,     ONLY: MaxImplBlk
    use ModUser,     ONLY: user_material_properties

    real, intent(out) :: ImplState_GVB(0:nI+1,0:nJ+1,0:nK+1,nw,MaxImplBlk)

    integer :: iImplBlock, iBlock, i, j, k
    real :: TeSi, Einternal
    !------------------------------------------------------------------------

    do iImplBlock = 1, nImplBLK
       iBlock = impl2iBLK(iImplBlock)
       do k = 0, nK+1; do j = 0, nJ+1; do i = 0, nI+1
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               TeSiOut = TeSi)

          ImplState_GVB(i,j,k,TeImpl_,iImplBlock) = &
               TeSi*Si2No_V(UnitTemperature_)

          ImplState_GVB(i,j,k,TradImpl_,iImplBlock) = &
               sqrt(sqrt(State_VGB(Eradiation_,i,j,k,iBlock) &
               /cRadiationNo))
       end do; end do; end do

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          Einternal = inv_gm1*State_VGB(p_,i,j,k,iBlock) &
               + State_VGB(EintExtra_,i,j,k,iBlock)

          Eint_VCB(TeImpl_,i,j,k,iBlock) = &
               Einternal &
               - SpecificHeat_VCB(TeImpl_,i,j,k,iBlock) &
               *ImplState_GVB(i,j,k,TeImpl_,iImplBlock)

          Eint_VCB(TradImpl_,i,j,k,iBlock) = &
               State_VGB(Eradiation_,i,j,k,iBlock) &
               - SpecificHeat_VCB(TradImpl_,i,j,k,iBlock) &
               *ImplState_GVB(i,j,k,TradImpl_,iImplBlock)

       end do; end do; end do

    end do

  end subroutine get_impl_gray_diff_state

  !==========================================================================

  subroutine get_gray_diffusion_rhs(iBlock, StateImpl_VG, Rhs_CV)

    use ModAdvance,  ONLY: Flux_VX, Flux_VY, Flux_VZ
    use ModGeometry, ONLY: dx_BLK, dy_BLK, dz_BLK
    use ModImplicit, ONLY: nw
    use ModMain,     ONLY: x_, y_, z_
    use ModParallel, ONLY: NOBLK, NeiLev

    integer, intent(in) :: iBlock
    real, intent(inout) :: StateImpl_VG(nw,-1:nI+2,-1:nJ+2,-1:nK+2)
    real, intent(out)   :: Rhs_CV(nI,nJ,nK,nw)

    real :: InvDx2, InvDy2, InvDz2
    real :: AbsorptionEmission
    integer :: i, j, k, iVar
    !------------------------------------------------------------------------

!!! set floating outer boundary
    if(NeiLev(1,iBlock) == NOBLK) StateImpl_VG(:,0   ,:,:) &
         =                        StateImpl_VG(:,1   ,:,:)
    if(NeiLev(2,iBlock) == NOBLK) StateImpl_VG(:,nI+1,:,:) &
         =                        StateImpl_VG(:,nI  ,:,:)
    if(NeiLev(3,iBlock) == NOBLK) StateImpl_VG(:,:,0   ,:) &
         =                        StateImpl_VG(:,:,1   ,:)
    if(NeiLev(4,iBlock) == NOBLK) StateImpl_VG(:,:,nJ+1,:) &
         =                        StateImpl_VG(:,:,nJ  ,:)
    if(NeiLev(5,iBlock) == NOBLK) StateImpl_VG(:,:,:,0   ) &
         =                        StateImpl_VG(:,:,:,1   )
    if(NeiLev(6,iBlock) == NOBLK) StateImpl_VG(:,:,:,nK+1) &
         =                        StateImpl_VG(:,:,:,nK  )

    ! Calculate radiative diffusion fluxes
    InvDx2 = 1./dx_BLK(iBlock)**2
    InvDy2 = 1./dy_BLK(iBlock)**2
    InvDz2 = 1./dz_BLK(iBlock)**2
    do k = 1, nK; do j = 1, nJ; do i = 1, nI+1
       Flux_VX(TradImpl_,i,j,k) = InvDx2*DiffusionRad_FDB(i,j,k,x_,iBlock) &
            *(StateImpl_VG(TradImpl_,i,j,k) - StateImpl_VG(TradImpl_,i-1,j,k))
    end do; end do; end do
    do k = 1, nK; do j = 1, nJ+1; do i = 1, nI
       Flux_VY(TradImpl_,i,j,k) = InvDy2*DiffusionRad_FDB(i,j,k,y_,iBlock) &
            *(StateImpl_VG(TradImpl_,i,j,k) - StateImpl_VG(TradImpl_,i,j-1,k))
    end do; end do; end do
    do k = 1, nK+1; do j = 1, nJ; do i = 1, nI
       Flux_VZ(TradImpl_,i,j,k) = InvDz2*DiffusionRad_FDB(i,j,k,z_,iBlock) &
            *(StateImpl_VG(TradImpl_,i,j,k) - StateImpl_VG(TradImpl_,i,j,k-1))
    end do; end do; end do

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       Rhs_CV(i,j,k,TradImpl_) = &
            + Flux_VX(TradImpl_,i+1,j,k) - Flux_VX(TradImpl_,i,j,k) &
            + Flux_VY(TradImpl_,i,j+1,k) - Flux_VY(TradImpl_,i,j,k) &
            + Flux_VZ(TradImpl_,i,j,k+1) - Flux_VZ(TradImpl_,i,j,k)
    end do; end do; end do

    ! Source term due to absorption and emission
    ! Sigma_a*(cRadiation*Te**4-Erad)
    do k = 1, nK; do j = 1, nJ; do i = 1, nI

       AbsorptionEmission = RelaxationCoef_CB(i,j,k,iBlock)&
            *(StateImpl_VG(TeImpl_,i,j,k) - StateImpl_VG(TradImpl_,i,j,k))

       ! dErad/dt = + AbsorptionEmission
       Rhs_CV(i,j,k,TradImpl_) = Rhs_CV(i,j,k,TradImpl_) + AbsorptionEmission

       ! dE/dt    = - AbsorptionEmission
       Rhs_CV(i,j,k,TeImpl_)   =                         - AbsorptionEmission

    end do; end do; end do

    do iVar = 1, nw; do k = 1, nK; do j = 1, nJ; do i = 1, nI
       Rhs_CV(i,j,k,iVar) = Rhs_CV(i,j,k,iVar) &
            /SpecificHeat_VCB(iVar,i,j,k,iBlock)
    end do; end do; end do; end do

  end subroutine get_gray_diffusion_rhs

  !==========================================================================
  subroutine update_impl_gray_diff(iBlock, StateImpl_GV)

    use ModAdvance,  ONLY: State_VGB, p_, Eradiation_
    use ModEnergy,   ONLY: calc_energy_cell
    use ModImplicit, ONLY: nw
    use ModPhysics,  ONLY: inv_gm1, No2Si_V, UnitEnergyDens_
    use ModUser,     ONLY: user_material_properties

    integer, intent(in) :: iBlock
    real, intent(in) :: StateImpl_GV(nI,nJ,nK,nw)

    integer :: i, j, k
    real :: TeSi, EinternalSi, Einternal, Gamma
    !------------------------------------------------------------------------

    do k = 1,nK; do j = 1,nJ; do i = 1,nI

       Eint_VCB(:,i,j,k,iBlock) = Eint_VCB(:,i,j,k,iBlock) &
            + SpecificHeat_VCB(:,i,j,k,iBlock)*StateImpl_GV(i,j,k,:)

       State_VGB(Eradiation_,i,j,k,iBlock) = Eint_VCB(TradImpl_,i,j,k,iBlock)

       Einternal = Eint_VCB(TeImpl_,i,j,k,iBlock)
       EinternalSi = Einternal*No2Si_V(UnitEnergyDens_)

       call user_material_properties(State_VGB(:,i,j,k,iBlock), &
            EinternalSiIn = EinternalSi, GammaOut = Gamma)

       State_VGB(p_,i,j,k,iBlock) = Einternal*(Gamma-1.0)

       State_VGB(EintExtra_,i,j,k,iBlock) = &
            (1.0/(Gamma-1.0) - inv_gm1)*State_VGB(p_,i,j,k,iBlock)

    end do; end do; end do

    call calc_energy_cell(iBlock)

  end subroutine update_impl_gray_diff

end module ModGrayDiffusion
