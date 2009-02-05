!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
!============================================================================
module ModGrayDiffusion

  use ModVarIndexes, ONLY: p_

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
  public :: get_gray_diffusion_rhs
  public :: get_gray_diff_jacobian
  public :: update_impl_gray_diff

  ! Logical for adding Gray Diffusion
  logical, public :: IsNewBlockGrayDiffusion = .true.
  logical, public :: IsNewTimestepGrayDiffusion = .true.

  ! Coefficients for two-temperature electron-radiation model
  real, allocatable, public :: DiffusionRad_FDB(:,:,:,:,:)
  real, allocatable         :: RelaxationCoef_CB(:,:,:,:)
  real, allocatable         :: RelaxSemiCoef_VCB(:,:,:,:,:)
  real, allocatable         :: DiffSemiCoef_VGB(:,:,:,:,:)

  ! radiation energy used for calculating radiative energy flux
  real, allocatable :: Erad_G(:,:,:)

  real, parameter :: GammaRel = 4.0/3.0

  ! local index for extra internal energy to keep the compiler happy
  integer, parameter :: EintExtra_ = p_ - 1

  ! variables needed if the semi-implicit gray-diffusion uses the
  ! generic implicit solver
  integer, parameter :: EradImpl_ = 1, EintImpl_ = 2, Planck_ = 2

  integer :: nDiffusion = 1 ! for now

  real, parameter :: TradMinSi = 300.0 ! K
  real            :: EradMin

contains

  !============================================================================

  subroutine init_gray_diffusion

    use ModAdvance,    ONLY: Eradiation_
    use ModProcMH,     ONLY: iProc
    use ModSize,       ONLY: nI, nJ, nK, MaxBlock, nDim
    use ModVarIndexes, ONLY: NameVar_V
    use ModImplicit,   ONLY: UseSemiImplicit
    use ModPhysics,    ONLY: Si2No_V, UnitTemperature_, cRadiationNo

    character(len=*), parameter :: NameSub = "init_gray_diffusion"
    !------------------------------------------------------------------------

    if(allocated(Erad_G)) RETURN

    ! Make sure that Eradiation_ is correct
    if(NameVar_V(Eradiation_) /= "Erad") call stop_mpi(NameSub// &
         ": incorrect index for Erad variable in ModEquation")

    allocate( &
         Erad_G(0:nI+1,0:nJ+1,0:nK+1), &
         DiffusionRad_FDB(1:nI+1,1:nJ+1,1:nK+1,nDim,MaxBlock) )
       
    allocate(RelaxationCoef_CB(1:nI,1:nJ,1:nK,MaxBlock))
       
    if(UseSemiImplicit)then
       allocate( &
            RelaxSemiCoef_VCB(2,nI,nJ,nK,MaxBlock), &
            DiffSemiCoef_VGB(nDiffusion,-1:nI+2,-1:nJ+2,-1:nK+2,MaxBlock))
       RelaxSemiCoef_VCB = 0.0
       DiffSemiCoef_VGB = 0.0
    end if

    EradMin = cRadiationNo*(TradMinSi*Si2No_V(UnitTemperature_))**4

  end subroutine init_gray_diffusion

  !==========================================================================

  subroutine get_radiation_energy_flux( &
       iDir, i, j, k, iBlock, State_V, EradFlux_D)

    !\
    ! Calculate the diffusion part of the radiation energy flux.
    !/
    use ModAdvance,    ONLY: State_VGB, Eradiation_
    use ModPhysics,    ONLY: Si2No_V, UnitX_, UnitU_, cRadiationNo, Clight
    use ModTemperature,ONLY: UseRadFluxLimiter, TypeRadFluxLimiter
    use ModUser,       ONLY: user_material_properties
    use ModVarIndexes, ONLY: nVar

    integer, intent(in) :: iDir, i, j, k, iBlock
    real,    intent(in) :: State_V(nVar)
    real,    intent(out):: EradFlux_D(3)

    real :: RosselandMeanOpacitySi, DiffRad
    real :: RosselandMeanOpacity
    real :: FaceGrad_D(3), Erad, Grad2ByErad2
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
             DiffRad = Clight &
                  /(3.0*RosselandMeanOpacity + sqrt(Grad2ByErad2))
          case("max")
             DiffRad = Clight &
                  /max(3.0*RosselandMeanOpacity,sqrt(Grad2ByErad2))
          case("larsen")
             DiffRad = Clight &
                  /sqrt((3.0*RosselandMeanOpacity)**2 + Grad2ByErad2)
          end select
       else
          DiffRad = Clight/(3.0*RosselandMeanOpacity)
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
    use ModSize,     ONLY: nI, nJ, nK

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
    use ModConst,      ONLY: cLightSpeed
    use ModGeometry,   ONLY: vInv_CB, y_BLK, TypeGeometry
    use ModImplicit,   ONLY: UseFullImplicit
    use ModNodes,      ONLY: NodeY_NB
    use ModPhysics,    ONLY: cRadiationNo, Si2No_V, UnitTemperature_, UnitT_
    use ModMain,       ONLY: nI, nJ, nK, UseGrayDiffusion
    use ModUser,       ONLY: user_material_properties
    use ModVarIndexes, ONLY: Energy_, RhoUy_

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: TeSi, Te, vInv, DivU
    real :: RadCompression, AbsorptionEmission, PlanckOpacitySi
    character(len=19), parameter:: NameSub = "calc_source_gray_diffusion"
    !------------------------------------------------------------------------

    do k=1,nK; do j=1,nJ; do i=1,nI

      
       vInv = vInv_CB(i,j,k,iBlock)
       DivU = uDotArea_XI(i+1,j,k,1) - uDotArea_XI(i,j,k,1) &
            + uDotArea_YI(i,j+1,k,1) - uDotArea_YI(i,j,k,1) &
            + uDotArea_ZI(i,j,k+1,1) - uDotArea_ZI(i,j,k,1) 
       

       ! Adiabatic compression of radiation energy by fluid velocity (fluid 1)
       ! (GammaRel-1)*Erad*Div(U)
       RadCompression = (GammaRel-1.0)*State_VGB(Eradiation_,i,j,k,iBlock) &
            *vInv*DivU

       ! dErad/dt = - adiabatic compression
       Source_VC(Eradiation_,i,j,k) = Source_VC(Eradiation_,i,j,k) &
            - RadCompression

       ! dE/dt    = + adiabatic compression
       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + RadCompression

       if(IsNewTimestepGrayDiffusion)then
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               AbsorptionOpacitySiOut = PlanckOpacitySi)

          RelaxationCoef_CB(i,j,k,iBlock) = &
               PlanckOpacitySi*cLightSpeed/Si2No_V(UnitT_)
       end if

       if(.not.UseFullImplicit) CYCLE

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

  subroutine get_impl_gray_diff_state(StateImpl_VGB)

    use ModAdvance,  ONLY: Eradiation_, State_VGB
    use ModImplicit, ONLY: nw, nImplBlk, impl2iBlk, TypeSemiImplicit, kr
    use ModMain,     ONLY: nI, nJ, nK, MaxImplBlk, dt, &
         iTest, jTest, kTest, BlkTest
    use ModPhysics,  ONLY: inv_gm1, Clight, cRadiationNo, &
         Si2No_V, UnitTemperature_, UnitEnergyDens_, UnitX_
    use ModUser,     ONLY: user_material_properties
    use ModTemperature, ONLY: UseRadFluxLimiter, TypeRadFluxLimiter
    use ModGeometry, ONLY: dx_BLK, dy_BLK, dz_BLK
    use ModParallel, ONLY: NOBLK, NeiLev

    real, intent(out) :: StateImpl_VGB(nw,0:nI+1,0:nJ+1,0:nK+1,MaxImplBlk)

    integer :: iImplBlock, iBlock, i, j, k
    real :: PlanckOpacitySi, PlanckOpacity, CvSi, Cv, TeSi, Te
    real :: RosselandMeanOpacitySi, RosselandMeanOpacity
    real :: Grad2ByErad2, DiffRad, InvDx2, InvDy2, InvDz2
    !--------------------------------------------------------------------------

    do iImplBlock = 1, nImplBLK

       iBlock = impl2iBLK(iImplBlock)
       do k = 0, nK+1; do j = 0, nJ+1; do i = 0, nI+1
          StateImpl_VGB(EradImpl_,i,j,k,iImplBlock) = &
               State_VGB(Eradiation_,i,j,k,iBlock)
          
          if(nw == 1) CYCLE

          StateImpl_VGB(EintImpl_,i,j,k,iImplBlock) = &
               inv_gm1*State_VGB(p_,i,j,k,iBlock) &
               +  State_VGB(EintExtra_,i,j,k,iBlock)
       end do; end do; end do

       InvDx2 = 0.5/dx_BLK(iBlock)
       InvDy2 = 0.5/dy_BLK(iBlock)
       InvDz2 = 0.5/dz_BLK(iBlock)

       ! calculate coefficients for linearized energy exchange and diffusion
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               AbsorptionOpacitySiOut = PlanckOpacitySi, &
               RosselandMeanOpacitySiOut = RosselandMeanOpacitySi, &
               CvSiOut = CvSi, TeSiOut = TeSi)

          PlanckOpacity = PlanckOpacitySi/Si2No_V(UnitX_)
          RosselandMeanOpacity = RosselandMeanOpacitySi/Si2No_V(UnitX_)
          Cv = CvSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
          Te = TeSi*Si2No_V(UnitTemperature_)

          if(nw==1)then
             ! This coefficient is cR'' = cR/(1+dt*cR*dPlanck/dEint)
             RelaxSemiCoef_VCB(1,i,j,k,iBlock) = Clight*PlanckOpacity  &
                  /(1 + dt*Clight*PlanckOpacity*4.0*cRadiationNo*Te**3 / Cv)

             ! This is just the Planck function at time level * saved
             RelaxSemiCoef_VCB(Planck_,i,j,k,iBlock) = cRadiationNo*Te**4
          else
             RelaxSemiCoef_VCB(EradImpl_,i,j,k,iBlock) = Clight*PlanckOpacity
             RelaxSemiCoef_VCB(EintImpl_,i,j,k,iBlock) = Clight*PlanckOpacity &
                  * 4.0*cRadiationNo*Te**3 / Cv
          end if

          ! Calculate the cell centered diffusion coefficients
          if(UseRadFluxLimiter)then
             Grad2ByErad2 = &
                  (((State_VGB(Eradiation_,i+1,j,k,iBlock) &
                  -  State_VGB(Eradiation_,i-1,j,k,iBlock))*InvDx2)**2 &
                  +((State_VGB(Eradiation_,i,j+1,k,iBlock) &
                  -  State_VGB(Eradiation_,i,j-1,k,iBlock))*InvDy2)**2 &
                  +((State_VGB(Eradiation_,i,j,k+1,iBlock) &
                  -  State_VGB(Eradiation_,i,j,k-1,iBlock))*InvDz2)**2 &
                  )/  State_VGB(Eradiation_,i,j,k,iBlock)**2
             select case(TypeRadFluxLimiter)
             case("sum")
                DiffRad = Clight/(3.0*RosselandMeanOpacity + sqrt(Grad2ByErad2))
             case("max")
                DiffRad = Clight/max(3.0*RosselandMeanOpacity,sqrt(Grad2ByErad2))
             case("larsen")
                DiffRad = Clight/sqrt((3.0*RosselandMeanOpacity)**2 + Grad2ByErad2)
             end select
          else
             DiffRad = Clight/(3.0*RosselandMeanOpacity)
          end if

          ! Since this will always be averaged to the face, divide by two right now
          DiffSemiCoef_VGB(EradImpl_,i,j,k,iBlock) = 0.5*DiffRad

          !if(i==iTest.and.j==jTest.and.k==kTest.and.iBlock==BlkTest)then
          !   write(*,*)'!!! State_V=',State_VGB(:,i,j,k,iBlock)
          !   write(*,*)'!!! UseRadFluxLimiter, TypeRadFluxLimiter=', &
          !        UseRadFluxLimiter, TypeRadFluxLimiter
          !   write(*,*)'!!! RosselandMeanOpacity=',RosselandMeanOpacity
          !   if(UseRadFluxLimiter)then
          !      write(*,*)'!!!InvDx2, Eradiation_,i-1:i+1)=', InvDx2, &
          !           State_VGB(Eradiation_,i-1:i+1,j,k,iBlock)
          !      write(*,*)'!!!InvDy2, Eradiation_,j-1:j+1)=', InvDy2, &
          !           State_VGB(Eradiation_,i,j-1:j+1,k,iBlock)
          !      write(*,*)'!!!InvDz2, Eradiation_,k-1:k+1)=', InvDz2, &
          !           State_VGB(Eradiation_,i,j,k-1:k+1,iBlock)
          !      write(*,*)'!!! Grad2ByErad2=',Grad2ByErad2
          !   end if
          !   write(*,*)'!!! DiffRad =',DiffRad
          !   
          !end if

       end do; end do; end do

       ! set floating outer boundary
       if(NeiLev(1,iBlock) == NOBLK) DiffSemiCoef_VGB(:,0   ,:,:,iBlock) &
            =                        DiffSemiCoef_VGB(:,1   ,:,:,iBlock)
       if(NeiLev(2,iBlock) == NOBLK) DiffSemiCoef_VGB(:,nI+1,:,:,iBlock) &
            =                        DiffSemiCoef_VGB(:,nI  ,:,:,iBlock)
       if(NeiLev(3,iBlock) == NOBLK) DiffSemiCoef_VGB(:,:,0   ,:,iBlock) &
            =                        DiffSemiCoef_VGB(:,:,1   ,:,iBlock)
       if(NeiLev(4,iBlock) == NOBLK) DiffSemiCoef_VGB(:,:,nJ+1,:,iBlock) &
            =                        DiffSemiCoef_VGB(:,:,nJ  ,:,iBlock)
       if(NeiLev(5,iBlock) == NOBLK) DiffSemiCoef_VGB(:,:,:,0   ,iBlock) &
            =                        DiffSemiCoef_VGB(:,:,:,1   ,iBlock)
       if(NeiLev(6,iBlock) == NOBLK) DiffSemiCoef_VGB(:,:,:,nK+1,iBlock) &
            =                        DiffSemiCoef_VGB(:,:,:,nK  ,iBlock)

    end do

    ! DoOneLayer, DoFacesOnly, No UseMonoteRestrict
    call message_pass_cells8(.true., .true., .false., nDiffusion, &
         DiffSemiCoef_VGB)

    !!! DiffusionRad_FDB = 0.0 !!!
    !!! RelaxSemiCoef_VCB = 0.0 !!!

  end subroutine get_impl_gray_diff_state

  !============================================================================

  subroutine get_gray_diffusion_rhs(iBlock, StateImpl_VG, Rhs_VC, IsLinear)

    use ModAdvance,  ONLY: Flux_VX, Flux_VY, Flux_VZ, State_VGB, Rho_
    use ModGeometry, ONLY: dx_BLK, dy_BLK, dz_BLK, y_Blk, &
         TypeGeometry, vInv_CB
    use ModImplicit, ONLY: nw
    use ModMain,     ONLY: nI, nJ, nK, x_, y_, z_, TypeBc_I, &
         iTest, jTest, kTest, BlkTest, Dt
    use ModNodes,    ONLY: NodeY_NB
    use ModParallel, ONLY: NOBLK, NeiLev
    use ModPhysics,  ONLY: gm1, cRadiationNo, No2Si_V, Si2No_V, &
         UnitEnergyDens_, UnitTemperature_
    use ModUser,     ONLY: user_material_properties

    integer, intent(in) :: iBlock
    real, intent(inout) :: StateImpl_VG(nw,-1:nI+2,-1:nJ+2,-1:nK+2)
    real, intent(out)   :: Rhs_VC(nw,nI,nJ,nK)
    logical, intent(in) :: IsLinear

    real :: InvDx2, InvDy2, InvDz2
    real :: EinternalSi, TeSi, Te, AbsorptionEmission
    integer :: i, j, k, iVar
    !--------------------------------------------------------------------------

    if(NeiLev(1,iBlock) == NOBLK)then 
       if(IsLinear .and. TypeBc_I(1) /= 'reflect')then
          StateImpl_VG(:,0   ,:,:) = 0.0
       else
          StateImpl_VG(:,0   ,:,:) = StateImpl_VG(:,1   ,:,:)
       end if
    end if
    if(NeiLev(2,iBlock) == NOBLK)then
       if(IsLinear .and. TypeBc_I(2) /= 'reflect')then
          StateImpl_VG(:,nI+1,:,:) = 0.0
       else
          StateImpl_VG(:,nI+1,:,:) = StateImpl_VG(:,nI  ,:,:)
       end if
    end if
    if(NeiLev(3,iBlock) == NOBLK)then
       if(IsLinear .and. TypeBc_I(3) /= 'reflect')then
          StateImpl_VG(:,:,0   ,:) = 0.0
       else
          StateImpl_VG(:,:,0   ,:) = StateImpl_VG(:,:,1   ,:)
       end if
    end if
    if(NeiLev(4,iBlock) == NOBLK) then
       if(IsLinear .and. TypeBc_I(4) /= 'reflect')then
          StateImpl_VG(:,:,nJ+1,:) = 0.0
       else
          StateImpl_VG(:,:,nJ+1,:) = StateImpl_VG(:,:,nJ  ,:)
       end if
    end if
    if(NeiLev(5,iBlock) == NOBLK) then
       if(IsLinear .and. TypeBc_I(5) /= 'reflect')then
          StateImpl_VG(:,:,:,0   ) = 0.0
       else
          StateImpl_VG(:,:,:,0   ) = StateImpl_VG(:,:,:,1   )
       end if
    end if
    if(NeiLev(6,iBlock) == NOBLK)then
       if(IsLinear .and. TypeBc_I(6) /= 'reflect')then
          StateImpl_VG(:,:,:,nK+1) = 0.0
       else
          StateImpl_VG(:,:,:,nK+1) = StateImpl_VG(:,:,:,nK  )
       end if
    end if

    ! Calculate radiative diffusion fluxes
    InvDx2 = 1./dx_BLK(iBlock)**2
    InvDy2 = 1./dy_BLK(iBlock)**2
    InvDz2 = 1./dz_BLK(iBlock)**2
    do k = 1, nK; do j = 1, nJ; do i = 1, nI+1
       Flux_VX(EradImpl_,i,j,k) = InvDx2 &
            *sum(DiffSemiCoef_VGB(EradImpl_,i-1:i,j,k,iBlock)) &
            *(StateImpl_VG(EradImpl_,i,j,k) - StateImpl_VG(EradImpl_,i-1,j,k))
    end do; end do; end do
    do k = 1, nK; do j = 1, nJ+1; do i = 1, nI
       Flux_VY(EradImpl_,i,j,k) = InvDy2 &
            *sum(DiffSemiCoef_VGB(EradImpl_,i,j-1:j,k,iBlock)) &
            *(StateImpl_VG(EradImpl_,i,j,k) - StateImpl_VG(EradImpl_,i,j-1,k))
    end do; end do; end do
    if(.not.TypeGeometry=='rz')then
       do k = 1, nK+1; do j = 1, nJ; do i = 1, nI
          Flux_VZ(EradImpl_,i,j,k) = InvDz2 &
               *sum(DiffSemiCoef_VGB(EradImpl_,i,j,k-1:k,iBlock)) &
               *(StateImpl_VG(EradImpl_,i,j,k) &
               - StateImpl_VG(EradImpl_,i,j,k-1))
       end do; end do; end do
    end if

    if(TypeGeometry=='rz')then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          Rhs_VC(EradImpl_,i,j,k) = &
               + Flux_VX(EradImpl_,i+1,j,k) - Flux_VX(EradImpl_,i,j,k) &
               + ( Flux_VY(EradImpl_,i,j+1,k)*abs(NodeY_NB(i,j+1,k,iBlock)) &
               -   Flux_VY(EradImpl_,i,j,k)*abs(NodeY_NB(i,j,k,iBlock)) ) &
               /   abs(y_Blk(i,j,k,iBlock))
       end do; end do; end do
    else
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          Rhs_VC(EradImpl_,i,j,k) = &
               + Flux_VX(EradImpl_,i+1,j,k) - Flux_VX(EradImpl_,i,j,k) &
               + Flux_VY(EradImpl_,i,j+1,k) - Flux_VY(EradImpl_,i,j,k) &
               + Flux_VZ(EradImpl_,i,j,k+1) - Flux_VZ(EradImpl_,i,j,k)
       end do; end do; end do
    end if

    if(.false. .and. iBlock==BlkTest)then
       i=iTest; j=jTest; k=kTest; iVar=EradImpl_
       write(*,*)'IsLinear=',IsLinear

       write(*,*)'Dt     =',Dt
       write(*,*)'Volume =',1/vInv_CB(iTest,jTest,kTest,iBlock)

       write(*,*)'!!! Flux(i+1,i)=',Flux_VX(iVar,i+1,j,k), Flux_VX(iVar,i,j,k)
       write(*,*)'!!! Flux(j+1,j)=',Flux_VY(iVar,i,j+1,k), Flux_VY(iVar,i,j,k)
       write(*,*)'!!! Flux(k+1,k)=',Flux_VZ(iVar,i,j,k+1), Flux_VZ(iVar,i,j,k)
       
       write(*,*)'!!! DiffCoef(j-1,j)=',2*DiffSemiCoef_VGB(iVar,i,j-1:j,k,iBlock)
       write(*,*)'!!! StateImp(j-1,j)=',StateImpl_VG(iVar,i,j-1:j,k)

       write(*,*)'!!! StateImpl(test) =',StateImpl_VG(iVar,i,j,k)
       write(*,*)'!!! Rhs(test)       =',Rhs_VC(iVar,i,j,k)
       write(*,*)'!!! sum(StateImpl)  =',sum(StateImpl_VG(iVar,1:nI,1:nJ,1:nK))
       write(*,*)'!!! sum(Rhs)        =',sum(Rhs_VC(iVar,:,:,:))

    end if

    ! Source term due to absorption and emission
    ! Sigma_a*(cRadiation*Te**4-Erad)
    if(IsLinear)then
       if(nw==1)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             Rhs_VC(EradImpl_,i,j,k) = Rhs_VC(EradImpl_,i,j,k) &
                  - RelaxSemiCoef_VCB(EradImpl_,i,j,k,iBlock) &
                  *StateImpl_VG(EradImpl_,i,j,k)
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             Rhs_VC(EintImpl_,i,j,k) = &
                  - RelaxSemiCoef_VCB(EintImpl_,i,j,k,iBlock) &
                  *StateImpl_VG(EintImpl_,i,j,k) &
                  + RelaxSemiCoef_VCB(EradImpl_,i,j,k,iBlock) &
                  *StateImpl_VG(EradImpl_,i,j,k)
             Rhs_VC(EradImpl_,i,j,k) = Rhs_VC(EradImpl_,i,j,k) &
                  - Rhs_VC(EintImpl_,i,j,k) 
          end do; end do; end do
       end if
    else
       if(nw==1)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             AbsorptionEmission = RelaxSemiCoef_VCB(EradImpl_,i,j,k,iBlock) &
               * (RelaxSemiCoef_VCB(Planck_,i,j,k,iBlock) &
               -  StateImpl_VG(EradImpl_,i,j,k))

             Rhs_VC(EradImpl_,i,j,k) = Rhs_VC(EradImpl_,i,j,k) + AbsorptionEmission

          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             EinternalSi = StateImpl_VG(EintImpl_,i,j,k)*No2Si_V(UnitEnergyDens_)
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  EinternalSiIn = EinternalSi, TeSiOut = TeSi)
             Te = TeSi*Si2No_V(UnitTemperature_)

             AbsorptionEmission = RelaxSemiCoef_VCB(EradImpl_,i,j,k,iBlock)&
                  *(cRadiationNo*Te**4 - StateImpl_VG(EradImpl_,i,j,k))

             ! dErad/dt = + AbsorptionEmission
             Rhs_VC(EradImpl_,i,j,k) = Rhs_VC(EradImpl_,i,j,k) &
                  + AbsorptionEmission

             ! dE/dt    = - AbsorptionEmission
             Rhs_VC(EintImpl_,i,j,k) = &
                  - AbsorptionEmission
          end do; end do; end do
       end if
    end if

  end subroutine get_gray_diffusion_rhs

  !===========================================================================

  subroutine get_gray_diff_jacobian(iBlock, nVar, Jacobian_VVCI)

    use ModAdvance,  ONLY: State_VGB
    use ModGeometry, ONLY: Dx_Blk, Dy_Blk, Dz_Blk, &
         vInv_CB, y_Blk, TypeGeometry
    use ModImplicit, ONLY: kr, nw
    use ModMain,     ONLY: nI, nJ, nK, nDim, iTest, jTest, kTest, BlkTest
    use ModNodes,    ONLY: NodeY_NB

    integer, parameter:: nStencil = 2*nDim + 1

    integer, intent(in) :: iBlock, nVar
    real, intent(out) :: Jacobian_VVCI(nVar,nVar,nI,nJ,nK,nStencil)

    integer :: iVar, i, j, k, iDim, Di, Dj, Dk, iCond
    real :: DiffLeft, DiffRight, Dxyz_D(nDim), Coeff
    
    real :: dSdErad, dSdEint
    !--------------------------------------------------------------------------

    ! All elements have to be set
    Jacobian_VVCI(:,:,:,:,:,:) = 0.0

    do k = 1, nK; do j = 1, nJ; do i = 1, nI

       ! energy exchange
       dSdErad = RelaxSemiCoef_VCB(EradImpl_,i,j,k,iBlock)

       ! dSrad/dErad (diagonal)
       Jacobian_VVCI(EradImpl_,EradImpl_,i,j,k,1) = -dSdErad

       if(nw == 1) CYCLE

       dSdEint = RelaxSemiCoef_VCB(EintImpl_,i,j,k,iBlock)

       ! dSint/dErad (off diagonal)
       Jacobian_VVCI(EintImpl_,EradImpl_,i,j,k,1) = +dSdErad

       ! dSint/dEint (diagonal)
       Jacobian_VVCI(EintImpl_,EintImpl_,i,j,k,1) = -dSdEint

       ! dSrad/dEint (off diagonal)
       Jacobian_VVCI(EradImpl_,EintImpl_,i,j,k,1) = +dSdEint

    end do; end do; end do

    Dxyz_D = (/dx_BLK(iBlock), dy_BLK(iBlock), dz_Blk(iBlock)/)
    iVar = EradImpl_
    do iDim = 1, nDim
       Coeff = 1.0/Dxyz_D(iDim)**2
       Di = kr(iDim,1)
       Dj = kr(iDim,2)
       Dk = kr(iDim,3)
       if(TypeGeometry=='rz'.and.iDim==2)then
          do k=1,nK; do j=1,nJ; do i=1,nI
             if(j==1)then
                DiffLeft = 0.0
             else
                DiffLeft = Coeff &
                     *sum(DiffSemiCoef_VGB(EradImpl_,i,j-1:j,k,iBlock)) &
                     *abs(NodeY_NB(i,j,k,iBlock)/y_Blk(i,j,k,iBlock))
             end if
             if(j==nJ)then
                DiffRight = 0.0
             else
                DiffRight = Coeff &
                     *sum(DiffSemiCoef_VGB(EradImpl_,i,j:j+1,k,iBlock)) &
                     *abs(NodeY_NB(i,j+1,k,iBlock)/y_Blk(i,j,k,iBlock))
             end if
             Jacobian_VVCI(iVar,iVar,i,j,k,1) = &
                  Jacobian_VVCI(iVar,iVar,i,j,k,1) - (DiffLeft + DiffRight)
             Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim)   = DiffLeft
             Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim+1) = DiffRight
          end do; end do; end do
          EXIT ! Done with cylindrical
       end if
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(iDim==1.and.i==1 .or. iDim==2.and.j==1 .or. iDim==3.and.k==1)then
             DiffLeft = 0.0
          else
             DiffLeft = Coeff &
                  *sum(DiffSemiCoef_VGB(EradImpl_,i-Di:i,j-Dj:j,k-Dk:k,iBlock))
          end if
          if(iDim==1.and.i==nI .or. iDim==2.and.j==nJ .or. &
               iDim==3.and.k==nK) then
             DiffRight = 0.0
          else
             DiffRight = Coeff &
                  *sum(DiffSemiCoef_VGB(EradImpl_,i:i+Di,j:j+Dj,k:k+Dk,iBlock))
          end if
          Jacobian_VVCI(iVar,iVar,i,j,k,1) = &
               Jacobian_VVCI(iVar,iVar,i,j,k,1) - (DiffLeft + DiffRight)
          Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim)   = DiffLeft
          Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim+1) = DiffRight
       end do; end do; end do
    end do

  end subroutine get_gray_diff_jacobian

  !============================================================================

  subroutine update_impl_gray_diff(iBlock, StateImpl_VG)

    use ModAdvance,  ONLY: State_VGB, p_, Eradiation_
    use ModEnergy,   ONLY: calc_energy_cell
    use ModImplicit, ONLY: nw
    use ModMain,     ONLY: nI, nJ, nK, dt
    use ModPhysics,  ONLY: inv_gm1, No2Si_V, Si2No_V, UnitEnergyDens_, UnitP_
    use ModUser,     ONLY: user_material_properties

    integer, intent(in) :: iBlock
    real, intent(in) :: StateImpl_VG(nw,nI,nJ,nK)

    integer :: i, j, k
    real :: Einternal, EinternalSi, PressureSi, AbsorptionEmission

    character(len=*), parameter :: NameSub = 'update_impl_gray_diff'
    !--------------------------------------------------------------------------

    do k = 1,nK; do j = 1,nJ; do i = 1,nI
       State_VGB(Eradiation_,i,j,k,iBlock) = &
            max(EradMin, StateImpl_VG(EradImpl_,i,j,k))

       if(nw == 1)then
          AbsorptionEmission = RelaxSemiCoef_VCB(EradImpl_,i,j,k,iBlock) &
               * (RelaxSemiCoef_VCB(Planck_,i,j,k,iBlock) &
               -  State_VGB(Eradiation_,i,j,k,iBlock))
          Einternal = inv_gm1*State_VGB(p_,i,j,k,iBlock) &
               +  State_VGB(EintExtra_,i,j,k,iBlock) &
               - dt*AbsorptionEmission
       else
          Einternal = StateImpl_VG(EintImpl_,i,j,k)
       end if
       EinternalSi = Einternal*No2Si_V(UnitEnergyDens_)

       if(State_VGB(Eradiation_,i,j,k,iBlock) < 0.0 .or. Einternal < 0.0)then
          write(*,*)NameSub,': ERROR negative EradOrig, Erad, or Eint=', &
               StateImpl_VG(EradImpl_,i,j,k), &
               State_VGB(Eradiation_,i,j,k,iBlock), Einternal
          write(*,*)NameSub,': ERROR at i,j,k,iBlock=', i, j, k, iBlock
          call stop_mpi(NameSub//' negative Erad')
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
