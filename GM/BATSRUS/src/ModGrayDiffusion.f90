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
  public :: advance_temperature

  ! Logical for adding Gray Diffusion
  logical, public :: IsNewBlockGrayDiffusion = .true.
  logical, public :: IsNewTimestepGrayDiffusion = .true.

  ! Parameters for radiation flux limiter
  logical,           public :: UseRadFluxLimiter  = .false.
  character(len=20), public :: TypeRadFluxLimiter = 'larsen'

  ! Coefficients for two-temperature electron-radiation model
  real, allocatable, public :: DiffusionRad_FDB(:,:,:,:,:)
  real, allocatable         :: RelaxationCoef_VCB(:,:,:,:,:)
  real, allocatable         :: SpecificHeat_VCB(:,:,:,:,:)
  real, allocatable         :: HeatConductionCoef_IGB(:,:,:,:,:)

  ! Store block Jacobi preconditioner
  real, allocatable :: JacPrec_VVCB(:,:,:,:,:,:)

  ! For the purpose of message passing the heat conduction coefficients,
  ! they are compactly stored without gaps. A pointer is used to indicate
  ! which temperature variables involve heat conduction
  integer, allocatable :: iCond_I(:)

  ! number of variables that uses heat conduction
  ! currently, only radiation needs "heat conduction"
  integer, parameter :: nCond = 1

  ! change in internal energy due to heat conduction/diffusion
  ! and energy exchange
  real, allocatable :: DelEint_VCB(:,:,:,:,:)

  ! Right-hand side of non-linear coupled system of temperature equations
  real, allocatable :: Source_VCB(:,:,:,:,:)

  ! radiation energy used for calculating radiative energy flux
  real, allocatable :: Erad_G(:,:,:)

  real, parameter :: GammaRel = 4.0/3.0

  ! local index for extra internal energy to keep the compiler happy
  integer, parameter :: EintExtra_ = p_ - 1

  ! the solution vector containing the independent temperature variables
  real, allocatable :: Temperature_VGB(:,:,:,:,:)

  ! number of independent temperature variables
  integer, parameter :: nTemperature = 2

  ! named indexes for temerature variables
  integer, parameter :: Te_ = 1, Trad_ = 2

contains

  !==========================================================================

  subroutine init_gray_diffusion

    use ModAdvance,    ONLY: Eradiation_
    use ModImplicit,   ONLY: UseFullImplicit
    use ModProcMH,     ONLY: iProc
    use ModSize,       ONLY: nI, nJ, nK, nBlk, nDim
    use ModVarIndexes, ONLY: NameVar_V

    character(len=*), parameter :: NameSub = "init_gray_diffusion"
    !------------------------------------------------------------------------

    ! Make sure that Eradiation_ is correct
    if(NameVar_V(Eradiation_) /= "Erad") call stop_mpi(NameSub// &
         ": incorrect index for Erad variable in ModEquation")


    if(UseFullImplicit)then

       if(.not.allocated(Erad_G)) then
          allocate( &
               Erad_G(0:nI+1,0:nJ+1,0:nK+1), &
               DiffusionRad_FDB(1:nI+1,1:nJ+1,1:nK+1,nDim,nBlk) )

          if(nTemperature>1) &
               allocate(RelaxationCoef_VCB(2:nTemperature,1:nI,1:nJ,1:nK,nBlk))
       
       end if

       RETURN
    end if


    ! Initializing the semi-implicit gray-diffusion

    if(iProc == 0)then
       write(*,*) "==================================================="
       write(*,*) "WARNING !!!"
       write(*,*) "The semi-implicit gray-diffusion is not yet working"
       write(*,*) "==================================================="
    end if

    if(.not.allocated(SpecificHeat_VCB))then
       allocate( &
            SpecificHeat_VCB(1:nTemperature,1:nI,1:nJ,1:nK,nBlk), &
            HeatConductionCoef_IGB(1:nCond,-1:nI+2,-1:nJ+2,-1:nK+2,nBlk), &
            iCond_I(1:nCond), &
            Temperature_VGB(1:nTemperature,-1:nI+2,-1:nJ+2,-1:nK+2,nBlk), &
            JacPrec_VVCB(nTemperature,nTemperature,nI,nJ,nK,nBLK), &
            DelEint_VCB(1:nTemperature,1:nI,1:nJ,1:nK,nBlk), &
            Source_VCB(1:nTemperature,1:nI,1:nJ,1:nK,nBlk) )

       if(nTemperature>1) &
            allocate(RelaxationCoef_VCB(2:nTemperature,1:nI,1:nJ,1:nK,nBlk))

    end if

    ! Radiation has diffusion now. 
    ! Heat conduction for electrons to be done later !!!
    iCond_I(1) = Trad_

  end subroutine init_gray_diffusion

  !==========================================================================

  subroutine set_frozen_coefficients

    use ModAdvance,  ONLY: State_VGB, Eradiation_
    use ModConst,    ONLY: cLightSpeed
    use ModMain,     ONLY: nBlock, unusedBlk
    use ModPhysics,  ONLY: Si2No_V, UnitT_, UnitX_, UnitU_, &
         UnitEnergyDens_, UnitTemperature_, cRadiationNo
    use ModSize,     ONLY: nI, nJ, nK
    use ModUser,     ONLY: user_material_properties

    integer :: i, j, k, iBlock
    real :: PlanckOpacitySi, RosselandMeanOpacitySi
    real :: CvSi, TeSi, Te, Trad, DiffRad
    real :: Grad_D(3), Grad2ByErad2
    !------------------------------------------------------------------------

    do iBlock = 1, nBlock
       if(unusedBlk(iBlock)) CYCLE

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               TeSiOut = TeSi)

          Te = TeSi*Si2No_V(UnitTemperature_)
          Trad = sqrt(sqrt(State_VGB(Eradiation_,i,j,k,iBlock)/cRadiationNo))

          Temperature_VGB(Te_,i,j,k,iBlock) = Te
          Temperature_VGB(Trad_,i,j,k,iBlock) = Trad
       end do; end do; end do
    end do

    if(UseRadFluxLimiter)then
       ! The radiation flux limiters use the gradient of the
       ! radiation temperature and therefor one layer of ghost cells
       ! needs to be correct. We massage pass all temperatures.

       ! DoOneLayer, DoFacesOnly, No UseMonoteRestrict
       call message_pass_cells8(.true., .true., .false., nTemperature, &
            Temperature_VGB)
    end if

    if(UseRadFluxLimiter) &
         call set_gray_outer_bcs(nTemperature, Temperature_VGB)

    do iBlock = 1, nBlock
       if(unusedBlk(iBlock)) CYCLE

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               AbsorptionOpacitySiOut = PlanckOpacitySi, &
               RosselandMeanOpacitySiOut = RosselandMeanOpacitySi, &
               CvSiOut = CvSi)

          Te = Temperature_VGB(Te_,i,j,k,iBlock)
          Trad = Temperature_VGB(Trad_,i,j,k,iBlock)

          SpecificHeat_VCB(Te_,i,j,k,iBlock) = CvSi &
               *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)

          SpecificHeat_VCB(Trad_,i,j,k,iBlock) = 4.0*cRadiationNo*Trad**3

          RelaxationCoef_VCB(Trad_,i,j,k,iBlock) = &
               PlanckOpacitySi*cLightSpeed/Si2No_V(UnitT_) &
               *cRadiationNo*(Te+Trad)*(Te**2+Trad**2)

          call get_heat_conduction_coef
       end do; end do; end do
    end do

    ! message pass one ghost cell layer of heat conduction coefficients

    ! DoOneLayer, DoFacesOnly, No UseMonoteRestrict
    call message_pass_cells8(.true., .true., .false., nCond, &
         HeatConductionCoef_IGB)

    call set_gray_outer_bcs(nCond, HeatConductionCoef_IGB)

  contains

    subroutine get_heat_conduction_coef

      !------------------------------------------------------------------------

      DiffRad = cLightSpeed/(3.0*RosselandMeanOpacitySi) &
           *Si2No_V(UnitU_)*Si2No_V(UnitX_)

      if(UseRadFluxLimiter)then
         call calc_cell_gradient
         Grad2ByErad2 = 16.0*sum(Grad_D**2)/Trad**2

         select case(TypeRadFluxLimiter)
         case("sum")
            DiffRad = 1.0/(1.0/DiffRad+sqrt(Grad2ByErad2))
         case("max")
            DiffRad = 1.0/max(1.0/DiffRad,sqrt(Grad2ByErad2))
         case("larsen")
            DiffRad = 1.0/sqrt(1.0/DiffRad**2+Grad2ByErad2)
         end select
      end if

      HeatConductionCoef_IGB(1,i,j,k,iBlock) = DiffRad &
           *4.0*cRadiationNo*Trad**3

    end subroutine get_heat_conduction_coef

    !==========================================================================

    subroutine calc_cell_gradient

      ! currently, only a cartesian version

      use ModGeometry, ONLY: Dx_Blk, Dy_Blk, Dz_Blk
      use ModMain,     ONLY: x_, y_, z_

      real :: InvDx2, InvDy2, InvDz2
      !------------------------------------------------------------------------

      InvDx2 = 0.5/Dx_Blk(iBlock)
      InvDy2 = 0.5/Dy_Blk(iBlock)
      InvDz2 = 0.5/Dz_Blk(iBlock)

      Grad_D(x_) = &
           ( Temperature_VGB(Trad_,i+1,j,  k,  iBlock) &
           - Temperature_VGB(Trad_,i-1,j,  k,  iBlock) )*InvDx2
      Grad_D(y_) = &
           ( Temperature_VGB(Trad_,i,  j+1,k,  iBlock) &
           - Temperature_VGB(Trad_,i,  j-1,k,  iBlock) )*InvDy2
      Grad_D(z_) = &
           ( Temperature_VGB(Trad_,i,  j,  k+1,iBlock) &
           - Temperature_VGB(Trad_,i,  j,  k-1,iBlock) )*InvDz2

    end subroutine calc_cell_gradient

  end subroutine set_frozen_coefficients

  !============================================================================

  subroutine get_radiation_energy_flux( &
       iDir, i, j, k, iBlock, State_V, EradFlux_D)

    !\
    ! Calculate the diffusion part of the radiation energy flux.
    !/
    use ModAdvance,    ONLY: State_VGB, Eradiation_
    use ModConst,      ONLY: cLightSpeed
    use ModPhysics,    ONLY: Si2No_V, UnitX_, UnitU_, cRadiationNo
    use ModUser,       ONLY: user_material_properties
    use ModVarIndexes, ONLY: nVar

    integer, intent(in) :: iDir, i, j, k, iBlock
    real,    intent(in) :: State_V(nVar)
    real,    intent(out):: EradFlux_D(3)

    real :: RosselandMeanOpacitySi, DiffRad
    real :: FaceGrad_D(3), Erad, Grad2ByErad2
    !------------------------------------------------------------------------

    call calc_face_gradient(iDir, i, j, k, iBlock, FaceGrad_D)


    if(IsNewTimestepGrayDiffusion)then

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
    use ModGeometry,   ONLY: vInv_CB, y_BLK
    use ModImplicit,   ONLY: UseFullImplicit
    use ModMain,       ONLY: IsCylindrical
    use ModPhysics,    ONLY: cRadiationNo, Si2No_V, UnitTemperature_, UnitT_
    use ModSize,       ONLY: nI, nJ, nK
    use ModUser,       ONLY: user_material_properties
    use ModVarIndexes, ONLY: Energy_

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: TeSi, Te, vInv
    real :: RadCompression, AbsorptionEmission, PlanckOpacitySi
    character(len=19), parameter:: NameSub = "calc_source_gray_diffusion"
    !------------------------------------------------------------------------

    do k=1,nK; do j=1,nJ; do i=1,nI

       if(IsCylindrical)then
          ! Multiply volume with radius (=Y) at cell center
          ! -> divide inverse volume
          vInv = vInv_CB(i,j,k,iBlock)/abs(y_BLK(i,j,k,iBlock))
       else
          vInv = vInv_CB(i,j,k,iBlock)
       end if

       ! Adiabatic compression of radiation energy by fluid velocity (fluid 1)
       ! (GammaRel-1)*Erad*Div(U)
       RadCompression = (GammaRel-1.0)*State_VGB(Eradiation_,i,j,k,iBlock) &
            *vInv &
            *(uDotArea_XI(i+1,j,k,1) - uDotArea_XI(i,j,k,1) &
             +uDotArea_YI(i,j+1,k,1) - uDotArea_YI(i,j,k,1) &
             +uDotArea_ZI(i,j,k+1,1) - uDotArea_ZI(i,j,k,1))

       ! dErad/dt = - adiabatic compression
       Source_VC(Eradiation_,i,j,k) = Source_VC(Eradiation_,i,j,k) &
            - RadCompression

       ! dE/dt    = + adiabatic compression
       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + RadCompression

       ! In the semi-implicit solver, the energy exchange is already added
       if(.not.UseFullImplicit) CYCLE

       if(IsNewTimestepGrayDiffusion)then
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               AbsorptionOpacitySiOut = PlanckOpacitySi)

          RelaxationCoef_VCB(Trad_,i,j,k,iBlock) = &
               PlanckOpacitySi*cLightSpeed/Si2No_V(UnitT_)
       end if

       call user_material_properties(State_VGB(:,i,j,k,iBlock), &
            TeSiOut = TeSi)

       Te = TeSi*Si2No_V(UnitTemperature_)

       ! Source term due to absorption and emission
       ! Sigma_a*(cRadiation*Te**4-Erad)
       AbsorptionEmission =  RelaxationCoef_VCB(Trad_,i,j,k,iBlock) &
            *(cRadiationNo*Te**4 - State_VGB(Eradiation_,i,j,k,iBlock))

       ! dErad/dt = + AbsorptionEmission
       Source_VC(Eradiation_,i,j,k) = Source_VC(Eradiation_,i,j,k) &
            + AbsorptionEmission

       ! dE/dt    = - AbsorptionEmission
       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
            - AbsorptionEmission

    end do; end do; end do

  end subroutine calc_source_gray_diffusion

  !============================================================================
  subroutine advance_temperature

    use ModGeometry, ONLY: vInv_CB, y_BLK
    use ModMain,     ONLY: nBlock, unusedBLK, Dt, IsCylindrical
    use ModSize,     ONLY: nI, nJ, nK

    integer :: i, j, k, iBlock
    integer :: Iter
    real, parameter :: Tolerance = 1e-5
    logical :: DoTest, DoTestMe

    character(len=*), parameter :: NameSub = 'advance_temperature'
    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    ! prevent doing useless updates
    if(Dt == 0.0)return

    call set_frozen_coefficients

    do iBlock = 1, nBlock
       if(unusedBlk(iBlock)) CYCLE

       Source_VCB(:,:,:,:,iBlock) = SpecificHeat_VCB(:,:,:,:,iBlock) &
            *Temperature_VGB(:,1:nI,1:nJ,1:nK,iBlock)
       DelEint_VCB(:,:,:,:,iBlock) = -Source_VCB(:,:,:,:,iBlock)

       do k = 1,nK; do j = 1,nJ; do i = 1,nI
          Source_VCB(:,i,j,k,iBlock) = Source_VCB(:,i,j,k,iBlock) &
               /(vInv_CB(i,j,k,iBlock)*Dt)
       end do; end do; end do

       if(IsCylindrical)then
          do k = 1,nK; do j = 1,nJ; do i = 1,nI
             Source_VCB(:,i,j,k,iBlock) = Source_VCB(:,i,j,k,iBlock) &
                  *abs(y_BLK(i,j,k,iBlock))
          end do; end do; end do
       end if

    end do

    ! do a preconditioned conjugate gradient
    call get_jacobi_preconditioner(nTemperature)

    call cg(heat_conduction, nTemperature, Source_VCB, Temperature_VGB, &
         .true., Tolerance, 'rel', Iter, jacobi_preconditioner)

    if(DoTestMe)write(*,*)NameSub,': Number of iterations =',Iter

    do iBlock = 1, nBlock
       if(unusedBlk(iBlock)) CYCLE

       DelEint_VCB(:,:,:,:,iBlock) = DelEint_VCB(:,:,:,:,iBlock) &
            + SpecificHeat_VCB(:,:,:,:,iBlock) &
            *Temperature_VGB(:,1:nI,1:nJ,1:nK,iBlock)

       call update_conservative_energy(iBlock)
    end do

  end subroutine advance_temperature

  !==========================================================================

  subroutine update_conservative_energy(iBlock)

    use ModAdvance,  ONLY: State_VGB, p_, Eradiation_
    use ModEnergy,   ONLY: calc_energy_cell
    use ModPhysics,  ONLY: inv_gm1, No2Si_V, Si2No_V, UnitEnergyDens_, UnitP_
    use ModSize,     ONLY: nI, nJ, nK
    use ModUser,     ONLY: user_material_properties

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: EinternalSi, Einternal, Gamma, PressureSi
    !------------------------------------------------------------------------

    do k = 1,nK; do j = 1,nJ; do i = 1,nI

       State_VGB(Eradiation_,i,j,k,iBlock) = &
            State_VGB(Eradiation_,i,j,k,iBlock) &
            + DelEint_VCB(Trad_,i,j,k,iBlock)

       Einternal = inv_gm1*State_VGB(p_,i,j,k,iBlock) &
            + State_VGB(EintExtra_,i,j,k,iBlock) &
            + DelEint_VCB(Te_,i,j,k,iBlock)

       EinternalSi = Einternal*No2Si_V(UnitEnergyDens_)

       call user_material_properties(State_VGB(:,i,j,k,iBlock), &
            EinternalSiIn = EinternalSi, PressureSiOut = PressureSi)

       State_VGB(p_,i,j,k,iBlock) = PressureSi*Si2No_V(UnitP_)

       State_VGB(EintExtra_,i,j,k,iBlock) = &
            Einternal - inv_gm1*State_VGB(p_,i,j,k,iBlock)

    end do; end do; end do

    call calc_energy_cell(iBlock)

  end subroutine update_conservative_energy

  !============================================================================

  subroutine heat_conduction(nVar, Temp_VGB, Rhs_VCB)

    use ModMain, ONLY: nBlock, unusedBLK
    use ModSize, ONLY: nI, nJ, nK, gcn, nBlk

    integer, intent(in) :: nVar
    real, intent(inout) :: &
         Temp_VGB(nVar,1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK)
    real, intent(out)   :: Rhs_VCB(nVar,1:nI,1:nJ,1:nK,nBLK)

    integer :: iBlock
    !--------------------------------------------------------------------------

    ! DoOneLayer, DoFacesOnly, No UseMonoteRestrict
    call message_pass_cells8(.true., .true., .false., nVar, Temp_VGB)

    call set_gray_outer_bcs(nVar, Temp_VGB)

    do iBlock = 1, nBlock
       if(unusedBlk(iBlock))CYCLE

       call get_heat_conduction(iBlock, nVar, &
            Temp_VGB(:,:,:,:,iBlock), Rhs_VCB(:,:,:,:,iBlock))
    end do

  end subroutine heat_conduction

  !============================================================================

  subroutine get_heat_conduction(iBlock, nVar, Temp_VG, Rhs_VC)

    use ModAdvance,  ONLY: Flux_VX, Flux_VY, Flux_VZ
    use ModGeometry, ONLY: dx_BLK, dy_BLK, dz_BLK, fAx_BLK, fAy_BLK, fAz_BLK, &
         vInv_CB, y_Blk
    use ModMain,     ONLY: x_, y_, z_, Dt, IsCylindrical
    use ModNodes,    ONLY: NodeY_NB
    use ModSize,     ONLY: nI, nJ, nK, gcn

    integer, intent(in) :: iBlock, nVar
    real, intent(inout) :: &
         Temp_VG(nVar,1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn)
    real, intent(out)   :: Rhs_VC(nVar,nI,nJ,nK)

    real :: AreaxInvDx, AreayInvDy, AreazInvDz, vInv
    integer :: i, j, k, iT, iCond
    !--------------------------------------------------------------------------

    ! Calculate heat conduction/diffusion fluxes
    AreaxInvDx = fAx_Blk(iBlock)/Dx_Blk(iBlock)
    AreayInvDy = fAy_Blk(iBlock)/Dy_Blk(iBlock)
    AreazInvDz = fAz_Blk(iBlock)/Dz_Blk(iBlock)
    do iCond = 1, nCond
       iT = iCond_I(iCond)
       do k = 1, nK; do j = 1, nJ; do i = 1, nI+1
          Flux_VX(iT,i,j,k) = -AreaxInvDx &
               *0.5*( HeatConductionCoef_IGB(iCond,i-1,j,k,iBlock) &
               +      HeatConductionCoef_IGB(iCond,i,  j,k,iBlock) ) &
               *(Temp_VG(iT,i,j,k) - Temp_VG(iT,i-1,j,k))
       end do; end do; end do
       do k = 1, nK; do j = 1, nJ+1; do i = 1, nI
          Flux_VY(iT,i,j,k) = -AreayInvDy &
               *0.5*( HeatConductionCoef_IGB(iCond,i,j-1,k,iBlock) &
               +      HeatConductionCoef_IGB(iCond,i,j,  k,iBlock) ) &
               *(Temp_VG(iT,i,j,k) - Temp_VG(iT,i,j-1,k))
       end do; end do; end do
       if(IsCylindrical)then
          ! Multiply fluxes with radius = abs(Y) at the X and Y faces
          do k = 1, nK; do j = 1, nJ; do i = 1, nI+1
             Flux_VX(iT,i,j,k) = Flux_VX(iT,i,j,k) &
                  *abs(y_BLK(i,j,k,iBlock))
          end do; end do; end do
          do k = 1, nK; do j = 1, nJ+1; do i = 1, nI
             Flux_VY(iT,i,j,k) = Flux_VY(iT,i,j,k) &
                  *abs(NodeY_NB(i,j,k,iBlock))
          end do; end do; end do

          ! There are no fluxes in the azimuthal direction
          do k = 1, nK+1; do j = 1, nJ; do i = 1, nI
             Flux_VZ(iT,i,j,k) = 0.0
          end do; end do; end do
       else
          do k = 1, nK+1; do j = 1, nJ; do i = 1, nI
             Flux_VZ(iT,i,j,k) = -AreazInvDz &
                  *0.5*( HeatConductionCoef_IGB(iCond,i,j,k-1,iBlock) &
                  +      HeatConductionCoef_IGB(iCond,i,j,k,  iBlock) ) &
                  *(Temp_VG(iT,i,j,k) - Temp_VG(iT,i,j,k-1))
          end do; end do; end do
       end if
    end do

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(IsCylindrical)then
          ! Multiply volume with radius (=Y) at cell center
          ! -> divide inverse volume
          vInv = vInv_CB(i,j,k,iBlock)/abs(y_BLK(i,j,k,iBlock))
       else
          vInv = vInv_CB(i,j,k,iBlock)
       end if

       ! Heat conduction
       Rhs_VC(:,i,j,k) = SpecificHeat_VCB(:,i,j,k,iBlock) &
            *Temp_VG(:,i,j,k)/Dt
       do iCond = 1, nCond
          iT = iCond_I(iCond)
          Rhs_VC(iT,i,j,k) = Rhs_VC(iT,i,j,k) + vInv*( &
               + Flux_VX(iT,i+1,j,k) - Flux_VX(iT,i,j,k) &
               + Flux_VY(iT,i,j+1,k) - Flux_VY(iT,i,j,k) &
               + Flux_VZ(iT,i,j,k+1) - Flux_VZ(iT,i,j,k) )
       end do

       if(nVar > 1)then
          ! Energy exchange
          Rhs_VC(Te_,i,j,k) = Rhs_VC(Te_,i,j,k) &
               + sum(RelaxationCoef_VCB(2:nVar,i,j,k,iBlock)) &
               *Temp_VG(Te_,i,j,k) &
               - sum( RelaxationCoef_VCB(2:nVar,i,j,k,iBlock) &
               *      Temp_VG(2:nVar,i,j,k) )

          Rhs_VC(2:nVar,i,j,k) = Rhs_VC(2:nVar,i,j,k) &
               + RelaxationCoef_VCB(2:nVar,i,j,k,iBlock) &
               *Temp_VG(2:nVar,i,j,k) &
               - RelaxationCoef_VCB(2:nVar,i,j,k,iBlock) &
               *Temp_VG(Te_,i,j,k)
       end if

       ! multiply by control volume
       Rhs_VC(:,i,j,k) = Rhs_VC(:,i,j,k)/vInv
    end do; end do; end do

  end subroutine get_heat_conduction

  !==========================================================================

  subroutine cg(matvec, nVar, Residual_VCB, Solution_VGB, IsInit, &
       Tolerance, TypeStop, Iter, preconditioner)

    ! conjugated gradient for symmetric and positive definite matrix A

    use ModMain, ONLY: nBlock, unusedBLK
    use ModSize, ONLY: nI, nJ, nK, gcn, nBlk

    ! subroutine for matrix vector multiplication 
    interface
       subroutine matvec(nVar, VectorIn_VGB, VectorOut_VCB)
         use ModSize, ONLY: nI, nJ, nK, gcn, nBlk

         ! Calculate VectorOut = M dot VectorIn,
         ! where M is the symmetric positive definite matrix

         ! To allow the usage of message_pass_cells,
         ! the input vector should be declared with intent(inout)
         ! and with ghost cells:
         integer, intent(in) :: nVar
         real, intent(inout) :: &
              VectorIn_VGB(nVar,1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK)
         real, intent(out)   :: VectorOut_VCB(nVar,1:nI,1:nJ,1:nK,nBLK)
       end subroutine matvec
    end interface

    optional :: preconditioner
    interface
       subroutine preconditioner(nVar, VectorIn_VCB, VectorOut_VCB)
         use ModSize

         ! Calculate VectorOut = M^{-1} \cdot VectorInwhere the preconditioner
         ! matrix, M^{-1}, should be symmetric positive definite
         integer, intent(in) :: nVar
         real, intent(inout) :: VectorIn_VCB(nVar,1:nI,1:nJ,1:nK,nBLK)
         real, intent(out)   :: VectorOut_VCB(nVar,1:nI,1:nJ,1:nK,nBLK)
       end subroutine preconditioner
    end interface


    integer, intent(in) :: nVar
    real, intent(inout) :: &      ! right hand side vector
         Residual_VCB(nVar,1:nI,1:nJ,1:nK,nBLK)
    real, intent(inout) :: &      ! initial guess / solution vector
         Solution_VGB(nVar,1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK)
    logical, intent(in) :: IsInit ! true if Solution_VGB contains initial guess
    real, intent(in) :: Tolerance ! required / achieved residual

    character (len=3), intent(in) :: TypeStop
    !      Determine stopping criterion (||.|| denotes the 2-norm):
    !      typestop='rel'  -- relative stopping crit.: ||res|| <= Tol*||res0||
    !      typestop='abs'  -- absolute stopping crit.: ||res|| <= Tol
 
    integer, intent(out) :: Iter ! actual number of iterations

    real :: rDotR, pDotADotP, rDotR0

    real, dimension(:,:,:,:,:), allocatable :: &
         P_VGB, ADotP_VCB, MMinusOneDotR_VCB

    integer :: i, j, k, iBlock
    !--------------------------------------------------------------------------

    allocate( &
         P_VGB(nVar,1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK), &
         ADotP_VCB(nVar,1:nI,1:nJ,1:nK,nBLK) )

    P_VGB = 0.0
    ADotP_VCB = 0.0

    if(present(preconditioner))then
       allocate(MMinusOneDotR_VCB(nVar,1:nI,1:nJ,1:nK,nBLK))
       MMinusOneDotR_VCB = 0.0
    end if

    !-------------- compute initial right hand side vector --------------

    if(IsInit)then
       call matvec(nVar, Solution_VGB, ADotP_VCB)
       do iBlock = 1, nBlock
          if(unusedBlk(iBlock))CYCLE
          Residual_VCB(:,:,:,:,iBlock) = Residual_VCB(:,:,:,:,iBlock) &
               - ADotP_VCB(:,:,:,:,iBlock)
       end do
    else
       do iBlock = 1, nBlock
          if(unusedBlk(iBlock))CYCLE
          Solution_VGB(:,:,:,:,iBlock) = 0.0
       end do
    end if

    Iter = 0

    LOOP: do

       if(present(preconditioner))then
          call preconditioner(nVar, Residual_VCB, MMinusOneDotR_VCB)
          rDotR = dot_product_mpi(0, Residual_VCB, MMinusOneDotR_VCB)
       else
          rDotR = dot_product_mpi(0, Residual_VCB, Residual_VCB)
       end if

       if(Iter == 0)then
          rDotR0 = rDotR
          if(rDotR0 == 0.0) EXIT
       end if

       if(TypeStop == 'abs' .and. rDotR <= Tolerance       ) EXIT
       if(TypeStop == 'rel' .and. rDotR <= Tolerance*rDotR0) EXIT

       do iBlock = 1, nBlock
          if(unusedBlk(iBlock))CYCLE
          if(present(preconditioner))then
             P_VGB(:,1:nI,1:nJ,1:nK,iBlock) = P_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
                  + MMinusOneDotR_VCB(:,:,:,:,iBlock)/rDotR
          else
             P_VGB(:,1:nI,1:nJ,1:nK,iBlock) = P_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
                  + Residual_VCB(:,:,:,:,iBlock)/rDotR
          end if
       end do

       call matvec(nVar, P_VGB, aDotP_VCB)

       pDotADotP = dot_product_mpi(gcn, P_VGB(:,:,:,:,:), aDotP_VCB)
       do iBlock = 1, nBlock
          if(unusedBlk(iBlock))CYCLE
          Residual_VCB(:,:,:,:,iBlock) = Residual_VCB(:,:,:,:,iBlock) &
               - aDotP_VCB(:,:,:,:,iBlock)/pDotADotP
          Solution_VGB(:,1:nI,1:nJ,1:nK,iBlock) = &
               Solution_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
               + P_VGB(:,1:nI,1:nJ,1:nK,iBlock)/pDotADotP
       end do

       Iter = Iter + 1
    end do LOOP

    deallocate(P_VGB, aDotP_VCB)
    if(allocated(MMinusOneDotR_VCB)) deallocate(MMinusOneDotR_VCB)

  contains

    real function dot_product_mpi(gcn, A_VGB, B_VCB)

      use ModMain,     ONLY: nBlock, unusedBLK
      use ModGeometry, ONLY: true_cell, body_blk
      use ModMpi
      use ModProcMH,   ONLY: iComm

      integer, intent(in) :: gcn
      real, intent(in) :: &
           A_VGB(nVar,1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBlk)
      real, intent(in) :: B_VCB(nVar,1:nI,1:nJ,1:nK,nBlk)
 
      real :: DotProduct, DotProductMpi
      integer :: iError, i, j, k, iBlock
      !------------------------------------------------------------------------

      DotProduct = 0.0

      do iBlock = 1, nBlock
         if(unusedBlk(iBlock))CYCLE
         if(body_blk(iBlock))then
            do k=1,nK; do j=1,nJ; do i=1,nI
               if(.not.true_cell(i,j,k,iBlock))CYCLE
               DotProduct = DotProduct &
                    + sum( A_VGB(:,i,j,k,iBlock)*B_VCB(:,i,j,k,iBlock) )
            end do; end do; end do
         else
            DotProduct = DotProduct &
                 + sum( A_VGB(:,1:nI,1:nJ,1:nK,iBlock)*B_VCB(:,:,:,:,iBlock) )
         end if
      end do

      call MPI_allreduce(DotProduct, DotProductMpi, 1, MPI_REAL, MPI_SUM, &
           iComm, iError)

      dot_product_mpi = DotProductMpi

    end function dot_product_mpi

  end subroutine cg

  !============================================================================

  subroutine jacobi_preconditioner(nVar, VectorIn_VCB, VectorOut_VCB)

    ! Calculate VectorOut = M^{-1} \cdot VectorIn where the stored
    ! preconditioner matrix, M^{-1}, should be symmetric positive definite

    use ModGeometry, ONLY: true_cell, body_blk
    use ModMain, ONLY: nBlock, unusedBLK
    use ModSize, ONLY: nI, nJ, nK, nBlk

    integer, intent(in) :: nVar
    real, intent(inout) :: VectorIn_VCB(nVar,1:nI,1:nJ,1:nK,nBLK)
    real, intent(out)   :: VectorOut_VCB(nVar,1:nI,1:nJ,1:nK,nBLK)

    integer :: iVar, i, j , k, iBlock

    character(len=*), parameter :: NameSub = 'jacobi_preconditioner'
    !--------------------------------------------------------------------------

    do iBlock = 1, nBlock
       if(unusedBlk(iBlock))CYCLE

       if(body_blk(iBlock))then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock))CYCLE
             do iVar = 1, nVar
                VectorOut_VCB(iVar,i,j,k,iBlock) = sum( &
                     VectorIn_VCB(:,i,j,k,iBlock)*JacPrec_VVCB(:,iVar,i,j,k,iBlock))
             end do
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             do iVar = 1, nVar
                VectorOut_VCB(iVar,i,j,k,iBlock) = sum( &
                     VectorIn_VCB(:,i,j,k,iBlock)*JacPrec_VVCB(:,iVar,i,j,k,iBlock))
             end do
          end do; end do; end do
       end if
    end do

  end subroutine jacobi_preconditioner

  !============================================================================

  subroutine get_jacobi_preconditioner(nVar)

    use ModGeometry, ONLY: true_cell, body_blk, Dx_Blk, Dy_Blk, Dz_Blk, &
         fAx_BLK, fAy_BLK, fAz_BLK, vInv_CB
    use ModMain, ONLY: nBlock, unusedBLK
    use ModSize, ONLY: nI, nJ, nK

    ! Inverts the nTemperature * nTemperature sub-matrices at the
    ! block diagonal of the operator of heat-conduction + relaxation

    integer, intent(in) :: nVar

    integer :: iVar, i, j , k, iBlock

    real, dimension(nVar) :: Diag_V
    real, allocatable, dimension(:) :: OffDiag_V

    real :: AreaxInvDx, AreayInvDy, AreazInvDz

    character(len=*), parameter :: NameSub = 'get_jacobi_preconditioner'
    !--------------------------------------------------------------------------

    if(nVar>1) allocate(OffDiag_V(2:nVar))

    do iBlock = 1, nBlock
       if(unusedBlk(iBlock))CYCLE

       AreaxInvDx = fAx_Blk(iBlock)/Dx_Blk(iBlock)
       AreayInvDy = fAy_Blk(iBlock)/Dy_Blk(iBlock)
       AreazInvDz = fAz_Blk(iBlock)/Dz_Blk(iBlock)

       select case(nVar)
       case(1)
          if(body_blk(iBlock))then
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                if(.not.true_cell(i,j,k,iBlock))CYCLE
                call block_jacobi_1t
             end do; end do; end do
          else
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                call block_jacobi_1t
             end do; end do; end do
          end if
       case(2)
          if(body_blk(iBlock))then
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                if(.not.true_cell(i,j,k,iBlock))CYCLE
                call block_jacobi_2t
             end do; end do; end do
          else
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                call block_jacobi_2t
             end do; end do; end do
          end if
       case(3)
          if(body_blk(iBlock))then
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                if(.not.true_cell(i,j,k,iBlock))CYCLE
                call block_jacobi_3t
             end do; end do; end do
          else
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                call block_jacobi_3t
             end do; end do; end do
          end if
       case default
          call stop_mpi(NameSub// &
               ': More than three temperatures is not allowed')
       end select
    end do

    if(allocated(OffDiag_V)) deallocate(OffDiag_V)

  contains

    subroutine block_jacobi_1t

      !------------------------------------------------------------------------

      call get_diagonal_subblock

      JacPrec_VVCB(1,1,i,j,k,iBlock) = 1.0/Diag_V(1)

    end subroutine block_jacobi_1t

    !==========================================================================

    subroutine block_jacobi_2t

      !------------------------------------------------------------------------

      call get_diagonal_subblock

      JacPrec_VVCB(:,:,i,j,k,iBlock) = reshape( (/ &
           Diag_V(2)    , -OffDiag_V(2), &
           -OffDiag_V(2), Diag_V(1)  /), (/2,2/) ) &
           /(Diag_V(1)*Diag_V(2) - OffDiag_V(2)**2)

    end subroutine block_jacobi_2t

    !==========================================================================

    subroutine block_jacobi_3t

      !------------------------------------------------------------------------

      call get_diagonal_subblock

      JacPrec_VVCB(:,:,i,j,k,iBlock) = reshape( (/ &
           Diag_V(2)*Diag_V(3), -Diag_V(3)*OffDiag_V(2), -Diag_V(2)*OffDiag_V(3), &
           -Diag_V(3)*OffDiag_V(2), Diag_V(1)*Diag_V(3) - OffDiag_V(3)**2, OffDiag_V(2)*OffDiag_V(3), &
           -Diag_V(2)*OffDiag_V(3), OffDiag_V(2)*OffDiag_V(3), Diag_V(1)*Diag_V(2) - OffDiag_V(2)**2 /), (/3,3/) ) &
           /( Diag_V(3)*(Diag_V(1)*Diag_V(2) - OffDiag_V(2)**2) - Diag_V(2)*OffDiag_V(3)**2 )

    end subroutine block_jacobi_3t

    !==========================================================================

    subroutine get_diagonal_subblock

      use ModGeometry, ONLY: y_Blk
      use ModMain,     ONLY: Dt, IsCylindrical
      use ModNodes,    ONLY: NodeY_NB

      integer :: iCond, iT
      real :: Volume
      !------------------------------------------------------------------------

      if(IsCylindrical)then

         ! Multiply volume with radius (=Y) at cell center
         Volume = abs(y_Blk(i,j,k,iBlock))/vInv_CB(i,j,k,iBlock)

         ! Heat conduction
         Diag_V(:) = SpecificHeat_VCB(:,i,j,k,iBlock)*Volume/Dt
         do iCond = 1, nCond
            iT = iCond_I(iCond)
            Diag_V(iT) = Diag_V(iT) &
                 + AreaxInvDx*abs(y_Blk(i,j,k,iBlock))*( &
                 +     0.5*HeatConductionCoef_IGB(iCond,i-1,j,  k,iBlock) &
                 +         HeatConductionCoef_IGB(iCond,i,  j,  k,iBlock) &
                 +     0.5*HeatConductionCoef_IGB(iCond,i+1,j,  k,iBlock) ) &
                 + AreayInvDy*0.5*( &
                 +     abs(NodeY_NB(i,j,k,iBlock))*( &
                 +         HeatConductionCoef_IGB(iCond,i,  j-1,k,iBlock) &
                 +         HeatConductionCoef_IGB(iCond,i,  j,  k,iBlock)) &
                 +     abs(NodeY_NB(i,j+1,k,iBlock))*( &
                 +         HeatConductionCoef_IGB(iCond,i,  j  ,k,iBlock) &
                 +         HeatConductionCoef_IGB(iCond,i,  j+1,k,iBlock)) )
         end do

      else

         Volume = 1.0/vInv_CB(i,j,k,iBlock)

         ! Heat conduction
         Diag_V(:) = SpecificHeat_VCB(:,i,j,k,iBlock)*Volume/Dt
         do iCond = 1, nCond
            iT = iCond_I(iCond)
            Diag_V(iT) = Diag_V(iT) &
              + AreaxInvDx*( &
              +     0.5*HeatConductionCoef_IGB(iCond,i-1,j,  k,  iBlock) &
              +         HeatConductionCoef_IGB(iCond,i,  j,  k,  iBlock) &
              +     0.5*HeatConductionCoef_IGB(iCond,i+1,j,  k,  iBlock) )&
              + AreayInvDy*( &
              +     0.5*HeatConductionCoef_IGB(iCond,i,  j-1,k,  iBlock) &
              +         HeatConductionCoef_IGB(iCond,i,  j,  k,  iBlock) &
              +     0.5*HeatConductionCoef_IGB(iCond,i,  j+1,k,  iBlock) )&
              + AreazInvDz*( &
              +     0.5*HeatConductionCoef_IGB(iCond,i,  j,  k-1,iBlock) &
              +         HeatConductionCoef_IGB(iCond,i,  j,  k,  iBlock) &
              +     0.5*HeatConductionCoef_IGB(iCond,i,  j,  k+1,iBlock) )
         end do

      end if

      if(nVar==1) return

      ! Energy exchange
      Diag_V(Te_) = Diag_V(Te_) &
           + Volume*sum(RelaxationCoef_VCB(2:nVar,i,j,k,iBlock))
      Diag_V(2:nVar) = Diag_V(2:nVar) &
           + Volume*RelaxationCoef_VCB(2:nVar,i,j,k,iBlock)
      OffDiag_V(2:nVar) = -Volume*RelaxationCoef_VCB(2:nVar,i,j,k,iBlock)

    end subroutine get_diagonal_subblock

  end subroutine get_jacobi_preconditioner

  !===========================================================================

  subroutine set_gray_outer_bcs(nVar, Var_VGB)

    use ModMain, ONLY: nI, nJ, nK, nBlk, nBlock, UnusedBlk
    use ModParallel, ONLY: NOBLK, NeiLev

    integer, intent(in):: nVar
    real, intent(inout):: Var_VGB(nVar,-1:nI+2,-1:nJ+2,-1:nK+2,nBlk)
    
    integer :: iBlock
    !----------------------------------------------------------------------

    do iBlock = 1, nBlock
       if(UnusedBlk(iBlock))CYCLE

       ! set floating outer boundary !!! shear to be added !!!
       if(NeiLev(1,iBlock) == NOBLK) Var_VGB(:,0   ,:,:,iBlock) &
            =                        Var_VGB(:,1   ,:,:,iBlock)
       if(NeiLev(2,iBlock) == NOBLK) Var_VGB(:,nI+1,:,:,iBlock) &
            =                        Var_VGB(:,nI  ,:,:,iBlock)
       if(NeiLev(3,iBlock) == NOBLK) Var_VGB(:,:,0   ,:,iBlock) &
            =                        Var_VGB(:,:,1   ,:,iBlock)
       if(NeiLev(4,iBlock) == NOBLK) Var_VGB(:,:,nJ+1,:,iBlock) &
            =                        Var_VGB(:,:,nJ  ,:,iBlock)
       if(NeiLev(5,iBlock) == NOBLK) Var_VGB(:,:,:,0   ,iBlock) &
            =                        Var_VGB(:,:,:,1   ,iBlock)
       if(NeiLev(6,iBlock) == NOBLK) Var_VGB(:,:,:,nK+1,iBlock) &
            =                        Var_VGB(:,:,:,nK  ,iBlock)
    end do

  end subroutine set_gray_outer_bcs

end module ModGrayDiffusion
