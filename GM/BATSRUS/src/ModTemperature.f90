!^CFG COPYRIGHT UM
!==============================================================================
module ModTemperature

  use ModVarIndexes, ONLY: p_

  implicit none
  save

  ! Solves the non-linear heat conduction problem in the coupled system
  ! of the energy equations. The internal energies involved can be any
  ! combination of electrons, ions, and/or radiation.
  !!! Currently, the only boundary conditions that are allowed for the
  !!! temperatures are floating.

  private !except

  ! Public methods
  public :: read_temperature_param
  public :: init_temperature_diffusion
  public :: calc_source_temperature_diff
  public :: advance_temperature

  Logical, public :: UseTemperatureDiffusion = .false.

  ! Parameters for radiation flux limiter
  logical,           public :: UseRadFluxLimiter  = .false.
  character(len=20), public :: TypeRadFluxLimiter = 'larsen'

  ! Coefficients for the temperature diffusion model
  real, allocatable :: RelaxationCoef_VCB(:,:,:,:,:)
  real, allocatable :: SpecificHeat_VCB(:,:,:,:,:)
  real, allocatable :: HeatConductionCoef_IGB(:,:,:,:,:)

  ! Face centered values for the heat conduction coefficients
  ! multiplied by face area divided by cell-to-cell
  ! center distance across the given face  
  real, allocatable, dimension(:,:,:,:) ::&
       Cond_VFX, Cond_VFY, Cond_VFZ

  ! For the purpose of message passing the heat conduction coefficients
  ! are compactly stored without gaps. A pointer is used to indicate
  ! which temperature variables involve heat conduction
  integer, allocatable :: iCond_I(:)

  ! number of variables that uses "heat conduction"
  integer :: nCond

  real, parameter :: GammaRel = 4.0/3.0

  ! local index for extra internal energy to keep the compiler happy
  integer, parameter :: EintExtra_ = p_ - 1

  ! the solution vector containing the independent temperature variables
  real, allocatable :: Temperature_VGB(:,:,:,:,:), &
       TemperatureOld_VCB(:,:,:,:,:)

  ! number of independent temperature variables
  integer :: nTemperature

  ! Named indexes for the temerature variables:
  ! The electron temperature is always the first temperature. If it
  ! is the only temperature in use, then it is effectively the material
  ! temperature.
  integer, parameter :: Te_ = 1
  integer, parameter :: aTe4_ = Te_
  integer :: Tion_, Trad_
  integer :: aTrad4_

  ! logicals if the ion and radiation temperature are in use
  ! These switches are also used to determine which energy exchange terms
  ! are needed: UseTion indicates exchange between electron and ions,
  ! UseTrad indicates exchange between electron and radiation.
  logical :: UseTion, UseTrad

  ! Indices to figure out where the energies variables are in State_VGB
  ! The energies are the electron and radiation internal energy
  ! The material energy is known, since this is Energy_
  integer :: iEe, iErad

  ! logical to select if temperature T (.true.) will be used or T^4 (.false.)
  logical :: UseTemperatureVariable = .true.

  ! variables for the accuracy of the conjugate gradient
  character(len=3) :: TypeStopCriterion = 'rel'
  real :: MaxErrorResidual = 1e-6

  ! Right-hand side of non-linear coupled system of temperature equations
  real, allocatable :: Source_VCB(:,:,:,:,:)

  character(len=20) :: TypePreconditioner = 'blockjacobi'

  ! Heptadiagonal mbilu preconditioner
  integer, parameter:: nStencil = 7
  real, allocatable :: Prec_VVCIB(:,:,:,:,:,:,:)

  ! Store block Jacobi preconditioner
  real, allocatable :: JacPrec_VVCB(:,:,:,:,:,:)

  !To switch the boundary condition at the outer boundary:
  logical :: IsFirstCgIteration

  ! To make sure that the dominant contribution to the norm
  ! TRad^2/Cv(TRad)\sim (1/T) is not due to infinitesimal temperatures
  real :: TradMin, EradMin

  real :: rDotRPe, pDotADotPPe

contains

  !============================================================================

  subroutine read_temperature_param(NameCommand)

    use ModMain, ONLY: UseGrayDiffusion, UseHeatConduction
    use ModPhysics,   ONLY: cRadiationNo, Si2No_V, UnitTemperature_
    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand

    real :: TradMinSi

    character(len=*), parameter :: NameSub = 'read_temperature_param'
    !--------------------------------------------------------------------------

    select case(NameCommand)
    case("#IMPLICITTEMPERATURE")
       call read_var('UseTemperatureDiffusion', UseTemperatureDiffusion)
       if(UseTemperatureDiffusion)then
          call read_var('UseTemperatureVariable', UseTemperatureVariable)
          call read_var('TypePreconditioner', TypePreconditioner, &
               IsLowerCase=.true.)
          select case(TypePreconditioner)
          case('blockjacobi','gs','dilu','mbilu')
          case default
             call stop_mpi(NameSub//': unknown TypePreconditioner='&
                  //TypePreconditioner)
          end select
          call read_var('TypeStopCriterion', TypeStopCriterion, &
               IsLowerCase=.true.)
          call read_var('MaxErrorResidual', MaxErrorResidual)
       end if
    case("#HEATCONDUCTION")
       call read_var('UseHeatConduction', UseHeatConduction)

    case("#RADIATION")
       call read_var('UseGrayDiffusion', UseGrayDiffusion)
       if(UseGrayDiffusion)then
          call read_var('UseRadFluxLimiter', UseRadFluxLimiter)
          if(UseRadFluxLimiter)then
             call read_var('TypeRadFluxLimiter', TypeRadFluxLimiter, &
                  IsLowerCase=.true.)

             select case(TypeRadFluxLimiter)
             case("larsen","sum","max")
             case default
                call stop_mpi(NameSub//': unknown TypeRadFluxLimiter='&
                     //TypeRadFluxLimiter)
             end select
          end if
          call read_var('TradMinSi', TradMinSi)
          TradMin = TradMinSi*Si2No_V(UnitTemperature_)
          EradMin = cRadiationNo*TradMin**4
       end if
    case default
       call stop_mpi(NameSub//' invalid NameCommand='//NameCommand)
    end select

  end subroutine read_temperature_param

  !============================================================================

  subroutine init_temperature_diffusion

    use ModMain,       ONLY: UseHeatConduction, UseGrayDiffusion
    use ModProcMH,     ONLY: iProc
    use ModSize,       ONLY: nI, nJ, nK, nBlk, nDim
    use ModVarIndexes, ONLY: NameVar_V, nVar
    integer :: iVar

    character(len=*), parameter :: NameSub = "init_temperature_diffusion"
    !--------------------------------------------------------------------------

    ! First default the additional temperatures and energies to not in use
    UseTrad = .false.
    UseTion = .false.
    iErad = -1
    iEe = -1

    ! Based on the energy variable, check if the radiation and ion
    ! temperatures are needed.
    do iVar = 1, nVar
       if(trim(NameVar_V(iVar)) == "Erad")then
          UseTrad = .true.
          iErad = iVar
       end if
       if(trim(NameVar_V(iVar)) == "Ee")then
          !!! Ee should be changed if a different NameVar will be used.
          ! The electron and ion temperature are not combined
          ! (material temperature), so we need besides Te also Tion.
          UseTion = .true.
          iEe = iVar
       end if
    end do

    ! The electron temperature is by default switched on.
    ! Do set the additional temperature temperature variables.
    nTemperature = 1
    if(UseTion)then
       nTemperature = nTemperature + 1
       Tion_ = nTemperature
    end if
    if(UseTrad)then
       nTemperature = nTemperature + 1
       Trad_ = nTemperature
       aTrad4_ = Trad_
    end if

    ! count the number of temperature variables that involve "heat conduction"
    nCond = 0
    if(UseHeatConduction) nCond = nCond + 1
    if(UseGrayDiffusion)  nCond = nCond + 1

    if(.not.allocated(SpecificHeat_VCB))then
       allocate( &
            SpecificHeat_VCB(nTemperature,nI,nJ,nK,nBlk), &
            TemperatureOld_VCB(nTemperature,nI,nJ,nK,nBlk), &
            Temperature_VGB(nTemperature,-1:nI+2,-1:nJ+2,-1:nK+2,nBlk), &
            Source_VCB(nTemperature,nI,nJ,nK,nBlk) )

       if(nCond > 0) then
            allocate( &
            HeatConductionCoef_IGB(nCond,-1:nI+2,-1:nJ+2,-1:nK+2,nBlk), &
            iCond_I(nCond) )
            
            allocate(Cond_VFX(nCond, nI+1, nJ, nK), &
                     Cond_VFY(nCond, nI, nJ+1, nK), &
                     Cond_VFZ(nCond, nI, nJ, nK+1))
         end if

       if(nTemperature > 1) &
            allocate(RelaxationCoef_VCB(2:nTemperature,nI,nJ,nK,nBlk))

       select case(TypePreconditioner)
       case('blockjacobi','gs','dilu')
          allocate(JacPrec_VVCB(nTemperature,nTemperature,nI,nJ,nK,nBLK))
       case('mbilu')
          allocate(Prec_VVCIB(nTemperature,nTemperature,nI,nJ,nK,nStencil,nBlk))
       end select

    end if

    ! Select which temperature variables involve "heat conduction"
    if(UseHeatConduction) iCond_I(1) = Te_
    if(UseGrayDiffusion)  iCond_I(nCond) = Trad_

  end subroutine init_temperature_diffusion

  !============================================================================

  subroutine set_temperature

    use ModAdvance, ONLY: State_VGB,p_
    use ModMain,    ONLY: nI, nJ, nK, nBlock, unusedBlk
    use ModPhysics, ONLY: cRadiationNo, Si2No_V, UnitTemperature_
    use ModUser,    ONLY: user_material_properties

    integer :: i, j, k, iBlock
    real :: Te, TeSi, Trad

    character(len=*), parameter:: NameSub = 'set_temperature'
    !--------------------------------------------------------------------------

    ! set the temperature variables
    !!! ion temperature can not yet be set
    do iBlock = 1, nBlock
       if(unusedBlk(iBlock)) CYCLE

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(State_VGB(p_,i,j,k,iBlock) <= 0.0)&
               call stop_mpi('Negative pressure in set_temperature')
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               TeSiOut = TeSi)

          Te = TeSi*Si2No_V(UnitTemperature_)
          if(Te <= 0.0)call stop_mpi('negative temperature in set_temperature')
          if(UseTemperatureVariable)then
             Temperature_VGB(Te_,i,j,k,iBlock) = Te
             if(UseTrad)then
                if(State_VGB(iERad,i,j,k,iBlock)<0.0) &
                     call stop_mpi('negative radiation energy')
                Trad = sqrt(sqrt(State_VGB(iErad,i,j,k,iBlock) &
                     /cRadiationNo))
                Temperature_VGB(Trad_,i,j,k,iBlock) = Trad
             end if
          else
             Temperature_VGB(aTe4_,i,j,k,iBlock) = cRadiationNo*Te**4
             if(UseTrad)then
                Temperature_VGB(aTrad4_,i,j,k,iBlock) = &
                     State_VGB(iErad,i,j,k,iBlock)
             end if
          end if

          TemperatureOld_VCB(:,i,j,k,iBlock) = Temperature_VGB(:,i,j,k,iBlock)
       end do; end do; end do
    end do

  end subroutine set_temperature

  !============================================================================

  subroutine set_frozen_coefficients

    use ModAdvance,  ONLY: State_VGB
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK
    use ModMain,     ONLY: nI, nJ, nK, nBlock, unusedBlk, &
         UseHeatConduction, UseGrayDiffusion
    use ModPhysics,  ONLY: Si2No_V, UnitX_, UnitEnergyDens_, &
         UnitTemperature_, cRadiationNo, Clight, UnitU_
    use ModProcMH,   ONLY: iProc
    use ModUser,     ONLY: user_material_properties

    integer :: i, j, k, iBlock
    real :: PlanckOpacitySi, RosselandMeanOpacitySi, HeatConductionCoefSi
    real :: PlanckOpacity, RosselandMeanOpacity, HeatConductionCoef
    real :: CvSi, Cv, TeSi, Te, Trad, DiffRad

    character(len=*), parameter:: NameSub = 'set_frozen_coefficients'
    !--------------------------------------------------------------------------

    if(UseRadFluxLimiter)then
       ! The radiation flux limiters use the gradient of the
       ! radiation temperature and therefor one layer of ghost cells
       ! needs to be correct. We massage pass all temperatures.

       ! DoOneLayer, DoFacesOnly, No UseMonoteRestrict
       call message_pass_cells8(.true., .true., .false., nTemperature, &
            Temperature_VGB)
       call set_gray_outer_bcs(nTemperature, Temperature_VGB,'float')
    end if

    ! fix the specific heat, heat conduction/diffusion,
    ! and energy exchange coeffients.
    !!! This part of the code is at the moment for
    !!! non-equilibrium gray-diffusion radiation only.
    do iBlock = 1, nBlock
       if(unusedBlk(iBlock)) CYCLE

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               AbsorptionOpacitySiOut = PlanckOpacitySi, &
               RosselandMeanOpacitySiOut = RosselandMeanOpacitySi, &
               CvSiOut = CvSi, TeSiOut = TeSi, &
               HeatConductionCoefSiOut = HeatConductionCoefSi)

          RosselandMeanOpacity = RosselandMeanOpacitySi/Si2No_V(UnitX_)
          PlanckOpacity = PlanckOpacitySi/Si2No_V(UnitX_)
          HeatConductionCoef = HeatConductionCoefSi &
               *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_) &
               *Si2No_V(UnitU_)*Si2No_V(UnitX_)
          Cv = CvSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
          Te = TeSi*Si2No_V(UnitTemperature_)

          if(UseTemperatureVariable)then
             SpecificHeat_VCB(Te_,i,j,k,iBlock) = Cv

             if(UseGrayDiffusion)then
                Trad = max(Temperature_VGB(Trad_,i,j,k,iBlock),TradMin)
                SpecificHeat_VCB(Trad_,i,j,k,iBlock) = 4.0*cRadiationNo*Trad**3
                RelaxationCoef_VCB(Trad_,i,j,k,iBlock) = Clight*PlanckOpacity &
                     *cRadiationNo*(Te+Trad)*(Te**2+Trad**2)

                call get_radiation_diffusion_coef
                HeatConductionCoef_IGB(nCond,i,j,k,iBlock) = &
                     DiffRad*4.0*cRadiationNo*Trad**3
             end if

             if(UseHeatConduction) &
                HeatConductionCoef_IGB(1,i,j,k,iBlock) = HeatConductionCoef

          else
             SpecificHeat_VCB(aTe4_,i,j,k,iBlock) = Cv/(4.0*cRadiationNo*Te**3)
             SpecificHeat_VCB(aTrad4_,i,j,k,iBlock) = 1.0
             RelaxationCoef_VCB(aTrad4_,i,j,k,iBlock) = Clight*PlanckOpacity

             call get_radiation_diffusion_coef
             HeatConductionCoef_IGB(1,i,j,k,iBlock) = DiffRad
          end if

       end do; end do; end do
    end do

    ! message pass one ghost cell layer of heat conduction coefficients.

    ! DoOneLayer, DoFacesOnly, No UseMonoteRestrict
    call message_pass_cells8(.true., .true., .false., nCond, &
         HeatConductionCoef_IGB)

    ! The ghost cell filling of the heat conduction coefficients
    ! is not really needed since the temperature boundary conditions are float
    call set_gray_outer_bcs(nCond, HeatConductionCoef_IGB,'float')

  contains

    subroutine get_radiation_diffusion_coef

      real :: Grad2ByErad2, Erad, Grad_D(3)
      !------------------------------------------------------------------------

      if(UseRadFluxLimiter)then
         if(UseTemperatureVariable)then
            call calc_cell_gradient(Trad_, Grad_D)
            Grad2ByErad2 = 16.0*sum(Grad_D**2)/Trad**2
         else
            call calc_cell_gradient(aTrad4_, Grad_D)
            Erad = Temperature_VGB(aTrad4_,i,j,k,iBlock)
            Grad2ByErad2 = sum(Grad_D**2)/Erad**2
         end if

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

    end subroutine get_radiation_diffusion_coef

    !==========================================================================

    subroutine calc_cell_gradient(iVar, Grad_D)

      ! currently, only a cartesian version

      use ModGeometry, ONLY: Dx_Blk, Dy_Blk, Dz_Blk
      use ModMain,     ONLY: x_, y_, z_

      integer, intent(in) :: iVar
      real, intent(out) :: Grad_D(3)

      real :: InvDx2, InvDy2, InvDz2
      !------------------------------------------------------------------------

      InvDx2 = 0.5/Dx_Blk(iBlock)
      InvDy2 = 0.5/Dy_Blk(iBlock)
      InvDz2 = 0.5/Dz_Blk(iBlock)

      Grad_D(x_) = &
           ( Temperature_VGB(iVar,i+1,j,  k,  iBlock) &
           - Temperature_VGB(iVar,i-1,j,  k,  iBlock) )*InvDx2
      Grad_D(y_) = &
           ( Temperature_VGB(iVar,i,  j+1,k,  iBlock) &
           - Temperature_VGB(iVar,i,  j-1,k,  iBlock) )*InvDy2
      Grad_D(z_) = &
           ( Temperature_VGB(iVar,i,  j,  k+1,iBlock) &
           - Temperature_VGB(iVar,i,  j,  k-1,iBlock) )*InvDz2

    end subroutine calc_cell_gradient

  end subroutine set_frozen_coefficients

  !============================================================================

  subroutine calc_source_temperature_diff(iBlock)

    use ModAdvance,    ONLY: State_VGB, Source_VC, &
         uDotArea_XI, uDotArea_YI, uDotArea_ZI
    use ModGeometry,   ONLY: vInv_CB, TypeGeometry, y_BLK
    use ModMain,       ONLY: nI, nJ, nK, UseGrayDiffusion
    use ModVarIndexes, ONLY: Energy_, RhoUy_

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: vInv, DivU, RadCompression
    character(len=*), parameter :: NameSub = "calc_source_temperature_diff"
    !--------------------------------------------------------------------------
    
    do k=1,nK; do j=1,nJ; do i=1,nI

       vInv = vInv_CB(i,j,k,iBlock)
       DivU = uDotArea_XI(i+1,j,k,1) - uDotArea_XI(i,j,k,1) &
            + uDotArea_YI(i,j+1,k,1) - uDotArea_YI(i,j,k,1) &
            + uDotArea_ZI(i,j,k+1,1) - uDotArea_ZI(i,j,k,1) 
 
       ! Adiabatic compression of radiation by fluid velocity
       if(UseTrad)then
          ! (GammaRel-1)*Erad*Div(U)
          RadCompression = &
               (GammaRel-1.0)*State_VGB(iErad,i,j,k,iBlock)*vInv*DivU

          ! dErad/dt = - adiabatic compression
          Source_VC(iErad,i,j,k) = Source_VC(iErad,i,j,k) - RadCompression

          ! dE/dt    = + adiabatic compression
          Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + RadCompression
       end if

    end do; end do; end do

    if(TypeGeometry=='rz' .and. UseGrayDiffusion)then
       ! Add "geometrical source term" p/r to the radial momentum equation
       ! The "radial" direction is along the Y axis
       ! NOTE: here we have to use signed radial distance!
       do k=1,nK; do j=1, nJ; do i=1, nI
          Source_VC(RhoUy_,i,j,k) = Source_VC(RhoUy_,i,j,k) &
               + (1./3.)*State_VGB(iErad,i,j,k,iBlock) &
               / y_BLK(i,j,k,iBlock)
       end do; end do; end do
    end if

  end subroutine calc_source_temperature_diff

  !============================================================================

  subroutine advance_temperature

    use ModAdvance,  ONLY: State_VGB
    use ModGeometry, ONLY: vInv_CB
    use ModMain,     ONLY: nI, nJ, nK, nBlock, unusedBLK, Dt

    integer :: i, j, k, iBlock, Iter
    real :: Error
    logical :: DoTestKrylov, DoTestKrylovMe

    character(len=*), parameter :: NameSub = 'advance_temperature'
    !--------------------------------------------------------------------------

    call set_oktest('krylov', DoTestKrylov, DoTestKrylovMe)

    ! prevent doing useless updates
    if(Dt == 0.0)return

    call set_temperature

    call set_frozen_coefficients

    do iBlock = 1, nBlock
       if(unusedBlk(iBlock)) CYCLE

       do k = 1,nK; do j = 1,nJ; do i = 1,nI
          Source_VCB(:,i,j,k,iBlock) = SpecificHeat_VCB(:,i,j,k,iBlock) &
               *Temperature_VGB(:,i,j,k,iBlock)/(vInv_CB(i,j,k,iBlock)*Dt)
       end do; end do; end do

    end do

    ! do a preconditioned conjugate gradient

    !To switch the boundary condition at zero iteration
    IsFirstCgIteration = .true. 

    Error = MaxErrorResidual

    select case(TypePreconditioner)
    case('blockjacobi')
       call get_jacobi_preconditioner(nTemperature)

       call cg(heat_conduction, nTemperature, Source_VCB, Temperature_VGB, &
            .true., Error, TypeStopCriterion, &
            Iter, DoTestKrylovMe, jacobi_preconditioner)
    case('gs')
       call get_jacobi_preconditioner(nTemperature)

       call cg(heat_conduction, nTemperature, Source_VCB, Temperature_VGB, &
            .true., Error, TypeStopCriterion, &
            Iter, DoTestKrylovMe, dilu_preconditioner)
    case('dilu')
       call get_dilu_preconditioner(nTemperature)

       call cg(heat_conduction, nTemperature, Source_VCB, Temperature_VGB, &
            .true., Error, TypeStopCriterion, &
            Iter, DoTestKrylovMe, dilu_preconditioner)
    case('mbilu')
       call get_mbilu_preconditioner(nTemperature)

       call cg(heat_conduction, nTemperature, Source_VCB, Temperature_VGB, &
            .true., Error, TypeStopCriterion, &
            Iter, DoTestKrylovMe, mbilu_preconditioner)
    end select

    if(DoTestKrylovMe)write(*,*)NameSub,': Number of iterations, Error =', &
         Iter, Error

    do iBlock = 1, nBlock
       if(unusedBlk(iBlock)) CYCLE
       if(nTemperature==2 .and. UseTemperatureVariable)call check_temperature
       call update_conservative_energy(iBlock)
    end do

  contains

    subroutine check_temperature

      real :: AveragedTemperature, SpecificHeatTotal, DiffT
      real, parameter :: TemperatureJumpAllowed = 0.9
      !------------------------------------------------------------------------

      do k = 1,nK; do j = 1,nJ; do i = 1,nI
         SpecificHeatTotal = sum(SpecificHeat_VCB(:,i,j,k,iBlock))
         AveragedTemperature = sum(Temperature_VGB(:,i,j,k,iBlock)*&
                                   SpecificHeat_VCB(:,i,j,k,iBlock))&
                                   /SpecificHeatTotal
                      
         if(AveragedTemperature<=0.0) &
              call stop_mpi('NEGATIVE AVERAGED TEMPERATURE')

         DiffT = Temperature_VGB(Te_,i,j,k,iBlock) &
              - Temperature_VGB(TRad_,i,j,k,iBlock)
         DiffT = min(max(DiffT,-AveragedTemperature*SpecificHeatTotal/&
              SpecificHeat_VCB(TRad_,i,j,k,iBlock)*TemperatureJumpAllowed),&
              AveragedTemperature*SpecificHeatTotal/&
              SpecificHeat_VCB(Te_,i,j,k,iBlock)*TemperatureJumpAllowed)
         Temperature_VGB(Te_,i,j,k,iBlock) = AveragedTemperature + &
              DiffT*&
              SpecificHeat_VCB(TRad_,i,j,k,iBlock)/SpecificHeatTotal
         Temperature_VGB(TRad_,i,j,k,iBlock) = AveragedTemperature - &
              DiffT*&
              SpecificHeat_VCB(Te_,i,j,k,iBlock)/SpecificHeatTotal
      end do; end do; end do
    end subroutine check_temperature

  end subroutine advance_temperature

  !============================================================================

  subroutine update_conservative_energy(iBlock)

    use ModAdvance,  ONLY: State_VGB, p_
    use ModEnergy,   ONLY: calc_energy_cell
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK
    use ModPhysics,  ONLY: inv_gm1, No2Si_V, Si2No_V, UnitEnergyDens_, &
         UnitP_, cRadiationNo,UnitTemperature_
    use ModProcMH,   ONLY: iProc
    use ModSize,     ONLY: nI, nJ, nK
    use ModUser,     ONLY: user_material_properties

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: EinternalSi, Einternal, Gamma, PressureSi, TeSi
    logical :: DoReport = .true.
    character(len=*), parameter:: NameSub = 'update_conservative_energy'
    !--------------------------------------------------------------------------

    !!! This routine can currently not treat split electron/ion temperatures

    do k = 1,nK; do j = 1,nJ; do i = 1,nI

       if(UseTrad)then
          State_VGB(iErad,i,j,k,iBlock) = State_VGB(iErad,i,j,k,iBlock) &
               + SpecificHeat_VCB(Trad_,i,j,k,iBlock) &
               *(Temperature_VGB(Trad_,i,j,k,iBlock) &
               - TemperatureOld_VCB(Trad_,i,j,k,iBlock))

          if(UseTemperatureVariable)then
             ! Fix energy if it goes negative
             if(State_VGB(iErad,i,j,k,iBlock) < 0.0) then

                if(DoReport .or. &
                     Temperature_VGB(Trad_,i,j,k,iBlock) < 0.0)then

                   write(*,*)NameSub, ' WARNING: negative Erad=', &
                        State_VGB(iErad,i,j,k,iBlock)
                   write(*,*)NameSub, ' i,j,k,iBlock,iProc=',i,j,k,iBlock,iProc
                   write(*,*)NameSub, ' x, y, z=', x_BLK(i,j,k,iBlock), &
                        y_BLK(i,j,k,iBlock), z_BLK(i,j,k,iBlock)
                   write(*,*)NameSub, ' fix energy using Trad^n+1 =', &
                        Temperature_VGB(Trad_,i,j,k,iBlock)

                   ! There is no way to fix this
                   if(Temperature_VGB(Trad_,i,j,k,iBlock) < 0) &
                        call stop_mpi(NameSub//' negative Trad!!!')

                   ! Report only once
                   DoReport = .false.
                end if

                State_VGB(iErad,i,j,k,iBlock) = &
                     cRadiationNo*Temperature_VGB(Trad_,i,j,k,iBlock)**4
             end if
          else
             State_VGB(iErad,i,j,k,iBlock) &
                  = max(State_VGB(iErad,i,j,k,iBlock),EradMin)
          end if
       end if

       ! update material internal energy
       Einternal = inv_gm1*State_VGB(p_,i,j,k,iBlock) &
            + State_VGB(EintExtra_,i,j,k,iBlock) &
            + SpecificHeat_VCB(Te_,i,j,k,iBlock) &
            *(Temperature_VGB(Te_,i,j,k,iBlock) &
            - TemperatureOld_VCB(Te_,i,j,k,iBlock))

       EinternalSi = Einternal*No2Si_V(UnitEnergyDens_)
       if(EInternalSi <= 0.0)then
          TeSi = Temperature_VGB(Te_,i,j,k,iBlock)*No2Si_V(UnitTemperature_)
          if(TeSi <= 0.0)call stop_mpi('Negative temperature')
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               TeSiIn = TeSi, EinternalSiOut = EinternalSi, &
               PressureSiOut = PressureSi)
          Einternal = EInternalSi*Si2No_V(UnitEnergyDens_)
       else

          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               EinternalSiIn = EinternalSi, PressureSiOut = PressureSi)
       end if
       State_VGB(p_,i,j,k,iBlock) = PressureSi*Si2No_V(UnitP_)

       State_VGB(EintExtra_,i,j,k,iBlock) = &
            Einternal - inv_gm1*State_VGB(p_,i,j,k,iBlock)

    end do; end do; end do

    call calc_energy_cell(iBlock)

  end subroutine update_conservative_energy

  !============================================================================

  subroutine get_face_coef(iBlock)

    ! Auxiliary routine, which calculates face centered values for the heat
    ! conduction coefficients multiplied by face area divided by cell-to-cell
    ! center distance across the given face

    use ModGeometry, ONLY: dx_BLK, dy_BLK, dz_BLK, fAx_BLK, fAy_BLK, fAz_BLK, &
         vInv_CB, y_Blk, TypeGeometry, UseCovariant
    use ModCovariant, ONLY: FaceAreaI_DFB, FaceAreaJ_DFB, FaceAreaK_DFB
    use ModMain,     ONLY: x_, y_, z_
    use ModSize,     ONLY: nI, nJ, nK

    integer, intent(in) :: iBlock

    real :: AreaxInvDx, AreayInvDy, AreazInvDz, Dxyz_D(3)
    integer :: i, j, k, iT, iCond
    !--------------------------------------------------------------------------

    if(.not.UseCovariant)then

       AreaxInvDx = fAx_Blk(iBlock)/Dx_Blk(iBlock)
       AreayInvDy = fAy_Blk(iBlock)/Dy_Blk(iBlock)
       AreazInvDz = fAz_Blk(iBlock)/Dz_Blk(iBlock)

       Cond_VFX(1:nCond,1:nI+1,1:nJ,1:nK) = AreaxInvDx * 0.50 * &
            (HeatConductionCoef_IGB(1:nCond,0:nI  ,1:nJ  ,1:nK  ,iBlock) &
            +HeatConductionCoef_IGB(1:nCond,1:nI+1,1:nJ  ,1:nK  ,iBlock) )
       Cond_VFY(1:nCond,1:nI,1:nJ+1,1:nK) = AreayInvDy * 0.50 * &
            (HeatConductionCoef_IGB(1:nCond,1:nI  ,0:nJ  ,1:nK  ,iBlock) &
            +HeatConductionCoef_IGB(1:nCond,1:nI  ,1:nJ+1,1:nK  ,iBlock) )
       Cond_VFZ(1:nCond,1:nI,1:nJ,1:nK+1) = AreazInvDz * 0.50 * &
            (HeatConductionCoef_IGB(1:nCond,1:nI  ,1:nJ  ,0:nK  ,iBlock) &
            +HeatConductionCoef_IGB(1:nCond,1:nI  ,1:nJ  ,1:nK+1,iBlock) )

    elseif(TypeGeometry=='cartesian')then

       AreaxInvDx = FaceAreaI_DFB(x_,1,1,1,iBlock)/Dx_Blk(iBlock)
       AreayInvDy = FaceAreaJ_DFB(y_,1,1,1,iBlock)/Dy_Blk(iBlock)
       AreazInvDz = FaceAreaK_DFB(z_,1,1,1,iBlock)/Dz_Blk(iBlock)

       Cond_VFX(1:nCond,1:nI+1,1:nJ,1:nK) = AreaxInvDx * 0.50 * &
            (HeatConductionCoef_IGB(1:nCond,0:nI  ,1:nJ  ,1:nK  ,iBlock) &
            +HeatConductionCoef_IGB(1:nCond,1:nI+1,1:nJ  ,1:nK  ,iBlock) )
       Cond_VFY(1:nCond,1:nI,1:nJ+1,1:nK) = AreayInvDy * 0.50 * &
            (HeatConductionCoef_IGB(1:nCond,1:nI  ,0:nJ  ,1:nK  ,iBlock) &
            +HeatConductionCoef_IGB(1:nCond,1:nI  ,1:nJ+1,1:nK  ,iBlock) )
       Cond_VFZ(1:nCond,1:nI,1:nJ,1:nK+1) = AreazInvDz * 0.50 * &
            (HeatConductionCoef_IGB(1:nCond,1:nI  ,1:nJ  ,0:nK  ,iBlock) &
            +HeatConductionCoef_IGB(1:nCond,1:nI  ,1:nJ  ,1:nK+1,iBlock) )
    
    elseif(TypeGeometry=='rz')then

       AreaxInvDx = 1/Dx_Blk(iBlock)
       AreayInvDy = 1/Dy_Blk(iBlock)
     
       do k=1,nK; do j=1,nJ; do i=1,nI+1
          Cond_VFX(1:nCond,i,j,k) = AreaxInvDx * 0.50 *           &
            (HeatConductionCoef_IGB(1:nCond,i-1,j  ,k,iBlock)     &
            +HeatConductionCoef_IGB(1:nCond,i  ,j  ,k,iBlock) ) * &
            FaceAreaI_DFB(x_,i,j,k,iBlock)
       end do; end do; end do

       do k=1,nK; do j=1,nJ+1; do i=1,nI
          Cond_VFY(1:nCond,i,j,k) = AreayInvDy * 0.50 *           &
            (HeatConductionCoef_IGB(1:nCond,i  ,j-1,k,iBlock)     &
            +HeatConductionCoef_IGB(1:nCond,i  ,j  ,k,iBlock) ) * &
            FaceAreaJ_DFB(y_,i,j,k,iBlock)
       end do; end do; end do

       Cond_VFZ(1:nCond,1:nI,1:nJ,1:nK+1) = 0.0
    else
       call stop_mpi('General covariant geometry is expected soon')
    end if
  end subroutine get_face_coef

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
    if(IsFirstCgIteration)then
       call set_gray_outer_bcs(nVar,Temp_VGB,'float')
       IsFirstCgIteration = .false.
    else
       call set_gray_outer_bcs(nVar, Temp_VGB,'zero')
    end if
    pDotADotPPe = 0.0
    do iBlock = 1, nBlock
       if(unusedBlk(iBlock))CYCLE

       call get_heat_conduction(iBlock, nVar, &
            Temp_VGB(:,:,:,:,iBlock), Rhs_VCB(:,:,:,:,iBlock))
    end do

  end subroutine heat_conduction

  !============================================================================

  subroutine get_heat_conduction(iBlock, nVar, Temp_VG, Rhs_VC)
    use ModSize, ONLY: nI, nJ, nK, gcn
    use ModMain, ONLY: Dt
    use ModGeometry, ONLY: vInv_CB
    use ModParallel, ONLY: NOBLK, NeiLev

    integer, intent(in) :: iBlock, nVar
    real, intent(inout) :: &
         Temp_VG(nVar,1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn)
    real, intent(out)   :: Rhs_VC(nVar,nI,nJ,nK)

    integer :: i, j, k, iVar, iCond, iT
    real    :: Volume
    real    :: Coef_FX(1:nI+1,1:nJ,1:nK)
    real    :: Coef_FY(1:nI,1:nJ+1,1:nK)
    real    :: Coef_FZ(1:nI,1:nJ,1:nK+1)
    !--------------------------------------------------------------------------

    call get_face_coef(iBlock)

    Coef_FX = 0.50
    if(NeiLev(1,iBlock)==NOBLK)Coef_FX(1   ,:,:) = 1.0
    if(NeiLev(2,iBlock)==NOBLK)Coef_FX(nI+1,:,:) = 1.0

    Coef_FY = 0.50
    if(NeiLev(3,iBlock)==NOBLK)Coef_FY(:,1   ,:) = 1.0
    if(NeiLev(4,iBlock)==NOBLK)Coef_FY(:,nJ+1,:) = 1.0

    Coef_FZ = 0.50
    if(NeiLev(5,iBlock)==NOBLK)Coef_FZ(:,:,1   ) = 1.0
    if(NeiLev(6,iBlock)==NOBLK)Coef_FZ(:,:,nK+1) = 1.0

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       Volume = 1/vInv_CB(i,j,k,iBlock)

       ! Specific heat
       Rhs_VC(:,i,j,k) = SpecificHeat_VCB(:,i,j,k,iBlock) &
            *Temp_VG(:,i,j,k)*Volume/Dt

       pDotADotPPe = pDotADotPPe + sum(Rhs_VC(:,i,j,k)*Temp_VG(:,i,j,k))

       ! Heat conduction
       do iCond = 1, nCond
          iT = iCond_I(iCond)
          Rhs_VC(iT,i,j,k) = Rhs_VC(iT,i,j,k) + ( &
               + Cond_VFX(iCond,i+1,j,k) &
               *(Temp_VG(iT,i,j,k) - Temp_VG(iT,i+1,j,k)) &
               + Cond_VFX(iCond,i,j,k) &
               *(Temp_VG(iT,i,j,k) - Temp_VG(iT,i-1,j,k)) &
               + Cond_VFY(iCond,i,j+1,k) &
               *(Temp_VG(iT,i,j,k) - Temp_VG(iT,i,j+1,k)) &
               + Cond_VFY(iCond,i,j,k) &
               *(Temp_VG(iT,i,j,k) - Temp_VG(iT,i,j-1,k)) &
               + Cond_VFZ(iCond,i,j,k+1) &
               *(Temp_VG(iT,i,j,k) - Temp_VG(iT,i,j,k+1)) &
               + Cond_VFZ(iCond,i,j,k) &
               *(Temp_VG(iT,i,j,k) - Temp_VG(iT,i,j,k-1)) )

          pDotADotPPe = pDotADotPPe + &
                  Cond_VFX(iCond,i+1,j,k) * Coef_FX(i+1,j,k) *&
                 (Temp_VG(iT,i,j,k) - Temp_VG(iT,i+1,j,k))**2 + &
                  Cond_VFX(iCond,i  ,j,k) * Coef_FX(i  ,j,k) *&
                 (Temp_VG(iT,i,j,k) - Temp_VG(iT,i-1,j,k))**2 + &
                  Cond_VFY(iCond,i,j+1,k) * Coef_FY(i,j+1,k) *&
                 (Temp_VG(iT,i,j,k) - Temp_VG(iT,i,j+1,k))**2 + &
                  Cond_VFY(iCond,i,j  ,k) * Coef_FY(i,j  ,k) *&
                 (Temp_VG(iT,i,j,k) - Temp_VG(iT,i,j-1,k))**2 + &
                  Cond_VFZ(iCond,i,j,k+1) * Coef_FZ(i,j,k+1) *&
                 (Temp_VG(iT,i,j,k) - Temp_VG(iT,i,j,k+1))**2 + &
                  Cond_VFZ(iCond,i,j,k  ) * Coef_FZ(i,j,k  ) *&
                 (Temp_VG(iT,i,j,k) - Temp_VG(iT,i,j,k-1))**2 

       end do

       do iVar=2,nVar
          ! Energy exchange
          Rhs_VC(Te_,i,j,k)  = Rhs_VC(Te_,i,j,k) &
               + RelaxationCoef_VCB(iVar,i,j,k,iBlock) * Volume &
               * (Temp_VG(Te_,i,j,k) - Temp_VG(iVar,i,j,k) )

          pDotADotPPe = pDotADotPPe + &
                 RelaxationCoef_VCB(iVar,i,j,k,iBlock) * Volume &
               * (Temp_VG(Te_,i,j,k) - Temp_VG(iVar,i,j,k) )**2

          Rhs_VC(iVar,i,j,k) = Rhs_VC(iVar,i,j,k) &
               + RelaxationCoef_VCB(iVar,i,j,k,iBlock) * Volume &
               * (Temp_VG(iVar,i,j,k) - Temp_VG(Te_,i,j,k) )
       end do

    end do; end do; end do

  end subroutine get_heat_conduction

  !==========================================================================

  subroutine cg(matvec, nVar, Residual_VCB, Solution_VGB, IsInit, &
       Tolerance, TypeStop, Iter, DoTest, preconditioner)

    ! conjugated gradient for symmetric and positive definite matrix A

    use ModMain, ONLY: nBlock, unusedBLK
    use ModSize, ONLY: nI, nJ, nK, gcn, nBlk
    use ModProcMH,ONLY: iComm
    use ModMpi

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
         real, intent(in)    :: VectorIn_VCB(nVar,1:nI,1:nJ,1:nK,nBLK)
         real, intent(out)   :: VectorOut_VCB(nVar,1:nI,1:nJ,1:nK,nBLK)
       end subroutine preconditioner
    end interface


    integer, intent(in) :: nVar
    real, intent(inout) :: &      ! right hand side vector
         Residual_VCB(nVar,1:nI,1:nJ,1:nK,nBLK)
    real, intent(inout) :: &      ! initial guess / solution vector
         Solution_VGB(nVar,1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK)
    logical, intent(in) :: IsInit ! true if Solution_VGB contains initial guess
    real, intent(inout) :: Tolerance ! required / achieved residual

    character (len=3), intent(in) :: TypeStop
    !      Determine stopping criterion (||.|| denotes the 2-norm):
    !      typestop='rel'  -- relative stopping crit.: ||res|| <= Tol*||res0||
    !      typestop='abs'  -- absolute stopping crit.: ||res|| <= Tol
 
    integer, intent(out) :: Iter ! actual number of iterations
    logical, intent(in) :: DoTest

    integer, parameter :: MaxIter = 30000

    real :: rDotR, pDotADotP, rDotR0, rDotRMax, rDotRInv

    real, dimension(:,:,:,:,:), allocatable :: &
         P_VGB, ADotP_VCB, MMinusOneDotR_VCB

    integer :: i, j, k, iBlock,iError
    !--------------------------------------------------------------------------

    allocate( &
         P_VGB(nVar,1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK), &
         ADotP_VCB(nVar,1:nI,1:nJ,1:nK,nBLK) )

    do iBlock = 1, nBlock
       if(unusedBlk(iBlock))CYCLE
       P_VGB(:,:,:,:,iBlock) = 0.0
    end do

    if(present(preconditioner)) &
       allocate(MMinusOneDotR_VCB(nVar,1:nI,1:nJ,1:nK,nBLK))

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

    do Iter = 0, MaxIter

       if(present(preconditioner))then
          call preconditioner(nVar, Residual_VCB, MMinusOneDotR_VCB)
          select case(TypePreconditioner)
          case('blockjacobi','gs','dilu')
             call MPI_ALLREDUCE(rDotRPe, rDotR, 1, MPI_REAL, MPI_SUM, iComm, iError)
          case default
             rDotR = dot_product_mpi(0, Residual_VCB, MMinusOneDotR_VCB)
          end select
       else
          rDotR = dot_product_mpi(0, Residual_VCB, Residual_VCB)
       end if

       if(Iter == 0)then
          rDotR0 = rDotR

          if(TypeStop=='abs')then
             rDotRMax = Tolerance**2
          else
             rDotRMax = Tolerance**2 * rDotR0
          end if

          if(DoTest)write(*,*)'CG rDotR0, rDotRMax=', rDotR0, rDotRMax

          if(rDotR0 == 0.0) EXIT
       end if
  
       if(DoTest)write(*,*)'CG Iter, rDotR = ', Iter, rDotR

       if(rDotR <= rDotRMax) EXIT
       rDotRInv = 1/rDotR
       do iBlock = 1, nBlock
          if(unusedBlk(iBlock))CYCLE
          if(present(preconditioner))then
             P_VGB(:,1:nI,1:nJ,1:nK,iBlock) = P_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
                  + MMinusOneDotR_VCB(:,:,:,:,iBlock) * rDotRInv
          else
             P_VGB(:,1:nI,1:nJ,1:nK,iBlock) = P_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
                  + Residual_VCB(:,:,:,:,iBlock) * rDotRInv
          end if
       end do

       call matvec(nVar, P_VGB, aDotP_VCB)
       call MPI_ALLREDUCE(pDotADotPPe, pDotADotP, 1, MPI_REAL, MPI_SUM, iComm, iError)
      ! write(*,*)pDotADotP , dot_product_mpi(gcn, P_VGB(:,:,:,:,:), aDotP_VCB)
      ! pDotADotP = dot_product_mpi(gcn,P_VGB,aDotP_VCB)
       do iBlock = 1, nBlock
          if(unusedBlk(iBlock))CYCLE
          Residual_VCB(:,:,:,:,iBlock) = Residual_VCB(:,:,:,:,iBlock) &
               - aDotP_VCB(:,:,:,:,iBlock)/pDotADotP
       end do

!!$       write(*,*)dot_product_mpi(gcn,P_VGB,Residual_VCB)

       do iBlock = 1, nBlock
          if(unusedBlk(iBlock))CYCLE
          Solution_VGB(:,1:nI,1:nJ,1:nK,iBlock) = &
               Solution_VGB(:,1:nI,1:nJ,1:nK,iBlock) &
               + P_VGB(:,1:nI,1:nJ,1:nK,iBlock)/pDotADotP
       end do

    end do

    if(iTer > MaxIter)call stop_mpi("ModTemperature::cg did not converge")

    ! Calculate the achieved tolerance
    if(TypeStop=='abs')then
       Tolerance = sqrt(rDotR)
    else
       Tolerance = sqrt(rDotR/rDotR0)
    end if

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
    real, intent(in)    :: VectorIn_VCB(nVar,1:nI,1:nJ,1:nK,nBLK)
    real, intent(out)   :: VectorOut_VCB(nVar,1:nI,1:nJ,1:nK,nBLK)

    integer :: iVar, i, j , k, iBlock

    character(len=*), parameter :: NameSub = 'jacobi_preconditioner'
    !--------------------------------------------------------------------------

    rDotRPe = 0.0

    do iBlock = 1, nBlock
       if(unusedBlk(iBlock))CYCLE

       if(body_blk(iBlock))then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock))CYCLE
             do iVar = 1, nVar
                VectorOut_VCB(iVar,i,j,k,iBlock) = sum( &
                     VectorIn_VCB(:,i,j,k,iBlock) &
                     *JacPrec_VVCB(:,iVar,i,j,k,iBlock))
                rDotRPe = rDotRPe + VectorOut_VCB(iVar,i,j,k,iBlock) * &
                     VectorIn_VCB(iVar,i,j,k,iBlock)
             end do
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             do iVar = 1, nVar
                VectorOut_VCB(iVar,i,j,k,iBlock) = sum( &
                     VectorIn_VCB(:,i,j,k,iBlock) &
                     *JacPrec_VVCB(:,iVar,i,j,k,iBlock))
                rDotRPe = rDotRPe + VectorOut_VCB(iVar,i,j,k,iBlock) * &
                     VectorIn_VCB(iVar,i,j,k,iBlock)
             end do
          end do; end do; end do
       end if
    end do

  end subroutine jacobi_preconditioner

  !============================================================================

  subroutine get_jacobi_preconditioner(nVar)

    use ModGeometry, ONLY: true_cell, body_blk, Dx_Blk, Dy_Blk, Dz_Blk, &
         vInv_CB
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
       call get_face_coef(iBlock)
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

      use ModGeometry, ONLY: vInv_CB
      use ModMain,     ONLY: Dt

      real :: Volume
      !------------------------------------------------------------------------


      Volume = 1/vInv_CB(i,j,k,iBlock)

      ! Specific heat
      Diag_V(:) = SpecificHeat_VCB(:,i,j,k,iBlock)*Volume/Dt

      ! Heat conduction
      Diag_V(iCond_I) = Diag_V(iCond_I) &
           + Cond_VFX(:,i+1,j  ,k  )  &
           + Cond_VFX(:,i  ,j  ,k  )  &
           + Cond_VFY(:,i  ,j+1,k  )  &
           + Cond_VFY(:,i  ,j  ,k  )  &
           + Cond_VFZ(:,i  ,j  ,k+1)  &
           + Cond_VFZ(:,i  ,j  ,k  )  

      if(nVar==1) return

      ! Energy exchange
      Diag_V(Te_) = Diag_V(Te_) &
           + Volume*sum(RelaxationCoef_VCB(2:nVar,i,j,k,iBlock))
      Diag_V(2:nVar) = Diag_V(2:nVar) &
           + Volume*RelaxationCoef_VCB(2:nVar,i,j,k,iBlock)
      OffDiag_V(2:nVar) = -Volume*RelaxationCoef_VCB(2:nVar,i,j,k,iBlock)

    end subroutine get_diagonal_subblock

  end subroutine get_jacobi_preconditioner

  !============================================================================

  subroutine get_dilu_preconditioner(nVar)

    use ModGeometry, ONLY: body_blk,vInv_CB
    use ModMain, ONLY: nBlock, unusedBLK, Dt
    use ModSize, ONLY: nI, nJ, nK

    ! Inverts the nTemperature * nTemperature sub-matrices at the
    ! block diagonal of the operator of heat-conduction + relaxation

    integer, intent(in) :: nVar

    integer :: iVar, i, j , k, iBlock, iCond, iT

    real, dimension(nVar) :: Diag_V
    real, allocatable, dimension(:) :: OffDiag_V

    real :: Volume

    character(len=*), parameter :: NameSub = 'get_jacobi_preconditioner'
    !--------------------------------------------------------------------------

    if(nVar>1) allocate(OffDiag_V(2:nVar))
  
    do iBlock = 1, nBlock
       if(unusedBlk(iBlock))CYCLE
       if(body_blk(iBlock))call stop_mpi(&
            'D_ILU preconditioner does not work for body block')
       call get_face_coef(iBlock)
       !Fluxes across the block boundary are not accounted for

       !Save non-inverted diagonal matrices
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          Volume = 1/vInv_CB(i,j,k,iBlock)

          ! Heat conduction
          Diag_V(:) = SpecificHeat_VCB(:,i,j,k,iBlock)*Volume/Dt
          do iCond = 1, nCond
             iT = iCond_I(iCond)
             Diag_V(iT) = Diag_V(iT) &
                  + Cond_VFX(iCond,i+1,j  ,k  )  &
                  + Cond_VFX(iCond,i  ,j  ,k  )  &
                  + Cond_VFY(iCond,i  ,j+1,k  )  &
                  + Cond_VFY(iCond,i  ,j  ,k  )  &
                  + Cond_VFZ(iCond,i  ,j  ,k+1)  &
                  + Cond_VFZ(iCond,i  ,j  ,k  )  
          end do
          JacPrec_VVCB(:,:,i,j,k,iBlock) = 0.0
          if(nVar/=1) then

             ! Energy exchange

             Diag_V(Te_) = Diag_V(Te_) &
                  + Volume*sum(RelaxationCoef_VCB(2:nVar,i,j,k,iBlock))
             Diag_V(2:nVar) = Diag_V(2:nVar) &
                  + Volume*RelaxationCoef_VCB(2:nVar,i,j,k,iBlock)
             JacPrec_VVCB(Te_,2:nVar,i,j,k,iBlock) = &
                  -Volume*RelaxationCoef_VCB(2:nVar,i,j,k,iBlock)
             JacPrec_VVCB(2:nVar,Te_,i,j,k,iBlock) = &
                  -Volume*RelaxationCoef_VCB(2:nVar,i,j,k,iBlock)
          end if
          do iVar=1,nVar
             JacPrec_VVCB(iVar,iVar,i,j,k,iBlock) = Diag_V(iVar)
          end do
       end do;end do;end do
      
       select case(nVar)
       case(1)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             call block_jacobi_1t
             call modify_upper_cell_jacobian
          end do; end do; end do
       case(2)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             call block_jacobi_2t
             call modify_upper_cell_jacobian
          end do; end do; end do
       case(3)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             call block_jacobi_3t
             call modify_upper_cell_jacobian
          end do; end do; end do
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
      do iVar = 1,nVar
         Diag_V(iVar) = JacPrec_VVCB(iVar,iVar,i,j,k,iBlock)
      end do
      if(nVar==1) return
      OffDiag_V(2:nVar) = JacPrec_VVCB(Te_,2:nVar,i,j,k,iBlock)
    end subroutine get_diagonal_subblock
    !=========================================================================!
    subroutine modify_upper_cell_jacobian
       integer :: jVar
       real :: LocI_V(nVar),LocJ_V(nVar),LocK_V(nVar)
       !------------------------------------------------------------------
       LocI_V = 0.0; LocJ_V = 0.0; LocK_V = 0.0
       !return
       do iCond = 1,nCond
          iT = iCond_I(iCond)
          LocI_V(iT) = Cond_VFX(iCond,i+1,j,k)
          LocJ_V(iT) = Cond_VFY(iCond,i,j+1,k)
          LocK_V(iT) = Cond_VFZ(iCond,i,j,k+1)
       end do
       if(i<nI)then
          do jVar=1,nVar; do iVar = 1,nVar
             JacPrec_VVCB(iVar,jVar,i+1,j,k,iBlock) = &
                  JacPrec_VVCB(iVar,jVar,i+1,j,k,iBlock) - &
                  JacPrec_VVCB(iVar,jVar,i,j,k,iBlock) &
                  *LocI_V(iVar)*LocI_V(jVar)
          end do;end do
       end if
       if(j<nJ)then
          do jVar=1,nVar; do iVar = 1,nVar
             JacPrec_VVCB(iVar,jVar,i,j+1,k,iBlock) = &
                  JacPrec_VVCB(iVar,jVar,i,j+1,k,iBlock) - &
                  JacPrec_VVCB(iVar,jVar,i,j,k,iBlock) &
                  *LocJ_V(iVar)*LocJ_V(jVar)
          end do;end do
       end if
       if(k<nK)then
          do jVar=1,nVar; do iVar = 1,nVar
             JacPrec_VVCB(iVar,jVar,i,j,k+1,iBlock) = &
                  JacPrec_VVCB(iVar,jVar,i,j,k+1,iBlock) - &
                  JacPrec_VVCB(iVar,jVar,i,j,k,iBlock) &
                  *LocK_V(iVar)*LocK_V(jVar)
          end do;end do
       end if
    end subroutine modify_upper_cell_jacobian
  end subroutine get_dilu_preconditioner
  !============================================================================
  subroutine dilu_preconditioner(nVar, VectorIn_VCB, VectorOut_VCB)

    ! Calculate VectorOut = M^{-1} \cdot VectorIn where the stored
    ! preconditioner matrix, M^{-1}, should be symmetric positive
    ! and is approximately decomposed as 
    ! (D+U)^{-1}\cdot D \cdot (D+L)^{-1)
    ! In symmetric Gauss-Seidel preconditioner the matrix D (in fact 
    ! we need D^{-1} is the same as in the block Jacobi preconditioner
    ! In D_ILU preconditioner the D^{-1} is somewhat modified 

    use ModGeometry, ONLY: true_cell, body_blk
    use ModMain, ONLY: nBlock, unusedBLK
    use ModSize, ONLY: nI, nJ, nK, nBlk

    integer, intent(in) :: nVar
    real, intent(in)    :: VectorIn_VCB(nVar,1:nI,1:nJ,1:nK,nBLK)
    real, intent(out)   :: VectorOut_VCB(nVar,1:nI,1:nJ,1:nK,nBLK)

    integer :: iVar, i, j , k, iBlock
    real :: Loc_V(nVar)
    real :: Sol_VG(nVar,0:nI+1,0:nJ+1,0:nK+1)

    character(len=*), parameter :: NameSub = 'dilu_preconditioner'
    !--------------------------------------------------------------------------
    Sol_VG = 0.0 ! To set zero values in the ghostcell
    rDotRPe = 0.0
    do iBlock = 1, nBlock
       if(unusedBlk(iBlock))CYCLE

       if(body_blk(iBlock))then
          call stop_mpi('DILU preconditioner does not work for body block')
       else
          call get_face_coef(iBlock)

          !Stage "UP", fluxes from below are accounted for
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             Loc_V = VectorIn_VCB(:,i,j,k,iBlock)
             Loc_V(iCond_I) = Loc_V(iCond_I) &
                  + Sol_VG(iCond_I,i-1,j,k)*Cond_VFX(:,i,j,k) &
                  + Sol_VG(iCond_I,i,j-1,k)*Cond_VFY(:,i,j,k) &
                  + Sol_VG(iCond_I,i,j,k-1)*Cond_VFZ(:,i,j,k)
             do iVar = 1, nVar
                Sol_VG(iVar,i,j,k) = sum( &
                      Loc_V*JacPrec_VVCB(:,iVar,i,j,k,iBlock))
                rDotRPe = rDotRPe + Loc_V(iVar) * Sol_VG(iVar,i,j,k)
             end do
          end do; end do; end do

          !Stage "DOWN", fluxes from above are accounted for
          do k = nK,1,-1; do j = nJ,1,-1; do i = nI,1,-1
             Loc_V = 0.0
             Loc_V(iCond_I) = Loc_V(iCond_I) &
                  + Sol_VG(iCond_I,i+1,j,k)*Cond_VFX(:,i+1,j,k) &
                  + Sol_VG(iCond_I,i,j+1,k)*Cond_VFY(:,i,j+1,k) &
                  + Sol_VG(iCond_I,i,j,k+1)*Cond_VFZ(:,i,j,k+1)
             do iVar = 1, nVar
                Sol_VG(iVar,i,j,k) = Sol_VG(iVar,i,j,k) &
                     + sum(Loc_V*JacPrec_VVCB(:,iVar,i,j,k,iBlock))
             end do
          end do; end do; end do

          VectorOut_VCB(:,:,:,:,iBlock) = Sol_VG(:,1:nI,1:nJ,1:nK)
       end if
    end do
  end subroutine dilu_preconditioner

  !============================================================================

  subroutine set_gray_outer_bcs(nVar, Var_VGB,TypeBc)

    use ModMain, ONLY: nI, nJ, nK, nBlk, nBlock, UnusedBlk,TypeBc_I
    use ModParallel, ONLY: NOBLK, NeiLev

    integer, intent(in) :: nVar
    real, intent(inout) :: Var_VGB(nVar,-1:nI+2,-1:nJ+2,-1:nK+2,nBlk)
    character(LEN=*), intent(in) :: TypeBc

    integer :: iBlock
    logical :: IsZero_S(1:6)
    !--------------------------------------------------------------------------

    IsZero_S = TypeBc=='zero' .and. TypeBc_I(1:6)/='reflect'

    do iBlock = 1, nBlock
       if(UnusedBlk(iBlock))CYCLE

       if(NeiLev(1,iBlock) == NOBLK)then
          if(IsZero_S(1))then
             Var_VGB(:,0,:,:,iBlock) = 0.0
          else
             Var_VGB(:,0,:,:,iBlock) = Var_VGB(:,1,:,:,iBlock)
          end if
       end if
       if(NeiLev(2,iBlock) == NOBLK)then
          if(IsZero_S(2))then
             Var_VGB(:,nI+1,:,:,iBlock) = 0.0
          else
             Var_VGB(:,nI+1,:,:,iBlock) = Var_VGB(:,nI,:,:,iBlock)
          end if
       end if
       if(NeiLev(3,iBlock) == NOBLK)then
          if(IsZero_S(3))then
             Var_VGB(:,:,0,:,iBlock) = 0.0
          else
             Var_VGB(:,:,0,:,iBlock) = Var_VGB(:,:,1,:,iBlock)
          end if
       end if
       if(NeiLev(4,iBlock) == NOBLK)then
          if(IsZero_S(4))then
             Var_VGB(:,:,nJ+1,:,iBlock) = 0.0
          else
             Var_VGB(:,:,nJ+1,:,iBlock) = Var_VGB(:,:,nJ,:,iBlock)
          end if
       end if
       if(NeiLev(5,iBlock) == NOBLK)then
          if(IsZero_S(5))then
             Var_VGB(:,:,:,0,iBlock) = 0.0
          else
             Var_VGB(:,:,:,0,iBlock) =  Var_VGB(:,:,:,1,iBlock)
          end if
       end if
       if(NeiLev(6,iBlock) == NOBLK)then
          if(IsZero_S(6))then
             Var_VGB(:,:,:,nK+1,iBlock) = 0.0
          else
             Var_VGB(:,:,:,nK+1,iBlock) = Var_VGB(:,:,:,nK,iBlock)
          end if
       end if
    end do

  end subroutine set_gray_outer_bcs

  !============================================================================

  subroutine mbilu_preconditioner(nVar, VectorIn_VCB, VectorOut_VCB)

    use ModLinearSolver, ONLY: Uhepta, Lhepta
    use ModMain, ONLY: nI, nJ, nK, nBlock, nBlk, unusedBlk, nIJK

    integer, intent(in) :: nVar
    real, intent(in)    :: VectorIn_VCB(nVar,nI,nJ,nK,nBlk)
    real, intent(out)   :: VectorOut_VCB(nVar,nI,nJ,nK,nBlk)

    integer :: iBlock
    !--------------------------------------------------------------------------

    do iBlock = 1, nBlock
       if(unusedBlk(iBlock))CYCLE

       VectorOut_VCB(:,:,:,:,iBlock) = VectorIn_VCB(:,:,:,:,iBlock)

       ! Rhs --> P_L.Rhs, where P_L=U^{-1}.L^{-1}
       call Lhepta(nIJK, nVar, nI, nI*nJ, &
            VectorOut_VCB(1,1,1,1,iBlock), &
            Prec_VVCIB(1,1,1,1,1,1,iBlock), &   ! Main diagonal
            Prec_VVCIB(1,1,1,1,1,2,iBlock), &   ! -i diagonal
            Prec_VVCIB(1,1,1,1,1,4,iBlock), &   ! -j
            Prec_VVCIB(1,1,1,1,1,6,iBlock) )    ! -k
       call Uhepta(.true., nIJK, nVar, nI, nI*nJ, &
            VectorOut_VCB(1,1,1,1,iBlock), &
            Prec_VVCIB(1,1,1,1,1,3,iBlock), &   ! +i diagonal
            Prec_VVCIB(1,1,1,1,1,5,iBlock), &   ! +j
            Prec_VVCIB(1,1,1,1,1,7,iBlock) )    ! +k
    end do

  end subroutine mbilu_preconditioner

  !============================================================================

  subroutine get_mbilu_preconditioner(nVar)

    use ModLinearSolver, ONLY: prehepta
    use ModMain, ONLY: nI, nJ, nK, nIJK, nBlock, unusedBlk

    integer, intent(in) :: nVar

    integer :: iBlock
    !--------------------------------------------------------------------------
    
    do iBlock = 1, nBlock
       if(unusedBlk(iBlock))CYCLE

       call get_jacobian(iBlock, nVar, Prec_VVCIB(:,:,:,:,:,:,iBlock))

       ! Preconditioning: Jacobian --> LU
       call prehepta(nIJK, nVar, nI, nI*nJ, 0.0, &
            Prec_VVCIB(1,1,1,1,1,1,iBlock), &
            Prec_VVCIB(1,1,1,1,1,2,iBlock), &
            Prec_VVCIB(1,1,1,1,1,3,iBlock), &
            Prec_VVCIB(1,1,1,1,1,4,iBlock), &
            Prec_VVCIB(1,1,1,1,1,5,iBlock), &
            Prec_VVCIB(1,1,1,1,1,6,iBlock), &
            Prec_VVCIB(1,1,1,1,1,7,iBlock) )
    end do

  end subroutine get_mbilu_preconditioner

  !============================================================================

  subroutine get_jacobian(iBlock, nVar, Jacobian_VVCI)

    use ModGeometry, ONLY: vInv_CB
    use ModMain,     ONLY: nI, nJ, nK, nDim, Dt

    integer, intent(in) :: iBlock, nVar
    real, intent(out) :: Jacobian_VVCI(nVar,nVar,nI,nJ,nK,nstencil)

    integer :: iVar, i, j, k, iDim, iCond
    real :: DiffLeft, DiffRight
    real :: Volume, AreaInvDxyz_FD(1:nI+1,1:nJ+1,1:nK+1,nDim)
    !--------------------------------------------------------------------------

    Jacobian_VVCI(:,:,:,:,:,:) = 0.0

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       Volume = 1.0/vInv_CB(i,j,k,iBlock)
    
       ! specific heat
       do iVar = 1, nVar
          Jacobian_VVCI(iVar,iVar,i,j,k,1) = &
               SpecificHeat_VCB(iVar,i,j,k,iBlock)*Volume/Dt
       end do

       if(nVar == 1)CYCLE

       ! energy exchange
       Jacobian_VVCI(Te_,Te_,i,j,k,1) = Jacobian_VVCI(Te_,Te_,i,j,k,1) &
            + Volume*sum(RelaxationCoef_VCB(2:nVar,i,j,k,iBlock))
       do iVar = 2, nVar
          Jacobian_VVCI(iVar,iVar,i,j,k,1) = &
               Jacobian_VVCI(iVar,iVar,i,j,k,1) &
               + Volume*RelaxationCoef_VCB(iVar,i,j,k,iBlock)
       end do
       Jacobian_VVCI(Te_,2:nVar,i,j,k,1) = &
            -Volume*RelaxationCoef_VCB(2:nVar,i,j,k,iBlock)
       Jacobian_VVCI(2:nVar,Te_,i,j,k,1) = Jacobian_VVCI(Te_,2:nVar,i,j,k,1)
    end do; end do; end do

    ! Heat conduction and diffusion contribution
    call get_face_coef(iBlock)

    do iDim = 1, nDim
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          do iCond = 1, nCond
             iVar = iCond_I(iCond)
             select case(iDim)
             case(1)
                ! ignore radiative flux at the boundaries of the AMR blocks
                DiffLeft  = Cond_VFX(iCond,i,j,k)
                if(i==1)DiffLeft = 0.0
                DiffRight = Cond_VFX(iCond,i+1,j,k) 
                if(i==nI)DiffRight = 0.0
             case(2) 
                DiffLeft  = Cond_VFY(iCond,i,j,k)
                if(j==1) DiffLeft = 0.0
                DiffRight = Cond_VFY(iCond,i,j+1,k) 
                if(j==nJ) DiffRight = 0.0
             case(3)
                DiffLeft  = Cond_VFZ(iCond,i,j,k)
                if(k==1) DiffLeft = 0.0
                DiffRight = Cond_VFZ(iCond,i,j,k+1) 
                if(k==nK) DiffRight = 0.0
             end select
             Jacobian_VVCI(iVar,iVar,i,j,k,1) = &
                  Jacobian_VVCI(iVar,iVar,i,j,k,1) + DiffLeft  + DiffRight
             Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim) = -DiffLeft
             Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim+1) = -DiffRight
          end do
       end do; end do; end do
    end do

  end subroutine get_jacobian

end module ModTemperature
