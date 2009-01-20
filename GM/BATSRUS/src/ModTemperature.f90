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
  public :: calc_source_compression
  public :: advance_temperature

  Logical, public :: UseTemperatureDiffusion = .false.

  ! Parameters for radiation flux limiter
  logical,           public :: UseRadFluxLimiter  = .false.
  character(len=20), public :: TypeRadFluxLimiter = 'larsen'

  ! Coefficients for the temperature diffusion model
  real, allocatable :: RelaxationCoef_VCB(:,:,:,:,:)
  real, allocatable :: SpecificHeat_VCB(:,:,:,:,:)
  real, allocatable :: HeatConductionCoef_IGB(:,:,:,:,:)

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
  real :: MaxErrorResidual = 1e-5

  ! Right-hand side of non-linear coupled system of temperature equations
  real, allocatable :: Source_VCB(:,:,:,:,:)

  ! Store block Jacobi preconditioner
  real, allocatable :: JacPrec_VVCB(:,:,:,:,:,:)

contains

  !============================================================================

  subroutine read_temperature_param(NameCommand)

    use ModMain, ONLY: UseGrayDiffusion, UseHeatConduction
    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand

    character(len=*), parameter :: NameSub = 'read_temperature_param'
    !--------------------------------------------------------------------------

    select case(NameCommand)
    case("#IMPLICITTEMPERATURE")
       call read_var('UseTemperatureDiffusion', UseTemperatureDiffusion)
       if(UseTemperatureDiffusion)then
          call read_var('UseTemperatureVariable', UseTemperatureVariable)
          call read_var('TypeStopCriterion', TypeStopCriterion)
          call read_var('MaxErrorResidual', MaxErrorResidual)
       end if
    case("#HEATCONDUCTION")
       call read_var('UseHeatConduction', UseHeatConduction)

       call stop_mpi(NameSub// &
            ': isotropic heat conduction can not yet be activated')
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
    use ModVarIndexes, ONLY: NameVar_V

    integer :: iVar, nVar

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
            JacPrec_VVCB(nTemperature,nTemperature,nI,nJ,nK,nBLK), &
            Source_VCB(nTemperature,nI,nJ,nK,nBlk) )

       if(nCond > 0) &
            allocate( &
            HeatConductionCoef_IGB(nCond,-1:nI+2,-1:nJ+2,-1:nK+2,nBlk), &
            iCond_I(nCond) )

       if(nTemperature > 1) &
            allocate(RelaxationCoef_VCB(2:nTemperature,nI,nJ,nK,nBlk))
    end if

    ! Select which temperature variables involve "heat conduction"
    if(UseHeatConduction) iCond_I(1) = Te_
    if(UseGrayDiffusion)  iCond_I(nCond) = Trad_

  end subroutine init_temperature_diffusion

  !============================================================================

  subroutine set_frozen_coefficients

    use ModAdvance,  ONLY: State_VGB
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK
    use ModMain,     ONLY: nI, nJ, nK, nBlock, unusedBlk
    use ModPhysics,  ONLY: Si2No_V, UnitX_, UnitEnergyDens_, &
         UnitTemperature_, cRadiationNo, Clight
    use ModProcMH,   ONLY: iProc
    use ModUser,     ONLY: user_material_properties

    integer :: i, j, k, iBlock
    real :: PlanckOpacitySi, RosselandMeanOpacitySi
    real :: PlanckOpacity, RosselandMeanOpacity
    real :: CvSi, Cv, TeSi, Te, Trad, DiffRad

    character(len=*), parameter:: NameSub = 'set_frozen_coefficients'
    !--------------------------------------------------------------------------

    ! set the temperature variables
    !!! ion temperature can not yet be set
    do iBlock = 1, nBlock
       if(unusedBlk(iBlock)) CYCLE

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               TeSiOut = TeSi)

          Te = TeSi*Si2No_V(UnitTemperature_)

          if(UseTemperatureVariable)then
             Temperature_VGB(Te_,i,j,k,iBlock) = Te
             if(UseTrad)then
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

    if(UseRadFluxLimiter)then
       ! The radiation flux limiters use the gradient of the
       ! radiation temperature and therefor one layer of ghost cells
       ! needs to be correct. We massage pass all temperatures.

       ! DoOneLayer, DoFacesOnly, No UseMonoteRestrict
       call message_pass_cells8(.true., .true., .false., nTemperature, &
            Temperature_VGB)

       call set_gray_outer_bcs(nTemperature, Temperature_VGB)
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
               CvSiOut = CvSi, TeSiOut = TeSi)

          RosselandMeanOpacity = RosselandMeanOpacitySi/Si2No_V(UnitX_)
          PlanckOpacity = PlanckOpacitySi/Si2No_V(UnitX_)
          Cv = CvSi*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)
          Te = TeSi*Si2No_V(UnitTemperature_)
          call get_radiation_diffusion_coef

          if(UseTemperatureVariable)then
             SpecificHeat_VCB(Te_,i,j,k,iBlock) = Cv

             Trad = Temperature_VGB(Trad_,i,j,k,iBlock)
             SpecificHeat_VCB(Trad_,i,j,k,iBlock) = 4.0*cRadiationNo*Trad**3
             RelaxationCoef_VCB(Trad_,i,j,k,iBlock) = Clight*PlanckOpacity &
                  *cRadiationNo*(Te+Trad)*(Te**2+Trad**2)
             HeatConductionCoef_IGB(1,i,j,k,iBlock) = &
                  DiffRad*4.0*cRadiationNo*Trad**3
          else
             SpecificHeat_VCB(aTe4_,i,j,k,iBlock) = Cv/(4.0*cRadiationNo*Te**3)
             SpecificHeat_VCB(aTrad4_,i,j,k,iBlock) = 1.0
             RelaxationCoef_VCB(aTrad4_,i,j,k,iBlock) = Clight*PlanckOpacity
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
    call set_gray_outer_bcs(nCond, HeatConductionCoef_IGB)

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

  subroutine calc_source_compression(iBlock)

    use ModAdvance,    ONLY: State_VGB, Source_VC, &
         uDotArea_XI, uDotArea_YI, uDotArea_ZI
    use ModGeometry,   ONLY: vInv_CB, y_BLK, IsCylindrical
    use ModNodes,      ONLY: NodeY_NB
    use ModSize,       ONLY: nI, nJ, nK
    use ModVarIndexes, ONLY: Energy_

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: vInv, DivU, RadCompression
    character(len=*), parameter :: NameSub = "calc_source_compression"
    !--------------------------------------------------------------------------

    do k=1,nK; do j=1,nJ; do i=1,nI

       if(IsCylindrical)then
          ! Multiply volume with radius (=Y) at cell center
          ! -> divide inverse volume
          vInv = vInv_CB(i,j,k,iBlock)/abs(y_BLK(i,j,k,iBlock))
          DivU = (uDotArea_XI(i+1,j,k,1) - uDotArea_XI(i,j,k,1)) &
               *abs(y_BLK(i,j,k,iBlock)) &
               + uDotArea_YI(i,j+1,k,1)*abs(NodeY_NB(i,j+1,k,iBlock)) &
               - uDotArea_YI(i,j,k,1)*abs(NodeY_NB(i,j,k,iBlock))
       else
          vInv = vInv_CB(i,j,k,iBlock)
          DivU = uDotArea_XI(i+1,j,k,1) - uDotArea_XI(i,j,k,1) &
               + uDotArea_YI(i,j+1,k,1) - uDotArea_YI(i,j,k,1) &
               + uDotArea_ZI(i,j,k+1,1) - uDotArea_ZI(i,j,k,1) 
       end if

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

  end subroutine calc_source_compression

  !============================================================================

  subroutine advance_temperature

    use ModAdvance,  ONLY: State_VGB
    use ModGeometry, ONLY: vInv_CB, y_BLK, IsCylindrical
    use ModMain,     ONLY: nI, nJ, nK, nBlock, unusedBLK, Dt

    integer :: i, j, k, iBlock
    integer :: Iter
    logical :: DoTestKrylov, DoTestKrylovMe

    character(len=*), parameter :: NameSub = 'advance_temperature'
    !--------------------------------------------------------------------------

    call set_oktest('krylov', DoTestKrylov, DoTestKrylovMe)

    ! prevent doing useless updates
    if(Dt == 0.0)return

    call set_frozen_coefficients

    do iBlock = 1, nBlock
       if(unusedBlk(iBlock)) CYCLE

       do k = 1,nK; do j = 1,nJ; do i = 1,nI
          Source_VCB(:,i,j,k,iBlock) = SpecificHeat_VCB(:,i,j,k,iBlock) &
               *Temperature_VGB(:,i,j,k,iBlock)/(vInv_CB(i,j,k,iBlock)*Dt)
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
         .true., MaxErrorResidual, TypeStopCriterion, &
         Iter, jacobi_preconditioner)

    if(DoTestKrylovMe)write(*,*)NameSub,': Number of iterations =',Iter

    do iBlock = 1, nBlock
       if(unusedBlk(iBlock)) CYCLE

       call update_conservative_energy(iBlock)
    end do

  end subroutine advance_temperature

  !============================================================================

  subroutine update_conservative_energy(iBlock)

    use ModAdvance,  ONLY: State_VGB, p_
    use ModEnergy,   ONLY: calc_energy_cell
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK
    use ModPhysics,  ONLY: inv_gm1, No2Si_V, Si2No_V, UnitEnergyDens_, &
         UnitP_, cRadiationNo
    use ModProcMH,   ONLY: iProc
    use ModSize,     ONLY: nI, nJ, nK
    use ModUser,     ONLY: user_material_properties

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: EinternalSi, Einternal, Gamma, PressureSi
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

          ! Fix energy if it goes negative
          if(State_VGB(iErad,i,j,k,iBlock) < 0.0) then

             if(DoReport .or. Temperature_VGB(Trad_,i,j,k,iBlock) < 0.0 &
                  .or. .not.UseTemperatureVariable)then

                write(*,*)NameSub, ' WARNING: negative Erad=', &
                     State_VGB(iErad,i,j,k,iBlock)
                write(*,*)NameSub, ' i,j,k,iBlock,iProc=', i, j, k,iBlock,iProc
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
       end if

       ! update material internal energy
       Einternal = inv_gm1*State_VGB(p_,i,j,k,iBlock) &
            + State_VGB(EintExtra_,i,j,k,iBlock) &
            + SpecificHeat_VCB(Te_,i,j,k,iBlock) &
            *(Temperature_VGB(Te_,i,j,k,iBlock) &
            - TemperatureOld_VCB(Te_,i,j,k,iBlock))

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
         vInv_CB, y_Blk, IsCylindrical
    use ModMain,     ONLY: x_, y_, z_, Dt
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

    do iBlock = 1, nBlock
       if(unusedBlk(iBlock))CYCLE
       P_VGB(:,:,:,:,iBlock) = 0.0
    end do
    do iBlock = 1, nBlock
       if(unusedBlk(iBlock))CYCLE
       ADotP_VCB(:,:,:,:,iBlock) = 0.0
    end do

    if(present(preconditioner))then
       allocate(MMinusOneDotR_VCB(nVar,1:nI,1:nJ,1:nK,nBLK))

       do iBlock = 1, nBlock
          if(unusedBlk(iBlock))CYCLE
          MMinusOneDotR_VCB(:,:,:,:,iBlock) = 0.0
       end do
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
       end do
       do iBlock = 1, nBlock
          if(unusedBlk(iBlock))CYCLE
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
                     VectorIn_VCB(:,i,j,k,iBlock) &
                     *JacPrec_VVCB(:,iVar,i,j,k,iBlock))
             end do
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             do iVar = 1, nVar
                VectorOut_VCB(iVar,i,j,k,iBlock) = sum( &
                     VectorIn_VCB(:,i,j,k,iBlock) &
                     *JacPrec_VVCB(:,iVar,i,j,k,iBlock))
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

      use ModGeometry, ONLY: y_Blk, IsCylindrical
      use ModMain,     ONLY: Dt
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

  !============================================================================

  subroutine set_gray_outer_bcs(nVar, Var_VGB)

    use ModMain, ONLY: nI, nJ, nK, nBlk, nBlock, UnusedBlk
    use ModParallel, ONLY: NOBLK, NeiLev

    integer, intent(in):: nVar
    real, intent(inout):: Var_VGB(nVar,-1:nI+2,-1:nJ+2,-1:nK+2,nBlk)
    
    integer :: iBlock
    !--------------------------------------------------------------------------

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

end module ModTemperature
