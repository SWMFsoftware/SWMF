!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModResistivity

  use BATL_lib, ONLY: &
       test_start, test_stop, iTest, jTest, kTest, iBlockTest

  ! Resistivity related variables and methods

  use ModSize,   ONLY: MaxBlock
  use BATL_size, ONLY: nDim, Maxdim, nI, nJ, nK, &
       MinI, MaxI, MinJ, MaxJ, MinK, MaxK

  use ModImplicit, ONLY: UseSemiHallResist, UseSemiResistivity

  implicit none
  save

  private ! except

  public:: read_resistivity_param
  public:: init_mod_resistivity
  public:: set_resistivity
  public:: calc_resistivity_source
  public:: calc_heat_exchange
  public:: init_impl_resistivity
  public:: init_impl_hall_resist
  public:: get_impl_resistivity_state
  public:: get_resistivity_rhs
  public:: add_jacobian_resistivity
  public:: add_jacobian_hall_resist
  public:: update_impl_resistivity

  logical, public           :: UseResistivity   = .false.
  logical, public           :: UseResistiveFlux = .false.
  logical, public           :: UseHeatExchange  = .true.
  character(len=30), public :: TypeResistivity='none'
  real, public, allocatable :: Eta_GB(:,:,:,:)
  real, public              :: Eta0=0.0, Eta0Si=0.0

  real, public              :: Si2NoEta

  ! Local variables
  logical :: UseJouleHeating  = .true.
  logical :: DoMessagePassResistivity = .false.
  logical :: UseCentralDifference = .false.            

  real :: EtaPerpSpitzerSi = 0.0
  real :: CoulombLogarithm = 20.0
  real :: Eta0AnomSi=0.0, Eta0Anom
  real :: EtaMaxAnomSi=0.0, EtaMaxAnom
  real :: JcritAnomSi=1.0, JcritInv
  real :: EtaCoeff = 0.0
  real :: JoverBCrit = 0.0

  ! Description of the region where resistivity is used
  character(len=200):: StringResistRegion ='none'

  ! Indexes of regions defined with the #REGION commands
  integer, allocatable:: iRegionResist_I(:)

  ! resistivity pre-multiplied by the face area vector
  real, allocatable :: Eta_DFDB(:,:,:,:,:,:)

  ! B/ne vector for Hall Resistivity
  real, allocatable :: Bne_DFDB(:,:,:,:,:,:)

  ! Whistler diffusion coefficients.
  real, allocatable :: WhistlerCoeff_FDB(:,:,:,:,:)

  ! Named indices for semi-implicit variables
  integer, public, parameter :: BxImpl_ = 1, ByImpl_ = 2, BzImpl_ = 3

  ! Include resistive fluxes at all?
  logical :: DoResistiveFlux = .true.

  character(len=*), private, parameter :: NameMod = 'ModResistivity'
contains
  !============================================================================
  subroutine read_resistivity_param(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'read_resistivity_param'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    select case(NameCommand)
    case('#MESSAGEPASSRESISTIVITY')
       call read_var("DoMessagePassResistivity", DoMessagePassResistivity)
    case("#RESISTIVITYOPTIONS")
       call read_var("UseResistiveFlux", DoResistiveFlux)
       call read_var("UseJouleHeating",  UseJouleHeating)
       call read_var("UseHeatExchange",  UseHeatExchange)
    case("#RESISTIVITY")
       call read_var('UseResistivity',UseResistivity)
       if(UseResistivity)then
          call read_var('TypeResistivity', TypeResistivity, IsLowerCase=.true.)
          select case(TypeResistivity)
          case('spitzer')
             call read_var('CoulombLogarithm', CoulombLogarithm)
          case('user')
             call read_var('Eta0Si',Eta0Si)
          case('constant')
             call read_var('Eta0Si',Eta0Si)
          case('anomalous')
             call read_var('Eta0Si',       Eta0Si)
             call read_var('Eta0AnomSi',   Eta0AnomSi)
             call read_var('EtaMaxAnomSi', EtaMaxAnomSi)
             call read_var('JcritAnomSi',  JcritAnomSi)
          case('raeder')
             call read_var('EtaCoeff', EtaCoeff)
             call read_var('JoverBCrit', JoverBCrit)
          case default
             call stop_mpi(NameSub//': unknown TypeResistivity='&
                  //TypeResistivity)
          end select
       end if
    case('#RESISTIVITYSCHEME')
       call read_var('UseCentralDifference', UseCentralDifference)
    case("#RESISTIVITYREGION", "#RESISTIVEREGION")
       call read_var('StringResistRegion', StringResistRegion)
    case default
       call stop_mpi(NameSub//': unknown command='//NameCommand)
    end select

    call test_stop(NameSub, DoTest)
  end subroutine read_resistivity_param
  !============================================================================

  subroutine init_mod_resistivity

    use ModPhysics, ONLY: Si2No_V, UnitX_, UnitT_, UnitJ_
    use ModConst,   ONLY: cLightSpeed, cElectronCharge, &
         cElectronMass, cEps, cBoltzmann, cTwoPi
    use BATL_lib,   ONLY: get_region_indexes

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'init_mod_resistivity'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    UseResistiveFlux = DoResistiveFlux

    Si2NoEta = Si2No_V(UnitX_)**2/Si2No_V(UnitT_)

    Eta0       = Eta0Si       * Si2NoEta
    Eta0Anom   = Eta0AnomSi   * Si2NoEta
    EtaMaxAnom = EtaMaxAnomSi * Si2NoEta
    JcritInv   = 1.0/(JcritAnomSi  * Si2No_V(UnitJ_))

    ! Spitzer resistivity coefficient with Coulomb logarithm in SI units
    EtaPerpSpitzerSi = 0.5*sqrt(cElectronMass)    * &
         (cLightSpeed*cElectronCharge)**2/ &
         (3*cEps*(cTwoPi*cBoltzmann)**1.5) &
         *CoulombLogarithm

    if(.not.allocated(Eta_GB)) &
         allocate(Eta_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    if(TypeResistivity == 'constant')then
       if(DoTest) write(*,*)NameSub, ': setting Eta_GB=Eta0=',Eta0
       Eta_GB = Eta0
    else
       if(DoTest) write(*,*)NameSub, ': setting Eta_GB=0.0 TypeResistivity=', &
            TypeResistivity
       Eta_GB = 0.0
    end if

    ! The following will ensure that the explicit evaluation of the
    ! resistive diffusion is switched off
    if(UseSemiResistivity) UseResistiveFlux = .false.

    ! Get signed indexes for resistive region(s)
    call get_region_indexes(StringResistRegion, iRegionResist_I)

    if(DoTest)then
       write(*,*)NameSub, ': DoResistiveFlux  = ', DoResistiveFlux
       write(*,*)NameSub, ': UseResistiveFlux = ', UseResistiveFlux
       write(*,*)NameSub, ': UseJouleHeating  = ', UseJouleHeating
       write(*,*)NameSub, ': UseHeatExchange  = ', UseHeatExchange
       write(*,*)NameSub, ': Si2NoEta = ',Si2NoEta
       write(*,*)NameSub, ': Si2NoJ   = ',Si2No_V(UnitJ_)
       write(*,*)NameSub, ': Eta0, Eta0Anom, EtaMaxAnom=', &
            Eta0, Eta0Anom, EtaMaxAnom
       write(*,*)NameSub, ': JcritInv = ', JcritInv
    end if

    call test_stop(NameSub, DoTest)
  end subroutine init_mod_resistivity
  !============================================================================

  subroutine set_resistivity

    use ModMain,    ONLY: nBlock, Unused_B
    use BATL_lib,   ONLY: message_pass_cell, block_inside_regions

    real, allocatable:: ResistFactor_G(:,:,:)

    integer :: iBlock

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'set_resistivity'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(allocated(iRegionResist_I)) &
         allocate(ResistFactor_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK))

    do iBlock = 1, nBlock
       if (Unused_B(iBlock)) CYCLE
       select case(TypeResistivity)
       case('constant')
          Eta_GB(:,:,:,iBlock) = Eta0
       case('spitzer')
          call spitzer_resistivity(iBlock, Eta_GB(:,:,:,iBlock))
       case('anomalous')
          call anomalous_resistivity(iBlock, Eta_GB(:,:,:,iBlock))
       case('raeder')
          call raeder_resistivity(iBlock, Eta_GB(:,:,:,iBlock))
       case('user')
          call user_set_resistivity(iBlock, Eta_GB(:,:,:,iBlock))
       case default
          call stop_mpi(NameSub//': invalid TypeResistivity='//TypeResistivity)
       end select

       ! Modify resistivity of resistive regions are defined
       if(allocated(iRegionResist_I))then
          call block_inside_regions(iRegionResist_I, iBlock, &
               size(ResistFactor_G), 'ghosts', Value_I=ResistFactor_G)
          Eta_GB(:,:,:,iBlock) = Eta_GB(:,:,:,iBlock)*ResistFactor_G
       end if
    end do

    if(allocated(iRegionResist_I)) deallocate(ResistFactor_G)

    if(DoMessagePassResistivity) &
         call message_pass_cell(Eta_GB, nWidthIn=1)

    call test_stop(NameSub, DoTest)
  end subroutine set_resistivity
  !============================================================================

  subroutine spitzer_resistivity(iBlock, Eta_G)

    use ModConst,   ONLY: cProtonMass, cElectronCharge
    use ModPhysics, ONLY: No2Si_V, UnitB_, UnitRho_, UnitTemperature_
    use ModAdvance, ONLY: State_VGB, Rho_, P_
    use ModB0,      ONLY: B0_DGB
    use ModMain,    ONLY: UseB0

    ! Compute Spitzer-type, classical resistivity

    integer, intent(in) :: iBlock
    real,    intent(out):: Eta_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK)

    real :: EtaSi, Coef, B0_DG(3,MinI:MaxI, MinJ:MaxJ, MinK:MaxK)
    integer :: i, j, k
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'spitzer_resistivity'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    Coef =((cProtonMass/cElectronCharge)*No2Si_V(UnitB_)/No2Si_V(UnitRho_))**2
    if(UseB0)then
       B0_DG = B0_DGB(:,:,:,:,iBlock)
    else
       B0_DG = 0.0
    end if
    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI

       ! EtaSi = EtaPerpSpitzer/Te^1.5 in SI units
       EtaSi = EtaPerpSpitzerSi &
            / ( No2Si_V(UnitTemperature_)* &
            State_VGB(P_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock))**1.5

       ! Normalize Eta_G
       Eta_G(i,j,k) = EtaSi*Si2NoEta
    end do; end do; end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine spitzer_resistivity
  !============================================================================

  subroutine anomalous_resistivity(iBlock, Eta_G)

    ! Compute current dependent anomalous resistivity

    use ModCurrent, ONLY: get_current

    integer, intent(in) :: iBlock
    real,    intent(out):: Eta_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK)

    real :: Current_D(3), AbsJ
    integer :: i, j, k

    ! Compute the magnitude of the current density |J|
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'anomalous_resistivity'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1;

       call get_current(i,j,k,iBlock,Current_D)
       AbsJ = sqrt(sum(current_D**2))

       ! Compute the anomalous resistivity::
       ! Eta = Eta0 + Eta0Anom*(|J|/Jcrit-1)

       Eta_G(i,j,k) = Eta0 + &
            min(EtaMaxAnom, max(0.0, Eta0Anom*(AbsJ*JcritInv - 1.0)))

    end do; end do; end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine anomalous_resistivity
  !============================================================================

  subroutine raeder_resistivity(iBlock, Eta_G)

    ! Compute resistivity based on Raeder's formula
    ! Eta = Alpha * j'^2    j' = |j|*gridspacing/(|B| + epsilon)
    ! Eta = 0 if j' < JoverBCrit

    use ModCurrent,  ONLY: get_current
    use ModAdvance,  ONLY: State_VGB
    use ModB0,       ONLY: B0_DGB
    use ModVarIndexes, ONLY: Bx_, Bz_
    use ModMain, ONLY: UseB0

    integer, intent(in) :: iBlock
    real,    intent(out):: Eta_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    real :: Current_D(3), AbsJ, AbsB, JoverB, b_D(3)
    integer :: i, j, k, l, m, n

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'raeder_resistivity'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    if(nDim /= 3) call stop_mpi(NameSub//' is implemented for 3D only')

    ! Compute |J| and |B|
    do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
       call get_current(i,j,k,iBlock,Current_D)
       AbsJ = sqrt(sum(Current_D**2))

       ! Calculate AbsB from the average of 26 neighboring cells and itself
       AbsB = 0.0
       do n=-1,1; do m=-1,1; do l=-1,1
          b_D = State_VGB(Bx_:Bz_,i+l,j+m,k+n,iBlock)
          if(UseB0) b_D = b_D + B0_DGB(:,i+l,j+m,k+n,iBlock)
          AbsB = AbsB + sqrt(sum(b_D**2))
       end do; end do; end do
       AbsB = AbsB/27.0

       ! No dx here, because that causes problems at resolution changes
       JoverB = min(AbsJ/(AbsB + 1e-8), 1.0)

       ! Compute Raeder resistivity
       if(JoverB >= JoverBCrit)then
          Eta_G(i,j,k) = EtaCoeff*(JoverB**2)
       else
          Eta_G(i,j,k) = 0.0
       end if
    end do; end do; end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine raeder_resistivity
  !============================================================================

  subroutine calc_resistivity_source(iBlock)

    use ModMain,       ONLY: x_
    use ModGeometry,   ONLY: true_cell
    use BATL_lib,      ONLY: IsRzGeometry, Xyz_DGB
    use ModCurrent,    ONLY: get_current
    use ModPhysics,    ONLY: GammaElectronMinus1, InvGammaElectronMinus1
    use ModVarIndexes, ONLY: p_, Pe_, Ppar_, Bz_, Energy_
    use ModAdvance,    ONLY: Source_VC, &
         UseElectronPressure, UseAnisoPressure, UseAnisoPe

    integer, intent(in):: iBlock

    ! Variables needed for Joule heating
    real :: Current_D(3), JouleHeating

    integer:: i, j, k

    logical :: DoTestCell = .false.

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_resistivity_source'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    JouleHeating = 0.0

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE

       DoTestCell = DoTest .and. &
            i == iTest .and. j == jTest .and. k == kTest

       if (DoTestCell) then
          write(*,*) NameSub, ': initial JouleHeating    =', &
               JouleHeating
          write(*,*) NameSub, ': initial P source        =', &
               Source_VC(P_,i,j,k)
          if (UseElectronPressure) &
               write(*,*) NameSub, ': initial Pe source       =', &
               Source_VC(Pe_,i,j,k)
          if (UseAnisoPressure) &
               write(*,*) NameSub, ': initial Ppar Source     =',&
               Source_VC(Ppar_,i,j,k)
          write(*,*) NameSub, ': initial energy source   =', &
               Source_VC(Energy_,i,j,k)
          write(*,*) NameSub, ': initial Bz source       =', &
               Source_VC(Bz_,i,j,k)
       end if

       if(UseAnisoPe) call stop_mpi(NameSub// &
            ' for UseAnisoPe has not been implemented yet.')

       if(UseJouleHeating)then
          call get_current(i,j,k,iBlock,Current_D)
          JouleHeating = GammaElectronMinus1 &
               * Eta_GB(i,j,k,iBlock) * sum(Current_D**2)
       end if
       if(UseElectronPressure) then
          Source_VC(Pe_,i,j,k) = Source_VC(Pe_,i,j,k) + JouleHeating
          ! Remove Joule heating from ion energy equation
          if(UseResistiveFlux) &
               Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
               - InvGammaElectronMinus1*JouleHeating
       else
          Source_VC(P_,i,j,k) = Source_VC(P_,i,j,k) + JouleHeating

          ! the same amount of Joule heating applies on Ppar
          if(UseAnisoPressure) &
               Source_VC(Ppar_,i,j,k)  = Source_VC(Ppar_,i,j,k) + JouleHeating

          if(.not.UseResistiveFlux) &
               Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
               + InvGammaElectronMinus1*JouleHeating
       end if

       ! rz-geometrical source terms
       if(IsRzGeometry .and. UseResistiveFlux)then
          ! calculate cell centered current if it is not yet done
          if(.not.UseJouleHeating) call get_current(i,j,k,iBlock,Current_D)

          ! Source[Bphi] = -eta*Jz / radius
          Source_VC(Bz_,i,j,k) = Source_VC(Bz_,i,j,k) &
               - Eta_GB(i,j,k,iBlock)*Current_D(x_)/Xyz_DGB(2,i,j,k,iBlock)
       end if

       if (DoTestCell) then
          write(*,*) NameSub, ': corrected JouleHeating  =', &
               JouleHeating
          write(*,*) NameSub, ': corrected P source      =', &
               Source_VC(P_,i,j,k)
          if (UseElectronPressure) &
               write(*,*) NameSub, ': corrected Pe source     =', &
               Source_VC(Pe_,i,j,k)
          if (UseAnisoPressure) &
               write(*,*) NameSub, ': corrected Ppar source     =', &
               Source_VC(Ppar_,i,j,k)
          write(*,*) NameSub, ': corrected energy source =', &
               Source_VC(Energy_,i,j,k)
          write(*,*) NameSub, ': initial Bz source       =', &
               Source_VC(Bz_,i,j,k)
       end if
    end do; end do; end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine calc_resistivity_source
  !============================================================================

  subroutine calc_heat_exchange

    use ModMain,       ONLY: Cfl, nBlock, Unused_B
    use ModGeometry,   ONLY: true_cell
    use ModPhysics,    ONLY: GammaMinus1, GammaElectronMinus1, IonMassPerCharge
    use ModVarIndexes, ONLY: Rho_, p_, Pe_, Ppar_
    use ModAdvance,    ONLY: time_blk, State_VGB, UseAnisoPressure, UseAnisoPe
    use ModEnergy,     ONLY: calc_energy_cell

    real :: DtLocal
    real :: HeatExchange, HeatExchangePeP, HeatExchangePePpar
    integer:: i, j, k, iBlock
    logical :: DoTestCell
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_heat_exchange'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    HeatExchange = 0.0
    HeatExchangePeP = 0.0
    HeatExchangePePpar = 0.0

    call set_resistivity

    do iBlock = 1, nBlock
       if (Unused_B(iBlock)) CYCLE

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          DoTestCell = DoTest .and. iBlock == iBlockTest .and. &
               i == iTest .and. j == jTest .and. k == kTest

          DtLocal = Cfl*time_BLK(i,j,k,iBlock)

          ! For single ion fluid the ion-electron collision results in a
          ! heat exchange term for the electron pressure
          ! See eq. 4.124c in Schunk and Nagy.

          ! Explicit heat exchange (energy)
          HeatExchange = Eta_GB(i,j,k,iBlock) * &
               3*State_VGB(Rho_,i,j,k,iBlock)*(1./IonMassPerCharge**2)

          if (DoTestCell) &
               write(*,*) NameSub, ': explicit HeatExchange    =', HeatExchange

          ! Point-implicit correction for stability: H' = H/(1+4./3*dt*H)
          HeatExchange = HeatExchange / (1 + 4.0/3.0*DtLocal*HeatExchange)

          HeatExchangePeP = HeatExchange &
               *(State_VGB(P_,i,j,k,iBlock) - State_VGB(Pe_,i,j,k,iBlock))

          if (DoTestCell) then
             write(*,*) NameSub, ': corrected HeatExchange  =', HeatExchange
             write(*,*) NameSub, ': heatExchangePeP         =', HeatExchangePeP
             write(*,*) NameSub, ': initial State_VGB(P_)   =', &
                  State_VGB(P_,i,j,k,iBlock)
             if (UseAnisoPressure) &
                  write(*,*) NameSub, ': initial State_VGB(Ppar_)=', &
                  State_VGB(Ppar_,i,j,k,iBlock)
             write(*,*) NameSub, ': initial State_VGB(Pe_)  =', &
                  State_VGB(Pe_,i,j,k,iBlock)
          end if

          if(UseAnisoPe) call stop_mpi(NameSub// &
               ' for UseAnisoPe has not been implemented yet.')

          ! Heat exchange for parallel ion pressure
          if(UseAnisoPressure)then
             HeatExchangePePpar = HeatExchange &
                  *(State_VGB(Ppar_,i,j,k,iBlock) &
                  - State_VGB(Pe_,i,j,k,iBlock))

             State_VGB(Ppar_,i,j,k,iBlock) = State_VGB(Ppar_,i,j,k,iBlock) &
                  - DtLocal*GammaMinus1*HeatExchangePePpar
          end if

          ! Heat exchange for the ions
          State_VGB(P_,i,j,k,iBlock) = State_VGB(P_,i,j,k,iBlock) &
               - DtLocal*GammaMinus1*HeatExchangePeP

          ! Heat exchange for the electrons
          State_VGB(Pe_,i,j,k,iBlock) = State_VGB(Pe_,i,j,k,iBlock) &
               + DtLocal*GammaElectronMinus1*HeatExchangePeP

          if (DoTestCell) then
             write(*,*) NameSub, ': corrected State_VGB(P_)  =', &
                  State_VGB(P_,i,j,k,iBlock)
             if (UseAnisoPressure) then
                write(*,*) NameSub, ': heatExchangePePpAr       =', &
                     HeatExchangePePpar
                write(*,*) NameSub, ': corrected State_VGB(Ppar_) =', &
                     State_VGB(Ppar_,i,j,k,iBlock)
             end if
             write(*,*) NameSub, ': corrected State_VGB(Pe_) =', &
                  State_VGB(Pe_,i,j,k,iBlock)
          end if

       end do; end do; end do

       call calc_energy_cell(iBlock)
    end do

    call test_stop(NameSub, DoTest)
  end subroutine calc_heat_exchange
  !============================================================================

  ! Interface for (Semi-)implicit collisional and Hall resistivity

  subroutine get_impl_resistivity_state(SemiAll_VCB)

    use ModAdvance,    ONLY: State_VGB
    use ModImplicit,   ONLY: nVarSemiAll, nBlockSemi, iBlockFromSemi_B
    use ModVarIndexes, ONLY: Bx_, Bz_
    use BATL_lib, ONLY: nBlock, Unused_B, iProc

    real, intent(out):: SemiAll_VCB(nVarSemiAll,nI,nJ,nK,nBlock)

    ! Initialize variables for RHS and Jacobian calculation
    ! and find the blocks that require semi-implicit solver

    integer:: i, j, k, iBlock

    logical:: IsSemiBlock

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'get_impl_resistivity_state'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(UseSemiResistivity) call init_impl_resistivity
    if(UseSemiHallResist)  call init_impl_hall_resist

    ! Collect the semi-implicit blocks
    nBlockSemi = 0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       ! Check if the block needs to be advanced with semi-implicit scheme
       IsSemiBlock = .false.
       if(UseSemiResistivity) &
            IsSemiBlock = any(Eta_DFDB(:,:,:,:,:,iBlock) /= 0.0)
       if(UseSemiHallResist .and. .not.IsSemiBlock) &
            IsSemiBlock = any(Bne_DFDB(:,:,:,:,:,iBlock) /= 0.0)
       if(.not.IsSemiBlock) CYCLE

       ! Add the semi-implicit block
       nBlockSemi = nBlockSemi + 1

       ! Store index mapping
       iBlockFromSemi_B(nBlockSemi) = iBlock

       ! Store the magnetic field in SemiAll_VCB
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          SemiAll_VCB(BxImpl_:BzImpl_,i,j,k,nBlockSemi) = &
               State_VGB(Bx_:Bz_,i,j,k,iBlock)
       end do; end do; end do
    end do

    if(DoTest) write(*,*) NameSub,' iProc, nBlockSemi=', iProc, nBlockSemi

    call test_stop(NameSub, DoTest)
  end subroutine get_impl_resistivity_state
  !============================================================================

  subroutine init_impl_resistivity

    use BATL_lib,        ONLY: IsCartesian, message_pass_cell, &
         CellFace_DB, FaceNormal_DDFB, nBlock, Unused_B
    use ModNumConst,     ONLY: i_DD

    integer:: iDim, i, j, k, Di, Dj, Dk, iBlock
    real:: Eta
    real:: FaceNormal_D(nDim)

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'init_impl_resistivity'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(.not.UseResistivity) RETURN

    ! Calculate Eta*FaceNormal_D to be used for resistive flux
    if(.not.allocated(Eta_DFDB)) &
         allocate(Eta_DFDB(nDim,nI+1,nJ+1,nK+1,nDim,MaxBlock))

    ! Message pass resistivity to fill in ghost cells
    call message_pass_cell(Eta_GB, DoSendCornerIn=.false.)

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       Eta_DFDB(:,:,:,:,:,iBlock) = 0.0

       do iDim = 1, nDim
          Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
          if(IsCartesian)then
             FaceNormal_D = 0; FaceNormal_D(iDim) = CellFace_DB(iDim,iBlock)
          end if
          do k = 1, nK+Dk; do j = 1, nJ+Dj; do i = 1, nI+Di
             if(.not.IsCartesian) &
                  FaceNormal_D = FaceNormal_DDFB(:,iDim,i,j,k,iBlock)

             Eta= 0.5*(Eta_GB(i,j,k,iBlock) + Eta_GB(i-Di,j-Dj,k-Dk,iBlock))

             Eta_DFDB(:,i,j,k,iDim,iBlock) = Eta*FaceNormal_D
          end do; end do; end do
       end do
    end do

    call test_stop(NameSub, DoTest)
  end subroutine init_impl_resistivity
  !============================================================================
  subroutine init_impl_hall_resist

    use BATL_lib,        ONLY: nBlock, Unused_B, IsCartesian, &
         Xyz_DGB, CellSize_DB, CellFace_DB, FaceNormal_DDFB
    use ModAdvance,      ONLY: State_VGB
    use ModMain,         ONLY: UseB0
    use ModNumConst,     ONLY: i_DD, cPi
    use ModVarIndexes,   ONLY: Bx_, Bz_
    use ModMultiFluid,   ONLY: iRhoIon_I
    use ModB0,           ONLY: B0_DGB
    use ModHallResist,   ONLY: UseHallResist, &
         set_hall_factor_face, HallFactor_DF, IsHallBlock, HallCmaxFactor

    integer:: iDim, i, j, k, Di, Dj, Dk, i1, j1, k1, iBlock

    real:: HallCoeff, Rho, b_D(MaxDim)

    real:: InvDxyz

    real:: FaceNormal_D(nDim)

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'init_impl_hall_resist'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(.not.UseHallResist) RETURN

    ! Calculate B/ne = m/e * B/rho on the face to be used for the Hall flux
    if(.not.allocated(Bne_DFDB)) &
         allocate(Bne_DFDB(MaxDim,nI+1,nJ+1,nK+1,nDim,MaxBlock))

    ! Whistler diffusion coefficient
    if(HallCmaxFactor > 0 .and. .not.allocated(WhistlerCoeff_FDB)) &
         allocate(WhistlerCoeff_FDB(nI+1,nJ+1,nK+1,nDim,MaxBlock))

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       Bne_DFDB(:,:,:,:,:,iBlock) = 0.0
       if(HallCmaxFactor > 0) WhistlerCoeff_FDB(:,:,:,:,iBlock) = 0.0

       call set_hall_factor_face(iBlock)
       if(.not.IsHallBlock) CYCLE

       do iDim = 1, nDim
          Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
          if(HallCmaxFactor > 0 .and.  IsCartesian) then
             InvDxyz = 1./CellSize_DB(iDim,iBlock)
             FaceNormal_D = 0.0; FaceNormal_D(iDim) = CellFace_DB(iDim,iBlock)
          end if
          do k = 1, nK+Dk; do j = 1, nJ+Dj; do i = 1, nI+Di

             ! Check if the Hall coefficient is positive for this face
             HallCoeff = HallFactor_DF(iDim,i,j,k)
             if(HallCoeff <= 0.0) CYCLE

             ! Cell center indexes for the cell on the left of the face
             i1 = i-Di; j1 = j-Dj; k1 = k - Dk

             ! Average rho and B to the face
             ! NOTE: HallCoeff is normalized to total ion mass density
             Rho = 0.5*(sum(State_VGB(iRhoIon_I,i1,j1,k1,iBlock)) &
                  +     sum(State_VGB(iRhoIon_I,i, j, k, iBlock)) )

             b_D = 0.5*(State_VGB(Bx_:Bz_,i1,j1,k1,iBlock) &
                  +     State_VGB(Bx_:Bz_,i,j,k,iBlock)    )

             ! B0 on the face is actually average of cell center B0
             ! It is easier to use B0_DGB then then the 3 face arrays
             ! The res change on the coarse side will get overwritten
             if(UseB0) b_D = b_D + &
                  0.5*(B0_DGB(:,i1,j1,k1,iBlock) + B0_DGB(:,i,j,k,iBlock))

             ! Calculate B/ne for the face
             Bne_DFDB(:,i,j,k,iDim,iBlock) = HallCoeff*b_D/Rho

             if(HallCmaxFactor > 0) then
                if(.not. IsCartesian) then
                   InvDxyz = 1.0/sqrt( sum( &
                        ( Xyz_DGB(:, i, j, k, iBlock)&
                        - Xyz_DGB(:, i1, j1, k1, iBlock))**2) )
                   FaceNormal_D = FaceNormal_DDFB(:,iDim,i,j,k,iBlock)
                endif
                ! Diffuse the linear stage by
                ! 0.5*HallCmaxFactor*c_whistler (like first order TVD).
                ! The whistler speed in the face normal direction is
                ! c_whistler = |B_n| * pi/(e*n_e*dx)
                ! Store the D coefficient for the diffusive flux
                ! F_i+1/2 = D_i+1/2 * (B_i+1 - B_i) where
                ! D = HallCmaxFactor*0.5*pi*|Area*B_n|/(e*n_e*dx)
                WhistlerCoeff_FDB(i,j,k,iDim,iBlock) = HallCmaxFactor*0.5*cPi*&
                     abs(sum(FaceNormal_D*Bne_DFDB(1:nDim,i,j,k,iDim,iBlock)))&
                     *InvDxyz
             end if
          end do; end do; end do
       end do
    end do

    call test_stop(NameSub, DoTest)
  end subroutine init_impl_hall_resist
  !============================================================================

  subroutine get_resistivity_rhs(iBlock, StateImpl_VG, Rhs_VC, IsLinear)

    use BATL_lib,        ONLY: store_face_flux, IsCartesian, IsRzGeometry, &
         Xyz_DGB, CellSize_DB, CellVolume_GB, CellFace_DB, FaceNormal_DDFB, &
         nG, i0_, j0_, k0_, nIp1_, nJp1_, nKp1_
    use ModFaceGradient, ONLY: get_face_curl
    use ModImplicit,     ONLY: nVarSemi, &
         FluxImpl_VXB, FluxImpl_VYB, FluxImpl_VZB
    use ModNumConst,     ONLY: i_DD
    use ModSize,         ONLY: x_, y_, z_
    use ModGeometry,     ONLY: true_cell, true_BLK
    use ModHallResist,   ONLY: HallCmaxFactor
    use ModCellGradient, ONLY: calc_cell_curl_ghost
    use ModCoordTransform, ONLY: determinant

    integer, intent(in) :: iBlock
    real, intent(inout) :: StateImpl_VG(nVarSemi,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(out)   :: Rhs_VC(nVarSemi,nI,nJ,nK)
    logical, intent(in) :: IsLinear

    ! resistive flux for operator split scheme
    real, allocatable, save :: FluxImpl_VFD(:,:,:,:,:)

    real :: DetJ, InvDxHalf, InvDyHalf, InvDzHalf
    real, allocatable :: Egen_DG(:,:,:,:)
    real :: dBgen_D(3)
    real :: InvJac_DD(3,3)
    
    integer :: iDim, i, j, k, Di, Dj, Dk
    real    :: Current_D(MaxDim), Jx, InvDy2, FaceNormal_D(nDim)
    real    :: Jnormal, BneNormal
    logical :: IsNewBlock
    logical:: DoTest, DoTestCell
    character(len=*), parameter:: NameSub = 'get_resistivity_rhs'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    if(DoTest)write(*,*)NameSub,' true cell=',true_cell(iTest,jTest,ktest,iBlockTest)

    IsNewBlock = .true.

    Rhs_VC = 0.0

    ! Should we keep this or not?
    !if(.not. true_BLK(iblock)) RETURN

    if(.not.allocated(FluxImpl_VFD)) allocate( &
         FluxImpl_VFD(nVarSemi,nI+1,nJ+1,nK+1,nDim))

    if(.not.UseCentralDifference)then
       ! Loop over face directions
       do iDim = 1, nDim
          Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
          
          if(UseSemiHallResist .and. IsCartesian) then
             FaceNormal_D = 0.0; FaceNormal_D(iDim) = CellFace_DB(iDim,iBlock)
          end if
          
          ! Loop over cell faces orthogonal to iDim
          do k = 1, nK+Dk; do j = 1, nJ+Dj; do i = 1, nI+Di
             if(  .not.true_cell(i,j,k,iBlock) .and. &
                  .not.true_cell(i-Di,j-Dj,k-Dk,iBlock)) CYCLE
             
             call get_face_curl(iDim, i, j, k, iBlock, IsNewBlock, &
                  StateImpl_VG, Current_D)
             
             FluxImpl_VFD(BxImpl_:BzImpl_,i,j,k,iDim) = 0.0
             
             if(UseSemiResistivity)then
                ! Resistive flux
                ! dB/dt = -curl E
                !       = -sum(FaceNormal x eta*J)
                !       = -sum(eta*FaceNormal x J)
                !       = -div(FluxImpl)
                if(nDim == 2) FluxImpl_VFD(BxImpl_,i,j,k,iDim) = &
                     + Eta_DFDB(y_,i,j,k,iDim,iBlock)*Current_D(z_)
                if(nDim == 3) FluxImpl_VFD(BxImpl_,i,j,k,iDim) = &
                     + Eta_DFDB(y_,i,j,k,iDim,iBlock)*Current_D(z_) &
                     - Eta_DFDB(z_,i,j,k,iDim,iBlock)*Current_D(y_)
                
                if(nDim < 3) FluxImpl_VFD(ByImpl_,i,j,k,iDim) = &
                     - Eta_DFDB(x_,i,j,k,iDim,iBlock)*Current_D(z_)
                if(nDim == 3)FluxImpl_VFD(ByImpl_,i,j,k,iDim) = &
                     + Eta_DFDB(z_,i,j,k,iDim,iBlock)*Current_D(x_) &
                     - Eta_DFDB(x_,i,j,k,iDim,iBlock)*Current_D(z_)
                
                if(nDim == 1) FluxImpl_VFD(BzImpl_,i,j,k,iDim) = &
                     + Eta_DFDB(x_,i,j,k,iDim,iBlock)*Current_D(y_)
                if(nDim > 1) FluxImpl_VFD(BzImpl_,i,j,k,iDim) = &
                     + Eta_DFDB(x_,i,j,k,iDim,iBlock)*Current_D(y_) &
                     - Eta_DFDB(y_,i,j,k,iDim,iBlock)*Current_D(x_)   
             end if
             
             if(UseSemiHallResist)then
                ! Hall MHD flux
                ! dB/dt = -curl E
                !       = -div(uH B - B uH)  where uH = -J/(n e)
                !       = -div(-J B/ne + B/ne J)
                !       = -sum(-A.J B/ne + A.B/ne J)
                !
                ! Note that B is frozen in, only J varies with StateImpl
                
                if(.not.IsCartesian) &
                  FaceNormal_D = FaceNormal_DDFB(:,iDim,i,j,k,iBlock)
                
                ! Normal component of current and B/ne vectors
                Jnormal   = sum(FaceNormal_D*Current_D(1:nDim))
                BneNormal = sum(FaceNormal_D*&
                     Bne_DFDB(1:nDim,i,j,k,iDim,iBlock))
                
                ! Flux = Bn/ne J - Jn B/ne
                FluxImpl_VFD(BxImpl_:BzImpl_,i,j,k,iDim) = &
                     FluxImpl_VFD(BxImpl_:BzImpl_,i,j,k,iDim) &
                     - Jnormal*Bne_DFDB(:,i,j,k,iDim,iBlock) &
                     + BneNormal*Current_D
                
                ! Add whistler diffusion in the linear phase.
                ! The diffusive flux is F_whistler=-WhistlerCoeff*(B_i+1 - B_i)
                if(HallCmaxFactor > 0 .and. IsLinear) &
                     FluxImpl_VFD(BxImpl_:BzImpl_,i,j,k,iDim) = &
                     FluxImpl_VFD(BxImpl_:BzImpl_,i,j,k,iDim) + &
                     WhistlerCoeff_FDB(i,j,k,iDim,iBlock)* &
                     (StateImpl_VG(BxImpl_:BzImpl_,i-Di,j-Dj,k-Dk) &
                     -StateImpl_VG(BxImpl_:BzImpl_,i,j,k) )
             end if
          end do; end do; end do
       end do

       ! Store the fluxes at resolution changes for restoring conservation
       call store_face_flux(iBlock, nVarSemi, FluxImpl_VFD, &
            FluxImpl_VXB, FluxImpl_VYB, FluxImpl_VZB)
       
       do iDim = 1, nDim
          Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             
             ! RHS_i += -(Flux_i+1/2 - Flux_i-1/2)                                                                                                          
             Rhs_VC(:,i,j,k) = Rhs_VC(:,i,j,k) &
                  -(FluxImpl_VFD(:,i+Di,j+Dj,k+Dk,iDim) &
                  - FluxImpl_VFD(:,i,j,k,iDim))/CellVolume_GB(i,j,k,iBlock)
          end do; end do; end do
       end do
          
    else ! Use central difference scheme for resistivity term

       ! Hall term
       if(UseSemiHallResist)then
          ! Loop over face directions                                     
          do iDim = 1, nDim
             Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
             
             if(UseSemiHallResist .and. IsCartesian) then
                FaceNormal_D = 0.0 
                FaceNormal_D(iDim) = CellFace_DB(iDim,iBlock)
             end if

             ! Loop over cell faces orthogonal to iDim  
             do k = 1, nK+Dk; do j = 1, nJ+Dj; do i = 1, nI+Di
                if(  .not.true_cell(i,j,k,iBlock) .and. &
                     .not.true_cell(i-Di,j-Dj,k-Dk,iBlock)) CYCLE
                
                call get_face_curl(iDim, i, j, k, iBlock, IsNewBlock, &
                     StateImpl_VG, Current_D)
                
                FluxImpl_VFD(BxImpl_:BzImpl_,i,j,k,iDim) = 0.0
                
                ! Hall MHD flux                                           
                ! dB/dt = -curl E                                    
                !       = -div(uH B - B uH)  where uH = -J/(n e)     
                !       = -div(-J B/ne + B/ne J)                    
                !       = -sum(-A.J B/ne + A.B/ne J)                  
                !                                                
                ! Note that B is frozen in, only J varies with StateImpl
                if(.not.IsCartesian) &
                     FaceNormal_D = FaceNormal_DDFB(:,iDim,i,j,k,iBlock)
                
                ! Normal component of current and B/ne vectors        
                Jnormal   = sum(FaceNormal_D*Current_D(1:nDim))
                BneNormal = sum(FaceNormal_D*&
                     Bne_DFDB(1:nDim,i,j,k,iDim,iBlock))
                
                ! Flux = Bn/ne J - Jn B/ne                            
                FluxImpl_VFD(BxImpl_:BzImpl_,i,j,k,iDim) = &
                     FluxImpl_VFD(BxImpl_:BzImpl_,i,j,k,iDim) &
                     - Jnormal*Bne_DFDB(:,i,j,k,iDim,iBlock) &
                     + BneNormal*Current_D
                
                ! Add whistler diffusion in the linear phase.                
                ! The diffusive flux is F_whistler=-WhistlerCoeff*(B_i+1 - B_i)
                if(HallCmaxFactor > 0 .and. IsLinear) &
                     FluxImpl_VFD(BxImpl_:BzImpl_,i,j,k,iDim) = &
                     FluxImpl_VFD(BxImpl_:BzImpl_,i,j,k,iDim) + &
                     WhistlerCoeff_FDB(i,j,k,iDim,iBlock)* &
                     (StateImpl_VG(BxImpl_:BzImpl_,i-Di,j-Dj,k-Dk) &
                     -StateImpl_VG(BxImpl_:BzImpl_,i,j,k) )
             end do; end do; end do
          end do
       end if

       ! Resistivity term
       ! Central difference scheme calculation of Rhs for contraining divB
       if(.not.allocated(Egen_DG)) &
            allocate(Egen_DG(3,MinI:MaxI,MinJ:MaxJ,MinK:MaxK))
       
       InvDxHalf = 0.5/CellSize_DB(1,iBlock)
       InvDyHalf = 0.5/CellSize_DB(2,iBlock)
       InvDzHalf = 0.5/CellSize_DB(3,iBlock)

       ! Calc E
       if(IsCartesian)then ! Cartesian coord.
          do k=k0_,nKp1_; do j=j0_,nJp1_; do i=i0_,nIp1_
             if(.not.true_cell(i,j,k,iBlock)) CYCLE

             ! Get current from curl B
             call calc_cell_curl_ghost(i,j,k,iBlock,StateImpl_VG,Current_D)

             Egen_DG(:,i,j,k) = Eta_GB(i,j,k,iBlock)*Current_D
          end do; end do; end do
       else ! Generalized coord.
          do k=k0_,nKp1_; do j=j0_,nJp1_; do i=i0_,nIp1_
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             ! Calculate inv. jacobian matrix from Cartesian to general Coord.
             InvJac_DD = get_InvJacobian(i,j,k,iBlock)

             ! Get current from curl B
             call calc_cell_curl_ghost(i,j,k,iBlock,StateImpl_VG,Current_D)

             ! E in generalized coord.
             Egen_DG(:,i,j,k) = matmul(transpose(InvJac_DD),&
                  Eta_GB(i,j,k,iBlock)*Current_D)
          end do; end do; end do
       end if

       ! Central difference update of RHS of induction equation
       if(IsCartesian)then
          do k=1,nK; do j=1,nJ; do i=1,nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             
             Rhs_VC(BxImpl_,i,j,k) = &
                  (Egen_DG(y_,i,j,k+1) - Egen_DG(y_,i,j,k-1)) * InvDzHalf - &
                  (Egen_DG(z_,i,j+1,k) - Egen_DG(z_,i,j-1,k)) * InvDyHalf
             
             Rhs_VC(ByImpl_,i,j,k) = &
                  (Egen_DG(z_,i+1,j,k) - Egen_DG(z_,i-1,j,k)) * InvDxHalf - &
                  (Egen_DG(x_,i,j,k+1) - Egen_DG(x_,i,j,k-1)) * InvDzHalf
             
             Rhs_VC(BzImpl_,i,j,k) = &
                  (Egen_DG(x_,i,j+1,k) - Egen_DG(x_,i,j-1,k)) * InvDyHalf - &
                  (Egen_DG(y_,i+1,j,k) - Egen_DG(y_,i-1,j,k)) * InvDxHalf
          end do; end do; end do
       else
          do k=1,nK; do j=1,nJ; do i=1,nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             
             InvJac_DD = get_InvJacobian(i,j,k,iBlock)

             ! Calculate determinant of Jacobian = 1/det(InvJac)
             DetJ = 1 / determinant(InvJac_DD)
             
             ! Calculate curl Egen
             dBgen_D(x_) = &
                  (Egen_DG(y_,i,j,k+1) - Egen_DG(y_,i,j,k-1)) * InvDzHalf - &
                  (Egen_DG(z_,i,j+1,k) - Egen_DG(z_,i,j-1,k)) * InvDyHalf
             
             dBgen_D(y_) = &
                  (Egen_DG(z_,i+1,j,k) - Egen_DG(z_,i-1,j,k)) * InvDxHalf - &
                  (Egen_DG(x_,i,j,k+1) - Egen_DG(x_,i,j,k-1)) * InvDzHalf
             
             dBgen_D(z_) = &
                  (Egen_DG(x_,i,j+1,k) - Egen_DG(x_,i,j-1,k)) * InvDyHalf - &
                  (Egen_DG(y_,i+1,j,k) - Egen_DG(y_,i-1,j,k)) * InvDxHalf
             
             ! Convert dBgen into dB
             Rhs_VC(BxImpl_:BzImpl_,i,j,k) = &
                  DetJ * matmul(InvJac_DD, dBgen_D)
          end do; end do; end do         
       end if

    end if
   
    if(IsRzGeometry)then
       InvDy2 = 0.5/CellSize_DB(y_,iBlock)
       
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE
          ! Jx = Dbz/Dy - Dby/Dz (Jz = Dbphi/Dr in rz-geometry)
          ! Note that set_block_field3 has already been called so that we can
          ! use central difference.
          Jx = ( StateImpl_VG(BzImpl_,i,j+1,k) &
               - StateImpl_VG(BzImpl_,i,j-1,k) )*InvDy2
          
          ! Correct current for rz-geometry: Jz = Jz + Bphi/radius
          Jx = Jx + StateImpl_VG(BzImpl_,i,j,k)/Xyz_DGB(y_,i,j,k,iBlock)
          
          ! in rz-geonetry: Rhs[Bphi] = -eta*Jz / radius
          Rhs_VC(BzImpl_,i,j,k) = Rhs_VC(BzImpl_,i,j,k) &
               - Eta_GB(i,j,k,iBlock)*Jx/Xyz_DGB(y_,i,j,k,iBlock)
       end do; end do; end do
       
       ! Anything to do here for Hall term? !!!
       
    end if
    
    call test_stop(NameSub, DoTest, iBlock)

contains
  function get_InvJacobian(i,j,k,iBlock) RESULT(InvJac_DD)
    
    integer, intent(in) :: i,j,k,iBlock
    real :: InvJac_DD(3,3)
    !-----------------------------------------------------

    ! Calculate the inverse Jacobian matrix
    InvJac_DD(:,1) = InvDxHalf *&
         (Xyz_DGB(:,i+1,j,k,iBlock) - Xyz_DGB(:,i-1,j,k,iBlock))
    InvJac_DD(:,2) = InvDyHalf *&
         (Xyz_DGB(:,i,j+1,k,iBlock) - Xyz_DGB(:,i,j-1,k,iBlock))
    InvJac_DD(:,3) = InvDzHalf *&
         (Xyz_DGB(:,i,j,k+1,iBlock) - Xyz_DGB(:,i,j,k-1,iBlock))
  end function get_InvJacobian

  end subroutine get_resistivity_rhs
  !============================================================================

  subroutine add_jacobian_resistivity(iBlock, nVarImpl, Jacobian_VVCI)

    ! Calculate the Jacobian for the preconditioning of
    ! collisional resistivity.

    use BATL_lib,        ONLY: IsCartesianGrid, CellSize_DB, CellVolume_GB
    use ModFaceGradient, ONLY: set_block_jacobian_face, DcoordDxyz_DDFD
    use ModImplicit,     ONLY: UseNoOverlap, nStencil
    use ModNumConst,     ONLY: i_DD
    use ModGeometry,     ONLY: true_cell
    use ModImplicit,     ONLY: nStencil
    use ModAdvance,      ONLY: B_

    integer, intent(in):: iBlock
    integer, intent(in):: nVarImpl
    real, intent(inout):: Jacobian_VVCI(nVarImpl,nVarImpl,nI,nJ,nK,nStencil)

    integer :: iB
    integer :: iDim, iDir, jDir, i, j, k, Di, Dj, Dk
    integer:: iVar, jVar
    real :: DiffLeft, DiffRight, InvDcoord_D(nDim), Coeff

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'add_jacobian_resistivity'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    ! Set the base index value for magnetic field variables
    if(UseSemiResistivity)then
       iB = 0
    else
       iB = B_
    end if

    InvDcoord_D = 1/CellSize_DB(1:nDim,iBlock)

    ! the transverse diffusion is ignored in the Jacobian
    if(IsCartesianGrid)then
       do iDim = 1, nDim
          Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             Coeff = InvDcoord_D(iDim)/CellVolume_GB(i,j,k,iBlock)

             do iDir = 1, MaxDim
                if(iDim == iDir) CYCLE

                ! Variable index in the Jacobian
                iVar = iB + iDir

                DiffLeft  = Coeff*Eta_DFDB(iDim,i,j,k,iDim,iBlock)
                DiffRight = Coeff*Eta_DFDB(iDim,i+Di,j+Dj,k+Dk,iDim,iBlock)

                Jacobian_VVCI(iVar,iVar,i,j,k,1) = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,1) - (DiffLeft + DiffRight)

                if(UseNoOverlap)then
                   if(  iDim==1.and.i==1  .or. &
                        iDim==2.and.j==1  .or. &
                        iDim==3.and.k==1)        DiffLeft = 0.0
                   if(  iDim==1.and.i==nI .or. &
                        iDim==2.and.j==nJ .or. &
                        iDim==3.and.k==nK)       DiffRight = 0.0
                end if

                Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim)   = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim) + DiffLeft
                Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim+1) = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim+1) + DiffRight
             end do
          end do; end do; end do
       end do
    else
       call set_block_jacobian_face(iBlock)

       do iDim = 1, nDim
          Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE
             Coeff = InvDcoord_D(iDim)/CellVolume_GB(i,j,k,iBlock)

             do iDir = 1, MaxDim; do jDir = 1, nDim
                if(iDir ==jDir) CYCLE

                ! Variable index in the Jacobian
                iVar = iB + iDir
                jVar = iB + jDir

                DiffLeft = Eta_DFDB(jDir,i,j,k,iDim,iBlock) &
                     *DcoordDxyz_DDFD(iDim,jDir,i,j,k,iDim)*Coeff
                DiffRight = Eta_DFDB(jDir,i+Di,j+Dj,k+Dk,iDim,iBlock) &
                     *DcoordDxyz_DDFD(iDim,jDir,i+Di,j+Dj,k+Dk,iDim)*Coeff

                Jacobian_VVCI(iVar,iVar,i,j,k,1) = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,1) - (DiffLeft + DiffRight)

                if(UseNoOverlap)then
                   if(  iDim==1.and.i==1  .or. &
                        iDim==2.and.j==1  .or. &
                        iDim==3.and.k==1)        DiffLeft = 0.0
                   if(  iDim==1.and.i==nI .or. &
                        iDim==2.and.j==nJ .or. &
                        iDim==3.and.k==nK)       DiffRight = 0.0
                end if

                Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim)   = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim) + DiffLeft
                Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim+1) = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,2*iDim+1) + DiffRight

                DiffLeft = -Eta_DFDB(jDir,i,j,k,iDim,iBlock) &
                     *DcoordDxyz_DDFD(iDim,iDir,i,j,k,iDim)*Coeff
                DiffRight = -Eta_DFDB(jDir,i+Di,j+Dj,k+Dk,iDim,iBlock) &
                     *DcoordDxyz_DDFD(iDim,iDir,i+Di,j+Dj,k+Dk,iDim)*Coeff

                Jacobian_VVCI(iVar,jVar,i,j,k,1) = &
                     Jacobian_VVCI(iVar,jVar,i,j,k,1) - (DiffLeft + DiffRight)

                if(UseNoOverlap)then
                   if(  iDim==1.and.i==1  .or. &
                        iDim==2.and.j==1  .or. &
                        iDim==3.and.k==1)        DiffLeft = 0.0
                   if(  iDim==1.and.i==nI .or. &
                        iDim==2.and.j==nJ .or. &
                        iDim==3.and.k==nK)       DiffRight = 0.0
                end if

                Jacobian_VVCI(iVar,jVar,i,j,k,2*iDim)   = &
                     Jacobian_VVCI(iVar,jVar,i,j,k,2*iDim) + DiffLeft
                Jacobian_VVCI(iVar,jVar,i,j,k,2*iDim+1) = &
                     Jacobian_VVCI(iVar,jVar,i,j,k,2*iDim+1) + DiffRight

             end do; end do
          end do; end do; end do
       end do
    end if

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine add_jacobian_resistivity
  !============================================================================

  subroutine add_jacobian_hall_resist(iBlock, nVarImpl, Jacobian_VVCI)

    ! Preconditioner for the induction equation in Hall MHD
    ! Based on Toth et al. "Hall MHD on Block Adaptive Grids", JCP 2008

    use BATL_lib,        ONLY: IsCartesianGrid, CellSize_DB, FaceNormal_DDFB,&
         CellVolume_GB
    use ModFaceGradient, ONLY: set_block_jacobian_face, DcoordDxyz_DDFD
    use ModImplicit,     ONLY: nStencil, UseNoOverlap
    use ModNumConst,     ONLY: i_DD, iLeviCivita_III
    use ModGeometry,     ONLY: true_cell
    use ModAdvance,      ONLY: B_
    use ModHallResist,   ONLY: HallCmaxFactor

    integer, intent(in):: iBlock
    integer, intent(in):: nVarImpl
    real, intent(inout):: Jacobian_VVCI(nVarImpl,nVarImpl,nI,nJ,nK,nStencil)

    integer:: iB
    integer:: iDim, iDir, jDir, i, j, k, i2, j2, k2, iSign
    integer:: iSub, iSup, iFace, kDim, lDir, jKlEpsilon, iKlEpsilon
    integer:: iVar, jVar
    real:: Term, TermSub, TermSup, InvDcoord2_D(nDim)

    real:: InvVolume

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'add_jacobian_hall_resist'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    ! Set the base index value for magnetic field variables
    if(UseSemiHallResist)then
       iB = 0
    else
       iB = B_
    end if

    if(HallCmaxFactor > 0 .and. UseSemiHallResist) then
       ! Jacobian = dR(B_iDir)/dB_iDir from the whistler diffusion
       ! R_i(B) = -[D_i+1/2 * (B_i+1 - B_i) - D_i-1/2 * (B_i - B_i-1)/V_i
       do iDim = 1, nDim

          iSub = 2*iDim    ! stencil index of subdiagonal elements
          iSup = iSub + 1  ! stencil index of superdiagonal elements

          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE

             ! Index of the "right" cell face
             i2 = i + i_DD(1,iDim)
             j2 = j + i_DD(2,iDim)
             k2 = k + i_DD(3,iDim)

             InvVolume = 1./CellVolume_GB(i,j,k,iBlock)
             ! dR/dB_i-1 = H_i-1/2 / V_i
             TermSub =  InvVolume*WhistlerCoeff_FDB(i,j,k,iDim,iBlock)
             ! dR/dB_i+1 = H_i+1/2 / V_i
             TermSup =  InvVolume*WhistlerCoeff_FDB(i2,j2,k2,iDim,iBlock)

             ! Main diagonal
             do iVar = iB+1, iB+MaxDim
                Jacobian_VVCI(iVar,iVar,i,j,k,1) = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,1) - (TermSub + TermSup)
             end do

             if(UseNoOverlap)then
                if(  iDim==1.and.i==1  .or. &
                     iDim==2.and.j==1  .or. &
                     iDim==3.and.k==1)        TermSub = 0.0
                if(  iDim==1.and.i==nI .or. &
                     iDim==2.and.j==nJ .or. &
                     iDim==3.and.k==nK)       TermSup = 0.0
             end if

             ! Off diagonals
             do iVar = iB+1, iB+MaxDim
                Jacobian_VVCI(iVar,iVar,i,j,k,iSub)   = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,iSub) + TermSub
                Jacobian_VVCI(iVar,iVar,i,j,k,iSup) = &
                     Jacobian_VVCI(iVar,iVar,i,j,k,iSup) + TermSup
             enddo
          enddo; enddo; enddo
       enddo
    end if

    ! the transverse part of curl(B) is ignored in the Jacobian
    if(IsCartesianGrid)then
       ! Is there something more to do for RZ geometry??? !!!
       InvDcoord2_D = 1/CellSize_DB(1:nDim,iBlock)**2

       ! Loop over face directions
       do iDim = 1, nDim

          iSub = 2*iDim    ! stencil index of subdiagonal elements
          iSup = iSub + 1  ! stencil index of superdiagonal elements

          ! Loop over cell centers
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.true_cell(i,j,k,iBlock)) CYCLE

             ! Index for the 'right' face of the cell
             i2 = i + i_DD(1,iDim)
             j2 = j + i_DD(2,iDim)
             k2 = k + i_DD(3,iDim)

             ! Loop over magnetic field flux components
             do iDir = 1, MaxDim
                ! There is a Levi_Civita symbol with iDim,iDir,jDir
                ! (see eqs. 52-53 with iDim=s, iDir=j, jDir=l)
                ! so all 3 have to be different. Their sum is 1+2+3=6.
                if(iDim == iDir) CYCLE

                ! Magnetic field component for partial derivative
                jDir = 6 - iDir - iDim

                ! Variable indexes for the Jacobian
                iVar = iDir + iB
                jVar = jDir + iB

                ! Get the sign
                iSign = iLeviCivita_III(iDim,iDir,jDir)

                ! B/(ne dCoord^2) from eqs 52-53
                TermSub  = iSign* &
                     InvDcoord2_D(iDim)*Bne_DFDB(iDim,i,j,k,iDim,iBlock)
                TermSup = iSign* &
                     InvDcoord2_D(iDim)*Bne_DFDB(iDim,i2,j2,k2,iDim,iBlock)

                ! Jacobian = dF(B_iDir)/dB_jDir
                ! Main diagonal (eq. 52)
                Jacobian_VVCI(iVar,jVar,i,j,k,1) = &
                     Jacobian_VVCI(iVar,jVar,i,j,k,1) - (TermSub + TermSup)

                if(UseNoOverlap)then
                   if(  iDim==1.and.i==1  .or. &
                        iDim==2.and.j==1  .or. &
                        iDim==3.and.k==1)        TermSub = 0.0
                   if(  iDim==1.and.i==nI .or. &
                        iDim==2.and.j==nJ .or. &
                        iDim==3.and.k==nK)       TermSup = 0.0
                end if

                ! Off diagonals (eq. 53)
                Jacobian_VVCI(iVar,jVar,i,j,k,iSub)   = &
                     Jacobian_VVCI(iVar,jVar,i,j,k,iSub) + TermSub
                Jacobian_VVCI(iVar,jVar,i,j,k,iSup) = &
                     Jacobian_VVCI(iVar,jVar,i,j,k,iSup) + TermSup

             end do
          end do; end do; end do
       end do
    else
       ! Jacobian for generalized coordinate (eqs. 50-51)

       ! Get dCoord/Dxyz
       call set_block_jacobian_face(iBlock)

       ! kDim is the index for gen.coord in DcoordDxyz(kDim,...)
       ! lDim is the index of the field component in dR(B_j)/d(B_l)
       do kDim = 1,nDim; do lDir = 1,3
          if(kDim == lDir) CYCLE

          ! jDir is the index for flux component in dR(B_j)/d(B_l)
          do jDir = 1, 3
             jklEpsilon = iLeviCivita_III(jDir,kDim,lDir)

             ! Variable indexes for the Jacobian
             iVar = jDir + iB
             jVar = lDir + iB

             ! Index for face normal components
             do iDim = 1, nDim
                if(iDim == jDir) CYCLE  ! Terms cancel out

                iklEpsilon = iLeviCivita_III(iDim, kDim, lDir)
                if(iklEpsilon == 0 .and. jklEpsilon == 0) CYCLE

                do iFace = 1, nDim  ! nDim directions of gen. coordinates

                   iSub = 2*iFace   ! stencil index of subdiagonal elements
                   iSup = iSub + 1  ! stencil index of superdiagonal elements

                   do k=1,nK; do j=1,nJ; do i=1,nI
                      if(.not.true_cell(i,j,k,iBlock)) CYCLE

                      ! Index for the opposite face of the cell
                      i2 = i + i_DD(1,iFace)
                      j2 = j + i_DD(2,iFace)
                      k2 = k + i_DD(3,iFace)

                      ! Eq. 51:
                      ! dR(Bj)/dBl = +- 1/(V*ds)*(Area_i*T_ks
                      !    *(Bi/ne*jklEpsilon - Bj/ne*iklEpsilon)^S

                      Term = -1.0/(CellVolume_GB(i,j,k,iBlock) &
                           * CellSize_DB(iFace,iBlock))

                      TermSub = Term &
                           *FaceNormal_DDFB(iDim,iFace,i,j,k,iBlock)&
                           *DcoordDxyz_DDFD(iFace,kDim,i,j,k,iFace) &
                           *(Bne_DFDB(iDim,i,j,k,iFace,iBlock)*jklEpsilon &
                           - Bne_DFDB(jDir,i,j,k,iFace,iBlock)*iklEpsilon)

                      TermSup = Term &
                           *FaceNormal_DDFB(iDim,iFace,i2,j2,k2,iBlock)&
                           *DcoordDxyz_DDFD(iFace,kDim,i2,j2,k2,iFace) &
                           *(Bne_DFDB(iDim,i2,j2,k2,iFace,iBlock)*jklEpsilon &
                           - Bne_DFDB(jDir,i2,j2,k2,iFace,iBlock)*iklEpsilon)

                      ! The main diagonal is -1 * the sum of off-diagonal term
                      Jacobian_VVCI(iVar,jVar,i,j,k,1) = &
                           Jacobian_VVCI(iVar,jVar,i,j,k,1) - TermSub - TermSup

                      if(UseNoOverlap)then
                         ! Exclude boundaries
                         if(  iFace==1.and.i==1  .or. &
                              iFace==2.and.j==1  .or. &
                              iFace==3.and.k==1)        TermSub = 0.0
                         if(  iFace==1.and.i==nI .or. &
                              iFace==2.and.j==nJ .or. &
                              iFace==3.and.k==nK)       TermSup = 0.0
                      end if

                      ! Off-diagonal terms
                      Jacobian_VVCI(iVar,jVar,i,j,k,iSub) = &
                           Jacobian_VVCI(iVar,jVar,i,j,k,iSub) + TermSub

                      Jacobian_VVCI(iVar,jVar,i,j,k,iSup) = &
                           Jacobian_VVCI(iVar,jVar,i,j,k,iSup) + TermSup

                   end do; end do; end do
                end do
             end do
          end do
       end do; end do

    end if

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine add_jacobian_hall_resist
  !============================================================================

  subroutine update_impl_resistivity(iBlock, NewSemiAll_VC)

    use ModAdvance,    ONLY: State_VGB
    use ModEnergy,     ONLY: calc_energy_cell
    use ModImplicit,   ONLY: nVarSemiAll
    use ModVarIndexes, ONLY: Bx_, Bz_
    use ModGeometry,   ONLY: true_cell

    integer, intent(in):: iBlock
    real,    intent(in):: NewSemiAll_VC(nVarSemiAll,nI,nJ,nK)

    integer :: i, j, k

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'update_impl_resistivity'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE
       State_VGB(Bx_:Bz_,i,j,k,iBlock) = NewSemiAll_VC(BxImpl_:BzImpl_,i,j,k)
    end do; end do; end do

    call calc_energy_cell(iBlock)

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine update_impl_resistivity
  !============================================================================

end module ModResistivity
!==============================================================================
