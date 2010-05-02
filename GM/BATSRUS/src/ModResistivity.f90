!^CFG FILE DISSFLUX
!^CFG COPYRIGHT UM
module ModResistivity

  ! Resistivity related variables and methods

  use ModSize, ONLY: nI, nJ, nK, MaxBlock

  implicit none

  logical            :: UseResistivity   = .false.
  logical            :: UseResistiveFlux = .true.
  logical            :: UseJouleHeating  = .true.
  logical            :: UseHeatExchange  = .true.

  character (len=30) :: TypeResistivity='none'
  real               :: Eta0Si=0.0, Eta0
  real               :: EtaPerpSpitzerSi = 0.0
  real               :: CoulombLogarithm = 20.0
  real               :: Eta0AnomSi=0.0, Eta0Anom
  real               :: EtaMaxAnomSi=0.0, EtaMaxAnom
  real               :: jCritAnomSi=1.0, jCritInv
  real               :: Si2NoEta

  real, allocatable :: Eta_GB(:,:,:,:)

  character(len=*), private, parameter :: NameMod = 'ModResistivity'

contains
  !===========================================================================
  subroutine read_resistivity_param(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand
    
    character(len=*), parameter:: NameSub = 'read_resistivity_param'
    !------------------------------------------------------------------------
    select case(NameCommand)
    case("#RESISTIVITYOPTIONS")
       call read_var("UseResistiveFlux", UseResistiveFlux)
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
             call read_var('jCritAnomSi',  jCritAnomSi)
          case default
             call stop_mpi(NameSub//': unknown TypeResistivity='&
                  //TypeResistivity)
          end select
       end if
    case default
       call stop_mpi(NameSub//': unknown command='//NameCommand)
    end select

  end subroutine read_resistivity_param

  !===========================================================================

  subroutine init_mod_resistivity

    use ModConst, ONLY: cLightSpeed, cElectronCharge, cElectronMass, cEps, &
         cBoltzmann, cTwoPi

    use ModPhysics, ONLY: Si2No_V, UnitX_, UnitT_, UnitJ_
    !------------------------------------------------------------------------
    Si2NoEta = Si2No_V(UnitX_)**2/Si2No_V(UnitT_)

    Eta0       = Eta0Si       * Si2NoEta
    Eta0Anom   = Eta0AnomSi   * Si2NoEta
    EtaMaxAnom = EtaMaxAnomSi * Si2NoEta
    jCritInv   = 1.0/(jCritAnomSi  * Si2No_V(UnitJ_))

    ! Spitzer resistivity coefficient with Coulomb logarithm in SI units
    EtaPerpSpitzerSi = 0.5*sqrt(cElectronMass)    * &
         (cLightSpeed*cElectronCharge)**2/ &
         (3*cEps*(cTwoPi*cBoltzmann)**1.5) &
         *CoulombLogarithm

    if(.not.allocated(Eta_GB))then
       allocate(Eta_GB(-1:nI+2,-1:nJ+2,-1:nK+2,MaxBlock))

       if(TypeResistivity == 'constant')then
          Eta_GB = Eta0
       else
          Eta_GB = 0.0
       end if
    end if

  end subroutine init_mod_resistivity

  !===========================================================================

  subroutine spitzer_resistivity(iBlock, Eta_G)

    use ModConst,   ONLY: cProtonMass, cElectronCharge
    use ModPhysics, ONLY: No2Si_V, UnitB_, UnitRho_, UnitTemperature_
    use ModAdvance, ONLY: State_VGB, Rho_, P_, Bx_, By_, Bz_, &
         B0_DGB
    use ModMain, ONLY: UseB0

    ! Compute Spitzer-type, classical resistivity 

    integer, intent(in) :: iBlock
    real,    intent(out):: Eta_G(-1:nI+2, -1:nJ+2, -1:nK+2)

    real :: EtaSi, Coef, B0_DG(3,-1:nI+2, -1:nJ+2, -1:nK+2)
    integer :: i, j, k
    !-------------------------------------------------------------------------

    Coef =((cProtonMass/cElectronCharge)*No2Si_V(UnitB_)/No2Si_V(UnitRho_))**2
    if(UseB0)then
       B0_DG=B0_DGB(:,:,:,:,iBlock)
    else
       B0_DG=0.00
    end if
    do k=-1, nK+2; do j=-1, nJ+2; do i=-1, nI+2

       ! EtaSi = EtaPerpSpitzer/Te^1.5 in SI units
       EtaSi = EtaPerpSpitzerSi &
            / ( No2Si_V(UnitTemperature_)* &
            State_VGB(P_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock))**1.5

       ! Take into account the dependence on the B field:
       !    Eta' = Eta*(1 + [B*mp/(rho*e*Eta)]^2)
       EtaSi = EtaSi * (1.0 + Coef*sum( &
            (State_VGB(Bx_:Bz_,i,j,k,iBlock)+B0_DG(:,i,j,k))**2) &
            / (State_VGB(Rho_,i,j,k,iBlock)*EtaSi)**2)

       ! Normalize Eta_G
       Eta_G(i,j,k) = EtaSi*Si2NoEta
    end do; end do; end do

  end subroutine spitzer_resistivity

  !============================================================================

  subroutine anomalous_resistivity(iBlock, Eta_G)

    ! Compute current dependent anomalous resistivity

    use ModCurrent, ONLY: get_current

    integer, intent(in) :: iBlock
    real,    intent(out):: Eta_G(-1:nI+2, -1:nJ+2, -1:nK+2)

    real :: Current_D(3), AbsJ
    integer :: i, j, k
    !-------------------------------------------------------------------------

    ! Compute the magnitude of the current density |J|
    do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1;

       call get_current(i,j,k,iBlock,Current_D)
       AbsJ = sqrt(sum(current_D**2))

       ! Compute the anomalous resistivity:: 
       ! Eta = Eta0 + Eta0Anom*(|J|/Jcrit-1) 

       Eta_G(i,j,k) = Eta0 + &
            min(EtaMaxAnom, max(0.0, Eta0Anom*(AbsJ*JcritInv - 1.0)))

    end do; end do; end do

  end subroutine anomalous_resistivity

  !===========================================================================
  subroutine calc_resistivity_source(iBlock)

    use ModMain,       ONLY: Cfl
    use ModGeometry,   ONLY: true_cell
    use ModCurrent,    ONLY: get_current
    use ModPhysics,    ONLY: gm1, inv_gm1, IonMassPerCharge
    use ModVarIndexes, ONLY: Rho_, p_, Pe_, Ppar_, Energy_
    use ModAdvance,    ONLY: time_blk, State_VGB, Source_VC, &
         UseElectronPressure, UseAnisoPressure

    integer, intent(in):: iBlock

    ! Variables needed for Joule heating
    real :: Current_D(3), JouleHeating, HeatExchange

    integer:: i, j, k
    !-----------------------------------------------------------------------
    JouleHeating = 0.0
    HeatExchange = 0.0

    do k=1,nK; do j=1,nJ; do i=1,nI

       if(.not.true_cell(i,j,k,iBlock)) CYCLE
       if(UseJouleHeating)then
          call get_current(i,j,k,iBlock,Current_D)
          JouleHeating = gm1 * Eta_GB(i,j,k,iBlock) * sum(Current_D**2)
       end if
       if(UseElectronPressure) then
          ! For single ion fluid the ion-electron collision results in a 
          ! heat exchange term for the electron pressure 
          ! See eq. 4.124c in Schunk and Nagy.
          if(UseHeatExchange)then
             ! Explicit heat exchange
             HeatExchange = gm1 * Eta_GB(i,j,k,iBlock) * &
                  3*State_VGB(Rho_,i,j,k,iBlock)*(1./IonMassPerCharge**2)

             ! Point-implicit correction for stability: H' = H/(1+dt*H)
             HeatExchange = HeatExchange / &
                  (1 + Cfl*HeatExchange*time_BLK(i,j,k,iBlock)) &
                  *(State_VGB(P_,i,j,k,iBlock) &
                  - State_VGB(Pe_,i,j,k,iBlock))
          end if

          Source_VC(Pe_,i,j,k) = Source_VC(Pe_,i,j,k) &
               + JouleHeating + HeatExchange

          ! Heat exchange applies to ions too
          Source_VC(P_,i,j,k) = Source_VC(P_,i,j,k) - HeatExchange

          if(UseAnisoPressure) then
             ! Heat exchange for parallel ion pressure
             Source_VC(Ppar_,i,j,k) = Source_VC(Ppar_,i,j,k) - HeatExchange

             ! Relaxation term due to collisions
             Source_VC(Ppar_,i,j,k)  = Source_VC(Ppar_,i,j,k)  &
                  - Eta_GB(i,j,k,iBlock) &
                  *State_VGB(Rho_,i,j,k,iBlock)/IonMassPerCharge**2 &
                  *(State_VGB(Ppar_,i,j,k,iBlock) - State_VGB(p_,i,j,k,iBlock))
          end if

          ! Remove Joule heating and apply heat exchange to ion energy
          Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
               - inv_gm1*(JouleHeating + HeatExchange)
       else
          Source_VC(P_,i,j,k) = Source_VC(P_,i,j,k) + JouleHeating
       end if
    end do; end do; end do

  end subroutine calc_resistivity_source

end module ModResistivity

