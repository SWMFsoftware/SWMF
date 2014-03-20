!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModResistivity

  ! Resistivity related variables and methods

  use ModSize,   ONLY: MaxBlock
  use BATL_size, ONLY: nDim, Maxdim, nI, nJ, nK, &
       MinI, MaxI, MinJ, MaxJ, MinK, MaxK

  implicit none
  save

  private ! except

  public:: read_resistivity_param
  public:: init_mod_resistivity
  public:: set_resistivity
  public:: mask_resistivity
  public:: calc_resistivity_source
  public:: calc_heat_exchange
  public:: init_impl_resistivity
  public:: get_impl_resistivity_state
  public:: get_resistivity_rhs
  public:: add_jacobian_resistivity
  public:: update_impl_resistivity

  logical, public           :: UseResistivity   = .false.
  logical, public           :: UseResistiveFlux = .false.
  logical, public           :: UseHeatExchange  = .true.
  character(len=30), public :: TypeResistivity='none'
  real, public, allocatable :: Eta_GB(:,:,:,:)
  real, public              :: Eta0, Eta0Si=0.0

  ! Local variables
  logical            :: UseJouleHeating  = .true.
  logical            :: DoMessagePassResistivity = .false.

  real               :: EtaPerpSpitzerSi = 0.0
  real               :: CoulombLogarithm = 20.0
  real               :: Eta0AnomSi=0.0, Eta0Anom
  real               :: EtaMaxAnomSi=0.0, EtaMaxAnom
  real               :: jCritAnomSi=1.0, jCritInv
  real, public       :: Si2NoEta
  real               :: EtaCoeff = 0.0
  real               :: jInvbCrit = 0.0

  real:: rZeroResist = -1.0
  real:: rFullResist = 1e30

  ! resistivity pre-multiplied by the face area vector
  real, allocatable :: Eta_DFDB(:,:,:,:,:,:)

  ! B/ne vector for Hall Resistivity
  real, allocatable :: Bne_DFDB(:,:,:,:,:,:)

  ! Named indices for semi-implicit variables
  integer, public, parameter :: BxImpl_ = 1, ByImpl_ = 2, BzImpl_ = 3

  ! Include resistive fluxes at all?
  logical :: DoResistiveFlux = .true.

  character(len=*), private, parameter :: NameMod = 'ModResistivity'
contains
  !===========================================================================
  subroutine read_resistivity_param(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'read_resistivity_param'
    !------------------------------------------------------------------------

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
             call read_var('jCritAnomSi',  jCritAnomSi)
          case('raeder')
             call read_var('EtaCoeff', EtaCoeff)
             call read_var('jInvbCrit', jInvbCrit)
          case default
             call stop_mpi(NameSub//': unknown TypeResistivity='&
                  //TypeResistivity)
          end select
       end if
    case("#RESISTIVITYREGION", "#RESISTIVEREGION")
       call read_var('rZeroResist', rZeroResist)
       call read_var('rFullResist', rFullResist)
    case default
       call stop_mpi(NameSub//': unknown command='//NameCommand)
    end select

  end subroutine read_resistivity_param

  !===========================================================================

  subroutine init_mod_resistivity

    use ModImplicit, ONLY: UseSemiImplicit, TypeSemiImplicit
    use ModPhysics,  ONLY: Si2No_V, UnitX_, UnitT_, UnitJ_
    use ModConst,    ONLY: cLightSpeed, cElectronCharge, &
         cElectronMass, cEps, cBoltzmann, cTwoPi

    logical:: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'init_mod_resistivity'
    !------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    UseResistiveFlux = DoResistiveFlux 

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

    if(.not.allocated(Eta_GB)) &
         allocate(Eta_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    if(TypeResistivity == 'constant')then
       Eta_GB = Eta0
    else
       Eta_GB = 0.0
    end if

    if(UseSemiImplicit .and. TypeSemiImplicit == 'resistivity')then
       ! The following will ensure that the explicit evaluation of the
       ! resistive diffusion is switched off
       UseResistiveFlux = .false.
    end if

    if(DoTestMe)then
       write(*,*)NameSub, ': DoResistiveFlux  = ', DoResistiveFlux
       write(*,*)NameSub, ': UseResistiveFlux = ', UseResistiveFlux
       write(*,*)NameSub, ': UseJouleHeating  = ', UseJouleHeating
       write(*,*)NameSub, ': UseHeatExchange  = ', UseHeatExchange
       write(*,*)NameSub, ': Si2NoEta = ',Si2NoEta
       write(*,*)NameSub, ': Si2NoJ   = ',Si2No_V(UnitJ_)
       write(*,*)NameSub, ': Eta0, Eta0Anom, EtaMaxAnom=', &
            Eta0, Eta0Anom, EtaMaxAnom
       write(*,*)NameSub, ': jCritInv = ', jCritInv
    end if

  end subroutine init_mod_resistivity

  !===========================================================================

  subroutine set_resistivity

    use ModMain,    ONLY: nBlock, Unused_B
    use BATL_lib,   ONLY: message_pass_cell

    character (len=*), parameter :: NameSub = 'set_resistivity'

    integer :: iBlock
    !--------------------------------------------------------------------------

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

       call mask_resistivity(iBlock, Eta_GB(:,:,:,iBlock))
    end do

    if(DoMessagePassResistivity) &
         call message_pass_cell(Eta_GB, nWidthIn=1)

  end subroutine set_resistivity

  !===========================================================================

  subroutine spitzer_resistivity(iBlock, Eta_G)

    use ModConst,   ONLY: cProtonMass, cElectronCharge
    use ModPhysics, ONLY: No2Si_V, UnitB_, UnitRho_, UnitTemperature_
    use ModAdvance, ONLY: State_VGB, Rho_, P_, Bx_, Bz_
    use ModB0,      ONLY: B0_DGB
    use ModMain,    ONLY: UseB0

    ! Compute Spitzer-type, classical resistivity 

    integer, intent(in) :: iBlock
    real,    intent(out):: Eta_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK)

    real :: EtaSi, Coef, B0_DG(3,MinI:MaxI, MinJ:MaxJ, MinK:MaxK)
    integer :: i, j, k
    !-------------------------------------------------------------------------

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
    real,    intent(out):: Eta_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK)

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

  subroutine mask_resistivity(iBlock, Eta_G)

    use ModGeometry, ONLY: rMin_BLK, r_BLK

    ! Mask Eta_G if required

    integer, intent(in)   :: iBlock
    real,    intent(inout):: Eta_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK)

    integer:: i, j, k
    real:: r
    !------------------------------------------------------------------------

    ! Check if there is any region with masked resistivity
    if(rZeroResist < 0.0) RETURN

    ! Check if the block is fully outside the masked region
    if(rMin_BLK(iBlock) > rFullResist) RETURN

    do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
       r = r_BLK(i,j,k,iBlock)
       if(r < rZeroResist)then
          Eta_G(i,j,k) = 0.0
       else if(r < rFullResist)then
          Eta_G(i,j,k) = Eta_G(i,j,k) &
               *(r - rZeroResist)/(rFullResist - rZeroResist)
       end if
    end do; end do; end do

  end subroutine mask_resistivity

  !===========================================================================
  subroutine raeder_resistivity(iBlock, Eta_G)

    ! Compute resistivity based on Raeder's formula
    ! Eta = Alpha * j'^2    j' = |j|*gridspacing/(|B| + epsilon)
    ! Eta = 0 if j' < jInvbCrit

    use ModCurrent,  ONLY: get_current
    use ModAdvance,  ONLY: State_VGB
    use ModB0,       ONLY: B0_DGB
    use ModVarIndexes, ONLY: Bx_, Bz_
    use ModMain, ONLY: UseB0

    integer, intent(in) :: iBlock
    real,    intent(out):: Eta_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    real :: Current_D(3), AbsJ, AbsB, JInvB, b_D(3)
    integer :: i, j, k, l, m, n
    !--------------------------------------------------------------------------
    ! Compute |J| and |B|
    do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
       call get_current(i,j,k,iBlock,Current_D)
       AbsJ = sqrt(sum(current_D**2))
       ! Calculate AbsB from the average of 26 neighboring cells and itself
       AbsB = 0.0
       do n=-1,1; do m=-1,1; do l=-1,1
          b_D = State_VGB(Bx_:Bz_,i+l,j+m,k+n,iBlock)
          if(UseB0) b_D = b_D + B0_DGB(:,i+l,j+m,k+n,iBlock)
          AbsB = AbsB + sqrt(sum(b_D**2))
       end do; end do; end do
       AbsB = AbsB/27.0
       ! Make J/B dimensionless by multiplying length scale 1Re
       JInvB = min(AbsJ*1.0/(AbsB + 1e-8), 1.0) 

       ! Compute Raeder resistivity
       if(JInvB >= jInvbCrit)then
          Eta_G(i,j,k) = EtaCoeff*(JInvB**2) 
       else
          Eta_G(i,j,k) = 0.0
       end if
    end do; end do; end do

  end subroutine raeder_resistivity

  !===========================================================================
  subroutine calc_resistivity_source(iBlock)

    use ModMain,       ONLY: x_
    use ModGeometry,   ONLY: true_cell
    use BATL_lib,      ONLY: IsRzGeometry, Xyz_DGB
    use ModCurrent,    ONLY: get_current
    use ModPhysics,    ONLY: Gm1, Inv_Gm1
    use ModVarIndexes, ONLY: p_, Pe_, Ppar_, Bz_, Energy_
    use ModAdvance,    ONLY: Source_VC, &
         UseElectronPressure, UseAnisoPressure

    integer, intent(in):: iBlock

    ! Variables needed for Joule heating
    real :: Current_D(3), JouleHeating

    integer:: i, j, k
    !-----------------------------------------------------------------------
    JouleHeating = 0.0

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE

       if(UseJouleHeating)then
          call get_current(i,j,k,iBlock,Current_D)
          JouleHeating = gm1 * Eta_GB(i,j,k,iBlock) * sum(Current_D**2)
       end if
       if(UseElectronPressure) then
          Source_VC(Pe_,i,j,k) = Source_VC(Pe_,i,j,k) + JouleHeating
       else
          Source_VC(P_,i,j,k) = Source_VC(P_,i,j,k) + JouleHeating

          ! the same amount of Joule heating applies on Ppar
          if(UseAnisoPressure) &
               Source_VC(Ppar_,i,j,k)  = Source_VC(Ppar_,i,j,k) + JouleHeating

          if(.not.UseResistiveFlux) &
               Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
               + Inv_Gm1*JouleHeating
       end if

       ! rz-geometrical source terms
       if(IsRzGeometry .and. UseResistiveFlux)then
          ! calculate cell centered current if it is not yet done
          if(.not.UseJouleHeating) call get_current(i,j,k,iBlock,Current_D)

          ! Source[Bphi] = -eta*Jz / radius
          Source_VC(Bz_,i,j,k) = Source_VC(Bz_,i,j,k) &
               - Eta_GB(i,j,k,iBlock)*Current_D(x_)/Xyz_DGB(2,i,j,k,iBlock)
       end if
    end do; end do; end do

  end subroutine calc_resistivity_source

  !============================================================================

  subroutine calc_heat_exchange

    use ModMain,       ONLY: Cfl, nBlock, Unused_B
    use ModGeometry,   ONLY: true_cell
    use ModPhysics,    ONLY: gm1, IonMassPerCharge
    use ModVarIndexes, ONLY: Rho_, p_, Pe_, Ppar_
    use ModAdvance,    ONLY: time_blk, State_VGB, &
         UseAnisoPressure

    real :: DtLocal
    real :: HeatExchange, HeatExchangePeP, HeatExchangePePpar
    integer:: i, j, k, iBlock
    !--------------------------------------------------------------------------
    HeatExchange = 0.0
    HeatExchangePeP = 0.0
    HeatExchangePePpar = 0.0

    call set_resistivity

    do iBlock = 1, nBlock
       if (Unused_B(iBlock)) CYCLE

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          DtLocal = Cfl*time_BLK(i,j,k,iBlock)

          ! For single ion fluid the ion-electron collision results in a 
          ! heat exchange term for the electron pressure 
          ! See eq. 4.124c in Schunk and Nagy.

          ! Explicit heat exchange
          HeatExchange = gm1 * Eta_GB(i,j,k,iBlock) * &
               3*State_VGB(Rho_,i,j,k,iBlock)*(1./IonMassPerCharge**2)

          ! Point-implicit correction for stability: H' = H/(1+2*dt*H)
          HeatExchange = HeatExchange / (1 + 2.0*DtLocal*HeatExchange)

          HeatExchangePeP = HeatExchange &
               *(State_VGB(P_,i,j,k,iBlock) - State_VGB(Pe_,i,j,k,iBlock))

          ! Heat exchange for parallel ion pressure
          if(UseAnisoPressure)then
             HeatExchangePePpar = HeatExchange &
                  *(State_VGB(Ppar_,i,j,k,iBlock) &
                  - State_VGB(Pe_,i,j,k,iBlock))

             State_VGB(Ppar_,i,j,k,iBlock) = State_VGB(Ppar_,i,j,k,iBlock) &
                  - DtLocal*HeatExchangePePpar
          end if

          ! Heat exchange for the ions
          State_VGB(P_,i,j,k,iBlock) = State_VGB(P_,i,j,k,iBlock) &
               - DtLocal*HeatExchangePeP

          ! Heat exchange for the electrons
          State_VGB(Pe_,i,j,k,iBlock) = State_VGB(Pe_,i,j,k,iBlock) &
               + DtLocal*HeatExchangePeP
       end do; end do; end do
    end do

  end subroutine calc_heat_exchange

  !============================================================================
  ! Interface for (Semi-)implicit collisional and Hall resistivity
  !============================================================================

  subroutine get_impl_resistivity_state(SemiAll_VGB)

    use ModSize,       ONLY: j0_, nJp1_, k0_, nKp1_
    use ModAdvance,    ONLY: State_VGB
    use ModImplicit,   ONLY: nVarSemiAll, nBlockSemi, iBlockFromSemi_I
    use ModVarIndexes, ONLY: Bx_, Bz_

    real, intent(out):: &
         SemiAll_VGB(nVarSemiAll,0:nI+1,j0_:nJp1_,k0_:nKp1_,nBlockSemi)

    integer:: i, j, k, iBlock, iBlockSemi

    character(len=*), parameter :: &
         NameSub = 'ModResistivity::get_impl_resistivity_state'
    !--------------------------------------------------------------------------

    ! Copy magnetic field into the implicit variable
    do iBlockSemi = 1, nBlockSemi
       iBlock = iBlockFromSemi_I(iBlockSemi)

       ! Store the magnetic field in SemiAll_VGB
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          SemiAll_VGB(BxImpl_:BzImpl_,i,j,k,iBlockSemi) = &
               State_VGB(Bx_:Bz_,i,j,k,iBlock)
       end do; end do; end do

    end do

    ! Initialize variables for RHS and Jacobian calculation
    call init_impl_resistivity

  end subroutine get_impl_resistivity_state

  !===========================================================================

  subroutine init_impl_resistivity

    use BATL_lib,        ONLY: IsCartesian, message_pass_cell, &
         CellFace_DB, FaceNormal_DDFB, nBlock, Unused_B
    use ModAdvance,      ONLY: State_VGB
    use ModMain,         ONLY: UseB0
    use ModNumConst,     ONLY: i_DD
    use ModVarIndexes,   ONLY: Rho_, Bx_, Bz_
    use ModHallResist,   ONLY: UseHallResist, hall_factor, &
         set_ion_mass_per_charge, IonMassPerCharge_G
    use ModB0,           ONLY: B0_DGB

    integer:: iDim, i, j, k, Di, Dj, Dk, i1, j1, k1, iBlock
    real:: Eta
    real:: FaceNormal_D(nDim)

    real:: HallCoeff, Rho, b_D(MaxDim)

    character(len=*), parameter :: &
         NameSub = 'ModResistivity::init_impl_resistivity'
    !--------------------------------------------------------------------------
    if(UseResistivity)then
       ! Calculate Eta*FaceNormal_D to be used for resistive flux
       if(.not.allocated(Eta_DFDB)) &
            allocate(Eta_DFDB(nDim,nI+1,nJ+1,nK+1,nDim,MaxBlock))

       ! Message pass resistivity to fill in ghost cells
       call message_pass_cell(Eta_GB, DoSendCornerIn=.false.)

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

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
    end if

    if(UseHallResist)then

       ! Calculate B/ne = m/e * B/rho on the face to be used for the Hall flux
       if(.not.allocated(Bne_DFDB)) &
            allocate(Bne_DFDB(MaxDim,nI+1,nJ+1,nK+1,nDim,MaxBlock))

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          call set_ion_mass_per_charge(iBlock)

          do iDim = 1, nDim
             Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
             do k = 1, nK+Dk; do j = 1, nJ+Dj; do i = 1, nI+Di

                ! Check if the Hall coefficient is positive for this face
                HallCoeff = hall_factor(iDim, i, j, k, iBlock) 
                if(HallCoeff <= 0.0)then
                   Bne_DFDB(:,i,j,k,iDim,iBlock) = 0.0
                else
                   ! Cell center indexes for the cell on the left of the face
                   i1 = i-Di; j1 = j-Dj; k1 = k - Dk

                   ! Average m/e, rho and B to the face
                   HallCoeff = HallCoeff* &
                        0.5*( IonMassPerCharge_G(i1,j1,k1) &
                        +     IonMassPerCharge_G(i,j,k)    )

                   Rho = 0.5*(State_VGB(Rho_,i1,j1,k1,iBlock) &
                        +     State_VGB(Rho_,i,j,k,iBlock)    )

                   b_D = 0.5*(State_VGB(Bx_:Bz_,i1,j1,k1,iBlock) &
                        +     State_VGB(Bx_:Bz_,i,j,k,iBlock)    )

                   ! B0 on the face is actually average of cell center B0
                   ! It is easier to use B0_DGB then then the 3 face arrays
                   ! The res change on the coarse side will get overwritten
                   if(UseB0) b_D = b_D + &
                        0.5*(B0_DGB(:,i1,j1,k1,iBlock)+ B0_DGB(:,i,j,k,iBlock))
                
                   ! Calculate B/ne for the face
                   Bne_DFDB(:,i,j,k,iDim,iBlock) = HallCoeff*b_D/Rho
                end if
             end do; end do; end do
          end do
       end do
    end if

  end subroutine init_impl_resistivity

  !============================================================================

  subroutine get_resistivity_rhs(iBlock, StateImpl_VG, Rhs_VC, IsLinear)

    use BATL_lib,        ONLY: store_face_flux, IsCartesian, IsRzGeometry, &
         Xyz_DGB, CellSize_DB, CellVolume_GB, CellFace_DB, FaceNormal_DDFB
    use ModFaceGradient, ONLY: get_face_curl
    use ModImplicit,     ONLY: nVarSemi, &
         FluxImpl_VXB, FluxImpl_VYB, FluxImpl_VZB
    use ModNumConst,     ONLY: i_DD
    use ModSize,         ONLY: x_, y_, z_
    use ModGeometry,     ONLY: true_cell, true_BLK
    use ModHallResist,   ONLY: UseHallResist

    integer, intent(in) :: iBlock
    real, intent(inout) :: StateImpl_VG(nVarSemi,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(out)   :: Rhs_VC(nVarSemi,nI,nJ,nK)
    logical, intent(in) :: IsLinear

    ! resistive flux for operator split scheme
    real, allocatable, save :: FluxImpl_VFD(:,:,:,:,:)

    integer :: iDim, i, j, k, Di, Dj, Dk
    real    :: Current_D(MaxDim), Jx, InvDy2, FaceNormal_D(nDim)
    real    :: jNormal, BneNormal
    logical :: IsNewBlock
    !--------------------------------------------------------------------------

    IsNewBlock = .true.

    Rhs_VC = 0.0

    if(.not. true_BLK(iblock)) RETURN

    if(.not.allocated(FluxImpl_VFD)) allocate( &
         FluxImpl_VFD(nVarSemi,nI+1,nJ+1,nK+1,nDim))

    ! Loop over face directions
    do iDim = 1, nDim
       Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)

       if(UseHallResist .and. IsCartesian) then
          FaceNormal_D = 0.0; FaceNormal_D(iDim) = CellFace_DB(iDim,iBlock)
       end if

       ! Loop over cell faces orthogonal to iDim 
       do k = 1, nK+Dk; do j = 1, nJ+Dj; do i = 1, nI+Di
          if(.not.true_cell(i,j,k,iBlock)) CYCLE

          call get_face_curl(iDim, i, j, k, iBlock, IsNewBlock, StateImpl_VG, &
               Current_D)

          FluxImpl_VFD(BxImpl_:BzImpl_,i,j,k,iDim) = 0.0

          if(UseResistivity)then
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
          if(UseHallResist)then
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
             BneNormal = sum(FaceNormal_D*Bne_DFDB(1:nDim,i,j,k,iDim,iBlock))

             ! Flux = Bn/ne J - Jn B/ne
             FluxImpl_VFD(BxImpl_:BzImpl_,i,j,k,iDim) = &
                  FluxImpl_VFD(BxImpl_:BzImpl_,i,j,k,iDim) &
                  - Jnormal*Bne_DFDB(:,i,j,k,iDim,iBlock) &
                  + BneNormal*Current_D
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

  end subroutine get_resistivity_rhs

  !============================================================================

  subroutine add_jacobian_resistivity(iBlock, nVarImpl, Jacobian_VVCI)

    ! Calculate the Jacobian for the preconditioning of 
    ! collisional and Hall resistivity.

    use ModProcMH,     ONLY: iProc
    use ModHallResist, ONLY: UseHallResist
    use ModImplicit,   ONLY: UseSemiImplicit, nStencil
    use BATL_lib,      ONLY: rot_to_cart
    use ModMain,       ONLY: iTest, jTest, kTest, BlkTest, ProcTest
    use ModAdvance,    ONLY: State_VGB, Bx_, Bz_, B_

    integer, intent(in):: iBlock
    integer, intent(in):: nVarImpl
    real, intent(inout):: Jacobian_VVCI(nVarImpl,nVarImpl,nI,nJ,nK,nStencil)

    integer:: iStencil, iB, i, j, k

    real, allocatable, save:: Jac_VVI(:,:,:)

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'add_jacobian_resistivity'
    !--------------------------------------------------------------------------
    if(.not.(UseResistivity .or. UseHallResist)) RETURN

    if(iProc == ProcTest .and. iBlock == BlkTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    ! Set the base index value for magnetic field variables
    if(UseSemiImplicit)then
       iB = 0
    else
       iB = B_
    end if

    if(DoTestMe)then
       write(*,*) NameSub,' starting with UseResistivity, UseHallResist=', &
            UseResistivity, UseHallResist
       do iStencil = 1, nStencil
          write(*,'(a,i1,a,10es13.5)') &
               NameSub//': unrot JAC(:,:,TestCell,', iStencil, ')=', &
               rot_to_cart( &
               Jacobian_VVCI(iB+1:iB+3,iB+1:iB+3,iTest,jTest,kTest,iStencil))
       end do
    end if

    if(UseResistivity) &
         call add_jacobian_eta_resist(iBlock, iB, nVarImpl, Jacobian_VVCI)

    if(UseHallResist) &
         call add_jacobian_hall_resist(iBlock, iB, nVarImpl, Jacobian_VVCI)

    if(DoTestMe)then
       write(*,*) NameSub, ' finished with:'
       i = iTest; j = jTest; k = kTest
       do iStencil = 1, nStencil
          write(*,'(a,i1,a,10es13.5)') &
               NameSub//': unrot JAC(:,:,TestCell,', iStencil, ')=', &
               rot_to_cart(Jacobian_VVCI(iB+1:iB+3,iB+1:iB+3,i,j,k,iStencil))
       end do

       if(.false.)then !!! UseSemiImplicit)then
          if(.not.allocated(Jac_VVI)) &
               allocate(Jac_VVI(nVarImpl,nVarImpl,nStencil))

          !!! Circular dependency to be resolved
          !call test_semi_impl_jacobian( &
          !     State_VGB(Bx_:Bz_,:,:,:,iBlock), &
          !     1e-4, get_resistivity_rhs, Jac_VVI)

          do iStencil = 1, nStencil
             write(*,'(a,i1,a,10es13.5)') &
                  NameSub//': numer JAC(:,:,TestCell,', iStencil, ')=', &
                  rot_to_cart(Jac_VVI(:,:,iStencil))
          end do

          do iStencil = 1, nStencil
             write(*,'(a,i1,a,10es13.5)') &
                  NameSub//': error JAC(:,:,TestCell,', iStencil, ')=', &
                  rot_to_cart(Jac_VVI(:,:,iStencil) &
                  - Jacobian_VVCI(iB+1:iB+3,iB+1:iB+3,i,j,k,iStencil))
          end do
       end if
    end if

  end subroutine add_jacobian_resistivity

  !========================================================================

  subroutine add_jacobian_eta_resist(iBlock, iB, nVarImpl, Jacobian_VVCI)

    ! Calculate the Jacobian for the preconditioning of 
    ! collisional resistivity.

    use BATL_lib,        ONLY: IsCartesianGrid, CellSize_DB, CellVolume_GB
    use ModFaceGradient, ONLY: set_block_jacobian_face, DcoordDxyz_DDFD
    use ModImplicit,     ONLY: UseNoOverlap, nStencil
    use ModNumConst,     ONLY: i_DD
    use ModGeometry,     ONLY: true_cell

    integer, intent(in):: iBlock
    integer, intent(in):: iB
    integer, intent(in):: nVarImpl
    real, intent(inout):: Jacobian_VVCI(nVarImpl,nVarImpl,nI,nJ,nK,nStencil)

    integer :: iDim, iDir, jDir, i, j, k, Di, Dj, Dk
    integer:: iVar, jVar
    real :: DiffLeft, DiffRight, InvDcoord_D(nDim), Coeff

    !--------------------------------------------------------------------------
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

  end subroutine add_jacobian_eta_resist

  !============================================================================

  subroutine add_jacobian_hall_resist(iBlock, iB, nVarImpl, Jacobian_VVCI)

    ! Preconditioner for the induction equation in Hall MHD
    ! Based on Toth et al. "Hall MHD on Block Adaptive Grids", JCP 2008

    use BATL_lib,        ONLY: IsCartesianGrid, CellSize_DB, FaceNormal_DDFB,&
         CellVolume_GB
    use ModFaceGradient, ONLY: set_block_jacobian_face, DcoordDxyz_DDFD
    use ModImplicit,     ONLY: nStencil, UseNoOverlap
    use ModNumConst,     ONLY: i_DD, iLeviCivita_III
    use ModGeometry,     ONLY: true_cell

    integer, intent(in):: iBlock
    integer, intent(in):: nVarImpl
    integer, intent(in):: iB
    real, intent(inout):: Jacobian_VVCI(nVarImpl,nVarImpl,nI,nJ,nK,nStencil)

    integer:: iDim, iDir, jDir, i, j, k, i2, j2, k2, iSign
    integer:: iSub, iSup, iFace, kDim, lDir, jklEpsilon, iklEpsilon
    integer:: iVar, jVar
    real:: Term, TermSub, TermSup, InvDcoord2_D(nDim)

    character(len=*), parameter:: NameSub = 'add_jacobian_hall_resist'
    !--------------------------------------------------------------------------
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

          ! jDri is the index for flux component in dR(B_j)/d(B_l)
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

  end subroutine add_jacobian_hall_resist

  !============================================================================

  subroutine update_impl_resistivity(iBlock, SemiAll_VC)

    use ModAdvance,    ONLY: State_VGB
    use ModEnergy,     ONLY: calc_energy_cell
    use ModImplicit,   ONLY: nVarSemiAll
    use ModVarIndexes, ONLY: Bx_, Bz_
    use ModGeometry,   ONLY: true_cell

    integer, intent(in):: iBlock
    real,    intent(in):: SemiAll_VC(nVarSemiAll,nI,nJ,nK)

    integer :: i, j, k
    !--------------------------------------------------------------------------

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE
       State_VGB(Bx_:Bz_,i,j,k,iBlock) = SemiAll_VC(BxImpl_:BzImpl_,i,j,k)
    end do; end do; end do

    call calc_energy_cell(iBlock)

  end subroutine update_impl_resistivity

end module ModResistivity
