!^CFG FILE DISSFLUX
!^CFG COPYRIGHT UM
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
  public:: spitzer_resistivity
  public:: anomalous_resistivity
  public:: mask_resistivity
  public:: calc_resistivity_source
  public:: get_impl_resistivity_state
  public:: get_resistivity_rhs
  public:: add_jacobian_resistivity
  public:: update_impl_resistivity

  logical, public           :: UseResistivity   = .false.
  logical, public           :: UseResistiveFlux = .true.
  character(len=30), public :: TypeResistivity='none'
  real, public, allocatable :: Eta_GB(:,:,:,:)
  real, public              :: Eta0, Eta0Si=0.0
  
  ! Local variables
  logical            :: UseJouleHeating  = .true.
  logical            :: UseHeatExchange  = .true.

  real               :: EtaPerpSpitzerSi = 0.0
  real               :: CoulombLogarithm = 20.0
  real               :: Eta0AnomSi=0.0, Eta0Anom
  real               :: EtaMaxAnomSi=0.0, EtaMaxAnom
  real               :: jCritAnomSi=1.0, jCritInv
  real               :: Si2NoEta

  real:: rZeroResist = -1.0
  real:: rFullResist = 1e30

  character(len=*), private, parameter :: NameMod = 'ModResistivity'

  ! Array needed for third-order interpolation of ghost cells
  real, allocatable :: State1_VG(:,:,:,:)

  ! resistivity pre-multiplied by the face area
  real, allocatable :: Eta_DFDB(:,:,:,:,:,:)

  ! resistive flux for operator split scheme
  real, allocatable :: FluxImpl_VFD(:,:,:,:,:)

  ! Named indices for semi-implicit variables
  integer, parameter :: BxImpl_ = 1, ByImpl_ = 2, BzImpl_ = 3

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
    case("#RESISTIVITYREGION", "#RESISTIVEREGION")
       call read_var('rZeroResist', rZeroResist)
       call read_var('rFullResist', rFullResist)
    case default
       call stop_mpi(NameSub//': unknown command='//NameCommand)
    end select

  end subroutine read_resistivity_param

  !===========================================================================

  subroutine init_mod_resistivity

    use BATL_lib,    ONLY: IsCartesian, IsRzGeometry
    use ModConst,    ONLY: cLightSpeed, cElectronCharge, cElectronMass, cEps, &
         cBoltzmann, cTwoPi
    use ModImplicit, ONLY: UseSemiImplicit, TypeSemiImplicit, nVarSemi
    use ModPhysics,  ONLY: Si2No_V, UnitX_, UnitT_, UnitJ_
    use ModVarIndexes, ONLY: nVar

    character(len=*), parameter :: &
         NameSub = 'ModResistivity::init_mod_resistivity'
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

    if(UseSemiImplicit .and. TypeSemiImplicit == 'resistivity')then
       ! The following will ensure that the explicit evaluation of the
       ! resistive diffusion is switched off
       UseResistiveFlux = .false.

       allocate( &
            State1_VG(nVar,-1:nI+2,-1:nJ+2,-1:nK+2), &
            FluxImpl_VFD(nVarSemi,nI+1,nJ+1,nK+1,nDim), &
            Eta_DFDB(MaxDim,nI+1,nJ+1,nK+1,nDim,MaxBlock) )

       if(.not.(IsCartesian.or.IsRzGeometry)) &
            call stop_mpi(NameSub//': semi-implicit resistivity not '// &
            'yet implemented for non-cartesian')
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

  subroutine mask_resistivity(iBlock, Eta_G)

    use ModGeometry, ONLY: rMin_BLK, r_BLK

    ! Mask Eta_G if required

    integer, intent(in)   :: iBlock
    real,    intent(inout):: Eta_G(-1:nI+2, -1:nJ+2, -1:nK+2)

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
  subroutine calc_resistivity_source(iBlock)

    use ModMain,       ONLY: Cfl, x_
    use ModGeometry,   ONLY: true_cell, IsRzGeometry, y_BLK
    use ModCurrent,    ONLY: get_current
    use ModPhysics,    ONLY: gm1, inv_gm1, IonMassPerCharge
    use ModVarIndexes, ONLY: Rho_, p_, Pe_, Ppar_, Energy_, Bz_
    use ModAdvance,    ONLY: time_blk, State_VGB, Source_VC, &
         UseElectronPressure, UseAnisoPressure

    integer, intent(in):: iBlock

    ! Variables needed for Joule heating
    real :: Current_D(3), JouleHeating, HeatExchange

    integer:: i, j, k
    !-----------------------------------------------------------------------
    JouleHeating = 0.0
    HeatExchange = 0.0

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
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
       else
          Source_VC(P_,i,j,k) = Source_VC(P_,i,j,k) + JouleHeating
       end if

       ! rz-geometrical source terms
       if(IsRzGeometry .and. UseResistiveFlux)then
          ! calculate cell centered current if it is not yet done
          if(.not.UseJouleHeating) call get_current(i,j,k,iBlock,Current_D)

          ! Source[Bphi] = -eta*Jz / radius
          Source_VC(Bz_,i,j,k) = Source_VC(Bz_,i,j,k) &
               - Eta_GB(i,j,k,iBlock)*Current_D(x_)/y_BLK(i,j,k,iBlock)
       end if
    end do; end do; end do

  end subroutine calc_resistivity_source

  !============================================================================
  ! Operator split, semi-implicit subroutines
  !============================================================================

  subroutine get_impl_resistivity_state(StateImpl_VGB, DconsDsemi_VCB)

    use BATL_lib,        ONLY: IsCartesian, message_pass_cell
    use ModAdvance,      ONLY: State_VGB
    use ModImplicit,     ONLY: nw, nImplBLK, impl2iBlk
    use ModMain,         ONLY: MaxImplBlk
    use ModNumConst,     ONLY: i_DD
    use ModVarIndexes,   ONLY: nVar, Bx_, Bz_

    real, intent(out) :: StateImpl_VGB(nw,0:nI+1,0:nJ+1,0:nK+1,MaxImplBlk)
    real, intent(inout) :: DconsDsemi_VCB(nw,nI,nJ,nK,MaxImplBlk)

    integer :: iDim, i, j, k, Di, Dj, Dk, iBlock, iImplBlock
    real :: Eta
    real :: area, Normal_D(MaxDim)

    character(len=*), parameter :: &
         NameSub = 'ModResistivity::get_impl_resistivity_state'
    !--------------------------------------------------------------------------

    do iImplBlock = 1, nImplBLK
       iBlock = impl2iBLK(iImplBlock)

       ! Store the magnetic field in StateImpl_VGB
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          StateImpl_VGB(BxImpl_:BzImpl_,i,j,k,iImplBlock) = &
               State_VGB(Bx_:Bz_,i,j,k,iBlock)
       end do; end do; end do

       ! Calculate the cell-centered resistivity
       call set_resistivity(iBlock)
    end do

    ! Message pass to fill in ghost cells
    call message_pass_cell(1, Eta_GB, DoSendCornerIn=.false.)

    do iImplBlock = 1, nImplBLK
       iBlock = impl2iBLK(iImplBlock)

       do iDim = 1, nDim
          Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
          if(IsCartesian) call set_cartesian_cell_face(iDim, iBlock)
          do k = 1, nK+Dk; do j = 1, nJ+Dj; do i = 1, nI+Di
             if(.not.IsCartesian) &
                  call set_covariant_cell_face(iDim,i,j,k,iBlock)

             Eta = 0.5*(Eta_GB(i,j,k,iBlock) + Eta_GB(i-Di,j-Dj,k-Dk,iBlock))

             Eta_DFDB(:,i,j,k,iDim,iBlock) = Area*Eta*Normal_D
          end do; end do; end do
       end do
    end do

  contains

    subroutine set_cartesian_cell_face(iDim, iBlock)

      use BATL_lib, ONLY: CellFace_DB

      integer, intent(in) :: iDim, iBlock
      !------------------------------------------------------------------------

      Area = CellFace_DB(iDim,iBlock)
      Normal_D = 0.0; Normal_D(iDim) = 1.0

    end subroutine set_cartesian_cell_face

    !==========================================================================

    subroutine set_covariant_cell_face(iDim, i, j, k, iBlock)

      use ModGeometry, ONLY: FaceAreaI_DFB, FaceAreaJ_DFB, FaceAreaK_DFB, &
           FaceArea2MinI_B, FaceArea2MinJ_B, FaceArea2MinK_B, &
           x_BLK, y_BLK, z_BLK
      use ModSize,     ONLY: x_, y_, z_

      integer, intent(in) :: iDim, i, j, k, iBlock

      integer :: iLeft, jLeft, kLeft
      real :: Area_D(3), Area2, Area2Min
      !------------------------------------------------------------------------

      select case(iDim)
      case(x_)
         iLeft = i-1; jLeft = j; kLeft = k
         Area_D = FaceAreaI_DFB(:,i,j,k,iBlock)
         Area2Min = FaceArea2MinI_B(iBlock)
      case(y_)
         iLeft = i; jLeft = j-1; kLeft = k
         Area_D = FaceAreaJ_DFB(:,i,j,k,iBlock)
         Area2Min = FaceArea2MinJ_B(iBlock)
      case(z_)
         iLeft = i; jLeft = j; kLeft = k-1
         Area_D = FaceAreaK_DFB(:,i,j,k,iBlock)
         Area2Min = FaceArea2MinK_B(iBlock)
      end select

      Area2 = sum(Area_D**2)

      if(Area2 < 0.5*Area2Min)then
         ! The face is at the pole
         Normal_D(x_) = x_BLK(i,j,k,iBlock) - x_BLK(iLeft,jLeft,kLeft,iBlock)
         Normal_D(y_) = y_BLK(i,j,k,iBlock) - y_BLK(iLeft,jLeft,kLeft,iBlock)
         Normal_D(z_) = z_BLK(i,j,k,iBlock) - z_BLK(iLeft,jLeft,kLeft,iBlock)
         Normal_D = Normal_D/sqrt(sum(Normal_D**2))
         Area = 0.0
      else
         Area = sqrt(Area2)
         Normal_D = Area_D/Area
      end if

    end subroutine set_covariant_cell_face

  end subroutine get_impl_resistivity_state

  !============================================================================

  subroutine get_resistivity_rhs(iBlock, StateImpl_VG, Rhs_VC, IsLinear)

    use BATL_lib,        ONLY: store_face_flux
    use ModFaceGradient, ONLY: get_face_curl
    use ModGeometry,     ONLY: vInv_CB
    use ModImplicit,     ONLY: nVarSemi, FluxImpl_VXB, FluxImpl_VYB, &
         FluxImpl_VZB
    use ModNumConst,     ONLY: i_DD
    use ModSize,         ONLY: x_, y_, z_

    integer, intent(in) :: iBlock
    real, intent(inout) :: StateImpl_VG(nVarSemi,-1:nI+2,-1:nJ+2,-1:nK+2)
    real, intent(out)   :: Rhs_VC(nVarSemi,nI,nJ,nK)
    logical, intent(in) :: IsLinear

    integer :: iDim, i, j, k, Di, Dj, Dk
    real :: Current_D(MaxDim)
    logical :: IsNewBlock
    !--------------------------------------------------------------------------

    IsNewBlock = .true.

    do iDim = 1, nDim
       Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
       do k = 1, nK+Dk; do j = 1, nJ+Dj; do i = 1, nI+Di

          call get_face_curl(iDim, i, j, k, iBlock, IsNewBlock, StateImpl_VG, &
               Current_D)

          FluxImpl_VFD(BxImpl_,i,j,k,iDim) = &
               + Eta_DFDB(y_,i,j,k,iDim,iBlock)*Current_D(z_) &
               - Eta_DFDB(z_,i,j,k,iDim,iBlock)*Current_D(y_)
          FluxImpl_VFD(ByImpl_,i,j,k,iDim) = &
               + Eta_DFDB(z_,i,j,k,iDim,iBlock)*Current_D(x_) &
               - Eta_DFDB(x_,i,j,k,iDim,iBlock)*Current_D(z_)
          FluxImpl_VFD(BzImpl_,i,j,k,iDim) = &
               + Eta_DFDB(x_,i,j,k,iDim,iBlock)*Current_D(y_) &
               - Eta_DFDB(y_,i,j,k,iDim,iBlock)*Current_D(x_)

       end do; end do; end do
    end do

    ! Store the fluxes at resolution changes for restoring conservation
    call store_face_flux(iBlock, nVarSemi, FluxImpl_VFD, &
         FluxImpl_VXB, FluxImpl_VYB, FluxImpl_VZB)

    Rhs_VC = 0.0

    do iDim = 1, nDim
       Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          Rhs_VC(:,i,j,k) = Rhs_VC(:,i,j,k) - vInv_CB(i,j,k,iBlock) &
               *(FluxImpl_VFD(:,i+Di,j+Dj,k+Dk,iDim) &
               - FluxImpl_VFD(:,i,j,k,iDim))
       end do; end do; end do
    end do

  end subroutine get_resistivity_rhs

  !============================================================================

  subroutine add_jacobian_resistivity(iBlock, Jacobian_VVCI)

    use BATL_lib,        ONLY: CellSize_DB
    use ModGeometry,     ONLY: vInv_CB
    use ModImplicit,     ONLY: nVarSemi
    use ModNumConst,     ONLY: i_DD

    integer, parameter:: nStencil = 2*nDim + 1

    integer, intent(in) :: iBlock
    real, intent(inout) :: Jacobian_VVCI(nVarSemi,nVarSemi,nI,nJ,nK,nStencil)

    integer :: iDim, iDir, i, j, k, Di, Dj, Dk, iB
    real :: DiffLeft, DiffRight, InvDcoord_D(MaxDim)
    !--------------------------------------------------------------------------

    InvDcoord_D = 1/CellSize_DB(:,iBlock)

    ! the transverse diffusion is ignored in the Jacobian
    do iDim = 1, nDim
       Di = i_DD(iDim,1); Dj = i_DD(iDim,2); Dk = i_DD(iDim,3)
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          do iDir = 1, MaxDim
             if(iDim == iDir) CYCLE

             iB = iDir

             DiffLeft = vInv_CB(i,j,k,iBlock) &
                  *Eta_DFDB(iDim,i,j,k,iDim,iBlock)*InvDcoord_D(iDim)
             DiffRight = vInv_CB(i,j,k,iBlock) &
                  *Eta_DFDB(iDim,i+Di,j+Dj,k+Dk,iDim,iBlock)*InvDcoord_D(iDim)

             Jacobian_VVCI(iB,iB,i,j,k,1) = &
                  Jacobian_VVCI(iB,iB,i,j,k,1) - (DiffLeft + DiffRight)

             if(iDim==1.and.i==1 .or. iDim==2.and.j==1 .or. iDim==3.and.k==1) &
                  DiffLeft = 0.0
             if(iDim==1.and.i==nI .or. iDim==2.and.j==nJ .or. &
                  iDim==3.and.k==nK) DiffRight = 0.0

             Jacobian_VVCI(iB,iB,i,j,k,2*iDim)   = &
                  Jacobian_VVCI(iB,iB,i,j,k,2*iDim) + DiffLeft
             Jacobian_VVCI(iB,iB,i,j,k,2*iDim+1) = &
                  Jacobian_VVCI(iB,iB,i,j,k,2*iDim+1) + DiffRight
          end do
       end do; end do; end do
    end do

  end subroutine add_jacobian_resistivity

  !============================================================================

  subroutine update_impl_resistivity(iBlock, iImplBlock, StateImpl_VG)

    use ModAdvance,    ONLY: State_VGB
    use ModEnergy,     ONLY: calc_energy_cell
    use ModImplicit,   ONLY: nw
    use ModVarIndexes, ONLY: Bx_, Bz_

    integer, intent(in) :: iBlock, iImplBlock
    real, intent(in) :: StateImpl_VG(nw,nI,nJ,nK)

    integer :: i, j, k
    !--------------------------------------------------------------------------

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       State_VGB(Bx_:Bz_,i,j,k,iBlock) = StateImpl_VG(BxImpl_:BzImpl_,i,j,k)
    end do; end do; end do

    call calc_energy_cell(iBlock)

  end subroutine update_impl_resistivity

end module ModResistivity
