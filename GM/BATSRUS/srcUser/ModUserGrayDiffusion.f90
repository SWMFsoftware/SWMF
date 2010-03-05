!^CFG COPYRIGHT UM
!============================================================================
module ModUser
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_init_session,               &
       IMPLEMENTED3 => user_normalization,              &
       IMPLEMENTED4 => user_set_ics,                    &
       IMPLEMENTED5 => user_set_plot_var,               &
       IMPLEMENTED6 => user_material_properties,        &
       IMPLEMENTED7 => user_update_states,              &
       IMPLEMENTED8 => user_set_outerbcs,               &
       IMPLEMENTED9 => user_amr_criteria

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'Hd with Gray-Diffusion'

  ! local variables for the reference solutions of Lowrie's
  ! radiative shock tube tests
  integer :: iLowrieTest
  character(len=100) :: NameLowrieFile
  integer :: nCellLowrie
  integer :: &
       iRhoLowrie  = 1, &
       iUxLowrie   = 2, &
       iTgasLowrie = 3, &
       iTradLowrie = 4, &
       nVarLowrie  = 4
  real, allocatable :: xLowrie_C(:), StateLowrie_VC(:,:)
  real :: U0, X0
  real :: EradBc1, EradBc2
  
  real, parameter :: Gamma = 5.0/3.0

contains
!============================================================================
  subroutine user_read_inputs

    use ModIO,          ONLY: write_prefix, write_myname, iUnitOut
    use ModMain,        ONLY: lVerbose
    use ModProcMH,      ONLY: iProc
    use ModReadParam,   ONLY: read_line, read_command, read_var

    character (len=100) :: NameCommand
    !------------------------------------------------------------------------
    if(iProc == 0 .and. lVerbose > 0)then
       call write_prefix;
       write(iUnitOut,*)'User read_input starts'
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
       case("#LOWRIETEST")
          call read_var('iLowrieTest',iLowrieTest)
          call read_var('NameLowrieFile',NameLowrieFile)

       case('#USERINPUTEND')
          if(iProc == 0 .and. lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input ends'
          endif
          EXIT

       case default
          if(iProc == 0) then
             call write_myname; write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) '  *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do

  end subroutine user_read_inputs

  !============================================================================

  subroutine user_init_session

    use ModIoUnit,  ONLY: UnitTmp_
    use ModPhysics, ONLY: cRadiationNo

    integer :: iError, iCell
    real :: Mach, Entropy

    character(len=*), parameter :: NameSub = "user_init_session"
    !--------------------------------------------------------------------------

    select case(iLowrieTest)
    case(1)
       ! Mach 1.05 test
       nCellLowrie = 6188
       U0 = -1.05               ! veloxity added to lowrie's solution
       X0 = 0.03675             ! shift in x-direction 
    case(2)
       ! Mach 2 test
       nCellLowrie = 4840
       U0 = -2.0
       X0 = 0.0125
    case(3)
       ! Mach 5 test with variable cross-sections
       nCellLowrie = 8904
       U0 = -5.0
       X0 = 0.02
    case default
       call stop_mpi(NameSub // " wrong iLowrieTest")
    end select

    allocate( xLowrie_C(nCellLowrie), StateLowrie_VC(nVarLowrie,nCellLowrie) )

    open(UnitTmp_, FILE=NameLowrieFile, STATUS="old", IOSTAT=iError)

    if(iError /= 0)call stop_mpi(NameSub // &
         " could not open Lowrie's file="//NameLowrieFile)

    do iCell = 1, nCellLowrie
       read(UnitTmp_,*) xLowrie_C(iCell), &
            StateLowrie_VC(:,iCell), Mach, Entropy
    end do

    close(UnitTmp_)

    ! The reference solutions use p=rho*T/gamma
    StateLowrie_VC(iTgasLowrie,:) = StateLowrie_VC(iTgasLowrie,:)/Gamma
    StateLowrie_VC(iTradLowrie,:) = StateLowrie_VC(iTradLowrie,:)/Gamma

    EradBc1 = cRadiationNo*StateLowrie_VC(iTradLowrie,1)**4
    EradBc2 = cRadiationNo*StateLowrie_VC(iTradLowrie,nCellLowrie)**4

  end subroutine user_init_session

  !==========================================================================

  subroutine user_normalization

    use ModConst,   ONLY: cRadiation, cProtonMass, cBoltzmann
    use ModPhysics, ONLY: No2Si_V, UnitRho_, UnitU_
    
    character (len=*), parameter :: NameSub = 'user_normalization'
    !------------------------------------------------------------------------

    No2Si_V = 1.0
    ! The following density unit is needed to get a normalized radiation
    ! constant  with value 1.0e-4. The gamma dependence is needed since
    ! the reference solution uses p=rho*T/gamma
    No2Si_V(UnitRho_) = 1.0e+4*cRadiation*(cProtonMass/cBoltzmann)**4 &
         *No2Si_V(UnitU_)**6/Gamma**4

  end subroutine user_normalization

  !==========================================================================

  subroutine user_set_ics

    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: x_Blk, y_Blk
    use ModMain,       ONLY: GlobalBlk, nI, nJ, nK, x_, y_
    use ModPhysics,    ONLY: ShockSlope, No2Si_V, Si2No_V, &
         UnitTemperature_, UnitEnergyDens_, inv_gm1, cRadiationNo
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUy_, RhoUz_, Erad_, &
         ExtraEint_, p_

    integer :: iBlock, i, j, k, iCell
    real :: x, Weight1, Weight2
    real :: Rho, Ux, Tgas, Trad, p, Erad, RhoU_D(3)
    real :: SinSlope, CosSlope, Rot_II(2,2)

    character(len=*), parameter :: NameSub = "user_set_ics"
    !------------------------------------------------------------------------

    iBlock = GlobalBlk

    ! Calculate sin and cos from the tangent = ShockSlope
    SinSlope = ShockSlope/sqrt(1.0+ShockSlope**2)
    CosSlope =        1.0/sqrt(1.0+ShockSlope**2)
    ! Set rotational matrix
    Rot_II = reshape( (/CosSlope, SinSlope, -SinSlope, CosSlope/), &
         (/2,2/) )

    do j = 1, nJ; do i = 1, nI
       
       x = x_Blk(i,j,0,iBlock)*CosSlope + y_Blk(i,j,0,iBlock)*SinSlope - X0

       do iCell = 1, nCellLowrie
          if(xLowrie_C(iCell) >= x) EXIT
       end do
       if(iCell == 1) call stop_mpi(NameSub // &
            " Lowrie solution does not cover the left boundary")

       if(iCell > nCellLowrie)then
          ! Cell is beyond the last point of Lowrie input: use last cell
          iCell   = nCellLowrie
          Weight1 = 0.0
          Weight2 = 1.0
       else
          ! Assign weights for linear interpolation between iCell-1 and iCell
          Weight1 = (xLowrie_C(iCell) - x) &
               /    (xLowrie_C(iCell) - xLowrie_C(iCell-1))
          Weight2 = 1.0 - Weight1
       end if

       Rho  = ( Weight1*StateLowrie_VC(iRhoLowrie, iCell-1) &
            +   Weight2*StateLowrie_VC(iRhoLowrie, iCell) )
       Ux   = ( Weight1*StateLowrie_VC(iUxLowrie, iCell-1) &
            +   Weight2*StateLowrie_VC(iUxLowrie, iCell) )
       Tgas = ( Weight1*StateLowrie_VC(iTgasLowrie, iCell-1) &
            +   Weight2*StateLowrie_VC(iTgasLowrie, iCell) )
       Trad = ( Weight1*StateLowrie_VC(iTradLowrie, iCell-1) &
            +   Weight2*StateLowrie_VC(iTradLowrie, iCell) )

       p = Rho*Tgas
       Erad = cRadiationNo*Trad**4
       RhoU_D(1) = Rho*(Ux+U0)
       RhoU_D(2) = 0.0
       RhoU_D(3) = 0.0

       do k = 1, nk
          State_VGB(Rho_,i,j,k,iBlock) = Rho
          State_VGB(RhoUx_:RhoUy_,i,j,k,iBlock) = matmul(Rot_II,RhoU_D(x_:y_))
          State_VGB(RhoUz_,i,j,k,iBlock) = 0.0
          State_VGB(Erad_,i,j,k,iBlock) = Erad
          State_VGB(ExtraEint_,i,j,k,iBlock) = &
               (1.0/(Gamma-1.0) - inv_gm1)*p
          State_VGB(p_,i,j,k,iBlock) = p
       end do

    end do; end do

  end subroutine user_set_ics

  !===========================================================================

  subroutine user_set_outerbcs(iBlock,iSide, TypeBc, IsFound)

    use ModAdvance,    ONLY: State_VGB
    use ModImplicit,   ONLY: StateSemi_VGB, iEradImpl
    use ModMain,       ONLY: nI
    use ModVarIndexes, ONLY: Erad_

    integer,          intent(in)  :: iBlock, iSide
    character(len=20),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer :: j, k

    character (len=*), parameter :: NameSub = 'user_set_outerbcs'
    !-------------------------------------------------------------------------

    if(.not. (iSide==1 .or. iSide==2) )then
       write(*,*) NameSub//' : user boundary not defined at iSide = ', iSide
       call stop_mpi(NameSub)
    end if

    select case(TypeBc)
    case('user')
       ! For shocktube slope < 1, float and shear bc's have the same effect
       ! on the boundaries in the x-direction. For convenience, float is used
       ! instead of shear.
       select case(iSide)
       case(1)
          ! float for all variables except radiation
          State_VGB(:, 0,:,:,iBlock) = State_VGB(:,1,:,:,iBlock)
          State_VGB(:,-1,:,:,iBlock) = State_VGB(:,1,:,:,iBlock)

          State_VGB(Erad_,-1:0,:,:,iBlock) = EradBc1
       case(2)
          ! float for all variables except radiation
          State_VGB(:,nI+1,:,:,iBlock) = State_VGB(:,nI,:,:,iBlock)
          State_VGB(:,nI+2,:,:,iBlock) = State_VGB(:,nI,:,:,iBlock)

          State_VGB(Erad_,nI+1:nI+2,:,:,iBlock) = EradBc2
       end select
    case('usersemi')
       select case(iSide)
       case(1)
          StateSemi_VGB(iEradImpl,0,:,:,iBlock) = EradBc1
       case(2)
          StateSemi_VGB(iEradImpl,nI+1,:,:,iBlock) = EradBc2
       end select
    end select

    IsFound = .true.

  end subroutine user_set_outerbcs

  !===========================================================================

  subroutine user_update_states(iStage,iBlock)

    use ModSize,    ONLY: nI, nJ, nK
    use ModAdvance, ONLY: State_VGB, p_, ExtraEint_, &
         UseNonConservative, IsConserv_CB, &
         Source_VC, uDotArea_XI, uDotArea_YI, uDotArea_ZI
    use ModGeometry, ONLY: vInv_CB
    use ModPhysics,  ONLY: g, inv_gm1, Si2No_V, No2Si_V, &
         UnitP_, UnitEnergyDens_
    use ModEnergy,  ONLY: calc_energy_cell

    implicit none

    integer, intent(in):: iStage,iBlock

    integer:: i, j, k
    real   :: PressureSi, EinternalSi, GammaEos, DivU
    logical:: IsConserv

    character(len=*), parameter :: NameSub = 'user_update_states'
    !------------------------------------------------------------------------
    ! Fix adiabatic compression source for pressure
    if(UseNonConservative)then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          DivU          =        uDotArea_XI(i+1,j,k,1) - uDotArea_XI(i,j,k,1)
          if(nJ>1) DivU = DivU + uDotArea_YI(i,j+1,k,1) - uDotArea_YI(i,j,k,1)
          if(nK>1) DivU = DivU + uDotArea_ZI(i,j,k+1,1) - uDotArea_ZI(i,j,k,1)
          DivU = vInv_CB(i,j,k,iBlock)*DivU

          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, GammaOut=GammaEos)

          Source_VC(p_,i,j,k) = Source_VC(p_,i,j,k) &
               -(GammaEos-g)*State_VGB(p_,i,j,k,iBlock)*DivU
       end do; end do; end do
    end if

    call update_states_MHD(iStage,iBlock)

    ! update of pressure, ionization and total energies
    do k=1,nK; do j=1,nJ; do i=1,nI
       ! Total internal energy ExtraEint + P/(\gamma -1) transformed to SI

       if(allocated(IsConserv_CB))then
          IsConserv = IsConserv_CB(i,j,k,iBlock)
       else
          IsConserv = .not. UseNonConservative
       end if

       if(IsConserv)then
          ! At this point p=(g-1)(e-rhov^2/2) with the ideal gamma g.
          ! Use this p to get total internal energy density.
          EinternalSi = No2Si_V(UnitEnergyDens_)*&
               (inv_gm1*State_VGB(P_,i,j,k,iBlock) + &
               State_VGB(ExtraEint_,i,j,k,iBlock))
          call user_material_properties(State_VGB(:,i,j,k,iBlock),&
               EinternalIn=EinternalSi, PressureOut=PressureSi)
      
          ! Set true pressure
          State_VGB(p_,i,j,k,iBlock) = PressureSi*Si2No_V(UnitP_)
       else
          call user_material_properties(State_VGB(:,i,j,k,iBlock),&
               EinternalOut=EinternalSi)
       end if

       ! Set ExtraEint = Total internal energy - P/(gamma -1)
       State_VGB(ExtraEint_,i,j,k,iBlock) = &
            Si2No_V(UnitEnergyDens_)*EinternalSi &
            - inv_gm1*State_VGB(p_,i,j,k,iBlock)

    end do; end do; end do

    call calc_energy_cell(iBlock)

  end subroutine user_update_states

  !===========================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use BATL_size,     ONLY: nI, nJ, nK, nG, MinI, MaxI
    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: x_Blk, y_Blk
    use ModMain,       ONLY: Time_Simulation
    use ModPhysics,    ONLY: Si2No_V, No2Si_V, UnitT_, &
         UnitTemperature_, UnitEnergyDens_, ShockSlope, cRadiationNo
    use ModVarIndexes, ONLY: p_, Rho_, Erad_

    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(-1:nI+2, -1:nJ+2, -1:nK+2)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    character (len=*), parameter :: NameSub = 'user_set_plot_var'

    integer :: i, j, k, iCell
    real :: x, Weight1, Weight2
    real :: Ux, Tgas, Trad
    real :: SinSlope, CosSlope
    integer, parameter:: jMin = 1 - 2*min(1,nJ-1), jMax = nJ + 2*min(1,nJ-1)
    integer, parameter:: kMin = 1 - 2*min(1,nK-1), kMax = nK + 2*min(1,nK-1)
    !------------------------------------------------------------------------

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

    IsFound = .true.
    select case(NameVar)
    case('tgas')
       do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
          PlotVar_G(i,j,k) = State_VGB(p_,i,j,k,iBlock) &
               /State_VGB(Rho_,i,j,k,iBlock)
       end do; end do; end do

    case('trad')
       do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI
          PlotVar_G(i,j,k) = (State_VGB(Erad_,i,j,k,iBlock)/cRadiationNo)**0.25
       end do; end do; end do

    case('rho0','ux0','uy0','tgas0','trad0')

       ! Calculate sin and cos from the tangent = ShockSlope
       SinSlope = ShockSlope/sqrt(1.0+ShockSlope**2)
       CosSlope =        1.0/sqrt(1.0+ShockSlope**2)

       do k = kMin, kMax; do j = jMin, jMax; do i = MinI, MaxI

          x = x_Blk(i,j,k,iBlock)*CosSlope + y_Blk(i,j,k,iBlock)*SinSlope
          x = x - U0*Time_Simulation*Si2No_V(UnitT_) - X0

          do iCell = 1, nCellLowrie
             if(xLowrie_C(iCell) >= x) EXIT
          end do
          if(iCell == 1) call stop_mpi(NameSub // &
               " Lowrie solution does not cover the left boundary")

          if(iCell > nCellLowrie)then
             ! Cell is beyond the last point of Lowrie input: use last cell
             iCell   = nCellLowrie
             Weight1 = 0.0
             Weight2 = 1.0
          else
             ! Assign weights for linear interpolation between
             ! iCell-1 and iCell
             Weight1 = (xLowrie_C(iCell) - x) &
                  /    (xLowrie_C(iCell) - xLowrie_C(iCell-1))
             Weight2 = 1.0 - Weight1
          end if

          select case(NameVar)
          case('rho0')

             PlotVar_G(i,j,k) = &
                  ( Weight1*StateLowrie_VC(iRhoLowrie, iCell-1) &
                  + Weight2*StateLowrie_VC(iRhoLowrie, iCell) )

          case('ux0','uy0')

             Ux   = ( Weight1*StateLowrie_VC(iUxLowrie, iCell-1) &
                  +   Weight2*StateLowrie_VC(iUxLowrie, iCell) )
             
             select case(NameVar)
             case('ux0')
                PlotVar_G(i,j,k) = (Ux+U0)*CosSlope
             case('uy0')
                PlotVar_G(i,j,k) = (Ux+U0)*SinSlope
             end select

          case('tgas0')

             Tgas = ( Weight1*StateLowrie_VC(iTgasLowrie, iCell-1) &
                  +   Weight2*StateLowrie_VC(iTgasLowrie, iCell) )

             PlotVar_G(i,j,k) = Tgas

          case('trad0')

             Trad = ( Weight1*StateLowrie_VC(iTradLowrie, iCell-1) &
                  +   Weight2*StateLowrie_VC(iTradLowrie, iCell) )

             PlotVar_G(i,j,k) = Trad

          end select
       end do; end do; end do

    case default
       IsFound = .false.
    end select

  end subroutine user_set_plot_var

  !==========================================================================

  subroutine user_material_properties(State_V, i, j, k, iBlock, iDir, &
       EinternalIn, TeIn, NatomicOut, &
       EinternalOut, TeOut, PressureOut, &
       CvOut, GammaOut, HeatCondOut, IonHeatCondOut, TeTiRelaxOut, &
       OpacityPlanckOut_W, OpacityRosselandOut_W, PlanckOut_W)

    ! The State_V vector is in normalized units

    use ModAdvance,    ONLY: nWave
    use ModConst,      ONLY: cLightSpeed
    use ModPhysics,    ONLY: No2Si_V, Si2No_V, &
         UnitTemperature_, UnitEnergyDens_, UnitT_, UnitU_, UnitX_, UnitP_
    use ModVarIndexes, ONLY: nVar, Rho_, p_

    real, intent(in) :: State_V(nVar)
    integer, optional, intent(in):: i, j, k, iBlock, iDir  ! cell/face index
    real, optional, intent(in)  :: EinternalIn             ! [J/m^3]
    real, optional, intent(in)  :: TeIn                    ! [K]
    real, optional, intent(out) :: NatomicOut              ! [1/m^3]
    real, optional, intent(out) :: EinternalOut            ! [J/m^3]
    real, optional, intent(out) :: TeOut                   ! [K]
    real, optional, intent(out) :: PressureOut             ! [Pa]
    real, optional, intent(out) :: CvOut                   ! [J/(K*m^3)]
    real, optional, intent(out) :: GammaOut                ! dimensionless
    real, optional, intent(out) :: HeatCondOut             ! [J/(m*K*s)]
    real, optional, intent(out) :: IonHeatCondOut          ! [J/(m*K*s)]
    real, optional, intent(out) :: TeTiRelaxOut            ! [1/s]
    real, optional, intent(out) :: &
         OpacityPlanckOut_W(nWave)                         ! [1/m]
    real, optional, intent(out) :: &
         OpacityRosselandOut_W(nWave)                      ! [1/m]

    ! Multi-group specific interface. The variables are respectively:
    !  Group Planckian spectral energy density
    real, optional, intent(out) :: PlanckOut_W(nWave)      ! [J/m^3]

    real :: Temperature, OpacityPlanck, DiffusionRad

    character (len=*), parameter :: NameSub = 'user_material_properties'
    !-------------------------------------------------------------------

    if(present(EinternalIn))then
       Temperature = EinternalIn*Si2No_V(UnitEnergyDens_) &
            *(Gamma - 1.0)/State_V(Rho_)
    else
       Temperature = State_V(p_)/State_V(Rho_)
    end if

    if(present(EinternalOut)) EinternalOut = &
         State_V(Rho_)*Temperature/(Gamma - 1.0) *No2Si_V(UnitEnergyDens_)

    if(present(PressureOut)) &
         PressureOut = State_V(Rho_)*Temperature*No2Si_V(UnitP_)

    if(present(CvOut)) CvOut = State_V(Rho_)/(Gamma - 1.0) &
         *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_)

    select case(iLowrieTest)
    case(1,2)
       DiffusionRad = 1.0
       OpacityPlanck = 1.0E6
    case(3)
       DiffusionRad = 0.00175*(Gamma*Temperature)**3.5/State_V(Rho_)
       OpacityPlanck = 1.0E6/DiffusionRad
    end select

    if(present(GammaOut)) GammaOut = Gamma

    if(present(TeOut)) TeOut = Temperature*No2Si_V(UnitTemperature_)

    if(present(OpacityPlanckOut_W)) OpacityPlanckOut_W = &
         OpacityPlanck/No2Si_V(UnitT_)/cLightSpeed

    if(present(OpacityRosselandOut_W)) OpacityRosselandOut_W = &
         cLightSpeed/(3.0*DiffusionRad*No2Si_V(UnitU_)*No2Si_V(UnitX_))

    if(present(HeatCondOut)) HeatCondOut = 0.0
    if(present(TeTiRelaxOut)) TeTiRelaxOut = 0.0

  end subroutine user_material_properties

  !============================================================================

  subroutine user_amr_criteria(iBlock, UserCriteria, TypeCriteria, IsFound)

    use ModSize,       ONLY: nI, nJ, nK
    use ModAdvance,    ONLY: State_VGB, Rho_, p_
    use ModAMR,        ONLY: RefineCritMin_I, CoarsenCritMax
    use ModGeometry,   ONLY: y_Blk, Dy_BLK, Dx_BLK, y1, y2
    use ModPhysics,    ONLY: cRadiationNo
    use ModVarIndexes, ONLY: Rho_, p_, Erad_

    ! Variables required by this user subroutine
    integer, intent(in)          :: iBlock
    real, intent(out)            :: UserCriteria
    character (len=*),intent(in) :: TypeCriteria
    logical ,intent(inout)       :: IsFound

    real, parameter:: TemperatureMin = 5.2, DTradDxMin = 1e3
    real:: Temperature, DTradDx, TradL, TradR, DyRefine
    integer:: i, j, k
    !--------------------------------------------------------------------------

    DyRefine = (y2 - y1)/8

    UserCriteria = 0.0

    ! capture the embedded hydro shock
    LOOPCELL: do k = 1, nK; do j=1, nJ; do i = -1, nI+2
       Temperature = State_VGB(p_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
       if(Temperature > TemperatureMin &
            .and. abs(y_Blk(i,j,k,iBlock)) < DyRefine)then
          UserCriteria = 1.0
          EXIT LOOPCELL
       end if
    end do; end do; end do LOOPCELL

    ! capture the precursor
    LOOPCELL2: do k = 1, nK; do j=1, nJ; do i = 0, nI+1
       TradL = sqrt(sqrt(State_VGB(Erad_,i-1,j,k,iBlock)/cRadiationNo))
       TradR = sqrt(sqrt(State_VGB(Erad_,i+1,j,k,iBlock)/cRadiationNo))
       DTradDx = (TradR - TradL)*0.5/Dx_Blk(iBlock)
       if(DTradDx > DTradDxMin .and. abs(y_Blk(i,j,k,iBlock)) < DyRefine)then
          UserCriteria = 1.0
          EXIT LOOPCELL2
       end if
    end do; end do; end do LOOPCELL2

    ! Do not refine blocks far from discontinuity (crit=0.0)
    ! Do not coarsen blocks near discontinuity    (crit=1.0)
    RefineCritMin_I = 0.5
    CoarsenCritMax  = 0.5

    IsFound = .true.

  end subroutine user_amr_criteria

end module ModUser
