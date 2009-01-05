!^CFG COPYRIGHT UM
!============================================================================
module ModUser
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_init_session,               &
       IMPLEMENTED3 => user_normalization,              &
       IMPLEMENTED4 => user_set_ics,                    &
       IMPLEMENTED5 => user_set_plot_var,               &
       IMPLEMENTED6 => user_material_properties

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
    use ModPhysics, ONLY: g

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
    StateLowrie_VC(iTgasLowrie,:) = StateLowrie_VC(iTgasLowrie,:)/g
    StateLowrie_VC(iTradLowrie,:) = StateLowrie_VC(iTradLowrie,:)/g

  end subroutine user_init_session

  !==========================================================================

  subroutine user_normalization

    use ModConst,   ONLY: cRadiation, cProtonMass, cBoltzmann
    use ModPhysics, ONLY: g, No2Si_V, UnitRho_, UnitU_
    
    character (len=*), parameter :: NameSub = 'user_normalization'
    !------------------------------------------------------------------------

    No2Si_V = 1.0
    ! The following density unit is needed to get a normalized radiation
    ! constant  with value 1.0e-4. The gamma dependence is needed since
    ! the reference solution uses p=rho*T/gamma
    No2Si_V(UnitRho_) = 1.0e+4*cRadiation*(cProtonMass/cBoltzmann)**4 &
         *No2Si_V(UnitU_)**6/g**4

  end subroutine user_normalization

  !==========================================================================

  subroutine user_set_ics

    use ModAdvance,    ONLY: State_VGB
    use ModConst,      ONLY: cRadiation
    use ModGeometry,   ONLY: x_Blk, y_Blk
    use ModMain,       ONLY: GlobalBlk, nI, nJ, nK, x_, y_
    use ModPhysics,    ONLY: ShockSlope, No2Si_V, Si2No_V, &
         UnitTemperature_, UnitEnergyDens_
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

    do j=-1,nJ+2; do i=-1,nI+2
       
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
       Erad = cRadiation*(Trad*No2Si_V(UnitTemperature_))**4 &
            *Si2No_V(UnitEnergyDens_)
       RhoU_D(1) = Rho*(Ux+U0)
       RhoU_D(2) = 0.0
       RhoU_D(3) = 0.0

       do k=-1,nk+2
          State_VGB(Rho_,i,j,k,iBlock) = Rho
          State_VGB(RhoUx_:RhoUy_,i,j,k,iBlock) = matmul(Rot_II,RhoU_D(x_:y_))
          State_VGB(RhoUz_,i,j,k,iBlock) = 0.0
          State_VGB(Erad_,i,j,k,iBlock) = Erad
          State_VGB(ExtraEint_,i,j,k,iBlock) = 0.0
          State_VGB(p_,i,j,k,iBlock) = p
       end do

    end do; end do

  end subroutine user_set_ics

  !==========================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModAdvance,    ONLY: State_VGB
    use ModConst,      ONLY: cRadiation
    use ModGeometry,   ONLY: x_Blk, y_Blk
    use ModMain,       ONLY: Time_Simulation
    use ModPhysics,    ONLY: Si2No_V, No2Si_V, UnitT_, &
         UnitTemperature_, UnitEnergyDens_, ShockSlope
    use ModSize,       ONLY: nI, nJ, nK
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
    !------------------------------------------------------------------------

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

    IsFound = .true.
    select case(NameVar)
    case('tgas')
       PlotVar_G = State_VGB(p_,:,:,:,iBlock)/State_VGB(Rho_,:,:,:,iBlock)

    case('trad')
       PlotVar_G = (State_VGB(Erad_,:,:,:,iBlock)*No2Si_V(UnitEnergyDens_) &
            /cRadiation)**0.25 *Si2No_V(UnitTemperature_)

    case('rho0','ux0','uy0','tgas0','trad0')

       ! Calculate sin and cos from the tangent = ShockSlope
       SinSlope = ShockSlope/sqrt(1.0+ShockSlope**2)
       CosSlope =        1.0/sqrt(1.0+ShockSlope**2)

       do k=-1,nK+2; do j=-1,nJ+2; do i=-1,nI+2

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

  subroutine user_material_properties(State_V, EinternalSiIn, &
       TeSiOut, AbsorptionOpacitySiOut, RosselandMeanOpacitySiOut, &
       CvSiOut, PressureSiOut)

    ! The State_V vector is in normalized units

    use ModConst,      ONLY: cLightSpeed
    use ModPhysics,    ONLY: g, gm1, inv_gm1, No2Si_V, Si2No_V, &
         UnitTemperature_, UnitEnergyDens_, UnitT_, UnitU_, UnitX_
    use ModVarIndexes, ONLY: nVar, Rho_, p_

    real, intent(in) :: State_V(nVar)
    real, optional, intent(in)  :: EinternalSiIn             ! [J/m^3]
    real, optional, intent(out) :: TeSiOut                   ! [K]
    real, optional, intent(out) :: AbsorptionOpacitySiOut    ! [1/m]
    real, optional, intent(out) :: RosselandMeanOpacitySiOut ! [1/m]
    real, optional, intent(out) :: CvSiOut                   ! [J/(K*m^3)]
    real, optional, intent(out) :: PressureSiOut             ! [Pa]

    real :: Temperature, AbsorptionOpacity, DiffusionRad

    character (len=*), parameter :: NameSub = 'user_material_properties'
    !-------------------------------------------------------------------

    if(present(EinternalSiIn))then
       Temperature = EinternalSiIn*Si2No_V(UnitEnergyDens_) &
            *gm1/State_V(Rho_)
    else
       Temperature = State_V(p_)/State_V(Rho_)
    end if

    if(present(PressureSiOut)) PressureSiOut = EinternalSiIn*(g - 1.0)

    if(present(CvSiOut)) CvSiOut = inv_gm1*State_V(Rho_) &
         *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_)

    select case(iLowrieTest)
    case(1,2)
       DiffusionRad = 1.0
       AbsorptionOpacity = 1.0E6
    case(3)
       DiffusionRad = 0.00175*(g*Temperature)**3.5/State_V(Rho_)
       AbsorptionOpacity = 1.0E6/DiffusionRad
    end select

    if(present(TeSiOut)) TeSiOut = Temperature*No2Si_V(UnitTemperature_)

    if(present(AbsorptionOpacitySiOut)) AbsorptionOpacitySiOut = &
         AbsorptionOpacity/No2Si_V(UnitT_)/cLightSpeed

    if(present(RosselandMeanOpacitySiOut)) RosselandMeanOpacitySiOut = &
         cLightSpeed/(3.0*DiffusionRad*No2Si_V(UnitU_)*No2Si_V(UnitX_))

  end subroutine user_material_properties

end module ModUser
