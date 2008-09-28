!^CFG COPYRIGHT UM
!========================================================================
module ModUser
  ! This is the default user module which contains empty methods defined
  ! in ModUserEmpty.f90

  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_update_states,              &
       IMPLEMENTED2 => user_calc_sources,               &
       IMPLEMENTED3 => user_initial_perturbation,       &
       IMPLEMENTED4 => user_read_inputs,                &
       IMPLEMENTED5 => user_set_plot_var

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'HYDRO + IONIZATION EQUILIBRIUM + LEVEL SETS'

  ! Some default values
  real :: xBe = 1000.0, yBe = 200., RhoWallDim = 50.0

contains

  subroutine user_read_inputs

    use ModReadParam
    character (len=100) :: NameCommand
    !------------------------------------------------------------------------
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#MATERIAL")
          call read_var('xBe', xBe)
          call read_var('yBe', yBe)
          call read_var('RhoWallDim', RhoWallDim)
       case('#USERINPUTEND')
          EXIT
       case default
          call stop_mpi('ERROR in ModUserCrash: unknown command='//NameCommand)
       end select
    end do

  end subroutine user_read_inputs

  !============================================================================
  subroutine user_update_states(iStage,iBlock)

    use ModVarIndexes
    use ModSize
    use ModAdvance, ONLY: State_VGB
    use ModMain,    ONLY: nStage
    use ModPhysics
    use ModEnergy,  ONLY: calc_energy_cell
    use ModEos

    implicit none

    integer, intent(in):: iStage,iBlock

    integer:: i, j, k, iMaterial, Loc_I(1)
    real   :: PressureSI, EInternal, EInternalSI, RhoSI
    !------------------------------------------------------------------------

    call update_states_MHD(iStage,iBlock)

    ! Do not use previous solution in Newton iteration to avoid dependency
    ! on number and order of blocks on a processor
    UsePreviousTe = .false.

    ! update of pressure and relaxation energy::
    do k=1,nK; do j=1,nJ; do i=1,nI
       ! Total internal energy ExtraEInt + P/(\gamma -1) transformed to SI
       EInternalSI = No2Si_V(UnitEnergyDens_)*&
            (inv_gm1*State_VGB(P_,i,j,k,iBlock) + &
            State_VGB(ExtraEInt_,i,j,k,iBlock))

       ! Density, transformed to SI
       RhoSI = No2Si_V(UnitRho_)*State_VGB(Rho_,i,j,k,iBlock)

       ! Find maximum level set value. 
       ! Note that for now plastic is now taken as carbon
       Loc_I = maxloc(State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock))
       iMaterial = Loc_I(1) - 1

       ! Apply the EOS, get pressure in SI
       call eos(&
            UDensityTotal=EInternalSI,& !Input total energy density SI,[J/m^3]
            Rho=RhoSI,                & !Input mass density, SI [kg/m^3] 
            iMaterial=iMaterial,      & !Input: sort of material
            PTotalOut=PressureSI      ) !Output, OPTIONAL, pressure, SI [Pa]

       ! Set pressure and ExtraEInt = Total internal energy - P/(gamma -1)
       State_VGB(P_,i,j,k,iBlock) = PressureSI*Si2No_V(UnitP_)
       State_VGB(ExtraEInt_,i,j,k,iBlock) = Si2No_V(UnitEnergyDens_)*&
            (EInternalSI - PressureSI*inv_gm1)

    end do; end do; end do

    call calc_energy_cell(iBlock)

  end subroutine user_update_states

  !===========================================================================

  subroutine user_calc_sources

    use ModMain,     ONLY: nI, nJ, nK, GlobalBlk
    use ModAdvance,  ONLY: State_VGB, LevelXe_, LevelPl_, &
         Source_VC, uDotArea_XI, uDotArea_YI, uDotArea_ZI
    use ModGeometry, ONLY: vInv_CB

    character (len=*), parameter :: NameSub = 'user_calc_sources'
    integer :: i, j, k, iBlock
    !-------------------------------------------------------------------

    iBlock = globalBlk

    ! Add Level*div(u) as a source term so level sets beome advected scalars
    ! Note that all levels use the velocity of the first (and only) fluid
    
    do k=1,nK; do j=1,nJ; do i=1,nI
       Source_VC(LevelXe_:LevelPl_,i,j,k) =                 &
            Source_VC(LevelXe_:LevelPl_,i,j,k)              &
            + State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock)     &
            * vInv_CB(i,j,k,iBlock)*                        &
            ( uDotArea_XI(i+1,j,k,1) - uDotArea_XI(i,j,k,1) &
            + uDotArea_YI(i,j+1,k,1) - uDotArea_YI(i,j,k,1) &
            + uDotArea_ZI(i,j,k+1,1) - uDotArea_ZI(i,j,k,1))
    end do; end do; end do


  end subroutine user_calc_sources

  !===========================================================================

  subroutine user_initial_perturbation

    use ModMain,    ONLY: nI, nJ, nK, nBlock, UnusedBlk, UseUserSource
    use ModPhysics, ONLY: ShockPosition, ShockSlope, ShockRightState_V, &
         Io2No_V, UnitRho_, UnitP_
    use ModAdvance, ONLY: State_VGB, Rho_, RhoUx_, RhoUz_, p_, &
         LevelBe_, LevelXe_, LevelPl_
    use ModGeometry,ONLY: x_BLK, y_BLK, y2
    use ModEnergy,  ONLY: calc_energy_ghost

    integer :: i, j, k, iBlock
    real    :: x, y, xBeSlope

    character (len=*), parameter :: NameSub = 'user_initial_perturbation'
    !-------------------------------------------------------------------
  
    do iBlock = 1, nBlock

       if (unusedBLK(iBlock)) CYCLE

       do k=-1, nK+2; do j=-1, nJ+2; do i=-1, nI+2 

          x = x_BLK(i,j,k,iBlock)
          y = y_BLK(i,j,k,iBlock)

          ! Be - Xe interface is at xBe + slope due to shockslope
          xBeSlope = xBe - ShockSlope*y

          ! Set plastic wall state
          if(abs(y) > yBe)then
             ! Use the density given by the #MATERIAL command
             State_VGB(Rho_,i,j,k,iBlock) = RhoWallDim*Io2No_V(UnitRho_)
             ! Assume pressure equilibrium with Xenon
             State_VGB(p_,i,j,k,iBlock) = ShockRightState_V(p_)*Io2No_V(UnitP_)
             ! Assume that plastic wall is at rest
             State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
          end if

          ! Plastic is present outside of yBe
          State_VGB(LevelPl_,i,j,k,iBlock) = abs(y) - yBe

          ! Berilium is inside plastic wall and left to xBe
          State_VGB(LevelBe_,i,j,k,iBlock) = min(xBeSlope - x, yBe - abs(y))

          ! Xenon is inside plastic wall and right to xBe
          State_VGB(LevelXe_,i,j,k,iBlock) = min(x - xBeSlope, yBe - abs(y))

          if(.not.UseUserSource) &
               State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock) = &
               State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock) &
               *State_VGB(Rho_,i,j,k,iBlock)

       end do; end do; end do

       call calc_energy_ghost(iBlock)

    end do

  end subroutine user_initial_perturbation

  !====================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModSize,    ONLY: nI, nJ, nK
    use ModAdvance, ONLY: State_VGB, LevelXe_, LevelPl_

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

    character (len=*), parameter :: Name='user_set_plot_var'

    integer :: i, j, k, Loc_I(1)
    !------------------------------------------------------------------------  
    IsFound = NameVar == 'level'
    if(.not. IsFound) RETURN

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0
    
    do k=-1, nK+1; do j=-1, nJ+1; do i=-1,nI+2
       Loc_I = maxloc(State_VGB(LevelXe_:LevelPl_,i,j,k,iBlock))
       PlotVar_G(i,j,k) = Loc_I(1)
    end do; end do; end do

  end subroutine user_set_plot_var

end module ModUser
