!^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_init_session,               &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_set_outerbcs,               &
       IMPLEMENTED5 => user_update_states,              &
       IMPLEMENTED6 => user_set_plot_var,               &
       IMPLEMENTED7 => user_material_properties,        &
       IMPLEMENTED8 => user_normalization

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'heat conduction'

  character(len=20) :: TypeProblem
  real :: HeatConductionCoef, AmplitudeTemperature, Tmin
  real :: Bx, By, Time0
  real :: u0, x0

  integer :: &
       nCellRef = -1, &
       nVarRef  = -1, &
       iRhoRef  = -1, &
       iUxRef   = -1, &
       iUrRef   = -1, &
       iTeRef   = -1, &
       iTiRef   = -1
  real, allocatable :: rRef_C(:), StateRef_VC(:,:)

  integer :: nCellFin
  integer :: &
       iRhoFin = 1, &
       iTeFin  = 2, &
       iUrFin  = 3, &
       nVarFin = 3
  real, allocatable :: rFin_C(:), StateFin_VC(:,:)

  real, parameter :: GammaRel = 4.0/3.0

contains

  !============================================================================

  subroutine user_init_session

    use ModIoUnit,  ONLY: UnitTmp_
    use ModMain,    ONLY: Time_Simulation
    use ModPhysics, ONLY: g
    use ModProcMH,  ONLY: iProc

    integer :: iCell, iError
    integer :: nStepRef, nDimRef, nParamRef, nVarRef
    real :: TimeRef, GammaRef
    character(len=100) :: NameRefFile
    character(len=500) :: StringHeaderRef, NameVarRef
    integer :: nStepFin, nDimFin, nParamFin, nVarFin
    real :: TimeFin, GammaFin
    character(len=100) :: NameFinFile
    character(len=500) :: StringHeaderFin, NameVarFin

    character(len=*), parameter :: NameSub = 'user_init_session'
    !--------------------------------------------------------------------------

    if(TypeProblem=='gaussian' .or. TypeProblem=='rz' &
         .or. TypeProblem=='parcond' .or. TypeProblem=='parcondsemi')then
       if(Time_Simulation <= 0.0)then
          if(iProc == 0) write(*,*) NameSub// &
               ' : starting simulation time should be larger than 0'
          call stop_mpi('reset time with #TIMESIMULATION')
       end if
    end if

    select case(TypeProblem)
    case('gaussian')
       HeatConductionCoef = 0.1
       AmplitudeTemperature = 10.0
       Tmin = 10.0

    case('rz')
       HeatConductionCoef = 0.1
       AmplitudeTemperature = 10.0
       Tmin = 3.0

    case('rmtv')
       Tmin = 1.0e-12   ! minimum temperature at far field

       NameRefFile = 'rmtv_initial.out'

       open(UnitTmp_, FILE=NameRefFile, STATUS="old", IOSTAT=iError)

       if(iError /= 0)call stop_mpi(NameSub // &
            " could not open reference file="//NameRefFile)

       read(UnitTmp_,'(a)') StringHeaderRef
       read(UnitTmp_,*) nStepRef, TimeRef, nDimRef, nParamRef, nVarRef
       read(UnitTmp_,*) nCellRef
       read(UnitTmp_,*) GammaRef
       read(UnitTmp_,'(a)') NameVarRef

       iRhoRef = 1; iTeRef = 2; iUrRef = 3
       nVarRef = 3

       allocate( rRef_C(nCellRef), StateRef_VC(nVarRef,nCellRef) )

       do iCell = 1, nCellRef
          read(UnitTmp_,*) rRef_C(iCell), StateRef_VC(:,iCell)
       end do

       close(UnitTmp_)

       do iCell = 1, nCellRef
          StateRef_VC(iTeRef,iCell) = &
               max(StateRef_VC(iTeRef,iCell),Tmin)
       end do

       NameFinFile = 'rmtv_final.out'

       open(UnitTmp_, FILE=NameFinFile, STATUS="old", IOSTAT=iError)

       if(iError /= 0)call stop_mpi(NameSub // &
            " could not open reference file="//NameFinFile)

       read(UnitTmp_,'(a)') StringHeaderFin
       read(UnitTmp_,*) nStepFin, TimeFin, nDimFin, nParamFin, nVarFin
       read(UnitTmp_,*) nCellFin
       read(UnitTmp_,*) GammaFin
       read(UnitTmp_,'(a)') NameVarFin

       allocate( rFin_C(nCellFin), StateFin_VC(nVarFin,nCellFin) )

       do iCell = 1, nCellFin
          read(UnitTmp_,*) rFin_C(iCell), StateFin_VC(:,iCell)
       end do

       close(UnitTmp_)

       do iCell = 1, nCellFin
          StateFin_VC(iTeFin,iCell) = &
               max(StateFin_VC(iTeFin,iCell),Tmin)
       end do

    case('lowrie')
       ! Mach 5 test with non-linear heat conduction and electron-ion
       ! interaction rate. This test is based one of Lowrie's non-equilibrium
       ! gray radiation diffusion tests.
       NameRefFile = 'lowrie3.out'

       open(UnitTmp_, FILE=NameRefFile, STATUS="old", IOSTAT=iError)

       if(iError /= 0)call stop_mpi(NameSub // &
            " could not open reference file="//NameRefFile)

       read(UnitTmp_,'(a)') StringHeaderRef
       read(UnitTmp_,*) nStepRef, TimeRef, nDimRef, nParamRef, nVarRef
       read(UnitTmp_,*) nCellRef
       read(UnitTmp_,*) GammaRef
       read(UnitTmp_,'(a)') NameVarRef

       iRhoRef = 1; iUxRef = 2; iTiRef = 3; iTeRef = 4
       nVarRef = 4

       allocate( rRef_C(nCellRef), StateRef_VC(nVarRef,nCellRef) )

       do iCell = 1, nCellRef
          read(UnitTmp_,*) rRef_C(iCell), StateRef_VC(:,iCell)
       end do

       close(UnitTmp_)

       ! The reference solutions use p=rho*T/gamma
       StateRef_VC(iTiRef,:) = StateRef_VC(iTiRef,:)/g
       StateRef_VC(iTeRef,:) = StateRef_VC(iTeRef,:)/g

       u0 = -5.0  ! veloxity added to lowrie's solution
       x0 = 0.02  ! shift in x-direction

    case('parcond')
       HeatConductionCoef = 1.0
       AmplitudeTemperature = 9.0
       Tmin = 1.0
       Bx = 1.7
       By = 1.0
       Time0 = Time_Simulation

    case('parcondsemi')
       AmplitudeTemperature = 100.0
       Tmin = 0.01
       Bx = 1.7
       By = 1.0
       Time0 = Time_Simulation

    case default
       call stop_mpi(NameSub//' : undefined problem type='//TypeProblem)
    end select

  end subroutine user_init_session

  !============================================================================

  subroutine user_normalization

    use ModConst,   ONLY: cRadiation, cProtonMass, cBoltzmann
    use ModPhysics, ONLY: No2Si_V, UnitRho_, UnitU_, g

    character (len=*), parameter :: NameSub = 'user_normalization'
    !--------------------------------------------------------------------------

    if(TypeProblem == 'lowrie')then
       No2Si_V = 1.0
       ! The following density unit is needed to get a normalized radiation
       ! constant  with value 1.0e-4. The gamma dependence is needed since
       ! the reference solution uses p=rho*T/gamma
       No2Si_V(UnitRho_) = 1.0e+4*cRadiation*(cProtonMass/cBoltzmann)**4 &
            *No2Si_V(UnitU_)**6/g**4
    end if

  end subroutine user_normalization

  !============================================================================

  subroutine user_update_states(iStage, iBlock)

    use ModAdvance,  ONLY: nVar, Flux_VX, Flux_VY, Flux_VZ, Source_VC, &
         UseElectronPressure
    use ModImplicit, ONLY: UseSemiImplicit
    use ModVarIndexes, ONLY: Pe_

    integer, intent(in) :: iStage, iBlock

    character(len=*), parameter :: NameSub = 'user_update_states'
    !--------------------------------------------------------------------------

    ! No call to update_states_MHD to nullify the effect of the hydro solver
    ! call update_states_MHD(iStage,iBlock)
    if(TypeProblem == 'parcond' .and. .not.UseSemiImplicit)then
       if(UseElectronPressure)then
          Flux_VX(1:Pe_-1,:,:,:) = 0.0; Flux_VX(Pe_+1:nVar+1,:,:,:) = 0.0
          Flux_VY(1:Pe_-1,:,:,:) = 0.0; Flux_VY(Pe_+1:nVar+1,:,:,:) = 0.0
          Flux_VZ(1:Pe_-1,:,:,:) = 0.0; Flux_VZ(Pe_+1:nVar+1,:,:,:) = 0.0
          Source_VC(1:Pe_-1,:,:,:) = 0.0; Source_VC(Pe_+1:nVar+1,:,:,:) = 0.0
       else
          Flux_VX(1:nVar,:,:,:) = 0.0
          Flux_VY(1:nVar,:,:,:) = 0.0
          Flux_VZ(1:nVar,:,:,:) = 0.0
          Source_VC(1:nVar,:,:,:) = 0.0
       end if
       call update_states_MHD(iStage,iBlock)
    end if

    if(TypeProblem == 'lowrie') call update_states_electron

  contains

    subroutine update_states_electron

      use ModAdvance,    ONLY: State_VGB
      use ModMain,       ONLY: nI, nJ, nK
      use ModPhysics,    ONLY: inv_gm1, Si2No_V, No2Si_V, UnitEnergyDens_, &
           UnitP_
      use ModVarIndexes, ONLY: Pe_, ExtraEint_

      integer :: i, j, k
      real :: PeSi, Ee, EeSi
      !------------------------------------------------------------------------

      call update_states_MHD(iStage,iBlock)

      do k = 1, nK; do j = 1, nJ; do i = 1, nI
         ! At this point Pe=(g-1)*Ee with the ideal gamma g.
         ! Use this Pe to get electron internal energy density.

         Ee = inv_gm1*State_VGB(Pe_,i,j,k,iBlock) &
              + State_VGB(ExtraEint_,i,j,k,iBlock)
         EeSi = Ee*No2Si_V(UnitEnergyDens_)

         call user_material_properties(State_VGB(:,i,j,k,iBlock), &
              i, j, k, iBlock, &
              EinternalIn=EeSi, PressureOut=PeSi)

         ! Set true electron pressure
         State_VGB(Pe_,i,j,k,iBlock) = PeSi*Si2No_V(UnitP_)

         ! Set ExtraEint = electron internal energy - Pe/(gamma -1)
         State_VGB(ExtraEint_,i,j,k,iBlock) = &
              Ee - inv_gm1*State_VGB(Pe_,i,j,k,iBlock)

         if(State_VGB(ExtraEint_,i,j,k,iBlock)<0.0)then
            write(*,*)NameSub,': ERROR extra internal energy =', &
                 State_VGB(ExtraEint_,i,j,k,iBlock)
            write(*,*)NameSub,': ERROR at i,j,k,iBlock=', i, j, k, iBlock
            call stop_mpi(NameSub//': ERROR negative extra internal energy')
         end if

      end do; end do; end do

    end subroutine update_states_electron

  end subroutine user_update_states

  !============================================================================

  subroutine user_read_inputs

    use ModIO,          ONLY: write_prefix, write_myname, iUnitOut
    use ModMain,        ONLY: lVerbose
    use ModProcMH,      ONLY: iProc
    use ModReadParam,   ONLY: read_line, read_command, read_var

    character (len=100) :: NameCommand
    !--------------------------------------------------------------------------
    if(iProc == 0 .and. lVerbose > 0)then
       call write_prefix;
       write(iUnitOut,*)'User read_input starts'
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
       case("#HEATCONDUCTIONTEST")
          call read_var('TypeProblem',TypeProblem)

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

  subroutine user_set_ics

    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: x_Blk, y_Blk
    use ModMain,       ONLY: GlobalBlk, nI, nJ, nK, Time_Simulation, x_, y_
    use ModPhysics,    ONLY: ShockSlope, cRadiationNo, inv_gm1
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUy_, RhoUz_, p_, ExtraEint_, Pe_

    integer :: iBlock, i, j, k, iCell
    real :: x, y
    real :: r, Weight1, Weight2, Rho, Te, Ur, p, RhoU_D(2)
    real :: SinSlope, CosSlope, Rot_II(2,2), Ee, Pe, Ti, Ux

    character(len=*), parameter :: NameSub = "user_set_ics"
    !--------------------------------------------------------------------------
    iBlock = GlobalBlk

    select case(TypeProblem)
    case('gaussian')
       do i = -1, nI+2
          call get_state_gaussian(i, iBlock)
       end do

    case('rz')
       do j=-1,nJ+2; do i=-1,nI+2
          call get_state_rz(i, j, iBlock)
       end do; end do

    case('rmtv')
       do j=-1,nJ+2; do i=-1,nI+2
          x = x_Blk(i,j,0,iBlock)
          y = y_Blk(i,j,0,iBlock)
          r = sqrt(x**2+y**2)

          do iCell = 1, nCellRef
             if(rRef_C(iCell) >= r) EXIT
          end do
          if(iCell == 1) call stop_mpi(NameSub // &
               " Reference solution does not cover the left boundary")

          if(iCell > nCellRef)then
             ! Cell is beyond the last point of Reference input: use last cell

             Rho = 1.0/r**(19.0/9.0)
             Te = StateRef_VC(iTeRef,nCellRef)
             Ur = StateRef_VC(iUrRef,nCellRef)
          else
             ! Assign weights for linear interpolation between
             ! iCell-1 and iCell
             Weight1 = (rRef_C(iCell) - r) &
                  /    (rRef_C(iCell) - rRef_C(iCell-1))
             Weight2 = 1.0 - Weight1

             Rho  = ( Weight1*StateRef_VC(iRhoRef, iCell-1) &
                  +   Weight2*StateRef_VC(iRhoRef, iCell) )
             Te =   ( Weight1*StateRef_VC(iTeRef, iCell-1) &
                  +   Weight2*StateRef_VC(iTeRef, iCell) )
             Ur   = ( Weight1*StateRef_VC(iUrRef, iCell-1) &
                  +   Weight2*StateRef_VC(iUrRef, iCell) )
          end if

          RhoU_D(1) = Rho*Ur*x/r
          RhoU_D(2) = Rho*Ur*y/r

          p = Rho*Te

          do k=-1,nK+2
             State_VGB(Rho_,i,j,k,iBlock) = Rho
             State_VGB(RhoUx_:RhoUy_,i,j,k,iBlock) = RhoU_D
             State_VGB(RhoUz_,i,j,k,iBlock) = 0.0
             State_VGB(ExtraEint_,i,j,k,iBlock) = 0.0
             State_VGB(p_,i,j,k,iBlock) = p
          end do
       end do; end do

    case('lowrie')
       ! Calculate sin and cos from the tangent = ShockSlope
       SinSlope = ShockSlope/sqrt(1.0+ShockSlope**2)
       CosSlope =        1.0/sqrt(1.0+ShockSlope**2)
       ! Set rotational matrix
       Rot_II = reshape( (/CosSlope, SinSlope, -SinSlope, CosSlope/), (/2,2/) )

       do j = -1, nJ + 2; do i = -1, nI+2

          x = x_Blk(i,j,0,iBlock)*CosSlope + y_Blk(i,j,0,iBlock)*SinSlope - x0

          do iCell = 1, nCellRef
             if(rRef_C(iCell) >= x) EXIT
          end do
          if(iCell == 1) call stop_mpi(NameSub // &
               " Reference solution does not cover the left boundary")

          if(iCell > nCellRef)then
             ! Cell is beyond the last point of reference input: use last cell
             iCell   = nCellRef
             Weight1 = 0.0
             Weight2 = 1.0
          else
             ! Assign weights for linear interpolation between
             ! iCell-1 and iCell
             Weight1 = (rRef_C(iCell) - x) &
                  /    (rRef_C(iCell) - rRef_C(iCell-1))
             Weight2 = 1.0 - Weight1
          end if

          Rho = ( Weight1*StateRef_VC(iRhoRef, iCell-1) &
               +  Weight2*StateRef_VC(iRhoRef, iCell) )
          Ux  = ( Weight1*StateRef_VC(iUxRef, iCell-1) &
               +  Weight2*StateRef_VC(iUxRef, iCell) )
          Ti  = ( Weight1*StateRef_VC(iTiRef, iCell-1) &
               +  Weight2*StateRef_VC(iTiRef, iCell) )
          Te  = ( Weight1*StateRef_VC(iTeRef, iCell-1) &
               +  Weight2*StateRef_VC(iTeRef, iCell) )

          p = Rho*Ti
          Ee = cRadiationNo*Te**4
          Pe = (GammaRel - 1)*Ee
          RhoU_D(1) = Rho*(Ux + U0)
          RhoU_D(2) = 0.0

          do k = -1, nk+2
             State_VGB(Rho_,i,j,k,iBlock) = Rho
             State_VGB(RhoUx_:RhoUy_,i,j,k,iBlock) = &
                  matmul(Rot_II,RhoU_D(x_:y_))
             State_VGB(RhoUz_,i,j,k,iBlock) = 0.0
             State_VGB(Pe_,i,j,k,iBlock) = Pe
             State_VGB(ExtraEint_,i,j,k,iBlock) = &
                  Ee - inv_gm1*State_VGB(Pe_,i,j,k,iBlock)
             State_VGB(p_,i,j,k,iBlock) = p
          end do

       end do; end do

    case('parcond', 'parcondsemi')
       do j=-1,nJ+2; do i=-1,nI+2
          call get_state_parcond(i, j, iBlock)
       end do; end do

    case default
       call stop_mpi(NameSub//' : undefined problem type='//TypeProblem)
    end select

  end subroutine user_set_ics

  !============================================================================

  subroutine user_set_outerbcs(iBlock, iSide, TypeBc, IsFound)

    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: x_Blk, y_Blk
    use ModImplicit,   ONLY: StateSemi_VGB
    use ModMain,       ONLY: nI, nJ
    use ModVarIndexes, ONLY: Rho_, RhoUx_, p_

    integer,          intent(in)  :: iBlock, iSide
    character(len=20),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer :: i, j
    real :: r, Temperature

    character (len=*), parameter :: NameSub = 'user_set_outerbcs'
    !--------------------------------------------------------------------------

    select case(TypeProblem)
    case('gaussian')
       select case(TypeBc)
       case('user')
          select case(iSide)
          case(1)
             do i = -1, 0
                call get_state_gaussian(i, iBlock)
             end do
          case(2)
             do i = nI+1, nI+2
                call get_state_gaussian(i, iBlock)
             end do
          end select
       case('usersemi')
          select case(iSide)
          case(1)
             call get_temperature_gaussian(0, iBlock, Temperature)
             StateSemi_VGB(1,0,:,:,iBlock) = Temperature
          case(2)
             call get_temperature_gaussian(nI+1, iBlock, Temperature)
             StateSemi_VGB(1,nI+1,:,:,iBlock) = Temperature
          end select
       end select

    case('rz')
       select case(TypeBc)
       case('user')
          select case(iSide)
          case(1)
             do j = -1, nJ+2; do i = -1, 0
                call get_state_rz(i, j, iBlock)
             end do; end do
          case(2)
             do j = -1, nJ+2; do i = nI+1, nI+2
                call get_state_rz(i, j, iBlock)
             end do; end do
          end select
       case('usersemi')
          select case(iSide)
          case(1)
             do j = 1, nJ
                call get_temperature_rz(0, j, iBlock, Temperature)
                StateSemi_VGB(1,0,j,:,iBlock) = Temperature
             end do
          case(2)
             do j = 1, nJ
                call get_temperature_rz(nI+1, j, iBlock, Temperature)
                StateSemi_VGB(1,nI+1,j,:,iBlock) = Temperature
             end do
          end select
       end select

    case('rmtv')
       if(.not. (iSide==2 .or. iSide==4) )then
          write(*,*) NameSub//' : user boundary not defined at iSide = ', iSide
          call stop_mpi(NameSub)
       end if

       select case(TypeBc)
       case('user')
          select case(iSide)
          case(2) ! z-direction in rz-geometry
             do j = -1, nJ+2; do i = nI+1, nI+2
                r = sqrt(x_Blk(i,j,1,iBlock)**2+y_Blk(i,j,1,iBlock)**2)
                State_VGB(Rho_,i,j,:,iBlock) = 1.0/r**(19.0/9.0)
                State_VGB(RhoUx_:p_-1,i,j,:,iBlock) = &
                     State_VGB(RhoUx_:p_-1,nI,j,:,iBlock)
                State_VGB(p_,i,j,:,iBlock) = State_VGB(Rho_,i,j,:,iBlock)*Tmin
             end do; end do
          case(4) ! r-direction in rz-geometry
             do j = nJ+1, nJ+2; do i = -1, nI+2
                r = sqrt(x_Blk(i,j,1,iBlock)**2+y_Blk(i,j,1,iBlock)**2)
                State_VGB(Rho_,i,j,:,iBlock) = 1.0/r**(19.0/9.0)
                State_VGB(RhoUx_:p_-1,i,j,:,iBlock) = &
                     State_VGB(RhoUx_:p_-1,i,nJ,:,iBlock)
                State_VGB(p_,i,j,:,iBlock) = State_VGB(Rho_,i,j,:,iBlock)*Tmin
             end do; end do
          end select
       case('usersemi')
          select case(iSide)
          case(2)
             StateSemi_VGB(1,nI+1,:,:,iBlock) = Tmin
          case(4)
             StateSemi_VGB(1,:,nJ+1,:,iBlock) = Tmin
          end select
       end select

    case('parcond', 'parcondsemi')
       select case(TypeBc)
       case('user')
          select case(iSide)
          case(1)
             do j = -1, nJ+2; do i = -1, 0
                call get_state_parcond(i, j, iBlock)
             end do; end do
          case(2)
             do j = -1, nJ+2; do i = nI+1, nI+2
                call get_state_parcond(i, j, iBlock)
             end do; end do
          case(3)
             do j = -1, 0; do i = -1, nI+2
                call get_state_parcond(i, j, iBlock)
             end do; end do
          case(4)
             do j = nJ+1, nJ+2; do i = -1, nI+2
                call get_state_parcond(i, j, iBlock)
             end do; end do
          end select
       case('usersemi')
          select case(iSide)
          case(1)
             do j = 0, nJ+1
                call get_temperature_parcond(0, j, iBlock, Temperature)
                StateSemi_VGB(1,0,j,:,iBlock) = Temperature
             end do
          case(2)
             do j = 0, nJ+1
                call get_temperature_parcond(nI+1, j, iBlock, Temperature)
                StateSemi_VGB(1,nI+1,j,:,iBlock) = Temperature
             end do
          case(3)
             do i = 0, nI+1
                call get_temperature_parcond(i, 0, iBlock, Temperature)
                StateSemi_VGB(1,i,0,:,iBlock) = Temperature
             end do
          case(4)
             do i = 0, nI+1
                call get_temperature_parcond(i, nJ+1, iBlock, Temperature)
                StateSemi_VGB(1,i,nI+1,:,iBlock) = Temperature
             end do
          end select
       end select

    case default
       call stop_mpi(NameSub//' : undefined problem type='//TypeProblem)
    end select

    IsFound = .true.

  end subroutine user_set_outerbcs

  !============================================================================

  subroutine get_state_gaussian(i, iBlock)

    use ModAdvance,    ONLY: State_VGB
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, p_

    integer, intent(in) :: i, iBlock

    real :: Temperature
    !--------------------------------------------------------------------------
    call get_temperature_gaussian(i, iBlock, Temperature)
    State_VGB(Rho_,i,:,:,iBlock) = 1.0
    State_VGB(RhoUx_:RhoUz_,i,:,:,iBlock) = 0.0
    State_VGB(p_,i,:,:,iBlock) = Temperature

  end subroutine get_state_gaussian

  !============================================================================

  subroutine get_temperature_gaussian(i, iBlock, Temperature)

    use ModGeometry, ONLY: x_Blk
    use ModMain,     ONLY: Time_Simulation
    use ModNumConst, ONLY: cPi

    integer, intent(in) :: i, iBlock
    real,    intent(out):: Temperature

    real :: Spread
    !--------------------------------------------------------------------------

    Spread = 4.0*HeatConductionCoef*Time_Simulation
    Temperature = Tmin + AmplitudeTemperature/(sqrt(cPi*Spread)) &
         *exp(-x_Blk(i,1,1,iBlock)**2/Spread)

  end subroutine get_temperature_gaussian

  !============================================================================

  subroutine get_state_rz(i, j, iBlock)

    use ModAdvance,    ONLY: State_VGB
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, p_

    integer, intent(in) :: i, j, iBlock

    real :: Temperature
    !--------------------------------------------------------------------------
    call get_temperature_rz(i, j, iBlock, Temperature)
    State_VGB(Rho_,i,j,:,iBlock) = 1.0
    State_VGB(RhoUx_:RhoUz_,i,j,:,iBlock) = 0.0
    State_VGB(p_,i,j,:,iBlock) = Temperature

  end subroutine get_state_rz

  !============================================================================

  subroutine get_temperature_rz(i, j, iBlock, Temperature)

    use ModGeometry, ONLY: x_Blk, y_Blk, y2
    use ModMain,     ONLY: Time_Simulation
    use ModNumConst, ONLY: cPi

    integer, intent(in) :: i, j, iBlock
    real,    intent(out):: Temperature

    real :: r, Lambda, Spread
    !--------------------------------------------------------------------------

    Lambda = -(3.831705970/y2)**2
    Spread = 4.0*HeatConductionCoef*Time_Simulation
    r = abs(y_BLK(i,j,1,iBlock))
    Temperature = Tmin + AmplitudeTemperature &
         *exp(Lambda*HeatConductionCoef*Time_Simulation) &
         *bessj0(sqrt(-Lambda)*r) &
         /(sqrt(cPi*Spread))*exp(-x_Blk(i,j,1,iBlock)**2/Spread)

  end subroutine get_temperature_rz

  !============================================================================

  subroutine get_state_parcond(i, j, iBlock)

    use ModAdvance,    ONLY: State_VGB, UseElectronPressure
    use ModPhysics,    ONLY: inv_gm1, ElectronTemperatureRatio
    use ModVarIndexes, ONLY: Rho_, Bx_, By_, p_, Pe_, ExtraEint_

    integer, intent(in) :: i, j, iBlock

    real :: Te
    !--------------------------------------------------------------------------
    call get_temperature_parcond(i, j, iBlock, Te)
    State_VGB(:,i,j,:,iBlock) = 0.0
    State_VGB(Rho_,i,j,:,iBlock) = 1.0
    State_VGB(Bx_,i,j,:,iBlock) = Bx
    State_VGB(By_,i,j,:,iBlock) = By
    if(TypeProblem == 'parcondsemi')then
       State_VGB(p_,i,j,:,iBlock) = Te
       State_VGB(ExtraEint_,i,j,:,iBlock) = (1.0/3.5)*Te**3.5 - inv_gm1*Te
    else
       if(UseElectronPressure)then
          State_VGB(Pe_,i,j,:,iBlock) = Te
       else
          State_VGB(p_,i,j,:,iBlock) = Te*(1 + ElectronTemperatureRatio)
       end if
    end if

  end subroutine get_state_parcond

  !============================================================================

  subroutine get_temperature_parcond(i, j, iBlock, Temperature)

    use ModGeometry, ONLY: x_Blk, y_Blk
    use ModMain,     ONLY: Time_Simulation
    use ModNumConst, ONLY: cPi

    integer, intent(in) :: i, j, iBlock
    real,    intent(out):: Temperature

    real :: x, y, xx, yy, Spread, Spread0
    !--------------------------------------------------------------------------

    x = x_Blk(i,j,0,iBlock)
    y = y_Blk(i,j,0,iBlock)
    xx = (Bx*x+By*y)/sqrt(Bx**2+By**2)
    yy = (Bx*y-By*x)/sqrt(Bx**2+By**2)

    Spread0 = 4.0*Time0
    Spread = 4.0*Time_Simulation

    Temperature = Tmin + AmplitudeTemperature/(sqrt(cPi*Spread)) &
         *exp(-xx**2/Spread-yy**2/Spread0)

    if(TypeProblem == 'parcondsemi') Temperature = Temperature**(1.0/3.5)

  end subroutine get_temperature_parcond

  !============================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModAdvance,    ONLY: State_VGB, UseIdealEos, UseElectronPressure
    use ModGeometry,   ONLY: x_Blk, y_Blk
    use ModMain,       ONLY: nI, nJ, nK, Time_Simulation
    use ModPhysics,    ONLY: Si2No_V, UnitTemperature_, UnitT_, ShockSlope, &
         ElectronTemperatureRatio
    use ModVarIndexes, ONLY: p_, Rho_, Pe_

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

    real :: r, Temperature, Weight1, Weight2, x, y
    real :: Rho, Ur, U_D(2)
    real :: SinSlope, CosSlope, Ux, Te, Ti, TeSi
    integer :: i, j, k, iCell

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    !--------------------------------------------------------------------------

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

    IsFound = .true.

    select case(TypeProblem)
    case('gaussian')
       select case(NameVar)
       case('t0','temp0')
          do i=-1,nI+2
             call get_temperature_gaussian(i, iBlock, Temperature)
             PlotVar_G(i,:,:) = Temperature
          end do
       case default
          IsFound = .false.
       end select

    case('rz')
       select case(NameVar)
       case('t0','temp0')
          do j=-1,nJ+2; do i=-1,nI+2
             call get_temperature_rz(i, j, iBlock, Temperature)
             PlotVar_G(i,j,:) = Temperature
          end do; end do
       case default
          IsFound = .false.
       end select

    case('rmtv')
       do j=-1,nJ+2; do i=-1,nI+2
          x = x_Blk(i,j,0,iBlock)
          y = y_Blk(i,j,0,iBlock)
          r = sqrt(x**2+y**2)

          do iCell = 1, nCellFin
             if(rFin_C(iCell) >= r) EXIT
          end do
          if(iCell == 1) call stop_mpi(NameSub // &
               " Reference solution does not cover the left boundary")

          if(iCell > nCellFin)then
             ! Cell is beyond the last point of Reference input: use last cell

             Rho = 1.0/r**(19.0/9.0)
             Temperature = StateFin_VC(iTeFin,nCellFin)
             Ur = StateFin_VC(iUrFin,nCellFin)
          else
             ! Assign weights for linear interpolation between
             ! iCell-1 and iCell
             Weight1 = (rFin_C(iCell) - r) &
                  /    (rFin_C(iCell) - rFin_C(iCell-1))
             Weight2 = 1.0 - Weight1

             Rho  = ( Weight1*StateFin_VC(iRhoFin, iCell-1) &
                  +   Weight2*StateFin_VC(iRhoFin, iCell) )
             Temperature = ( Weight1*StateFin_VC(iTeFin, iCell-1) &
                  +          Weight2*StateFin_VC(iTeFin, iCell) )
             Ur   = ( Weight1*StateFin_VC(iUrFin, iCell-1) &
                  +   Weight2*StateFin_VC(iUrFin, iCell) )
          end if

          U_D(1) = Ur*x/r
          U_D(2) = Ur*y/r

          select case(NameVar)
          case('rho0')
             PlotVar_G(i,j,-1:nK+2) = Rho
          case('t0','temp0')
             PlotVar_G(i,j,-1:nK+2) = Temperature
          case('ux0')
             PlotVar_G(i,j,-1:nK+2) = U_D(1)
          case('uy0')
             PlotVar_G(i,j,-1:nK+2) = U_D(2)
          case default
             IsFound = .false.
          end select
       end do; end do

    case('lowrie')
       select case(NameVar)
       case('ti')
          PlotVar_G = State_VGB(p_,:,:,:,iBlock)/State_VGB(Rho_,:,:,:,iBlock)

       case('te')
          do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  TeOut=TeSi)
             PlotVar_G(i,j,k) = TeSi*Si2No_V(UnitTemperature_)
          end do; end do; end do

       case('rho0','ux0','uy0','ti0','te0')
          ! Calculate sin and cos from the tangent = ShockSlope
          SinSlope = ShockSlope/sqrt(1.0+ShockSlope**2)
          CosSlope =        1.0/sqrt(1.0+ShockSlope**2)

          do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2

             x = x_Blk(i,j,k,iBlock)*CosSlope + y_Blk(i,j,k,iBlock)*SinSlope
             x = x - u0*Time_Simulation*Si2No_V(UnitT_) - x0

             do iCell = 1, nCellRef
                if(rRef_C(iCell) >= x) EXIT
             end do
             if(iCell == 1) call stop_mpi(NameSub // &
                  " Reference solution does not cover the left boundary")

             if(iCell > nCellRef)then
                ! Cell is beyond the last point of reference input:
                ! use last cell
                iCell   = nCellRef
                Weight1 = 0.0
                Weight2 = 1.0
             else
                ! Assign weights for linear interpolation between
                ! iCell-1 and iCell
                Weight1 = (rRef_C(iCell) - x) &
                     /    (rRef_C(iCell) - rRef_C(iCell-1))
                Weight2 = 1.0 - Weight1
             end if

             select case(NameVar)
             case('rho0')
                PlotVar_G(i,j,k) = &
                     ( Weight1*StateRef_VC(iRhoRef, iCell-1) &
                     + Weight2*StateRef_VC(iRhoRef, iCell) )

             case('ux0','uy0')
                Ux = ( Weight1*StateRef_VC(iUxRef, iCell-1) &
                     + Weight2*StateRef_VC(iUxRef, iCell) )

                select case(NameVar)
                case('ux0')
                   PlotVar_G(i,j,k) = (Ux+u0)*CosSlope
                case('uy0')
                   PlotVar_G(i,j,k) = (Ux+u0)*SinSlope
                end select

             case('ti0')
                Ti = ( Weight1*StateRef_VC(iTiRef, iCell-1) &
                     + Weight2*StateRef_VC(iTiRef, iCell) )

                PlotVar_G(i,j,k) = Ti

             case('te0')
                Te = ( Weight1*StateRef_VC(iTeRef, iCell-1) &
                     + Weight2*StateRef_VC(iTeRef, iCell) )

                PlotVar_G(i,j,k) = Te

             end select
          end do; end do; end do

       end select

    case('parcond', 'parcondsemi')
       select case(NameVar)
       case('te0')
          do j=-1,nJ+2; do i=-1,nI+2
             call get_temperature_parcond(i, j, iBlock, Te)
             PlotVar_G(i,j,:) = Te
          end do; end do
       case('te')
          if(UseIdealEos)then
             do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
                if(UseElectronPressure)then
                   PlotVar_G(i,j,k) = State_VGB(Pe_,i,j,k,iBlock) &
                        /State_VGB(Rho_,i,j,k,iBlock)
                else
                   PlotVar_G(i,j,k) = State_VGB(p_,i,j,k,iBlock) &
                        /State_VGB(Rho_,i,j,k,iBlock) &
                        /(1 + ElectronTemperatureRatio)
                end if
             end do; end do; end do
          else
             do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
                call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                     TeOut = PlotVar_G(i,j,k))
                PlotVar_G(i,j,:) = PlotVar_G(i,j,k)*Si2No_V(UnitTemperature_)
             end do; end do; end do
          end if
       case default
          IsFound = .false.
       end select

    end select

  end subroutine user_set_plot_var

  !============================================================================

  subroutine user_material_properties(State_V, i, j, k, iBlock, iDir, &
       EinternalIn, TeIn, NatomicOut, &
       EinternalOut, TeOut, PressureOut, &
       CvOut, GammaOut, HeatCondOut, TeTiRelaxOut, &
       OpacityPlanckOut_W, OpacityRosselandOut_W, &
       PlanckOut_W, CgTeOut_W, CgTgOut_W, TgOut_W)

    ! The State_V vector is in normalized units

    use ModAdvance,    ONLY: nWave, UseElectronPressure
    use ModConst,      ONLY: cBoltzmann
    use ModPhysics,    ONLY: gm1, inv_gm1, No2Si_V, Si2No_V, &
         UnitRho_, UnitP_, UnitEnergyDens_, UnitTemperature_, &
         UnitX_, UnitT_, UnitU_, UnitN_, cRadiationNo, g, Clight
    use ModVarIndexes, ONLY: nVar, Rho_, p_, Pe_, ExtraEint_

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
    real, optional, intent(out) :: TeTiRelaxOut            ! [1/s]
    real, optional, intent(out) :: &
         OpacityPlanckOut_W(nWave)                         ! [1/m]
    real, optional, intent(out) :: &
         OpacityRosselandOut_W(nWave)                      ! [1/m]

    ! Multi-group specific interface. The variables are respectively:
    !  Group Planckian spectral energy density
    !  Derivative of group Planckian by electron temperature
    !  Group specific heat of the radiation
    !  Group radiation temperature
    real, optional, intent(out) :: PlanckOut_W(nWave)      ! [J/m^3]
    real, optional, intent(out) :: CgTeOut_W(nWave)        ! [J/(m^3*K)]
    real, optional, intent(out) :: CgTgOut_W(nWave)        ! [J/(m^3*K)]
    real, optional, intent(out) :: TgOut_W(nWave)          ! [K]

    real :: Rho, Pressure, Te, Ti
    real :: RhoSi, pSi, TeSi
    real :: Cv, HeatCond
    real :: EeSi, Ee, NatomicSi

    character(len=*), parameter :: NameSub = 'user_material_properties'
    !--------------------------------------------------------------------------

    Rho = State_V(Rho_)
    RhoSi = Rho*No2Si_V(Rho_)

    if(present(EinternalIn))then
       if(TypeProblem == 'parcondsemi')then
          Te = (3.5*EinternalIn*Si2No_V(UnitEnergyDens_))**(1.0/3.5)
          Pressure = Rho*Te
          pSi = Pressure*No2Si_V(UnitP_)
       elseif(TypeProblem == 'lowrie')then
          Ee = EinternalIn*Si2No_V(UnitEnergyDens_)
          Te = sqrt(sqrt(Ee/cRadiationNo))
          pSi = (GammaRel - 1)*EinternalIn
       else
          pSi = EinternalIn*gm1
          Pressure = pSi*Si2No_V(UnitP_)
          Te = Pressure/Rho
       end if
       TeSi = Te*No2Si_V(UnitTemperature_)
    elseif(present(TeIn))then
       TeSi = TeIn
       Te = TeSi*Si2No_V(UnitTemperature_)
       if(TypeProblem == 'lowrie')then
          Ee = cRadiationNo*Te**4
          pSi = (GammaRel - 1)*Ee*No2Si_V(UnitP_)
       else
          Pressure = Rho*Te
          pSi = Pressure*No2Si_V(UnitP_)
       end if
    else
       if(TypeProblem == 'lowrie')then
          pSi = State_V(Pe_)*No2Si_V(UnitP_)
          Ee = inv_gm1*State_V(Pe_) + State_V(ExtraEint_)
          Te = sqrt(sqrt(Ee/cRadiationNo))
       else
          if(UseElectronPressure)then
             Pressure = State_V(Pe_)
          else
             Pressure = State_V(p_)
          end if
          pSi = Pressure*No2Si_V(UnitP_)
          Te = Pressure/Rho
       end if
       TeSi = Te*No2Si_V(UnitTemperature_)
    end if

    if(present(EinternalOut))then
       if(TypeProblem=='parcondsemi')then
          EinternalOut = Te**3.5*No2Si_V(UnitEnergyDens_)/3.5
       elseif(TypeProblem == 'lowrie')then
          EinternalOut = cRadiationNo*Te**4*No2Si_V(UnitEnergyDens_)
       else
          EinternalOut = pSi*inv_gm1
       end if
    end if

    if(present(TeOut)) TeOut = TeSi
    if(present(PressureOut)) PressureOut = pSi

    if(present(CvOut))then
       if(TypeProblem == 'parcondsemi')then
          Cv = Te**2.5
       elseif(TypeProblem == 'lowrie')then
          Cv = 4.0*cRadiationNo*Te**3
       else
          Cv = inv_gm1*Rho
       end if
       CvOut = Cv*No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_)
    end if

    if(present(HeatCondOut))then
       select case(TypeProblem)
       case('rmtv')
          HeatCond = Te**6.5/Rho**2
       case('lowrie')
          Ti = State_V(p_)/State_V(Rho_)
          HeatCond = 0.00175*(g*Ti)**3.5/State_V(Rho_)*4.0*cRadiationNo*Te**3
       case default
          HeatCond = HeatConductionCoef
       end select
       HeatCondOut = HeatCond &
            *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_) &
            *No2Si_V(UnitU_)*No2Si_V(UnitX_)
    end if

    if(present(TeTiRelaxOut))then
       if(TypeProblem == 'lowrie')then
          NatomicSi = State_V(Rho_)*No2Si_V(UnitN_)
          Ti = State_V(p_)/State_V(Rho_)
          TeTiRelaxOut = 1.0E6/(0.00175*(g*Ti)**3.5/State_V(Rho_)) &
               *cRadiationNo*(Te + Ti)*(Te**2 + Ti**2) &
               /(cBoltzmann*NatomicSi)*No2Si_V(UnitEnergyDens_) &
               /(No2Si_V(UnitTemperature_)*No2Si_V(UnitT_))
       else
          TeTiRelaxOut = 0.0
       end if
    end if

    if(present(NatomicOut)) NatomicOut = State_V(Rho_)*No2Si_V(UnitN_)
    if(present(OpacityPlanckOut_W)) OpacityPlanckOut_W = 0.0
    if(present(OpacityRosselandOut_W)) OpacityRosselandOut_W = 0.0
    if(present(CgTeOut_W)) CgTeOut_W = 0.0

  end subroutine user_material_properties

  !============================================================================

  real function bessj0(x)

    ! bessj0 from numerical recipes (changed to free format)

    real :: x, z, ax, xx

    real y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6
    data p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4, &
         -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1, &
         .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
    data r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,651619640.7d0, &
         -11214424.18d0,77392.33017d0,-184.9052456d0/, &
         s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0, &
         9494680.718d0,59272.64853d0,267.8532712d0,1.d0/

    !--------------------------------------------------------------------------

    if(abs(x).lt.8.)then
       y=x**2
       bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) &
            /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
    else
       ax=abs(x)
       z=8./ax
       y=z**2
       xx=ax-.785398164
       bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y &
            *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
    endif

  end function bessj0

end module ModUser
