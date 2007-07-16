! You must set UseUserInitSession, UseUserSetIcs and UseUserOuterBcs 
! to .true. for this user module to be effective.

!========================================================================
Module ModUser
  use ModNumConst, ONLY: cHalf,cTwo,cThree,&
       cFour,cE1,cHundred,cHundredth,cZero,&
       cOne,cTiny
  use ModSize,     ONLY: nI,nJ,nK,gcn,nBLK
  use ModVarIndexes, ONLY: nVar
  use ModUserEmpty,               &
       IMPLEMENTED1 => user_init_session,               &
       IMPLEMENTED2 => user_set_ics,                    &
       IMPLEMENTED3 => user_set_outerbcs,               &
       IMPLEMENTED4 => user_set_plot_var,               &
       IMPLEMENTED5 => user_update_states

  include 'user_module.h' !list of public methods

  !summed MHD quantities
  integer, public :: MaxSumMhdVar=nVar !(8+2)
  real :: tSumStart, tSumEnd
  real :: StateSum_VC(nVar,nI, nJ, nK)
 
  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'Earth Mag X-ray (EarthXray), Hansen, July, 2007'

contains

  !=====================================================================
  subroutine user_init_session

  use ModVarIndexes
  use ModPhysics, ONLY: FaceState_VI, CellState_VI, SW_rho, BodyRho_I
  use ModNumConst, ONLY: cTiny
  use ModSize, ONLY: east_, west_, south_, north_, bot_, top_
  use ModMain, ONLY: body1_
  integer :: iBoundary

    !-------------------------------------------------------------------

    !\
    ! We are using this routine to initialize the arrays that control the
    ! default value for the inner body and hence the boundary condition.
    ! Note these values are typically set in set_physics and they are set to
    ! the BodyRho_I value read from the PARAM.in file.  We want to use 
    ! this same strategy for multi-species but have to do it here to avoid 
    ! modifying the core source code.
    !/
  
    ! FaceState_VI is used to set the inner boundary condition.  Setting
    ! the correct values here for the extra species will assure that 
    ! the inner boundary is done correctly.
    FaceState_VI(rhosw_,body1_) =cTiny*BodyRho_I(1)
    FaceState_VI(rhoion_,body1_)=BodyRho_I(1)

    ! We set the following array for the outer boundaries.  Although
    ! only CellState_VI is used we set both.  Not that these are 
    ! used for only some outerboundary cases (fixed, for example) and
    ! are ignored for vary and other types.  We code them as in set_physics
    ! just to be safe.
    do iBoundary=east_,top_
       FaceState_VI(rhosw_, iBoundary)  = SW_rho
       FaceState_VI(rhoion_, iBoundary) = cTiny*sw_rho
    end do
    CellState_VI=FaceState_VI

  end subroutine user_init_session

  !=====================================================================
  subroutine user_set_ics
    use ModMain,     ONLY: globalBLK, time_simulation, dt
    use ModGeometry, ONLY: r_BLK
    use ModAdvance,  ONLY: State_VGB, rhoion_, rhosw_
    use ModPhysics,  ONLY: BodyRho_I, sw_rho, rBody
    use ModNumConst, ONLY: cTiny
    use ModBlockData,ONLY: put_block_data
    implicit none

    integer :: iBlock,iBlockLast = -1


    !--------------------------------------------------------------------------
    iBlock = globalBLK

    where(r_BLK(:,:,:,iBlock)<2.0*Rbody)
       State_VGB(rhoion_,:,:,:,iBlock) = BodyRho_I(1)
       State_VGB(rhosw_,:,:,:,iBlock)  = cTiny*sw_rho
    elsewhere
       State_VGB(rhoion_,:,:,:,iBlock) = cTiny*BodyRho_I(1)
       State_VGB(rhosw_,:,:,:,iBlock)  = sw_rho
    end where


    !\
    ! Initiallize the arrays that contain the time average of the State 
    ! variables with the initial conditions of the cells
    !/
    tSumStart = Time_Simulation
    StateSum_VC = 0.0
    if(iBlock /= iBlockLast)then
       iBlockLast = iBlock
       call put_block_data(iBlock, MaxSumMhdVar, nI, nJ, nK, StateSum_VC)
    end if


  end subroutine user_set_ics


  !=====================================================================
  subroutine user_set_outerbcs(iBlock,iSide, TypeBc,found)

  use ModMain,      ONLY : time_simulation
  use ModMain,      ONLY : time_accurate
  use ModVarIndexes
  use ModAdvance,   ONLY : State_VGB
  use ModPhysics,   ONLY : CellState_VI

    integer,intent(in)::iBlock, iSide
    logical,intent(out) :: found
    character (len=20),intent(in) :: TypeBc
    real :: time_now

    character (len=*), parameter :: Name='user_set_outerbcs'

    time_now = time_simulation

    !-------------------------------------------------------------------

    if(TypeBc=='vary'.and.time_accurate)then
       call BC_solar_wind(time_now)
    else
       call BC_fixed(1,nVar,CellState_VI(:,iSide))
       call BC_fixed_B
    end if

    ! Note that the above code does not set the extra density species (vary)
    ! or sets them to undefined values (fixed). 
    !
    ! The solar wind species is the only one at the upstream boundary.
    ! The ionosphere species is zero.
    State_VGB(rhosw_,:,:,:,iBlock)   = State_VGB(rho_,:,:,:,iBlock)
    State_VGB(rhoion_,:,:,:,iBlock)  = cTiny*State_VGB(rho_,:,:,:,iBlock)

    found = .true.

  end subroutine user_set_outerbcs

  !===========================================================================
  subroutine user_update_states(iStage,iBlock)
    use ModVarIndexes
    use ModSize
    use ModAdvance, ONLY: State_VGB
    use ModMain,    ONLY: nStage, time_simulation, dt
    use ModBlockData,ONLY: get_block_data, put_block_data, use_block_data

    implicit none
    integer,intent(in):: iStage,iBlock

    !--------------------------------------------------------------------------

    !\
    ! do the normal update states
    !/
    call update_states_MHD(iStage,iBlock)

    !\
    ! Now compute the sum of the state variables (we will divide by the number of
    ! time steps later to get the average state value
    !/
    tSumEnd = Time_Simulation
    if(use_block_data(iBlock))then
       call get_block_data(iBlock, MaxSumMhdVar, nI, nJ, nK, &
            StateSum_VC, DoNotAdvance=.true.)
       StateSum_VC = StateSum_VC + State_VGB(:,1:nI,1:nJ,1:nK,iBlock)*dt
!       StateSum_VC = StateSum_VC
       call put_block_data(iBlock, MaxSumMhdVar, nI, nJ, nK, &
            StateSum_VC, DoAllowReplace=.true.)
    end if

  end subroutine user_update_states

  !====================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModPhysics
    use ModMain, ONLY: Body1_, time_simulation
    use ModAdvance
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, r_BLK, IsBoundaryBlock_IB
    use ModMain, ONLY: iTest, jTest, kTest, ProcTest, BlkTest, GLOBALBLK
    use ModProcMH,   ONLY: iProc
    use ModBlockData,ONLY: get_block_data, put_block_data, use_block_data

    implicit none

    integer,          intent(in) :: iBlock
    character(len=*), intent(in) :: NameVar
    logical,          intent(in) :: IsDimensional
    real,             intent(out):: PlotVar_G(-1:nI+2, -1:nJ+2, -1:nK+2)
    real,             intent(out):: PlotVarBody
    logical,          intent(out):: UsePlotVarBody
    character(len=*), intent(out):: NameTecVar
    character(len=*), intent(out):: NameTecUnit
    character(len=*), intent(out):: NameIdlUnit
    logical,          intent(out):: IsFound

    integer :: iUnitVar

    character (len=*), parameter :: Name='user_set_plot_var'

    logical :: oktest,oktest_me

    !------------------------------------------------------------------------  
    if(iProc==PROCtest .and. iBlock==BLKtest)then
       call set_oktest('user_set_plot_var',oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    end if

    ! Get the averaged data from the Block storage arrays
    call get_block_data(iBlock, MaxSumMhdVar, nI, nJ, nK, &
         StateSum_VC, DoNotAdvance=.true.)

    select case(NameVar)
    case('rhoave')
       NameTecVar = '<`r>'
       iUnitVar   = UnitRho_
       PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(rho_,:,:,:)
    case('rhoswave')
       NameTecVar = '<`r>'
       iUnitVar   = UnitRho_
       PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(rhosw_,:,:,:)
    case('rhoionave')
       NameTecVar = '<`r>'
       iUnitVar   = UnitRho_
       PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(rhoion_,:,:,:)
    case('bx') 
       NameTecVar = '<B_x>'
       iUnitVar   = UnitB_
       PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(Bx_,:,:,:)
    case('by') 
       NameTecVar = '<B_y>'
       iUnitVar   = UnitB_
       PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(By_,:,:,:)
    case('bz') 
       NameTecVar = '<B_z>'
       iUnitVar   = UnitB_
       PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(Bz_,:,:,:)
    case('ux') 
       NameTecVar = '<U_x>'
       iUnitVar   = UnitU_
       PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(RhoUx_,:,:,:)/ &
            StateSum_VC(Rho_,:,:,:)
    case('uy') 
       NameTecVar = '<U_y>'
       iUnitVar   = UnitU_
       PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(RhoUy_,:,:,:)/ &
            StateSum_VC(Rho_,:,:,:)
    case('uz') 
       NameTecVar = '<U_z>'
       iUnitVar   = UnitU_
       PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(RhoUz_,:,:,:)/ &
            StateSum_VC(Rho_,:,:,:)
    case('p','pth')
       NameTecVar = '<p>'
       iUnitVar   = UnitP_
       PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(p_,:,:,:)
    case default
       call stop_mpi(Name//': unimplemented variable='//NameVar)
    end select

    UsePlotVarBody = .true.
    PlotVarBody    = 0.0

    ! The Sum variables store the summation over time.  We want averages so divide by
    ! the elapsed time.
    PlotVar_G = PlotVar_G/(tSumEnd-tSumStart)

    if(IsDimensional) PlotVar_G = PlotVar_G*No2Io_V(iUnitVar)

    !reset the block data for summing over the next period
    StateSum_VC = 0.0
    tSumStart = Time_Simulation
    tSumEnd = tSumStart
    call put_block_data(iBlock, MaxSumMhdVar, nI, nJ, nK, &
         StateSum_VC, DoAllowReplace=.true.)
    
  end subroutine user_set_plot_var


end module ModUser

