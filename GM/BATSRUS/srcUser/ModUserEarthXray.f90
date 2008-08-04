! Note that for some reason that I don't understand yet, this user file 
! will only work with the normalization type set to "solarwind".  You 
! should set this in the PARAM.in file.



!========================================================================
Module ModUser
  use ModNumConst, ONLY: cHalf,cTwo,cThree,&
       cFour,cE1,cHundred,cHundredth,cZero,&
       cOne,cTiny
  use ModSize,     ONLY: nI,nJ,nK,gcn,nBLK,MaxBlock
  use ModVarIndexes, ONLY: nVar
  use ModMain, ONLY: iTest, jTest, kTest, BlkTest, ProcTest, n_step
  use ModProcMH, ONLY: iProc
  use ModUserEmpty,               &
       IMPLEMENTED1 => user_init_session,               &
       IMPLEMENTED2 => user_set_ics,                    &
       IMPLEMENTED3 => user_set_outerbcs,               &
       IMPLEMENTED4 => user_set_plot_var,               &
       IMPLEMENTED5 => user_update_states

  include 'user_module.h' !list of public methods

  !summed MHD quantities
  integer, parameter, public :: MaxSumMhdVar=nVar+2 !(8+2+2)
  integer, parameter, public :: Umag_=nVar+1
  integer, parameter, public :: Tsw_ =nVar+2 
  real :: tSumStart, tSumEnd
  real :: StateSum_VC(MaxSumMhdVar,nI, nJ, nK)
  logical :: IsRestartSum(MaxBlock)
 
  real, parameter :: VersionUserModule = 1.1
  character (len=*), parameter :: &
       NameUserModule = 'Earth Mag X-ray (EarthXray), Hansen, Jan, 2008'

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


!    !\
!    ! Initiallize the arrays that contain the time average of the State 
!    ! variables with the initial conditions of the cells
!    !/
!    tSumStart = Time_Simulation
!    StateSum_VC = 0.0
!    if(iBlock /= iBlockLast)then
!       iBlockLast = iBlock
!       call put_block_data(iBlock, MaxSumMhdVar, nI, nJ, nK, StateSum_VC)
!    end if
!
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

    if(TypeBc=='uservary'.and.time_accurate)then
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
    integer :: iBlockLast = -1

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
    if(iBlock /= iBlockLast)then
       iBlockLast = iBlock
       if(use_block_data(iBlock))then
          call get_block_data(iBlock, MaxSumMhdVar, nI, nJ, nK, &
               StateSum_VC, DoNotAdvance=.true.)

          !reset the block data for summing over the next period if it was just printed
          if (IsRestartSum(iBlock)) then
             IsRestartSum(iBlock) = .false.
             StateSum_VC = 0.0
             tSumStart = Time_Simulation
             tSumEnd = tSumStart
          else
             StateSum_VC(rho_,:,:,:) = StateSum_VC(rho_,:,:,:) + &
                  State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock)*dt
             StateSum_VC(rhosw_,:,:,:) = StateSum_VC(rhosw_,:,:,:) + &
                  State_VGB(rhosw_,1:nI,1:nJ,1:nK,iBlock)*dt
             StateSum_VC(rhoion_,:,:,:) = StateSum_VC(rhoion_,:,:,:) + &
                  State_VGB(rhoion_,1:nI,1:nJ,1:nK,iBlock)*dt
             StateSum_VC(Bx_,:,:,:) = StateSum_VC(Bx_,:,:,:) + &
                  State_VGB(Bx_,1:nI,1:nJ,1:nK,iBlock)*dt
             StateSum_VC(By_,:,:,:) = StateSum_VC(By_,:,:,:) + &
                  State_VGB(By_,1:nI,1:nJ,1:nK,iBlock)*dt
             StateSum_VC(Bz_,:,:,:) = StateSum_VC(Bz_,:,:,:) + &
                  State_VGB(Bz_,1:nI,1:nJ,1:nK,iBlock)*dt
             StateSum_VC(p_,:,:,:) = StateSum_VC(p_,:,:,:) + &
                  State_VGB(p_,1:nI,1:nJ,1:nK,iBlock)*dt
             StateSum_VC(Ux_,:,:,:) = StateSum_VC(Ux_,:,:,:) + &
                  State_VGB(rhoux_,1:nI,1:nJ,1:nK,iBlock) /   &
                  State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock)*dt
             StateSum_VC(Uy_,:,:,:) = StateSum_VC(Uy_,:,:,:) + &
                  State_VGB(rhouy_,1:nI,1:nJ,1:nK,iBlock) /   &
                  State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock)*dt
             StateSum_VC(Uz_,:,:,:) = StateSum_VC(Uz_,:,:,:) + &
                  State_VGB(rhouz_,1:nI,1:nJ,1:nK,iBlock) /   &
                  State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock)*dt
             StateSum_VC(Umag_,:,:,:) = StateSum_VC(Umag_,:,:,:) + &
                  sqrt(State_VGB(RhoUx_,1:nI,1:nJ,1:nK,iBlock)**2 + &
                       State_VGB(RhoUy_,1:nI,1:nJ,1:nK,iBlock)**2 + &
                       State_VGB(RhoUz_,1:nI,1:nJ,1:nK,iBlock)**2 )/ &
                  State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock)*dt
             ! Note that the temperature of the solar wind plasma is only
             ! going to work here for cells were the ionospheric plasma is 
             ! small since we are using the full pressure and the partial 
             ! density.
             StateSum_VC(Tsw_,:,:,:) = StateSum_VC(Tsw_,:,:,:) + &
                  State_VGB(p_,1:nI,1:nJ,1:nK,iBlock)/           &
                  State_VGB(rhosw_,1:nI,1:nJ,1:nK,iBlock)*dt

          end if

          call put_block_data(iBlock, MaxSumMhdVar, nI, nJ, nK, &
               StateSum_VC, DoAllowReplace=.true.)

       else

          IsRestartSum(iBlock) = .false.
          tSumStart = Time_Simulation
          StateSum_VC = 0.0
          call put_block_data(iBlock, MaxSumMhdVar, nI, nJ, nK, StateSum_VC)

       end if
    end if

  end subroutine user_update_states

  !====================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModPhysics
    use ModMain, ONLY: Body1_, time_simulation,x_,y_,z_
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
    integer :: iBlockLast = -1

    character (len=*), parameter :: Name='user_set_plot_var'
    logical :: oktest,oktest_me

    !------------------------------------------------------------------------  
    if(iProc==PROCtest .and. iBlock==BLKtest)then
       call set_oktest('user_set_plot_var',oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    end if

    !If we are in this routine then reset the logical that tell the sum to start over

    IsRestartSum = .true.

    ! Get the averaged data from the Block storage arrays if this block has not been read yet
    ! if it has just use the already saved StateSum_VC which should have the right info in it
    if(iBlock /= iBlockLast) then
       iBlockLast = iBlock
       if(use_block_data(iBlock)) &
          call get_block_data(iBlock, MaxSumMhdVar, nI, nJ, nK, StateSum_VC)
    end if

    ! Now load the plot arrays
    if(use_block_data(iBlock))then

       isFound = .true.
    
       select case(NameVar)
       case('rhoave')
          NameTecVar = '<`r>'
          iUnitVar   = UnitRho_
          PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(rho_,:,:,:)

       case('rhoswave')
          NameTecVar = '<`rsw>'
          iUnitVar   = UnitRho_
          PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(rhosw_,:,:,:)
       case('rhoionave')
          NameTecVar = '<`rion>'
          iUnitVar   = UnitRho_
          PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(rhoion_,:,:,:)
       case('bxave') 
          NameTecVar = '<B_x>'
          iUnitVar   = UnitB_
          ! Note: here we add B0x to the summed B1.  We have to multiply by the time because the
          ! summed variables have been multiplied by dt in the summation.
          PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(Bx_,:,:,:)+ &
               B0_DGB(x_,1:nI,1:nJ,1:nK,iBlock)*(tSumEnd-tSumStart)
       case('byave') 
          NameTecVar = '<B_y>'
          iUnitVar   = UnitB_
          ! Note: here we add B0y to the summed B1.  We have to multiply by the time because the
          ! summed variables have been multiplied by dt in the summation.
          PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(By_,:,:,:)+ &
               B0_DGB(y_,1:nI,1:nJ,1:nK,iBlock)*(tSumEnd-tSumStart)
       case('bzave') 
          NameTecVar = '<B_z>'
          iUnitVar   = UnitB_
          ! Note: here we add B0z to the summed B1.  We have to multiply by the time because the
          ! summed variables have been multiplied by dt in the summation.
          PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(Bz_,:,:,:)+ &
               B0_DGB(z_,1:nI,1:nJ,1:nK,iBlock)*(tSumEnd-tSumStart)
       case('uxave') 
          NameTecVar = '<U_x>'
          iUnitVar   = UnitU_
          PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(RhoUx_,:,:,:)
       case('uyave') 
          NameTecVar = '<U_y>'
          iUnitVar   = UnitU_
          PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(RhoUy_,:,:,:)
       case('uzave') 
          NameTecVar = '<U_z>'
          iUnitVar   = UnitU_
          PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(RhoUz_,:,:,:)
       case('uave') 
          NameTecVar = '<|U|>'
          iUnitVar   = UnitU_
          PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(Umag_,:,:,:)
       case('pave','pthave')
          NameTecVar = '<p>'
          iUnitVar   = UnitP_
          PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(p_,:,:,:)
       case('tswave')
          NameTecVar = '<Tsw>'
          iUnitVar   = UnitTemperature_
          PlotVar_G(1:nI,1:nJ,1:nK)  = StateSum_VC(Tsw_,:,:,:)
       case default
          IsFound = .false.
          call stop_mpi(Name//': unimplemented variable='//NameVar)
       end select
    
       UsePlotVarBody = .true.
       PlotVarBody    = 0.0
    
       ! The Sum variables store the summation over time.  We want averages so divide by
       ! the elapsed time.
       PlotVar_G = PlotVar_G/((tSumEnd-tSumStart)*Io2No_V(UnitT_))
    
       if(IsDimensional) PlotVar_G = PlotVar_G*No2Io_V(iUnitVar)
    
    end if
          
  end subroutine user_set_plot_var


end module ModUser

