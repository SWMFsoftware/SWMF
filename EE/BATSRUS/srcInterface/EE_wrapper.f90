! Wrapper for Eruptive Event (EE) component
!==============================================================================
subroutine EE_set_param(CompInfo, TypeAction)

  use CON_comp_info
  use CON_physics, ONLY: get_time
  use EE_ModProcMH
  use EE_ModIO, ONLY: iUnitOut, StringPrefix, STDOUT_, NamePlotDir
  use EE_ModRestartFile, ONLY: NameRestartInDir, NameRestartOutDir
  use EE_ModMain, ONLY : CodeVersion, NameThisComp, &
       time_accurate, StartTime, iStartTime_I
  use ModTimeConvert, ONLY: time_real_to_int
  implicit none

  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do

  logical :: DoTest, DoTestMe

  character (len=*), parameter :: NameSub='EE_set_param'
  !----------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  if(DoTest)write(*,*)NameSub,' called with TypeAction, iProc=', &
       TypeAction,iProc

  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,&
          Use        =.true.,                        &
          NameVersion='EE_BATSRUS (Univ. of Michigan)', &
          Version    =CodeVersion)
  case('MPI')
     call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc,&
          Name=NameThisComp)

     NamePlotDir(1:2)       = NameThisComp
     NameRestartInDir(1:2)  = NameThisComp
     NameRestartOutDir(1:2) = NameThisComp
  case('READ')
     call EE_set_parameters('READ')
  case('CHECK')
     call get_time( &
          DoTimeAccurateOut = time_accurate, &
          tStartOut         = StartTime)
     call time_real_to_int(StartTime, iStartTime_I)

     call EE_set_parameters('CHECK')
  case('STDOUT')
     iUnitOut=STDOUT_
     if(iProc==0)then
        StringPrefix = NameThisComp//':'
     else
        write(StringPrefix,'(a,i4.4,a)')NameThisComp,iProc,':'
     end if
  case('FILEOUT')
     call get(CompInfo, iUnitOut=iUnitOut)
     StringPrefix=''
  case('GRID')
     call EE_set_grid
  case default
     call CON_stop(NameSub//' SWMF_ERROR: invalid TypeAction='//TypeAction)
  end select

end subroutine EE_set_param

!==============================================================================
subroutine EE_set_grid

  use CON_coupler
  use CON_comp_param, ONLY: EE_
  use EE_domain_decomposition
  use EE_ModGeometry, ONLY: TypeGeometry
  use EE_ModMain, ONLY: TypeCoordSystem, NameVarCouple
  use EE_ModPhysics, ONLY: No2Si_V, UnitX_
  implicit none

  character(len=*), parameter :: NameSub='EE_set_grid'
  !----------------------------------------------------------------------------

  if(done_dd_init(EE_))RETURN

  call init_decomposition( &
       GridID_ = EE_, &
       CompID_ = EE_, &
       nDim = 3,      &
       IsTreeDecomposition = .true.)

  call set_coord_system( &
       GridID_ = EE_, &
       TypeCoord = TypeCoordSystem, &
       UnitX = No2Si_V(UnitX_), &
       NameVar = NameVarCouple)

  if(is_proc(EE_))then
     ! Initialize the local grid
     call init_decomposition(&
          DomainDecomposition = MH_DomainDecomposition, &
          CompID_ = EE_, &
          nDim = 3, &
          IsTreeDecomposition = .true.)

     ! Get the octree root array
     call MH_get_root_decomposition(MH_DomainDecomposition)

     ! Get the whole octree after the initial refinement
     call MH_update_local_decomposition(MH_DomainDecomposition)
     MH_DomainDecomposition%IsLocal=.true.
  end if

  ! Repeat the initialization at the global grid level:
  ! Octree root array:
  if(is_proc0(EE_))call MH_get_root_decomposition(EE_)

  ! Broadcast root array:
  call bcast_decomposition(EE_)

  ! Synchronize global and local grids:
  call synchronize_refinement( &
       GridID_ = EE_, &
       localDD = MH_domaindecomposition)

end subroutine EE_set_grid

!==============================================================================
subroutine EE_init_session(iSession, TimeSimulation)

  use EE_ModMain,  ONLY: Time_Simulation
  use CON_physics, ONLY: get_time

  implicit none

  integer, intent(in) :: iSession         ! session number (starting from 1)
  real,    intent(in) :: TimeSimulation   ! seconds from start time

  logical :: IsUninitialized = .true.
  logical :: DoTest, DoTestMe

  character(len=*), parameter :: NameSub='EE_init_session'
  !----------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  if(IsUninitialized)then

     call get_time(tSimulationOut=Time_Simulation)

     call EE_BATS_setup
     IsUninitialized = .false.
  end if
  call EE_BATS_init_session

  if(DoTest)write(*,*)NameSub,' finished for session ',iSession

end subroutine EE_init_session

!==============================================================================
subroutine EE_finalize(TimeSimulation)

  use EE_ModMain, ONLY: time_loop
  implicit none

  real, intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='EE_finalize'
  !----------------------------------------------------------------------------
  ! We are not advancing in time any longer
  time_loop = .false.

  call EE_BATS_save_files('FINAL')

  call EE_BATSRUS_finalize

end subroutine EE_finalize

!==============================================================================
subroutine EE_save_restart(TimeSimulation)

  implicit none

  real, intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='EE_save_restart'
  !----------------------------------------------------------------------------
  call EE_BATS_save_files('RESTART')

end subroutine EE_save_restart

!==============================================================================
subroutine EE_run(TimeSimulation, TimeSimulationLimit)

  use EE_ModProcMH, ONLY: iProc
  use EE_ModMain,   ONLY: Time_Simulation
  implicit none

  real, intent(inout) :: TimeSimulation   ! current time of component
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  logical :: DoTest, DoTestMe

  character(len=*), parameter :: NameSub='EE_run'
  !----------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  if(DoTest)write(*,*)NameSub,' called with tSim, tSimLimit, iProc=',&
       TimeSimulation, TimeSimulationLimit, iProc

  if(abs(Time_Simulation-TimeSimulation)>0.0001) then
     write(*,*)NameSub,' EE time=',Time_Simulation,' SWMF time=',TimeSimulation
     call CON_stop(NameSub//' SWMF_ERROR: EE and SWMF simulation times differ')
  end if

  call EE_BATS_advance(TimeSimulationLimit)

  ! Return time after the time step
  TimeSimulation = Time_Simulation

end subroutine EE_run
