! Wrapper for Global Magnetosphere (GM) component
!==========================================================================
subroutine GM_set_param(CompInfo, TypeAction)

  use CON_comp_info
  use ModProcMH
  use ModIO, ONLY: iUnitOut, StringPrefix, STDOUT_, NamePlotDir
  use ModRestartFile, ONLY: NameRestartInDir, NameRestartOutDir
  use ModMain, ONLY : CodeVersion, NameThisComp, &
       time_accurate, StartTime, iStartTime_I, &
       dt_UpdateB0, DoUpdateB0, UseRotatingBc
  use CON_physics, ONLY: get_time, get_planet
  use ModTimeConvert, ONLY: time_real_to_int

  implicit none

  character (len=*), parameter :: NameSub='GM_set_param'

  ! Arguments
  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do

  logical :: DoTest,DoTestMe
  !-------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest,DoTestMe)

  if(DoTest)write(*,*)NameSub,' called with TypeAction,iProc=',&
       TypeAction,iProc

  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,&
          Use        =.true.,                        &
          NameVersion='BATSRUS (Univ. of Michigan)', &
          Version    =CodeVersion)
  case('MPI')
     call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc,&
          Name=NameThisComp)

     NamePlotDir(1:2)       = NameThisComp
     NameRestartInDir(1:2)  = NameThisComp
     NameRestartOutDir(1:2) = NameThisComp
  case('READ')
     call MH_set_parameters('READ')
  case('CHECK')
     call get_time( &
          DoTimeAccurateOut = time_accurate, &
          tStartOut         = StartTime)
     call get_planet( &
          DtUpdateB0Out  = dt_updateB0,   &
          DoUpdateB0Out  = DoUpdateB0,    &
          UseRotationOut = UseRotatingBc)
     call time_real_to_int(StartTime,iStartTime_I)

     call MH_set_parameters('CHECK')
  case('STDOUT')
     iUnitOut=STDOUT_
     if(iProc==0)then
        StringPrefix = NameThisComp//':'
     else
        write(StringPrefix,'(a,i4.4,a)')NameThisComp,iProc,':'
     end if
  case('FILEOUT')
     call get(CompInfo,iUnitOut=iUnitOut)
     StringPrefix=''
  case('GRID')
     call GM_set_grid
  case default
     call CON_stop(NameSub//' SWMF_ERROR: invalid TypeAction='//TypeAction)
  end select

end subroutine GM_set_param
!======================================================================
!BOP
!ROUTINE: GM_set_grid - intialize, set and broadcast adaptive block grid
!INTERFACE:
subroutine GM_set_grid
  !USES:
  use MH_domain_decomposition
  use CON_coupler
  use CON_test_global_message_pass
  use ModMain,ONLY:TypeCoordSystem
  use ModPhysics,ONLY:No2Si_V, UnitX_
  use CON_comp_param,ONLY:GM_

  implicit none
  !REVISION HISTORY:
  !23Aug03 I.Sokolov <igorsok@umich.edu> - initial prototype/prolog/code
  !03Sep03 G.Toth    <gtoth@umich.edu> - removed test_message_pass
  !                                      call synchronize_refinement directly
  !EOP
  logical ::DoTest,DoTestMe
  DoTest=.false.;DoTestMe=.false.
  if(done_dd_init(GM_))return
  call init_decomposition(GM_,GM_,3,.true.)
  call set_coord_system(GM_,TypeCoordSystem,No2Si_V(UnitX_))
  
  if(is_proc(GM_))then
     call init_decomposition(&
          MH_DomainDecomposition,GM_,3,.true.)
     call MH_get_root_decomposition(MH_DomainDecomposition)
     call MH_update_local_decomposition(MH_DomainDecomposition)
     MH_DomainDecomposition%IsLocal=.true.
  end if
  call CON_set_do_test('test_grids',DoTest,DoTestMe)


  if(is_proc0(GM_))call MH_get_root_decomposition(GM_)

  call bcast_decomposition(GM_)

  call synchronize_refinement(GM_,MH_domaindecomposition)

  if(DoTest) call test_global_message_pass(GM_)
end subroutine GM_set_grid
!===================================================================!
!BOP
!ROUTINE: GM_synchronize_refinement - synchronize global grid for GM_
!INTERFACE:
subroutine GM_synchronize_refinement(iProc0,iCommUnion)

  !USES:
  use ModProcMH
  use MH_domain_decomposition
  use CON_comp_param,ONLY:GM_
  !INPUT ARGUMENTS:
  integer,intent(in) :: iProc0,iCommUnion
  !REVISION HISTORY:
  !23AUG03  I.Sokolov <igorsok@umich.edu> - initial prototype/code/prolog
  !03SEP03  G.Toth    <gtoth@umich.edu> - arguments are not optional now
  !EOP

  if(is_proc(GM_)) &
       call MH_update_local_decomposition(MH_DomainDecomposition)

  call synchronize_refinement(&
       GM_,MH_domaindecomposition,iProc0,iCommUnion)

end subroutine GM_synchronize_refinement
!==============================================================================
subroutine GM_print_variables(NameSource)

  use ModMain, ONLY: NameThisComp
  use ModNumConst
  use ModImPressure                                  !^CFG IF RCM
  use ModIonoPotential,ONLY:IonoPotential_II
  use ModIeGrid,       ONLY:nThetaIono, nPhiIono, ThetaIono_I, PhiIono_I
  use ModIoUnit, ONLY: UNITTMP_
  implicit none
  character(len=*), parameter :: NameSub='GM_print_variables'

  character(len=*),intent(in) :: NameSource
  integer            :: nFile=0
  character(len=100) :: NameFile
  character(len=100) :: NameVar
  integer            :: i,j
  !--------------------------------------------------------------------------

  if(NameThisComp/='GM') RETURN

  select case(NameSource)
  case('IM')                        !^CFG IF RCM
     NameVar='j i lon lat density pressure' !^CFG IF RCM
  case('IE','IE_swmf')
     NameVar='i j theta phi pot'
  case default
     write(*,*)NameSub,': incorrect NameSource=',NameSource
     RETURN
  end select

  nFile=nFile+1
  write(NameFile,'(a,i1,a)')'GM_from_'//NameSource//'_',nFile,'.dat'
  open(UNITTMP_,file=NameFile)
  write(UNITTMP_,'(a)')trim(NameVar)

  select case(NameSource)
  case('IM')                             !^CFG IF RCM BEGIN
     do i=1,iSize
        do j=1,jSize
           write(UNITTMP_,'(2i4,4G14.6)')j,i,RCM_lon(j),RCM_lat(i), &
                RCM_dens(i,j),RCM_p(i,j)
        end do
     end do                              !^CFG END RCM
  case('IE')
     do j=1,nPhiIono
        do i=1,nThetaIono
           write(UNITTMP_,'(2i4,3G14.6)')i,j,&
                ThetaIono_I(i),PhiIono_I(j),IonoPotential_II(i,j)
        end do
     end do
  end select
  close(UNITTMP_)

end subroutine GM_print_variables

!=============================================================================

subroutine GM_init_session(iSession, TimeSimulation)

  use ModProcMH,   ONLY: iProc
  use ModMain,     ONLY: Time_Simulation, UseIe, UsePw, TypeBC_I, west_, body1_
  use ModMain,     ONLY: UseIM                            !^CFG IF RCM
  use CON_physics, ONLY: get_time
  use CON_coupler, ONLY: Couple_CC, IE_, IM_, GM_, IH_, PW_
  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='GM_init_session'

  logical :: IsUninitialized = .true.
  logical :: DoTest, DoTestMe
  !----------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest, DoTestMe)

  UseIm = Couple_CC(IM_,GM_) % DoThis !^CFG IF RCM
  UsePw = Couple_CC(PW_,GM_) % DoThis
  UseIe = Couple_CC(IE_,GM_) % DoThis

  ! Check if the boundary condition is properly set
  if(UsePw) TypeBC_I(body1_) = 'polarwind'

  if(Couple_CC(IH_,GM_) % DoThis .neqv. (TypeBc_I(west_)=='coupled'))then
     if(Couple_CC(IH_,GM_) % DoThis) then
        TypeBc_I(west_)='coupled'
     else
        if(iProc==0)write(*,*)NameSub,' WARNING: IH and GM are not coupled,',&
             ' changing west boundary type from "coupled" to "vary"'
        TypeBc_I(west_)='vary'
     end if
  end if

  if(IsUninitialized)then

     call get_time(tSimulationOut=Time_Simulation)
     call BATS_setup
     IsUninitialized = .false.
  end if
  call BATS_init_session

  if(DoTest)write(*,*)NameSub,' finished for session ',iSession

end subroutine GM_init_session

!==============================================================================

subroutine GM_finalize(TimeSimulation)

  use ModMain, ONLY: UseIe, time_loop
  use ModFieldAlignedCurrent, ONLY: &     !^CFG IF IONOSPHERE
       clean_mod_field_aligned_current    !^CFG IF IONOSPHERE
  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='GM_finalize'

  integer :: iError
  !----------------------------------------------------------------------------
  ! We are not advancing in time any longer
  time_loop = .false.

  call BATS_save_files('FINAL')

  if (UseIe) call clean_mod_field_aligned_current  !^CFG IF IONOSPHERE

  call error_report('PRINT',0.,iError,.true.)

end subroutine GM_finalize

!==============================================================================

subroutine GM_save_restart(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='GM_save_restart'

  call BATS_save_files('RESTART')

end subroutine GM_save_restart

!=============================================================================

subroutine GM_run(TimeSimulation,TimeSimulationLimit)

  use ModProcMH, ONLY: iProc
  use ModMain,   ONLY: Time_Simulation

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='GM_run'

  logical :: DoTest, DoTestMe
  !----------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest,DoTestMe)

  if(DoTest)write(*,*)NameSub,' called with tSim, tSimLimit, iProc=',&
       TimeSimulation, TimeSimulationLimit, iProc

  if(abs(Time_Simulation-TimeSimulation)>0.0001) then
     write(*,*)NameSub,' GM time=',Time_Simulation,' SWMF time=',TimeSimulation
     call CON_stop(NameSub//' SWMF_ERROR: GM and SWMF simulation times differ')
  end if

  call BATS_advance(TimeSimulationLimit)

  ! Return time after the time step
  TimeSimulation = Time_Simulation

end subroutine GM_run

