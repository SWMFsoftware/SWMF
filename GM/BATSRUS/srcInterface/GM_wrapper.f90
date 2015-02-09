!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module GM_wrapper

  ! Wrapper for BATSRUS Global Magnetosphere (GM) component

  use GM_couple_pt
  use GM_couple_ie
  use GM_couple_im
  use GM_couple_ih
  use GM_couple_rb
  use GM_couple_pw
  use GM_couple_pc

  implicit none

  private ! except

  ! CON wrapper
  public:: GM_set_param
  public:: GM_init_session
  public:: GM_run
  public:: GM_save_restart
  public:: GM_finalize
  public:: GM_print_variables

  ! Coupling toolkit
  public:: GM_synchronize_refinement

  ! Point coupling
  public:: GM_get_grid_info
  public:: GM_find_points

  ! IE coupling
  public:: GM_get_for_ie
  public:: GM_put_from_ie
  public:: GM_put_mag_from_ie

  ! IH coupling
  public:: GM_put_from_ih           ! coupling toolkit based coupler
  public:: GM_put_from_ih_buffer    ! buffer grid based coupler

  ! IM coupling
  public:: GM_get_for_im_trace_crcm ! for CRCM
  public:: GM_get_for_im_crcm       ! for CRCM
  public:: GM_get_sat_for_im_crcm   ! for CRCM
  public:: GM_get_for_im_trace      ! for RAM
  public:: GM_get_for_im_line       ! for RAM
  public:: GM_get_for_im            ! for RCM
  public:: GM_satinit_for_im        ! initialize satellite
  public:: GM_get_sat_for_im        ! get satellite info
  public:: GM_get_multi_for_im      ! check if multifluid is used
  public:: GM_put_from_im           ! from IM

  ! PT coupling
  public:: GM_get_for_pt

  ! PW coupling
  public:: GM_get_for_pw
  public:: GM_put_from_pw

  ! RB coupling
  public:: GM_get_for_rb_trace
  public:: GM_get_for_rb
  public:: GM_satinit_for_rb
  public:: GM_get_sat_for_rb

  ! PC coupling
  public:: GM_get_for_pc_dt
  public:: GM_get_for_pc_init
  public:: GM_get_for_pc
  public:: GM_put_from_pc

contains
  !==========================================================================

  subroutine GM_set_param(CompInfo, TypeAction)

    use CON_comp_info
    use ModProcMH
    use ModIO, ONLY: iUnitOut, StringPrefix, STDOUT_, NamePlotDir
    use ModRestartFile, ONLY: NameRestartInDir, NameRestartOutDir
    use ModMain, ONLY : CodeVersion, NameThisComp, &
         time_accurate, StartTime, iStartTime_I, UseRotatingBc
    use ModB0, ONLY: DtUpdateB0, DoUpdateB0
    use CON_physics, ONLY: get_time, get_planet
    use ModTimeConvert, ONLY: time_real_to_int

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
            DtUpdateB0Out  = DtUpdateB0,    &
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
    use ModMain,ONLY:TypeCoordSystem, NameVarCouple
    use ModPhysics,ONLY:No2Si_V, UnitX_
    use ModVarIndexes,ONLY: nVar
    use CON_comp_param,ONLY:GM_

    !REVISION HISTORY:
    !23Aug03 I.Sokolov <igorsok@umich.edu> - initial prototype/prolog/code
    !03Sep03 G.Toth    <gtoth@umich.edu> - removed test_message_pass
    !                                      call synchronize_refinement directly
    !EOP
    logical ::DoTest,DoTestMe
    DoTest=.false.;DoTestMe=.false.
    if(done_dd_init(GM_))return
    call init_decomposition(GM_,GM_,3,.true.)
    call set_coord_system(GM_,TypeCoordSystem,No2Si_V(UnitX_), &
         NameVar = NameVarCouple, nVar=nVar)

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
  subroutine GM_get_grid_info(nDimOut, iGridOut, iDecompOut)

    use BATL_lib, ONLY: nDim
    use ModMain,  ONLY: iNewGrid, iNewDecomposition

    integer, intent(out):: nDimOut    ! grid dimensionality
    integer, intent(out):: iGridOut   ! grid index (increases with AMR)
    integer, intent(out):: iDecompOut ! decomposition index 

    character(len=*), parameter :: NameSub='GM_get_grid_info'

    ! Return basic grid information useful for model coupling.
    ! The decomposition index increases with load balance and AMR.
    !---------------------------------------------------------------------------

    nDimOut    = nDim
    iGridOut   = iNewGrid
    iDecompOut = iNewDecomposition

  end subroutine GM_get_grid_info
  !==============================================================================
  subroutine GM_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

    use BATL_lib,   ONLY: MaxDim, find_grid_block
    use ModPhysics, ONLY: Si2No_V, UnitX_

    integer, intent(in) :: nDimIn                ! dimension of position vectors
    integer, intent(in) :: nPoint                ! number of positions
    real,    intent(in) :: Xyz_DI(nDimIn,nPoint) ! positions
    integer, intent(out):: iProc_I(nPoint)       ! processor owning position

    ! Find array of points and return processor indexes owning them
    ! Could be generalized to return multiple processors...

    real:: Xyz_D(MaxDim) = 0.0
    integer:: iPoint, iBlock

    character(len=*), parameter:: NameSub = 'GM_find_points'
    !--------------------------------------------------------------------------
    do iPoint = 1, nPoint
       Xyz_D(1:nDimIn) = Xyz_DI(:,iPoint)*Si2No_V(UnitX_)
       call find_grid_block(Xyz_D, iProc_I(iPoint), iBlock)
    end do

  end subroutine GM_find_points
  !============================================================================
  subroutine GM_print_variables(NameSource)

    use ModMain, ONLY: NameThisComp
    use ModNumConst
    use ModImPressure
    use ModIoUnit, ONLY: UNITTMP_

    character(len=*), parameter :: NameSub='GM_print_variables'

    character(len=*),intent(in) :: NameSource
    integer            :: nFile=0
    character(len=100) :: NameFile
    character(len=100) :: NameVar
    integer            :: i,j
    !--------------------------------------------------------------------------

    if(NameThisComp/='GM') RETURN

    select case(NameSource)
    case('IM')
       NameVar='j i lon lat density pressure'
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
    case('IM')
       do i=1,iSize
          do j=1,jSize
             write(UNITTMP_,'(2i4,4G14.6)')j,i,IM_lon(j),IM_lat(i), &
                  IM_dens(i,j),IM_p(i,j)
          end do
       end do
    case('IE')
       call print_iono_potential
    end select
    close(UNITTMP_)

  end subroutine GM_print_variables

  !============================================================================

  subroutine GM_init_session(iSession, TimeSimulation)

    use ModProcMH,   ONLY: iProc
    use ModMain,     ONLY: Time_Simulation, UseIe, UsePw, TypeBC_I, body1_
    use ModMain,     ONLY: UseIM
    use CON_physics, ONLY: get_time
    use CON_coupler, ONLY: Couple_CC, IE_, IM_, GM_, IH_, PW_

    !INPUT PARAMETERS:
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='GM_init_session'

    logical :: IsUninitialized = .true.
    logical :: DoTest, DoTestMe
    !----------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest, DoTestMe)

    UseIm = Couple_CC(IM_,GM_) % DoThis
    UsePw = Couple_CC(PW_,GM_) % DoThis
    UseIe = Couple_CC(IE_,GM_) % DoThis

    ! Check if the boundary condition is properly set
    if(UsePw) TypeBC_I(body1_) = 'polarwind'

    if(Couple_CC(IH_,GM_) % DoThis .neqv. (TypeBc_I(2)=='coupled'))then
       if(Couple_CC(IH_,GM_) % DoThis) then
          TypeBc_I(2)='coupled'
       else
          if(iProc==0)write(*,*)NameSub,' WARNING: IH and GM are not coupled,',&
               ' changing west boundary type from "coupled" to "vary"'
          TypeBc_I(2)='vary'
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

  !============================================================================

  subroutine GM_finalize(TimeSimulation)

    use ModMain, ONLY: UseIe, time_loop

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='GM_finalize'

    !--------------------------------------------------------------------------
    ! We are not advancing in time any longer
    time_loop = .false.

    call BATS_save_files('FINAL')

    if (UseIe) call clean_mod_field_aligned_current

    call BATSRUS_finalize

  end subroutine GM_finalize

  !============================================================================

  subroutine GM_save_restart(TimeSimulation)

    use CON_coupler, ONLY: NameRestartOutDirComp
    use ModRestartFile, ONLY: NameRestartOutDir

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='GM_save_restart'
    !--------------------------------------------------------------------------
    if( NameRestartOutDirComp /= '') NameRestartOutDir = NameRestartOutDirComp

    call BATS_save_files('RESTART')

  end subroutine GM_save_restart

  !============================================================================

  subroutine GM_run(TimeSimulation,TimeSimulationLimit)

    use ModProcMH, ONLY: iProc
    use ModMain,   ONLY: Time_Simulation

    !INPUT/OUTPUT ARGUMENTS:
    real, intent(inout) :: TimeSimulation   ! current time of component

    !INPUT ARGUMENTS:
    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    character(len=*), parameter :: NameSub='GM_run'

    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(DoTestMe)write(*,*)NameSub,' called with tSim, tSimLimit=',&
         TimeSimulation, TimeSimulationLimit

    if(abs(Time_Simulation-TimeSimulation)>0.0001) then
       write(*,*)NameSub,' GM time=', Time_Simulation, &
            ' SWMF time=', TimeSimulation
       call CON_stop(NameSub// &
            ' SWMF_ERROR: GM and SWMF simulation times differ')
    end if

    call BATS_advance(TimeSimulationLimit)

    ! Return time after the time step
    TimeSimulation = Time_Simulation

    if(DoTestMe)write(*,*)NameSub,' finished with tSim=', TimeSimulation

  end subroutine GM_run

end module GM_wrapper
