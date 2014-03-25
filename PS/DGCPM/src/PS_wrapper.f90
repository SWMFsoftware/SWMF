! Wrapper for DGCPM Plasmasphere Model
!==============================================================================
subroutine PS_set_param(CompInfo, TypeAction)

  use ModProcPS
  use ModMainDGCPM
  use ModIoDGCPM
  use ModTimeDGCPM,  ONLY: MaxTime
  use CON_physics,   ONLY: get_time
  use CON_TIME,      ONLY: tSimulationMax
  use ModIoUnit
  use CON_comp_info

  implicit none

  character (len=*), parameter :: NameSub='PS_set_param'

  ! Arguments
  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do

  integer :: iError

  !-------------------------------------------------------------------------
    demo =2

  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,&
          Use=.true.,                                    &
          NameVersion='Dynamic Global Core Plasma Model (DGCPM)',&
          Version=0.1)
  case('MPI')
     call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
!     if( NameOutputDir(1:3) /= 'PS/' ) NameOutputDir = 'PS/'//NameOutputDir
  case('READ','CHECK')
     call read_param
  case('STDOUT')
     iUnitOut=STDOUT_
     if(nProc==1)then
        StringPrefix='PS:'
     else
        write(StringPrefix,'(a,i1,a)')'PS',iProc,':'
     end if
  case('FILEOUT')
     call get(CompInfo,iUnitOut=iUnitOut)
     StringPrefix=''
  case('GRID')
     call PS_set_grid
  case default
     call CON_stop(NameSub//' PS_ERROR: invalid TypeAction='//TypeAction)
  end select

contains
  subroutine read_param
    
    use ModReadParam
    use ModUtilities, Only: fix_dir_name, check_dir, lower_case
    ! The name of the command
    character (len=100) :: NameCommand
    
    ! Read parameters
    logical :: DoEcho=.false., UseStrict=.true.
    
    ! Plot file parameters
    integer :: iFile, i, iError, iDebugProc
    character (len=50) :: plot_string
    character (len=100) :: imffilename, kpfilename
    character (len=100), dimension(100) :: cTempLines
    real :: tmax
    
    !-------------------------------------------------------------------------
    select case(TypeAction)
    case('CHECK')
       ! We should check and correct parameters here
       if(iProc==0)write(*,*) NameSub,': CHECK iSession =',i_session_read()
       RETURN
       
    case('READ')
       if(iProc==0)write(*,*) NameSub,': READ iSession =',i_session_read(),&
            ' iLine=',i_line_read(),' nLine =',n_line_read()
    end select
    
    ! Read Input Data from Text via ModReadParam
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       
       select case(NameCommand)
       case("#TIMESTEP")
          call read_var('DtStep', DT)
       case("#TIMING")
          call read_var('DtStep',DT)
          call read_var('TMAX', Tmax)
       case("#KP")
          call read_var('IKP',IKP)
          if (ikp.eq.0) then
             call read_var('KP',KP)
          else if (ikp.eq.6) then
             cTempLines(1) = '#NGDC_INDICES'
             call read_var('kpFileName', kpFileName)
             cTempLines(2) = kpFileName
             cTempLines(3) = " "
             cTempLines(4) = "#END"
             
             call IO_set_Inputs(cTempLines)
             call read_NGDC_Indices(iError)
             
             if (iError /= 0) then
                write(*,*) "Read indices was NOT successful (NGDC KP File)"
                EXIT
             endif
          endif
          
       case("#NAME")
          call read_var('Name', Name)
       case("#SHUE")
          call read_var('UseShue', UseShue)
       case("#OUTPUT")
          call read_var('WriteStatic',WriteStatic)
          call read_var('WriteDynamic',WriteDynamic)
          call read_var('OutputInterval',TINT)
          call read_var('OutputType', OutputType)
          call read_var('MagneticType', MagneticType)
       case("#MLTSLICE")
          DoMltSlice=.true.
          call read_var('nMltSlice', nMltSlice)
          call read_var('DtMltSlice', DtMltSlice)
          
       case("#RESTART")
          call read_var('WriteRestart', WriteRestart)
       case("#LOG")
          call read_var('WriteLogFile', WriteLogFile)
       case("#FILLING")
          call read_var('EmptyPeriodClosed', EmptyPeriodClosed)
          call read_var('EmptyPeriodOpen', EmptyPeriodOpen)
          call read_var('FillRate', FillDays) 
          call read_var('FluxMax', FluxMax)
       case("#TESTS")
          call read_var('TestFill', TestFill)
       case default
          if(iProc==0) then
             write(*,'(a,i4,a)')NameSub//' IE_ERROR at line ',i_line_read(),&
                  ' invalid command '//trim(NameCommand)
             if(UseStrict)call CON_stop('Correct PARAM.in!')
          end if
       end select
    end do
    
    DT = DT * 2.
    !  Tmax = tSimulationMax
    
    ithermfirst=1		! So we do the setup routines in THERMAL
    
    LAMGAM=2.       ! Empirical E-field shielding parameter
    
  end subroutine read_param
  
end subroutine PS_set_param

!=============================================================================
subroutine PS_set_grid

  ! Set the grid descriptor for PS
  ! Since PS has a static grid the descriptor has to be set once.
  ! There can be many couplers that attempt to set the descriptor,
  ! so we must check IsInitialized.
  use ModProcPS
  use ModSizeDGCPM
  use ModMainDGCPM
  use CON_coupler

  implicit none
  character (len=*), parameter :: NameSub='PS_set_grid'
  logical :: IsInitialized=.false.
  real    :: zTmp(1)

  logical :: DoTest, DoTestMe

  !------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest, DoTestMe)
  if(DoTest)write(*,*)NameSub,' IsInitialized=',IsInitialized
  if(IsInitialized) return
  IsInitialized=.true.

  zTmp = 0.0

  if (debug .gt. 0) write(*,*) "begin PS_set_grid"

  call set_grid_descriptor(                        &
       PS_,                                        &! component index
       nDim=2,                                     &! dimensionality
       nRootBlock_D=(/1,1/),                       &! radius and MLT
       nCell_D =(/nthetacells,nphicells/),         &! size of node based grid
       XyzMin_D=(/cHalf, cHalf/),            &! min colat and longitude indexes
       XyzMax_D=(/nthetacells+cHalf,         &
       nphicells+cHalf/),                    &! max colat and longitude indexes
       TypeCoord='SMG',                            &! solar magnetic coord.
       Coord1_I=90.0-vthetacells,                       &! radius
       Coord2_I=vphicells,                         &! longitudes/MLTs
       Coord3_I=zTmp)

  if (debug .gt. 0) write(*,*) "end PS_set_grid"

end subroutine PS_set_grid

!==============================================================================
subroutine PS_get_for_ie(nPoint,iPointStart,Index,Weight,Buff_V,nVar)

  ! Provide current for IE
  ! The value should be interpolated from nPoints with
  ! indexes stored in Index and weights stored in Weight
  ! The variables should be put into Buff_V(??)

  use CON_router,   ONLY: IndexPtrType, WeightPtrType
  use ModIonoDGCPM, ONLY: IONO_NORTH_RCM_JR,IONO_SOUTH_RCM_JR, IONO_nTheta, IONO_nPsi

  implicit none
  character(len=*), parameter :: NameSub='PS_get_for_ie'

  integer,intent(in)            :: nPoint, iPointStart, nVar
  real,intent(out)              :: Buff_V(nVar)
  type(IndexPtrType),intent(in) :: Index
  type(WeightPtrType),intent(in):: Weight

  integer :: iLat, iLon, iBlock, iPoint
  real    :: w

  !---------------------------------------------------------------------------
  Buff_V = 0.0

  do iPoint = iPointStart, iPointStart + nPoint - 1

     iLat   = Index % iCB_II(1,iPoint)
     iLon   = Index % iCB_II(2,iPoint)
     iBlock = Index % iCB_II(3,iPoint)
     w      = Weight % Weight_I(iPoint)

     if(iBlock/=1)then
        write(*,*)NameSub,': iPoint,Index % iCB_II=',&
             iPoint,Index%iCB_II(:,iPoint)
        call CON_stop(NameSub//&
             ' SWMF_ERROR iBlock should be 1=North in PS-IE coupling')
     end if

     if(iLat<1 .or. iLat>IONO_nTheta*2 .or. iLon<1 .or. iLon>IONO_nPsi+1)then
        write(*,*)'iLat,iLon=',iLat, IONO_nTheta*2, iLon, IONO_nPsi
        call CON_stop(NameSub//' SWMF_ERROR index out of range')
     end if

     ! Only worry about the northern hemisphere....  IE can fix the southern hemisphere.
     if (iLat <= IONO_nTheta .and. iLon <= IONO_nPsi) &
          Buff_V(1) = Buff_V(1) + w * IONO_NORTH_RCM_JR(iLat,iLon)

     if (iLat > IONO_nTheta .and. iLon <= IONO_nPsi) &
          Buff_V(1) = Buff_V(1) + w * IONO_SOUTH_RCM_JR(2*IONO_nTheta-iLat+1,iLon)

  end do

end subroutine PS_get_for_ie

!============================================================================
subroutine PS_put_from_ie(ntheta, nPhi, Potential_Out, FieldModel)

  use ModCoupleDGCPM
  use ModMainDGCPM
  
  implicit none
  
  integer, intent(in):: nTheta, nPhi
  real,    intent(in):: Potential_out(nTheta, nPhi)

  integer :: i,j
  character(len=100):: NameFile 
  character(len=100):: FieldModel

  !-------------------------------------------------------------------------
    coupled_potential = Potential_out
!    EFieldModel = FieldModel
    isCoupled = .true.
   
       RETURN

end subroutine PS_put_from_ie

!==============================================================================
subroutine PS_put_from_ie_complete

  implicit none
  
  !--------------------------------------------------------------------------
  
  write(*,*) "Don't know what this is really supposed to do.  I think that it is"
  write(*,*) "Supposed to be applying periodic boundaries...?"

end subroutine PS_put_from_ie_complete

!==============================================================================

subroutine PS_init_session(iSession, tSimulation)

  ! Initialize the Plasmasphere (PS) module for session iSession

  use CON_physics,   ONLY: get_time, get_planet, get_axes
  use CON_time,   ONLY: tSimulationMax
  use ModIoDGCPM
  use ModMainDGCPM
  use ModTimeDGCPM

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession      ! session number (starting from 1)
  real,     intent(in) :: tSimulation   ! seconds from start time
  
  !DESCRIPTION:
  ! Initialize the Plasmasphere (PS) module for session iSession

  ! Local variables:
  real :: Tmax = 0.0
  
  character(len=*), parameter :: NameSub='PS_init_session'
  logical :: DoTest,DoTestMe

  IsUninitialized = .true.
  !---------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest,DoTestMe)

  if(IsUninitialized)then

    if (debug .gt. 0) write(*,*) "PS_init_session"

     ! Use the SWMF time max, set max number of steps.
     Tmax = tSimulationMax
     NSTEP=NINT(TMAX/DT/2.)
     ! Setup Time variables
     time = tSimulation
     t=time
     nst=nint(time/dt/2.) + 1
     nkp=nint(10800./dt/2.)
     call get_time(tCurrentOut = StartTime)

     ! Read KP Input File via GETKPA
     if (debug .gt. 0) write(*,*) "GetKPA"
     call GetKPA(0)

     i2=(nst-1)/nkp + 1
     if (ikp >= 3) f107=f107r(i2)

     if (debug.gt.0) write(*,*) "thermal"
     call thermal   ! setup only

     if (TestFill.gt.0) call TestFilling(1.0)

     ! Finish Initialization     
     IsUninitialized = .false.
     WriteStatic = .true.
  end if

  if (debug .gt. 0) write(*,*) "Done with PS_init_session"

end subroutine PS_init_session

!==============================================================================
subroutine PS_finalize(tSimulation)

  use ModProcPS
  use CON_physics, ONLY: get_time
  use ModIoDGCPM,  ONLY: iUnitSlice, DoMltSlice, iUnitMlt, nMltSlice
  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: tSimulation   ! seconds from start time

  ! Other variables:
  integer :: i

  character(len=*), parameter :: NameSub='PS_finalize'

  !---------------------------------------------------------------------------

  call wresult(0)

  ! Close files:
  close(iUnitSlice)
  if(DoMltSlice) then
     do i=1, nMltSlice
        close(iUnitMlt(i))
     end do
     deallocate(iUnitMlt)
  end if

end subroutine PS_finalize

!==============================================================================

subroutine PS_save_restart(tSimulation)

  use ModProcPS, ONLY:  nProc

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: tSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PS_save_restart'

  ! I don't know how to restart...

  RETURN

end subroutine PS_save_restart

!==============================================================================

subroutine PS_run(tSimulation,tSimulationLimit)

  use ModProcPS
  use ModMainDGCPM
  use ModIoDGCPM
  use CON_physics, ONLY: get_time, get_axes, time_real_to_int
  use ModKind
  use ModTimeDGCPM
  use ModCoupleDGCPM

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: tSimulation   ! current time of component
  
  !INPUT ARGUMENTS:
  real, intent(in) :: tSimulationLimit ! simulation time not to be exceeded

  real(Real8_) :: tStart
  integer      :: i3
  real         :: dt_requested

  character(len=*), parameter :: NameSub='PS_run'

  logical :: DoTest,DoTestMe
  !----------------------------------------------------------------------------

  call CON_set_do_test(NameSub,DoTest,DoTestMe)

  if(DoTest)write(*,*)NameSub,': iProc,tSimulation,tSimulationLimit=',&
       iProc,tSimulation,tSimulationLimit

!  if (debug .gt. 0) write(*,*) "PS_run"
 
  dt_requested = tsimulationlimit-tsimulation
  if (dt_requested < dt*2) dt = dt_requested/2
 
  CurrentTime = StartTime + tSimulation
    

  t = tSimulation
  i3 = tSimulation/(2.0*dt)

  nst=nint(t/dt/2.) + 1
  nkp=nint(10800./dt/2.)
  i2=(nst-1)/nkp + 1

!  if (debug .gt. 0) write(*,*) "getkpa"
  call getkpa(i3)

!  if (debug .gt. 0) write(*,*) "magconv"
  if (.not.(isCoupled)) call magconv()
!  if (debug .gt. 0) write(*,*) "thermal"
  call thermal

!  if (debug .gt. 0) write(*,*) "wresult"
  if (i3.eq.nst) call wresult(1)

  ! Update timing.
  tSimulation = tSimulation+2.*dt
  CurrentTime = StartTime + tSimulation
  
  ! Log File Writing
  if (WriteLogFile .and. (mod(tSimulation, 300.0)<0.001)) then
     call LogFileDGCPM(cOutputDir, i3)
  endif


  ! General Output Writing
  if (mod(tSimulation, tint) < 1E-5) then
     call wresult(0)
  end if

  ! Slice file writing.
  if (mod(tSimulation, 300.0)<0.001) call write_lslice
  if (DoMltSlice .and. mod(tSimulation,DtMltSlice)<0.001) call write_mltslice
  return

end subroutine PS_run

!==============================================================================

!subroutine PS_put_from_ie(Buffer_IIV, iSize, jSize, nVar)

!  implicit none
!  character (len=*),parameter :: NameSub='PS_put_from_ie'

!  integer, intent(in)           :: iSize, jSize, nVar
!  real, intent(out)             :: Buffer_IIV(iSize,jSize,nVar)

  !NOTE: The Buffer variables have been pushed to all PS processors already.

!  write(*,*) NameSub,' -- called but not yet implemented.'

!end subroutine PS_put_from_ie

!=================================================================

subroutine PS_get_grid_size(iSize,jSize)

  implicit none
  integer, intent(out) :: iSize, jSize

!FIX
  iSize=1
  jSize=1

end subroutine PS_get_grid_size

!==============================================================================

subroutine PS_get_grid(MLTs,Lats)

  implicit none
  real :: MLTs(:,:), Lats(:,:)

!FIX


end subroutine PS_get_grid

!==============================================================================
