! Wrapper for DGCPM Plasmasphere Model
!==============================================================================
subroutine PS_set_param(CompInfo, TypeAction)

  use ModProcPS
  use ModMainDGCPM
  use ModHeidiDGCPM
  use ModIoDGCPM

  use ModIoUnit
  use ModReadParam, only: read_text
  use CON_comp_info

  implicit none

  character (len=*), parameter :: NameSub='PS_set_param'

  ! Arguments
  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do

  integer :: iError

  !-------------------------------------------------------------------------

  write(*,*) "typeaction : ",typeaction

  select case(TypeAction)

  case('VERSION')

     call put(CompInfo,&
          Use=.true.,                                    &
          NameVersion='DGCPM (Ober, Liemohn)',&
          Version=1.0)

  case('MPI')

     call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)

     if( nProc>1 )call CON_stop(NameSub//' PS_ERROR '//&
          'this version can run on 1 PE only!')

!     call check_dir("PS/Input")
!     call check_dir("PS/Output")
!     call check_dir("PS/RestartOUT")

     IsFramework = .true.

  case('READ','CHECK')

!     call DGCPM_read_inputs(cInputText)
     call read_text(cInputText)
     call readpara

  case('STDOUT')
     iOutputUnit_=STDOUT_
!     StringPrefix='PS:'
  case('FILEOUT')
     call get(CompInfo,iUnitOut=iOutputUnit_)
!     StringPrefix=''
  case('GRID')
     call PS_set_grid
  case default
     call CON_stop(NameSub//' PS_ERROR: invalid TypeAction='//TypeAction)
  end select

end subroutine PS_set_param

!=============================================================================
subroutine PS_set_grid

  ! Set the grid descriptor for PS
  ! Since PS has a static grid the descriptor has to be set once.
  ! There can be many couplers that attempt to set the descriptor,
  ! so we must check IsInitialized.
  use ModProcPS
  use ModSizeDGCPM
  use ModHeidiDGCPM
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

  write(*,*) "begin PS_set_grid"

  call set_grid_descriptor(                        &
       PS_,                                        &! component index
       nDim=2,                                     &! dimensionality
       nRootBlock_D=(/1,1/),                       &! radius and MLT
       nCell_D =(/nthetacells,nphicells/), &! size of node based grid
       XyzMin_D=(/cHalf, cHalf/),            &! min colat and longitude indexes
       XyzMax_D=(/nthetacells+cHalf,         &
       nphicells+cHalf/),                    &! max colat and longitude indexes
       TypeCoord='SMG',                            &! solar magnetic coord.
       Coord1_I=vthetacells,                       &! radius
       Coord2_I=vphicells,                         &! longitudes/MLTs
       Coord3_I=zTmp)

  write(*,*) "end PS_set_grid"

end subroutine PS_set_grid

!==============================================================================

subroutine PS_init_session(iSession, tSimulation)

  ! Initialize the Plasmasphere (PS) module for session iSession

  use CON_physics,   ONLY: get_time, get_planet, get_axes

  use ModIoDGCPM
  use ModMainDGCPM

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession      ! session number (starting from 1)
  real,     intent(in) :: tSimulation   ! seconds from start time

  !DESCRIPTION:
  ! Initialize the Plasmasphere (PS) module for session iSession

  character(len=*), parameter :: NameSub='PS_init_session'

  logical :: IsUninitialized=.true.

  logical :: DoTest,DoTestMe
  !---------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest,DoTestMe)

  if(IsUninitialized)then

     write(*,*) "PS_init_session"

     time = tSimulation
     t=time
     nst=nint(time/dt/2.) + 1
     nkp=nint(10800./dt/2.)
!     nibc=nint(tinj/dt/2.)

     write(*,*) "constant"
     call constant

     i2=(nst-1)/nkp + 1
     if (ikp >= 3) f107=f107r(i2)

     write(*,*) "arrays"
     call arrays
     !  call otherpara
     write(*,*) "thermal"
     call thermal   ! setup only

     !.......start the calculation

!     npr=nint(tint/dt/2.)
!     write(*,*) 'times:',nst,nstep,npr,nkp,i2,dt,nstep*2.*dt
     write(*,*) 'times:',nst,nstep,nkp,i2,dt,nstep*2.*dt

     IsUninitialized = .false.
  end if

  write(*,*) "Done with PS_init_session"

!  call get_time(  DoTimeAccurateOut = time_accurate)
!  call get_planet(DipoleStrengthOut = IONO_Bdp)
!  call get_axes(tSimulation, MagAxisTiltGsmOut = ThetaTilt)

!  IONO_Bdp = IONO_Bdp*1.0e9 ! Tesla -> nT

end subroutine PS_init_session

!==============================================================================
subroutine PS_finalize(tSimulation)

  use ModProcPS
  use CON_physics, ONLY: get_time
  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: tSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PS_finalize'

  !---------------------------------------------------------------------------

  call wresult(0)

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
  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: tSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: tSimulationLimit ! simulation time not to be exceeded

  real(Real8_) :: tStart
  integer      :: i3, npr

  character(len=*), parameter :: NameSub='PS_run'

  logical :: DoTest,DoTestMe
  !----------------------------------------------------------------------------

  call CON_set_do_test(NameSub,DoTest,DoTestMe)

  if(DoTest)write(*,*)NameSub,': iProc,tSimulation,tSimulationLimit=',&
       iProc,tSimulation,tSimulationLimit


  write(*,*) "PS_run"

  t = tSimulation
  i3 = tSimulation/(2.0*dt)

  write(*,*) "getkpa"
  call getkpa(i3)

  write(*,*) 'i3:',i3,t,kp

  write(*,*) "magconv"
  call magconv(i3)
  write(*,*) "thermal"
  call thermal

  write(*,*) "wresult"
  if (i3.eq.nst) call wresult(1)

  tSimulation = tSimulation+2.*dt

  npr=nint(tint/dt/2.)

  if (mod(i3,npr) == 0) then
     call wresult(0)
  end if

  return

end subroutine PS_run

!==============================================================================

subroutine PS_put_from_ie(Buffer_IIV, iSize, jSize, nVar)

  implicit none
  character (len=*),parameter :: NameSub='PS_put_from_ie'

  integer, intent(in)           :: iSize, jSize, nVar
  real, intent(out)             :: Buffer_IIV(iSize,jSize,nVar)

  !NOTE: The Buffer variables have been pushed to all PS processors already.

  write(*,*) NameSub,' -- called but not yet implemented.'

end subroutine PS_put_from_ie

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
  real, allocatable, intent(out) :: MLTs(:,:), Lats(:,:)

!FIX


end subroutine PS_get_grid

!==============================================================================
