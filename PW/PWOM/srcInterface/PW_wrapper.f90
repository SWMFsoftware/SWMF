!^CFG COPYRIGHT UM
! Wrapper for the empty PWOM (PW) component
!==========================================================================
subroutine PW_set_param(CompInfo, TypeAction)

  use CON_comp_info
  use ModIoUnit, only: STDOUT_
  use Mod_PW, only: IsFramework, iUnitOut, iProc, nProc, iComm, StringPrefix

  implicit none

  character (len=*), parameter :: NameSub='PW_set_param'

  ! Arguments
  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do
  !-------------------------------------------------------------------------
  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,&
          Use        =.true., &
          NameVersion='PWOM (A. Glocer et al.)', &
          Version    =1.0)
  case('MPI')
     call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
     IsFramework = .true.
  case('READ')
     call PW_set_parameters('READ')

  case('CHECK')
     call PW_set_parameters('CHECK')

  case('STDOUT')

     iUnitOut=STDOUT_
     if(nProc==1)then
        StringPrefix='PW:'
     else
        write(StringPrefix,'(a,i3.3,a)')'PW',iProc,':'
     end if

  case('FILEOUT')

     call get(CompInfo,iUnitOut=iUnitOut)
     StringPrefix=''

  case('GRID')
     ! Do nothing
  case default
     call CON_stop(NameSub//': PW_ERROR: empty version cannot be used!')
  end select

end subroutine PW_set_param

!==============================================================================

subroutine PW_init_session(iSession, TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PW_init_session'

  call CON_stop(NameSub//': PW_ERROR: empty version cannot be used!')

end subroutine PW_init_session

!==============================================================================

subroutine PW_finalize(TimeSimulation)

  use ModPwom, ONLY: nLine, iUnitGraphics

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PW_finalize'

  integer :: iLine

  do iLine=1,nLine
     CLOSE(UNIT=iUnitGraphics(iLine))
  enddo

end subroutine PW_finalize

!==============================================================================

subroutine PW_save_restart(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PW_save_restart'

  call CON_stop(NameSub//': PW_ERROR: empty version cannot be used!')

end subroutine PW_save_restart

!==============================================================================

subroutine PW_run(TimeSimulation,TimeSimulationLimit)

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='PW_run'

  call CON_stop(NameSub//': PW_ERROR: empty version cannot be used!')

end subroutine PW_run

