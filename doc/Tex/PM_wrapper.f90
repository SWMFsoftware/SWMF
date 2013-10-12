!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP =========================================================================
!ROUTINE: PM_set_param - set parameters for the PM module
!INTERFACE:
subroutine PM_set_param(CompInfo, TypeAction)

  !USES:
  use CON_comp_info, ONLY: CompInfoType, get, put   ! component info, access
  use CON_coupler,   ONLY: PM_, set_grid_descriptor ! component ID, descriptor
  use ModConst,      ONLY: cPi                      ! Pi = 3.1415...
  use ModIoUnit,     ONLY: STDOUT_                  ! unit for STDOUT

  use PM_ModProc, ONLY: PM_iComm, PM_iProc, PM_nProc ! MPI parameters
  use PM_ModMain, ONLY: PM_Dt, PM_TimeAccurate       ! Some variables to set
  use PM_ModSize, ONLY: PM_n, PM_m                   ! Grid size
  use PM_ModIo,   ONLY: PM_iUnit, PM_Prefix  ! I/O unit and prefix for output

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp

  !INPUT ARGUMENTS:
  character (len=*), intent(in)     :: TypeAction ! What to do

  !LOCAL VARIABLES:
  character (len=*), parameter :: NameSub = 'PM_set_param'

  logical :: DoTest, DoTestMe

  !DESCRIPTION:
  ! This subroutine serves multiple purposes. It is called multiple times
  ! with various values for the TypeAction argument. The CompInfo argument
  ! is a derived type that contains basic information about this component.
  ! The data in this object should be accessed with the get and put methods
  ! of the CON\_comp\_info class. This subroutine is first called with
  ! TypeAction='VERSION' so that the version name is provided by PM to CON.
  ! Next it is called with TypeAction='MPI' and CON provides the MPI 
  ! communication group parameters to PM. The next actions are 'STDOUT' and
  ! possibly 'FILEOUT', which set the prefix and the I/O unit for 
  ! the verbose output of PM. The 'STDOUT' and 'FILEOUT' actions may be 
  ! called multiple times during the run.
  ! In each session the physics module reads its input parameters
  ! (TypeAction='READ') and checks the parameters ('CHECK').
  ! The action 'READ' may be called multiple times before the 'CHECK' action.
  ! After the module is initialized for the session (with the 
  ! PM\_init\_session method), the physics module provides a 
  ! description of its grid to CON (TypeAction='GRID').
  ! 
  !EOP
  !---------------------------------------------------------------------------
  !BOC
  ! Set DoTest and DoTesMe logicals. Check for the name of this subroutine
  ! in the StringTest parameter of the #TEST command in PARAM.in.
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  if(DoTest)write(*,*)NameSub,' called with TypeAction, iProc=',&
       TypeAction, PM_iProc

  select case(TypeAction)
  case('VERSION')
     ! Provide module name and version number for CON
     call put(CompInfo,                          &
          Use        = .true.,                   &
          NameVersion= 'CODE NAME (developers)', &
          Version    = 1.4 )

  case('MPI')
     ! Set the MPI communicator, the processor rank, and number of PE-s
     call get(CompInfo, iComm=PM_iComm, iProc=PM_iProc, PM_nProc=PM_nProc)

  case('STDOUT')
     ! Set PM_iUnit to STDOUT and prefix to 'PMx:' where x is the PE rank
     PM_iUnit=STDOUT_
     write(PM_Prefix,'(a,i1,a)') 'PM', PM_iProc, ':'

  case('FILEOUT')
     ! Set PM_iUnit to the value provided by CON. No need for prefix.
     call get(CompInfo, iUnitOut=PM_iUnit)
     PM_Prefix=''

  case('READ')
     ! Read input parameters for this component
     call read_param

  case('CHECK')
     ! Get some basic parameters from CON
     call get_time(DoTimeAccurateOut = PM_TimeAccurate)

     ! Check interdependent parameters here for consistency

  case('GRID')
     ! Provide grid description for the coupling toolkit.
     ! The grid is a uniform spherical grid divided into two blocks
     ! along the Theta coordinate: northern and southern hemispheres.
     ! Both hemispheres have a PM_n by PM_m grid.
     ! The module can run on 1 or 2 PE-s,
     ! so the processor array is (/0,0/) or (/0,1/)

     call set_grid_descriptor(          &
          PM_,                          &! component index
          nDim        =2,               &! dimensionality of the grid
          nRootBlock_D=(/2,1/),         &! block array for the 2 hemispheres
          nCell_D     =(/PM_n, PM_m/),  &! grid size for each block
          XyzMin_D    =(/0., 0./),      &! min colatitude and longitude indexes
          XyzMax_D    =(/cPi, 2*cPi/),  &! max colatitude and longitude indexes
          TypeCoord   ='SMG',           &! solar magnetic coordinate system
          iProc_A     =(/0,PM_nProc-1/)) ! processor assigment for the blocks

  case default
     call CON_stop(NameSub//': invalid TypeAction='//TypeAction)
  end select

contains
  !========================================================================
  subroutine read_param
    ! Read parameter via ModReadParam
    use ModReadParam, ONLY: &
         read_line, read_command, read_var, i_line_read, lStringLine

    ! String for the name of the command
    character (len=lStringLine) :: NameCommand
    !---------------------------------------------------------------------
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
       case("#TIMESTEP")
          call read_var('Dt',PM_Dt)
       ! Read other component parameters
       case default
          if(PM_iProc==0) then
             write(*,'(a,i4,a)')NameSub//': ERROR at line ',i_line_read(),&
                  ' invalid command '//trim(NameCommand)
             if(UseStrict)call CON_stop('Correct PARAM.in!')
          end if
       end select
    end do

  end subroutine read_param
  !EOC
end subroutine PM_set_param

!BOP =========================================================================
!ROUTINE: PM_init_session - initialize PM at the beginning of the session
!INTERFACE:
subroutine PM_init_session(iSession, TimeSimulation)
  !USES:
  use PM_ModIo, ONLY: PM_iUnit, PM_Prefix ! I/O unit and prefix for output

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession       ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation ! seconds from start time

  !LOCAL VARIABLES:
  character(len=*), parameter :: NameSub = 'PM_init_session'

  logical :: DoInitialize = .true.
  !EOP
  !---------------------------------------------------------------------------
  !BOC

  ! Initialize module for the first time
  if(DoInitialize)then
     call PM_setup
     DoInitialize = .false.
  end if

  ! Initialization for all sessions comes here

  ! Here is an example of using PM_iUnit and PM_Prefix
  write(PM_iUnit,*) PM_Prefix, NameSub, ' finished for session ', iSession

  !EOC
end subroutine PM_init_session

!BOP ==========================================================================
!ROUTINE: PM_run - do one time step with PM
!INTERFACE:
subroutine PM_run(TimeSimulation, TimeSimulationLimit)

  !USES:
  use PM_ModMain, ONLY: PM_Time ! Simulation time of PM

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  !LOCAL VARIABLES:
  character(len=*), parameter :: NameSub = 'PM_run'
  real                        :: Dt                 ! Time step
  !EOP
  !---------------------------------------------------------------------------
  !BOC
  ! Check if simulation times agree (e.g. for restart)
  if(abs(PM_Time-TimeSimulation)>0.0001) then
     write(*,*)NameSub,' PM time=',PM_Time,' SWMF time=',TimeSimulation
     call CON_stop(NameSub//' ERROR: PM and SWMF simulation times differ')
  end if

  ! Limit time step
  Dt = min(PM_Dt, TimeSimulationLimit - TimeSimulation)

  ! Call the appropriate subroutine of the physics module to advance by Dt
  call PM_advance(Dt)

  ! Return the new simulation time after the time step
  TimeSimulation = PM_Time

  !EOC
end subroutine PM_run

!BOP ==========================================================================
!ROUTINE: PM_save_restart - save restart files for PM
!INTERFACE:
subroutine PM_save_restart(TimeSimulation)
  implicit none
  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time
  !EOP
  !---------------------------------------------------------------------------
  !BOC
  ! Call the appropriate subroutine of the physics module
  call PM_save_restart_files(TimeSimulation)
  !EOC
end subroutine PM_save_restart

!BOP ==========================================================================
!ROUTINE: PM_finalize - finalize PM at the end of the run
!INTERFACE:
subroutine PM_finalize(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time
  !EOP
  !---------------------------------------------------------------------------
  !BOC
  ! Call the appropriate subroutines of the physics module
  call PM_save_files_final(TimeSimulation)
  call PM_error_report
  !EOC
end subroutine PM_finalize
