!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! Wrapper for the empty PTOM (PT) component
!==========================================================================
subroutine PT_set_param(CompInfo, TypeAction)

  use CON_comp_info

  implicit none

  character (len=*), parameter :: NameSub='PT_set_param'
  integer :: iComm,iProc,nProc

  ! Arguments
  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do
  !-------------------------------------------------------------------------

!call AMPS_TimeStep()


  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,&
          Use        =.true., &
          NameVersion='AMPS', &
          Version    =1.0)

  case('MPI')
     call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
     call AMPS_SetMpiCommunicator(iComm, iProc, nProc)

  case('READ', 'CHECK')
     ! No input parameters

  case('STDOUT')
     ! call AMPS_set_stdout

  case('FILEOUT')
     ! call AMPS_set_fileout

  case('GRID')
     ! Grid info depends on BATSRUS

  case default
     call CON_stop(NameSub//': PT_ERROR: empty version cannot be used!')
  end select

end subroutine PT_set_param

!==============================================================================

subroutine PT_init_session(iSession, TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PT_init_session'

! call AMPS_init(iSession, TimeSimulation)

end subroutine PT_init_session

!==============================================================================

subroutine PT_finalize(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PT_finalize'

  ! call AMPS_finalize

end subroutine PT_finalize

!==============================================================================

subroutine PT_save_restart(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PT_save_restart'

  call CON_stop(NameSub//': PT_ERROR: empty version cannot be used!')

end subroutine PT_save_restart

!==============================================================================

subroutine PT_run(TimeSimulation,TimeSimulationLimit)

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='PT_run'

! call AMPS_run(TimeSimulation, TimeSimulationLimit)
call AMPS_TimeStep(TimeSimulation, TimeSimulationLimit) 

end subroutine PT_run

!==============================================================================
subroutine PT_get_grid_info(nDimOut, iGridOut, iDecompOut)

  ! Provide information about AMPS grid

  implicit none

  integer, intent(out):: nDimOut    ! grid dimensionality
  integer, intent(out):: iGridOut   ! grid index (increases with AMR)
  integer, intent(out):: iDecompOut ! decomposition index
  
  character(len=*), parameter :: NameSub = 'PT_get_grid_info'
  !---------------------------------------------------------------------------
  nDimOut    = 3
  iGridOut   = 1
  iDecompOut = 1

end subroutine PT_get_grid_info
!==============================================================================
subroutine PT_put_from_gm(UseData, &
     NameVar, nVar, nPoint, Pos_DI, Data_VI, iPoint_I)

  implicit none

  logical,          intent(in)   :: UseData ! true when data is transferred
                                            ! false if positions are asked
  character(len=*), intent(inout):: NameVar ! List of variables
  integer,          intent(inout):: nVar    ! Number of variables in Data_VI
  integer,          intent(inout):: nPoint  ! Number of points in Pos_DI

  real, pointer:: Pos_DI(:,:)               ! Position vectors

  real,    intent(in), optional:: Data_VI(nVar,nPoint)! Recv data array
  integer, intent(in), optional:: iPoint_I(nPoint)    ! Order of data

  integer:: iPoint, i, iProc=0

  character(len=*), parameter :: NameSub='PT_put_from_gm'
  !--------------------------------------------------------------------------
  if(.not.UseData)then
     ! Set variable names
     NameVar = 'Bx By Bz'
     ! Set number of variables needed
     nVar = 3

     ! set number of grid points on this processor
     call amps_get_center_point_number(nPoint)  !kkkkjjj10
     ! allocate array
     allocate(Pos_DI(3,nPoint))

     ! Fill in values
!    do iPoint = 1, nPoint
!       Pos_DI(:,iPoint) = (/ 10.0*iPoint, 20.0*iPoint, 30.0*iPoint/)
!    end do

     call amps_get_center_point_coordinates(Pos_DI) 

     RETURN
  end if

  write(*,*)NameSub,': iProc, iPoint, i, Pos, Data'


! do iPoint = 1, nPoint
!    i = iPoint_I(iPoint)
!    write(*,*)NameSub, iProc, iPoint, i, Data_VI(:,i)
!    ! Here should convert from SI to AMPS units and store data
! end do

  call amps_recieve_gm2amps_center_point_data(Data_VI,iPoint_I)

end subroutine PT_put_from_gm

