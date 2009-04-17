! File name: heidi_main.f90
!
! Contains: the main driver program for HEIDI
!	HEIDI
!
! Last Modified: March 2006, Mike Liemohn
!
!************************************************************************
! PROGRAM HEIDI - 
!   ***  Hot Electron and Ion Drift Integrator  ***
! A program to calculate the growth and decay of a given population of 
! ring current ions and electrons solving the bounce-averaged 
! Boltzmann equation considering drifts, charge exchange, atmospheric 
! loss, wave-particle interactions, electric and magnetic field effects
! and feedback, and Coulomb collisions.
!
! Numerical scheme for cons terms: Lax-Wendroff + superbee flux limiter
! Converted from mram05.f into heidi010.f90, March 2006, Mike Liemohn
!***********************************************************************

!.......NR=no. grids in radial direction, NT=no. grids in azimuth,
!	NE=no. of energy grids, NS=no. of species (e-, H+, He+, O+),
!	NPA=no. of grids in equatorial pitch angle

program heidi_main

  use ModHeidiIO 
  use ModInit
  use ModProcIM
  use ModHeidiMain
  implicit none 

  !****************************************************************************
  ! Initiallize MPI and get number of processors and rank of given processor
  !****************************************************************************

  !---------------------------------------------------------------------------
  call MPI_INIT(iError)
  iComm= MPI_COMM_WORLD

  call MPI_COMM_RANK(iComm, iProc, iError)
  call MPI_COMM_SIZE(iComm, nProc, iError)   

  !***************************************************************************
  ! Read the input file
  !***************************************************************************     

  IsFramework = .false.

  call heidi_read
  call heidi_check
  call IM_init_session(1, 0.0)

  do I3=NST,NSTEP
     call heidi_run
  end do			! End time loop

  call IM_finalize(0.0)

contains
  !==========================================================================
  subroutine IM_init_session(iSession, TimeSimulation)

    implicit none
    integer,  intent(in) :: iSession       ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time
    logical :: IsUninitialized = .true.
    !-----------------------------------------------------------------------

    if(IsUninitialized)then
       call heidi_init
       IsUninitialized = .false.
    end if
  end subroutine IM_init_session

  !==========================================================================  

  subroutine IM_finalize(TimeSimulation)

    use ModProcIM
    use ModInit, ONLY: nS

    implicit none
    real,     intent(in) :: TimeSimulation   ! seconds from start time
    integer:: iSpecies
    !-----------------------------------------------------------------------
    do iSpecies=1,NS
       CLOSE(15+iSpecies)          ! Closes continuous output file
    end do

    CLOSE(13)	            ! Closes sw1 input file
    CLOSE(15)		    ! Closes sw2 input file
    CLOSE(14)               ! Closes MPA input file
    CLOSE(16)               ! Closes SOPA input file
    CLOSE(18)               ! Closes FPOT input file

    call MPI_BARRIER(iComm,iError) 
    call MPI_finalize(iError)

  end subroutine IM_finalize
  !==========================================================================  

  subroutine IM_run(TimeSimulation,TimeSimulationLimit)

    implicit none
    real, intent(inout) :: TimeSimulation   ! current time of component
    real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded
    !-----------------------------------------------------------------------

    call heidi_run 

  end subroutine IM_run
  !==========================================================================  

end program heidi_main

subroutine CON_stop(StringError)
  
  use ModProcIM, ONLY: iProc, iComm
  use ModMpi
  
  implicit none
  character (len=*), intent(in) :: StringError
  integer :: iError,nError
  !-----------------------------------------------------------------------
  
  write(*,*)'Stopping execution! iProc=',iProc,' with msg:'
  write(*,*)StringError
  call MPI_abort(iComm, nError, iError)
  stop

end subroutine CON_stop

!==========================================================================

subroutine CON_set_do_test(String,DoTest,DoTestMe)
  
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
  !-----------------------------------------------------------------------
 
  DoTest = .false.; DoTestMe = .false.
  
end subroutine CON_set_do_test
!==========================================================================
subroutine CON_io_unit_new(iUnit)

  use ModIoUnit, ONLY: io_unit_new
  
  implicit none
  integer, intent(out) :: iUnit
  !-----------------------------------------------------------------------
  iUnit = io_unit_new()

end subroutine CON_io_unit_new

