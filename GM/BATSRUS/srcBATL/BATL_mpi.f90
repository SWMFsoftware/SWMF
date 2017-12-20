!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_mpi

  use ModMpi

  implicit none

  SAVE

  private ! except

  public:: init_mpi       ! initialize this module
  public:: clean_mpi      ! finalize MPI
  public:: barrier_mpi    ! use an MPI barrier to synchronize processors

  integer, public:: iComm       ! MPI communicator for the group of processors
  integer, public:: nProc = 1   ! number of processors in this group
  integer, public:: iProc = 0   ! processor rank from 0 to nProc-1
  integer, public:: nThread = 1 ! number of threads per MPI process
  integer, public:: iThread = 0 ! thread rank from 0 to nThread-1

contains
  !============================================================================
  subroutine init_mpi(iCommIn)

    ! Initialize iComm, nProc and iProc. If iCommIn is not present, set it
    ! to MPI_COMM_WORLD and also call MPI_init
    
    integer, optional, intent(in):: iCommIn
    integer :: iError

    ! OpenMP support levels
    integer, parameter :: lSupportRequired = MPI_THREAD_SINGLE
    integer            :: lSupportProvided = MPI_THREAD_SINGLE 
    !--------------------------------------------------------------------------
    if(.not.present(iCommIn))then
       iComm = MPI_COMM_WORLD
       call MPI_init_thread(lSupportRequired, lSupportProvided, iError)
       ! Check the threading support level
       if (lSupportProvided < lSupportRequired) then
          if (iProc == 0) write(*,*) &
               "Warning:  This MPI implementation provides ",   &
               "insufficient threading support. Switching to pure MPI..."
          !$ call omp_set_num_threads(1)
       end if
    else
       iComm = iCommIn
    end if
    call MPI_comm_rank(iComm, iProc, iError)
    call MPI_comm_size(iComm, nProc, iError)
    
  end subroutine init_mpi
  !============================================================================

  subroutine clean_mpi

    ! This should only be called if the whole application is finished

    integer :: iError
    !--------------------------------------------------------------------------
    call MPI_finalize(iError)

  end subroutine clean_mpi
  !============================================================================
  subroutine barrier_mpi

    use ModUtilities, ONLY: flush_unit
    use ModIoUnit,    ONLY: STDOUT_

    integer:: iError
    !--------------------------------------------------------------------------
    call flush_unit(STDOUT_)
    call MPI_barrier(iComm, iError)

  end subroutine barrier_mpi
  !============================================================================

end module BATL_mpi
!==============================================================================
