module BATL_mpi

  use ModMpi

  implicit none

  SAVE

  private ! except

  public:: init_mpi       ! initialize this module
  public:: clean_mpi      ! finalize MPI
  public:: barrier_mpi    ! use an MPI barrier to synchronize processors

  integer, public:: iComm ! MPI communicator for the group of processors
  integer, public:: nProc ! number of processors in this group
  integer, public:: iProc ! processor rank from 0 to nProc-1

contains

  !==========================================================================
  subroutine init_mpi(iCommIn)

    ! Initialize iComm, nProc and iProc. If iCommIn is not present, set it
    ! to MPI_COMM_WORLD and also call MPI_init

    integer, optional, intent(in):: iCommIn
    integer :: iError
    !------------------------------------------------------------------------
    if(.not.present(iCommIn))then
       call MPI_init(iError)
       iComm = MPI_COMM_WORLD
    else
       iComm = iCommIn
    end if
    call MPI_COMM_RANK (iComm, iProc, iError)
    call MPI_COMM_SIZE (iComm, nProc, iError)

  end subroutine init_mpi
  !==========================================================================
  subroutine clean_mpi

    ! This should only be called if the whole application is finished

    integer :: iError
    !------------------------------------------------------------------------
    call MPI_finalize(iError)

  end subroutine clean_mpi
  !==========================================================================
  subroutine barrier_mpi

    use ModUtilities, ONLY: flush_unit
    use ModIoUnit,    ONLY: STDOUT_

    integer:: iError
    !-----------------------------------------------------------------------
    call flush_unit(STDOUT_)
    call MPI_barrier(iComm, iError)

  end subroutine barrier_mpi

end module BATL_mpi
