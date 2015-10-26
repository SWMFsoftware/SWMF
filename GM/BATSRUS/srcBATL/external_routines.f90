! Collection of external routines called by BATL
!=============================================================================
subroutine CON_stop(String)
  use BATL_mpi, ONLY: iProc
  use ModMpi, ONLY: MPI_abort, MPI_COMM_WORLD
  implicit none
  integer:: iError, nError
  character (len=*), intent(in) :: String
  !--------------------------------------------------------------------------
  write(*,*)'CON_stop called on processor ',iProc,' with String='
  write(*,*) String

  call MPI_abort(MPI_COMM_WORLD, nError, iError)
  stop
end subroutine CON_stop
!=============================================================================
subroutine user_specify_refinement(iBlock, iArea, DoRefine)
  
  implicit none

  integer, intent(in) :: iBlock, iArea
  logical, intent(out):: DoRefine

  character(len=*), parameter :: NameSub = 'user_specify_refinement'
  !-------------------------------------------------------------------
  call CON_stop(NameSub//' is not implemented')

end subroutine user_specify_refinement
!=============================================================================
