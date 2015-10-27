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
subroutine user_block_inside_region(iArea, iBlock, nValue, NameLocation, &
     IsInside, IsInside_I, Value_I)

  implicit none

  integer,   intent(in):: iArea        ! area index in BATL_region
  integer,   intent(in):: iBlock       ! block index
  integer,   intent(in):: nValue       ! number of output values
  character, intent(in):: NameLocation ! c, g, x, y, z, or n

  logical, optional, intent(out) :: IsInside
  logical, optional, intent(out) :: IsInside_I(nValue)
  real,    optional, intent(out) :: Value_I(nValue)

  character(len=*), parameter :: NameSub = 'user_block_inside_region'
  !-------------------------------------------------------------------
  call CON_stop(NameSub//' is not implemented')

end subroutine user_block_inside_region
!=====================================================================
