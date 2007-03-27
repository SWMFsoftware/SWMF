subroutine GM_get_for_pw(nTotalLine, Buffer_I)
  
  use ModVarIndexes, ONLY: p_
  use ModPwGrid,     ONLY: CoordXyzPw_DI
  use ModMain,       ONLY: nBlock
  use ModProcMH,     ONLY: iComm
  use ModMpi
  use ModPhysics, ONLY: No2Si_V, UnitP_
  implicit none
  
  integer,intent(in)  :: nTotalLine
  real,   intent(out) :: Buffer_I(nTotalLine)
  real, allocatable   :: LocalBuffer_I(:)
  integer             :: iLine,iError
  !----------------------------------------------------------------------------

  call allocate(LocalBuffer_I(nTotalLine))

  !Get the pressure at each field line location.
  do iLine=1,nTotalLine
     call get_point_data(&
          0.0, CoordXyzPW_DI(:,iLine), 1, nBlock, p_, p_, LocalBuffer_I(iLine))
  enddo
  
  ! Sum contributions from all PE-s to processor 0
  call MPI_reduce(LocalBuffer_I, Buffer_I, nTotalLine,&
       MPI_REAL, MPI_SUM, 0, iComm, iError)
  
  ! Convert pressure to SI units
  Buffer_I = Buffer_I * No2Si_V(UnitP_)
   
  call deallocate(LocalBuffer_I)

 end subroutine GM_get_for_pw
