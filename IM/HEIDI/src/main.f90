
program GetIE

  use ModKind
  use ModErrors

  implicit none

  character (len=100), dimension(100) :: Lines
  integer :: iError, i
  integer, dimension(7) :: itime
  real (kind=dblprec) :: rtime
  real :: value
  real, dimension(25,1) :: mlts, lats, TempPotential

  Lines(1) = "#BACKGROUND"
  Lines(2) = "../srcData/"
  Lines(3) = "izmem"
  Lines(4) = "ihp"
  Lines(5) = "idontknow"
  Lines(6) = ""
  Lines(7) = "#DEBUG"
  Lines(8) = "2"
  Lines(9) = "0"
  Lines(10) = ""
  Lines(11) = "#END"

  write(*,*) "Calling IE_set_inputs"

  call IE_set_inputs(Lines)

  itime(1) = 1998
  itime(2) = 5
  itime(3) = 1
  itime(4) = 21
  itime(5) = 30
  itime(6) = 0
  itime(7) = 0

  write(*,*) "Calling time convertion"

  call time_int_to_real(itime, rtime)

  write(*,*) "=> Setting up Ionospheric Electrodynamics"

  call IE_Initialize(iError)

  write(*,*) "Setting nmlts and nlats"

  call IO_SetnMLTs(25)
  call IO_SetnLats(1)

  call IO_SetTime(rtime)

  call IO_SetIMFBz( -5.0)
  call IO_SetIMFBy( 10.0)
  call IO_SetSWV  (400.0)
  call IO_SetKp   (  3.0)
  call IO_SetNorth

  lats = 75.0
  do i=1,25
     mlts(i,1) = float(i)
  enddo

  write(*,*) "Setting grid"

  call IO_SetGrid(mlts,lats, iError)

  call IO_GetPotential(TempPotential, iError)

  if (iError /= 0) then
     write(*,*) "Error : ",cErrorCodes(iError)
  else
     write(*,*) TempPotential
  endif

  call IE_End(iError)

end program GetIE
