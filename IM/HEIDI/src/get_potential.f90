subroutine get_potential(iModel, imf_by, imf_bz, sw_v, &
     Year, Day, Hour, Minute, Second, &
     Lats1, MLTs1, OutPotential1, &
     Lats2, MLTs2, OutPotential2)

  use ModKind
  use ModErrors
  use ModHeidiSize

  implicit none

  integer, intent(in) :: iModel
  real, intent(in)    :: imf_by, imf_bz, sw_v
  real, intent(in)    :: Year, Day, Hour, Minute, Second

  real, dimension(NT, NR+3), intent(in)  :: Lats1, MLTs1
  real, dimension(NT, NR+3), intent(out) :: OutPotential1

  real, dimension(nphicells,nthetacells), intent(in)  :: Lats2, MLTs2
  real, dimension(nphicells,nthetacells), intent(out) :: OutPotential2

  integer                        :: iOutputError

  integer :: iError

  character (len=100), dimension(100) :: Lines

  real (kind=dblprec) :: rtime
  integer, dimension(7) :: itime

  logical :: IsFirstTime = .true.

  iOutputError = 0

  if (IsFirstTime) then

     Lines(1) = "#BACKGROUND"
     Lines(2) = "../src/"
     Lines(3) = "weimer96"
     Lines(4) = "ihp"
     Lines(5) = "idontknow"
     Lines(6) = ""
     Lines(7) = "#DEBUG"
     Lines(8) = "2"
     Lines(9) = "0"
     Lines(10) = ""
     Lines(11) = "#END"

     call IE_set_inputs(Lines)

     call IE_Initialize(iError)

     if (iError /=0) then
        write(*,*) "Error in IE_Initialize, called from get_potential"
        iOutputError = iError
        return
     endif

     IsFirstTime = .false.

  endif

  itime(1) = Year
  itime(2) = 1
  itime(3) = Day
  itime(4) = Hour
  itime(5) = Minute
  itime(6) = Second
  itime(7) = 0

  call time_int_to_real(itime, rtime)
  call IO_SetTime(rtime)

  call IO_SetIMFBz(imf_bz)
  call IO_SetIMFBy(imf_by)
  call IO_SetSWV  (abs(sw_v)/1000.0)

  call IO_SetnMLTs(NT)
  call IO_SetnLats(NR+3)

  call IO_SetGrid(MLTs1,Lats1, iError)

  call IO_GetPotential(OutPotential1, iError)

  if (iError /= 0) then
     write(*,*) "Error : ",cErrorCodes(iError)
     iOutputError = iError
     return
  endif

  call IO_SetnMLTs(nphicells)
  call IO_SetnLats(nthetacells)

  call IO_SetGrid(MLTs2,Lats2, iError)

  call IO_GetPotential(OutPotential2, iError)

  if (iError /= 0) then
     write(*,*) "Error : ",cErrorCodes(iError)
     iOutputError = iError
     return
  endif

end subroutine get_potential
