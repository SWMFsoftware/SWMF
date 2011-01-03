module IE_ModIo

  ! Unit number for temporary files
  use ModIoUnit, ONLY: iUnit => UNITTMP_, STDOUT_

  implicit none
  save

  integer            :: iUnitOut=STDOUT_
  integer, parameter :: lStringPrefix=6
  character (len=lStringPrefix) :: StringPrefix = ''

  ! Maximum number of plot files
  integer, parameter :: MaxFile = 10

  ! Actual number of plot files
  integer :: nFile=0

  character (LEN=3)   :: plot_form(maxfile)
  character (LEN=100) :: plot_vars(maxfile)
  character (LEN=80)  :: filename
  character (len=100) :: NameIonoDir='ionosphere/'
  character (len=100) :: NameFile

  integer             :: dn_output(MaxFile)=-1, dn_magoutput = -1
  real                :: dt_output(MaxFile)=-1.0, dt_magoutput = -1.0
  integer             :: t_output_last(MaxFile)=-1, t_magoutput_last = -1
  
contains

  !===========================================================================
  subroutine write_prefix

    if(iUnitOut==STDOUT_)write(*,'(a)',ADVANCE='NO')trim(StringPrefix)

  end subroutine write_prefix


end module IE_ModIo
