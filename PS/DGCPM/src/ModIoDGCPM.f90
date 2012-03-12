Module ModIoDGCPM
  !\
  ! Input/output variable definition module for the DGCPM program.
  ! Updated: Aron Dodger, January 2012 
  !/

  use ModSizeDGCPM
  use ModIoUnit, Only: iUnit => UNITTMP_, STDOUT_

  logical :: IsFramework

  ! Define a few time and geophysical index input variables
  ! Formerly: Common block PARAM2
  integer nstep,ikp,iwpi
  real :: ut,r,kp=0.0,f107,tint,time

  integer :: nst, nkp, nibc, i2

  ! Define some variables set with the input file
  ! Formerly: Common block PARAM3
  real TINJ
  character*5 name

  ! Define convection input parameters
  ! Formerly: Common block PARAM5
  integer ilame,ilambe,ippcm,ippc
  parameter (ilame=450, ippcm=5000)
  real lamgam,lambe(ilame),tlame(ilame),ppc(ippcm),tppc(ippcm)

  !----------------------------------------------------------------
  ! Stuff added by Aaron
  !----------------------------------------------------------------

  integer, parameter        :: iCharLen_     = 100

  integer                   :: iOutputUnit_  != UnitTmp_
  integer                   :: iInputUnit_   != UnitTmp_
  integer                   :: iRestartUnit_ != UnitTmp_

  integer, parameter        :: nInputMaxLines = 10000
  integer                   :: nInputLines
  character (len=iCharLen_) :: cInputText(nInputMaxLines)

  character (len=iCharLen_) :: cInputFile = "input.dgcpm"
  character (len=10) :: cOutputDir = "PS/Output/"
  character (len= 9) :: cInputDir  = "PS/Input/"

  integer :: iUnitOut=STDOUT_
  integer, parameter :: lStringPrefix=6
  character (len=lStringPrefix) :: StringPrefix = '' 

contains

  !===========================================================================
 
  subroutine write_prefix

    if(iUnitOut==STDOUT_)write(*,'(a)',ADVANCE='NO')trim(StringPrefix)

  end subroutine write_prefix

end Module ModIoDGCPM
