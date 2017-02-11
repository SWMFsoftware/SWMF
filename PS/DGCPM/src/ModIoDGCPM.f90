!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
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
  integer :: nstep,iwpi
  real    :: ut,r,tint,time

  integer :: nst, nkp, nibc
  
  ! Define variables for setting Kp:
  character(len=5) :: NameSourceKp='const'
  real             :: kpConst=0.0, Kp=0.0
  
  ! Define some variables set with the input file
  ! Formerly: Common block PARAM3
  real :: TINJ
  character(len=5) ::  name

  ! Define convection input parameters
  ! Formerly: Common block PARAM5
  integer, parameter :: ilame=450, ippcm=5000
  real :: lamgam,lambe(ilame),tlame(ilame),ppc(ippcm),tppc(ippcm)

  integer, parameter        :: iCharLen_     = 100

  integer, parameter        :: nInputMaxLines = 10000
  integer                   :: nInputLines
  character (len=iCharLen_) :: cInputText(nInputMaxLines)

  character (len=iCharLen_) :: cInputFile = "input.dgcpm"
  character(len=*), parameter :: cOutputDir  = "PS/Output/"
  character(len=*), parameter :: cInputDir   = "PS/Input/"
  character(len=*), parameter :: cRestartIn  = "PS/restartIN/"
  character(len=*), parameter :: cRestartOut = "PS/restartOUT/"
  

  integer :: iUnitOut=STDOUT_
  integer, parameter :: lStringPrefix=6
  character (len=lStringPrefix) :: StringPrefix = '' 

  ! MLTSlice file parameters.
  logical :: DoMltSlice = .false.
  integer :: nMltSlice = 4
  real    :: DtMltSlice = 300.0
  integer, allocatable :: iUnitMlt(:)

  ! LSlice file parameters.
  integer :: iUnitSlice

contains

  !===========================================================================
 
  subroutine write_prefix

    if(iUnitOut==STDOUT_)write(*,'(a)',ADVANCE='NO')trim(StringPrefix)

  end subroutine write_prefix

end Module ModIoDGCPM
