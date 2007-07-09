Module ModIoDGCPM
  !\
  ! Input/output variable definition module for the HEIDI program.
  ! Mike Liemohn, March 2006
  !/

  use ModSizeDGCPM
  use ModIoUnit, only : UnitTmp_

  logical :: IsFramework

  ! Define a few time and geophysical index input variables
  ! Formerly: Common block PARAM2
  integer year,day,nstep,ikp,iwpi
  real ut,r,ap,kp,f107,tint,time

  integer :: nst, nkp, nibc, i2

  ! Define some variables set with the input file
  ! Formerly: Common block PARAM3
  real TINJ
  character*5 name
  integer INI(NS),IRES(15),ISTORM,IST,IBC(NS),IA

  ! Define solar wind input variables
  ! Formerly: Common block CSWIND
  integer ISW
  real BYSW,BZSW,MDSW,USW,DPSW

  ! Define SW-dependent nightside plasma input variables
  ! Formerly: Common block CSWBND
  integer ISWB
  real NSWB,USWB,Ninj,Einj,Kinj

  ! Define input parameters for independent variable definition
  ! Formerly: Common block PARAM4
  integer ipa,ifac
  real elb,swe,rw,hmin

  ! Define variables for continuous output stream (source-loss numbers)
  ! Formerly: Common block CDNDE
  REAL RNS,RNL,RES,REL,ESN,ELN,ESE,ELE,ECN,ECE,ALN,ALE,CEN,CEE
  REAL  CONSL(NE,NS),LMP(NT)
  INTEGER ILMP(NT),Ilold(NT)

  ! Define input variables for initial condition setup
  ! Formerly: Common block PARAMB
  real FINI(NE),CHI(NR,NT)

  ! Define input variables for source cone boundary condition
  ! Formerly: Common block PARAMA
  integer Ib
  real Ab, Eob

  ! Define convection input parameters
  ! Formerly: Common block PARAM5
  integer ilame,ilambe,ippcm,ippc
  parameter (ilame=450, ippcm=5000)
  real lamgam,lambe(ilame),tlame(ilame),ppc(ippcm),tppc(ippcm)

  !----------------------------------------------------------------
  ! Stuff added by Aaron
  !----------------------------------------------------------------

  integer, parameter        :: iCharLen_     = 100

  integer                   :: iOutputUnit_  = UnitTmp_
  integer                   :: iInputUnit_   = UnitTmp_
  integer                   :: iRestartUnit_ = UnitTmp_

  integer, parameter        :: nInputMaxLines = 10000
  integer                   :: nInputLines
  character (len=iCharLen_) :: cInputText(nInputMaxLines)

  character (len=iCharLen_) :: cInputFile = "input.dgcpm"
  character (len=10) :: cOutputDir = "PS/Output/"
  character (len= 9) :: cInputDir  = "PS/Input/"

end Module ModIoDGCPM
