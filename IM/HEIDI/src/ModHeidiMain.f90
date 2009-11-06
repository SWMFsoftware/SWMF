Module ModHeidiMain
  !\
  ! The main variable definition module for the HEIDI program.
  ! Mike Liemohn, March 2006
  !/
  ! Include the parameter settings for MPI parallel processing
  use ModMpi
  use ModHeidiSize
  
  ! Define "constants" of the simulation
  ! Formerly: Common block CONST
  real :: q,me,mp,re,dkp,fluxfact(NS),pi
  real ::dayr(48),rkph(48),f107r(48),apr(48),rsunr(48)
  
  ! Define independent variables and grid-related factors
  ! Formerly: Common block CAR
  integer :: upa(NR)
  real    :: DL1,DR,LZ(NR),Z(NR),BE(NR,Slen),PHI(NT),DPHI,MLT(NT)
  real    :: MAS(NS),M1(NS),WE(NE),DE(NE),EKEV(NE),V(NE,NS)
  real    :: VBND(NE,NS),MU(NPA),DMU(NPA),WMU(NPA),EBND(NE)
  real    :: CONMU1,CONMU2,FFACTOR(NR,NT,NE,NPA),FACMU(NPA,NR,NT),CONF1,CONF2
  real    :: CEDR(NR,NT,NE,NPA,NS),CIDR(NR,NT,NE,NPA,NS)
  
  ! Define flux variable, and a few others
  ! Formerly: Common block CF2
  real ::F2(NR,NT,NE,NPA,NS),A,T,FGEOS(NT,NE,NPA,NS)
  
  ! Define parameters based on grid variables
  ! formerly: Common block CINIT
  real ::ENER(NR,NS),FACTOR(NS),LEC(NR,NS),ECOF(NR),WCD(NR,NT)
  
  ! Define thermal plasma variables
  ! Formerly:  Common block CXNE
  real    :: xne(NR,NT)
  integer :: itherminit,ithermfirst
  
  
  ! Define parallel computing variables
  ! Formerly: Common block PAR
  integer :: nParallelSpecies(NS)
  integer :: nSpecies, iSpecies

  real :: funt(nPa,nR,nT),funi(nPa,nR,nT)
end Module ModHeidiMain

