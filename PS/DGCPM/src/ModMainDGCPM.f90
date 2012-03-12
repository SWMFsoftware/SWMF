Module ModMainDGCPM
  !\
  ! The main variable defination program for DGCPM 
  ! Aron Dodger, May 2010
  !/

  use ModSizeDGCPM

  integer demo
  logical :: IsUninitialized

  ! Define "constants" of the simulation
  ! Formerly: Common block CONST
  real dayr(250),rkph(250),f107r(250),apr(250),rsunr(250), t, A, AP
  integer :: itherminit=1
  integer :: ithermfirst

  ! Define Debug Variable
  ! Controls write(*,*) output
  integer :: debug = -1

  ! Define Output Variable structure
  ! Formerly ModHeidiDGCPM
  real :: vthetacells(nthetacells),vphicells(nphicells)
  real :: vlzcells(nthetacells),vmltcells(nphicells)
  real :: vrcells(nthetacells)
  real :: potdgcpm(nthetacells,nphicells)
  real :: dendgcpm(nthetacells,nphicells)
  real :: mgridx(nthetacells,nphicells)
  real :: mgridy(nthetacells,nphicells)
  real :: mgridoc(nthetacells,nphicells)
  real :: mgridcoro(nthetacells,nphicells)
  real :: mgridvr(nthetacells,nphicells)
  real :: mgridvp(nthetacells,nphicells)
  real :: mgridsource(nthetacells,nphicells)
  real :: mgridfluxa(nthetacells,nphicells)
  real :: mgridfluxr(nthetacells,nphicells)
  real :: mgridn(nthetacells,nphicells)
  real :: mgridvol(nthetacells,nphicells)
  real :: mgridden(nthetacells,nphicells)
  real :: mgridpot(nthetacells,nphicells)
  real :: mgridb(nthetacells,nphicells)
  real :: mgridbi(nthetacells,nphicells)
  real :: mgrider(nthetacells,nphicells)
  real :: mgridep(nthetacells,nphicells)
  real :: mgridhalf(nthetacells,nphicells)
  real delr, delphi
  real :: CrossPotential=0.0
  real :: PotentialRatio=0.0
  
  ! Test Variables for WEIMER Errors.    
  real test1_den, test1_pot, test1_er, test1_ep, test1_vr, test1_vp
  real test1_fluxp, test1_fluxr, test1_sourcep, test1_sourcer
  real test2_den, test2_pot, test2_er, test2_ep, test2_vr, test2_vp
  real test2_fluxp, test2_fluxr, test2_sourcep, test2_sourcer


  ! Define electric field related variables
  character (len=100) :: EFieldModel ='VS'
  logical :: UseStaticKP = .false.

  ! Define Shue Magenetopause Variables
  real :: IMF_Bz
  real :: SW_V
  real :: SW_N
  logical :: useShue = .false.

  ! Define Output Variables
  logical :: writeStatic = .false.
  logical :: WriteLogFile = .false.
  logical :: WriteDynamic = .true.
  logical :: WriteRestart = .false.

  ! Define Testing Variables
  integer :: TestFill = 0

  ! Define Testing Variables for Filling/Emptying Tests
  real :: EmptyPeriodClosed = 3.0
  real :: EmptyPeriodOpen   = 1.0
  real :: FluxMax           = 2.0E12   

  ! Define Output Variable Type
  character(10) :: OutputType = 'OLD'
  character(10) :: MagneticType = 'DIPOLE'
 
end Module ModMainDGCPM

