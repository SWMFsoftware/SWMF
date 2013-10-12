!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModHeidiIO
  !\
  ! Input/output variable definition module for the HEIDI program.
  !/
  
  use ModHeidiSize, ONLY: nR, nT, nE, nS
  implicit none

  ! Define a few time and geophysical index input variables
  integer:: year,day,nstep,ikp,iwpi,month
  real   ::ut,r,ap,kp,f107,tint,time
  integer, dimension(1:7)::TimeArray
  
  ! Define some variables set with the input file
  real :: TINJ
  character(len=20) :: NameRun="test1"
  character(len=*),parameter :: NameOutputDir      = "IM/plots/"
  character(len=*),parameter :: NameRestartOutDir  = "IM/restartOUT/"
  character(len=*),parameter :: NameRestartInDir   = "IM/restartIN/"
  character(len=*),parameter :: NameInputDirectory = "IM/input/"
  integer :: INI(NS),IRES(16),ISTORM,IST,IBC(NS),IA
  
  ! Define solar wind input variables
  integer :: ISW
  real    :: BYSW,BZSW,MDSW,USW,DPSW
  
  ! Define SW-dependent nightside plasma input variables
  integer :: ISWB
  real    :: NSWB,USWB,Ninj,Einj,Kinj
  
  ! Define input parameters for independent variable definition
  integer :: ipa,ifac
  real    :: elb,swe,rw,hmin
  
  ! Define variables for continuous output stream (source-loss numbers)
  REAL    :: RNS,RNL,RES,REL,ESN,ELN,ESE,ELE,ECN,ECE,ALN,ALE,CEN,CEE
  REAL    :: CONSL(NE,NS),LMP(NT)
  INTEGER :: ILMP(NT),Ilold(NT)
  
  ! Define input variables for initial condition setup
  real :: FINI(NE),CHI(NR,NT)
  
  ! Define input variables for source cone boundary condition
  integer :: Ib
  real    :: Ab, Eob
  
  ! Define convection input parameters
  integer :: ilambe,ippc
  integer ,parameter :: ilame=450 
  integer, parameter :: ippcm=5000
  real :: lamgam,lambe(ilame),tlame(ilame),ppc(ippcm),tppc(ippcm)

  !Define new variables needed for the framework version  
  logical :: IsFramework
  integer :: iUnitStdOut =  6 ! Change to 9 for ABSOFT compiler
  character (len=7) :: StringPrefix = ''
  integer :: iUnitSw1,iUnitSw2,iUnitSopa,iUnitMpa,iUnitPot
  integer ::iUnitSal
  integer           :: iUnitSal1,iUnitSal2,iUnitSal3,iUnitSal4
contains
  ! ===========================================================================
  subroutine write_prefix
    
    write(iUnitStdOut,'(a)',ADVANCE='NO')trim(StringPrefix)
    
  end subroutine write_prefix
  ! ===========================================================================
end Module ModHeidiIO
