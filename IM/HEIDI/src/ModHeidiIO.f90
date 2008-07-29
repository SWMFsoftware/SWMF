Module ModHeidiIO
  !\
  ! Input/output variable definition module for the HEIDI program.
  ! Mike Liemohn, March 2006
  !/

	use ModHeidiSize

! Define a few time and geophysical index input variables
! Formerly: Common block PARAM2
	integer year,day,nstep,ikp,iwpi
	real ut,r,ap,kp,f107,tint,time	

! Define some variables set with the input file
! Formerly: Common block PARAM3
	real TINJ
	character*5 name
        character*10 :: cOutputDir="IM/Output/"
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

end Module ModHeidiIO
