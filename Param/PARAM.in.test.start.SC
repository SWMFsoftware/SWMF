!- PARAM.in: This is an input file for the 
!- NEW heliosphere

#ECHO
T

#TIMEACCURATE
F

#STARTTIME
1998		year
5		month
1		day
0		hour
0		minute
0		second
0.0 		FracSecond

#BEGIN_COMP SC


#PROBLEMTYPE
3                       problem_type

#USER_FLAGS
F			UseUserInnerBCs
F			UseUserSource
T			UseUserPerturbation
F			UseUserOuterBcs
F			UseUserICs
T			UseUserSpecifyRefinement
T			UseUserLogFiles
F			UseUserWritePlot
F			UseUserAMR
F			UseUserEchoInput
T			UseUserB0
F			UseUserSetPhysConst
T			UseUserUpdateStates

#USERINPUTBEGIN

#USEUSERHEATING
F			UseUserHeating

#PFSSM
T			UseUserB0
1.00E+00		Ro_PFSSM
2.50E+00		Rs_PFSSM
SC/CR1935_WSO.dat       File_PFSSM
12			iHead_PFSSM
-1.00E+00		dt_UpdateB0

#AWHEAT
1.50E+00		Bnot	[G]
1.6304543E+00		Tnot
9			DegFrm1
2.10E+00		TemRatio
F     			UseStaticICs
F                       UseZeroOmega

#USERINPUTEND

#HELIOSPHERE
1.80E+06		Tsun    [K]      1.0E6 
5.00E+13		NDsun [n/m-3]
10.00			Qsun
1.76			Theat
1.00			Rheat
4.50			SIGMAheat
F			UseFluxRope

#INERTIAL
T			UseInertial, if false then corotating frame
F			UseCorotation

#GAMMA
1.50000E+00		g

#GRID
1			proc_dims(1)
1			proc_dims(2)
1			proc_dims(3),
-24.0			x1  
 24.0			x2  
-24.0			y1 
 24.0			y2 
-24.0			z1 
 24.0			z2 

#FACEOUTERBC
6                       MaxBoundary
F                       DoFixOuterBoundary

#OUTERBOUNDARY
heliofloat		TypeBcOuter_E(East_)
heliofloat		TypeBcOuter_E(West_)
heliofloat		TypeBcOuter_E(South_)
heliofloat		TypeBcOuter_E(North_)
heliofloat		TypeBcOuter_E(Bot_)
heliofloat		TypeBcOuter_E(Top_)

#AMRINIT
helio_init
9			initial_refinement_levels

#OPTIMIZE
allopt			OptimizeMessagePass

#SCHEME
2			nORDER (1 or 2)
4			fluxfcn_type: 1=Roe 2=Rusanov 3=Linde 4=Sokolov
minmod			limiter_type mc or minmod

#TIMESTEPPING
2                       nStage
0.80                    CflExpl

#DIVB
T			UseDivbSource
F			UseDivbDiffusion
F			UseProjection
F			UseConstrainB

#SAVERESTART
T			save_restartfile
2000			dn_output(restart_)
-1.0			dt_output(restart_)


#SAVELOGFILE
T                       DoSaveLogfile
VAR                     StringLog
10                      dn_output(logfile_)
-1.0			dt_output(logfile_)
em_r ek_r et_r
4.0  10.0               rLog

#SAVEPLOT
3			nplotfiles
x=0 MHD tec             plot_string
10000                   dn
-1.0			dt
y=0 MHD tec		plot_string
10000			dn
-1.0			dt
z=0 MHD tec		plot_string
10000			dn
-1.0			dt

#END_COMP SC

#STOP
70000			nITER
-1.0			t_max

#END
