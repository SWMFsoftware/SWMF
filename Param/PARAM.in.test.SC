#ECHO
T			DoEcho

#DESCRIPTION
Input file for creating a steady state solar corona with Threaded-Field-Line LC

#TIMEACCURATE
F                       DoTimeAccurate

#STARTTIME
2011			iYear
2			iMonth
16			iDay
17			iHour
26			iMinute
0			iSecond
0.0			FracSecond

#SAVERESTART
T			DoSaveRestart
50000			DnSaveRestart
-1			DtSaveRestart


#BEGIN_COMP SC --------------------------------------------------------------

#COORDSYSTEM
HGR			TypeCoordSystem

#PLASMA
1.0			IonMass [amu]
1.0			AverageIonCharge [e]
1.0			ElectronTemperatureRatio

#FIELDLINETHREAD
T                       UseFieldLineThreads
45                      nPointInThreadMax
3.0e-3                  DsThreadMin

#THREADEDBC
T	  UseAlignedVelocity
T         DoConvergenceCheck
second    TypeBc

#BODY
T			UseBody
1.0			rBody
1.5e8			BodyNDim
1.5e6			BodyTDim

#BODY
F			UseBody

#MAGNETOGRAM
T			UseMagnetogram
1.0			rMagnetogram
2.5			rSourceSurface
0.0			HeightInnerBc
SC/Param/CR2077_GNG.dat		NameMagnetogramFile
12			nHeaderLine
0.0			PhiShift
2.0			UnitB

#USECURLB0
T			UseCurlB0
2.5			rCurrentFreeB0

#USER_FLAGS
F			UseUserInnerBcs
F			UseUserSource
F			UseUserPerturbation
F			UseUserOuterBcs
T			UseUserICs
F			UseUserSpecifyRefinement
T			UseUserLogFiles
T			UseUserWritePlot
F			UseUserAMR
F			UseUserEchoInput
F			UseUserB0
T			UseUserInitSession
F			UseUserUpdateStates

#POYNTINGFLUX
1.1e6			PoyntingFluxPerBSi [J/m^2/s/T]

#CORONALHEATING
turbulentcascade        TypeCoronalHeating
T                       UseWaveReflection
1.5e5                   LperpTimesSqrtBSi

#HEATPARTITIONING
stochasticheating	TypeHeatPartitioning
0.34			StochasticExponent
0.18			StochasticAmplitude


#HEATCONDUCTION
T			UseHeatConduction
spitzer			TypeHeatConduction

#HEATFLUXREGION
T			UseHeatFluxRegion
5.0			rCollisional
-8.0			rCollisionless

#HEATFLUXCOLLISIONLESS
T            UseHeatFluxCollisionless
1.05            CollisionlessAlpha  

#TEST
krylov radiative_cooling

#SEMIIMPLICIT
T			UseSemiImplicit
parcond			TypeSemiImplicit

#SEMIKRYLOV
GMRES			TypeKrylov
1.0e-5			ErrorMaxKrylov
10			MaxMatvecKrylov

#PRECONDITIONER
left            TypePrecondSide (left, symmetric, right)
bilu            TypePrecond (MBILU)

#RADIATIVECOOLING
T			UseRadCooling

#LOOKUPTABLE
radcool			NameTable
load			NameCommand
SC/Param/RadCoolCorona.dat
ascii			TypeFile

#LOOKUPTABLE
TR			NameTable
load			NameCommand
SC/Param/TR.dat
ascii			TypeFile

#LOOKUPTABLE
AW_TR			NameTable
load			NameCommand
SC/Param/AW_TR.dat
ascii			TypeFile

---------------------------------------------------
#GRIDGEOMETRY
spherical_lnr

! Use a true sphere, so box is made much larger
#GRID
1			nRootBlockX
2			nRootBlockY
1			nRootBlockZ
-100.0			xMin
 100.0			xMax
-100.0			yMin
 100.0			yMax
-100.0			zMin
 100.0			zMax

#LIMITRADIUS
1.05			rMin
24.0			rMax

----------------------------------------------------
#AMRREGION
InnerShell		NameRegion
box_gen			TypeRegion
1.05			Radius1
0.0
-70.0
1.7			Radius2
360.0
70.0

#AMRCRITERIALEVEL
3			nCriteria
Level			TypeCriteria
2			RefineTo
2			CoarsenTo
Level +InnerShell	TypeCriteria
5     			RefineTo
5			CoarsenTo
currentsheet		TypeCriteria
0.5			CoarsenLimit
0.5			RefineLimit
5			MaxLevel

#GRIDLEVEL
2			nLevelArea
initial			NameArea

#DOAMR
T			DoAmr
60			DnAmr
-1.0			DtAmr
T			IsStrictAmr

#TIMESTEPPING
2			nStage
0.8			CflExpl

#MINIMUMPRESSURE
1.0E-9			pMinDim
1.0e-9			peMinDim

#MINIMUMTEMPERATURE
5.0e4
5.0e4


#SCHEME
2			nOrder (1 or 2)
Sokolov			TypeFlux (Roe, Rusanov, Linde, Sokolov
mc3			TypeLimiter
1.2			LimiterBeta

#LIMITER
T
T
F

#NONCONSERVATIVE
T			UseNonConservative


#MESSAGEPASS
all			TypeMessagePass

#TVDRESCHANGE
T			UseTVDResChange

#SAVELOGFILE
T			DoSaveLogfile
RAW			StringLogfile
1			DnSaveLogfile
-1.			DtSaveLogfile

#COARSEAXIS
F			UsePoleDiffusion
T                       UseCoarseAxis
2			nCoarseLayer

#OUTERBOUNDARY
fieldlinethreads	TypeBcEast   user for spherical and user_outerbcs
float			TypeBcWest
float			TypeBcSouth
float			TypeBcNorth
float			TypeBcBottom
float			TypeBcTop

#SAVEPLOT
3			nPlotFiles
x=0 VAR idl		StringPlot
10			DnSavePlot
-1.			DtSavePlot
-1.0			DxSavePlot
rho ux uy uz ti te bx by bz I01 I02 p
g rbody
y=0 VAR idl		StringPlot
10			DnSavePlot
-1.			DtSavePlot
-1.0			DxSavePlot
rho ux uy uz ti te bx by bz I01 I02 p
g rbody
z=0 VAR idl		StringPlot
10			DnSavePlot
-1.			DtSavePlot
-1.0			DxSavePlot
rho ux uy uz ti te bx by bz I01 I02 p
g rbody
3d VAR tec
100000			DnSavePlot
-1.			DtSavePlot
rho ux uy uz ti te bx by bz I01 I02 p
g rbody

#SAVEINITIAL
T			DoSaveInitial

#END_COMP SC ----------------------------------------------------------------

#COMPONENT
EE			NameComp
F			UseComp

#STOP
70			MaxIteration
-1.0			tSimulationMax
