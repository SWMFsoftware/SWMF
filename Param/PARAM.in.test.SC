#COMPONENTMAP
SC  0 -1  1		CompMap runs on all PEs

#DESCRIPTION
Param/PARAM.in.test.SC: Input file for creating a steady state solar corona with Threaded-Field-Line LC		StringDescription

#TIMEACCURATE
F			IsTimeAccurate

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


#BEGIN_COMP SC ---------------------------------------------------------------

#GRIDGEOMETRY
spherical_lnr		TypeGeometry

! Use a true sphere, so box is made much larger
#GRID
1			nRootBlock1
2			nRootBlock2
1			nRootBlock3
-100.0			xMin
 100.0			xMax
-100.0			yMin
 100.0			yMax
-100.0			zMin
 100.0			zMax

#GRIDBLOCKALL
1200			MaxBlock

#LIMITRADIUS
1.05			rMin
24.0			rMax

#DOAMR
T			DoAmr
60			DnAmr
-1.0			DtAmr
T			IsStrictAmr

#AMRREGION
InnerShell		NameRegion
box_gen			StringShape
1.05			Coord1MinBox
0.0			Coord2MinBox
-70.0			Coord3MinBox
1.7			Coord1MaxBox
360.0			Coord2MaxBox
70.0			Coord3MaxBox

#AMRCRITERIALEVEL
3			nRefineCrit
Level			StringRefine
2			RefineTo
2			CoarsenFrom
Level +InnerShell		StringRefine
5			RefineTo
5			CoarsenFrom
currentsheet		StringRefine
0.5			CoarsenLimit
0.5			RefineLimit
5			MaxLevel

#GRIDLEVEL
2			nLevel
initial			StringShape

#COORDSYSTEM
HGR			TypeCoordSystem

#PLASMA
1.0			FluidMass [amu]
1.0			IonCharge [e]
1.0			ElectronTemperatureRatio

#FIELDLINETHREAD
T			UseFieldLineThreads
45			nPointThreadMax
3.0e-3			DsThreadMin

#THREADEDBC
T			UseAlignedVelocity
T			DoConvergenceCheck
second			TypeBc
1e-6			Tolerance
20			MaxIter

#HARMONICSFILE
SC/Param/CR2077_GNG.dat		NameHarmonicsFile

#CURLB0
T			UseCurlB0
2.5			rCurrentFreeB0
F			UseB0MomentumFlux

#USERSWITCH
+init +ic		StringSwitch

#POYNTINGFLUX
1.1e6			PoyntingFluxPerBSi [J/m^2/s/T]

#CORONALHEATING
turbulentcascade		TypeCoronalHeating
T			UseReflectionParam
1.5e5			LperpTimesSqrtBSi
0.0			rMinWaveReflection
F			UseReynoldsDecomposition

#HEATPARTITIONING
stochasticheating		TypeHeatPartitioning
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
T			UseHeatFluxCollisionless
1.05			CollisionlessAlpha

#TEST
krylov radiative_cooling		StringTest

#SEMIIMPLICIT
T			UseSemiImplicit
parcond			TypeSemiImplicit

#SEMIKRYLOV
GMRES			TypeKrylov
1.0e-5			ErrorMaxKrylov
10			MaxMatvecKrylov

#RADIATIVECOOLING
T			UseRadCooling

#LOOKUPTABLE
radcool			NameTable
load			NameCommand
SC/Param/RadCoolCorona.dat		NameFile
ascii			TypeFile

#LOOKUPTABLE
TR			NameTable
load			NameCommand
SC/Param/TR.dat		NameFile
ascii			TypeFile


#TIMESTEPPING
2			nStage
0.8			CflExpl

#MINIMUMPRESSURE
1.0E-9			pMinDim
1.0e-9			PeMinDim

#MINIMUMTEMPERATURE
5.0e4			TminDim
5.0e4			TeMinDim


#SCHEME
2			nOrder (1 or 2)
Sokolov			TypeFlux (Roe, Rusanov, Linde, Sokolov
mc3			TypeLimiter
1.2			LimiterBeta

#LIMITER
T			UseLogRhoLimiter
T			UseLogPLimiter
F			UseRhoRatioLimiter

#NONCONSERVATIVE
T			UseNonConservative


#MESSAGEPASS
all			TypeMessagePass

#TVDRESCHANGE
T			UseTvdReschange

#SAVELOGFILE
T			DoSaveLogfile
RAW			StringLog
1			DnSaveLogfile
-1.			DtSaveLogfile

#COARSEAXIS
T			UseCoarseAxis
2			nCoarseLayer

#OUTERBOUNDARY
fieldlinethreads		TypeBc1   user for spherical and user_outerbcs
float			TypeBc2
float			TypeBc3
float			TypeBc4
float			TypeBc5
float			TypeBc6

#SAVEPLOT
3			nPlotFile
x=0 VAR idl		StringPlot
10			DnSavePlot
-1.			DtSavePlot
-1.0			DxSavePlot
rho ux uy uz ti te bx by bz I01 I02 p		NameVars
g rbody			NamePars
y=0 VAR idl		StringPlot
10			DnSavePlot
-1.			DtSavePlot
-1.0			DxSavePlot
rho ux uy uz ti te bx by bz I01 I02 p		NameVars
g rbody			NamePars
z=0 VAR idl		StringPlot
10			DnSavePlot
-1.			DtSavePlot
-1.0			DxSavePlot
rho ux uy uz ti te bx by bz I01 I02 p		NameVars
g rbody			NamePars

#SAVEINITIAL
T			DoSaveInitial

#END_COMP SC -----------------------------------------------------------------

#STOP
70			MaxIter
-1.0			TimeMax
