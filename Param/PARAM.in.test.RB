#ECHO
T			DoEcho

#TEST
RB_put_from_gm

#DESCRIPTION
Run GM and RB/RBE

#PLANET
EARTH

#SAVERESTART
F

#TIMEACCURATE
T			DoTimeAccurate

#STARTTIME
2000			iYear
3			iMonth
21			iDay
10			iHour
45			iMinute
0			iSecond
0.0			FracSecond

#IDEALAXES
! align rotation and magnetic axes with z axis

#ROTATION
F			UseRotation

#COUPLE1
GM			NameComp1
RB			NameComp2
-1                      DnCouple
3.0                     Dtcouple

! Allow GM to step through coupling time
#COUPLETIME
GM			NameComp
F			DoCoupleOnTime

! This is needed so that RBE reads its input file
#BEGIN_COMP RB #############################################################

#SAVEPLOT
3.0		DtSavePlot
F		UseSeparatePlotFiles
2002f296	OutName

#RESTART
F		IsRestart

#TIMESTEP
3.0		Dt

#TIMESIMULATION
0.0		tStart

#SPECIES
e		NameSpecies

#STARTUPTIME
0.0		trans

#BMODEL
MHD		NameModel
F		UseFixedB

#IEMODEL
1		iConvect

#INPUTDATA
2002_296

#PLASMASPHERE
F		UsePlasmaSphere

#END_COMP RB ###############################################################

#BEGIN_COMP GM #############################################################

#MAGNETOSPHERE
T                       UseBody
2.5                     rBody
3.5                     rCurrents
10.0			BodyRhoDim
25000.0			BodyTDim

#GRID
1			nRootBlockX
1			nRootBlockY
1			nRootBlockZ
-224.			xMin
  32.			xMax
-128.			yMin
 128.			yMax
-128.			zMin
 128.			zMax

#GRIDRESOLUTION
1.0                     Resolution
initial                 NameArea

#INCLUDE
GM/Param/EARTH/Grid

#OUTERBOUNDARY
outflow                 TypeBcEast
vary                    TypeBcWest
float                   TypeBcSouth
float                   TypeBcNorth
float                   TypeBcBottom
float                   TypeBcTop

#INNERBOUNDARY
ionosphere		TypeBcInner

#SAVELOGFILE
T
RAW			StringLogfile
1			DnSaveLogfile
-1.			DtSaveLogfile

Avoid exactly zero By so that ray tracing is more reproducible
#UPSTREAM_INPUT_FILE
T                                                          UseUpstreamInputFile
GM/Param/TESTSUITE/Inputfiles/IMF_NSturning_1nT_tilt.dat   NameUpstreamFile

#SCHEME
2			nOrder
Rusanov			TypeFlux
minmod			TypeLimiter

#TVDRESCHANGE
F

#TIMESTEPPING
2			nStage
0.8			CflExpl

#SAVEPLOT
3                       nPlotFile
x=0 FUL idl             StringPlot
-1                      DnSavePlot
10.			DtSavePlot
-1.                     DxSavePlot (resolution, 0. maximum, -1. unstructured)
y=0 FUL idl             StringPlot
-1                      DnSavePlot
10.			DtSavePlot
-1.                     DxSavePlot (resolution, 0. maximum, -1. unstructured)
z=0 FUL idl             StringPlot
-1                      DnSavePlot
10.			DtSavePlot
-1.                     DxSavePlot (resolution, 0. maximum, -1. unstructured)

#END_COMP GM ###############################################################

#STOP
-1			MaxIteration
9.0			tSimulationMax
