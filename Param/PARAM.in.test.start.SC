#ECHO
T

#DESCRIPTION
Input file for creating a steady state solar corona

#TIMEACCURATE
F

#STARTTIME
1997			iYear
4			iMonth
24			iDay
17			iHour
2			iMinute
37			iSecond
0.0 			FracSecond

#SAVERESTART
T			DoSaveRestart
2000			DnSaveRestart
-1.0			DtSaveRestart

#COMPONENT
IH
F

#BEGIN_COMP SC --------------------------------------------------------------



#COORDSYSTEM
HGR                     TypeCoordSystem

#USER_FLAGS
T			UseUserInnerBCs
F			UseUserSource
T			UseUserPerturbation
F			UseUserOuterBcs
F			UseUserICs
T			UseUserSpecifyRefinement
F			UseUserLogFiles
F			UseUserWritePlot
F			UseUserAMR
F			UseUserEchoInput
T			UseUserB0
F			UseUserSetPhysConst
T			UseUserUpdateStates

#USERINPUTBEGIN ------------------------

#PFSSM
T			UseUserB0
1.0			Ro_PFSSM
2.5			Rs_PFSSM
0.0			H_PFSSM
SC/CR1922_MDI.dat       File_PFSSM
12			iHead_PFSSM
-1.0			PhiShift
1.0			UnitB
-1.0			dt_UpdateB0

#EMPIRICALSW
WSA			SWModelType

#USERINPUTEND ----------------------------

#BODY
T			Usebody
1.0			rBody
3.0E+8			Body_Rho_Dim
2.0E+6			Body_T_Dim

#GAMMA
1.5000E+00		g

#GRID
1			nProcX
1			nProcY
1			nProcZ
-24.0			xMin
 24.0			xMax 
-24.0			yMin
 24.0			yMax
-24.0			zMin
 24.0			zMax

#OUTERBOUNDARY
float	                TypeBcOuter_E(East_)
float   	        TypeBcOuter_E(West_)
float	                TypeBcOuter_E(South_)
float	  		TypeBcOuter_E(North_)
float	 	        TypeBcOuter_E(Bot_)
float	   		TypeBcOuter_E(Top_)

#AMRINIT
helio_init		TypeAmrInit
9			nRefineLevelInit

#AMR
3000			DnRefine
F			DoAutomaticRefinement 

#SCHEME
2			nOrder (1 or 2)
Sokolov			TypeFlux (Roe, Rusanov, Linde, Sokolov
mc                      TypeLimiter
1.2			BetaLimiter

#TVDRESCHANGE
T                       UseTVDatResChange

#RESCHANGEBOUNDARY
T

#DIVB
T
T
F
F

#TIMESTEPPING
2                       nStage
0.80                    CflExpl

#SAVELOGFILE
T			DoSaveLogfile
RAW			StringLogfile
10			DnSaveLogfile
-1.			DtSaveLogfile

#SAVEPLOT
4			nPlotFiles
y=0 MHD tec		StringPlot
1000                   	DnSavePlot
-1.			DtSavePlot
y=0 MHD tec		StringPlot
1000			DnSavePlot
-1			DtSavePlot
z=0 MHD tec		StringPlot
1000			DnSavePlot
-1.			DtSavePlot
3d MHD tec		StringPlot
14000                   DnSavePlot
-1.			DtSavePlot

#END_COMP SC ----------------------------------------------------------------

#STOP
6000
-1.0


#RUN ########################################################################

#BEGIN_COMP SC --------------------------------------------------------------

#AMR
-1			DnRefine
F			DoAutomaticRefinement 

#END_COMP SC -----------------------------------------------------------------

#STOP
14000
-1.0

#END
