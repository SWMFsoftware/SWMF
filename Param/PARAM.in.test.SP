Begin session: 1

#COMPONENTMAP
SP  0 -1  1		CompMap

#DESCRIPTION
Param/PARAM.in.test.SP: Input file for testing SP model for SEPs, uses input files with background data		StringDescription

#STARTTIME
2011		iYear
   2		iMonth
  16		iDay
  17		iHour
  26		iMinute
   0		iSecond
   0		FracSecond

#TIMEACCURATE
T			IsTimeAccurate

#FIELDLINE
SP			NameTarget
0			nSource

#BEGIN_COMP SP ---------------------------------------------------------------

#DORUN
F			DoRun

#TRACESHOCK
F			DoTraceShock

#READMHDATA
T			DoReadMhData
MH_data			NameInputDir

#INCLUDE
MH_data/MH_data.H

#SAVEPLOT
1			nFileOut
mh1d flux idl		StringPlot

#END_COMP SP -----------------------------------------------------------------

#STOP
-1			MaxIter
30			TimeMax

End session: 1
#END ##########################################################################
