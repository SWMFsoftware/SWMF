Begin session: 1

#ECHO
T			DoEcho

#DESCRIPTION
Input file for testing SP model for SEPs, uses input files with background data


#TIMEACCURATE
T

#BEGIN_COMP SP ----------------------------------------------------------------

#DORUN
F

#TEST
F			DoTraceShock
F			UseDiffusion

#GRID
1.05			RMin
19			RBufferMin
21			RBufferMax
220			RMax


#READMHDATA
T
MH_data
ascii

#INCLUDE
MH_data/MH_data.H

#SAVEPLOT
1
mh1d idl

#END_COMP SP   ----------------------------------------------------------------

#STOP
-1			MaxIter
30			TimeMax

End session: 1
#END ##########################################################################
