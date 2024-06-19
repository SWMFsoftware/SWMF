#COMPONENTMAP
SP 0 0 1

#DESCRIPTION
Param/PARAM.in.test.start.SP: Uses flat files for MHD data and starts SP

#TIMEACCURATE
T			DoTimeAccurate


#BEGIN_COMP SP	-------------------------------------------------------------

#RTRANSIENT
2.0

#NSMOOTH
3

#DORUN
T			DoRun

#DOREADMHDATA
T			DoReadMhData

#VERBOSE
1			lVerbose

#PLOT
1			DnPlot

#END_COMP SP ---------------------------------------------------------------

#STOP
-1                    MaxIteration
5999.9                tSimulationMax
