#ECHO
T

#DESCRIPTION
Uses flat files for MHD data and starts SP


#TIMEACCURATE
T			DoTimeAccurate


#BEGIN_COMP SP



#RTRANSIENT
2.0

#NSMOOTH
3

#DORUN
T

#DOREADMHDATA
T

#VERBOSE
1

#PLOT
1                     DnPlot

#END_COMP SP



#STOP
-1                    MaxIteration
5999.9                tSimulationMax

#END
