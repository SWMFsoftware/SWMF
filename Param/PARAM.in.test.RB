#ECHO
T

#VERBOSE
100

#DESCRIPTION
Test RB/RBE only

#SAVERESTART
F			DoSaveRestart

#STOP
-1			MaxIteration
60.0			tSimulationMax


! This is needed so that RBE reads its input file
#BEGIN_COMP RB

#END_COMP RB
