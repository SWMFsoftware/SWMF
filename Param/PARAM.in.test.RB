#ECHO
T			DoEcho


! Since RBE runs serial only, empty processors would stop SWMF in strict mode
#STRICT
F			UseStrict

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
