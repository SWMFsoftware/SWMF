SHELL =/bin/sh

DEFAULT_TARGET : lib  

include ../../Makefile.conf

MAINLIBOBJ = RosinaMeasurements.o Gravity.o Comet.o CometData.o main_lib.o OrbitCalculation.o ColumnIntegral.o SurfaceDataCG.o SampleSizeDistribution.o RosinaMeasurements_Liouville.o
MAINOBJ = main.o

lib : $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'main' is done." 

amps : $(MAINOBJ) 
	ar -rc main.a $(MAINOBJ)

