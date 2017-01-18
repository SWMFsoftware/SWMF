SHELL =/bin/sh

DEFAULT_TARGET : lib  

include ../../Makefile.conf

MAINLIBOBJ = Earth_Sampling.o main_lib.o Earth.o BoundaryInjection.o BoundaryInjection_SEP.o BoundaryInjection_GCR.o GCR_Badavi2011ASR.o CutoffRigidity.o ImpulseSource.o
MAINOBJ = main.o

lib : $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'main' is done." 

amps : $(MAINOBJ) 
	ar -rc main.a $(MAINOBJ)

