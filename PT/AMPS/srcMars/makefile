SHELL =/bin/sh

DEFAULT_TARGET : lib  

include ../../Makefile.conf

MAINLIBOBJ = newMars.o MarsBackgroundAtmosphereFox.o BackgroundAtmosphere.o 
MAINOBJ = mars.o

lib : $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'main' is done." 

amps : $(MAINOBJ) 
	ar -rc main.a $(MAINOBJ)

