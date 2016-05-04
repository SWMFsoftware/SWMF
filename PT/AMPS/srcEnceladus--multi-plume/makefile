SHELL =/bin/sh

DEFAULT_TARGET : lib  

include ../../Makefile.conf

MAINLIBOBJ = EnceladusMultiPlume_SourceLocation.o main_lib.o \
	EnceladusMultiPlume_TigerStripes3D.o EnceladusMultiPlume.o
MAINOBJ = main.o

lib : $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'main' is done." 

amps : $(MAINOBJ) 
	ar -rc main.a $(MAINOBJ)

