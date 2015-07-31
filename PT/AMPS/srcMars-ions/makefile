SHELL =/bin/sh

DEFAULT_TARGET : lib  

include ../../Makefile.conf

MAINLIBOBJ = main_lib.o mars-ions_init.o mars-ions_boundary.o mars-ions_exosphere.o mars-ions_source.o
MAINOBJ = main.o

lib : $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'main' is done." 

amps : $(MAINOBJ) 
	ar -rc main.a $(MAINOBJ)

