SHELL =/bin/sh

DEFAULT_TARGET : lib  

include ../../Makefile.conf

MAINLIBOBJ = MOP.o Enceladus.o SaturnianSystem.o main_lib.o 
MAINOBJ = main.o

lib : $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'MOP' is done." 

amps : $(MAINOBJ) 
	ar -rc main.a $(MAINOBJ)

