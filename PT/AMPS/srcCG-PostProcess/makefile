SHELL =/bin/sh

DEFAULT_TARGET : lib  

include ../../Makefile.conf

#the following is a marker that is used by ./Config.pl to add extra compiler parameters that are custom to the application
#FLAGC

MAINOBJ=main.o
MAINLIBOBJ=main.o \
	DustScatteringEfficientcy.o \
	VIRTIS-M.o \

lib:	$(MAINLIBOBJ)  
	ar -rc mainlib.a $(MAINLIBOBJ) 

amps:	$(MAINOBJ)	
	ar -rc main.a $(MAINOBJ)
