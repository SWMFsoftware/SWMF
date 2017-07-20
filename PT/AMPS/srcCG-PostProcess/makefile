SHELL =/bin/sh

DEFAULT_TARGET : lib  


include ../../Makefile.conf
#FLAGC
MAINOBJ=main.o
MAINLIBOBJ=main.o \
	DustScatteringEfficientcy.o \
	VIRTIS-M.o \

lib:	$(MAINLIBOBJ)  
	ar -rc mainlib.a $(MAINLIBOBJ) 

amps:	$(MAINOBJ)	
	ar -rc main.a $(MAINOBJ)
