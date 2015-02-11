SHELL =/bin/sh

DEFAULT_TARGET : sputtering

include ../../../Makefile.conf

SPUTTERINGOBJ = Sputtering.o Sputtering_Ice.o

sputtering : $(SPUTTERINGOBJ)
	ar -src sputtering.a $(SPUTTERINGOBJ)
	@echo "Target 'sputtering' is done." 
