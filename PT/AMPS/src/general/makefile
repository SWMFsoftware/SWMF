SHELL =/bin/sh

DEFAULT_TARGET : general  

include ../../Makefile.def
include ../../Makefile.conf

SEARCH=

GENERALOBJ = specfunc.o rnd.o global.o 

general : $(GENERALOBJ)
	ar -src general.a $(GENERALOBJ)
	@echo "Target 'general' is done." 
