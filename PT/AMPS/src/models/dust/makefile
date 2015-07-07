SHELL =/bin/sh

DEFAULT_TARGET : dust  

include ../../../Makefile.conf

DUSTOBJ = Dust_sampling.o  

dust : $(DUSTOBJ)
	ar -src dust.a $(DUSTOBJ)
	@echo "Target 'dust' is done." 
