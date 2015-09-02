SHELL =/bin/sh

DEFAULT_TARGET : dust  

include ../../../Makefile.conf

DUSTOBJ = Dust_sampling.o Dust.o  Dust_charging.o

dust : $(DUSTOBJ)
	ar -src dust.a $(DUSTOBJ)
	@echo "Target 'dust' is done." 
