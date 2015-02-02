SHELL =/bin/sh

DEFAULT_TARGET : electronimpact 

include ../../../Makefile.conf

ELECTRONIMPACTOBJ = ElectronImpact_Burger2010SSR.o  

electronimpact : $(ELECTRONIMPACTOBJ)
	ar -src electronimpact.a $(ELECTRONIMPACTOBJ)
	@echo "Target 'electronimpact' is done." 
