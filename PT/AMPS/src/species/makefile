SHELL =/bin/sh

DEFAULT_TARGET : species  

include ../../Makefile.conf

SPECIESOBJ = Na.o He.o 

species : $(SPECIESOBJ)
	ar -src species.a $(SPECIESOBJ)
	@echo "Target 'species' is done." 
