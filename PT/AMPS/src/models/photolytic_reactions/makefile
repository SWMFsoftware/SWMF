SHELL =/bin/sh

DEFAULT_TARGET : photolyticreactions 

include ../../../Makefile.conf

PHOTOLITICREACTIONSOBJ = PhotolyticReactions_Huebner1992ASS.o  

photolyticreactions : $(PHOTOLITICREACTIONSOBJ)
	ar -src photolyticreactions.a $(PHOTOLITICREACTIONSOBJ)
	@echo "Target 'photolyticreactions' is done." 
