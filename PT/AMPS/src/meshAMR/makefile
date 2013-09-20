SHELL =/bin/sh

DEFAULT_TARGET : mesh  

include ../../Makefile.conf

MESHOBJ = meshAMRgeneric.o  meshAMRcutcell.o
mesh : $(MESHOBJ)
	ar -src mesh.a $(MESHOBJ)
	@echo "Target 'mesh' is done." 
