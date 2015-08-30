SHELL =/bin/sh

DEFAULT_TARGET : interface

include ../../Makefile.conf
include ../../Makefile.def
include makefile.cell_centered_linear_interpolation


INTERFACES=

interface :
	$(foreach INTERFACE, ${INTERFACES}, ${MAKE} ${INTERFACE})
	${AR} interface.a ./*.o
	@echo "Target 'interface' is done." 
