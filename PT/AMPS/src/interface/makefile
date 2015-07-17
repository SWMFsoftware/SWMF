SHELL =/bin/sh

DEFAULT_TARGET : interface

include ../../Makefile.conf
include makefile.interpolation_amr


INTERFACES=
INTERFACEOBJ=interface.o

interface : ${INTERFACEOBJ}
	$(foreach INTERFACE, ${INTERFACES},$(MAKE) $(INTERFACE);)
	${AR} interface.a *.o
	@echo "Target 'interface' is done." 