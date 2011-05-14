#^CFG COPYRIGHT UM
SHELL=/bin/sh

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.RULES

help:
	@echo make targets:
	@echo '   libNOMPI.a'                             
	@echo '   distclean'

MY_LIB = ${LIBDIR}/libNOMPI.a

${SHAREDIR}/ModMpiOrig.o:
	cd ${SHAREDIR}; make ModMpiOrig.o

LIB:	${SHAREDIR}/ModMpiOrig.o
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: NOMPI.o
	${AR} ${MY_LIB} NOMPI.o

distclean: clean

# keep this line

