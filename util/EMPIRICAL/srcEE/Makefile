#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libEMPIRICALEE.a

#include Makefile.DEPEND
#include Makefile.RULES

OBJECTS = \
	EEE_ModMain.o


LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

distclean: clean
