#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.RULES

MY_LIB = ${LIBDIR}/libEMPIRICALUA.a

OBJECTS =\
	ModIri90.o\
        ModMsis90.o

LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}:  ${MODULES}  ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS} ${MODULES}

distclean: clean

