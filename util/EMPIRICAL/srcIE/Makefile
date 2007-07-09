#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libEMPIRICALIE.a

include Makefile.DEPEND
#include Makefile.RULES

MODULES = \
	ModEIEConductance.o\
	ModAMIE_Interface.o\
	ModErrors.o\
	ModEIE_Interface.o\
	ModMHD_Interface.o\
	ModExtras.o\
	ModEIEFiles.o

OBJECTS = \
	AMIE_Library.o\
	EIE_End.o\
	EIE_Initialize.o\
	EIE_Library.o\
	IO_Library.o\
	EIE_set_inputs.o\
	readAMIEoutput.o\
	readMHDoutput.o\
	wei01.o\
	wei96.o\
	weicom.o\
	mh86.o\
	hmr89.o\
	iz94.o\
	ihp.o\
	MHD_Library.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl -I${SHAREDIR} ${MODULES}

LIB:    DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES} ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES} ${OBJECTS}

distclean: clean
	rm -f Makefile.DEPEND
