
SHELL =/bin/sh

include ../../../Makefile.def
include ../../../Makefile.conf
-include Makefile.DEPEND

OBJECTS = \
	ModGridDescriptor.o\
	SC_wrapper.o

SEARCHDIR = -I${SHAREDIR} -I${COUPLERDIR} -I${LIBRARYDIR} -I${MAGNETOGRAMDIR} \
	-I${DEMTDIR} -I${EMPIRICALEEDIR} -I../src -I../srcBATL

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${OBJECTS}

#
#	Making executables and libraries
#
MY_LIB = ${LIBDIR}/libSC.a

LIB:	DEPEND
	$(MAKE) ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

BATSRUS_LIB = ../src/libBATSRUS.a 

${MY_LIB}: ${BATSRUS_LIB} ${OBJECTS}
	cp -f ${BATSRUS_LIB} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

clean: cleanfiles

distclean: clean
	rm -f SC_wrapper.f90  ModGridDescriptor.f90 ModBuffer.f90 
