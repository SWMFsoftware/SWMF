
SHELL =/bin/sh

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND

SEARCHDIR = -I${SHAREDIR} -I${COUPLERDIR} -I${LIBRARYDIR} -I../src

OBJECTS = SP_wrapper.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${OBJECTS}

#
#	Making executables and libraries
#
MY_LIB = ${LIBDIR}/libSP.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

MFLAMPA_LIB = ../src/libMFLAMPA.a 

${MY_LIB}: ${MFLAMPA_LIB} ${OBJECTS}
	cp -f ${MFLAMPA_LIB} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

clean: cleanfiles

distclean: clean
	rm -f Makefile.DEPEND

# keep this line
