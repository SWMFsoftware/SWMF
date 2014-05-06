
SHELL =/bin/sh

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND

SEARCH_EXTRA = -I${COUPLERDIR} -I${LIBRARYDIR} -I../src

OBJECTS = SP_wrapper.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

#
#	Making executables and libraries
#
MY_LIB = ${LIBDIR}/libSP.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

FLAMPA_LIB = ../src/libFLAMPA.a 

${MY_LIB}: ${FLAMPA_LIB} ${OBJECTS}
	cp -f ${FLAMPA_LIB} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

distclean: clean
	rm -f Makefile.DEPEND

# keep this line
