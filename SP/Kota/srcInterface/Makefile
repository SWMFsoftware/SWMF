
SHELL =/bin/sh

include ../../../Makefile.def
include ../../../Makefile.conf

-include Makefile.DEPEND

# Extra path for include files
INCL_EXTRA = -I../src

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

KOTA_LIB = ../src/libKota.a 

${MY_LIB}: ${KOTA_LIB} ${OBJECTS}
	cp -f ${KOTA_LIB} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

clean: cleanfiles

distclean: clean
	rm -f Makefile.DEPEND

# keep this line
