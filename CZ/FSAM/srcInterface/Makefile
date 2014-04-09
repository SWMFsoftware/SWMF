#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND

# Configurable objects

OBJECTS = \
	CZ_wrapper.o

SEARCH_EXTRA = -I${COUPLERDIR} -I${LIBRARYDIR} -I../src 

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

#
#	Making executables and libraries
#
MY_LIB = ${LIBDIR}/libCZ.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

FSAM_LIB = ../src/libFSAM.a 

${MY_LIB}: ${FSAM_LIB} ${OBJECTS} ${MAKEFILE_COMP_SELECT}
	cp -f ${FSAM_LIB} ${MY_LIB}
	cp ../input/libfishpak.a .
	ar -x libfishpak.a
	${AR} ${MY_LIB} ${OBJECTS} *.o
	rm -f libfishpak.a *.o

distclean: clean
	rm -f Makefile.DEPEND

