#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../../../Makefile.def

SEARCH_EXTRA = -I${LIBRARYDIR} -I${COUPLERDIR}

include ../../../Makefile.conf

# OBJECT FILE FOR THE WRAPPER

OBJECTS = \
	PT_wrapper.o \
	amps2swmf.o

#
#	Making library
#
MY_LIB    = ${LIBDIR}/libPT.a

LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

AMPS_LIB = ../srcTemp/libAMPS.a

${MY_LIB}: ${AMPS_LIB} ${OBJECTS} ${MAKEFILE_COMP_SELECT}
	cp -f ${AMPS_LIB} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

PT_wrapper.o: ${LIBRARYDIR}/CON_comp_info.o

#
#	cleaning
#

distclean: clean
