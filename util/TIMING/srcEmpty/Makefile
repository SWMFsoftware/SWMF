#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../../../Makefile.def
include ../../../Makefile.conf

# OBJECT FILES FOR THE EMPTY TIMING LIBRARY

OBJECTS = \
	timing_empty.o

#
#	Making library
#
MY_LIB = ${LIBDIR}/libTIMING.a

LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS} ${MAKEFILE_COMP_SELECT}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

#
#	cleaning
#

distclean: clean

# keep this line
