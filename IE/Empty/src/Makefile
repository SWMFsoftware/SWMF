#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

include ../../../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${LIBRARYDIR} -I${COUPLERDIR}

include ../../../Makefile.conf

# OBJECT FILE FOR THE WRAPPER

OBJ_WRAPPER = \
	IE_wrapper.o

#
#	Making library
#
MY_LIB = ${LIBDIR}/libIE.a


LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJ_WRAPPER}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJ_WRAPPER}

IE_wrapper.o: ${LIBRARYDIR}/CON_comp_info.o

#
#	cleaning
#

clean: cleanfiles

distclean: clean

# keep this line
