#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

include ../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${LIBRARYDIR} -I${COUPLERDIR}

include ../Makefile.conf

# OBJECT FILE FOR THE WRAPPER

OBJECTS = \
	PT_wrapper.o \
	amps2swmf.o  \
	amps_param_interface.o

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

${MY_LIB}: ${AMPS_LIB} ${OBJECTS}
	cp -f ${AMPS_LIB} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

PT_wrapper.o: ${LIBRARYDIR}/CON_comp_info.o

amps_interface.o:

#
#	cleaning
#

distclean: clean
