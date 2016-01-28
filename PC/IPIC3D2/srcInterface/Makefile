#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

include ../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${LIBRARYDIR} -I${COUPLERDIR} 

SEARCH_C = -I../src/include -I./

FLAGC_EXTRA = -DBATSRUS

include ../Makefile.conf

# OBJECT FILE FOR THE WRAPPER

OBJECTS = \
	PC_wrapper.o \
	ipic3d_interface.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH_C} ${OBJECTS} ${FLAGC_EXTRA}


LIBIPIC3D = ../libIPIC3D/libIPIC3D.a

#	Making library
#
MY_LIB    = ${LIBDIR}/libPC.a

LIB:    DEPEND ${LIBIPIC3D} 
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS} ${LIBIPIC3D}
	cp -f ${LIBIPIC3D} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

PC_wrapper.o: ${LIBRARYDIR}/CON_comp_info.o

#
#	cleaning
#

clean: cleanfiles

distclean: clean
	rm -f Makefile.DEPEND Makefile.RULES
