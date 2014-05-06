#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

include ../../../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${LIBRARYDIR} -I${COUPLERDIR}

include ../../../Makefile.conf

SC_wrapper.f90:
	cp -f ../../../IH/Empty/src/IH_wrapper.f90 SC_wrapper.f90
	perl -i -pe \
	's/IH/SC/g;s/BATSRUS/SC_BATSRUS/;s/Inner/Solar/;s/Heliosphere/Corona/'\
		SC_wrapper.f90

# OBJECT FILE FOR THE WRAPPER

OBJ_WRAPPER = \
	SC_wrapper.o

#
#	Making library
#
MY_LIB = ${LIBDIR}/libSC.a


LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJ_WRAPPER}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJ_WRAPPER}

SC_wrapper.o: ${LIBRARYDIR}/CON_comp_info.o

#
#	cleaning
#

distclean: clean
	   rm -f SC_wrapper.f90

# keep this line
