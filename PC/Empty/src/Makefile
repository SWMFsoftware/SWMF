#This code is a copyright protected software (c) 2002- University of Michigan

SHELL =/bin/sh

include ../../../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${LIBRARYDIR} -I${COUPLERDIR}

include ../../../Makefile.conf

PC_wrapper.f90:
	cp -f ../../../IH/Empty/src/IH_wrapper.f90 PC_wrapper.f90
	perl -i -pe \
	's/IH/PC/g;s/BATSRUS/PC_BATSRUS/;s/Inner/Solar/;s/Heliosphere/Corona/'\
		PC_wrapper.f90	

# OBJECT FILE FOR THE WRAPPER

OBJ_WRAPPER = \
	PC_wrapper.o


#
#	Making library
#
MY_LIB = ${LIBDIR}/libPC.a


LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJ_WRAPPER}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJ_WRAPPER}

PC_wrapper.o: ${LIBRARYDIR}/CON_comp_info.o

#
#	cleaning
#

distclean: clean
