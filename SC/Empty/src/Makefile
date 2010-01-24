#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../../../Makefile.def

SEARCH_EXTRA = -I${LIBRARYDIR} -I${COUPLERDIR}

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

${MY_LIB}: ${OBJ_WRAPPER} ${MAKEFILE_COMP_SELECT}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJ_WRAPPER}

SC_wrapper.o: ${LIBRARYDIR}/CON_comp_info.o

#
#	cleaning
#

distclean: clean
	   rm -f SC_wrapper.f90

# keep this line
