#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../../../Makefile.def

SEARCH_EXTRA = -I${LIBRARYDIR} -I${COUPLERDIR}

include ../../../Makefile.conf

OH_wrapper.f90:
	cp -f ../../../IH/Empty/src/IH_wrapper.f90 OH_wrapper.f90
	 perl -i -pe \
	's/IH/OH/g;s/_sc/_ih/;s/BATSRUS/OH_BATSRUS/;s/Inner/Outer/;'\
		OH_wrapper.f90

# OBJECT FILE FOR THE WRAPPER

OBJ_WRAPPER = \
	OH_wrapper.o

#
#	Making library
#
MY_LIB = ${LIBDIR}/libOH.a


LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJ_WRAPPER} ${MAKEFILE_COMP_SELECT}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJ_WRAPPER}

OH_wrapper.o: ${LIBRARYDIR}/CON_comp_info.o

#
#	cleaning
#

distclean: clean
	   rm -f OH_wrapper.f90

# keep this line
