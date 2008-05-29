#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../../../Makefile.def

SEARCH_EXTRA = -I${LIBRARYDIR} -I${COUPLERDIR}

include ../../../Makefile.conf

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

# keep this line
