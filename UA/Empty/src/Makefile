#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../../../Makefile.def

SEARCH_EXTRA = -I${LIBRARYDIR}

include ../../../Makefile.conf

# OBJECT FILE FOR THE WRAPPER

OBJ_WRAPPER = \
	UA_wrapper.o

#
#	Making library
#
MY_LIB  = ${LIBDIR}/libUA.a
MY_LIB2 = ${LIBDIR}/libIndices.a


LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} and ${MY_LIB2} have been brought up to date.
	@echo

${MY_LIB}: ${OBJ_WRAPPER} ${MAKEFILE_COMP_SELECT}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJ_WRAPPER}
	rm -f ${MY_LIB2}
	${AR} ${MY_LIB2} ${OBJ_WRAPPER}

UA_wrapper.o: ${LIBRARYDIR}/CON_comp_info.o

#
#	cleaning
#

distclean: clean

# keep this line
