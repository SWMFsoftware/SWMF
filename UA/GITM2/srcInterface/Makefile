#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND

# Configurable objects

OBJECTS = \
	ModProcUA.o\
	UA_finalize.o\
	UA_run.o\
	UA_init_session.o\
	UA_save_restart.o\
	UA_set_param.o\
	UA_extras.o

SEARCH_EXTRA = -I${COUPLERDIR} -I${LIBRARYDIR} -I../src -I../srcSphereAB

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

#
#	Making executables and libraries
#
MY_LIB = ${LIBDIR}/libUA.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

GITM_LIB = ../src/libGITM.a 

${MY_LIB}: ${GITM_LIB} ${OBJECTS} ${MAKEFILE_COMP_SELECT}
	cp -f ${GITM_LIB} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

distclean: clean
	rm -f Makefile.DEPEND

# keep this line
