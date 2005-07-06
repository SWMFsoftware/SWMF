#^CFG COPYRIGHT UM

SHELL =/bin/sh
include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND


OBJECTS = SP_main.o \
	read_ihdata_for_sp.o \
	sp_diffusive_acceleration.o \
	write_logfile_SP.o \
	write_plotfile_SP.o 


SEARCH_EXTRA = -I${LIBRARYDIR} -I${COUPLERDIR}

MY_LIB = libFLAMPA.a

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} \
		${MODULES} ${OBJECTS}

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES} ${OBJECTS} \
	${MAKEFILE_COMP_SELECT}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

distclean: clean
	rm -f Makefile.DEPEND
