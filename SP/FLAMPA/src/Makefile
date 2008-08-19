#^CFG COPYRIGHT UM

SHELL =/bin/sh
include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND


OBJECTS = ModLogAdvection.o \
	ModLinAdvection.o \
	SP_ModTurbulence.o \
	SP_main.o \
	read_ihdata_for_sp.o \
	sp_diffusive_acceleration.o \
	write_logfile_SP.o \
	write_plotfile_SP.o 

SEPEXE = ${BINDIR}/SEP.exe

SEP: ${SEPEXE}

SEARCH_EXTRA = -I${LIBRARYDIR} -I${COUPLERDIR}

${SEPEXE}: ${OBJECTS} SP_stand_alone.o
	${LINK.f90} ${SEARCH} ${SEARCH_EXTRA} -o ${SEPEXE} \
	${LIBDIR}/libSHARE.a ${OBJECTS} SP_stand_alone.o



MY_LIB = libFLAMPA.a

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} \
		${MODULES} ${OBJECTS} SP_stand_alone.o

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
