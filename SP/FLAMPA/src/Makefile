#^CFG COPYRIGHT UM

SHELL =/bin/sh
include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND

OBJECTS = sep19.o main.o read_ihdata.o

SEARCH_EXTRA = -I${LIBRARYDIR} -I${COUPLERDIR}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${SEARCH_EXTRA} ${OBJECTS}

SEPEXE = ${BINDIR}/SEP.exe

SEP: ${SEPEXE}

${SEPEXE}:	${OBJECTS} sp_stand_alone.o
	${LINK.f90} ${SEARCH} -o ${SEPEXE} ${OBJECTS} sp_stand_alone.o\
	${SHAREDIR}/io_unit_for77.o ${SHAREDIR}/ModIoUnit.o 

io_unit_for77.f: 
	cp  ${SHAREDIR}/io_unit_for77.f .

stdout.h:
	cp  ${SHAREDIR}/stdout.h .

SEP77:	io_unit_for77.f stdout.h
	f77 read_ihdata.f sp_stand_alone.f io_unit_for77.f main.f sep19.f

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
