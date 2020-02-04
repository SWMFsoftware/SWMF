#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
-include Makefile.DEPEND

OBJECTS = sep19.o main.o 

SEPEXE = ${BINDIR}/SEP.exe

SEP: ${SEPEXE}

${SEPEXE}:	${OBJECTS} sp_stand_alone.o io_unit_for77.o
	${LINK.f90} ${SEARCH} -o ${SEPEXE} ${OBJECTS}  \
	sp_stand_alone.o io_unit_for77.o

SEP77:	stdout.h
	f77 sp_stand_alone.f main.f sep19.f io_unit_for77.f

MY_LIB = libKota.a

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} \
		${OBJECTS}

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

clean: cleanfiles

distclean: clean
	rm -f Makefile.DEPEND
