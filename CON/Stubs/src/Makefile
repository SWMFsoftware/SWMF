#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

include ../../../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${LIBRARYDIR} -I${COUPLERDIR}

include ../../../Makefile.conf
include Makefile.DEPEND

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${MODULES}

MODULES = \
	CON_wrapper.o \
	CON_couple_all.o

MY_LIB = ${LIBDIR}/libINTERFACE.a

LIB:	DEPEND
	$(MAKE) ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES} ${MAKEFILE_COMP_SELECT}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES}

distclean: clean
	rm -f Makefile.DEPEND
