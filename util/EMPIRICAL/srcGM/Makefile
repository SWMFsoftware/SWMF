#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libEMPIRICALGM.a

#include Makefile.DEPEND
#include Makefile.RULES

OBJECTS = \
	t96_01.o\
	t04_s.o


#DEPEND:
#	@perl ${SCRIPTDIR}/depend.pl -I${SHAREDIR} ${MODULES}

.f.o:
	${COMPILE.f77} ${CFLAGS} ${PRECISION} $<

LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES} ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES} ${OBJECTS}

distclean: clean
	rm -f Makefile.DEPEND
