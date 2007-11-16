#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libEMPIRICALUA.a

#include Makefile.DEPEND
#include Makefile.RULES
	 
OBJECTS =\
	ModIri90.o\
	msis90.o

#	iri90.o\
#	ModMsis90.o\

#DEPEND:
#	@perl ${SCRIPTDIR}/depend.pl -I${SHAREDIR} ${MODULES}

#.f.o:
#	${COMPILE.f77} ${CFLAGS} ${PRECISION} $<
 
LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}:  ${MODULES}  ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS} ${MODULES}

distclean: clean
	rm -f Makefile.DEPEND
