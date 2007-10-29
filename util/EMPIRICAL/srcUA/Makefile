#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libEMPIRICALUA.a

#include Makefile.DEPEND
#include Makefile.RULES
	 
OBJECTS =\
	msis90.o\
	iri90.o

#       msis_library.o         [* this is the same file as msis90
#       init_iri.o & init_msis.o   [* this files in model's root dir] 


#DEPEND:
#	@perl ${SCRIPTDIR}/depend.pl -I${SHAREDIR} ${MODULES}
 

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
