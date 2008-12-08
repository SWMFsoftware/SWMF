#^CFG COPYRIGHT UM
SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND

#			
#	Menu of make options
#
help:	
	@echo ' '
	@echo '  You can "make" the following:'
	@echo ' '
	@echo '    LIB          (libIONO.a - ionosphere library)'
	@echo ' '
	@echo '    help         (Makefile option list)'
	@echo '    clean        (rm -f *~ *.o *.kmo *.mod *.T *.lst core)'
	@echo '    distclean    (make clean; rm -f Makefile.DEPEND)'
	@echo ' '


MODULES = \
	ModProcIE.o\
	ModSizeRIM.o\
	ModRIM.o\
	ModIoRIM.o\
	ModParamRIM.o\
	ModFiles.o\
	ModIEInterface.o

#
#	Object Lists
#

OBJECTS = \
	IE_wrapper.o\
	init_RIM.o\
	stop_RIM.o\
	advance_RIM.o\
	get_conductance.o\
	conductance_gradients.o\
	solve.o\
	set_imf.o\
	gather.o\
	distribute.o\
	write_output_RIM.o\
	calc_aurora.o\
	IE_couplers.o\
	logfileRIM.o\
	set_test.o

#
#	Libraries
#

SEARCH_EXTRA = -I${LIBRARYDIR} -I${COUPLERDIR} -I${DATAREADINDICESDIR} -I${EMPIRICALIEDIR} 

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${MODULES} ${OBJECTS}

MY_LIB = ${LIBDIR}/libIE.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}:	${MF} ${MODULES} ${OBJECTS} ${MAKEFILE_COMP_SELECT}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES} ${OBJECTS}

#       Cleaning

distclean: clean


# keep this line
