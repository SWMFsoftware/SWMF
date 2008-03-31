#^CFG COPYRIGHT UM
SHELL =/bin/sh

include ../../../Makefile.def
include ../../../Makefile.conf

#			
#	Menu of make options
#
help:	
	@echo ' '
	@echo '  You can "make" the following:'
	@echo ' '
	@echo '    LIB          (libPS.a - ionosphere library)'
	@echo '    dgcpm        (Stand Alone Code)'
	@echo ' '
	@echo '    help         (Makefile option list)'
	@echo '    clean        (rm -f *~ *.o *.kmo *.mod *.T *.lst core)'
	@echo '    distclean    (make clean; rm -f Makefile.DEPEND)'
	@echo '    backup       (tar -cf IONO.`date +%m_%d_%y_%H%M.tar` ...)'
	@echo ' '


MODULES = \
	ModConstants.o\
	ModSizeDGCPM.o\
	ModMainDGCPM.o\
	ModIoDGCPM.o\
	ModHeidiDGCPM.o\
	ModDGCPMWaves.o\
	ModProcPS.o

OBJECTS = \
	dgcpm_setup_011.o\
	dgcpm_coefficients_010.o\
	dgcpm_output_010.o\
	DGCPM_read_inputs.o\
	pbo.o\
	stop_dgcpm.o\
	DGCPM_read_inputs.o\
	PS_wrapper.o


#
#	Libraries
#

SEARCH_EXTRA = -I${LIBRARYDIR} -I${COUPLERDIR}

#DEPEND:
#	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${MODULES} ${OBJECTS}

MY_LIB = ${LIBDIR}/libPS.a

LIB:	
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}:	${MODULES} ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES} ${OBJECTS}

#
#	Stand-alone executables
#

dgcpm:
	make dgcpm.exe

dgcpm010:
	make dgcpm.exe

dgcpm.exe:	${MODULES} ${OBJECTS_010} dgcpm_main_010.o
	$(COMPILE.f90) $(LFLAGS) -o dgcpm.exe ${MODULES} ${OBJECTS_010} dgcpm_main_010.o

distclean: clean
	rm -f Makefile.RULES Makefile.DEPEND
