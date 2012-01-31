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
    ModSizeDGCPM.o\
    ModConstants.o\
    ModMainDGCPM.o\
    ModIoDGCPM.o\
    ModProcPS.o\
    ModTimeDGCPM.o\
    ModCoupleDGCPM.o\
    ModFunctionsDGCPM.o\
    ModIonoDGCPM.o

OBJECTS = \
    dgcpm_setup.o\
    dgcpm_coefficients.o\
    dgcpm_output.o\
    pbo.o\
    dgcpm_stop.o\
    dgcpm_test.o\
    PS_wrapper.o\
    dgcpm_log.o


#
#	Libraries
#

SEARCH_EXTRA = -I${LIBRARYDIR} -I${COUPLERDIR} -I${DATAREADINDICESDIR}

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
