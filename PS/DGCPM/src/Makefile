#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf
SHELL =/bin/sh

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND

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
	@echo ' '


OBJECTS = \
    ModSizeDGCPM.o\
    ModConstants.o\
    ModMainDGCPM.o\
    ModIoDGCPM.o\
    ModProcPS.o\
    ModTimeDGCPM.o\
    ModCoupleDGCPM.o\
    ModFunctionsDGCPM.o\
    ModIonoDGCPM.o\
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

SEARCHDIR = -I${SHAREDIR} -I${LIBRARYDIR} -I${COUPLERDIR} -I${DATAREADINDICESDIR}

DEPEND:
	@${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${OBJECTS}

MY_LIB = ${LIBDIR}/libPS.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS} ${MAKEFILE_COMP_SELECT}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

#
#	Stand-alone executable
#

dgcpm: DEPEND
	make dgcpm.exe

dgcpm.exe: ${OBJECTS} dgcpm_main.o
	$(COMPILE.f90) $(LFLAGS) -o dgcpm.exe ${OBJECTS} dgcpm_main.o

distclean: clean
	rm -f Makefile.RULES Makefile.DEPEND
