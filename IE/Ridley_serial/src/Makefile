#^CFG COPYRIGHT UM
SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND
include Makefile.RULES

#			
#	Menu of make options
#
help:	
	@echo ' '
	@echo '  You can "make" the following:'
	@echo ' '
	@echo '    LIB          (libIONO.a - ionosphere library)'
	@echo '    IONO         (Stand-Alone Ionosphere Model for Earth)'
	@echo '    PostIONO.exe (convert IDL files to TecPlot format)'
	@echo ' '
	@echo '    help         (Makefile option list)'
	@echo '    clean        (rm -f *~ *.o *.kmo *.mod *.T *.lst core)'
	@echo '    distclean    (make clean; rm -f Makefile.DEPEND)'
	@echo '    backup       (tar -cf IONO.`date +%m_%d_%y_%H%M.tar` ...)'
	@echo ' '


MODULES = \
	IE_ModMain.o\
	IE_ModIo.o\
	ModIonosphere.o\
	ModProcIE.o\
	precision.o \
	io.o \
	CSR.o \
	MSR.o \
	LUCSRMSR.o \
	BandLU.o \
	matrix-arithmetic.o \
	flop_module.o \
	iteration-defaults.o \
	GMRES.o \
	ilu.o \
	ModAMIE_Interface.o\
	ModErrors.o\
	ModIE_Interface.o\
	ModExtras.o\
	ModFiles.o\
	ModConductance.o

#
#	Object Lists
#

OBJECTS = \
	IE_solve.o\
	IE_interpolate.o\
	IE_wrapper.o\
	IE_End.o\
	IE_Initialize.o\
	IE_Library.o\
	IO_Library.o\
	AMIE_Library.o\
	readAMIEoutput.o\
        iono_conductances.o\
        iono_solver.o\
        iono_coupling.o\
	Create_Region2_Currents.o\
        ionosphere.o\
	wei01.o\
	wei96.o\
	weicom.o\
	mh86.o\
	hmr89.o\
	iz94.o\
	ihp.o

#
#	Libraries
#

SEARCH_EXTRA = -I${LIBRARYDIR} -I${COUPLERDIR}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${SEARCH_EXTRA} \
		${MODULES} ${OBJECTS}

MY_LIB = ${LIBDIR}/libIE.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}:	${MF} ${MODULES} ${OBJECTS} ${MAKEFILE_COMP_SELECT}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES} ${OBJECTS}

#
#	Stand-alone executables
#

PostIONO.exe:	PostIONO.o
	$(LINK.f90) -o ${BINDIR}/PostIONO.exe $(Lflag2) PostIONO.o

IONO: ${EXE2}
	@echo Program IONO has been brought up to date.

backup:	
	tar -cf IONO.`date +%m_%d_%y_%H%M.tar` Make* *.f90 *.h *.sub *.script *.llscript.* *.in mpi* perf_mjp.a BATio BATioPP* IOcat*

distclean: clean
	rm -f Makefile.DEPEND Makefile.RULES

# keep this line
