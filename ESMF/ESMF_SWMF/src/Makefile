# by default, this makefile builds the executable 
ALL:  EXE

# Get the SWMF definitions
-include ../../../Makefile.def
-include ../../../Makefile.conf

# Get the dependencies
-include Makefile.DEPEND

# Add ESMF include directories
-include ${ESMF_DIR}/esmf.mk
INCL_EXTRA=${ESMF_F90COMPILEPATHS}

# The executable
EXE = ${BINDIR}/ESMF_SWMF.exe

# A link to the executable in the run directory
EXELINK = ./ESMF_SWMF.exe

# object files which need to be compiled into the executable.
OBJECTS = \
	ESMF_SWMF_Coupler.o \
	SWMF_GridComp.o \
	ESMF_GridComp.o \
	ESMF_SWMF_GridComp.o \
	ESMF_SWMF_Driver.o \
	ESMF_SWMF_Mod.o

# Dependencies:
DEPEND:
	perl ${SCRIPTDIR}/depend.pl ${OBJECTS}

EXE: DEPEND
	make ${EXE}
	@echo ' '
	@echo '${EXE} has been brought up to date'
	@echo ' '

# build the executable.
$(EXE): $(OBJECTS)
	${LINK.f90} -o ${EXE} $(OBJECTS) ${ESMF_F90LINKPATHS} $(ESMF_F90ESMFLINKLIBS) \
		-L${LIBDIR} -lSWMF

# run the executable in parallel. 
parallelrun: ${EXE}
	@cd ../run; rm -f *.ESMF_LogFile; \
	${MPIRUN} ${EXELINK}; \
	perl -ne 'print if /error/i' PET*.ESMF_LogFile

# run only a single copy.
serialrun: ${EXE}
	@cd ../run; rm -f *.ESMF_LogFile; \
	${SERIAL} ${EXELINK}; \
	perl -ne 'print if /error/i' PET0.ESMF_LogFile

clean: cleanfiles

distclean: clean

