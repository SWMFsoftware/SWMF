# by default, this makefile builds the executable 
ALL:  EXE

# Get the SWMF definitions
-include ../../../Makefile.def
-include ../../../Makefile.conf

# Get the dependencies
-include Makefile.DEPEND

# Add ESMF include directories
-include ${ESMFMKFILE}
INCL_EXTRA = ${ESMF_F90COMPILEPATHS}

# The executable
EXE = ${BINDIR}/ESMF_SWMF.exe

# object files which need to be compiled into the executable.
OBJECTS = \
	ESMF_driver.o \
	ESMFSWMF_variables.o \
	ESMF_grid_comp.o \
	SWMF_grid_comp.o \
	IPE_grid_comp.o \
	RIM_grid_comp.o

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
	${LINK.f90} -o ${EXE} $(OBJECTS) \
		${ESMF_F90LINKPATHS} $(ESMF_F90ESMFLINKLIBS) \
		-L${LIBDIR} -lSWMF ${LflagMpi}

clean: cleanfiles

distclean: clean

