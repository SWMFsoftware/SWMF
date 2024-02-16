# by default, this makefile builds the executable 
ALL:  EXE

# Get the SWMF definitions
-include ../../../Makefile.def
-include ../../../Makefile.conf

# Get the dependencies
-include Makefile.DEPEND

# Add ESMF include directories
-include ${ESMF_DIR}/esmf.mk
INCL_EXTRA = ${ESMF_F90COMPILEPATHS}

# The executable
EXE = ${BINDIR}/ESMF_SWMF.exe

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
	${LINK.f90} -o ${EXE} $(OBJECTS) \
		${ESMF_F90LINKPATHS} $(ESMF_F90ESMFLINKLIBS) \
		-L${LIBDIR} -lSWMF

clean: cleanfiles

distclean: clean

