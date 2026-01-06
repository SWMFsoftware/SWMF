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
ifeq ($(ipe_ROOT),)
OBJECTS = \
	ESMF_driver.o \
	ESMFSWMF_variables.o \
	ESMF_grid_comp.o \
	SWMF_grid_comp.o \
	IPE_grid_comp.o \
	RIM_grid_comp.o \
	CON_cpl_comp.o
else
OBJECTS = \
	ESMF_driver.o \
	ESMFSWMF_variables.o \
	ESMF_grid_comp.o \
	SWMF_grid_comp.o \
	RIM_grid_comp.o \
	CON_cpl_comp.o
INCL_EXTRA += -I$(ipe_ROOT)/include
IPELIBS = -L$(ipe_ROOT)/lib -lipe -lipm -lshared \
          -L$(comio_ROOT)/lib -lcomio \
          -L$(parallel_netcdf_ROOT)/lib -lpnetcdf \
          -L$(hdf5_ROOT)/lib -lhdf5_fortran
Cflag2 += -DEXTERNAL_IPE
endif

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
		${ESMF_F90LINKOPTS} -L${LIBDIR} \
		-lSWMF ${LflagMpi} ${IPELIBS}

info:
	@echo ESMF_F90LINKPATHS=${ESMF_F90LINKPATHS}
	@echo ESMF_F90ESMFLINKLIBS ${ESMF_F90ESMFLINKLIBS}
	@echo ESMF_F90LINKOPTS=${ESMF_F90LINKOPTS}
	@echo LflagMpi=${LflagMpi}
	@echo IPELIBS=${IPELIBS}

clean: cleanfiles

distclean: clean

