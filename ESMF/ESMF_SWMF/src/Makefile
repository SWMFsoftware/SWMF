# by default, this makefile builds the executable 
ALL:  EXE

# Get the SWMF definitions
include ../../../Makefile.def

# The executable
EXE = ${BINDIR}/ESMF_SWMF.exe

# A link to the executable in the run directory
EXELINK = ./ESMF_SWMF.exe

### LOCDIR	  = src/Superstructure/AppDriver/seq_pairwise_fdriver_spmd

# The ESMF definitions
include $(ESMF_DIR)/makefile

DIRS = 

# Define these for the generic ESMF clean target
CLEANDIRS   =
CLEANFILES  = ${EXE}
CLOBBERDIRS =

# object files which need to be compiled into the executable.
APPOBJS = ESMF_SWMF_Coupler.o SWMF_GridComp.o ESMF_GridComp.o \
          ESMF_SWMF_GridComp.o ESMF_SWMF_Driver.o

EXE: ${EXE}
	@echo ' '
	@echo '${EXE} has been brought up to date'
	@echo ' '

# build the executable.
${EXE}: $(APPOBJS) ${LIBDIR}/libSWMF.a ${ESMFLIB}
	$(FLINKER) $(LINKOPTS) -o ${EXE} $(APPOBJS) $(FLINKLIBS) \
		-L${LIBDIR} -lSWMF

# run the executable in parallel. 
NP = 2

mpirun:	${EXE}
	@cd ../run; rm -f *.ESMF_LogFile; \
	${MPIRUN} -np ${NP} ${EXELINK}; \
	perl -ne 'print if /error/i' PET*.ESMF_LogFile

# run only a single copy.
run:	${EXE}
	@cd ../run; rm -f *.ESMF_LogFile; \
	${MPIRUN} -np 1 ${EXELINK}; \
	perl -ne 'print if /error/i' PET0.ESMF_LogFile
