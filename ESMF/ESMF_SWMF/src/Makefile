# by default, this makefile builds the executable 
ALL:  EXE

# Get the SWMF definitions
include ../../../Makefile.def

# Get the dependencies
include Makefile.DEPEND

# The executable
EXE = ${BINDIR}/ESMF_SWMF.exe

# A link to the executable in the run directory
EXELINK = ./ESMF_SWMF.exe

# The ESMF definitions
include $(ESMF_DIR)/makefile

# Overwrite the incorrect .f90.o rule (for ESMF 2.*):
.f90.o:
	$(FC) -c $(FC_MOD)$(ESMF_MODDIR) $(FOPTFLAGS) $(FFLAGS) \
		$(F_FREENOCPP) $<

DIRS = 

# Define these for the generic ESMF clean target
CLEANDIRS   =
CLEANFILES  = ${EXE}
CLOBBERDIRS =

# object files which need to be compiled into the executable.
APPOBJS = \
	ESMF_SWMF_Coupler.o \
	SWMF_GridComp.o \
	ESMF_GridComp.o \
	ESMF_SWMF_GridComp.o \
	ESMF_SWMF_Driver.o \
	ESMF_SWMF_Mod.o

# Dependencies:
DEPEND:
	perl ${SCRIPTDIR}/depend.pl ${APPOBJS}

EXE: DEPEND
	make ${EXE}
	@echo ' '
	@echo '${EXE} has been brought up to date'
	@echo ' '

# build the executable.
${EXE}: $(APPOBJS) ${LIBDIR}/libSWMF.a ${ESMFLIB}
	@(if [ "${OS}" = "Darwin" ];         then make EXE_Darwin;\
	else                                      make EXE_general; fi)

EXE_general:
	echo $(FLINKER) and  $(LINKOPTS) and  $(FLINKLIBS)
	$(FLINKER) $(LINKOPTS) -o ${EXE} $(APPOBJS) $(FLINKLIBS) \
		-L${LIBDIR} -lSWMF

# Link rule for Darwin
EXE_Darwin:
	ld -all_load -r -o LIBRARIES.o $(LIBDIR)/libSWMF.a \
	    $(ESMF_F90LINKPATHS) $(ESMF_F90ESMFLINKLIBS)
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) \
	                  $(LINKOPTS) -o ${EXE} LIBRARIES.o $(APPOBJS)  \


TESTOBJS = \
	test.o

test.exe: ${TESTOBJS} ${ESMFLIB}
	$(FLINKER) $(LINKOPTS) -o test.exe ${TESTOBJS} $(FLINKLIBS) \
		-L${LIBDIR}

test: test.exe
	${MPIRUN} test.exe

ESMFTESTOBJS = \
	user_model1.o \
	user_model2.o \
	user_coupler.o \
	ESMFTEST.o

esmftest.exe: ${ESMFTESTOBJS} ${ESMFLIB}
	$(FLINKER) $(LINKOPTS) -o esmftest.exe ${ESMFTESTOBJS} $(FLINKLIBS) \
		-L${LIBDIR}

esmftest: esmftest.exe
	${MPIRUN} esmftest.exe

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
	rm -f Makefile.DEPEND
