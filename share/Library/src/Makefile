#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND
include Makefile.RULES

# These modules provide the "infrastructure" for SWMF and components.
# A component does not have to use the whole library.
# The minimum reqguirement is ModIoUnit.o

# This rule is used to regenerate ModMpiInterface.f90 from a modified SWMF_MPI_routines.dat
MPIINTERFACE:
	${SCRIPTDIR}/CreateModMpi.pl

OBJECTS = \
	${BLAS}\
	ModBlasLapack.o\
	ModConst.o\
	ModCoordTransform.o\
	ModExactRS.o\
	ModFreq.o\
	ModHyperGeometric.o\
	ModInitialState.o\
	ModInterpolate.o\
	ModInterpolateAMR.o\
	ModIoUnit.o\
	ModHdf5Utils.o\
	ModKind.o\
	ModLinearAdvection.o\
	ModLinearSolver.o\
	ModLookupTable.o\
	ModMpiInterfaces.o\
	ModMpi.o\
	ModMpiOrig.o\
        ModNumConst.o\
	ModPlanetConst.o\
	ModPlotFile.o\
	ModProcessVarName.o\
	ModRandomNumber.o\
	ModReadParam.o\
	ModSort.o\
	ModSpice.o\
	ModTimeConvert.o\
	ModTriangulate.o\
	ModTriangulateSpherical.o\
	ModUtility.o\
	CON_axes.o\
	CON_geopack.o\
	CON_line_extract.o\
	CON_planet.o\
	CON_planet_field.o\
	CON_ray_trace.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${OBJECTS} PostIDL.o

MY_LIB = ${LIBDIR}/libSHARE.a

LIB:	DEPEND
	${MAKE} ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS} ModUtility_c.o
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS} ModUtility_c.o

LIB_NO_C:DEPEND
	perl -pi -e 's/(iError = make_dir_c)/\!$$1/' ModUtility.f90
	${MAKE} LIB_NO_C_
	perl -pi -e 's/.(iError = make_dir_c)/$$1/' ModUtility.f90
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

LIB_NO_C_: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

${BINDIR}/PostIDL.exe: PostIDL.o ${MY_LIB}
	${LINK.f90} -o ${BINDIR}/PostIDL.exe PostIDL.o -L${LIBDIR} -lSHARE ${Lflag2}

PIDL: DEPEND
	make ${BINDIR}/PostIDL.exe

clean: cleanfiles

distclean: clean
	rm -f Makefile.DEPEND ModHdf5Utils.f90 ModSpice.f90
