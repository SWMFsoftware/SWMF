#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND
include Makefile.RULES

# These modules provide the "infrastructure" for SWMF and components.
# A component does not have to use the whole library.
# The minimum reqguirement is ModIoUnit.o

# This rule is useful while modifying SWMF_MPI_routines.dat
#ModMpiInterfaces.f90: SWMF_MPI_routines.dat
#	${SCRIPTDIR}/CreateModMpi.pl

MODULES = \
	${BLAS}\
	ModBlasLapack.o\
	ModConst.o\
	ModCoordTransform.o\
	ModFreq.o\
	ModInterpolate.o\
	ModIoUnit.o\
	ModKind.o\
	ModLinearSolver.o\
	ModMpiConstants.o\
	ModMpiInterfaces.o\
	ModMpi.o\
        ModNumConst.o\
	ModPlanetConst.o\
	ModReadParam.o\
	ModSort.o\
	ModTimeConvert.o\
	ModTriangulate.o\
	ModUtilities.o\
	CON_axes.o\
	CON_geopack.o\
	CON_line_extract.o\
	CON_planet.o\
	CON_planet_field.o\
	CON_ray_trace.o

OBJECTS=\
	heat_capacity.o \
	get_ionization.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${MODULES} ${OBJECTS}

MY_LIB = ${LIBDIR}/libSHARE.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES} ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES}


distclean: clean
	rm -f Makefile.DEPEND


