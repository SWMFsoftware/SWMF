#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND

# These modules provide the "infrastructure" for SWMF and components.
# A component does not have to use the whole library.
# The minimum requirement is ModIoUnit.o

MODULES = \
	ModKind.o\
	ModTimeConvert.o\
	ModFreq.o\
	ModCoordTransform.o\
	ModSort.o\
        ModNumConst.o\
	ModConst.o\
	ModIoUnit.o\
	ModMpi.o\
	ModReadParam.o\
	ModUtilities.o\
	CON_axes.o\
	CON_geopack.o\
	CON_planet.o\
	CON_planet_field.o\
	CON_ray_trace.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${MODULES}

MY_LIB = ${LIBDIR}/libSHARE.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES}

install: 

distclean: clean
	rm -f Makefile.DEPEND 

# keep this line

