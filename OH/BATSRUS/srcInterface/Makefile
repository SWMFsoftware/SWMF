
SHELL =/bin/sh

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND



OBJECTS = \
	ModGridDescriptor.o\
        ModBuffer.o\
        update_lagrangian_grid.o\
        ModRadioWaveImage.o\
        ModRadioWaveRaytracing.o\
        ModDensityAndGradient.o\
        OH_set_buffer_grid.o\
        OH_get_for_mh_with_xyz.o\
	OH_wrapper.o

SEARCH_EXTRA = -I${COUPLERDIR} -I${LIBRARYDIR} -I../src -I../srcBATL

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

#
#	Making executables and libraries
#
MY_LIB = ${LIBDIR}/libOH.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

BATSRUS_LIB = ../src/libBATSRUS.a 

${MY_LIB}: ${BATSRUS_LIB} ${OBJECTS} ${MAKEFILE_COMP_SELECT}
	cp -f ${BATSRUS_LIB} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

distclean: clean
	rm -f  	ModGridDescriptor.f90\
	ModBuffer.f90 \
	update_lagrangian_grid.f90 \
	ModRadioWaveImage.f90 \
	ModRadioWaveRaytracing.f90 \
	ModDensityAndGradient.f90 \
	OH_get_for_mh_with_xyz.f90 \
	OH_wrapper.f90 Makefile.DEPEND  

# keep this line
