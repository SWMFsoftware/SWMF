
SHELL =/bin/sh

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND



OBJECTS = \
	EE_wrapper.o \
	ModBuffer.o \
	update_lagrangian_grid.o \
	ModRadioWaveImage.o \
	ModGridDescriptor.o \
	ModRadioWaveRaytracing.o \
	ModDensityAndGradient.o

SEARCH_EXTRA = -I${COUPLERDIR} -I${LIBRARYDIR} -I${MAGNETOGRAMDIR} \
	-I${EMPIRICALEEDIR} -I../src -I../srcBATL

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

#
#	Making executables and libraries
#
MY_LIB = ${LIBDIR}/libEE.a

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
	rm -f   Makefile.DEPEND ModBuffer.f90 update_lagrangian_grid.f90 \
	ModRadioWaveImage.f90 ModGridDescriptor.f90 \
	ModRadioWaveRaytracing.f90 ModDensityAndGradient.f90

# keep this line
