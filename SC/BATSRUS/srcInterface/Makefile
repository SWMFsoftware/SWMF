
SHELL =/bin/sh

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND



OBJECTS = \
	ModBuffer.o\
	ModGridDescriptor.o\
	SC_get_for_mh.o\
	SC_get_for_mh_with_xyz.o\
	SC_get_for_sp.o\
	SC_put_from_mh.o\
	SC_wrapper.o\
	update_lagrangian_grid.o\
	ModDensityAndGradient.o\
	ModRadioWaveImage.o\
	ModRadioWaveRaytracing.o \
	SC_set_buffer_grid.o

SEARCH_EXTRA = -I${COUPLERDIR} -I${LIBRARYDIR} -I${MAGNETOGRAMDIR} \
	-I${EMPIRICALEEDIR} -I../src -I../srcBATL

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

#
#	Making executables and libraries
#
MY_LIB = ${LIBDIR}/libSC.a

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
	rm -f   SC_wrapper.f90 SC_get_for_sp.f90 SC_set_buffer_grid.f90\
		ModGridDescriptor.f90 update_lagrangian_grid.f90 \
		Makefile.DEPEND ModBuffer.f90 ModDensityAndGradient.f90 \
		ModRadioWaveImage.f90 ModRadioWaveRaytracing.f90 \
		SC_get_for_mh.f90 SC_get_for_mh_with_xyz.f90 \
		SC_put_from_mh.f90

# keep this line
