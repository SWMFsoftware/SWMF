 #^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND

# Configurable objects

OBJECTS = \
	ModGridDescriptor.o\
	ModBuffer.o\
	update_lagrangian_grid.o\
	ModRadioWaveImage.o\
	ModRadioWaveRaytracing.o\
	ModDensityAndGradient.o\
	IH_set_buffer_grid.o\
	IH_get_for_gm.o\
	IH_get_for_mh.o\
	IH_get_for_mh_with_xyz.o\
	IH_get_for_sp.o\
	IH_put_from_mh.o\
	IH_wrapper.o

SEARCH_EXTRA = -I${COUPLERDIR} -I${LIBRARYDIR} -I../src -I../srcBATL -I${MAGNETOGRAMDIR} -I${DEMTDIR}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

#
#	Making executables and libraries
#
MY_LIB = ${LIBDIR}/libIH.a

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
	rm -f IH_*.f90 ModGridDescriptor.f90 ModBuffer.f90 \
		update_lagrangian_grid.f90
	rm -f Makefile.DEPEND ModDensityAndGradient.f90 \
		ModRadioWaveImage.f90 ModRadioWaveRaytracing.f90

# keep this line
