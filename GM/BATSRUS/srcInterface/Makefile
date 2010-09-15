#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND

# Configurable objects

#^CMP IF IH BEGIN
OBJ_IH = GM_put_from_ih.o
#^CMP END IH
#^CMP IF IM BEGIN
OBJ1   = ModGmImCoupling.o  #^CFG IF RAYTRACE
OBJ_IM = $(OBJ1) GM_get_for_im.o GM_put_from_im.o
#^CMP END IM
#^CMP IF RB BEGIN
OBJ2   = ModGmRbCoupling.o  #^CFG IF RAYTRACE
OBJ_RB = $(OBJ2) GM_get_for_rb.o
#^CMP END RB
#^CMP IF IE BEGIN
OBJ_IE =GM_get_for_ie.o GM_put_from_ie.o ModIeGrid.o
#^CMP END IE
#^CMP IF PW BEGIN
OBJ_PW = GM_get_for_pw.o GM_put_from_pw.o ModPwGrid.o
#^CMP END PW

OBJECTS = \
	ModGridDescriptor.o\
	ModBuffer.o\
	GM_put_from_ih_buffer.o\
	GM_wrapper.o\
	update_lagrangian_grid.o\
	ModDensityAndGradient.o\
	ModRadioWaveImage.o\
	ModRadioWaveRaytracing.o\
	ModAbsorption.o\
	$(OBJ_IH) $(OBJ_IM) $(OBJ_RB) $(OBJ_IE) $(OBJ_PW)

SEARCH_EXTRA = -I${COUPLERDIR} -I${LIBRARYDIR} -I../src -I../srcBATL -I${CRASHDIR}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

#
#	Making executables and libraries
#
MY_LIB = ${LIBDIR}/libGM.a

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
	rm -f Makefile.DEPEND

# keep this line
