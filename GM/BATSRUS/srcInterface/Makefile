#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND

# Configurable objects

OBJ_IH = GM_put_from_ih.o                                   #^CMP IF IH
OBJ_IM = ModGmImCoupling.o GM_get_for_im.o GM_put_from_im.o #^CMP IF IM
OBJ_RB = ModGmRbCoupling.o GM_get_for_rb.o                  #^CMP IF RB
OBJ_IE = GM_get_for_ie.o GM_put_from_ie.o ModIeGrid.o       #^CMP IF IE
OBJ_PT = GM_get_for_pt.o                                    #^CMP IF PT
OBJ_PW = GM_get_for_pw.o GM_put_from_pw.o ModPwGrid.o       #^CMP IF PW

OBJECTS = \
	ModGridDescriptor.o\
	ModBuffer.o\
	GM_put_from_ih_buffer.o\
	GM_wrapper.o\
	update_lagrangian_grid.o\
	ModDensityAndGradient.o\
	ModRadioWaveImage.o\
	ModRadioWaveRaytracing.o\
	$(OBJ_IH) $(OBJ_IM) $(OBJ_RB) $(OBJ_IE) $(OBJ_PT) $(OBJ_PW)

SEARCHDIR = -I${SHAREDIR} -I${COUPLERDIR} -I${LIBRARYDIR} -I${MAGNETOGRAMDIR} \
	-I../src -I../srcBATL -I${CRASHDIR}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${OBJECTS}

#
#	Making executables and libraries
#
MY_LIB = ${LIBDIR}/libGM.a

LIB:	DEPEND
	$(MAKE) ${MY_LIB}
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
