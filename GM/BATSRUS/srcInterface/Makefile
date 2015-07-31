#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND

# Configurable objects

OBJ_IH = GM_couple_ih.o read_ih_buffer.o                    #^CMP IF IH
OBJ_IM = GM_couple_im.o	 				    #^CMP IF IM
OBJ_RB = GM_couple_rb.o                                     #^CMP IF RB
OBJ_IE = GM_couple_ie.o calc_inner_bc_velocity.o logvar_ionosphere.o    #^CMP IF IE
OBJ_PT = GM_couple_pt.o                                     #^CMP IF PT
OBJ_PC = GM_couple_pc.o                                     #^CMP IF PC
OBJ_PW = GM_couple_pw.o read_pw_buffer.o                    #^CMP IF PW

OBJECTS = \
	ModGridDescriptor.o\
	ModBuffer.o\
	GM_wrapper.o\
	update_lagrangian_grid.o\
	ModDensityAndGradient.o\
	ModRadioWaveImage.o\
	ModRadioWaveRaytracing.o\
	$(OBJ_IH) $(OBJ_IM) $(OBJ_RB) $(OBJ_IE) $(OBJ_PT) $(OBJ_PW) $(OBJ_PC)

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

${MY_LIB}: ${BATSRUS_LIB} ${OBJECTS}
	cp -f ${BATSRUS_LIB} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

clean: cleanfiles

distclean: clean
	rm -f Makefile.DEPEND

# keep this line
