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
OBJ_IM = GM_get_for_im.o GM_put_from_im.o
#^CMP END IM
#^CMP IF IE BEGIN
OBJ_IE =GM_get_for_ie.o GM_put_from_ie.o \
        GM_get_for_ie_swmf.o GM_put_from_ie_swmf.o \
	ModInnerBC.o \
	ionosphere_bc.o map_ionosphere_bc.o
#^CMP END IE

OBJECTS = \
	ModGridDescriptor.o\
	GM_wrapper.o\
	$(OBJ_IH) $(OBJ_IM) $(OBJ_IE)

SEARCH_EXTRA = -I${COUPLERDIR} -I${LIBRARYDIR} -I../src

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${SEARCH_EXTRA} \
		${MODULES} ${OBJECTS} ${OBJ_WRAPPER}

#
#	Making executables and libraries
#
MY_LIB = ${LIBDIR}/libGM.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

BATSRUS_LIB = ${LIBDIR}/libBATSRUS.a 

${MY_LIB}: ${BATSRUS_LIB} ${OBJECTS} ${MAKEFILE_COMP_SELECT}
	cp -f ${BATSRUS_LIB} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

distclean: clean
	rm -f Makefile.DEPEND

# keep this line
