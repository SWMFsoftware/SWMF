#^CMP COPYRIGHT UM

SHELL =/bin/sh

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND


SEARCH_EXTRA = -I${COUPLERDIR} -I${LIBRARYDIR}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${SEARCH_EXTRA} \
		${MODULES} ${OBJECTS}

#^CMP IF GM BEGIN
MOD1 = CON_couple_gm_ie.o	#^CMP IF IE
MOD2 = CON_couple_gm_ie_swmf.o	#^CMP IF IE
MOD3 = CON_couple_ih_gm.o	#^CMP IF IH
MOD4 = CON_couple_gm_im.o	#^CMP IF IM
#^CMP END GM
#^CMP IF IE BEGIN
MOD5 = CON_couple_ie_im.o	#^CMP IF IM
MOD6 = CON_couple_ie_ua.o	#^CMP IF UA
#^CMP END IE

MODULES = ${MOD1} ${MOD2} ${MOD3} ${MOD4} ${MOD5} ${MOD6} \
	CON_wrapper.o \
	CON_couple_all.o

MY_LIB = ${LIBDIR}/libINTERFACE.a

LIB:	${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES} ${MAKEFILE_COMP_DEF}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES}

distclean: clean
	rm -f Makefile.DEPEND
