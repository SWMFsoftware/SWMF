#^CMP COPYRIGHT UM

SHELL =/bin/sh

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND


SEARCH_EXTRA = -I${COUPLERDIR} -I${LIBRARYDIR}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${MODULES}

#^CMP IF GM BEGIN
MOD01 = CON_couple_gm_ie.o	#^CMP IF IE
MOD03 = CON_couple_ih_gm.o	#^CMP IF IH
MOD04 = CON_couple_gm_im.o	#^CMP IF IM
MOD09 = CON_couple_gm_rb.o	#^CMP IF RB
MOD11 = CON_couple_gm_pw.o      #^CMP IF PW
#^CMP END GM
#^CMP IF IE BEGIN
MOD05 = CON_couple_ie_im.o	#^CMP IF IM
MOD06 = CON_couple_ie_ua.o	#^CMP IF UA
MOD10 = CON_couple_ie_pw.o	#^CMP IF PW
MOD12 = CON_couple_ie_ps.o	#^CMP IF PS
MOD13 = CON_couple_ie_rb.o	#^CMP IF RB
#^CMP END IE
#^CMP IF IH BEGIN
MOD07 = CON_couple_ih_sc.o       #^CMP IF SC
MOD14 = CON_couple_ih_oh.o       #^CMP IF OH
#^CMP END IH

MOD08 = CON_couple_mh_sp.o       #^CMP IF SP
MOD15 = CON_couple_sc_lc.o       #^CMP IF LC 

MODULES = ${MOD01} ${MOD03} ${MOD04} ${MOD05} ${MOD06} ${MOD14} ${MOD15}\
	${MOD07} ${MOD08} ${MOD09} ${MOD10} ${MOD11} ${MOD12} ${MOD13}\
	CON_wrapper.o \
	CON_couple_all.o

MY_LIB = ${LIBDIR}/libINTERFACE.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES} ${MAKEFILE_COMP_SELECT}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES}

distclean: clean
	rm -f Makefile.DEPEND
