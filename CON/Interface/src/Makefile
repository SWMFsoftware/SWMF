#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

include ../../../Makefile.def
include ../../../Makefile.conf
-include Makefile.DEPEND

SEARCHDIR = -I${SHAREDIR} -I${COUPLERDIR} -I${LIBRARYDIR}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${MODULES}

#^CMP IF GM BEGIN
MOD01 = CON_couple_gm_ie.o	#^CMP IF IE
MOD02 = CON_couple_gm_ih.o	#^CMP IF IH
MOD03 = CON_couple_gm_im.o	#^CMP IF IM
MOD04 = CON_couple_gm_pc.o	#^CMP IF PC
MOD05 = CON_couple_gm_ps.o	#^CMP IF PS
MOD06 = CON_couple_gm_pt.o	#^CMP IF PT
MOD07 = CON_couple_gm_pw.o      #^CMP IF PW
MOD08 = CON_couple_gm_rb.o	#^CMP IF RB
MOD09 = CON_couple_gm_sc.o	#^CMP IF SC
MOD10 = CON_couple_gm_ua.o	#^CMP IF UA
#^CMP END GM
#^CMP IF IE BEGIN
MOD11 = CON_couple_ie_im.o	#^CMP IF IM
MOD12 = CON_couple_ie_ps.o	#^CMP IF PS
MOD13 = CON_couple_ie_pw.o	#^CMP IF PW
MOD14 = CON_couple_ie_rb.o	#^CMP IF RB
MOD15 = CON_couple_ie_ua.o	#^CMP IF UA
#^CMP END IE
#^CMP IF IH BEGIN
MOD16 = CON_couple_ih_sc.o	#^CMP IF SC
MOD17 = CON_couple_ih_oh.o	#^CMP IF OH
MOD18 = CON_couple_ih_pt.o      #^CMP IF PT 
#^CMP END IH
#^CMP IF SC BEGIN
MOD19 = CON_couple_ee_sc.o	#^CMP IF EE
MOD20 = CON_couple_sc_pt.o	#^CMP IF PT
#^CMP END SC
#^CMP IF OH BEGIN
MOD21 = CON_couple_oh_pt.o	#^CMP IF PT
#^CMP END OH

MOD22 = CON_couple_mh_sp.o	#^CMP IF SP

MODULES = \
	${MOD01} ${MOD02} ${MOD03} ${MOD04} ${MOD05} ${MOD06} ${MOD07} \
	${MOD08} ${MOD09} ${MOD10} ${MOD11} ${MOD12} ${MOD13} ${MOD14} \
	${MOD15} ${MOD16} ${MOD17} ${MOD18} ${MOD19} ${MOD20} ${MOD21} \
	${MOD22} \
	CON_wrapper.o \
	CON_couple_all.o

MY_LIB = ${LIBDIR}/libINTERFACE.a

LIB:	DEPEND
	$(MAKE) ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES}

clean: cleanfiles

distclean: clean
