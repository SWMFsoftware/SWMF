include ../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${DATAREADINDICESDIR} -I${EMPIRICALIEDIR} -I${EMPIRICALUADIR}

include ../Makefile.conf

include Makefile.DEPEND

include Makefile.RULES

MODULES = \
	ModGrid.o\
	geomag.o\
	iri90.o\
	ModSeProduction.o\
	ModSeCross.o\
	ModBackground.o\
	ModElecTrans.o\
	twostream_init.o\
	ModCouplePWOMtoTWOSTREAM.o



METHODS = \

LIBOBJECTS = ${MODULES} ${METHODS}

EXEOBJECTS = main_twostream.o ${LIBOBJECTS}

DEPEND:
	@${SCRIPTDIR}/depend.pl ${SEARCH} ${EXEOBJECTS} 


MY_LIB = ${LIBDIR}/libSE.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo



${MY_LIB}: ${LIBOBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${LIBOBJECTS}

EXE = ${BINDIR}/TWOSTREAM.exe

TWOSTREAM:	DEPEND
	@make ${EXE}
	@echo "TWOSTREAM.exe has been created"

${EXE}:	${EXEOBJECTS}
	@(if [ ${OS} = "Darwin" ];         then make TWOSTREAM_for_Darwin; \
	else                                    make TWOSTREAM_general; fi)

TWOSTREAM_general:
	${LINK.f90} -o ${EXE} ${EXEOBJECTS} \
		-L${LIBDIR}  -lEMPIRICALIE -lEMPIRICALUA -lINDICES -lSHARE -lTIMING \
		${Lflag1}

TWOSTREAM_for_Darwin:
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ${LIBDIR}/libEMPIRICALIE.a; \
		ar -x ${LIBDIR}/libEMPIRICALUA.a; \
		ar -x ${LIBDIR}/libINDICES.a; \
		ar -x ${LIBDIR}/libSHARE.a
	${LINK.f90} -o ${EXE} ${EXEOBJECTS} Tmp_/*.o ${LBLAS} ${Lflag1}
	rm -rf Tmp_

clean: cleanfiles

distclean: clean

