#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

include ../../../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${MAGNETOGRAMDIR}

include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libEMPIRICALEE.a

include Makefile.DEPEND
include Makefile.RULES

OBJECTS = \
	EEE_ModCommonVariables.o \
	EEE_ModMain.o \
	EEE_ModTD99.o \
	EEE_ModGL98.o \
	EEE_ModArch.o \
	EEE_ModGetB0.o \
	EEE_ModShearFlow.o \
	EEE_ModCms.o

OBJECTS_EXE= flux_rope_magnetogram.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${OBJECTS} ${OBJECTS_EXE}

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

LIBSHARE  = ${LIBDIR}/libSHARE.a
LIBTIMING = ${LIBDIR}/libTIMING.a
LIBMAGNETOGRAM = ${LIBDIR}/libMAGNETOGRAM.a

flux_rope_magnetogram.o : ${LIBSHARE} ${LIBTIMING} ${LIBMAGNETOGRAM} ${MY_LIB}

FRM:	DEPEND
	make ${BINDIR}/FRMAGNETOGRAM.exe

${BINDIR}/FRMAGNETOGRAM.exe: flux_rope_magnetogram.o
			rm -rf Tmp_; mkdir Tmp_ ;
			cd Tmp_; ar -x ${MY_LIB}; ar -x ${LIBTIMING}; \
			ar -x ${LIBSHARE}; ar -x ${LIBMAGNETOGRAM}; 
			${LINK.f90} -o ${BINDIR}/FRMAGNETOGRAM.exe \
			flux_rope_magnetogram.o Tmp_/*.o ${Lflag1}    ; 
			rm -rf Tmp_

clean: cleanfiles

distclean: clean
	rm -f Makefile.DEPEND
