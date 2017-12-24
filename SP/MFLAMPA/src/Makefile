#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh
include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND

# Object files included into the library
OBJECTS = \
	ModAdvance.o \
	ModDistribution.o \
	ModDiffusion.o \
	ModGrid.o \
	ModLogAdvection.o \
	ModMain.o \
	ModProc.o \
	ModReadMhData.o \
	ModRestart.o \
	ModSize.o \
	ModUnit.o \
	ModWrite.o   

# Object files not included into the library
OBJECTS_EXE = SP_stand_alone.o

EXE = ${BINDIR}/MFLAMPA.exe

SEARCH_EXTRA = -I${LIBRARYDIR} -I${COUPLERDIR}

MY_LIB = libMFLAMPA.a

# Other required libraries
LIBSHARE  = ${LIBDIR}/libSHARE.a
LIBTIMING = ${LIBDIR}/libTIMING.a

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} \
		${MODULES} ${OBJECTS}

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES} ${OBJECTS} \
	${MAKEFILE_COMP_SELECT}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

# Additional ibraries should be compiled first
${OBJECTS_EXE}: ${LIBSHARE} ${LIBTIMING} ${MY_LIB}


MFLAMPA:
	make ${EXE}
	@echo ' '
	@echo Program MFLAMPA has been brought up to date.
	@echo ' '  

${EXE}: ${OBJECTS_EXE}
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ../${MY_LIB}; \
		ar -x ${LIBTIMING}; \
		ar -x ${LIBSHARE}
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} Tmp_/*.o \
		${Lflag1}
	rm -rf Tmp_ 

clean: cleanfiles

distclean: clean
	rm -f ModSize.f90
