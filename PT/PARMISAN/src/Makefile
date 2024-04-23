#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh
include ../Makefile.def
include ../Makefile.conf
-include Makefile.DEPEND
-include Makefile.RULES

# Object files included into the library
OBJECTS = \
	ModGrid.o \
	ModMain.o \
	ModOriginPoints.o \
	ModProc.o \
	ModReadMhData.o \
	ModSize.o \
	ModTime.o

ModSize.f90: ModSize_orig.f90
	cp -f ModSize_orig.f90 ModSize.f90

# Object files not included into the library
OBJECTS_EXE = PT_stand_alone.o

EXE = ${BINDIR}/PARMISAN.exe

SEARCH_EXTRA = -I${LIBRARYDIR} -I${COUPLERDIR} -I${EMPIRICALCRDIR}

MY_LIB = libPARMISAN.a

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


PARMISAN:
	make ${EXE}
	@echo ' '
	@echo Program PARMISAN has been brought up to date.
	@echo ' '  

${EXE}: ${OBJECTS_EXE}
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ../${MY_LIB}; \
		ar -x ${LIBTIMING}; \
		ar -x ${LIBSHARE}
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} Tmp_/*.o \
		${LflagMpi}
	rm -rf Tmp_ 

clean: cleanfiles

distclean: clean
	rm -f ModSize.f90
