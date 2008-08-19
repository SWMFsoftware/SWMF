#^CFG COPYRIGHT UM

SHELL =/bin/sh
include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND


OBJECTS = ModLogAdvection.o \
	ModLinAdvection.o \
	SP_ModTurbulence.o \
	SP_main.o \
	read_ihdata_for_sp.o \
	sp_diffusive_acceleration.o \
	write_logfile_SP.o \
	write_plotfile_SP.o 

SEPEXE = ${BINDIR}/SEP.exe

SEP: ${SEPEXE}

SEARCH_EXTRA = -I${LIBRARYDIR} -I${COUPLERDIR}



MY_LIB = libFLAMPA.a

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} \
		${MODULES} ${OBJECTS} SP_stand_alone.o

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES} ${OBJECTS} \
	${MAKEFILE_COMP_SELECT}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}


${SEPEXE}: ${MY_LIB} SP_stand_alone.o ${LIBDIR}/libSHARE.a
	@(if [ ${OS} = "Darwin" ];         then make SEP_for_Darwin; \
	else                                    make SEP_general; fi)

SEP_general:
	${LINK.f90} -o ${SEPEXE} ${OBJECTS_EXE} \
		-L. -lFLAMPA -L${LIBDIR} -lSHARE ${Lflag1}

SEP_for_Darwin:
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ../${MY_LIB}; \
		ar -x ${LIBDIR}/libSHARE.a
	${LINK.f90} -o ${SEPEXE} SP_stand_alone.o Tmp_/*.o \
		${Lflag1}
	rm -rf Tmp_




distclean: clean
	rm -f Makefile.DEPEND
