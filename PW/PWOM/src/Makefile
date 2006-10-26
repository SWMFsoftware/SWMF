include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND
include Makefile.RULES

include Makefile.planet

MODULES = \
	ModParameters.o\
	ModCommonPlanet.o\
	ModCommonVariables.o\
	ModPW.o\
	ModNeutralPW.o\
	ModPass.o

METHODS = \
	polar_wind.o\
	set_parameters.o\
	solver.o\
	cvmgp.o\
	prntel.o\
	prntef.o\
	prnthp.o\
	prntop.o\
	prntsm.o\
	elfldw.o\
	clfme2w.o\
	calcdt.o\
	diffusion.o\
	PRNT_Sources.o\
	calc_cfl.o\
	calc_diffusion.o\
	solver_rusanov.o\
	Get_Neutral.o\
	${OBJECTS_PLANET}

LIBOBJECTS = ${MODULES} ${METHODS}

EXEOBJECTS = pw.o ${LIBOBJECTS}

DEPEND:
	@${SCRIPTDIR}/depend.pl ${SEARCH} ${EXEOBJECTS}

MY_LIB = libPWOM.a

LIB:    DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${LIBOBJECTS} ${MAKEFILE_COMP_SELECT}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${LIBOBJECTS}

pw.o: ${LIBSHARE} ${MY_LIB}

EXE = ${BINDIR}/PWOM.exe

PWOM: DEPEND
	@make ${EXE}
	@echo "PWOM.exe has been created"

${EXE}: ${EXEOBJECTS}
	${LINK.f90} ${Lflag} -o ${EXE} ${EXEOBJECTS} -L${LIBDIR} -lSHARE

distclean: clean
	rm -f PLANET *_planet.f ModCommonPlanet.f90 Makefile.planet 
	rm -f Makefile.DEPEND Makefile.RULES

