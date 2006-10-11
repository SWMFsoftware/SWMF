include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND
# include Makefile.RULES


EXE = ../bin/pw3.exe

default: ${EXE}

include Makefile.planet

MODULES = \
	ModParameters.o\
	ModCommonPlanet.o\
	ModCommonVariables.o\
	Mod_PW.o\
	ModNeutralPW.o\
	ModPass.o

METHODS = \
	polar_wind.o\
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

EXE = ${BINDIR}/PWOM.exe

PWOM:
	@make ${EXE}
	@echo "PWOM.exe has been created"

${EXE}: ${EXEOBJECTS}
	${LINK.f90} ${Lflag} -o ${EXE} ${EXEOBJECTS}

distclean: clean
	rm -f PLANET *_planet.f ModCommonPlanet.f90 Makefile.planet 
	rm -f Makefile.DEPEND Makefile.RULES

