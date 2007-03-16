include ../Makefile.def
include ../Makefile.conf
include Makefile.planet
include Makefile.DEPEND
include Makefile.RULES

MODULES = \
	ModParameters.o\
	ModCommonPlanet.o\
	ModCommonVariables.o\
	ModPW.o\
	ModNeutralPW.o\
	ModPass.o

METHODS = \
	pw.o\
	polar_wind.o\
	set_parameters.o\
	initialize.o\
	solver.o\
	cvmgp.o\
	print_log.o\
	print_plot.o\
	EfieldPW.o\
	eHeatFluxPW.o\
	calcdt.o\
	diffusion.o\
	PRNT_Sources.o\
	calc_cfl.o\
	calc_diffusion.o\
	solver_rusanov.o\
	Get_Neutral.o\
	get_electrodynamic.o\
	print_electrodynamics.o\
	write_restart.o\
	${OBJECTS_PLANET}

LIBOBJECTS = ${MODULES} ${METHODS}

EXEOBJECTS = main.o ${LIBOBJECTS}

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

main.o: ${LIBSHARE}

EXE = ${BINDIR}/PWOM.exe

PWOM: DEPEND
	@make ${EXE}
	@echo "PWOM.exe has been created"

${EXE}: ${EXEOBJECTS}
	${LINK.f90} ${Lflag1} -o ${EXE} ${EXEOBJECTS} -L${LIBDIR} -lSHARE

distclean: clean
	rm -f PLANET *_planet.f ModCommonPlanet.f90 Makefile.planet 
