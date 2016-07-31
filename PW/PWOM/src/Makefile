include ../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${DATAREADINDICESDIR} -I${EMPIRICALIEDIR} -I${EMPIRICALUADIR}\
	    -I../srcSTET -I../srcTWOSTREAM

include ../Makefile.conf

include Makefile.planet
include Makefile.DEPEND

include Makefile.RULES

MODULES = \
	ModOvation.o\
	ModPW.o\
	ModParameters.o\
	ModCommonPlanet_planet.o\
	ModCommonVariables.o\
	ModNeutralPW.o\
	ModPass.o\
	ModPwImplicit.o\
	ModGmPressure.o\
	ModAurora.o\
	ModPwTime.o\
	ModLatLon.o\
	ModPwPlot.o\
	ModWaves.o\
	ModPhotoElectron.o\
	${MODULES_PLANET}

METHODS = \
	pw.o\
	polar_wind.o\
	set_parameters.o\
	initialize.o\
	solver.o\
	cvmgp.o\
	Efield.o\
	eheat_flux.o\
	eheat_flux_explicit.o\
	iheat_flux.o\
	calcdt.o\
	diffusion.o\
	calc_cfl.o\
	calc_diffusion.o\
	solver_rusanov.o\
	Get_Neutral.o\
	get_electrodynamic.o\
	print_electrodynamics.o\
	write_restart.o\
	set_upper_bc.o\
	get_weimer_potential.o\
	read_restart.o\
	set_vertical_grid.o\
	calc_dt.o\
	calc_centrifugal.o\
	${OBJECTS_PLANET}

LIBOBJECTS = ${MODULES} ${METHODS}

EXEOBJECTS = main.o ${LIBOBJECTS}

DEPEND:
	@${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${EXEOBJECTS} 


MY_LIB = libPWOM.a
LIBSE  = ${LIBDIR}/libSE.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${LIBSE} ${LIBOBJECTS}
	cp -f ${LIBSE} ${MY_LIB}
	${AR} ${MY_LIB} ${LIBOBJECTS}

main.o: ${LIBSHARE}

EXE = ${BINDIR}/PWOM.exe

PWOM:	DEPEND
	@make ${EXE}
	@echo "PWOM.exe has been created"

${EXE}:	${EXEOBJECTS}
	${LINK.f90} -o ${EXE} ${EXEOBJECTS} \
		-L${LIBDIR} -lEMPIRICALIE -lEMPIRICALUA -lINDICES  \
		-lTIMING  -lSE -lSHARE\
		${Lflag1}

PWOM_for_Darwin:
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ${LIBDIR}/libEMPIRICALIE.a; \
		ar -x ${LIBDIR}/libEMPIRICALUA.a; \
		ar -x ${LIBDIR}/libINDICES.a; \
		ar -x ${LIBDIR}/libSHARE.a
		ar -x ${LIBDIR}/libSE.a
	${LINK.f90} -o ${EXE} ${EXEOBJECTS} Tmp_/*.o ${LBLAS} ${Lflag1}
	rm -rf Tmp_

clean: cleanfiles

distclean: clean
	rm -f PLANET *_planet.f *_planet.f90
