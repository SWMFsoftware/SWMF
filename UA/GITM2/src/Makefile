
ABDIR   = ../srcSphereAB
EDDIR   = ../srcRudy
IEDIR   = ../srcIE
IODIR   = ../srcIO

SEARCH_EXTRA = -I${ABDIR} -I${IODIR}

include ../srcMake/Makefile.conf
include Makefile.DEPEND
include Makefile.RULES

MODULES = \
	ModSize.o\
	ModTime.o\
	ModConstants.o\
	ModPlanet.o\
	ModGITM.o\
	ModSources.o\
	ModVertical.o\
	ModSphereInterface.o\
	ModInputs.o\
	ModTimingGITM.o\
	ModRates.o\
	ModEUV.o\
	ModElectrodynamics.o\
	ModSatellites.o

OBJECTS = \
	planet.o\
	satellites.o\
	calc_electron_temperature.o\
	calc_etemp_sources.o\
	aurora.o\
	get_potential.o\
	calc_efield.o\
	set_inputs.o\
	init_energy_deposition.o\
	init_grid.o\
	time_routines.o\
	vertical_solver.o\
	library.o\
	stretch_grid.o\
	init_mpi.o\
	restart.o\
	initialize.o\
	init_msis.o\
	init_iri.o\
	init_b0.o\
	init_altitude.o\
	finalize.o\
	calc_neutral_friction.o\
	lubksb.o\
	ludcmp.o\
	calc_electrodynamics.o\
	divergence.o\
	calc_pressure.o\
	calc_rates.o\
	calc_timestep.o\
	exchange_messages_sphere.o\
	timing.o\
	advance.o\
	advance_vertical.o\
	advance_horizontal.o\
	set_vertical_bcs.o\
	output.o\
	logfile.o\
	calc_conduction.o\
	calc_physics.o\
	calc_euv.o\
	calc_neutral_density.o\
	calc_ion_density.o\
	calc_chemistry.Earth.o\
	calc_ion_v.o\
	chapman_new.o\
	stop_file.o\
	read_inputs.o\
	add_sources.o\
	calc_sources.o\
	msis_library.o\
	iri90.o\
	apex.o\
	write_output.o\
	gradient.o

DEPEND:
	@${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS} ${MODULES}

GITM:	DEPEND
	@make GITM.exe
	@echo "GITM.exe has been created"

GITM.exe:	${OBJECTS} ${MODULES} main.o
	${LINK.f90} ${Lflag1} -o GITM.exe main.o ${OBJECTS} ${MODULES}\
	-L${ABDIR} -lSphere -L${IODIR} -lindices -L${IEDIR} -lIE\
	-L${EDDIR} -lED -L. -lSHARE ${IEDIR}/ModConductance.o

STATIC:
	cp -f ModGITM_static.f90    ModGITM.f90
	cp -f ModEUV_static.f90     ModEUV.f90
	cp -f ModSources_static.f90 ModSources.f90

S2D = ${SCRIPTDIR}/StaticToDynamic.pl
DYNAMIC:	
	${S2D} -l=nBlocksMax -d=4 ModGITM_static.f90    > ModGITM.f90
	${S2D} -l=nBlocksMax -d=4 ModEUV_static.f90     > ModEUV.f90
	${S2D} -l=nBlocksMax -d=4 ModSources_static.f90 > ModSources.f90
	# Replaces nBlocksMax with nBlocks to make allocation flexible
	perl -pi -e 's/nBlocksMax/nBlocks/ if /allocate/' \
		ModGITM.f90 ModEUV.f90 ModSources.f90

MY_LIB = libGITM.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${LIBPREV} ${MF} ${MODULES} ${OBJECTS}
	cp -f ${LIBPREV} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS} ${MODULES}

distclean:
	rm -f 	ModGITM.f90 \
		ModEUV.f90 ModSources.f90 \
		ModPlanet.f90 \
		calc_electron_temperature.f90 \
		calc_ion_density.f90 \
		calc_neutral_density.f90 \
		calc_rates.f90 \
		init_iri.f90 \
		init_msis.f90 \
		output.f90 \
		planet.f90 \
		Makefile.DEPEND Makefile.RULES
