#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

include ../Makefile.def

SEARCHDIR = -I${SHAREDIR} \
	-I../srcBATL -I${MAGNETOGRAMDIR} -I${DEMTDIR} -I${EMPIRICALEEDIR} \
	-I${CRASHDIR}

include ../Makefile.conf
include Makefile.DEPEND
include Makefile.RULES

# OBJECT FILES FOR LIBRARY / EXECUTABLE

OBJECTS = \
	ModAdvance.o \
	ModAMR.o \
	ModB0.o \
	ModBatlInterface.o \
	ModBlockData.o \
	ModBoundaryCells.o \
	ModCalcSource.o \
	ModCellBoundary.o \
	ModCellGradient.o \
	ModCharacteristic.o \
	ModChromosphere.o \
	ModCoarseAxis.o \
	ModConserveFlux.o \
	ModCoronalHeating.o \
	ModCT.o \
	ModCurrent.o \
	ModFaceBoundary.o \
	ModFieldLineThread.o \
	ModImplicit.o \
	ModImplHypre.o \
	ModExtraVariables.o \
	ModHeatConduction.o \
	ModEnergy.o \
	ModElectricField.o \
	ModEquation.o \
	ModExpansionFactors.o \
	ModFaceFlux.o \
	ModFaceGradient.o \
	ModFaceValue.o \
	ModParticleFieldLine.o \
	ModGeometry.o \
	ModGroundMagPerturb.o \
	ModHallResist.o \
	ModHdf5.o \
	ModHeatFluxCollisionless.o \
	ModIeCoupling.o \
	ModIO.o \
	ModIonElectron.o \
	ModLaserHeating.o \
	ModLoadBalance.o \
	ModLocalTimeStep.o \
	ModMain.o \
	ModMessagePass.o \
	ModMultiFluid.o \
	ModMultiIon.o \
	ModNodes.o \
	ModParallel.o \
	ModPartImplicit.o \
	ModPartSteady.o \
	ModPIC.o \
	ModPhysics.o \
	ModPlotShell.o \
	ModPlotBox.o \
	ModPointImplicit.o \
	ModProcMH.o \
	ModProject.o \
	ModRadDiffusion.o \
	ModRadiativeCooling.o \
	ModRaytrace.o \
	ModRestartFile.o \
	ModResistivity.o \
	ModSatelliteFile.o \
	ModSemiImplicit.o \
	ModSingleFluid.o \
	ModSize.o \
	ModSolarwind.o \
	ModThreadedLC.o \
	ModTimeStepControl.o \
	ModUserEmpty.o \
	ModUserInterface.o \
	ModUser.o \
	ModViscosity.o \
	ModWaves.o \
	BATS_methods.o \
	MH_set_parameters.o \
	advect_points.o \
	clean_divb.o \
	constrain_B.o \
	explicit.o \
	fix_axis_cells.o \
	get_im_pressure.o \
	library.o \
	project_B.o \
	proj_bicgstab.o \
	proj_cg.o \
	ray_trace.o \
	ray_trace_new.o \
	ray_pass.o \
	rotate.o \
	set_block_geometry.o \
	set_ICs.o \
	update_states.o \
	update_states_MHD.o \
	user_interface.o \
	write_logfile.o \
	write_plot_common.o \
	write_plot_idl.o \
	write_plot_tec.o \
	write_plot_los.o \
	write_plot_particle.o \
	write_plot_radiowave.o \
	write_progress.o

# Default user module
ModUser.f90:
	cp ModUserDefault.f90 ModUser.f90

############################################################################

DEPEND:
	perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${OBJECTS}

#
#	Making executables and libraries
#

#
# Library for framework and also for stand alone code
#

MY_LIB  = libBATSRUS.a
LIBBATL = ../srcBATL/libBATL.a

LIB:	DEPEND
	$(MAKE) ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${LIBBATL} ${OBJECTS}
	cp -f ${LIBBATL} ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

# Object files not included into the library
OBJECTS_EXE = main.o

# Other requireed libraries
LIBSHARE  = ${LIBDIR}/libSHARE.a
LIBTIMING = ${LIBDIR}/libTIMING.a

# Libraries should be compiled first, because modules are used in main.
${OBJECTS_EXE}: ${LIBSHARE} ${LIBTIMING} ${MY_LIB}

EXE = ${BINDIR}/BATSRUS.exe

BATSRUS:
	$(MAKE) ${EXE}
	@echo ' '
	@echo Program BATSRUS has been brought up to date.
	@echo ' '	

${EXE}: ${OBJECTS_EXE}
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ../${MY_LIB}; \
		ar -x ${LIBDIR}/libMAGNETOGRAM.a; \
		ar -x ${LIBDIR}/libDEMT.a; \
		ar -x ${LIBDIR}/libEMPIRICALEE.a; \
		ar -x ${LIBTIMING}; \
		ar -x ${LIBSHARE}
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} Tmp_/*.o \
		${LBLAS} ${HDFLIB} ${Lflag1} ${HYPRELIB} ${SPICELIB}

# Compilation of CRASH ############

# Other required libraries
LIBCRASH  = ${LIBDIR}/libCRASH.a

# Libraries should be compiled first, because modules are used in main.
${CRASH_OBJECT_EXE}: ${LIBSHARE} ${LIBTIMING} ${LIBCRASH} ${MY_LIB}

CRASHEXE = ${BINDIR}/CRASH.exe

CRASH:
	$(MAKE) ${CRASHEXE}
	@echo ' '
	@echo Program CRASH has been brought up to date.
	@echo ' '	

${CRASHEXE}: ${OBJECTS_EXE} ${LIBCRASH}
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ../${MY_LIB}; \
		ar -x ${LIBDIR}/libMAGNETOGRAM.a; \
		ar -x ${LIBDIR}/libDEMT.a; \
		ar -x ${LIBDIR}/libEMPIRICALEE.a; \
		ar -x ${LIBTIMING}; \
		ar -x ${LIBCRASH}; \
		ar -x ${LIBSHARE}
	${LINK.f90} -o ${CRASHEXE} ${OBJECTS_EXE} Tmp_/*.o \
		${LBLAS} ${HDFLIB} ${Lflag1} ${HYPRELIB}

#########################################################################

#
#	cleaning
#

clean: cleanfiles
	rm -rf Tmp_

distclean: clean
	rm -f ModUser.f90.safe
	mv ModUser.f90 ModUser.f90.safe
	rm -f ModSize.f90 ModEquation.f90 ModHdf5.f90 ModImplHypre.f90
