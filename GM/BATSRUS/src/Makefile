#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../Makefile.def

SEARCH_EXTRA = \
	-I../srcBATL -I${MAGNETOGRAMDIR} -I${DEMTDIR} -I${EMPIRICALEEDIR} \
	-I${CRASHDIR}

include ../Makefile.conf
include Makefile.DEPEND
include Makefile.RULES

# CONFIGURABLE OBJECT FILES

OBJ1 = ModCT.o constrain_B.o #				  ^CFG IF CONSTRAINB

#                                                     ^CFG IF IMPLICIT BEGIN
OBJ2=   ModImplicit.o \
        ModImplHypre.o \
        ModRadDiffusion.o \
        ModHeatConduction.o \
        implicit.o\
	impl_interface.o\
	impl_jacobian.o\
	impl_matvec.o\
	impl_newton.o
#                                                     ^CFG END IMPLICIT

OBJ3= ModRaytrace.o ray_trace_new.o ray_trace.o ray_pass.o #^CFG IF RAYTRACE

OBJ4= get_im_pressure.o #			            ^CFG IF RCM

OBJ5= ModProject.o project_B.o proj_cg.o proj_bicgstab.o #  ^CFG IF PROJECTION

OBJ6= ModResistivity.o #                                  ^CFG IF DISSFLUX

OBJ7= clean_divb.o #					  ^CFG IF DIVBDIFFUSE

# OBJECT FILES FOR LIBRARY / EXECUTABLE

OBJECTS = \
	ModAdvance.o\
	ModAMR.o\
	ModB0.o\
	ModBatlInterface.o\
	ModBlockData.o\
	ModBoundaryCells.o\
	ModCalcSource.o\
	ModCellBoundary.o\
	ModConserveFlux.o\
	ModChromosphere.o \
	ModCoronalHeating.o \
	ModFaceBoundary.o \
	ModIonoVelocity.o\
	ModExtraVariables.o\
	ModSingleFluid.o\
	ModMultiFluid.o\
	ModMultiIon.o\
	ModEnergy.o\
	ModEquation.o\
	ModCharacteristic.o\
	ModCurrent.o\
	ModExpansionFactors.o \
	ModFaceFlux.o\
	ModFaceGradient.o\
	ModFaceValue.o\
	ModGeometry.o\
	ModGroundMagPerturb.o\
	ModGmGeoindices.o \
	ModHallResist.o\
	ModHeatFluxCollisionless.o\
	ModIO.o\
	ModLaserHeating.o\
	ModMain.o\
	ModMessagePass.o\
	ModNodes.o\
	ModParallel.o\
	ModPartSteady.o\
	ModPIC.o\
	ModPhysics.o\
	ModPointImplicit.o\
	ModProcMH.o\
	ModRadiativeCooling.o \
	ModRestartFile.o\
	ModSatelliteFile.o\
	ModSize.o\
	ModSolarwind.o \
	ModTimeStepControl.o \
	ModUserEmpty.o \
	ModUser.o \
        ModViscosity.o \
	user_interface.o \
	ModWaves.o \
	ModHdf5.o \
	BATS_methods.o \
	MH_set_parameters.o\
	advect_points.o\
	explicit.o\
	fix_axis_cells.o\
	library.o\
	load_balance.o\
	rotate.o\
	set_block_geometry.o\
	set_ICs.o\
	set_physics.o\
	update_states.o\
	update_states_MHD.o\
	write_logfile.o\
	write_plot_common.o\
	write_plot_idl.o\
	write_plot_tec.o\
	write_plot_los.o \
	write_plot_sph.o \
	write_plot_radiowave.o \
	write_progress.o\
	${OBJ1} ${OBJ2} ${OBJ3} ${OBJ4} ${OBJ5} ${OBJ6} ${OBJ7}

# Default user module
ModUser.f90:
	cp ModUserDefault.f90 ModUser.f90

############################################################################

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

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

${MY_LIB}: ${LIBBATL} ${OBJECTS} ${MAKEFILE_COMP_SELECT}
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
		${LBLAS} ${HDFLIB} ${Lflag1} ${HYPRELIB}
	rm -rf Tmp_

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
	rm -rf Tmp_

#########################################################################

#
#	cleaning
#

distclean: clean
	rm -f ModUser.f90.safe
	mv ModUser.f90 ModUser.f90.safe
	rm -f ModSize.f90 ModEquation.f90 ModHdf5.f90 ModImplHypre.f90
