#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../Makefile.def

SEARCH_EXTRA = -I${MAGNETOGRAMDIR} -I${EMPIRICALEEDIR}

include ../Makefile.conf
include Makefile.DEPEND
include Makefile.RULES

# CONFIGURABLE MODULES

MOD1 = ModCT.o       #^CFG IF CONSTRAINB
MOD2 = ModImplicit.o #^CFG IF IMPLICIT
MOD3 = ModRaytrace.o #^CFG IF RAYTRACE
MOD4 = ModProject.o  #^CFG IF PROJECTION


MODULES = \
	ModAdvance.o\
	ModAMR.o\
	ModBlockData.o\
	ModBoundaryCells.o\
	ModConserveFlux.o\
	ModCube.o\
	ModSingleFluid.o\
	ModMultiFluid.o\
	ModMultiIon.o\
	ModEnergy.o\
	ModEquation.o\
	ModCharacteristic.o\
	ModCovariant.o\
	ModFaceFlux.o\
	ModFaceValue.o\
	ModGeometry.o\
	ModHallResist.o\
	ModIO.o\
	ModMain.o\
	ModMessagePass.o\
	ModNodes.o\
	ModOctree.o\
	ModParallel.o\
	ModPartSteady.o\
	ModPhysics.o\
	ModPointImplicit.o\
	ModPolarNeighbor.o\
	ModProcMH.o\
	ModRestartFile.o\
	ModSize.o\
	ModSolarwind.o \
	ModUserEmpty.o \
	ModUser.o \
	ModExpansionFactors.o \
	${MOD1}	${MOD2}	${MOD3} ${MOD4} 

# OBJECT FILES FOR MAIN EXECUTABLE

# CONFIGURABLE OBJECTS				^CFG IF CONFIGURE

OBJ1=constrain_B.o #				^CFG IF CONSTRAINB

#^CFG IF IMPLICIT BEGIN
OBJ2=   implicit.o\
	impl_interface.o\
	impl_jacobian.o\
	impl_matvec.o\
	impl_newton.o
#^CFG END IMPLICIT

OBJ3=ray_trace_new.o ray_trace.o ray_pass.o #	^CFG IF RAYTRACE
OBJ3B=get_im_pressure.o #			^CFG IF RCM
OBJ4=project_B.o proj_cg.o proj_bicgstab.o #	^CFG IF PROJECTION

OBJ10=  calc_heat_flux.o set_resistivity.o ModResistivity.o #^CFG IF DISSFLUX

OBJ15=clean_divb.o                          #^CFG IF DIVBDIFFUSE

OBJECTS = \
	BATS_methods.o \
	MH_set_parameters.o\
	advect_points.o\
	allocate_vars.o\
	amr.o\
	amr_criteria.o\
	amr_physics.o\
	calc_gradients.o\
	calc_sources.o\
	calc_timestep.o\
	coarsen_grid.o\
	covariant.o\
	create_soln_blocks.o\
	exchange_messages.o\
	explicit.o\
	find_neighbors.o\
	fix_axis_cells.o\
	library.o\
	load_balance.o\
	message_pass_cells.o\
	message_pass_faces.o\
	message_pass_nodes.o\
	number_blocks.o\
	octree.o\
	potential_force_averages.o\
	refine_grid.o\
	rotate.o\
	satellites.o\
	set_b0.o\
	set_BCs.o\
	set_block_geometry.o\
	set_ICs.o\
	set_outer_BCs.o\
	set_physics.o\
	specify_refinement.o\
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
	${OBJ1} ${OBJ2} ${OBJ3} ${OBJ3B} ${OBJ4} \
	${OBJ10} ${OBJ15}


# Default user module
ModUser.f90:
	cp ModUserDefault.f90 ModUser.f90

NAME = Advance

STATIC1:
	cp -f Mod${NAME}_static.f90 Mod${NAME}.f90

STATIC:
	make -s STATIC1 NAME=Advance
	make -s STATIC1 NAME=Nodes
	make -s STATIC1 NAME=Geometry
	make -s STATIC1 NAME=Implicit	#^CFG IF IMPLICIT
	make -s STATIC1 NAME=Raytrace	#^CFG IF RAYTRACE
	make -s STATIC1 NAME=CT		#^CFG IF CONSTRAINB

# Transform modules with static declaration of large variables into 
# modules with dynamic allocation as used modules
DYNAMIC1: 
	${SCRIPTDIR}/StaticToDynamic.pl ${LARGE} \
		Mod${NAME}_static.f90 > Mod${NAME}.f90

DYNAMIC:
	make -s DYNAMIC1 NAME=Advance
	make -s DYNAMIC1 NAME=Nodes
	make -s DYNAMIC1 NAME=Geometry
	@#^CFG IF IMPLICIT BEGIN
	make -s DYNAMIC1 NAME=Implicit \
		LARGE="-d=1 -l='MaxImplVar|,nBLK|,MaxImplBLK'" 
	@#^CFG END IMPLICIT
	make -s DYNAMIC1 NAME=Raytrace	#^CFG IF RAYTRACE
	make -s DYNAMIC1 NAME=CT	#^CFG IF CONSTRAINB

# Dependency rules for static/dynamic modules

# Default variable declaration type for 1 module (select one)
STATIC_OR_DYNAMIC = STATIC1
#STATIC_OR_DYNAMIC = DYNAMIC1

ModAdvance.f90: ModAdvance_static.f90
	make -s ${STATIC_OR_DYNAMIC} NAME=Advance

ModNodes.f90: ModNodes_static.f90
	make -s ${STATIC_OR_DYNAMIC} NAME=Nodes

ModGeometry.f90: ModGeometry_static.f90
	make -s ${STATIC_OR_DYNAMIC} NAME=Geometry

#^CFG IF CONSTRAINB BEGIN
ModCT.f90: ModCT_static.f90
	make -s ${STATIC_OR_DYNAMIC} NAME=CT
#^CFG END CONSTRAINB

#^CFG IF RAYTRACE BEGIN
ModRaytrace.f90: ModRaytrace_static.f90
	make -s ${STATIC_OR_DYNAMIC} NAME=Raytrace
#^CFG END RAYTRACE

#^CFG IF IMPLICIT BEGIN
ModImplicit.f90: ModImplicit_static.f90
	make -s ${STATIC_OR_DYNAMIC} NAME=Implicit \
		LARGE="-d=1 -l='MaxImplVar|,nBLK|,MaxImplBLK'"
#^CFG END IMPLICIT
############################################################################

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} \
		${MODULES} ${OBJECTS}

#
#	Making executables and libraries
#

#
# Library for framework and also for stand alone code
#

MY_LIB = libBATSRUS.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES} ${OBJECTS} ${MAKEFILE_COMP_SELECT}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES} ${OBJECTS}

# Object files not included into the library
OBJECTS_EXE = main.o

# Other requireed libraries
LIBSHARE  = ${LIBDIR}/libSHARE.a
LIBTIMING = ${LIBDIR}/libTIMING.a

# Libraries should be compiled first, because modules are used in main.
${OBJECTS_EXE}: ${LIBSHARE} ${LIBTIMING} ${MY_LIB}

EXE = ${BINDIR}/BATSRUS.exe

BATSRUS:
	make ${EXE}
	@echo ' '
	@echo Program BATSRUS has been brought up to date.
	@echo ' '	

${EXE}: ${OBJECTS_EXE}
	@(if [ ${OS} = "Darwin" ];         then make BATSRUS_for_Darwin; \
	else                                    make BATSRUS_general; fi)

BATSRUS_general:
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} \
		-L. -lBATSRUS -L${LIBDIR} -lMAGNETOGRAM -lEMPIRICALEE -lTIMING -lSHARE \
		${LBLAS} ${Lflag1}

BATSRUS_for_Darwin:
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ../${MY_LIB}; \
		ar -x ${LIBDIR}/libMAGNETOGRAM.a; \
		ar -x ${LIBDIR}/libEMPIRICALEE.a; \
		ar -x ${LIBTIMING}; \
		ar -x ${LIBSHARE}
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} Tmp_/*.o \
		${LBLAS} ${Lflag1}
	rm -rf Tmp_
#
#	cleaning
#

distclean: clean
	rm -f ModUser.f90.safe
	mv ModUser.f90 ModUser.f90.safe
	rm -f ModSize.f90 ModEquation.f90
	rm -f ModAdvance.f90 ModCT.f90 ModGeometry.f90 ModImplicit.f90 \
		ModNodes.f90 ModRaytrace.f90
