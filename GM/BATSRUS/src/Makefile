#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND
include Makefile.RULES

# Temporary rules to convert from user routines to user modules and back
SRCUSER = amr_criteria.f90 BATS_methods.f90 calc_sources_covar.f90 \
	calc_sources.f90 create_soln_blocks.f90 MH_set_parameters.f90 \
	set_b0.f90 set_BCs.f90 set_block_geometry.f90 set_ICs.f90 \
	specify_refinement.f90 set_outer_BCs.f90 specify_refinement.f90 \
	update_states.f90 write_logfile.f90 write_progress.f90

USEUSER = use ModUser, ONLY:

USEUSERMODULE:
	perl -pi -e 's/ \!\!${USEUSER}/ ${USEUSER}/' ${SRCUSER}
	perl -pi -e 's/^OBJ13 = user/##OBJ13 = user/' Makefile
	perl -pi -e 's/^##OBJ13 = ModUser/OBJ13 = ModUser/' Makefile

USEUSERFILE:
	perl -pi -e 's/ ${USEUSER}/ \!\!${USEUSER}/' ${SRCUSER}
	perl -pi -e 's/^##OBJ13 = user/OBJ13 = user/' Makefile
	perl -pi -e 's/^OBJ13 = ModUser/##OBJ13 = ModUser/' Makefile

# end of temporary rules

# CONFIGURABLE MODULES

MOD1 = ModCT.o       #^CFG IF CONSTRAINB
MOD2 = ModImplicit.o #^CFG IF IMPLICIT
MOD3 = ModRaytrace.o #^CFG IF RAYTRACE
MOD4 = ModProject.o  #^CFG IF PROJECTION

#^CFG IF NOT MULTISP BEGIN
MODFLCOV=ModFlux_covar.o
#^CFG END MULTISP
#MODFLCOV=ModFlux_covar_mars.o         #^CFG UNCOMMENT IF MULTISP
#MODFLCOV=ModFlux_modif_covar.o        #^CFG UNCOMMENT IF ALWAVES

MOD7 = ModCovariant.o ${MODFLCOV}   #^CFG IF NOT CARTESIAN

#^CFG IF NOT GLOBALHELIOSPHERE BEGIN
#^CFG IF NOT MULTISP BEGIN
#^CFG IF NOT ALWAVES BEGIN
MODEL = ModVarIndexes.o 
#^CFG END ALWAVES
#^CFG END MULTISP
#^CFG END GLOBALHELIOSPHERE
#MODEL = ModVarIndexes_modif.o    #^CFG UNCOMMENT IF ALWAVES 
#MODEL = ModVarIndexes_mars.o     #^CFG UNCOMMENT IF MULTISP
#MODEL = ModVarIndexes_globhelio.o #^CFG UNCOMMENT IF GLOBALHELIOSPHERE

MODULES = \
	ModBlockData.o\
	ModInterface.o\
	ModProcMH.o\
	ModSize.o\
	ModMain.o\
	ModAdvance.o\
	ModGeometry.o\
	ModIO.o\
	ModPhysics.o\
	ModAMR.o\
	ModPartSteady.o\
	ModParallel.o\
	ModOctree.o\
	ModNodes.o\
	${MOD1}	${MOD2}	${MOD3} ${MOD4} ${MOD7} ${MODEL}

# OBJECT FILES FOR MAIN EXECUTABLE

# CONFIGURABLE OBJECTS				^CFG IF CONFIGURE

OBJ1=constrain_B.o #				^CFG IF CONSTRAINB

#^CFG IF IMPLICIT BEGIN
OBJ2=   bicgstab.o\
	gmres.o\
	implicit.o\
	impl_interface.o\
	impl_jacobian.o\
	impl_matvec.o\
	impl_newton.o\
	ModPrecond.o precond.o ${BLAS}
#^CFG END IMPLICIT

OBJ3=ray_trace_new.o ray_trace.o ray_pass.o #	^CFG IF RAYTRACE
OBJ3B=get_im_pressure.o #			^CFG IF RCM
OBJ4=project_B.o proj_cg.o proj_bicgstab.o #	^CFG IF PROJECTION

#^CFG IF RUSANOVFLUX BEGIN
RUSCART=calc_flux_Rusanov.o                 #^CFG IF CARTESIAN
RUSCOVAR=calc_flux_rusanov_covar.o          #^CFG IF NOT CARTESIAN
OBJ5=${RUSCART} ${RUSCOVAR}
#^CFG END RUSANOVFLUX

#^CFG IF LINDEFLUX BEGIN  
                                            #^CFG IF CARTESIAN BEGIN
LINDECART=calc_flux_Linde.o                 #^CFG IF NOT GLOBALHELIOSPHERE
#LINDECART=calc_flux_Linde_globhelio.o      #^CFG UNCOMMENT IF GLOBALHELIOSPHERE
                                            #^CFG END CARTESIAN
                                            #^CFG IF NOT CARTESIAN BEGIN
#LINDECOVAR=calc_flux_linde_covar_mars.o    #^CFG UNCOMMENT IF MULTISP
LINDECOVAR=calc_flux_linde_covar.o          #^CFG IF NOT MULTISP
                                            #^CFG END CARTESIAN 
OBJ6=${LINDECOVAR} ${LINDECART}     
#^CFG END LINDEFLUX                  

#^CFG IF AWFLUX BEGIN
                                            #^CFG IF NOT CARTESIAN BEGIN
AWCOVAR=calc_flux_aw_covar.o                #^CFG IF NOT ALWAVES
#AWCOVAR=calc_flux_aw_modif_covar.o         #^CFG UNCOMMENT IF ALWAVES
                                            #^CFG END CARTESIAN

                                            #^CFG IF CARTESIAN BEGIN
#AWCART=calc_flux_aw_mars.o                 #^CFG UNCOMMENT IF MULTISP
#AWCART=calc_flux_AW_modif.o                #^CFG UNCOMMENT IF ALWAVES
                                            #^CFG IF NOT MULTISP BEGIN
AWCART = calc_flux_AW.o                     #^CFG IF NOT ALWAVES
                                            #^CFG END MULTISP       
                                            #^CFG END CARTESIAN
OBJ7=${AWCOVAR} ${AWCART}   
#^CFG END AWFLUX

#^CFG IF ROEFLUX BEGIN
OBJ8=calc_flux_Roe.o                        #^CFG IF CARTESIAN
#^CFG END ROEFLUX
 
CARTSOURCES=calc_sources.o                  #^CFG IF CARTESIAN
COVARCOURCES=calc_sources_covar.o           #^CFG IF NOT CARTESIAN
OBJ9=${CARTSOURCES} ${COVARCOURCES}

#^CFG IF DISSFLUX BEGIN
OBJ10=\
	calc_heat_flux.o\
	calc_resistive_flux.o\
	set_ICs_diss_test.o
#^CFG END DISSFLUX

OBJ11=set_outer_BCs.o #				^CFG IF CELLOUTERBC

OBJ12= write_plot_los.o write_plot_sph.o \
		     set_ICs_arc.o set_ICs_cme.o      

#^CFG IF USERFILES BEGIN
OBJ13 = user_routines.o 
##OBJ13 = ModUserEmpty.o ModUser.o

user_routines.f90:
	cp user_routines_orig.f90 user_routines.f90

ModUser.f90:
	cp ModUserDefault.f90 ModUser.f90

#^CFG END USERFILES


OBJ14=covariant.o                           #^CFG IF NOT CARTESIAN



OBJECTS = \
	BATS_methods.o \
	MH_set_parameters.o\
	advect_points.o\
	allocate_vars.o\
	amr.o\
	amr_criteria.o\
	amr_physics.o\
	calc_facefluxes.o\
	calc_facevalues_opt.o\
	calc_gradients.o\
	calc_timestep.o\
	coarsen_grid.o\
	conserve_flux.o\
	create_soln_blocks.o\
	exchange_messages.o\
	explicit.o\
	find_neighbors.o\
	get_solar_wind_point.o\
	index_var.o\
	library.o\
	load_balance.o\
	message_pass_dir.o\
	message_pass_cells.o\
	message_pass_faces.o\
	message_pass_nodes.o\
	number_blocks.o\
	octree.o\
	potential_force_averages.o\
	refine_grid.o\
	restart.o\
	rotate.o\
	satellites.o\
	set_b0.o\
	set_BCs.o\
	set_block_geometry.o\
	set_ICs.o\
	set_physics.o\
	specify_refinement.o\
	update_states.o\
	update_states_MHD.o\
	write_logfile.o\
	write_plot_common.o\
	write_plot_idl.o\
	write_plot_tec.o\
	write_progress.o\
	${OBJ1} ${OBJ2} ${OBJ3} ${OBJ3B} ${OBJ4} \
	${OBJ5} ${OBJ6} ${OBJ7} ${OBJ8} ${OBJ9}\
	${OBJ10} ${OBJ11} ${OBJ12} ${OBJ13} ${OBJ14}

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

#
# Stand alone executable
#

# Object files not included into the library
OBJECTS_EXE = main.o stand_alone.o

# Other requireed libraries
LIBSHARE  = ${LIBDIR}/libSHARE.a
LIBTIMING = ${LIBDIR}/libTIMING.a

# Libraries should be compiled first, because modules are used in main.
${OBJECTS_EXE}: ${LIBSHARE} ${LIBTIMING} ${MY_LIB}

BATSRUS.exe: ${OBJECTS_EXE}
	@(if [ ${COMPILE.f90} = "xlf90" ]; then make BATSRUS_for_xlf; \
	elif [ ${OS} = "Darwin" ];         then make BATSRUS_for_Darwin; \
	else                                    make BATSRUS_general; fi)

BATSRUS_general:
	${LINK.f90} -o BATSRUS.exe ${OBJECTS_EXE} \
		-L. -lBATSRUS -L${LIBDIR} -lTIMING -lSHARE \
		${Lflag1}

BATSRUS_for_Darwin:
	ld -all_load -r -o LIBRARIES.o \
		${LIBSHARE} ${LIBTIMING} ${MY_LIB}
	${LINK.f90} -o BATSRUS.exe ${OBJECTS_EXE} \
		LIBRARIES.o ${LBLAS} ${Lflag1}
	rm -f LIBRARIES.o

BATSRUS_for_xlf:
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; ar -x ../${MY_LIB}; ar -x ${LIBTIMING}; ar -x ${LIBSHARE}
	${LINK.f90} -o BATSRUS.exe ${OBJECTS_EXE} Tmp_/*.o \
		${LBLAS} ${Lflag1}
	rm -rf Tmp_
#
#	cleaning
#

distclean: clean
	rm -f Makefile.RULES Makefile.DEPEND ModSize.f90 ModUser.f90
	rm -f ModAdvance.f90 ModCT.f90 ModGeometry.f90 ModImplicit.f90 \
		ModNodes.f90 ModRaytrace.f90
	rm -f stand_alone.f90

# keep this line
