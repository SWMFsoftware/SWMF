#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
include Makefile.OPTIONS
include Makefile.DEPEND
include Makefile.RULES

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

#MOD7 = ModCovariant.o ${MODFLCOV}   #^CFG IF NOT CARTESIAN

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
	${IH_WRAPPER_MOD}\
	ModCompatibility.o\
	ModInterface.o\
	ModProcMH.o\
	ModSize.o\
	ModMain.o\
	ModAdvance.o\
	ModGeometry.o\
	ModIO.o\
	ModPhysics.o\
	ModAMR.o\
	ModParallel.o\
	ModOctree.o\
	ModNodes.o\
	${MOD1}	${MOD2}	${MOD3} ${MOD4} ${MOD7} ${MODEL}

# OBJECT FILES FOR MAIN EXECUTABLE

# CONFIGURABLE OBJECTS				^CFG IF CONFIGURE

OBJ1=constrain_B.o #				^CFG IF CONSTRAINB

#^CFG IF IMPLICIT BEGIN
OBJ2A=bicgstab.o #				^CFG IF NOT SIMPLE
OBJ2=   ${OBJ2A}\
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
OBJ5=calc_flux_Rusanov.o                    #^CFG IF CARTESIAN
#OBJ5=calc_flux_rusanov_covar.o             #^CFG IF NOT CARTESIAN
#^CFG END RUSANOVFLUX

#^CFG IF LINDEFLUX BEGIN  
LINDECART=calc_flux_Linde.o                 #^CFG IF NOT GLOBALHELIOSPHERE
#LINDECART=calc_flux_Linde_globhelio.o      #^CFG UNCOMMENT IF GLOBALHELIOSPHERE    
OBJ6=${LINDECART}                           #^CFG IF CARTESIAN

#LINDECOVAR=calc_flux_linde_covar_mars.o    #^CFG UNCOMMENT IF MULTISP
LINDECOVAR=calc_flux_linde_covar.o          #^CFG IF NOT MULTISP
#OBJ6=${LINDECOVAR}                         #^CFG IF NOT CARTESIAN 
#^CFG END LINDEFLUX                  

#^CFG IF AWFLUX BEGIN
#OBJ7=calc_flux_aw_covar.o                  #^CFG IF NOT CARTESIAN
#OBJ7=calc_flux_aw_modif_covar.o            #^CFG UNCOMMENT IF ALWAVES
#^CFG IF CARTESIAN BEGIN
#OBJ7=calc_flux_aw_mars.o          #^CFG UNCOMMENT IF MULTISP
#OBJ7=calc_flux_AW_modif.o         #^CFG UNCOMMENT IF ALWAVES
#^CFG IF NOT MULTISP BEGIN
OBJ7= calc_flux_AW.o                 	    #^CFG IF NOT ALWAVES
#^CFG END MULTISP       
#^CFG END CARTESIAN
#^CFG END AWFLUX

#^CFG IF ROEFLUX BEGIN
OBJ8=calc_flux_Roe.o                        #^CFG IF CARTESIAN
#^CFG END ROEFLUX
 
OBJ9=calc_sources.o                        #^CFG IF CARTESIAN
#OBJ9=calc_sources_covar.o covariant.o     #^CFG IF NOT CARTESIAN

#^CFG IF DISSFLUX BEGIN
OBJ10=\
	calc_heat_flux.o\
	calc_resistive_flux.o\
	set_ICs_diss_test.o
#^CFG END DISSFLUX

OBJ11=set_outer_BCs.o #				^CFG IF CELLOUTERBC

#^CFG IF NOT SIMPLE BEGIN
OBJ12= write_plot_los.o write_plot_sph.o \
		     set_ICs_arc.o set_ICs_cme.o      
#^CFG END SIMPLE

#^CFG IF USERFILES BEGIN
OBJ13 = user_routines.o 
user_routines.f90:
	cp user_routines_orig.f90 user_routines.f90
#^CFG END USERFILES

OBJECTS = \
	${IH_WRAPPER_OBJ}\
	BATS_methods.o \
	MH_set_parameters.o\
	conserve_flux.o\
	allocate_vars.o\
	amr.o\
	amr_criteria.o\
	amr_physics.o\
	calc_facefluxes.o\
	calc_facevalues_opt.o\
	calc_gradients.o\
	calc_timestep.o\
	coarsen_grid.o\
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
	perfctr_cray.o\
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
	${OPTION_FORCE}\
	${OBJ1} ${OBJ2} ${OBJ3} ${OBJ3B} ${OBJ4} \
	${OBJ5} ${OBJ6} ${OBJ7} ${OBJ8} ${OBJ9}\
	${OBJ10} ${OBJ11} ${OBJ12} ${OBJ13}

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

${MY_LIB}: ${MODULES} ${OBJECTS} \
	Makefile.OPTIONS ${MAKEFILE_COMP_SELECT}
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
	@(if [ ${OS} != "Darwin" ]; then make BATSRUS_general; \
	else                             make BATSRUS_Darwin; fi)

BATSRUS_general:
	${LINK.f90} -o BATSRUS.exe ${OBJECTS_EXE} \
		-L. -lBATSRUS -L${LIBDIR} -lTIMING -lSHARE \
		${Lflag1}

BATSRUS_Darwin:
	ld -all_load -r -o LIBRARIES.o \
		${LIBSHARE} ${LIBTIMING} ${MY_LIB}
	${LINK.f90} -o BATSRUS.exe ${OBJECTS_EXE} \
		LIBRARIES.o ${LBLAS} \
		${Lflag1}
	rm -f LIBRARIES.o

#
#	cleaning
#

distclean: clean
	rm -f Makefile.RULES Makefile.DEPEND Makefile.OPTIONS ModSize.f90
	rm -f stand_alone.f90

# keep this line
