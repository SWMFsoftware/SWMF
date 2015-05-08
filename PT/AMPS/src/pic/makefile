SHELL =/bin/sh

DEFAULT_TARGET : amps  

BATL=nobatl

include ../../Makefile.conf
include ../../Makefile.local

PICOBJ = pic.o \
	pic_parallel.o \
	pic_bc.o \
	pic_init_const.o \
	pic_mesh.o \
	pic_mol.o \
	pic_mol_parser.o \
	pic_mover.o \
	pic_mover_boris.o \
	pic_parser.o \
	pic_pbuffer.o \
	pic_weight_time.o \
	pic_internal_boundary_sphere.o \
	pic_sample_distribution_function.o \
	pic_sample_flux_distribution.o \
	pic_vel_distribution.o \
	pic_ices.o \
	pic_photolytic_reactions.o \
	pic_volume_particle_injection.o \
	pic_background_atmosphere.o \
	pic_column_integration.o \
	pic_particle_collision_models.o \
	pic_run_time_system_state.o \
	pic_lb.o \
	pic_qlb.o \
	pic_debugger.o \
	pic__model__electrically_charged_dust.o \
	pic_internal_boundary_rotation_body.o \
	pic_internal_boundary_nastran_surface.o \
	pic_swmf.o \
	pic__ray_tracing.o \
	pic_cut_cells.o \
	pic_pt.o \
	pic_sample_pitch_angle_distribution.o \
	pic_kameleon.o \
	pic_restart.o \
	pic_arms.o \
	pic_tecplot.o \
	pic_initial_conditions.o \
	pic_datafile.o
	
ifneq ($(BATL),nobatl)	
PICOBJ+=pic_batsrus.o \
	pic_batsrus_wrapper.o
endif

amps : $(PICOBJ)
	ar -src amps.a $(PICOBJ)
	@echo "Target 'amps' is done." 
