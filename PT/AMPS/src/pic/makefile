SHELL =/bin/sh

DEFAULT_TARGET : amps  

clean:
	rm -f *.o *.a

CC=mpicxx 
Flags=-g 
IncludeList=-I../general -I../meshAMR 

.SUFFIXES: .cpp .o 
.cpp.o: 
	${CC} -c ${Flags} ${IncludeList} $< 


PICOBJ = pic.o \
	pic_parallel.o \
	pic_bc.o \
	pic_init_const.o \
	pic_mesh.o \
	pic_mol.o \
	pic_mol_parser.o \
	pic_mover.o \
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
	pic_run_time_system_state.o 
amps : $(PICOBJ)
	ar -src amps.a $(PICOBJ)
	@echo "Target 'amps' is done." 
