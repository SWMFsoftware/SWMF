include ../../Makefile.def

SEARCH_EXTRA = -I../ -I${DATAREADINDICESDIR} -I${EMPIRICALIEDIR} -I${EMPIRICALUADIR} \
	-I${LIBDIR}

include ../../Makefile.conf

MPIRUN = mpirun -np 2

help:
	@echo "Tests for Particle module: unit tests and combinations"
	@echo
	@echo "   test_axes             - test CON_axes"
	@echo "   test_coord_transform  - test ModTimeConvert"
	@echo "   test                  - test all the above"
	@echo
	@echo "To test everything with easy to read output try"
	@echo
	@echo "   make -s test"
	@echo
	@echo "The parallel tests can be executed like this"
	@echo
	@echo "   make MPIRUN='mpirun -np 4' test_line_extract"
	@echo
	@echo "The default value for MPIRUN is mpirun -np 2."

test:   unit_test_sample
	@echo ---------------------------------------------------------
	@ls -l *.diff

distclean: clean
	rm -rf *.tmp *.diff test*.out test_initial_state.ref \
		test_process_var_name.log

#######################################################################
OBJ_SAMPLE = \
	../ModParticle.o 


unit_test_sample_exe:
	@(cd .. ;make ModParticle.o)
	@make unit_test_sample.o 
	${LINK.f90} -o unit_test_sample.exe ${OBJ_SAMPLE} \
		unit_test_sample.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_sample: 
	@echo unit_test_sample starting
	@make unit_test_sample_exe 
	@echo ./unit_test_sample_exe
	-@(./unit_test_sample_exe )


#######################################################################
OBJ_SPLIT_JOIN = \
	../ModParticle.o 


unit_test_split_join_exe:
	@(cd .. ;make ModParticle.o)
	@make unit_test_split_join.o 
	${LINK.f90} -o unit_test_split_join.exe ${OBJ_SPLIT_JOIN} \
		unit_test_split_join.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_split_join: 
	@echo unit_test_split_join starting
	@make unit_test_split_join_exe 
	@echo ./unit_test_split_join_exe
	-@(./unit_test_split_join_exe )

#######################################################################

OBJ_PUSHER = \
	../ModParticle.o 


unit_test_pusher_exe:
	@(cd .. ;make ModParticle.o)
	@make unit_test_pusher.o 
	${LINK.f90} -o unit_test_pusher.exe ${OBJ_PUSHER} \
		unit_test_pusher.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_pusher: 
	@echo unit_test_pusher starting
	@make unit_test_pusher_exe 
	@echo ./unit_test_pusher_exe
	-@(./unit_test_pusher_exe )

#######################################################################

#######################################################################

OBJ_COULOMB_COLLISION = \
	../ModParticle.o 


unit_test_coulomb_collision_exe:
	@(cd .. ;make ModParticle.o)
	@make unit_test_coulomb_collision.o 
	${LINK.f90} -o unit_test_coulomb_collision.exe ${OBJ_COULOMB_COLLISION} \
		unit_test_coulomb_collision.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_coulomb_collision: 
	@echo unit_test_coulomb_collision starting
	@make unit_test_coulomb_collision_exe 
	@echo ./unit_test_coulomb_collision_exe
	-@(./unit_test_coulomb_collision_exe )

#######################################################################

#######################################################################

OBJ_WPI = \
	../ModParticle.o 


unit_test_wpi_exe:
	@(cd .. ;make ModParticle.o)
	@make unit_test_wpi.o 
	${LINK.f90} -o unit_test_wpi.exe ${OBJ_WPI} \
		unit_test_wpi.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_wpi: 
	@echo unit_test_wpi starting
	@make unit_test_wpi_exe 
	@echo ./unit_test_wpi_exe
	-@(./unit_test_wpi_exe )

#######################################################################

#######################################################################
OBJ_COMBINE_FLUID_PARTICLE = \
	../ModParticle.o 


unit_test_combine_fluid_particle_exe:
	@(cd .. ;make ModParticle.o)
	@make unit_test_combine_fluid_particle.o 
	${LINK.f90} -o unit_test_combine_fluid_particle.exe ${OBJ_COMBINE_FLUID_PARTICLE} \
		unit_test_combine_fluid_particle.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_combine_fluid_particle: 
	@echo unit_test_combine_fluid_particle starting
	@make unit_test_combine_fluid_particle_exe 
	@echo ./unit_test_combine_fluid_particle_exe
	-@(./unit_test_combine_fluid_particle_exe )


