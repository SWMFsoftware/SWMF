#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf

help:
	@echo "You can test the following modules in ../src:"
	@echo
	@echo "   test_time_convert  - test ModTimeConvert"
	@echo "   test_ray_trace     - test CON_ray_trace"
	@echo "   test_planet_field  - test CON_planet_field"
	@echo "   tests              - test all the above"
	@echo
	@echo "To test everything with easy to read output try"
	@echo
	@echo "   make -s tests"
	@echo

tests: \
	test_time_convert \
	test_ray_trace \
	test_planet_field

#######################################################################
OBJ_RAY_TEST = \
	test_ray_trace.o \
	../src/ModMpi.o \
	../src/CON_ray_trace.o

test_ray_trace.exe: ${OBJ_RAY_TEST}
	cd ../src; make ModMpi.o CON_ray_trace.o
	${LINK.f90} ${Lflag1} -o test_ray_trace.exe ${OBJ_RAY_TEST}

test_ray_trace: test_ray_trace.exe
	@echo -----------------------------------------------------------------
	@echo test_ray_trace starting
	@echo
	mpirun -np 2 test_ray_trace.exe > test_ray_trace.tmp
	@echo 'test_ray_trace: there should be no difference:'
	diff test_ray_trace.tmp test_ray_trace.log
	@echo
	@echo test_ray_trace done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_TIME_TEST = \
	test_time_convert.o \
	../src/ModTimeConvert.o

test_time_convert.exe: ${OBJ_TIME_TEST}
	cd ../src; make ModTimeConvert.o
	${LINK.f90} ${Lflag2} -o test_time_convert.exe ${OBJ_TIME_TEST}

test_time_convert: test_time_convert.exe
	@echo -----------------------------------------------------------------
	@echo test_time_convert starting
	@echo
	@echo 'test_time_convert: there should be no output'
	./test_time_convert.exe
	@echo
	@echo test_time_convert done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_FIELD_TEST = \
	test_planet_field.o \
	../src/ModCoordTransform.o \
	../src/ModIoUnit.o \
	../src/ModReadParam.o \
	../src/ModTimeConvert.o \
	../src/ModUtilities.o \
	../src/CON_axes.o \
	../src/CON_planet.o \
	../src/CON_planet_field.o

test_planet_field.exe: ${OBJ_FIELD_TEST}
	${LINK.f90} ${Lflag2} -o test_planet_field.exe ${OBJ_FIELD_TEST}

test_planet_field: 
	@echo -----------------------------------------------------------------
	@echo test_planet_field starting
	@echo
	@(cd ../src; make LIB)
	@make test_planet_field.exe
	./test_planet_field.exe > test_planet_field.tmp
	@echo 'test_planet_field: there should be no difference:'
	diff test_planet_field.tmp test_planet_field.log
	@echo
	@echo test_planet_field done
	@echo -----------------------------------------------------------------
#######################################################################

distclean: clean
	rm -f *.tmp

# keep this line

