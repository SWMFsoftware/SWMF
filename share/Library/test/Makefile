#^CFG COPYRIGHT UM

include ../../../Makefile.def

SEARCH_EXTRA = -I../src

include ../../../Makefile.conf

NP = 2

help:
	@echo "You can test the following modules in ../src:"
	@echo
	@echo "   test_time_convert     - test ModTimeConvert"
	@echo "   test_coord_transform  - test ModTimeConvert"
	@echo "   test_freq             - test ModFreq"
	@echo "   test_sort             - test ModSort"
	@echo "   test_io_unit          - test ModIoUnit"
	@echo "   test_ray_trace        - test CON_ray_trace"
	@echo "   test_line_extract     - test CON_line_extract"
	@echo "   test_geopack          - test CON_geopack"
	@echo "   test_axes             - test CON_axes"
	@echo "   test_planet_field     - test CON_planet_field"
	@echo "   tests                 - test all the above"
	@echo
	@echo "To test everything with easy to read output try"
	@echo
	@echo "   make -s tests"
	@echo
	@echo "The number of PE-s for the parallel tests can be set by"
	@echo
	@echo "   make NP=4 test_line_extract"
	@echo
	@echo "The default value is 2."

tests: \
	test_time_convert \
	test_coord_transform \
	test_freq \
	test_sort \
	test_io_unit \
	test_ray_trace \
	test_line_extract \
	test_geopack \
	test_axes \
	test_planet_field

#######################################################################
OBJ_TIME_CONVERT = \
	test_time_convert.o \
	../src/ModTimeConvert.o

test_time_convert.exe: ${OBJ_TIME_CONVERT}
	${LINK.f90} -o test_time_convert.exe ${OBJ_TIME_CONVERT} ${Lflag2}

test_time_convert:
	@echo -----------------------------------------------------------------
	@echo test_time_convert starting
	@echo
	@(cd ../src; make ModTimeConvert.o)
	@make test_time_convert.exe
	@echo
	@echo './test_time_convert.exe'
	@./test_time_convert.exe
	@echo
	@echo test_time_convert done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_COORD_TRANSFORM = \
	test_coord_transform.o \
	../src/ModCoordTransform.o

test_coord_transform.exe: ${OBJ_COORD_TRANSFORM}
	${LINK.f90} -o test_coord_transform.exe ${OBJ_COORD_TRANSFORM} \
		${Lflag2}

test_coord_transform:
	@echo -----------------------------------------------------------------
	@echo test_coord_transform starting
	@echo
	@(cd ../src; make ModCoordTransform.o)
	@make test_coord_transform.exe
	@echo
	@echo './test_coord_transform.exe > test_coord_transform.tmp'
	@./test_coord_transform.exe > test_coord_transform.tmp
	@echo 'diff test_coord_transform.tmp test_coord_transform.log'
	-@(diff test_coord_transform.tmp test_coord_transform.log)
	@echo
	@echo test_coord_transform done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_FREQ = \
	test_freq.o \
	../src/ModFreq.o

test_freq.exe: ${OBJ_FREQ}
	${LINK.f90} -o test_freq.exe ${OBJ_FREQ} ${Lflag2}

test_freq:
	@echo -----------------------------------------------------------------
	@echo test_freq starting
	@echo
	@(cd ../src; make ModFreq.o)
	@make test_freq.exe
	@echo
	@echo './test_freq.exe'
	@./test_freq.exe
	@echo
	@echo test_freq done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_SORT = \
	test_sort.o \
	../src/ModSort.o

test_sort.exe: ${OBJ_SORT}
	${LINK.f90} -o test_sort.exe ${OBJ_SORT} ${Lflag2}

test_sort:
	@echo -----------------------------------------------------------------
	@echo test_sort starting
	@echo
	@(cd ../src; make ModSort.o)
	@make test_sort.exe
	@echo
	@echo './test_sort.exe'
	@./test_sort.exe
	@echo
	@echo test_sort done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_IO_UNIT = \
	test_io_unit.o \
	../src/ModIoUnit.o

test_io_unit.exe: ${OBJ_IO_UNIT}
	${LINK.f90} -o test_io_unit.exe ${OBJ_IO_UNIT} ${Lflag2}

test_io_unit:
	@echo -----------------------------------------------------------------
	@echo test_io_unit starting
	@echo
	@(cd ../src; make ModIoUnit.o)
	@make test_io_unit.exe
	@echo
	@echo './test_io_unit.exe'
	@./test_io_unit.exe
	@echo
	@echo test_io_unit done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_RAY_TRACE = \
	test_ray_trace.o \
	../src/ModMpi.o \
	../src/CON_ray_trace.o

test_ray_trace.exe: ${OBJ_RAY_TRACE}
	${LINK.f90} -o test_ray_trace.exe ${OBJ_RAY_TRACE} ${Lflag1}

test_ray_trace:
	@echo -----------------------------------------------------------------
	@echo test_ray_trace starting
	@echo
	@(cd ../src; make ModMpi.o CON_ray_trace.o)
	@make test_ray_trace.exe
	@echo
	@echo 'mpirun -np ${NP} test_ray_trace.exe > test_ray_trace.tmp'
	@mpirun -np ${NP} test_ray_trace.exe > test_ray_trace.tmp
	@echo 'diff test_ray_trace.tmp test_ray_trace.log'
	-@(diff test_ray_trace.tmp test_ray_trace.log)
	@echo
	@echo test_ray_trace done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_LINE_EXTRACT = \
	test_line_extract.o \
	../src/ModMpi.o \
	../src/CON_line_extract.o

test_line_extract.exe: ${OBJ_LINE_EXTRACT}
	${LINK.f90} -o test_line_extract.exe ${OBJ_LINE_EXTRACT} ${Lflag1}

test_line_extract:
	@echo -----------------------------------------------------------------
	@echo test_line_extract starting
	@echo
	@(cd ../src; make ModMpi.o CON_line_extract.o)
	@make test_line_extract.exe
	@echo
	@echo 'mpirun -np ${NP} test_line_extract.exe'
	@mpirun -np ${NP} test_line_extract.exe 
	@echo
	@echo test_line_extract done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_GEOPACK = \
	test_geopack.o \
	../src/CON_geopack.o \
	../src/ModCoordTransform.o

test_geopack.exe: ${OBJ_GEOPACK}
	${LINK.f90} -o test_geopack.exe ${OBJ_GEOPACK} ${Lflag2}

test_geopack:
	@echo -----------------------------------------------------------------
	@echo test_geopack starting
	@echo
	@(cd ../src; make ModCoordTransform.o CON_geopack.o)
	@make test_geopack.exe
	@echo
	@echo './test_geopack.exe > test_geopack.tmp'
	@./test_geopack.exe > test_geopack.tmp
	@echo 'diff test_geopack.tmp test_geopack.log'
	-@(diff test_geopack.tmp test_geopack.log)
	@echo
	@echo test_geopack done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_AXES = \
	test_axes.o \
	../src/CON_axes.o \
	../src/CON_planet.o \
	../src/ModCoordTransform.o \
	../src/ModIoUnit.o \
	../src/ModReadParam.o \
	../src/ModTimeConvert.o \
	../src/ModUtilities.o

test_axes.exe: ${OBJ_AXES}
	${LINK.f90} -o test_axes.exe ${OBJ_AXES} ${Lflag1}

test_axes: 
	@echo -----------------------------------------------------------------
	@echo test_axes starting
	@echo
	@(cd ../src; make CON_axes.o)
	@make test_axes.exe
	@echo
	@echo './test_axes.exe'
	@./test_axes.exe
	@echo
	@echo test_axes done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_PLANET_FIELD = \
	test_planet_field.o \
	../src/CON_planet_field.o \
	../src/CON_planet.o \
	../src/CON_axes.o \
	../src/ModCoordTransform.o \
	../src/ModIoUnit.o \
	../src/ModReadParam.o \
	../src/ModTimeConvert.o \
	../src/ModUtilities.o

test_planet_field.exe: ${OBJ_PLANET_FIELD}
	${LINK.f90} -o test_planet_field.exe ${OBJ_PLANET_FIELD} ${Lflag1}

test_planet_field: 
	@echo -----------------------------------------------------------------
	@echo test_planet_field starting
	@echo
	@(cd ../src; make CON_planet_field.o)
	@make test_planet_field.exe
	@echo
	@echo './test_planet_field.exe > test_planet_field.tmp'
	@./test_planet_field.exe > test_planet_field.tmp
	@echo 'diff test_planet_field.tmp test_planet_field.log'
	-@(diff test_planet_field.tmp test_planet_field.log)
	@echo
	@echo test_planet_field done
	@echo -----------------------------------------------------------------
#######################################################################

distclean: clean
	rm -f *.tmp

# keep this line

