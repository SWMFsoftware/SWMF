#^CFG COPYRIGHT UM

include ../../../Makefile.def

SEARCH_EXTRA = -I../src

include ../../../Makefile.conf

MPIRUN = mpirun -np 2

help:
	@echo "You can test the following modules in ../src:"
	@echo
	@echo "   test_time_convert     - test ModTimeConvert"
	@echo "   test_coord_transform  - test ModTimeConvert"
	@echo "   test_freq             - test ModFreq"
	@echo "   test_linear_solver    - test ModLinearSolver"
	@echo "   test_sort             - test ModSort"
	@echo "   test_io_unit          - test ModIoUnit"
	@echo "   test_interpolate      - test ModInterpolate"
	@echo "   test_lookup_table     - test ModLookupTable"
	@echo "   test_plot_file        - test ModPlotFile"
	@echo "   test_ray_trace        - test CON_ray_trace"
	@echo "   test_line_extract     - test CON_line_extract"
	@echo "   test_geopack          - test CON_geopack"
	@echo "   test_axes             - test CON_axes"
	@echo "   test_planet_field     - test CON_planet_field"
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

test: \
	test_time_convert \
	test_coord_transform \
	test_freq \
	test_sort \
	test_io_unit \
	test_interpolate \
	test_lookup_table \
	test_ray_trace \
	test_line_extract \
	test_geopack \
	test_axes \
	test_planet_field \
	test_plot_file \
	test_linear_advection

#######################################################################
OBJ_TIME_CONVERT = \
	test_time_convert.o \
	../src/ModTimeConvert.o

test_time_convert_exe:
	@(cd ../src; make ModTimeConvert.o)
	@make test_time_convert.o
	${LINK.f90} -o test_time_convert.exe ${OBJ_TIME_CONVERT} ${Lflag2}

test_time_convert:
	@echo -----------------------------------------------------------------
	@echo test_time_convert starting
	@echo
	@make test_time_convert_exe
	@echo
	@echo ./test_time_convert.exe
	-@(./test_time_convert.exe)
	@echo
	@echo test_time_convert done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_COORD_TRANSFORM = \
	test_coord_transform.o \
	../src/ModCoordTransform.o

test_coord_transform_exe:
	@(cd ../src; make ModCoordTransform.o)
	@make test_coord_transform.o
	${LINK.f90} -o test_coord_transform.exe ${OBJ_COORD_TRANSFORM} \
		${Lflag2}

test_coord_transform: 
	@echo -----------------------------------------------------------------
	@echo test_coord_transform starting
	@echo
	@make test_coord_transform_exe
	@echo
	@echo './test_coord_transform.exe > test_coord_transform.tmp'
	-@(./test_coord_transform.exe > test_coord_transform.tmp)
	@echo 'DiffNum.pl test_coord_transform.tmp test_coord_transform.log'
	-@(${SCRIPTDIR}/DiffNum.pl -b -a=1e-12 \
		test_coord_transform.tmp test_coord_transform.log)
	@echo
	@echo test_coord_transform done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_FREQ = \
	test_freq.o \
	../src/ModFreq.o

test_freq_exe:
	@(cd ../src; make ModFreq.o)
	@make test_freq.o
	${LINK.f90} -o test_freq.exe ${OBJ_FREQ} ${Lflag2}

test_freq:
	@echo -----------------------------------------------------------------
	@echo test_freq starting
	@echo
	@make test_freq_exe
	@echo
	@echo './test_freq.exe'
	-@(./test_freq.exe)
	@echo
	@echo test_freq done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_LINEAR_SOLVER = \
	test_linear_solver.o \
	../src/blas.o \
	../src/lapack.o \
	../src/ModMpi.o \
	../src/ModMpiInterfaces.o \
	../src/ModBlasLapack.o \
	../src/ModLinearSolver.o

test_linear_solver_exe:
	@(cd ../src; make blas.o lapack.o ModLinearSolver.o)
	@make test_linear_solver.o
	${LINK.f90} -o test_linear_solver.exe ${OBJ_LINEAR_SOLVER} ${Lflag1}

test_linear_solver:
	@echo -----------------------------------------------------------------
	@echo test_linear_solver starting
	@echo
	@make test_linear_solver_exe
	@echo
	@echo './test_linear_solver.exe'
	-@(./test_linear_solver.exe)
	@echo
	@echo test_linear_solver done
	@echo -----------------------------------------------------------------
#######################################################################
OBJ_LINEAR_ADVECTION = \
	test_linear_advection.o \
	../src/ModLinearAdvection.o 

test_linear_advection_exe:
	@(cd ../src; make ModLinearAdvection.o)
	@make test_linear_advection.o
	${LINK.f90} -o test_linear_advection.exe ${OBJ_LINEAR_ADVECTION} ${Lflag1}

test_linear_advection:
	@echo -----------------------------------------------------------------
	@echo test_linear_advection starting
	@echo
	@make test_linear_advection_exe
	@echo
	@echo './test_linear_advection.exe'
	-@(./test_linear_advection.exe)
	@echo
	@echo test_linear_advection done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_SORT = \
	test_sort.o \
	../src/ModSort.o

test_sort_exe:
	@(cd ../src; make ModSort.o)
	@make test_sort.o
	${LINK.f90} -o test_sort.exe ${OBJ_SORT} ${Lflag2}

test_sort:
	@echo -----------------------------------------------------------------
	@echo test_sort starting
	@echo
	@make test_sort_exe
	@echo
	@echo './test_sort.exe'
	-@(./test_sort.exe)
	@echo
	@echo test_sort done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_IO_UNIT = \
	test_io_unit.o \
	../src/ModIoUnit.o

test_io_unit_exe:
	@(cd ../src; make ModIoUnit.o)
	@make test_io_unit.o
	${LINK.f90} -o test_io_unit.exe ${OBJ_IO_UNIT} ${Lflag2}

test_io_unit:
	@echo -----------------------------------------------------------------
	@echo test_io_unit starting
	@echo
	@make test_io_unit_exe
	@echo
	@echo ./test_io_unit.exe
	-@(./test_io_unit.exe)
	@echo
	@echo test_io_unit done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_INTERPOLATE = \
	test_interpolate.o \
	../src/ModInterpolate.o

test_interpolate_exe:
	@(cd ../src; make ModInterpolate.o)
	@make test_interpolate.o
	${LINK.f90} -o test_interpolate.exe ${OBJ_INTERPOLATE} ${Lflag2}

test_interpolate:
	@echo -----------------------------------------------------------------
	@echo test_interpolate starting
	@echo
	@make test_interpolate_exe
	@echo
	@echo './test_interpolate_exe'
	-@(./test_interpolate.exe)
	@echo
	@echo test_interpolate done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_RAY_TRACE = \
	../src/ModMpiInterfaces.o \
	../src/ModMpi.o \
	../src/CON_ray_trace.o

test_ray_trace_exe:
	@(cd ../src; make CON_ray_trace.o)
	@make test_ray_trace.o
	${LINK.f90} -o test_ray_trace.exe \
		test_ray_trace.o ${OBJ_RAY_TRACE} ${Lflag1}

test_ray_trace:
	@echo -----------------------------------------------------------------
	@echo test_ray_trace starting
	@echo
	@make test_ray_trace_exe
	@echo
	@echo '${MPIRUN} test_ray_trace.exe > test_ray_trace.tmp'
	-@${MPIRUN} test_ray_trace.exe > test_ray_trace.tmp
	@echo 'diff test_ray_trace.tmp test_ray_trace.log'
	-@(diff test_ray_trace.tmp test_ray_trace.log)
	@echo
	@echo test_ray_trace done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_LINE_EXTRACT = \
	test_line_extract.o \
	../src/ModMpiInterfaces.o \
	../src/ModMpi.o \
	../src/ModSort.o \
	../src/CON_line_extract.o

test_line_extract_exe: ${OBJ_LINE_EXTRACT}
	@(cd ../src; make CON_line_extract.o)
	@make test_line_extract.o
	${LINK.f90} -o test_line_extract.exe ${OBJ_LINE_EXTRACT} ${Lflag1}

test_line_extract:
	@echo -----------------------------------------------------------------
	@echo test_line_extract starting
	@echo
	@make test_line_extract_exe
	@echo
	@echo '${MPIRUN} test_line_extract.exe'
	-@(${MPIRUN} test_line_extract.exe)
	@echo
	@echo test_line_extract done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_GEOPACK = \
	test_geopack.o \
	../src/CON_geopack.o \
	../src/ModCoordTransform.o

test_geopack_exe:
	@(cd ../src; make ModCoordTransform.o CON_geopack.o)
	@make test_geopack.o
	${LINK.f90} -o test_geopack.exe ${OBJ_GEOPACK} ${Lflag2}

test_geopack:
	@echo -----------------------------------------------------------------
	@echo test_geopack starting
	@echo
	@make test_geopack_exe
	@echo
	@echo './test_geopack.exe > test_geopack.tmp'
	-@(./test_geopack.exe > test_geopack.tmp)
	@echo DiffNum.pl -b -r=1e-8 test_geopack.tmp test_geopack.log
	-@(${SCRIPTDIR}/DiffNum.pl -b -r=1e-8 \
		test_geopack.tmp test_geopack.log)
	@echo
	@echo test_geopack done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_AXES = \
	test_axes.o \
	../src/CON_axes.o \
	../src/CON_geopack.o \
	../src/CON_planet.o \
	../src/ModMpi.o \
	../src/ModMpiInterfaces.o \
	../src/ModCoordTransform.o \
	../src/ModIoUnit.o \
	../src/ModReadParam.o \
	../src/ModTimeConvert.o \
	../src/ModUtilities.o \
	../src/ModKind.o \
	../src/ModConst.o \
	../src/ModPlanetConst.o 

test_axes_exe:
	@(cd ../src; make CON_axes.o)
	@make test_axes.o
	${LINK.f90} -o test_axes.exe ${OBJ_AXES} ${Lflag1}

test_axes: 
	@echo -----------------------------------------------------------------
	@echo test_axes starting
	@echo
	@make test_axes_exe
	@echo
	@echo './test_axes.exe'
	-@(./test_axes.exe)
	@echo
	@echo test_axes done
	@echo -----------------------------------------------------------------

#######################################################################
OBJ_PLANET_FIELD = \
	test_planet_field.o \
	../src/CON_planet_field.o \
	../src/CON_planet.o \
	../src/CON_axes.o \
	../src/CON_geopack.o \
	../src/ModMpi.o \
	../src/ModMpiInterfaces.o \
	../src/ModPlanetConst.o \
	../src/ModCoordTransform.o \
	../src/ModIoUnit.o \
	../src/ModReadParam.o \
	../src/ModTimeConvert.o \
	../src/ModUtilities.o

test_planet_field_exe:
	@(cd ../src; make CON_planet_field.o)
	@make test_planet_field.o
	${LINK.f90} -o test_planet_field.exe ${OBJ_PLANET_FIELD} ${Lflag1}

test_planet_field: 
	@echo -----------------------------------------------------------------
	@echo test_planet_field starting
	@echo
	make test_planet_field_exe
	@echo
	@echo './test_planet_field.exe > test_planet_field.tmp'
	-(./test_planet_field.exe > test_planet_field.tmp)
	@echo DiffNum.pl -b -r=1e-6 test_planet_field.tmp test_planet_field.log
	-(${SCRIPTDIR}/DiffNum.pl -b -r=1e-6 \
		 test_planet_field.tmp test_planet_field.log)
	@echo
	@echo test_planet_field done
	@echo -----------------------------------------------------------------
#######################################################################

OBJ_PLOT_FILE = \
	test_plot_file.o \
	../src/ModIoUnit.o \
	../src/ModUtilities.o \
	../src/ModPlotFile.o

test_plot_file_exe:
	@(cd ../src; make ModPlotFile.o)
	@make test_plot_file.o
	${LINK.f90} -o test_plot_file.exe ${OBJ_PLOT_FILE} ${Lflag2}

test_plot_file:
	@echo -----------------------------------------------------------------
	@echo test_plot_file starting
	@echo
	@make test_plot_file_exe
	@echo
	@echo './test_plot_file.exe'
	-@(./test_plot_file.exe)
	@echo
	@echo test_plot_file done
	@echo -----------------------------------------------------------------
#######################################################################

OBJ_LOOKUP_TABLE = \
	test_lookup_table.o \
	../src/ModIoUnit.o \
	../src/ModMpi.o \
	../src/ModMpiInterfaces.o \
	../src/ModReadParam.o \
	../src/ModInterpolate.o \
	../src/ModUtilities.o \
	../src/ModPlotFile.o \
	../src/ModLookupTable.o

test_lookup_table_exe:
	@(cd ../src; make ModLookupTable.o)
	@make test_lookup_table.o
	${LINK.f90} -o test_lookup_table.exe ${OBJ_LOOKUP_TABLE} ${Lflag1}

test_lookup_table:
	@echo -----------------------------------------------------------------
	@echo test_lookup_table starting
	@echo
	@make test_lookup_table_exe
	@echo
	@echo '${MPIRUN} test_lookup_table.exe'
	-@(${MPIRUN} test_lookup_table.exe) | grep -i -C2 'error' \
		> test_lookup_table.diff
	@echo
	@echo test_lookup_table done
	@echo -----------------------------------------------------------------
#######################################################################

OBJ_UTILITIES = \
	test_utilities.o \
	../src/ModUtilities.o 

test_utilities_exe:
	@(cd ../src; make ModUtilities.o)
	@make test_utilities.o
	${LINK.f90} -o test_utilities.exe ${OBJ_UTILITIES} ${Lflag2}

test_utilities:
	@echo -----------------------------------------------------------------
	@echo test_utilities starting
	@echo
	@make test_utilities_exe
	@echo
	@echo './test_utilities.exe'
	-@(./test_utilities.exe)
	@echo
	@echo test_utilities done
	@echo -----------------------------------------------------------------

distclean: clean
	rm -f *.tmp
