#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf

OBJ_RAY_TEST = \
	test_ray_trace.o \
	../src/ModMpi.o \
	../src/CON_ray_trace.o

test_ray_trace.exe: ${OBJ_RAY_TEST}
	cd ../src; make ModMpi.o CON_ray_trace.o
	${LINK.f90} ${Lflag1} -o test_ray_trace.exe ${OBJ_RAY_TEST}

test_ray_trace: test_ray_trace.exe
	mpirun -np 2 test_ray_trace.exe > test_ray_trace.tmp
	@echo 'There should be no difference:'
	diff test_ray_trace.tmp test_ray_trace.log

distclean: clean
	rm -f *.tmp

# keep this line

