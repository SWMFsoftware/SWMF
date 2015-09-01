
INCLUDE = -I../include

OBJECTS = \
  test_arrays.o \
  asserts.o \
  errors.o \
  MPIdata.o \
  debug.o

# if compiling on Xeon
#
#PERFORMANCE_FLAGS = -O3 -fno-exceptions -xHost -vec-report
#MPI_LIBS =  -I/opt/intel/impi/4.1.2.040/intel64/include
#MPI_LINK = \
#  /opt/intel/impi/4.1.2.040/intel64/lib/libmpigc4.so \
#  /opt/intel/impi/4.1.2.040/intel64/lib/libmpigf.so \
#  /opt/intel/impi/4.1.2.040/intel64/lib/libmpi.so \
#  /opt/intel/impi/4.1.2.040/intel64/lib/libmpigi.a
#MPIRUN = mpirun

# if compiling on my macbook
#
PERFORMANCE_FLAGS = -O3 -fno-exceptions
MPI_LIBS =  -I/opt/local/include/openmpi
MPI_LINK = \
  /opt/local/lib/libmpi_cxx.dylib \
  /opt/local/lib/libmpi.dylib
MPIRUN = mpiexec

# need -lrt to handle clock_gettime()
ARRAY_FLAGS = #-DCHAINED_ARRAYS #-DFLAT_ARRAYS #-DCHECK_BOUNDS
FLAGS = $(PERFORMANCE_FLAGS) -DNO_MPI $(ARRAY_FLAGS) #-ggdb 

COMPILER = g++ # icpc # c++

	#$(COMPILER) -E testcomm.cpp $(INCLUDE) $(MPI_LIBS)
test: test.cpp
	$(COMPILER) test.cpp $(INCLUDE) $(PERFORMANCE_FLAGS) $(MPI_LIBS) $(MPI_LINK) -o test

run.test: test
	$(MPIRUN) -n 1 ./test

testcomm: testcomm.cpp
	$(COMPILER) testcomm.cpp $(INCLUDE) $(PERFORMANCE_FLAGS) $(MPI_LIBS) $(MPI_LINK) -o testcomm

run.testcomm: testcomm
	$(MPIRUN) -n 2 ./testcomm

testvbasic.mic: testvbasic.cpp
	icpc -mmic testvbasic.cpp $(INCLUDE) -o testvbasic.mic

run.testvbasic.mic: testvbasic.mic
	ssh `hostname`-mic1 tests/testvbasic.mic

testvrinke.mic: testvrinke.cpp
	KMP_AFFINITY=balanced,granularity=thread icpc -mmic $(FLAGS) -fopenmp testvrinke.cpp $(INCLUDE) -o testvrinke.mic

run.testvrinke.mic: testvrinke.mic
	OMP_NUM_THREADS=50 mpirun -n 1 -host `hostname`-mic1 -env OMP_NUM_THREADS=50 ./testvrinke.mic

test_mic_particles.mic: test_mic_particles.cpp Makefile
	KMP_AFFINITY=compact icpc -mmic $(FLAGS) -fopenmp test_mic_particles.cpp \
           $(INCLUDE) -o test_mic_particles.mic

run.test_mic_particles.mic: test_mic_particles.mic
	mpirun -n 1 -host `hostname`-mic1 -env OMP_NUM_THREADS=2 ./test_mic_particles.mic

testaa: testaa.cpp Makefile
	$(COMPILER) $(FLAGS) -fopenmp testaa.cpp $(INCLUDE) -o testaa

testv2: testv2.cpp Makefile
	$(COMPILER) $(FLAGS) -fopenmp testv2.cpp $(INCLUDE) -o testv2

#icpc -mmic -S -O3 testv2.cpp -o testv2.mic.s
testv2.mic: testv2.cpp Makefile
	KMP_AFFINITY=compact icpc -mmic $(FLAGS) -fopenmp testv2.cpp $(INCLUDE) -o testv2.mic

#icc -mmic -S -O3 testv2.cpp -std=c99
testv2.xeon: testv2.cpp Makefile
	KMP_AFFINITY=compact icpc $(FLAGS) -fopenmp testv2.cpp $(INCLUDE) -o testv2.xeon

run.testv2.mic: testv2.mic
	mpirun -n 1 -host `hostname`-mic1 -env OMP_NUM_THREADS=2 ./testv2.mic

run.testv2.xeon: testv2.xeon
	mpirun -n 1 -env OMP_NUM_THREADS=2 ./testv2.xeon

testv.mic: testv.cpp Makefile
	icpc testv.cpp -g -O3 -fopenmp -fp-model fast=2 -vec-report6 -mmic -o testv.mic

testv.xeon: testv.cpp Makefile
	icpc testv.cpp -g -O3 -fopenmp -vec-report6 -o testv.xeon

vtune.mic: testv.mic Makefile
	amplxe-cl -collect knc-general-exploration -knob enable-vpu-metrics=true -knob enable-tlb-metrics=true -- mpirun -host knc1-mic0 -np 1 -env OMP_NUM_THREADS 1 ./testv.mic > vtune.txt

test_arrays: clean_test_arrays $(OBJECTS)
	$(COMPILER) $(FLAGS) $(INCLUDE) $(OBJECTS) -o test_arrays

test_arrays.o:
	$(COMPILER) $(FLAGS) -c test_arrays.cpp $(INCLUDE) -o test_arrays.o


asserts.o:
	$(COMPILER) $(FLAGS) -c ../utility/asserts.cpp $(INCLUDE) -o asserts.o

debug.o:
	$(COMPILER) $(FLAGS) -c ../utility/debug.cpp $(INCLUDE) -o debug.o

errors.o:
	$(COMPILER) $(FLAGS) -c ../utility/errors.cpp $(INCLUDE) -o errors.o

MPIdata.o:
	$(COMPILER) $(FLAGS) -c ../utility/MPIdata.cpp $(INCLUDE) -o MPIdata.o

clean: clean_old_test_arrays clean_test_arrays clean_executables
	rm -f test_arrays $(OBJECTS)

clean_executables:
	rm -f *.mic *.xeon

clean_test_arrays:
	rm -f test_arrays test_arrays.o

old_test_arrays: clean_old_test_arrays
	$(COMPILER) $(FLAGS) old_test_arrays.cpp  -o old_test_arrays

clean_old_test_arrays:
	rm -f old_test_arrays old_test_arrays.o
