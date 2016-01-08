SHELL =/bin/sh

DEFAULT_TARGET : pproc  

include ../../Makefile.conf

MAINLIBOBJ = PostProcess.o PostProcess3D_ColumnIntegration.o main.o 

lib : $(MAINLIBOBJ)
	echo "fff"

pproc :
	rm -f *.o
	mpicxx -g -c -I../../src/meshAMR -I../../src/general PostProcess.cpp
	mpicxx -g -c -I../../src/meshAMR -I../../src/general PostProcess3D_ColumnIntegration.cpp
	mpicxx -g -c -I../../src/meshAMR -I../../src/general main.cpp
	mpicxx -g -c -I../../src/meshAMR -I../../src/general ../../src/meshAMR/meshAMRcutcell.cpp
	mpicxx -g -c -I../../src/meshAMR -I../../src/general ../../src/general/global.cpp
	mpicxx -g -c -I../../src/meshAMR -I../../src/general ../../src/general/specfunc.cpp
	mpicxx -g -c -I../../src/meshAMR -I../../src/general ../../src/general/rnd.cpp
	mpicxx *.o 


