SHELL =/bin/sh

DEFAULT_TARGET : pproc  

MAINLIBOBJ = PostProcess.o PostProcess3D_ColumnIntegration.o main.o 

pproc :
	rm -f *.o
	mpicxx -g -c -I/Users/vtenishe/SPICE/Toolkit/cspice/include -I../../../utility/PostProcess -I../../../src/meshAMR -I../../../src/general main.cpp
	cd ../../../utility/PostProcess; make
	mpicxx main.o ../../../utility/PostProcess/pproc.a /Users/vtenishe/SPICE/Toolkit/cspice/lib/cspice.a 


