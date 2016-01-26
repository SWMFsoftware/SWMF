SHELL =/bin/sh

DEFAULT_TARGET : pproc  

pproc :
	rm -f *.o
	./Config.pl
	mpicxx -fopenmp  -g -c -I/Users/vtenishe/SPICE/Toolkit/cspice/include -I../utility/PostProcess -I../src/meshAMR -I../src/general main.cpp DustScatteringEfficientcy.cpp VIRTIS-M.cpp 
	cd ../utility/PostProcess; make
	mpicxx main.o DustScatteringEfficientcy.o VIRTIS-M.o ../utility/PostProcess/pproc.a /Users/vtenishe/SPICE/Toolkit/cspice/lib/cspice.a 


