SHELL =/bin/sh

DEFAULT_TARGET : pproc  
FLAGS=-O3 -fopenmp

include ../../Makefile.conf

pproc :
	rm -f *.o
	mpicxx  ${FLAGS} -c -I../../src/meshAMR -I../../src/general -I../../src/pic PostProcess3D.cpp
	mpicxx  ${FLAGS} -c -I../../src/meshAMR -I../../src/general -I../../src/pic PostProcess3D_ColumnIntegration.cpp
	mpicxx  ${FLAGS} -c -I../../src/meshAMR -I../../src/general -I../../src/pic PostProcess3D_ParticleTrajectory.cpp
	mpicxx  ${FLAGS} -c -I../../src/meshAMR -I../../src/general -I../../src/pic PostProcess3D_SurfaceData.cpp
	mpicxx  ${FLAGS} -c -I../../src/meshAMR -I../../src/general ../../src/general/global.cpp
	mpicxx  ${FLAGS} -c -I../../src/meshAMR -I../../src/general ../../src/general/specfunc.cpp
	mpicxx  ${FLAGS} -c -I../../src/meshAMR -I../../src/general ../../src/general/rnd.cpp

	ar -rc pproc.a PostProcess3D.o PostProcess3D_ColumnIntegration.o PostProcess3D_SurfaceData.o  global.o specfunc.o rnd.o PostProcess3D_ParticleTrajectory.o 

