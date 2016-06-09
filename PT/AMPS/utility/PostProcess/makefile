SHELL =/bin/sh

DEFAULT_TARGET : pproc  
FLAGS=

include ../../Makefile.conf

pproc :
	rm -f *.o
	mpicxx -g ${FLAGS} -c -I../../src/meshAMR -I../../src/general PostProcess3D.cpp
	mpicxx -g ${FLAGS} -c -I../../src/meshAMR -I../../src/general PostProcess3D_ColumnIntegration.cpp
	mpicxx -g ${FLAGS} -c -I../../src/meshAMR -I../../src/general PostProcess3D_ParticleTrajectory.cpp
	mpicxx -g ${FLAGS} -c -I../../src/meshAMR -I../../src/general PostProcess3D_SurfaceData.cpp
	mpicxx -g ${FLAGS} -c -I../../src/meshAMR -I../../src/general ../../src/meshAMR/meshAMRcutcell.cpp
	mpicxx -g ${FLAGS} -c -I../../src/meshAMR -I../../src/general ../../src/general/global.cpp
	mpicxx -g ${FLAGS} -c -I../../src/meshAMR -I../../src/general ../../src/general/specfunc.cpp
	mpicxx -g ${FLAGS} -c -I../../src/meshAMR -I../../src/general ../../src/general/rnd.cpp

	ar -rc pproc.a PostProcess3D.o PostProcess3D_ColumnIntegration.o PostProcess3D_SurfaceData.o meshAMRcutcell.o global.o specfunc.o rnd.o PostProcess3D_ParticleTrajectory.o


