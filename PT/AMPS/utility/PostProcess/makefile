SHELL =/bin/sh

DEFAULT_TARGET : pproc  

include ../../Makefile.conf

FLAGC+= -I../../src/meshAMR -I../../src/general -I../../src/pic 
MAINOBJ=PostProcess3D.o PostProcess3D_ColumnIntegration.o PostProcess3D_SurfaceData.o  PostProcess3D_ParticleTrajectory.o 

pproc : $(MAINOBJ)
	ar -rc pproc.a $(MAINOBJ) 

