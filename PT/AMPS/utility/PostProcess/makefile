SHELL =/bin/sh

DEFAULT_TARGET : pproc  
GENERALOBJ=PostProcess3D.o PostProcess3D_ColumnIntegration.o PostProcess3D_SurfaceData.o PostProcess3D_ParticleTrajectory.o 

include ../../../Makefile.conf

pproc : $(GENERALOBJ) 
	ar -src pproc.a $(GENERALOBJ)


