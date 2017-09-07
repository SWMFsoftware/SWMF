SHELL =/bin/sh

DEFAULT_TARGET : pproc  

include ../../Makefile.conf
include ../../Makefile.local

FLAGC+= -I../../src/meshAMR -I../../src/general -I../../src/pic  ${EXTRACOMPILEROPTIONS}
ifeq ($(COMPILE.c),pgcc)
FLAGC+=-mp
else
FLAGC+=-fopenmp
endif


MAINOBJ=PostProcess3D.o PostProcess3D_ColumnIntegration.o PostProcess3D_SurfaceData.o  PostProcess3D_ParticleTrajectory.o 

pproc : $(MAINOBJ)
	ar -rc pproc.a $(MAINOBJ) 

clean:
	rm -rf *.o
	rm -rf *.a