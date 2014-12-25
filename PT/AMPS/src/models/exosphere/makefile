SHELL =/bin/sh

DEFAULT_TARGET : exosphere

include ../../../Makefile.conf

EXOSPHEREOBJ = Exosphere.o \
	Exosphere_OrbitalMotion.o \
	Exosphere_SourceProcesses.o \
	Exosphere_Parallel.o \
	Exosphere_ReferenceGroundBasedObservations.o \
	Exosphere_ChamberlainExosphereModel.o \
	Exosphere_BackgroundPlasmaSource.o

exosphere : $(EXOSPHEREOBJ)
	ar -src exosphere.a $(EXOSPHEREOBJ)
	@echo "Target 'exosphere' is done." 
