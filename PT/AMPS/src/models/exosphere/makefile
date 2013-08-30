SHELL =/bin/sh

DEFAULT_TARGET : exosphere

clean:
	rm -f *.o *.a

CC=gcc
Flags=-g   
IncludeList=

.SUFFIXES: .cpp .o 
.cpp.o: 
	${CC} -c ${Flags} ${IncludeList} $< 


EXOSPHEREOBJ = Exosphere.o \
	Exosphere_OrbitalMotion.o \
	Exosphere_SourceProcesses.o \
	Exosphere_Parallel.o \
	Exosphere_ReferenceGroundBasedObservations.o \
	Exosphere_ChamberlainExosphereModel.o 

exosphere : $(EXOSPHEREOBJ)
	ar -src exosphere.a $(EXOSPHEREOBJ)
	@echo "Target 'exosphere' is done." 
