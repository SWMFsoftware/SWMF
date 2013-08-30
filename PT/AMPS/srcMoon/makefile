SHELL =/bin/sh

DEFAULT_TARGET : lib  

clean:
	rm -f *.o *.a

CC=gcc
Flags=-g   
IncludeList=

.SUFFIXES: .cpp .o 
.cpp.o: 
	${CC} -c ${Flags} ${IncludeList} $< 


MAINLIBOBJ = Exosphere_Helium.o Moon.o Moon_SampleVelocityDistribution.o Moon_SampleSubSolarLimbColumnIntegrals.o tvis.Kaguya.o main_lib.o 
MAINOBJ = main.o


lib : $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'main' is done." 

amps : $(MAINOBJ) 
	ar -rc main.a $(MAINOBJ)

