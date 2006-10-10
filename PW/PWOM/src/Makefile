#FTN = g77
#DEBUG = -g -C -O3 -ffixed-line-length-132

#FTN = f95
#DEBUG = -w -gline -nan -C=all -O0 -r8
#DEBUG = -w -O3 -nan -r8 -ieee=nonstd
#DEBUG = -w -nan -C=all -O3 -r8
#DEBUG = -w -nan -C=all -O3 -r8 -ieee=stop

FTN = mpif90
DEBUG = -w -gline -nan -C=all -O0 -r8
#DEBUG = -w -nan -C=all -O3 -r8 -ieee=stop
#FTN = ifort
#DEBUG = -132

Cflag = -c ${DEBUG}
Lflag = ${DEBUG}

EXE = ../bin/pw3.exe

default: ${EXE}

include Makefile.planet

MODULES = \
	ModParameters.o\
	ModCommonPlanet.o\
	ModCommonVariables.o\
	Mod_PW.o\
	ModNeutralPW.o\
	ModPass.o
METHODS = \
	pw.o\
	polar_wind.o\
	solver.o\
	cvmgp.o\
	prntel.o\
	prntef.o\
	prnthp.o\
	prntop.o\
	prntsm.o\
	elfldw.o\
	clfme2w.o\
	calcdt.o\
	diffusion.o\
	PRNT_Sources.o\
	calc_cfl.o\
	calc_diffusion.o\
	solver_rusanov.o\
	Get_Neutral.o\
	${OBJECTS_PLANET}

ModCommonPlanet.o: ModParameters.o

ModCommonVariables.o: ModParameters.o ModCommonPlanet.o

${METHODS}: ${MODULES}

OBJECTS = ${MODULES} ${METHODS}

${EXE}: ${OBJECTS}
	${FTN} ${Lflag} -o ${EXE} ${OBJECTS}

UNSPLIT: pw3.o
	${FTN} ${Lflag} -o ../bin/pw3.exe pw3.o

pw3.o: common_variables.f Makefile

.SUFFIXES:
.SUFFIXES: .f .f90 .o

.f.o:
	${FTN} ${Cflag} -132 $<

.f90.o:
	${FTN} ${Cflag} $<

clean:
	rm -f *.o *~ *.mod  

distclean: clean
	rm -f *_planet.f Makefile.planet PLANET *.exe fort.*
