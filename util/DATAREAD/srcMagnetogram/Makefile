#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libINDICES.a

OBJECTS = \
	main.o \
	ModMagHarmonics.o

main.o: ModMagHarmonics.o

HARMONICS:
	make ${EXE}
	@echo ${EXE} has been compiled.

EXE = ${BINDIR}/HARMONICS.exe

${EXE}: ${OBJECTS}
	${LINK.f90} -o ${EXE} ${OBJECTS} ${Lflag1}

distclean: clean
