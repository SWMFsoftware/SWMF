#^CFG COPYRIGHT UM

SHELL =/bin/sh
include ../../../Makefile.def
include ../../../Makefile.conf

OBJECTS = sep19.o

SEPEXE = ${BINDIR}/SEP.exe

SEP: ${SEPEXE}

${SEPEXE}:	${OBJECTS} main.o
	${LINK.f90} -o ${SEPEXE} ${OBJECTS} main.o

clean:
	rm -f *.o
