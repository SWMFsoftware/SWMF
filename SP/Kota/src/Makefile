#^CFG COPYRIGHT UM

SHELL =/bin/sh
include ../Makefile.def
include ../Makefile.conf

OBJECTS = sep19.o

SEPEXE = ${BINDIR}/SEP.exe

SEP: ${SEPEXE}

${SEPEXE}:	${OBJECTS} main.o
	${LINK.f90} ${SEARCH} -o ${SEPEXE} ${OBJECTS} main.o ${SHAREDIR}/io_unit_for77.o ${SHAREDIR}/ModIoUnit.o 

io_unit_for77.f: 
	cp  ${SHAREDIR}/io_unit_for77.f .

stdout.h:
	cp  ${SHAREDIR}/stdout.h .

SEP77:	io_unit_for77.f stdout.h
	f77 io_unit_for77.f main.f sep19.f
