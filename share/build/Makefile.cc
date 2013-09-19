# C language related part of Makefile.conf

COMPILE.c     = cc
COMPILE.mpicc = mpicc
COMPILE.mpicxx= mpicxx

DEBUGC = 
#DEBUGC = -g

.SUFFIXES: .c

FLAGC = ${SEARCH} -c ${OPT3} ${DEBUGC}

.c.o:
	${COMPILE.c} ${FLAGC} ${SEARCH} $<

.cpp.o:
	${COMPILE.mpicxx} ${FLAGC} ${SEARCH} $<
