# C language related part of Makefile.conf

COMPILE.c     = cc
COMPILE.mpicc = mpicc
COMPILE.mpicxx= mpicxx

DEBUGC = 
#DEBUGC = -g

.SUFFIXES: .c .cpp

FLAGC = ${SEARCH_C} -c ${OPT3} ${DEBUGC}

.c.o:
	${COMPILE.c} ${FLAGC} $<

.cpp.o:
	${COMPILE.mpicxx} ${FLAGC} $<
