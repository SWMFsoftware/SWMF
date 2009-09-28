#^CFG COPYRIGHT UM

SHELL=/bin/sh

#
#	Tru64 Unix Operating System, f90 - Compaq Fortran 90 Compiler
#	Space Weather Modeling Framework (SWMF) main Makefile
#

COMPILE.f77	= ${CUSTOMPATH_F}f90
COMPILE.f90     = ${CUSTOMPATH_F}f90
LINK.f90	= ${CUSTOMPATH_MPI}f90
AR = ar -rs

SINGLEPREC = -real_size 32
DOUBLEPREC = -real_size 64
PRECISION  = ${DOUBLEPREC}

LOCAL_LIBS = /usr/local/lib
MPILIB = -L$(LOCAL_LIBS) -lmpi -lelan	#For 64 bit compilation
#MPILIB = -L${LIBDIR} -lNOMPI

# This is the search path for used modules
# SEARCH_EXTRA should be given in the individual Makefiles

SEARCH = -I${SHAREDIR} ${SEARCH_EXTRA}

DEBUGFLAG = -C -g
DEBUG     =

OPT0 = -O0
OPT1 = -O1
OPT2 = -O2
OPT3 = -O3
OPT4 = -fast

ARCH	= ev6 
TUNE	= -tune ${ARCH} 
BITFLAG = -arch ${ARCH}

# For some parts (GITM, RCM) of the code the -fpe1 flag is needed 
# to avoid run time errors due to floating point underflows
CFLAG = ${SEARCH} -c -align dcommons ${TUNE} ${BITFLAG} ${DEBUG} -fpe1

Cflag0  = ${CFLAG} ${PRECISION} ${OPT0}
Cflag1  = ${CFLAG} ${PRECISION} ${OPT1}
Cflag2  = ${CFLAG} ${PRECISION} ${OPT2}
Cflag3  = ${CFLAG} ${PRECISION} ${OPT3}
Cflag4  = ${CFLAG} ${PRECISION} ${OPT4}

# RCM compilation flags
# To allow RCM to compile as double precision, add PRECISION flag
CFLAGS = ${SEARCH} -c -warn argument_checking -fpe1 ${DEBUG}

Lflag1  = ${OPTSAFE} ${PRECISION} ${BITFLAG} ${MPILIB} 
Lflag2  = ${OPT} ${PRECISION} ${BITFLAG} 

# BLAS and LAPACK libraries
LBLAS =
BLAS  = lapack.o blas.o


#
#       General rules
#

.SUFFIXES:
.SUFFIXES: .f90 .F90 .f .for .ftn .o

.f90.o:
	${COMPILE.f90} ${Cflag2} $<

.F90.o:
	${COMPILE.f90}  -DsysOSF1 ${Cflag2} $<


# Compiling the F77 source with Clfag3 results in a runtime error for the
# UA component. Do not change this unless you know what you are doing...

.f.o:
	${COMPILE.f90} ${Cflag2} -extend_source $<

.for.o:
	${COMPILE.f90} ${Cflag2} -extend_source $<

.ftn.o:
	${COMPILE.f90} ${Cflag2} -extend_source $<

clean:	
	rm -f *~ core *.o *.mod fort.* a.out *.exe *.a *.so *.protex

distclean: clean

# keep this line
