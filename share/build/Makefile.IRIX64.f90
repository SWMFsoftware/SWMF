#^CFG COPYRIGHT UM

SHELL = /bin/sh

#
#	Space Weather Modeling Framework (SWMF) 
#	IRIX64 Operating System, MIPSpro Fortran 90 compiler
#           IRIX64 specific part of Makefile
#

#MF = Makefile.${OS}
MF = 

COMPILE.f77     = ${CUSTOMPATH_F}f90
COMPILE.f90     = ${CUSTOMPATH_F}f90
LINK.f90        = ${CUSTOMPATH_F}f90
AR = ar -rs

SINGLEPREC =
DOUBLEPREC = -r8
PRECISION  = ${DOUBLEPREC}

MPILIB = -lmpi
#MPILIB = -L${LIBDIR} -lNOMPI

# This is the search path for used modules
# SEARCH_EXTRA should be given in the individual Makefiles

SEARCH = -I${SHAREDIR} ${SEARCH_EXTRA}

#
# DEBUG flags
#
# -------------------------------------------------------------------
# PLEASE NOTE: In order for subscript_check to work, you need to do:
# setenv F90_BOUNDS_CHECK_ABORT YES
# -------------------------------------------------------------------

DEBUGFLAG = -C -DEBUG:trap_uninitialized=ON
DEBUG     = 

BITLIB = -L/usr/lib32
BITFLAG = -n32
#BITLIB = -L/usr/lib64
#BITFLAG = -64

# Option        Action
#----------------------------------------------------------------------------
# -O0           No optimization.
# -O1           Local optimization.
# -O2, -O       Extensive optimization.  
# -O3           Aggressive optimization.  
# -Ofast[=ipxx] Maximizes performance for the target platform ipxx (use hinv)

# -------------------------------------------------------------------
# PLEASE NOTE: use -O2 or less for double PRECISION = -r8
# -------------------------------------------------------------------

OPT0 = -O0
OPT1 = -O1
OPT2 = -O2
OPT3 = -O3
OPT4 = -O4

CFLAG = ${SEARCH} -c ${BITFLAG} ${DEBUG}

Cflag0  = ${CFLAG} ${PRECISION} ${OPT0}
Cflag1  = ${CFLAG} ${PRECISION} ${OPT1}
Cflag2  = ${CFLAG} ${PRECISION} ${OPT2}
Cflag3  = ${CFLAG} ${PRECISION} ${OPT3}
Cflag4  = ${CFLAG} ${PRECISION} ${OPT4}

# RCM compilation flags
# To allow RCM to compile as double precision, add PRECISION flag
CFLAGS = ${CFLAG} ${OPT2} -w

Lflag1  = ${BITFLAG} ${MPILIB} ${BITLIB} -woff136
Lflag2  = ${BITFLAG} ${BITLIB} -woff136

# BLAS and LAPACK libraries
LBLAS =
BLAS  = blas.o lapack.o

# Using the SGI blas library is incompatible with the GM->IH renaming 
# right now. Our timings did not show any measurable speed up relative
# to the use of our own blas.f and lapack.f subroutines. So the SGI blas is 
# commented out for now. We know how to make this work if required.
# LBLAS = -lblas_mp
# BLAS  = lapack.o

#
#       General rules
#

.SUFFIXES:
.SUFFIXES: .f90 .F90 .f .for .ftn .o

.f90.o:
	${COMPILE.f90} ${Cflag3} $<

.F90.o:
	${COMPILE.f90} -DsysIRIX64 ${Cflag3} $<

.f.o:
	${COMPILE.f77} ${Cflag3} -extend_source $<

.for.o:
	${COMPILE.f77} ${Cflag2} -extend_source $<

.ftn.o:
	${COMPILE.f77} ${Cflag2} -extend_source $<

clean:	
	rm -f *~ core *.o *.mod fort.* a.out *.exe *.a *.so *.protex


# keep this line
