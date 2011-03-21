# This is the SWMF style Makefile for Randall LeVeque's CLAWPACK.

SHELL = /bin/sh

include ../../Makefile.def
include ../../Makefile.conf

# Force it to be single precision
#PRECISION =

OBJECTS = \
  qinit.o \
  rpn2ad1.o \
  rpt2ad1.o \
  b4step2.o \
  setaux.o \
  psi.o \
  claw2ez.o \
  bc2.o \
  claw2.o \
  step2.o \
  step2ds.o \
  dimsp2.o \
  flux2.o \
  copyq2.o \
  limiter.o \
  philim.o \
  src2.o 

# claw needs its own rules

.SUFFIXES:
.SUFFIXES: .f .o

.f.o:
	${COMPILE.f77} ${CFLAGS} ${PRECISION} $<

lib: $(OBJECTS)

xclaw: $(OBJECTS) driver.o
	$(FTN) $(LFLAGS) $(OBJECTS) driver.o -o xclaw

distclean: clean
	/bin/rm -f *~ fort.* xclaw

### DO NOT remove this line - make depends on it ###
