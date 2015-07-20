SHELL =/bin/sh

DEFAULT_TARGET : interface

include ../../Makefile.conf
include ../../Makefile.def

SHAREDOBJ=${SHAREDIR}/ModKind.o \
	${SHAREDIR}/ModInterpolateAMR.o

CELL_CENTERED_LINEAR_INTERPOLATION_OBJ=cell_centered_linear_interpolation_cpp.o \ 
	cell_centered_linear_interpolation_f90.o 

SHARED : $(SHAREDOBJ)

CELL_CENTERED_LINEAR_INTERPOLATION : $(CELL_CENTERED_LINEAR_INTERPOLATION_OBJ)

interface :
	echo $(SHAREDIR)
	make CELL_CENTERED_LINEAR_INTERPOLATION
	make SHARED
	${AR} -src interface.a ./*.o
	@echo "Target 'interface' is done." 
