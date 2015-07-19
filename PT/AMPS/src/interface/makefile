SHELL =/bin/sh

DEFAULT_TARGET : interface

include ../../Makefile.conf

CELL_CENTERED_LINEAR_INTERPOLATION_OBJ=cell_centered_linear_interpolation_cpp.o \
	cell_centered_linear_interpolation_f90.o \
	../../share/Library/src/ModKind.o \
	../../share/Library/src/ModInterpolateAMR.o

CELL_CENTERED_LINEAR_INTERPOLATION : $(CELL_CENTERED_LINEAR_INTERPOLATION_OBJ)

interface :
	make CELL_CENTERED_LINEAR_INTERPOLATION
	${AR} -src interface.a ./*.o
	@echo "Target 'interface' is done." 