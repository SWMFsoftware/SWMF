# Copyright (C) 2002 Regents of the University of Michigan, 
# portions used with permission 
# For more information, see http://csem.engin.umich.edu/tools/swmf

DEFAULT_TARGET = PARMISAN
DEFAULT_EXE    = PARMISAN.exe

default : ${DEFAULT_TARGET}

include Makefile.def
include Makefile.conf

# Menu of make options
help:
	@echo ' '
	@echo '  You can "make" the following:'
	@echo ' '
	@echo '    <default> ${DEFAULT_TARGET} in stand alone mode'
	@echo ' '
	@echo '    help          (show makefile option list)'
	@echo '    install       (install PARMISAN)'
	@echo ' '
	@echo '    LIB           (Component library libSP for SWMF)'
	@echo '    PARMISAN       (make PARMISAN.exe)'
	@echo '    NOMPI         (NOMPI library for compilation without MPI)'
	@echo ' '
	@echo '    rundir        (create run directory for standalone or SWMF)'
	@echo '    rundir RUNDIR=run_test (create run directory run_test)'
	@echo ' '
	@echo '    test          (run all tests)'
	@echo ' '
	@echo '    clean         (remove temp files like: *~ *.o etc)'
	@echo '    distclean     (equivalent to ./Config.pl -uninstall)'


install: src/ModSize.f90

src/ModSize.f90: src/ModSize_orig.f90
	cp -f src/ModSize_orig.f90 src/ModSize.f90

LIB:    install
	cd src;          make LIB
	cd srcInterface; make LIB

PARMISAN:
	cd ${SHAREDIR}; ${MAKE} LIB
	cd ${EMPIRICALCRDIR}; ${MAKE} LIB
	cd ${TIMINGDIR}; ${MAKE} LIB
	cd src; ${MAKE} LIB
	cd src; ${MAKE} PARMISAN
