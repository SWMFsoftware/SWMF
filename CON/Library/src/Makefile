#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND

# These modules provide the "infrastructure" for SWMF.

MODULES = \
	CON_comp_info.o\
	CON_comp_param.o\
	CON_world.o\
	CON_time.o\
	CON_physics.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${MODULES}

MY_LIB = ${LIBDIR}/libLIBRARY.a

LIB:	DEPEND
	$(MAKE) ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES}

test: test_registry

test_registry: DEPEND
	cd ${SHAREDIR}; $(MAKE) LIB
	$(MAKE) LIB
	$(MAKE) test_registry.exe
	mpirun -np 4 test_registry.exe
	cat LAYOUT.in

TEST_REGISTRY_OBJ= \
	test_registry.o

test_registry.exe: DEPEND ${TEST_REGISTRY_OBJ}
	cd ${SHAREDIR}; $(MAKE) LIB
	cd ${LIBRARYDIR}; $(MAKE) LIB
	${LINK.f90} -o test_registry.exe ${SEARCH} \
		${TEST_REGISTRY_OBJ} -L${LIBDIR} -lLIBRARY -lSHARE

install: 

distclean: clean
	rm -f Makefile.DEPEND 

# keep this line

