#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND

# These modules contain the coupling toolkit of SWMF.

MODULES = \
	CON_domain_decomposition.o\
	ModInitGridStorage.o\
	CON_grid_storage.o\
	CON_grid_descriptor.o\
	CON_global_vector.o\
	CON_router.o\
	CON_global_message_pass.o\
	CON_integrator.o\
	CON_buffer_grid.o\
	CON_test_global_message_pass.o\
	CON_coupler.o\
	CON_couple_points.o\
	CON_transfer_data.o

SEARCHDIR = -I${SHAREDIR} -I${LIBRARYDIR}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${MODULES}

MY_LIB = ${LIBDIR}/libCOUPLER.a

LIB:	DEPEND
	$(MAKE) ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES}

install: 

distclean: clean
	rm -f Makefile.DEPEND 

# keep this line

