#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND

# These modules contain the coupling toolkit of SWMF.

MODULES = \
	ModInitGridStorage.o\
	CON_grid_storage.o\
	CON_grid_descriptor.o\
	CON_test_global_message_pass.o\
	CON_coupler.o\
	CON_router.o\
	CON_global_message_pass.o\
	CON_domain_decomposition.o

SEARCH_EXTRA = -I${LIBRARYDIR}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${SEARCH_EXTRA} ${MODULES}

MY_LIB = ${LIBDIR}/libCOUPLER.a

LIB:	DEPEND
	make ${MY_LIB}
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

