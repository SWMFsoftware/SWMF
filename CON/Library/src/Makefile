#^CFG COPYRIGHT UM

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

