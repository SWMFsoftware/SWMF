#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libINDICES.a

include Makefile.DEPEND
#include Makefile.RULES

MODULES = \
	ModIndicesInterfaces.o\
	ModIndices.o

OBJECTS = \
	read_NGDC_indices.o\
	read_SWPC_indices.o\
	read_NOAAHPI_indices.o\
	read_MHDIMF_indices.o\
	insert_into_indices_array.o\
	indices_set_inputs.o\
	indices_library.o\
	IO_library.o\
	IO_set_inputs.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl -I${SHAREDIR} ${MODULES}

LIB:    DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${MODULES} ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES} ${OBJECTS}

distclean: clean
	rm -f Makefile.DEPEND
