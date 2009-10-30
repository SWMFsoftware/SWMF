
SHELL =/bin/sh

include ../Makefile.def

include ../Makefile.conf

include Makefile.DEPEND

install:
	touch Makefile.DEPEND

OBJECTS = \
	BATL_amr.o \
	BATL_geometry.o \
	BATL_grid.o \
	BATL_lib.o  \
	BATL_mpi.o  \
	BATL_pass_cell.o \
	BATL_pass_face.o \
	BATL_size.o \
	BATL_tree.o

ALLOBJECTS = \
	${OBJECTS} \
	main.o \
	advect_main.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${ALLOBJECTS}

BATL_size.f90: BATL_size_orig.f90
	cp -f BATL_size_orig.f90 BATL_size.f90

MY_LIB = libBATL.a

LIB: DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

BATL:
	make DEPEND
	make ${BINDIR}/BATL.exe

${BINDIR}/BATL.exe: main.o ${OBJECTS}
	${LINK.f90} -o ${BINDIR}/BATL.exe main.o ${OBJECTS} \
		-L${LIBDIR} -lTIMING -lSHARE ${Lflag1}

ADVECT:
	make DEPEND
	make ${BINDIR}/ADVECT.exe

${BINDIR}/ADVECT.exe: advect_main.o ${OBJECTS}
	${LINK.f90} -o ${BINDIR}/ADVECT.exe advect_main.o ${OBJECTS} \
		-L${LIBDIR} -lTIMING -lSHARE ${Lflag1}

distclean: clean
	rm -f BATL_size.f90 Makefile.DEPEND

