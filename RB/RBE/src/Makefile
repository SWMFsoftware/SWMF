
include ../Makefile.def

SEARCH_EXTRA = -I${EMPIRICALGMDIR} -I${EMPIRICALIEDIR}

include ../Makefile.conf
include Makefile.RULES

OBJECTS = \
	ModRbe.o\
	ModRbSat.o\
	ModWriteTec.o\
	ModPrerunField.o\
	rbe_subs.o\
	plot_potential.o\
	trace_line.o\
	trace_dipole.o\
	plot_fieldline.o\
	set_parameters.o\
	pbo_2.o\
	geopack.o

DEPEND:
	@${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

.f.o:
	${COMPILE.f77} ${CFLAGS} ${PRECISION} $<

.f90.o:
	${COMPILE.f90} ${CFLAGS} ${PRECISION} $<

MY_LIB = libRBE.a

LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS} 

# Object files not included into the library
OBJECTS_EXE = rbe_main.o

# Other requireed libraries
LIBSHARE  = ${LIBDIR}/libSHARE.a
LIBTIMING = ${LIBDIR}/libTIMING.a
LIBNOMPI  = ${LIBDIR}/libNOMPI.a
LIBEMPIRICALIE   = ${LIBDIR}/libEMPIRICALIE.a
LIBEMPIRICALGM = ${LIBDIR}/libEMPIRICALGM.a

# Libraries should be compiled first, because modules are used in main
${OBJECTS_EXE}: ${LIBSHARE} ${LIBTIMING} ${MY_LIB}   

EXE = ${BINDIR}/rbe.exe

RBE:
	@echo ' '
	@make ${EXE}
	@echo "Program ${EXE} has been brought up to date."
	@echo ' '

${EXE}: ${OBJECTS_EXE}
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ../${MY_LIB}; \
		ar -x ${LIBTIMING}; \
		ar -x ${LIBSHARE}; \
		ar -x ${LIBNOMPI}; \
		ar -x ${LIBEMPIRICALGM};\
		ar -x ${LIBEMPIRICALIE};
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} Tmp_/*.o ${LBLAS} ${Lflag2}
	rm -rf Tmp_

RBE_orig:
	${COMPILE.f90} ${Lflag2} -o ${EXE} rbe_main.f90 ${OBJECTS} \
	-L${LIBDIR} -lTIMING -lSHARE -lNOMPI

distclean: clean
	rm -f Makefile.DEPEND Makefile.RULES
