include ../Makefile.def
include ../Makefile.conf
# include Makefile.DEPEND
# include Makefile.RULES

OBJECTS = \
	rbe_subs.o\
	pbo_2.o\
	t04_s.o\
	t96_01.o\
	w2k.o\
	geopack.o

DEPEND:
	@${SCRIPTDIR}/depend.pl ${SEARCH} rbe_main.o ${OBJECTS}

RBE:
	@make rbe.exe
	@echo "rbe.exe has been created"

rbe_main.o: rbe_subs.o ${LIBDIR}/libSHARE.a ${LIBDIR}/libTIMING.a

.f.o:
	${COMPILE.f77} ${CFLAGS} ${PRECISION} $<

.f90.o:
	${COMPILE.f90} ${CFLAGS} ${PRECISION} $<

rbe.exe: rbe_main.o ${OBJECTS}
	${COMPILE.f90} ${Lflag2} -o rbe.exe rbe_main.f90 ${OBJECTS} \
	-L${LIBDIR} -lTIMING -lSHARE

MY_LIB = ${LIBDIR}/libRB.a

LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS} RB_wrapper.o
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS} RB_wrapper.o

distclean: clean
	rm -f Makefile.DEPEND Makefile.RULES
