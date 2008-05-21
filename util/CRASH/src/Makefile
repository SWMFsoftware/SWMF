#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND  

OBJECTS=\
	heat_capacity.o \
	get_ionization.o \
        ModIonizPotential.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

MY_LIB = libCRASH.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

LIBSHARE  = ${LIBDIR}/libSHARE.a
LIBTIMING = ${LIBDIR}/libTIMING.a

OBJECTS_EXE= test_ionization.o

${OBJECTS_EXE} : ${MY_LIB} ${LIBSHARE} ${LIBTIMING}

IONIZ : ${OBJECTS_EXE}
	@(if [ ${OS} = "Darwin" ];         then make IONIZ_for_Darwin; \
	else                                    make IONIZ_general; fi)

IONIZ_general:
	${LINK.f90} -o IONIZ.exe ${OBJECTS_EXE} \
		-L. -lCRASH -L${LIBDIR} -lSHARE -lTIMING \
		${Lflag1}

IONIZ_for_Darwin:
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; ar -x ../${MY_LIB}; ar -x ${LIBTIMING}; ar -x ${LIBSHARE}
	${LINK.f90} -o IONIZ.exe ${OBJECTS_EXE} Tmp_/*.o \
		${Lflag1}
	rm -rf Tmp_


distclean: clean
	rm -f Makefile.DEPEND


