#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND  

MY_LIB = ${LIBDIR}/libCRASH.a

OBJECTS=\
	heat_capacity.o \
	ModStatSum.o\
	ModAtomicMass.o\
        ModIonizPotential.o 

OBJECTS_EXE= test_ionization.o



DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS} ${OBJECTS_EXE} 



LIB:    DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: DEPEND ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

LIBSHARE  = ${LIBDIR}/libSHARE.a
LIBTIMING = ${LIBDIR}/libTIMING.a


${OBJECTS_EXE} : DEPEND ${MY_LIB} ${LIBSHARE} ${LIBTIMING}

IONIZ : ${OBJECTS_EXE}
	@(if [ ${OS} = "Darwin" ];         then make IONIZ_for_Darwin; \
	else                                    make IONIZ_general; fi)

IONIZ_general:
	${LINK.f90} -o IONIZ.exe ${OBJECTS_EXE} \
		-L${LIBDIR} -lCRASH  -lSHARE -lTIMING \
		${Lflag1}

IONIZ_for_Darwin:
	rm -rf Tmp_; mkdir Tmp_ ;
	cd Tmp_; ar -x ${MY_LIB}; ar -x ${LIBTIMING}; ar -x ${LIBSHARE};
	${LINK.f90} -o IONIZ.exe ${OBJECTS_EXE} Tmp_/*.o \
		${Lflag1}    ; 
	rm -rf Tmp_ ;


distclean: clean 
	rm -f Makefile.DEPEND; \
        rm -rf Tmp_;




