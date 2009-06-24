#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND  

MY_LIB = ${LIBDIR}/libCRASH.a

OBJECTS=\
	ModAtomicMass.o\
	ModEos.o\
	ModPowerLawEos.o\
	ModFermiGas.o\
        ModIonizPotential.o \
	ModPolyimide.o \
	ModStatSum.o \
	ModTransport.o


OBJECTS_EXE= test_ionization.o test_godunov.o test_Fermi.o



DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS} ${OBJECTS_EXE} 



LIB:    DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

LIBSHARE  = ${LIBDIR}/libSHARE.a
LIBTIMING = ${LIBDIR}/libTIMING.a


test_ionization.o : ${LIBSHARE} ${LIBTIMING} ${MY_LIB}

IONIZ:  DEPEND 
	make IONIZ.exe

IONIZ.exe: test_ionization.o
	@(if [ ${OS} = "Darwin" ];         then make IONIZ_for_Darwin; \
	else                                    make IONIZ_general; fi)

IONIZ_general:
	${LINK.f90} -o IONIZ.exe test_ionization.o \
		-L${LIBDIR} -lCRASH  -lSHARE -lTIMING \
		${Lflag1}

IONIZ_for_Darwin:
	rm -rf Tmp_; mkdir Tmp_ ;
	cd Tmp_; ar -x ${MY_LIB}; ar -x ${LIBTIMING}; ar -x ${LIBSHARE};
	${LINK.f90} -o IONIZ.exe test_ionization.o Tmp_/*.o \
		${Lflag1}    ; 
	rm -rf Tmp_ 

test_godunov.o : ${LIBSHARE} ${LIBTIMING} ${MY_LIB}

GODUNOV: DEPEND 
	 make GODUNOV.exe

GODUNOV.exe: test_godunov.o
	@(if [ ${OS} = "Darwin" ];         then make GODUNOV_for_Darwin; \
	else                                    make GODUNOV_general; fi)

GODUNOV_general:
	${LINK.f90} -o GODUNOV.exe test_godunov.o \
		-L${LIBDIR} -lCRASH  -lSHARE -lTIMING \
		${Lflag1}

GODUNOV_for_Darwin:
	rm -rf Tmp_; mkdir Tmp_ ;
	cd Tmp_; ar -x ${MY_LIB}; ar -x ${LIBTIMING}; ar -x ${LIBSHARE};
	${LINK.f90} -o GODUNOV.exe test_godunov.o Tmp_/*.o \
		${Lflag1}    ; 
	rm -rf Tmp_ 

test_eos: IONIZ.exe
	  IONIZ.exe;
	  rm -f *.diff test_re*;
	  cat ../doc/Table*.tex >test_result;
	  cat ../doc/Tex/ionization_table*.tex>test_ref;
	  diff test_result test_ref > test_eos.diff;
	  ls -l *.diff; rm -f test_re*

distclean: clean 
	rm -f Makefile.DEPEND; \
        rm -rf Tmp_;




test_Fermi.o : ${LIBSHARE} ${LIBTIMING} ${MY_LIB}

test_Fermi: DEPEND
	 make test_Fermi.exe

test_Fermi.exe: test_Fermi.o
	@(if [ ${OS} = "Darwin" ];         then make test_Fermi_for_Darwin; \
	else                                    make test_Fermi_general; fi)

test_Fermi_general:
	${LINK.f90} -o test_Fermi.exe test_Fermi.o \
		-L${LIBDIR} -lCRASH  -lSHARE -lTIMING \
		${Lflag1}

test_Fermi_for_Darwin:
	rm -rf Tmp_; mkdir Tmp_ ;
	cd Tmp_; ar -x ${MY_LIB}; ar -x ${LIBTIMING}; ar -x ${LIBSHARE};
	${LINK.f90} -o test_Fermi.exe test_Fermi.o Tmp_/*.o \
		${Lflag1}    ; 
	rm -rf Tmp_ 

