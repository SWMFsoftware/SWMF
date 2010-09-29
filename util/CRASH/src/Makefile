#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND
include Makefile.RULES  


MY_LIB = ${LIBDIR}/libCRASH.a

#DEBUG = ${DEBUGFLAG}

IOMIX=\
	ModIonMix.o \
	ModOpacityVoigt.o \
	ModMultiGroup.o	


OBJECTS=\
	ModAtomicMass.o \
	ModAtomicDataMix.o \
	ModEos.o\
	ModExcitation.o \
	ModExcitationData.o \
	ModPowerLawEos.o\
	ModFermiGas.o\
	ModHydrogenicModel.o \
        ModIonizPotential.o \
	ModMolecularData.o \
	ModPartition.o \
	ModPlastic.o \
	ModStatSum.o \
	ModTransport.o \
	${IOMIX}

MPIRUN = 

OBJECTS_EXE= test_ionization.o test_godunov.o test_Fermi.o test_absorption.o \
	     save_eos_table.o



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

test_absorption.o : ${LIBSHARE} ${LIBTIMING} ${MY_LIB}

ABS:  DEPEND 
	make ABS.exe

ABS.exe: test_absorption.o
	@(if [ ${OS} = "Darwin" ];         then make ABS_for_Darwin; \
	else                                    make ABS_general; fi)

ABS_general:
	${LINK.f90} -o ABS.exe test_absorption.o \
		-L${LIBDIR} -lCRASH  -lSHARE -lTIMING \
		${Lflag1}

ABS_for_Darwin:
	rm -rf Tmp_; mkdir Tmp_ ;
	cd Tmp_; ar -x ${MY_LIB}; ar -x ${LIBTIMING}; ar -x ${LIBSHARE};
	${LINK.f90} -o ABS.exe test_absorption.o Tmp_/*.o \
		${Lflag1}    ; 
	rm -rf Tmp_ 

opacity.o : ${LIBSHARE} ${LIBTIMING} ${MY_LIB}

OPAC:  DEPEND 
	make OPAC.exe
	OPAC.exe
	tar -czvf opacities.tgz Xe_*.dat Be_*.dat Pl_*.dat
	rm -f Xe_*.dat Be_*.dat Pl_*.dat

OPAC.exe: opacity.o
	@(if [ ${OS} = "Darwin" ];         then make OPAC_for_Darwin; \
	else                                    make OPAC_general; fi)

OPAC_general:
	${LINK.f90} -o OPAC.exe opacity.o \
		-L${LIBDIR} -lCRASH  -lSHARE -lTIMING \
		${Lflag1}

OPAC_for_Darwin:
	rm -rf Tmp_; mkdir Tmp_ ;
	cd Tmp_; ar -x ${MY_LIB}; ar -x ${LIBTIMING}; ar -x ${LIBSHARE};
	${LINK.f90} -o OPAC.exe opacity.o Tmp_/*.o \
		${Lflag1}    ; 
	rm -rf Tmp_ 

save_eos_table.o : ${LIBSHARE} ${LIBTIMING} ${MY_LIB}

EOS:  DEPEND 
	make EOS.exe
	${MPIRUN} EOS.exe

EOS.exe: save_eos_table.o
	@(if [ ${OS} = "Darwin" ];         then make EOS_for_Darwin; \
	else                                    make EOS_general; fi)

EOS_general:
	${LINK.f90} -o EOS.exe save_eos_table.o \
		-L${LIBDIR} -lCRASH  -lSHARE -lTIMING \
		${Lflag1}

EOS_for_Darwin:
	rm -rf Tmp_; mkdir Tmp_ ;
	cd Tmp_; ar -x ${MY_LIB}; ar -x ${LIBTIMING}; ar -x ${LIBSHARE};
	${LINK.f90} -o EOS.exe save_eos_table.o Tmp_/*.o \
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

test_eos: 
	  make IONIZ
	  IONIZ.exe
	  rm -f test_eos.diff test_re*
	  cat ../doc/Table*.tex >test_result
	  cat ../doc/Tex/ionization_table*.tex>test_ref
	  ${SCRIPTDIR}/DiffNum.pl test_result test_ref > test_eos.diff
	  ls -l *.diff
	  rm -f test_re*

test_abs: 
	  make ABS
	  ABS.exe
	  rm -f test_abs.diff test_re*
	  cat ../doc/polyimide_absorption.dat ../doc/excited_levels.tex>test_result
	   cat ../doc/Tex/polyimide_absorption.dat ../doc/Tex/excited_levels.tex>test_ref
	  ${SCRIPTDIR}/DiffNum.pl -t -r=1.0e-8 -a=1.0e-13 test_result test_ref > test_abs.diff
	  ls -l *.diff
	  rm -f test_re*

distclean: clean 
	rm -f Makefile.DEPEND Makefile.RULES; \
        rm -rf Tmp_;
	rm -f test*diff *report*



test_Fermi.o : ${LIBSHARE} ${LIBTIMING} ${MY_LIB}

test_Fermi: DEPEND
	 make test_Fermi.exe
	 test_Fermi.exe

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
convert.o : ${LIBSHARE} ${LIBTIMING} ${MY_LIB}

CONVERT: convert.o
	rm -rf Tmp_; mkdir Tmp_ ;
	cd Tmp_; ar -x ${MY_LIB}; ar -x ${LIBTIMING}; ar -x ${LIBSHARE};
	${LINK.f90} -o CONVERT.exe convert.o Tmp_/*.o \
		${Lflag1}    ; 
	rm -rf Tmp_ 

