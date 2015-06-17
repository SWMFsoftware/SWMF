SHELL=/bin/sh

DEFAULT_TARGET : amps

# These definitions may be overwritten by Makefile.def
SOURCES=src
WSD=srcTemp
InputFileAMPS=moon.input
SPICE=nospice

#Compiling with the CCMC's Kameleon
KAMELEON=nokameleon
BOOST=noboost
BATL=nobatl
SWMF=noswmf
TESTMODE=off

include Makefile.def
include Makefile.conf


#include the local makefile (defined the AMPS' compiling variables)  
include Makefile.local

#the default value of the c++ compiler flags
SEARCH_C=-DMPI_ON -LANG:std -I${CWD}/${WSD}/pic -I${CWD}/${WSD}/main  -I${CWD}/${WSD}/meshAMR -I${CWD}/${WSD}/general -I${CWD}/${WSD}/models/electron_impact -I${CWD}/${WSD}/models/sputtering -I${CWD}/${WSD}/models/charge_exchange -I${CWD}/${WSD}/models/photolytic_reactions -I${CWD}/${WSD}/species -I${CWD}/${WSD}/models/exosphere -I${SPICE}/include -I${BOOST}/include -I${KAMELEON}/src -I${CWD}

#the additional argument string for the fortran compiler
SEARCH_F= 

# These definitions are inherited from Makefile.def and Makefile.conf
CC=${COMPILE.mpicxx}
CWD=${MYDIR}
AMPSLINKER=${CC}


	AMPSLINKLIB= 
	
ifneq ($(BATL),nobatl)	
	AMPSLINKLIB+=${BATL}/lib/libREADAMR.a
	AMPSLINKER=${LINK.f90}

ifeq ($(COMPILE.mpicxx),gfortran)  
	SEARCH_F+= -J${BATL}/share/include 
else ifeq ($(COMPILE.mpicxx),mpif90)
	SEARCH_F+= -I${BATL}/share/include
else ifeq  ($(COMPILE.mpicxx),ifort)
	SEARCH_F+= -module ${BATL}/share/include
endif


endif
	
ifneq ($(SWMF),noswmf)	
	AMPSLINKLIB+=${SWMF}/lib/*
	AMPSLINKER=${LINK.f90}
endif

ifneq ($(KAMELEON),nokameleon)
	AMPSLINKLIB+=${KAMELEON}/lib/ccmc/libccmc.dylib 
endif

ifneq ($(SPICE),nospice)
	AMPSLINKLIB+=${SPICE}/lib/cspice.a
endif

ifeq ($(TESTMODE),on)  
	SEARCH_C+=-D _PIC_NIGHTLY_TEST_MODE_=_PIC_MODE_ON_
endif


#SPICE=/Users/vtenishe/SPICE/Toolkit/cspice/include

install:
	@echo " " > Makefile.local
	@rm -f .ampsConfig.Settings
	@./Config.pl -application=Moon
	@echo "AMPS installed"

distclean:
	./Config.pl -uninstall

allclean: clean
	rm -rf main srcTemp *.input* amps Makefile.local .ampsConfig.Settings

rundir:
	mkdir -p ${RUNDIR}/PT
	cd ${RUNDIR}/PT; mkdir restartIN restartOUT plots


# /Users/vtenishe/Debugger/eclipse-workspace/pic-input-preprocess

EXE=amps

#Lib=  -lm

#MPIRUN=mpirun -np 4
#RUNDIR=run

LIB_AMPS = ${WSD}/libAMPS.a

clean:
	rm -f ${LIB_AMPS}
	@(if [ -d ${WSD} ]; then cd ${WSD}/general;                    $(MAKE) clean; fi);
	@(if [ -d ${WSD} ]; then cd ${WSD}/meshAMR;                    $(MAKE) clean; fi);
	@(if [ -d ${WSD} ]; then cd ${WSD}/pic;                        $(MAKE) clean; fi);
	@(if [ -d ${WSD} ]; then cd ${WSD}/species;                    $(MAKE) clean; fi);
	@(if [ -d ${WSD} ]; then cd ${WSD}/models/exosphere;           $(MAKE) clean; fi);
	@(if [ -d ${WSD} ]; then cd ${WSD}/main;                       $(MAKE) clean; fi);
	@(if [ -d srcInterface ]; then cd srcInterface;                $(MAKE) clean; fi);
	@(if [ -d ${WSD} ]; then cd ${WSD}/models/sputtering;          $(MAKE) clean; fi);
	@(if [ -d ${WSD} ]; then cd ${WSD}/models/charge_exchange;     $(MAKE) clean; fi);
	@(if [ -d ${WSD} ]; then cd ${WSD}/models/electron_impact;     $(MAKE) clean; fi);
	@(if [ -d ${WSD} ]; then cd ${WSD}/models/photolytic_reactions;$(MAKE) clean; fi);

tar:
	cd ../pic-tower/sources/general; rm -f *.o *.a
	cd ../pic-tower/sources/dsmc; rm -f *.o *.a
	tar -cvf sources.tar sources

${WSD}:
	./ampsConfig.pl -input ${InputFileAMPS} -no-compile

${LIB_AMPS}: 
	(if [ -d ${WSD} ]; then rm -rf ${WSD}; fi);
	make ${WSD}
	cd ${WSD}/general;                     make SEARCH_C=
	cd ${WSD}/meshAMR;                     make SEARCH_C="${SEARCH_C}" 
	cd ${WSD}/pic;                         make SEARCH_C="${SEARCH_C}" SEARCH="${SEARCH_F}" 
	cd ${WSD}/species;                     make SEARCH_C="${SEARCH_C}"
	cd ${WSD}/models/electron_impact;      make SEARCH_C="${SEARCH_C}"
	cd ${WSD}/models/sputtering;           make SEARCH_C="${SEARCH_C}"
	cd ${WSD}/models/charge_exchange;      make SEARCH_C="${SEARCH_C}"
	cd ${WSD}/models/photolytic_reactions; make SEARCH_C="${SEARCH_C}" 

#compile external modules
	$(foreach src, $(ExternalModules), (cd ${WSD}/$(src); make SEARCH_C="${SEARCH_C}")) 
	cd ${WSD}/main; make SEARCH_C="${SEARCH_C}"
	cp -f ${WSD}/main/mainlib.a ${WSD}/libAMPS.a
ifeq ($(SPICE),nospice)
	cd ${WSD}; ${AR} libAMPS.a general/*.o meshAMR/*.o pic/*.o species/*.o models/electron_impact/*.o models/sputtering/*.o models/charge_exchange/*.o models/photolytic_reactions/*.o
	$(foreach src, $(ExternalModules), (cd ${WSD}; ${AR} libAMPS.a $(src)/*.o))
else
	rm -rf ${WSD}/tmpSPICE
	mkdir ${WSD}/tmpSPICE
	cp ${SPICE}/lib/cspice.a ${WSD}/tmpSPICE
	cd ${WSD}/tmpSPICE; ar -x cspice.a

	cd ${WSD}; ${AR} libAMPS.a general/*.o meshAMR/*.o pic/*.o species/*.o models/electron_impact/*.o models/sputtering/*.o models/charge_exchange/*.o models/photolytic_reactions/*.o tmpSPICE/*.o 
	$(foreach src, $(ExternalModules), (cd ${WSD}; ${AR} libAMPS.a $(src)/*.o))
endif

LIB: 
	(if [ -d ${WSD} ]; then rm -rf ${WSD}; fi);
	make ${LIB_AMPS} 
	cd srcInterface; make LIB SEARCH_C="${SEARCH_C}"

amps_link:
	${AMPSLINKER} -o amps srcTemp/main/main.a srcTemp/libAMPS.a \
		${CPPLIB} ${AMPSLINKLIB}

amps: ${LIB_AMPS}
	@rm -f amps
	cd ${WSD}/main; make amps SEARCH_C="${SEARCH_C}"
	make amps_link


TESTDIR = run_test

test:
	(if [ -d ${WSD} ]; then rm -rf ${WSD}; fi); 	
	./Config.pl -application=Moon -amps-test=on 
	rm -f *.diff
	-@($(MAKE) test_amps)
	#@ls -l *.diff

test_amps:
	@echo "test_amps_compile..." > test_amps.diff
	$(MAKE) test_amps_compile
	@echo "test_amps_rundir..." >> test_amps.diff
	$(MAKE) test_amps_rundir
	@echo "test_amps_run..." >> test_amps.diff
	$(MAKE) test_amps_run
	@echo "test_amps_check..." >> test_amps.diff
	$(MAKE) test_amps_check
	#if([ "${KEEP}" ]); then rm -rf run_$@; mv ${TESTDIR} run_$@; fi

test_amps_compile:
	rm -rf ${TESTDIR}
	./ampsConfig.pl -no-compile 
	$(MAKE) amps

test_amps_rundir:
	rm -rf ${TESTDIR}
	mkdir -p ${TESTDIR}
	mv amps ${TESTDIR}

test_amps_run:
	cd ${TESTDIR}; ${MPIRUN} ./amps

test_amps_check:
	-(${SCRIPTDIR}/DiffNum.pl ${TESTDIR}/PT/plots/amps.dat \
	output/test_amps.ref_np`ls ${TESTDIR}/PT/thread* | wc -l | tr -d ' '` \
	> test_amps.diff)
	@ls -l test_amps.diff



t:
	@(cd ${RUNDIR}; perl ${SCRIPTDIR}/DiffNum.pl amps.test.dat ../amps.test.reference.dat > amps.diff)
