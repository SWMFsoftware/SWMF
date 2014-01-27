SHELL=/bin/sh

DEFAULT_TARGET : amps

# These definitions may be overwritten by Makefile.def
SOURCES=src
WSD=srcTemp
InputFileAMPS=moon.input
#SPICE=nospice

include Makefile.def
include Makefile.conf


#include the local makefile (defined the AMPS' compiling varisbles) if exists  
include Makefile.local


# These definitions are inherited from Makefile.def and Makefile.conf
CC=${COMPILE.mpicxx}
CWD=${MYDIR}

install:
	echo " " > Makefile.local
	-rm .ampsConfig.Settings
	./Config.pl -application=Moon
	@echo "AMPS installed"

distclean:
	-rm Makefile.local .ampsConfig.Settings
	./Config.pl -uninstall

allclean:
	rm -rf main srcTemp *.input* amps Makefile.local .ampsConfig.Settings
	$(MAKE) clean

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
	@(if [ -d ${WSD} ]; then cd ${WSD}/general;         $(MAKE) clean; fi);
	@(if [ -d ${WSD} ]; then cd ${WSD}/meshAMR;         $(MAKE) clean; fi);
	@(if [ -d ${WSD} ]; then cd ${WSD}/pic;             $(MAKE) clean; fi);
	@(if [ -d ${WSD} ]; then cd ${WSD}/species;         $(MAKE) clean; fi);
	@(if [ -d ${WSD} ]; then cd ${WSD}/models/exosphere;$(MAKE) clean; fi);
	@(if [ -d ${WSD} ]; then cd ${WSD}/main;            $(MAKE) clean; fi);
	@(if [ -d srcInterface ]; then cd srcInterface;     $(MAKE) clean; fi);

tar:
	cd ../pic-tower/sources/general; rm -f *.o *.a
	cd ../pic-tower/sources/dsmc; rm -f *.o *.a
	tar -cvf sources.tar sources

#Flags=-O3 -fno-inline -ffinite-math-only  -ftrapping-math -fsignaling-nans -wd383 -wd981 -wd869 -wd1418 -LANG:std -Wall -DMPI_ON -DParticleSolver=dsmc -DCompilationTarget=0
#Flags=-O3  -fasm-blocks -use-asm  -fprefetch-loop-arrays -funroll-loops -unroll=3  -mmmx  -wd383 -wd981 -wd869 -wd1418 -LANG:std -Wall -DMPI_ON -DParticleSolver=dsmc -DCompilationTarget=0

#Flags=-g   -use-asm   -fprefetch-loop-arrays -funroll-loops -unroll=3    -LANG:std -Wall -DMPI_ON -DParticleSolver=dsmc -DCompilationTarget=0
export Flags

#Flags=-O3 -fasm-blocks  -wd383 -wd981 -wd869 -wd1418 -LANG:std -Wall -DMPI_ON -DParticleSolver=dsmc -DCompilationTarget=0

SEARCH=-DMPI_ON -LANG:std -I${CWD}/${WSD}/pic -I${CWD}/${WSD}/main  -I${CWD}/${WSD}/meshAMR -I${CWD}/${WSD}/general -I${CWD}/${WSD}/species -I${CWD}/${WSD}/models/exosphere -I${SPICE}/include -I${CWD}

${WSD}:
	./ampsConfig.pl -input ${InputFileAMPS} -no-compile

${LIB_AMPS}: 
	make ${WSD}
	cd ${WSD}/general; make SEARCH_C=
	cd ${WSD}/meshAMR; make SEARCH_C="${SEARCH}" 
	cd ${WSD}/pic; make SEARCH_C="${SEARCH}"
	cd ${WSD}/species; make SEARCH_C="${SEARCH}"

#compile external modules
	$(foreach src, $(ExternalModules), (cd ${WSD}/$(src); make SEARCH_C="${SEARCH}")) 

	cd ${WSD}/main; make SEARCH_C="${SEARCH}"
	cp -f ${WSD}/main/mainlib.a ${WSD}/libAMPS.a


ifeq ($(SPICE),nospice)
	cd ${WSD}; ${AR} libAMPS.a general/*.o meshAMR/*.o pic/*.o species/*.o 
	$(foreach src, $(ExternalModules), (cd ${WSD}; ${AR} libAMPS.a $(src)/*.o))
else
	rm -rf ${WSD}/tmpSPICE
	mkdir ${WSD}/tmpSPICE
	cp ${SPICE}/lib/cspice.a ${WSD}/tmpSPICE
	cd ${WSD}/tmpSPICE; ar -x cspice.a

	cd ${WSD}; ${AR} libAMPS.a general/*.o meshAMR/*.o pic/*.o species/*.o  tmpSPICE/*.o 
	$(foreach src, $(ExternalModules), (cd ${WSD}; ${AR} libAMPS.a $(src)/*.o))
endif

LIB: ${LIB_AMPS}
	cd srcInterface; make LIB SEARCH_C="${SEARCH}"

amps: ${LIB_AMPS}
	@rm -f amps
	cd ${WSD}/main; make amps SEARCH_C="${SEARCH}"

ifeq ($(SPICE),nospice)
	${CC} -o ${EXE} ${WSD}/main/main.a ${LIB_AMPS} ${Lib} ${MPILIB}  
else 
	${CC} -o ${EXE} ${WSD}/main/main.a ${LIB_AMPS} ${Lib} ${MPILIB} ${SPICE}/lib/cspice.a
endif

TESTDIR = run_test

test:
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
