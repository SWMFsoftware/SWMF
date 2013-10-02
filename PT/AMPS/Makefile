SHELL=/bin/sh

DEFAULT_TARGET : LIB

# These definitions may be overwritten by Makefile.def
SOURCES=src
WSD=srcTemp
#SPICE=nospice

include Makefile.def
include Makefile.conf

# These definitions are inherited from Makefile.def and Makefile.conf
CC=${COMPILE.mpicxx}
CWD=${MYDIR}

install:
	Config.pl -application=Moon
	@echo "AMPS installed"

distclean:
	./Config.pl -uninstall

allclean: cleansrc
	rm -rf main srcTemp *.input* amps

rundir:
	mkdir -p ${RUNDIR}/PT
	cd ${RUNDIR}/PT; mkdir restartIN restartOUT plots


# /Users/vtenishe/Debugger/eclipse-workspace/pic-input-preprocess

EXE=amps

#Lib=  -lm

#MPIRUN=mpirun -np 4
#RUNDIR=run

cleansrc:
	cd ${WSD}/general; make clean
	cd ${WSD}/meshAMR; make clean
	cd ${WSD}/pic; make clean
	cd ${WSD}/species; make clean
	cd ${WSD}/models/exosphere; make clean
	cd ${WSD}/main; make clean

#	cd srcInterface; rm -f *.o

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
	./ampsConfig.pl -no-compile

LIB_AMPS = srcTemp/libAMPS.a

${LIB_AMPS}: 
	make ${WSD}
	cd ${WSD}/general; make SEARCH_C=
	cd ${WSD}/meshAMR; make SEARCH_C="${SEARCH}" 
	cd ${WSD}/pic; make SEARCH_C="${SEARCH}"
	cd ${WSD}/species; make SEARCH_C="${SEARCH}"
	cd ${WSD}/models/exosphere; make SEARCH_C="${SEARCH}"
	cd ${WSD}/main; make SEARCH_C="${SEARCH}"
	cp -f ${WSD}/main/mainlib.a ${WSD}/libAMPS.a
	cd ${WSD}; ${AR} libAMPS.a general/*.o meshAMR/*.o pic/*.o \
		species/*.o models/exosphere/*.o main/main_lib.o

LIB: ${LIB_AMPS}
	cd srcInterface; make LIB SEARCH_C="${SEARCH}"

amps: ${LIB_AMPS}
	@rm -f amps
	cd ${WSD}/main; make amps SEARCH_C="${SEARCH}"
	${CC} -o ${EXE} ${WSD}/main/main.a ${LIB_AMPS} ${Lib} ${MPILIB}

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
