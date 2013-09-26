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

LIB:
	make cleansrc
	cd ${WSD}/main; rm -f *.o *.a

	cd ${WSD}/general; make SEARCH=
	cd ${WSD}/meshAMR; make SEARCH="${SEARCH}" 
	cd ${WSD}/pic; make SEARCH="${SEARCH}"
	cd ${WSD}/species; make SEARCH="${SEARCH}"
	cd ${WSD}/models/exosphere; make SEARCH="${SEARCH}"
	cd ${WSD}/main; make SEARCH="${SEARCH}"

#	cd srcInterface; make LIB 
#	cd srcInterface; ar -scr amps2swmf.a amps2swmf.o

#	ar -src ../../lib/libPT.a ${WSD}/general/*.o ${WSD}/meshAMR/*.o ${WSD}/pic/*.o ${WSD}/species/*.o ${WSD}/models/exosphere/exosphere.a srcInterface/PT_wrapper.o srcInterface/amps2swmf.o 
#	ar -rc /Users/vtenishe/SWMF/SWMF/lib//libPT.a ${WSD}/general/general.a ${WSD}/meshAMR/mesh.a ${WSD}/main/mainlib.a ${WSD}/pic/amps.a ${WSD}/species/species.a ${WSD}/models/exosphere/exosphere.a srcInterface/amps2swmf.a 
#srcInterface/amps2swmf.o 
#srcInterface/amps2swmf.o

amps:
	@rm -f amps
	make LIB
	cd ${WSD}/main; make amps SEARCH="${SEARCH}"

	${CC} -o ${EXE} ${WSD}/main/main.a ${WSD}/main/mainlib.a ${WSD}/general/general.a \
			${WSD}/pic/amps.a ${WSD}/species/species.a ${WSD}/models/exosphere/exosphere.a \
		 	${WSD}/meshAMR/mesh.a ${WSD}/pic/amps.a ${Lib} ${MPILIB}

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
