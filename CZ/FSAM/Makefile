default: xfsam

include Makefile.def

help:
	@echo Makefile targets:
	@echo 
	@echo 'make            - compile xfsam.exe'
	@echo 'make install    - install FSAM'
	@echo 'make LIB        - compile libCZ.a for SWMF'
	@echo 'make rundir     - create run directory'
	@echo 'make test       - test FSAM in stand alone mode'
	@echo 'make clean      - remove object files'
	@echo 'make distclean  - uninstall FSAM'
	@echo

INSTALLFILES =  src/Makefile.DEPEND \
		src/Makefile.RULES \
		srcInterface/Makefile.DEPEND

install: 
	touch ${INSTALLFILES}
	cp srcProblem/ModUserSetup.test.f90 src/ModUserSetup.f90
	cp srcProblem/ModPar.test.f90 src/ModPar.f90

xfsam:
	cd ${SHAREDIR}; $(MAKE) LIB
	cd ${UTILDIR};  $(MAKE) install
	cd src; 	$(MAKE) xfsam

nompirun: xfsam
	cd ${RUNDIR}; ./xfsam.exe

LIB:	
	cd src; $(MAKE) LIB
	cd srcInterface; $(MAKE) LIB

# default problem
PROBS=test

rundir: 
	mkdir -p ${RUNDIR}/CZ
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR}; \
			ln -s ${MYDIR}/src/xfsam.exe .; \
			cp ${MYDIR}/input/PARAM.in.${PROBS} PARAM.in; \
			cp ${MYDIR}/input/gong.l4b.14 .; \
	fi)

#####################################RUN########################################

# run directory for test
TESTDIR = run_test

test: 
	rm -f test_fsam.diff
	@echo "compile..." > test_fsam.diff
	${MAKE} test_compile
	@echo "rundir..." >> test_fsam.diff
	${MAKE} test_rundir
	@echo "run..." >> test_fsam.diff
	${MAKE} test_run
	@echo "check..." >> test_fsam.diff
	${MAKE} test_check

# This should be done with Config.pl -problem=test
# It should not recompile if the files did not change

test_compile:
	rm -f src/xfsam.exe
	cp srcProblem/ModUserSetup.test.f90 src/ModUserSetup.f90;
	cp srcProblem/ModPar.test.f90 src/ModPar.f90;
	$(MAKE) xfsam

test_rundir:
	rm -rf ${TESTDIR}
	${MAKE} rundir PROBS=test RUNDIR=${TESTDIR} STANDALONE="YES"

test_run:
	@if [ "$(MPIRUN)" = "mpiexec" ]; then \
		cd ${TESTDIR}; mpiexec -np 6 xfsam.exe > runlog; \
	else \
		cd ${TESTDIR}; mpirun -np 6 xfsam.exe > runlog; \
	fi;

test_check:
	grep -v 'Grid cell updates' ${TESTDIR}/runlog > ${TESTDIR}/runlog.new
	grep -v 'Grid cell updates' output/runlog.ref > ${TESTDIR}/runlog.ref
	-${SCRIPTDIR}/DiffNum.pl -a=1e-5 -t \
		${TESTDIR}/runlog.new ${TESTDIR}/runlog.ref \
		> test_fsam.diff
	ls -l test_fsam.diff

################################################################################

clean:	
	@touch ${INSTALLFILES}
	cd src; $(MAKE) clean
	cd srcInterface; $(MAKE) clean
	@(if [ -d util  ]; then cd util;  $(MAKE) clean; fi);
	@(if [ -d share ]; then cd share; $(MAKE) clean; fi);

distclean:
	./Config.pl -uninstall

allclean:
	@touch ${INSTALLFILES}
	cd src; $(MAKE) distclean
	cd srcInterface; $(MAKE) distclean
	rm -f *~
