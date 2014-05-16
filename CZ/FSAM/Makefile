default: xfsam

include Makefile.def

help:
	@echo Makefile targets:
	@echo 
	@echo 'make                          		- compile xfsam.exe'
	@echo 'make LIB 				- compile LibCZ.a for SWMF'
	@echo 'make run 	                        - create run directory'
	@echo 'make test 				- test FSAM in stand alone mode'
	@echo 'make clean				- remove object files'
	@echo

INSTALLFILES =  src/Makefile.DEPEND \
		src/Makefile.RULES \
		srcInterface/Makefile.DEPEND

install: 
	touch ${INSTALLFILES}

xfsam:
	cd ${SHAREDIR};      make LIB
	cd ${UTILDIR};       make install
	cd src; 	     make xfsam

nompirun: xfsam
	cd ${RUNDIR}; ./xfsam.exe

LIB:	
	cd ${UTILDIR};   make install
	cd src; make LIB
	cd srcInterface; make LIB

TESTDIR = run_test

#####################################RUN########################################

test: 
	@echo "starting..." 
	@echo "compiling..."
	${MAKE} compile PROBS=test
	@echo "making rundir..."
	${MAKE} rundir PROBS=test RUNDIR=${TESTDIR} STANDALONE="YES" PWDIR=`pwd`
	@echo "running..."
	${MAKE} run RUNDIR=${TESTDIR}

compile:
	rm -f src/xfsam.exe
	cp srcProblem/ModUserSetup.${PROBS}.f90 src/ModUserSetup.f90;
	cp srcProblem/ModPar.${PROBS}.f90 src/ModPar.f90;
	make xfsam

rundir: 
	mkdir -p ${RUNDIR}/CZ
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR}; \
			cp ${CZDIR}/src/xfsam.exe .; \
			cp ${CZDIR}/input/PARAM.in.${PROBS} PARAM.in; \
			cp ${CZDIR}/input/gong.l4b.14 .; \
			cp ${CZDIR}/input/jobscript.${PROBS} ./jobscript; \
	fi)

run:
	cd ${RUNDIR} ; \
	mpirun -n 6 xfsam.exe > runlog

################################################################################

clean:	
	@touch ${INSTALLFILES}
	cd src; $(MAKE) clean
	cd srcInterface; $(MAKE) clean

distclean: clean
	@touch ${INSTALLFILES}
	cd src; make distclean
	cd srcInterface; make distclean
	rm -f *~

allclean: distclean
