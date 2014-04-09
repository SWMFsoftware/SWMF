default: xfsam

include Makefile.def

help:
	@echo Makefile targets:
	@echo 
	@echo 'make                          		- compile xfsam.exe'
	@echo 'make LIB 				- compile LibCZ.a for SWMF'
	@echo 'make run 	                        - create run directory'
	@echo 'make test 				- test FSAM in stand alone mode'
	@echo 'make prob PROBS=dynamo RUNDIR=rundir	- convective dynamo run from seed fields'
	@echo 'make prob_rst PROBS=dynamo RUNDIR=rundir - convective dynamo run from a restart file'
	@echo 'make prob PROBS=risetube_qcz RUNDIR=rundir - a run of isolated rising tube in stable cz'
	@echo 'make prob PROBS=benchA RUNDIR=rundir 	- a hydro run for comparison with ASH'
	@echo 'make clean				- remove object files'
	@echo

INSTALLFILES =  src/Makefile.DEPEND \
		src/Makefile.RULES \
		srcInterface/Makefile.DEPEND

install: 
	touch ${INSTALLFILES}

xfsam:
	cd ${SHAREDIR};      make LIB
	cd src;              make xfsam

nompirun: xfsam
	cd ${RUNDIR}; ./xfsam.exe

LIB: 
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

prob:
	@echo "starting..."
	@echo "compiling..."
	${MAKE} compile
	@echo "making rundir..."
	${MAKE} rundir PROBS=${PROBS} RUNDIR=${RUNDIR} STANDALONE="YES" PWDIR=`pwd`
	@echo "running..."
	${MAKE} run RUNDIR=${RUNDIR}

prob_rst:
	@echo "starting..."
	@echo "compiling..."
	${MAKE} compile
	@echo "making rundir..."
	${MAKE} rundir_rst PROBS=${PROBS} RUNDIR=${RUNDIR} STANDALONE="YES" PWDIR=`pwd`
	@echo "running..."
	${MAKE} run	

compile:
	rm -f src/xfsam.exe
	cp src/problems/${PROBS}/ModUserSetup.${PROBS}.f90 src/ModUserSetup.f90; 
	cp src/problems/${PROBS}/ModPar.${PROBS}.f90 src/ModPar.f90; 
	make xfsam

rundir: 
	mkdir -p ${RUNDIR}/CZ
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR}; \
			cp ${CZDIR}/src/xfsam.exe .; \
			cp ${CZDIR}/src/problems/${PROBS}/PARAM.in .; \
			cp ${CZDIR}/input/gong.l4b.14 .; \
			cp ${CZDIR}/src/problems/${PROBS}/jobscript .; \
	fi)

rundir_rst:
	make rundir
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR}; cp ${CZDIR}/src/problems/${PROBS}/PARAM.in.rst PARAM.in;\
		cd ${RUNDIR}; cp ${CZDIR}/src/problems/${PROBS}/*rst.dat .; \
	fi)

run:
	cd ${RUNDIR} ; \
	bsub < jobscript

################################################################################

clean:
	@touch ${INSTALLFILES}
	cd src; make clean
	cd srcInterface; make clean

distclean:
	./Config.pl -uninstall

allclean:
	@touch ${INSTALLFILES}
	cd src; make distclean
	cd srcInterface; make distclean
	rm -f *~

