include Makefile.def

help:
	@echo ' '
	@echo '  You can "make" the following (among others):'
	@echo ' '
	@echo '    RCM         stand alone executable'
	@echo '    rundirSA    run directory for stand alone code'
	@echo ' '

install:
	touch src/Makefile.DEPEND

LIB:
	cd src; make LIB

RCM:
	cd ${SHAREDIR}; make LIB
	cd src; make RCM

test:
	@echo "There is no test for RCM2" > notest.diff

rundir:
	mkdir ${RUNDIR}/IM
	cd ${RUNDIR}/IM; mkdir restartIN restartOUT input plots
	cd ${IMDIR}/input; cp rcm* *.dat dktable ${RUNDIR}/IM/
	cd ${RUNDIR}/IM; \
		mv rcmpcp_inp rcmkp_inp input/; rm -f *_inp*; \
		mv dktable trf.dat elecoef.dat rcmcond rcmcrd* rcmlas1 input/;\
		touch rcm.printout rcm.index
	cp ${SCRIPTDIR}/Preplot.pl ${RUNDIR}/IM/

rundirSA:
	mkdir ${RUNDIR}/IM
	cd ${RUNDIR}/IM; mkdir restartIN restartOUT input plots
	cd ${IMDIR}/input; cp rcm* *.dat dktable ${RUNDIR}/IM/
	cd ${RUNDIR}/IM; mv *_inp* input/
	cd ${RUNDIR}/IM; mv dktable trf.dat elecoef.dat rcmcond rcmcrd* rcmlas1 input/
	cd ${RUNDIR}/IM; touch rcm.printout rcm.index

clean:	install
	cd src; make clean
	cd src/claw; make clean

distclean: 
	./Config.pl -uninstall

allclean: install
	cd src; make distclean
	cd src/claw; make distclean
