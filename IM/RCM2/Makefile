include ../../Makefile.def

install:
	touch src/Makefile.DEPEND
	@(if [ -f src/Makefile.RULES.${OS}${COMPILER} ]; then                \
		cp -f src/Makefile.RULES.${OS}${COMPILER} src/Makefile.RULES;\
	else \
		rm -f src/Makefile.RULES; touch src/Makefile.RULES; \
	fi);

LIB:
	cd src; make LIB

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
	@touch src/Makefile.RULES
	cd src; make clean
	cd src/claw; make clean

distclean: install
	cd src; make distclean
	cd src/claw; make distclean
	rm -f Makefile.def Makefile.conf
