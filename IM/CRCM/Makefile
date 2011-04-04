default : CRCM

include Makefile.def

INSTALLFILES =  src/Makefile.DEPEND \
		src/Makefile.RULES \
		srcInterface/Makefile.DEPEND


install: 
	touch ${INSTALLFILES}

#
#       General Housekeeping
#

CRCM:
	@cd ${SHAREDIR};  	make LIB
	@cd ${NOMPIDIR};	make LIB
	@cd ${TIMINGDIR}; 	make LIB 
	@cd ${EMPIRICALIEDIR};	make LIB
	@cd ${EMPIRICALGMDIR};	make LIB
	@cd src;	make LIB
	@cd src;	make CRCM

LIB:
	cd src; make LIB
	cd srcInterface; make LIB

TESTDIR = run_test

test:
	@echo "test_compile..." > test.diff
	make   test_compile
	@echo "test_rundir..." >> test.diff
	make   test_rundir
	@echo "test_run..."    >> test.diff
	make   test_run
	@echo "test_check..."  >> test.diff
	make   test_check


test_compile:
	make CRCM

test_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES IMDIR=`pwd`

test_run:
	cd ${TESTDIR}; ./crcm.exe > runlog 
#v10 > runlog


test_check:
#	gunzip -c data/output/2002f296_e.fls.standalone.gz \
#		> ${TESTDIR}/RB/2002f296_e.fls.ref
#	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
#		${TESTDIR}/RB/plots/2002f296_e.fls \
#		${TESTDIR}/RB/2002f296_e.fls.ref \
#		> test.diff
#	ls -l test.diff

clean:
	@touch ${INSTALLFILES}
	@cd src; make clean
	@cd srcInterface; make clean
	@(if [ -d util ];  then cd util;  make clean; fi);
	@(if [ -d share ]; then cd share; make clean; fi);

distclean: 
	rm -f config.log
	cp src/ModEarthHO.f90 src/ModPlanet.f90
	cp src/ModGrid_default.f90 src/ModGrid.f90
	./Config.pl -uninstall

allclean:
	@touch ${INSTALLFILES}
	@cd src; make distclean
	rm -f config.log
	cd srcInterface; make distclean
	rm -f *~

#
#       Create run directories
#
rundir:
	mkdir -p ${RUNDIR}/IM
	@(cd ${RUNDIR}; \
		if [ ! -e "EIE/README" ]; then \
			ln -s ${EMPIRICALIEDIR}/data EIE;\
		fi;)
	cd ${RUNDIR}/IM; \
		mkdir plots restartIN restartOUT
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR} ; \
		ln -s ${BINDIR}/crcm.exe .   ; \
		touch core ; chmod 444 core;\
	fi);
