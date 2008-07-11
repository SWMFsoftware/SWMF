default: HEIDI

include Makefile.def

INSTALLFILES =  src/Makefile.DEPEND \
		src/Makefile.RULES \
		srcInterface/Makefile.DEPEND \
		srcInterface/srcIONO/Makefile.DEPEND \
		srcInterface/srcIE/Makefile.DEPEND

install: 
	touch ${INSTALLFILES}

#
#       General Housekeeping
#

HEIDI:  install
	@cd ${NOMPIDIR};	make LIB
	@cd ${SHAREDIR};  	make LIB
	@cd ${TIMINGDIR}; 	make LIB 
	@cd ${EMPIRICALIEDIR};	make LIB
	@cd ${EMPIRICALGMDIR};	make LIB
	@cd src;	        make HEIDI
OUTDIR =output
TESTDIR = run_test

test:
	@echo "There is no test for HEIDI" > notest.diff
#	@echo "test_rundir..." >> test.diff
#	make   test_rundir
#	@echo "test_run..."    >> test.diff
#	make   test_run

test1:
	@echo "test_compile..." > test.diff
	make   test_compile
	@echo "test_rundir..." >> test.diff
	make   test_rundir
	@echo "test_run..."    >> test.diff
	make   test_run
	@echo "test_check..."  >> test.diff
	make   test_check

test_compile:
	make HEIDI

test_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR}
	cp input/*.in ${TESTDIR}
	cp input/*.dat ${TESTDIR}
	cp input/*.glo ${TESTDIR}
	cp input/*.initial ${TESTDIR}
	cp input/*.unff ${TESTDIR}
	cd ${TESTDIR} ;\
	mkdir test1_ionosphere

test_run:
	cd ${TESTDIR}; 
	cp ../../bin/HEIDI.exe ${TESTDIR}/HEIDI.exe; \
	pwd
	cd ${TESTDIR};./HEIDI.exe


test_check:
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/test1_h_prs.002  ${OUTDIR}/test1_h_prs.002\
			> test.diff
	ls -l test.diff



clean: install
	@(if [ -r "Makefile.conf" ]; then  \
	cd src;                      make clean;\
	cd ../srcInterface;             make clean;\
	cd srcIONO;     make clean;\
	cd ../srcIE;       make clean;\
	fi)

distclean: install
	@(if [ -r "Makefile.conf" ]; then \
		cd src;                      make distclean;\
		cd ../srcInterface;             make distclean;\
		cd srcIONO;     make distclean;\
		cd ../srcIE;       make distclean;\
	fi)
	rm -f *~
#
#       Create run directories
#
rundir:
	mkdir -p ${RUNDIR}/IM
	cd ${RUNDIR}/IM; \
		mkdir restartOUT restartIN plots
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR} ; \
		ln -s ${BINDIR}/HEIDI.exe .;\
	fi);
