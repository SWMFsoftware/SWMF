default : HEIDI

include Makefile.def

INSTALLFILES =  src/Makefile.DEPEND \
		src/Makefile.RULES 

install: 
	touch ${INSTALLFILES}

#
#       General Housekeeping
#

HEIDI:
	@cd ${NOMPIDIR};	make LIB
	@cd ${SHAREDIR};  	make LIB
	@cd ${TIMINGDIR}; 	make LIB 
	@cd ${EMPIRICALIEDIR};	make LIB
	@cd ${EMPIRICALGMDIR};	make LIB
	@cd src;	make LIB
	@cd src;	make HEIDI

LIB:
	cd src; make LIB
	cd srcInterface; make LIB

TESTDIR = run_test

test:
	@echo "There is no test for HEIDI" > notest.diff
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
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES IMDIR=`pwd`

test_run:
	cd ${TESTDIR}; ./HEIDI.exe > runlog 

test_check:
	gunzip -c output/2000f223_e.fls.standalone.gz > ${TESTDIR}/IM/2000f223_e.fls.ref
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/plots/2000f223_e.fls ${TESTDIR}/IM/2000f223_e.fls.ref \
		> test.diff
	ls -l test.diff

clean:
	@touch ${INSTALLFILES}
	@cd src; make clean
	# @cd srcInterface; make clean
	@(if [ -d util ];  then cd util;  make clean; fi);
	@(if [ -d share ]; then cd share; make clean; fi);

distclean: clean
	@touch ${INSTALLFILES}
	#@cd src; make distclean
	#@cd srcInterface; make distclean
	rm -f *~ ${INSTALLFILES}

#
#       Create run directories
#
rundir:
	#mkdir -p ${RUNDIR}/EIE
	mkdir -p ${RUNDIR}/IM
	cd ${RUNDIR}/IM; \
		mkdir restartOUT restartIN plots
	#cp ${EMPIRICALIEDIR}/w2k.dat ${RUNDIR}/EIE/
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR} ; \
		ln -s ${BINDIR}/HEIDI.exe .
		touch core ; chmod 444 core;\
	fi);
