default: HEIDI

include Makefile.def

INSTALLFILES =  src/Makefile.DEPEND \
		src/Makefile.RULES \
		srcInterface/Makefile.DEPEND

install: 
	touch ${INSTALLFILES}


HEIDI:  install
	@cd ${SHAREDIR};  	make LIB
	@cd ${NOMPIDIR};	make LIB
	@cd ${TIMINGDIR}; 	make LIB 
	@cd ${EMPIRICALIEDIR};	make LIB
	@cd ${EMPIRICALGMDIR};	make LIB
	@cd src;	        make HEIDI

LIB:  install
	@cd ${SHAREDIR};  	make LIB
	@cd ${NOMPIDIR};	make LIB
	@cd ${TIMINGDIR}; 	make LIB 
	@cd ${EMPIRICALIEDIR};	make LIB
	@cd ${EMPIRICALGMDIR};	make LIB
	@cd src;	        make LIB
	@cd srcInterface;	make LIB

OUTDIR    = plots
TESTDIR1  = run_test_analytic
TESTDIR2  = run_test_numeric
HEIDIDIR  = ${DIR}/IM/HEIDI
CHECKDIR  = output

test:
	@rm -f *.diff
	-@(make test_analytic)	
	-@(make test_numeric)


test_analytic:
	@echo "test_compile..." > test_analytic.diff
	make   test_compile
	@echo "test_rundir..." >> test_analytic.diff
	make   test_analytic_rundir
	@echo "test_run..."    >> test_analytic.diff
	make   test_analytic_run
	@echo "test_check..."  >> test_analytic.diff
	make   test_analytic_check

test_numeric:
	@echo "test_compile..." > test_numeric.diff
	make   test_compile
	@echo "test_rundir..." >> test_numeric.diff
	make   test_numeric_rundir
	@echo "test_run..."    >> test_numeric.diff
	make   test_numeric_run
	@echo "test_check..."  >> test_numeric.diff
	make   test_numeric_check


test_compile:
	make HEIDI


test_analytic_rundir:
	rm -rf ${TESTDIR1}
	make rundir RUNDIR=${TESTDIR1} STANDALONE="YES" PWDIR=`pwd`
	cd input; cp PARAM.analytic.in ../${TESTDIR1}/PARAM.in

test_numeric_rundir:
	rm -rf ${TESTDIR2}
	make rundir RUNDIR=${TESTDIR2} STANDALONE="YES" PWDIR=`pwd`
	cd input; cp PARAM.numeric.in ../${TESTDIR2}/PARAM.in

test_analytic_run: 
	cd ${TESTDIR1}; ${MPIRUN} ./HEIDI.exe

test_numeric_run: 
	cd ${TESTDIR2}; ${MPIRUN} ./HEIDI.exe

test_analytic_check:
	${SCRIPTDIR}/DiffNum.pl -t -r=0.001 -a=1e-10 \
		${TESTDIR1}/IM/${OUTDIR}/hydrogen/test1_h_prs.004  ${CHECKDIR}/test1_h_prs_analytic.004 \
			> test_analytic.diff
	ls -l test_analytic.diff

test_numeric_check:
	${SCRIPTDIR}/DiffNum.pl -t -r=0.001 -a=1e-10 \
		${TESTDIR2}/IM/${OUTDIR}/hydrogen/test1_h_prs.004  ${CHECKDIR}/test1_h_prs_numeric.004 \
			> test_numeric.diff
	ls -l test_numeric.diff



rundir:
	mkdir -p ${RUNDIR}/IM
	@(cd ${RUNDIR}; \
		if [ ! -e "EIE/README" ]; then \
			ln -s ${EMPIRICALIEDIR}/data EIE;\
		fi;)
	cd ${RUNDIR}/IM; \
		mkdir input plots restartIN restartOUT
	cd ${RUNDIR}/IM/plots; \
		mkdir electron hydrogen helium oxygen ionosphere
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR} ; \
		ln -s ${BINDIR}/HEIDI.exe .;\
	fi);
	ln -s ${HEIDIDIR}/input/* ${RUNDIR}/IM/input
	cp ${HEIDIDIR}/data/input/*.gz ${RUNDIR}/IM/restartIN
	cd ${RUNDIR}/IM/restartIN; gunzip *.gz

clean: install
	@(if [ -r "Makefile.conf" ]; then  \
	cd src;                      make clean;\
	cd ../srcInterface;             make clean;\
	fi)

distclean: 
	./Config.pl -uninstall

allclean: install
	@(if [ -r "Makefile.conf" ]; then \
		cd src;                      make distclean;\
		cd ../srcInterface;          make distclean;\
	fi)
	rm -f *~

