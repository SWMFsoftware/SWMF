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

OUTDIR   = plots
TESTDIR  = run_test
HEIDIDIR = ${DIR}/IM/HEIDI
CHECKDIR = output

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
	make HEIDI


test_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE="YES" PWDIR=`pwd`
	cd input; cp PARAM.in ../${TESTDIR}
test_run: 
	cd ${TESTDIR}; ${MPIRUN} ./HEIDI.exe

test_check:
	${SCRIPTDIR}/DiffNum.pl -t -r=0.001 -a=1e-10 \
		${TESTDIR}/IM/${OUTDIR}/hydrogen/test1_h_prs.004  ${CHECKDIR}/test1_h_prs.004 \
			> test.diff
	ls -l test.diff

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
	ln -s ${HEIDIDIR}/input/* ${RUNDIR}/IM/input;\
	cp ${HEIDIDIR}/input/*.gz ${RUNDIR}/IM/restartIN;\
	cd  ${RUNDIR}/IM/restartIN;\
	gzip -d *.gz


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

