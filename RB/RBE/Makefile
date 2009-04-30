default : RBE

include Makefile.def

INSTALLFILES =  src/Makefile.DEPEND \
		src/Makefile.RULES \
		srcInterface/Makefile.DEPEND


install: 
	touch ${INSTALLFILES}

#
#       General Housekeeping
#

RBE:
	@cd ${SHAREDIR};  	make LIB
	@cd ${NOMPIDIR};	make LIB
	@cd ${TIMINGDIR}; 	make LIB 
	@cd ${EMPIRICALIEDIR};	make LIB
	@cd ${EMPIRICALGMDIR};	make LIB
	@cd src;	make LIB
	@cd src;	make RBE

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
	make RBE

test_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES RBDIR=`pwd`

test_run:
	cd ${TESTDIR}; ./rbe.exe > runlog 
#v10 > runlog


test_check:
	gunzip -c output/2002f296_e.fls.standalone.gz > ${TESTDIR}/RB/2002f296_e.fls.ref
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/RB/plots/2002f296_e.fls ${TESTDIR}/RB/2002f296_e.fls.ref \
		> test.diff
	ls -l test.diff

clean:
	@touch ${INSTALLFILES}
	@cd src; make clean
	@cd srcInterface; make clean
	@(if [ -d util ];  then cd util;  make clean; fi);
	@(if [ -d share ]; then cd share; make clean; fi);

distclean: 
	./Config.pl -uninstall

allclean:
	@touch ${INSTALLFILES}
	@cd src; make distclean
	cd srcInterface; make distclean
	rm -f *~

#
#       Create run directories
#
rundir:
	mkdir -p ${RUNDIR}/RB
	@(cd ${RUNDIR}; \
		if [ ! -e "EIE/README" ]; then \
			ln -s ${EMPIRICALIEDIR}/data EIE;\
		fi;)
	cd ${RUNDIR}/RB; \
		mkdir restartOUT restartIN plots
	cp ${RBDIR}/input/rbe_e.fin ${RUNDIR}/RB/
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR} ; \
		ln -s ${BINDIR}/rbe.exe .   ; \
		cp ../input/2002_296.* RB/  ; \
		cp ../input/PARAM.in . ; \
		cp ../input/B_wave_eq.dat RB/   ; \
		cp ../input/Horne_chorus.tgz RB/ ; \
		touch core ; chmod 444 core;\
		cd RB ; \
			tar xfz Horne_chorus.tgz ; \
	fi);
