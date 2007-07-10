default : RBE

include Makefile.def

install: 
	touch src/Makefile.DEPEND srcInterface/Makefile.DEPEND

#
#       General Housekeeping
#

RBE:
	@cd ${NOMPIDIR};  make LIB
	@cd ${SHAREDIR};  make LIB
	@cd ${TIMINGDIR}; make LIB
	@cd src;          make LIB
	@cd src;          make RBE

LIB:
	cd src; make LIB
	cd srcInterface; make LIB

TESTDIR = run_test

test:
	@echo "test_compile..." > test.diff
	make test_compile
	@echo "test_rundir..." >> test.diff
	make test_rundir
	@echo "test_run..." >> test.diff
	make test_run
	@echo "test_check..." >> test.diff
	make test_check

test_compile:
	make RBE

test_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES RBDIR=`pwd`

test_run:
	cd ${TESTDIR}; ./rbe.exe > runlog

test_check:
	gunzip -c output/2000f223_e.fls.standalone.gz > ${TESTDIR}/2000f223_e.fls.ref
	${SCRIPTDIR}/DiffNum.pl -r=0.001 -a=1e-10 \
		${TESTDIR}/2000f223_e.fls ${TESTDIR}/2000f223_e.fls.ref \
		> test.diff
	ls -l test.diff

clean:
	@touch src/Makefile.DEPEND src/Makefile.RULES
	@cd src; make clean
	@(if [ -d util ];  then cd util;  make clean; fi);
	@(if [ -d share ]; then cd share; make clean; fi);

distclean: clean
	@cd src; make distclean
	rm -f *~

#
#       Create run directories
#
run:
	make rundir

rundir:
	mkdir -p ${RUNDIR}/RB
	cd ${RUNDIR}/RB; \
		mkdir restartOUT output; \
		ln -s restartOUT restartIN;
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR} ; \
		ln -s ${BINDIR}/rbe.exe .   ; \
		ln -s ../input/2000_223.* . ; \
		ln -s ../input/w2k.dat .    ;\
		ln -s ../input/rbe_e.fin .  ; \
		cp ../input/rbe_swmf.dat .  ; \
		touch core ; chmod 444 core ; \
		ln -s RB/* .; \
	fi);

