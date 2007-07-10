default : GITM

include Makefile.def

ABDIR   = srcSphereAB
EDDIR   = srcRudy
IEDIR   = ${EMPIRICALIEDIR}
IODIR   = ${DATAREADINDICESDIR}
MAINDIR = src

PLANET=earth

src/Makefile:
	cp src/Makefile.orig src/Makefile

src/ModSize.f90:
	cp src/ModSize.f90.orig src/ModSize.f90

install: src/Makefile src/ModSize.f90
	touch src/Makefile.DEPEND srcInterface/Makefile.DEPEND
	cd src; make DYNAMIC
#
#       General Housekeeping
#

GITM:
	@cd ${SHAREDIR}; make LIB
	@cd $(ABDIR);    make LIB
	@cd $(IEDIR);    make LIB
	@cd $(EDDIR);    make LIB
	@cd $(IODIR);    make LIB
	@cd $(MAINDIR);  make GITM

POST:
	@cd $(MAINDIR);  make POST

GITM2 = ${DIR}/UA/GITM2

LIB:
	cd $(ABDIR)     ; make                                         LIB
	cd $(EDDIR)     ; make LIBPREV=${GITM2}/${ABDIR}/libSphere.a   LIBADD
	cd $(MAINDIR)   ; make LIBPREV=${GITM2}/${EDDIR}/libABIEED.a   LIB
	cd srcInterface ; make LIBPREV=${GITM2}/${MAINDIR}/libUA.a     LIB

nompirun:
	make GITM
	cd run; ./GITM.exe

clean:
	@touch src/Makefile.DEPEND src/Makefile.RULES
	@cd $(ABDIR);    make clean
	@cd $(IEDIR);    make clean
	@cd $(EDDIR);    make clean
	@cd $(IODIR);    make clean
	@cd $(MAINDIR);  make clean
	@cd srcInterface;make clean
	@(if [ -d share ]; then cd share; make clean; fi);

distclean:
	@touch src/Makefile.DEPEND src/Makefile.RULES
	@cd $(ABDIR);    make clean
	@cd $(IEDIR);    make clean
	@cd $(EDDIR);    make clean
	@cd $(IODIR);    make clean
	@cd $(MAINDIR);  make distclean
	@cd srcInterface;make distclean
	rm -f src/Makefile src/ModSize.f90 *~
#
#       Create run directories
#
rundir:
	mkdir -p ${RUNDIR}/UA
	cd ${RUNDIR} ; ln -s ${BINDIR}/PostProcess.exe ./PostGITM.exe
	cd ${RUNDIR}/UA; \
		mkdir restartOUT data; \
		ln -s restartOUT restartIN; \
		ln -s ${BINDIR}/pGITM .; \
		ln -s ${UADIR}/srcData DataIn
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR} ; \
		ln -s ${BINDIR}/GITM.exe . ; \
		cp UA/DataIn/UAM.in . ; \
		touch core ; chmod 444 core ; \
		ln -s UA/* .; \
	fi);

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
	make clean
	./Config.pl -Earth -g=9,9,50,4
	make GITM

TESTDIR = run_test

MPIRUN = mpirun -np 2

test_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES UADIR=`pwd`
	cd ${TESTDIR}; cp UA/DataIn/UAM.in.test.noAPEX UAM.in

test_run:
	cd ${TESTDIR}; ${MPIRUN} ./GITM.exe > runlog

test_check:
	-@(${SCRIPTDIR}/DiffNum.pl -b -r=1e-5 \
		${TESTDIR}/UA/data/log00000002.dat \
		srcData/log00000002.dat.noAPEX > test.diff)
	ls -l test.diff

dist:
	make distclean
	tar cvzf gitm_`date "+%y%m%d"`.tgz Makefile* config.pl get_info.pl \
	    share src srcData srcIDL srcIE srcIO srcInterface \
	    srcMake srcRudy srcSphereAB
	make install
