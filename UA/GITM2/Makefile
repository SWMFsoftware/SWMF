
default : GITM

include Makefile.def

ABDIR   = srcSphereAB
EIEDIR  = ${EMPIRICALIEDIR}
EUADIR  = ${EMPIRICALUADIR}
IODIR   = ${DATAREADINDICESDIR}
MAINDIR = src
GLDIR   = srcGlow

PLANET=earth

src/ModSize.f90:
	cp src/ModSize.f90.orig src/ModSize.f90

INSTALLFILES =  src/Makefile.DEPEND \
		src/Makefile.RULES \
		srcInterface/Makefile.DEPEND

install: src/ModSize.f90
	touch ${INSTALLFILES}
	cd src; make DYNAMIC
#
#       General Housekeeping
#

GITM:
	@cd ${SHAREDIR}; make LIB
	@cd $(ABDIR);    make LIB
	@cd $(EIEDIR);   make LIB
	@cd ${EUADIR};   make LIB
	@cd $(IODIR);    make LIB
	@cd $(GLDIR);	 make LIB
	@cd $(MAINDIR);  make GITM

POST:
	@cd $(MAINDIR);  make POST

GITM2 = ${DIR}/UA/GITM2

LIB:
	cd $(ABDIR)     ; make                                         LIB
	cd $(GLDIR)     ; make LIBPREV=${GITM2}/${IODIR}/libABIEEDIO.a LIBADD
	cd $(MAINDIR)   ; make LIBPREV=${GITM2}/${ABDIR}/libSphere.a   LIBADD
	cd srcInterface ; make LIBPREV=${GITM2}/${MAINDIR}/libUA.a     LIB

nompirun:
	make GITM
	cd ${RUNDIR}; ./GITM.exe

clean:
	@touch ${INSTALLFILES}
	@cd $(ABDIR);    make clean
	@cd $(MAINDIR);  make clean
	@cd $(GLDIR);    make clean
	@cd srcInterface;make clean
	@(if [ -d share ]; then cd share; make clean; fi);
	@(if [ -d util ];  then cd util;  make clean; fi);

distclean: 
	./Config.pl -uninstall

allclean:
	@touch ${INSTALLFILES}
	@cd $(ABDIR);    make clean
	@cd $(MAINDIR);  make distclean
	@cd srcInterface;make distclean
	rm -f *~ srcData/UAM.in
#
#       Create run directories
#
rundir:
	mkdir -p ${RUNDIR}/UA
	@(cd ${RUNDIR}; \
		if [ ! -e "EIE/README" ]; then \
			ln -s ${EMPIRICALIEDIR}/data EIE;\
		fi;)
	cd ${RUNDIR}; ln -s ${BINDIR}/PostProcess.exe ./PostGITM.exe
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


TESTDIR = run_test

MPIRUN = mpirun -np 2

test:
	-@(make test_earth)
	-@(make test_mars)
	ls -l *.diff

test_earth:
	@echo "test_earth_compile..." > test_earth.diff
	make test_earth_compile
	@echo "test_earth_rundir..." >> test_earth.diff
	make test_earth_rundir
	@echo "test_earth_run..." >> test_earth.diff
	make test_earth_run
	@echo "test_earth_check..." >> test_earth.diff
	make test_earth_check

test_earth_compile:
	./Config.pl -Earth
	./Config.pl -g=9,9,50,4
	make GITM

test_earth_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES UADIR=`pwd`
	cd ${TESTDIR}; cp UA/DataIn/UAM.in.test.noAPEX UAM.in

test_earth_run:
	cd ${TESTDIR}; ${MPIRUN} ./GITM.exe > runlog

test_earth_check:
	-(${SCRIPTDIR}/DiffNum.pl -b -r=1e-5 \
		${TESTDIR}/UA/data/log00000002.dat \
		srcData/log00000002.dat.noAPEX >& test_earth.diff)
	ls -l test_earth.diff

test_mars:
	@echo "test_mars_compile..." > test_mars.diff
	make test_mars_compile
	@echo "test_mars_rundir..." >> test_mars.diff
	make test_mars_rundir
	@echo "test_mars_run..." >> test_mars.diff
	make test_mars_run
	@echo "test_mars_check..." >> test_mars.diff
	make test_mars_check

test_mars_compile:
	./Config.pl -Mars
	./Config.pl -g=1,1,90,1
	make GITM

test_mars_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES UADIR=`pwd`
	cd ${TESTDIR}; cp UA/DataIn/UAM.in.Mars UAM.in

test_mars_run:
	cd ${TESTDIR}; ${MPIRUN} ./GITM.exe > runlog

test_mars_check:
	-(${SCRIPTDIR}/DiffNum.pl -b -r=1e-5 \
		${TESTDIR}/UA/data/log00000002.dat \
		srcData/log00000002.dat.Mars >& test_mars.diff)
	ls -l test_mars.diff

dist:
	make distclean
	tar cvzf gitm_`date "+%y%m%d"`.tgz Makefile* config.pl get_info.pl \
	    share src srcData srcIDL srcIE srcIO srcInterface \
	    srcMake srcSphereAB
	make install
