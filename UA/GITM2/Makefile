default : GITM

include Makefile.def

ABDIR   = srcSphereAB
EDDIR   = srcRudy
IEDIR   = srcIE
IODIR   = srcIO
MAINDIR = src

install: Makefile.def.orig MAKEFILE_DEF
	@make install_cont;

Makefile.def.orig:
	mv Makefile.def Makefile.def.orig
	cp Makefile.def.orig Makefile.def

MAKEFILE_DEF:
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		echo UADIR=`pwd`                        >  Makefile.def; \
		echo OS=`uname`                         >> Makefile.def; \
		echo COMPILER=`./get_info.pl -c`        >> Makefile.def; \
		echo MPIVERSION=`./get_info.pl -m`      >> Makefile.def; \
		echo STANDALONE=${STANDALONE}           >> Makefile.def; \
		cat srcMake/Makefile.def                >> Makefile.def; \
	fi);

PLANET=earth

install_cont: 
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cp -f share/build/Makefile.${OS}${COMPILER} Makefile.conf; \
		cd share; make COMPILER=${COMPILER} MPIVERSION=${MPIVERSION} install;\
	else \
		echo include $(DIR)/Makefile.conf > Makefile.conf; \
	fi);
	@(if [ -f src/Makefile.RULES.${OS}${COMPILER} ]; then                \
		cp -f src/Makefile.RULES.${OS}${COMPILER} src/Makefile.RULES;\
	else \
		rm -f src/Makefile.RULES; touch src/Makefile.RULES; \
	fi);
	touch src/Makefile.DEPEND srcInterface/Makefile.DEPEND
	cd src; make DYNAMIC
	./config.pl -${PLANET} -compiler=${COMPILER}

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

GITM2 = ${DIR}/UA/GITM2

LIB:
	cd $(ABDIR)     ; make                                         LIB
	cd $(EDDIR)     ; make LIBPREV=${GITM2}/${ABDIR}/libSphere.a   LIBADD
	cd $(IODIR)     ; make LIBPREV=${GITM2}/${EDDIR}/libABIEED.a   LIBADD
	cd $(MAINDIR)   ; make LIBPREV=${GITM2}/${IODIR}/libABIEEDIO.a LIB
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

distclean: clean
	@(if [ -d share ]; then cd share; make distclean; fi);
	@cd $(MAINDIR); make distclean
	@(if [ -e src/Makefile.orig ]; then mv src/Makefile.orig src/Makefile; fi);
	@cd srcInterface; make distclean
	rm -f Makefile.conf Makefile.def *~
	mv Makefile.def.orig Makefile.def

#
#       Create run directories
#
rundir:
	mkdir -p ${RUNDIR}/UA
	cd ${RUNDIR}/UA; \
		mkdir restartOUT data; \
		ln -s restartOUT restartIN; \
		ln -s ${UADIR}/srcData DataIn
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR} ; \
		ln -s ${BINDIR}/GITM.exe . ; \
		cp UA/DataIn/UAM.in . ; \
		touch core ; chmod 444 core ; \
		ln -s UA/* .; \
	fi);

test:
	make test_compile
	make test_rundir
	make test_run
	make test_check

test_compile:
	make clean
	./GridSize.pl -g=9,9,50,4
	make GITM IEDIR=`pwd`/srcIE

TESTDIR = run_test

MPIRUN = mpirun -np 2

test_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE=YES UADIR=`pwd`
	cd ${TESTDIR}; cp UA/DataIn/UAM.in.test UAM.in

test_run:
	cd ${TESTDIR}; ${MPIRUN} ./GITM.exe > runlog

test_check:
	-@(${SCRIPTDIR}/DiffNum.pl -b \
		${TESTDIR}/UA/data/log00000002.dat \
		srcData/log00000002.dat > test.diff)
	ls -l test.diff

dist:
	make distclean
	tar cvzf gitm_`date "+%y%m%d"`.tgz Makefile* config.pl get_info.pl \
	    share src srcData srcIDL srcIE srcIO srcInterface \
	    srcMake srcRudy srcSphereAB
	make install
