default: PWOM

include Makefile.def
include Makefile.planet

help:
	@echo Makefile targets:
	@echo
	@echo 'make                        - compile PWOM.exe'
	@echo 'make LIB                    - compile libPW.a for SWMF'
	@echo 'make run                    - create run directory'
	@echo 'make test                   - test PWOM in stand alone mode'
	@echo 'make clean                  - remove object files'
	@echo 'make distclean              - remove all files not part of CVS'
	@echo

bin:
	mkdir bin

INSTALLFILES =  src/Makefile.DEPEND \
		src/Makefile.RULES \
		srcInterface/Makefile.DEPEND

install: bin
	touch ${INSTALLFILES}

PWOM:
	cd ${SHAREDIR};           make LIB
	cd ${TIMINGDIR};          make LIB
	cd ${EMPIRICALIEDIR};     make LIB
	cd ${EMPIRICALUADIR};    make  LIB
	cd ${DATAREADINDICESDIR}; make LIB
	cd src;                  make PWOM

nompirun: PWOM
	cd ${RUNDIR}; ./PWOM.exe

LIB:
	cd src; make LIB
	cd srcInterface; make LIB

TESTDIR = run_test

MPIRUN = mpirun -np 2

test:
	-@(make test_saturn)
	-@(make test_restart)
	-@(make test_earth)
	-@(make test_restart)

test_earth:	
	@echo "starting..." > test_Earth_plots.diff
	@echo "test_compile..." > test_Earth.diff
	make test_compile PLANET=Earth
	@echo "test_rundir..." >> test_Earth.diff
	make test_rundir
	@echo "test_run..." >> test_Earth.diff
	make test_run
	@echo "test_check..." >> test_Earth.diff
	make test_check

test_saturn:
	@echo "starting..." > test_Saturn_plots.diff
	@echo "test_compile..." > test_Saturn.diff
	make test_compile PLANET=Saturn
	@echo "test_rundir..." >> test_Saturn.diff
	make test_rundir
	@echo "test_run..." >> test_Saturn.diff
	make test_run
	@echo "test_check..." >> test_Saturn.diff
	make test_check

test_compile:
	./Config.pl -${PLANET}
	make PWOM

test_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE="YES" PWDIR=`pwd`

test_run: PWOM
	cd ${TESTDIR}; ${MPIRUN} ./PWOM.exe

test_check:
	-@(${SCRIPTDIR}/DiffNum.pl -b -r=1e-9 \
		${TESTDIR}/PW/restartOUT/restart_iline0001.dat \
		data/output/${PLANET}/restart_iline0001.dat \
		> test_${PLANET}${MYTEST}.diff)
	-@(${SCRIPTDIR}/DiffNum.pl -b -r=1e-9 \
		${TESTDIR}/PW/restartOUT/restart_iline0002.dat \
		data/output/${PLANET}/restart_iline0002.dat \
		>> test_${PLANET}${MYTEST}.diff)
	-@(${SCRIPTDIR}/DiffNum.pl -b -r=1e-9 \
		${TESTDIR}/PW/restartOUT/restart_iline0003.dat \
		data/output/${PLANET}/restart_iline0003.dat \
		>> test_${PLANET}${MYTEST}.diff)
	-@(${SCRIPTDIR}/DiffNum.pl -b -r=1e-9 \
		${TESTDIR}/PW/restartOUT/restart_iline0004.dat \
		data/output/${PLANET}/restart_iline0004.dat \
		>> test_${PLANET}${MYTEST}.diff)
	-@(${SCRIPTDIR}/DiffNum.pl -b -r=1e-8 \
		${TESTDIR}/PW/plots/north_plots_iline0001.out \
		data/output/${PLANET}/north_plots_iline0001.out \
		> test_${PLANET}${MYTEST}_plots.diff)
	ls -l *.diff

test_restart:
	@echo "test_restart_save..." > test_${PLANET}_restart.diff
	make   test_restart_save
	@echo "test_restart_read..." >> test_${PLANET}_restart.diff
	make   test_restart_read
	@echo "test_restart_check..." >> test_${PLANET}_restart.diff
	make   test_check MYTEST=_restart

test_restart_save:
	cp data/input/${PLANET}/restart_iline* ${TESTDIR}/PW/restartIN/
	cp input/${PLANET}/PARAM.in.restartsave ${TESTDIR}/PARAM.in
	cd ${TESTDIR}; ${MPIRUN} ./PWOM.exe
	cd ${TESTDIR}; mv PW/restartOUT/* PW/restartIN/
	cd ${TESTDIR}/PW; rm -rf plot_save; mv plots plots_save; mkdir plots

test_restart_read:
	cp input/${PLANET}/PARAM.in.restartread ${TESTDIR}/PARAM.in
	cd ${TESTDIR}; ${MPIRUN} ./PWOM.exe
	cd ${TESTDIR}/PW; \
	rm -rf plots_read; mv plots plots_read; mv plots_save plots
	cd ${TESTDIR}/PW; \
	cat plots_read/north_plots_iline0001.out >> plots/north_plots_iline0001.out;\
	cat plots_read/north_plots_iline0002.out >> plots/north_plots_iline0002.out;\
	cat plots_read/north_plots_iline0003.out >> plots/north_plots_iline0003.out;\
	cat plots_read/north_plots_iline0004.out >> plots/north_plots_iline0004.out

rundir:
	mkdir -p ${RUNDIR}/PW
	@(cd ${RUNDIR}; \
		if [ ! -e "EIE/README" ]; then \
			ln -s ${EMPIRICALIEDIR}/data EIE;\
		fi;)
	cd ${RUNDIR}/PW; \
		mkdir restartIN restartOUT plots; \
		cp ${PWDIR}/data/input/${PLANET}/restart_iline* restartIN/ ;\
		cp ${PWDIR}/data/input/${PLANET}/North.dat .
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR}; \
			ln -s ${BINDIR}/PWOM.exe .; \
			cp ${PWDIR}/input/${PLANET}/PARAM.in .; \
			touch core ; chmod 444 core ; \
			ln -s PW/* .; \
	fi)

PDF:
	@cd doc/Tex; make PDF

clean:
	@touch ${INSTALLFILES}
	cd src; make clean
	cd srcInterface; make clean
	cd doc/Tex; make clean
	@(if [ -d share ]; then cd share; make clean; fi);
	@(if [ -d util ];  then cd util;  make clean; fi);

distclean: 
	./Config.pl -uninstall

allclean:
	@touch ${INSTALLFILES}
	cd src; make distclean
	cd srcInterface; make distclean
	cd doc/Tex; make distclean
	rm -rf Makefile.planet *~ bin
