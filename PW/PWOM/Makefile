default: PWOM

include Makefile.def
include Makefile.planet

help:
	@echo Makefile targets:
	@echo
	@echo 'make install PLANET=Saturn COMPILER=g95 MPIVERSION=Altix'
	@echo
	@echo '      installs PWOM for Saturn with g95 compiler and MPI for Altix'
	@echo '      defaults are Earth, the default compiler for the OS and MPICH'
	@echo
	@echo 'make                        - compile PWOM.exe'
	@echo 'make LIB                    - compile libPW.a for SWMF'
	@echo 'make run                    - create run directory'
	@echo 'make test                   - test PWOM in stand alone mode'
	@echo 'make clean                  - remove object files'
	@echo 'make distclean              - remove all files not part of CVS'
	@echo

install: Makefile.def.orig MAKEFILE_DEF bin src/PLANET
	@make install_cont;

Makefile.def.orig:
	mv Makefile.def Makefile.def.orig
	cp Makefile.def.orig Makefile.def

MAKEFILE_DEF:
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		echo PWDIR=`pwd`                        >  Makefile.def; \
		echo OS=`uname`                         >> Makefile.def; \
		echo STANDALONE=${STANDALONE}           >> Makefile.def; \
		cat src/Makefile.def                    >> Makefile.def; \
	fi);

bin:
	mkdir bin

src/PLANET:
	mv Makefile.planet Makefile.planet.orig
	./config.pl -p=${PLANET}

install_cont:
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cp -f share/build/Makefile.${OS}${COMPILER} Makefile.conf; \
		cd share; make install;\
	else \
		echo include $(DIR)/Makefile.conf > Makefile.conf; \
	fi);
	@(if [ -f src/Makefile.RULES.${OS}${COMPILER} ]; then                \
		cp -f src/Makefile.RULES.${OS}${COMPILER} src/Makefile.RULES;\
	else \
		rm -f src/Makefile.RULES; touch src/Makefile.RULES; \
	fi);
	touch src/Makefile.DEPEND srcInterface/Makefile.DEPEND

PWOM:
	@cd ${SHAREDIR}; make LIB
	@cd src;         make PWOM

LIB:
	cd src; make LIB
	cd srcInterface; make LIB

TESTDIR = run_test

MPIRUN = mpirun -np 2

test:
	-@(make test_saturn)
	-@(make test_earth)

test_earth:	
	@echo "starting..." > test_Earth_plots.diff
	@echo "starting..." > test_Earth_restart.diff
	@echo "test_compile..." > test_Earth.diff
	make test_compile PLANET=Earth
	@echo "test_rundir..." >> test_Earth.diff
	make test_rundir
	@echo "test_run..." >> test_Earth.diff
	make test_run
	@echo "test_check..." >> test_Earth.diff
	make test_check

test_earth_compile:
	./config.pl -p=Earth
	make PWOM

test_saturn:
	@echo "starting..." > test_Saturn_plots.diff
	@echo "starting..." > test_Saturn_restart.diff
	@echo "test_compile..." > test_Saturn.diff
	make test_compile PLANET=Saturn
	@echo "test_rundir..." >> test_Saturn.diff
	make test_rundir
	@echo "test_run..." >> test_Saturn.diff
	make test_run
	@echo "test_check..." >> test_Saturn.diff
	make test_check

test_compile:
	./config.pl -p=${PLANET}
	make PWOM

test_rundir:
	rm -rf ${TESTDIR}
	make rundir RUNDIR=${TESTDIR} STANDALONE="YES" PWDIR=`pwd`

test_run: PWOM
	cd ${TESTDIR}; ${MPIRUN} ./PWOM.exe

test_check:
	-@(${SCRIPTDIR}/DiffNum.pl -b \
		${TESTDIR}/PW/log.out \
		output/${PLANET}/log.out \
		> test_${PLANET}.diff)
	-@(${SCRIPTDIR}/DiffNum.pl -b -r=1e-9 \
		${TESTDIR}/PW/restartOUT/restart_iline0001.dat \
		output/${PLANET}/restart_iline0001.dat \
		> test_${PLANET}_restart.diff)
	-@(${SCRIPTDIR}/DiffNum.pl -b -r=1e-9 \
		${TESTDIR}/PW/restartOUT/restart_iline0002.dat \
		output/${PLANET}/restart_iline0002.dat \
		>> test_${PLANET}_restart.diff)
	-@(${SCRIPTDIR}/DiffNum.pl -b -r=1e-9 \
		${TESTDIR}/PW/restartOUT/restart_iline0003.dat \
		output/${PLANET}/restart_iline0003.dat \
		>> test_${PLANET}_restart.diff)
	-@(${SCRIPTDIR}/DiffNum.pl -b -r=1e-9 \
		${TESTDIR}/PW/restartOUT/restart_iline0004.dat \
		output/${PLANET}/restart_iline0004.dat \
		>> test_${PLANET}_restart.diff)
	-@(${SCRIPTDIR}/DiffNum.pl -b -r=1e-8 \
		${TESTDIR}/PW/plots/plots_iline0001.out \
		output/${PLANET}/plots_iline0001.out \
		> test_${PLANET}_plots.diff)
	ls -l *.diff

rundir: run

run:
	mkdir -p ${RUNDIR}/PW
	cd ${RUNDIR}/PW; \
		mkdir restartIN restartOUT plots; \
		cp ${PWDIR}/input/${PLANET}/restart_iline* restartIN/ ;\
		cp ${PWDIR}/input/${PLANET}/North.dat .
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR}; \
		ln -s ${BINDIR}/PWOM.exe .; \
		cp ${PWDIR}/input/${PLANET}/pw.input .; \
		touch core ; chmod 444 core ; \
		ln -s PW/* .; \
	fi)

clean:
	@touch src/Makefile.DEPEND src/Makefile.RULES
	cd src; make clean
	@touch srcInterface/Makefile.DEPEND
	cd srcInterface; make clean
	@(if [ -d share ]; then cd share; make clean; fi);

distclean: uninstall

uninstall:
	@touch src/Makefile.DEPEND src/Makefile.RULES
	cd src; make distclean
	@touch srcInterface/Makefile.DEPEND
	cd srcInterface; make distclean
	@(if [ -d share ]; then cd share; make distclean; fi);
	rm -rf Makefile.conf Makefile.def Makefile.planet *~ bin
	mv Makefile.planet.orig Makefile.planet
	mv Makefile.def.orig Makefile.def
