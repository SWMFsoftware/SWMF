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
	touch src/Makefile.DEPEND

PWOM:
	@cd ${SHAREDIR}; make LIB
	@cd src;         make PWOM

LIB:
	cd src; make LIB

test:	run PWOM
	cd run; PWOM.exe
	diff -b run/PW/restartOUT/restart_iline0001.dat \
		output/${PLANET}/restart_iline0001.dat
	diff -b run/PW/plots/plots_iline0001.dat output/${PLANET}/plots_iline0001.dat
	diff -b run/log.out output/${PLANET}/log.out

testrun: PWOM
	cd run; PWOM.exe

rundir: run

run:
	mkdir -p ${RUNDIR}/PW
	cd ${RUNDIR}/PW; \
		mkdir restartIN restartOUT plots; \
		cp ${PWDIR}/input/${PLANET}/restart_iline0001.dat restartIN/ ;\
		cp ${PWDIR}/input/${PLANET}/North.dat .
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR}; \
		ln -s ${BINDIR}/PWOM.exe .; \
		cp ${PWDIR}/input/${PLANET}/pw.input .; \
		touch core ; chmod 444 core ; \
		ln -s PW/* .; \
	fi)

clean:
	cd src; make clean

distclean: uninstall

uninstall:
	@(if [ -d share ]; then cd share; make distclean; fi);
	@cd src; make distclean
	rm -f Makefile.conf Makefile.def Makefile.planet *~
	mv Makefile.planet.orig Makefile.planet
	mv Makefile.def.orig Makefile.def
	@(if [ "$(STANDALONE)" != "NO" ]; then rm -rf bin; fi)
