#
# List the default target first for stand alone mode
#
DEFAULT_TARGET = BATSRUS
DEFAULT_EXE    = BATSRUS.exe

default : ${DEFAULT_TARGET}

include Makefile.def

#
# Menu of make options
#
help:
	@echo ' '
	@echo '  You can "make" the following:'
	@echo ' '
	@echo '    <default> BATSRUS in stand alone mode, help in SWMF'
	@echo ' '
	@echo '    install   (create MAKEFILES, src/ModSize.f90, src/mpif90.h)'	
	@echo '    MAKEFILE_DEF (create/correct Makefile.def)'
	@echo ' '
	@echo '    LIB     (Component library libGM for SWMF)'
	@echo '    BATSRUS (Block Adaptive Tree Solar-Wind Roe Upwind Scheme)'
	@echo '    NOMPI   (NOMPI library for compilation without MPI)'
	@echo '    PIDL    (PostIDL program creates 1 .out file from local .idl files)'
	@echo '    PSPH    (PostSPH program creates spherical tec file from sph*.tec files)'
	@echo ' '
	@echo '    help      (makefile option list)'
	@echo '    clean     (rm -f *~ *.o *.kmo *.mod *.T *.lst core)'
	@echo '    distclean (make clean; rm -f *exe Makefile Makefile.DEPEND)'
	@echo '    dist      (create source distribution tar file)'
	@#^CFG IF DOC BEGIN
	@#	^CFG IF NOT REMOVEDOCTEX BEGIN
	@echo ' '
	@echo '    PDF       (Make PDF version of the documentation)'
	@#		^CFG IF DOCHTML BEGIN
	@echo '    HTML      (Make HTML version of the documentation)'
	@#		^CFG END DOCHTML
	@#	^CFG END REMOVEDOCTEX
	@#^CFG END DOC
	@echo ' '
	@echo '    rundir      (create run directory for standalone or SWMF)'
	@echo ' '
	@echo "    mpirun      (make BATSRUS and mpirun BATSRUS.exe on 8 PEs)"
	@echo "    mpirun NP=7 (make BATSRUS and mpirun BATSRUS.exe on 7 PEs)"
	@echo "    mprun NP=5  (make BATSRUS and mprun  BATSRUS.exe on 5 PEs)"
	@echo ' '
	@#^CFG IF CONFIGURE BEGIN
	@#    ^CFG IF COVARIANT BEGIN
	@echo "    spherical_src  (Make SPHERICAL directory with BATSRUS on spherical grid)"
	@echo "    spherical_conf (Make SPHERICAL directory and link it to BATSRUS_conf)"
	@echo "    covariant_src  (Make COVARIANT directory with BATSRUS on covariant grid)"
	@echo "    corelax_src    (Make CORELAX directory for the covariant version "
	@echo "                of the magnetogram-driven solar wind)"
	@#    ^CFG END COVARIANT
	@#    ^CFG IF NOT COVARIANT BEGIN
	@echo "    cartesian_src  (removes source code for covariant grid)"
	@echo "    relax_src      (Make RELAX directory for the Cartesian "
	@echo "                    version of the magnetogram-driven solar wind)"
	@#    ^CFG END COVARIANT
	@echo ' '
	@#^CFG END CONFIGURE


install: Makefile.def.orig MAKEFILE_DEF
	@make install_cont;

Makefile.def.orig:
	mv Makefile.def Makefile.def.orig
	cp Makefile.def.orig Makefile.def

MAKEFILE_DEF:
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		echo DIR=`pwd`                          >  Makefile.def; \
		echo GMDIR=`pwd`                        >> Makefile.def; \
		echo OS=`uname`                         >> Makefile.def; \
		echo STANDALONE=${STANDALONE}           >> Makefile.def; \
		cat src/Makefile.def                    >> Makefile.def; \
	fi);

install_cont: src/ModSize.f90
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cp -f share/build/Makefile.${OS}${COMPILER} Makefile.conf; \
		cp -f src/stand_alone_${STANDALONE}.f90 src/stand_alone.f90;\
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
	./Options.pl -u=Default -e=Mhd
	cd src; make STATIC

src/ModSize.f90:
	cp -f src/ModSize_orig.f90 src/ModSize.f90

LIB:
	cd src; make LIB
	cd srcInterface; make LIB

BATSRUS:
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd src; make LIB
	cd src; make BATSRUS.exe
	@echo ' '
	@echo Program BATSRUS has been brought up to date.
	@echo ' '	

NOMPI:
	cd util/NOMPI/src; make LIB

PIDL:
	cd srcPostProc; make PIDL
	@echo ' '
	@echo Program PostIDL has been brought up to date.
	@echo ' '

PSPH:
	cd srcPostProc; make PSPH
	@echo ' '
	@echo Program PostSPH has been brought up to date.
	@echo ' '

# The MACHINE variable holds the machine name for which scripts should
# be copied to the run directory when it is created.  This is used mostly
# when several different machines have the same operating system,
# but they require different batch queue scripts.
# If MACHINE is empty or not defined, all scripts for the current OS will
# be copied.
#
# The default is the short name of the current machine
MACHINE = `hostname | sed -e 's/\..*//'`

rundir:
	mkdir -p ${RUNDIR}/GM
	cd ${RUNDIR}/GM; \
		mkdir restartIN restartOUT IO2; \
		ln -s ${BINDIR}/PostIDL.exe .; \
		ln -s ${BINDIR}/PostSPH.exe .; \
		cp    ${GMDIR}/Scripts/IDL/pIDL .; \
		cp    ${GMDIR}/Scripts/TEC/pTEC .; \
		ln -s ${GMDIR}/Param .
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cp -f Param/PARAM.DEFAULT run/PARAM.in; \
		touch run/core; chmod 444 run/core; \
		touch Scripts/Run/${OS}/TMP_${MACHINE}; \
		cp Scripts/Run/${OS}/*${MACHINE}* run/; \
		rm -f run/TMP_${MACHINE} Scripts/Run/${OS}/TMP_${MACHINE}; \
		cp share/Scripts/PostProc.pl run/; \
		cp share/Scripts/Restart.pl run/; \
		cd run; ln -s ${BINDIR}/BATSRUS.exe .; \
		ln -s GM/* .;                          \
	fi);

#
#       Run the default code on NP processors
#

NP=8

mpirun: ${DEFAULT_TARGET}
	cd run; mpirun -np ${NP} ./${DEFAULT_EXE}

mprun: ${DEFAULT_TARGET}
	cd run; mprun -np ${NP} ./${DEFAULT_EXE}

nompirun: ${DEFAULT_TARGET}
	cd run; ./${DEFAULT_EXE}

#					^CFG IF DOC BEGIN
#	Create the documentation files      ^CFG IF NOT REMOVEDOCTEX BEGIN
#	
PDF:
	@cd Doc/Tex; make cleanpdf; make PDF

CLEAN1 = cleanpdf #				^CFG IF NOT MAKEPDF

#	Create HTML documentation		^CFG IF DOCHTML BEGIN
HTML:
	@cd Doc/Tex; make cleanhtml; make HTML

CLEAN2 = cleanhtml #				    ^CFG IF NOT MAKEHTML
#						^CFG END DOCHTML
#					    ^CFG END REMOVEDOCTEX
#					^CFG END DOC

#
# Cleaning
#

clean:
	@touch src/Makefile.DEPEND src/Makefile.RULES
	cd src; make clean
	@touch srcInterface/Makefile.DEPEND
	cd srcInterface; make clean
	cd srcPostProc;  make clean
	@(if [ -d util  ]; then cd util;  make clean; fi);
	@(if [ -d share ]; then cd share; make clean; fi);

distclean:
	@touch src/Makefile.DEPEND src/Makefile.RULES
	cd src; make distclean
	@touch srcInterface/Makefile.DEPEND
	cd srcInterface; make distclean
	cd srcPostProc;  make distclean
	@				#^CFG IF DOC BEGIN
	@					#^CFG IF NOT REMOVEDOCTEX BEGIN
	cd Doc/Tex; make clean ${CLEAN1} ${CLEAN2}
	@					#^CFG END REMOVEDOCTEX
	@				#^CFG END DOC
	@(if [ -d util  ]; then cd util;  make distclean; fi);
	@(if [ -d share ]; then cd share; make distclean; fi);
	rm -f Makefile.conf Makefile.def *~
	mv Makefile.def.orig Makefile.def

dist:	distclean
	@echo ' '
	@echo ' NOTE: All "run" or other created directories not included!'
	@echo ' '
	tar -cf tmp.tar  Makefile Makefile_CONFIGURE Makefile.def
	tar -rf tmp.tar  Copyrights
	tar -rf tmp.tar  CVS* .cvsignore	#^CFG IF CONFIGURE
	tar -rf tmp.tar  Configure.options	#^CFG IF CONFIGURE
	tar -rf tmp.tar  Configure.pl		#^CFG IF CONFIGURE
	tar -rf tmp.tar  Test*.pl TestCovariant	#^CFG IF TESTING
	tar -rf tmp.tar  Doc			#^CFG IF DOC
	tar -rf tmp.tar  PARAM.XML PARAM.pl
	tar -rf tmp.tar  Options.pl GridSize.pl
	tar -rf tmp.tar  Idl
	tar -rf tmp.tar  Param
	tar -rf tmp.tar  Scripts
	tar -rf tmp.tar  src srcInterface srcPostProc srcUser
	@(if [ -d util  ]; then tar -rf tmp.tar util; fi);
	@(if [ -d share ]; then tar -rf tmp.tar share; fi);
	@echo ' '
	gzip tmp.tar
	mv tmp.tar.gz BATSRUS_v${VERSION}_`date +%Y%b%d_%H%M.tgz`
	@echo ' '
	@ls -l BATSRUS_v*.tgz

include Makefile_CONFIGURE #^CFG IF CONFIGURE


test_titan:
	make test_titan_compile
	make test_titan_rundir
	make test_titan_run
	make test_titan_check
	make test_titan_restart

test_titan_compile:
	make clean
	Options.pl -u=Titan -e=MhdTitan
	GridSize.pl -g=4,4,4,2000,200
	make
	make PIDL

test_titan_rundir:
	make rundir
	cd run; \
		cp GM/Param/TITAN/PARAM.in .; \
		tar xzf GM/Param/TITAN/TitanInput.tgz

test_titan_run:
	cd run; ./BATSRUS.exe > runlog; PostProc.pl -M -o TitanTest

test_titan_check:
	-(share/Scripts/DiffNum.pl \
		Param/TITAN/TestOutput/log_n000001.log \
		run/TitanTest/GM/log_n000001.log > titan_diff.log)
	ls -l titan_diff.log titan_diff.runlog

test_titan_restart:
	make test_titan_restart_save
	make test_titan_restart_read
	make test_titan_restart_check

test_titan_restart_save:
	cd run; \
		cp GM/Param/TITAN/PARAM.in.restartsave PARAM.in; \
		./BATSRUS.exe > runlog; \
		PostProc.pl -M -o TitanTest/RestartSave; \
		Restart.pl -o RESTART_titan; 

test_titan_restart_read:
	cd run; \
		cp GM/Param/TITAN/PARAM.in.restartread PARAM.in; \
		Restart.pl -i RESTART_titan; 
		./BATSRUS.exe > runlog; \
		PostProc.pl -M -o TitanTest/RestartRead; \

test_titan_restart_check:
	cd run/TitanTest; \
		cp RestartSave/GM/log_n000001.log log_all.log; \
		tail -25 RestartRead/GM/log_n000001.log >> log_all.log
	-(share/Scripts/DiffNum.pl \
		Param/TITAN/TestOutput/log_n000001.log \
		run/TitanTest/log_all.log > titan_diff.log)
