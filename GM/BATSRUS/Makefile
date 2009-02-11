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
	@echo '    <default> ${DEFAULT_TARGET} in stand alone mode, help in SWMF'
	@echo ' '
	@echo '    help         (makefile option list)'
	@echo '    install      (install BATSRUS)'
	@#^CFG IF DOC BEGIN
	@#	^CFG IF NOT REMOVEDOCTEX BEGIN
	@echo ' '
	@echo '    PDF          (Make PDF version of the documentation)'
	@#		^CFG IF DOCHTML BEGIN
	@echo '    HTML         (Make HTML version of the documentation)'
	@#		^CFG END DOCHTML
	@#	^CFG END REMOVEDOCTEX
	@#^CFG END DOC
	@#^CFG IF TESTING BEGIN
	@echo '    test         (run all tests for BATSRUS)'
	@echo '    test_help    (show all options for running the tests)'
	@#^CFG END TESTING
	@#^CFG IF CONFIGURE BEGIN
	@echo '    config_help  (show all targets in Makefile_CONFIGURE'
	@#^CFG END CONFIGURE
	@echo ' '
	@echo '    LIB     (Component library libGM for SWMF)'
	@echo '    BATSRUS (Block Adaptive Tree Solar-Wind Roe Upwind Scheme)'
	@echo '    CRASH   (Code for Radiative Shock Hydrodynamics)'
	@echo '    NOMPI   (NOMPI library for compilation without MPI)'
	@echo '    PIDL    (PostIDL.exe creates 1 .out file from local .idl files)'
	@echo '    PSPH    (PostSPH.exe creates spherical tec file from sph*.tec files)'
	@echo '    EARTH_TRAJ (EARTH_TRAJ.exe creates Earth trajectory file for heliosphere)'
	@echo ' '
	@echo '    rundir      (create run directory for standalone or SWMF)'
	@echo '    rundir RUNDIR=run_test (create run directory run_test)'
	@echo ' '
	@echo "    nompirun    (make and run ${DEFAULT_EXE} on 1 PE)"
	@echo "    mpirun      (make and mpirun ${DEFAULT_EXE} on 8 PEs)"
	@echo "    mpirun NP=7 RUNDIR=run_test (run on 7 PEs in run_test)"
	@echo "    mprun NP=5  (make and mprun ${DEFAULT_EXE} on 5 PEs)"
	@echo ' '	
	@echo '    clean     (remove temp files like: *~ *.o *.kmo *.mod *.T *.lst core)'
	@echo '    distclean (equivalent to ./Config.pl -uninstall)'
	@echo '    dist      (create source distribution tar file)'

INSTALLFILES =	src/Makefile.DEPEND \
		src/Makefile.RULES \
		srcInterface/Makefile.DEPEND \
		srcPostProc/Makefile.RULES

install: src/ModSize.f90
	touch ${INSTALLFILES}
	./Config.pl -u=Default -e=Mhd
	cd src; make STATIC

src/ModSize.f90:
	cp -f src/ModSize_orig.f90 src/ModSize.f90

LIB:
	cd src; make LIB
	cd srcInterface; make LIB

BATSRUS:
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd ${MAGNETOGRAMDIR}; make LIB
	cd ${EMPIRICALEEDIR}; make LIB
	cd src; make LIB
	cd src; make BATSRUS

CRASH:
	cd ${SHAREDIR}; make LIB
	cd ${TIMINGDIR}; make LIB
	cd ${MAGNETOGRAMDIR}; make LIB
	cd ${EMPIRICALEEDIR}; make LIB
	cd ${CRASHDIR}; make LIB
	cd src; make LIB
	cd src; make CRASH

NOMPI:
	cd util/NOMPI/src; make LIB

PIDL:
	cd ${SHAREDIR}; make LIB
	cd srcPostProc; make PIDL
	@echo ' '
	@echo Program PostIDL has been brought up to date.
	@echo ' '

PSPH:
	cd srcPostProc; make PSPH
	@echo ' '
	@echo Program PostSPH has been brought up to date.
	@echo ' '

EARTH_TRAJ:
	cd srcPostProc; make EARTH_TRAJ
	@echo ' '
	@echo Program EARTH_TRAJ has been brought up to date.
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

COMPONENT = GM

rundir:
	mkdir -p ${RUNDIR}/${COMPONENT}
	cd ${RUNDIR}/${COMPONENT}; \
		mkdir restartIN restartOUT IO2; \
		ln -s ${BINDIR}/PostIDL.exe .; \
		ln -s ${BINDIR}/PostSPH.exe .; \
		cp    ${GMDIR}/Scripts/IDL/pIDL .; \
		cp    ${GMDIR}/Scripts/TEC/pTEC .; \
		ln -s ${GMDIR}/Param .
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cp -f Param/PARAM.DEFAULT ${RUNDIR}/PARAM.in; \
		touch ${RUNDIR}/core; chmod 444 ${RUNDIR}/core; \
		touch Scripts/Run/${OS}/TMP_${MACHINE}; \
		cp Scripts/Run/${OS}/*${MACHINE}* ${RUNDIR}/; \
		rm -f ${RUNDIR}/TMP_${MACHINE}; \
		rm -f Scripts/Run/${OS}/TMP_${MACHINE}; \
		cp ${SCRIPTDIR}/PostProc.pl ${RUNDIR}/; \
		cp ${SCRIPTDIR}/Restart.pl ${RUNDIR}/; \
		cd ${RUNDIR}; ln -s ${BINDIR}/${DEFAULT_EXE} .; \
		ln -s ${COMPONENT}/* .;                          \
	fi);

rundir_rh:
	make rundir RUNDIR=run STANDALONE=YES GMDIR=`pwd`
	cd run; ln -s ${BINDIR}/CRASH.exe .; \
	cp -f Param/CRASH/PARAM.in .                          
#
#       Run the default code on NP processors
#

NP=8

mpirun: ${DEFAULT_TARGET}
	cd ${RUNDIR}; mpirun -np ${NP} ./${DEFAULT_EXE}

mprun: ${DEFAULT_TARGET}
	cd ${RUNDIR}; mprun -np ${NP} ./${DEFAULT_EXE}

nompirun: ${DEFAULT_TARGET}
	cd ${RUNDIR}; ./${DEFAULT_EXE}

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
	@touch ${INSTALLFILES}
	cd src; make clean
	cd srcInterface; make clean
	cd srcPostProc;  make clean
	@(if [ -d util  ]; then cd util;  make clean; fi);
	@(if [ -d share ]; then cd share; make clean; fi);

distclean: 
	./Config.pl -uninstall

allclean:
	@touch ${INSTALLFILES}
	cd src; make distclean
	cd srcInterface; make distclean
	cd srcPostProc;  make distclean
	@				#^CFG IF DOC BEGIN
	@					#^CFG IF NOT REMOVEDOCTEX BEGIN
	cd Doc/Tex; make clean ${CLEAN1} ${CLEAN2}
	@					#^CFG END REMOVEDOCTEX
	@				#^CFG END DOC

dist:
	./Config.pl -uninstall
	@echo ' '
	@echo ' NOTE: All "run" or other created directories not included!'
	@echo ' '
	tar -cf tmp.tar  Makefile Makefile_CONFIGURE Makefile.test
	tar -rf tmp.tar  README Copyrights
	tar -rf tmp.tar  CVS* .cvsignore	#^CFG IF CONFIGURE
	tar -rf tmp.tar  Configure.options	#^CFG IF CONFIGURE
	tar -rf tmp.tar  Configure.pl		#^CFG IF CONFIGURE
	tar -rf tmp.tar  Test*.pl TestCovariant	#^CFG IF TESTING
	tar -rf tmp.tar  Doc			#^CFG IF DOC
	tar -rf tmp.tar  PARAM.XML
	tar -rf tmp.tar  Config.pl
	tar -rf tmp.tar  Idl
	tar -rf tmp.tar  Param
	tar -rf tmp.tar  Scripts
	tar -rf tmp.tar  src srcInterface srcPostProc srcTest srcUser
	@(if [ -d util  ]; then tar -rf tmp.tar util; fi);
	@(if [ -d share ]; then tar -rf tmp.tar share; fi);
	@echo ' '
	gzip tmp.tar
	mv tmp.tar.gz BATSRUS_v${VERSION}_`date +%Y%b%d_%H%M.tgz`
	@echo ' '
	@ls -l BATSRUS_v*.tgz

include Makefile_CONFIGURE #^CFG IF CONFIGURE

include Makefile.test #^CFG IF TESTING
