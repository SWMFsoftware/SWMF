#^CMP COPYRIGHT UM

#!QUOTE: \clearpage
#BOP
#!QUOTE: \section{Makefiles}
#!ROUTINE: Makefile - Space Weather Modeling Framework (SWMF) main Makefile
#!DESCRIPTION:
# This is the main Makefile to (re)compile executables, run the SWMF, 
# produce documentation, make distribution clean distribution, etc.
#EOP

SHELL=/bin/sh

#BOC
VERSION = 2.3

#
# The default target is SWMF so it is listed first
#
default : SWMF

#
# Make all common used executables.
#
ALL : SWMF PIDL PSPH EARTH_TRAJ

#
#       Definition of OS, component versions and directory structure
#
include Makefile.def

help:
	@echo ' '
	@echo '  You can "make" the following:'
	@echo ' '
	@echo '    [default]   SWMF'
	@echo ' '
	@echo '    help        (makefile option list)'
	@echo '    install     (to be used via Config.pl only)'
	@#^CMP IF NOT REMOVEDOCTEX BEGIN
	@echo ' '
	@echo '    PDF         (make PDF  version of the documentation)' #^CMP IF DOC
	@echo '    HTML        (make HTML version of the documentation)' #^CMP IF DOCHTML
	@#^CMP END  REMOVEDOCTEX
	@#^CFG IF TESTING BEGIN
	@echo '    test        (run all tests for the SWMF)'
	@echo '    test_help   (show all options for running the tests)'
	@echo ' '
	@#^CFG END TESTING
	@#^CMP IF IH BEGIN
	@echo '    IHBATSRUS   (copy and rename GM/BATSRUS source into IH/BATSRUS)'
	@#^CMP END IH
	@#^CMP IF SC BEGIN
	@echo '    SCBATSRUS   (configure and rename GM/BATSRUS source into SC/BATSRUS)'
	@#^CMP END SC
	@echo ' '
	@echo '    LIB         (lib/libSWMF.a the SWMF library)'
	@echo '    NOMPI       (lib/libNOMPI.a for single node execution with no MPI)'
	@echo ' '
	@echo '    SWMF        (bin/SWMF.exe the main executable for SWMF)'
	@echo '    ESMF_SWMF   (bin/ESMF_SWMF.exe the ESMF compatible executable)'
	@#^CMP IF IE BEGIN
	@echo '    PIONO       (bin/PostIONO.exe creates ionosphere tec file from idl files)'
	@#^CMP END IE
	@#^CMP IF GM BEGIN
	@echo '    PIDL        (bin/PostIDL.exe post-processes *.idl files)'
	@echo '    PSPH        (bin/PostSPH.exe post-processes sph*.tec files)'
	@echo '    EARTH_TRAJ  (bin/EARTH_TRAJ.exe generates Earth trajectory)'
	@echo '    ALL         (make SWMF PIDL PSPH EARTH_TRAJ)'
	@#^CMP END GM
	@echo ' '
	@echo ' '
	@echo '    rundir      (create run directory "run/")'
	@echo '    rundir RUNDIR=run_test  (create run directory "run_test/")'
	@echo '    rundir MACHINE=columbia (select job scripts for columbia)'
	@echo '    rundir_code (saves code archive in run directory)'
	@echo ' '
	@echo '    mpirun      (make SWMF and mpirun SWMF.exe on 2 PEs)'
	@echo '    mpirun NP=7 (make SWMF and mpirun SWMF.exe on 7 PEs)'
	@echo '    mprun  NP=5 (make SWMF and mprun  SWMF.exe on 5 PEs)'
	@echo '    nompirun    (make SWMF and run it without MPI)'
	@echo ' '
	@echo '    clean       (remove files: *~ *.o *.kmo *.mod *.T *.lst core)'
	@echo '    distlean    (remove all files which are not in distribution)'
	@echo '    dist        (create source distribution tar file)'
	@echo ' '
	@echo '    tags        (create etags for emacs for easy look up in source code)'
	@echo ' '
#EOC
#
# Check the variables DIR=`pwd` and OS = `uname`
#
ENV_CHECK:
	@(if([ "${DIR}" != `/bin/pwd` ]) || ([ "${OS}" != `uname` ]); \
	then \
	  echo; \
	  echo "DIR='${DIR}' and/or OS='${OS}' are incorrect!";\
	  echo "Correcting Makefile.def";\
	  perl -pi share/Scripts/FixMakefileDef.pl Makefile.def; \
	  echo "Type 'Config.pl -s' and retry the previous make command!"; \
	  exit 1; \
	fi);

#
# JOB SUBMISSION SCRIPTS TO COPY TO RUN
#
# The MACHINE variable holds the machine name for which scripts should
# be copied to the run directory when it is created.  This is used mostly
# when several different machines have the same operating system,
# but they require different batch queue scripts. 
# If MACHINE is empty or not defined, all scripts for the current OS will
# be copied.
#
# The default is the short name of the current machine
MACHINE = `hostname | sed -e 's/\..*//'`

#
# install for the 1st time
#

install: ENV_CHECK mkdir
	@echo VERSION ${VERSION}
	cd CON;			make install
	cd ESMF/ESMF_SWMF;	make install

	if([ -d "SC/BATSRUS" ]); \
		then cp GM/BATSRUS/Config.pl SC/BATSRUS; \
		     perl -i -pe 's/GM/SC/' SC/BATSRUS/Config.pl; \
	fi
	if([ -d "IH/BATSRUS" ]); \
		then cp GM/BATSRUS/Config.pl IH/BATSRUS; \
		     perl -i -pe 's/GM/IH/' IH/BATSRUS/Config.pl; \
	fi
	if([ -d "IH/BATSRUS_share" ]); \
		then cp GM/BATSRUS/Config.pl IH/BATSRUS_share; \
		     perl -i -pe "s/GM/IH/; s#'src#'../../GM/BATSRUS/src#" \
			IH/BATSRUS_share/Config.pl; \
	fi
	for i in `ls -d [A-Z][A-Z]/*/ | grep -v /CVS/ | grep -v /Empty/`; \
		do (echo Installing $$i; cd $$i; ./Config.pl -install=c); done
	@echo
	@echo Installation succeeded
	@echo

mkdir: ./lib ./bin

./lib:
	-@mkdir lib

./bin:
	-@mkdir bin

rmdir:
	-@rm -rf lib
	-@rm -rf bin

#
#	SWMF
#
SWMF:	ENV_CHECK
	@cd ${CONTROLDIR}; make SWMFEXE  
	@echo ' '

#
#       SWMF library (for external code)
#
LIB:	ENV_CHECK
	@cd ${CONTROLDIR}; make SWMFLIB
	@echo ' '

#
#       ESMF driver for the SWMF and a simple ESMF component
#
ESMF_SWMF: LIB
	@cd ESMF/ESMF_SWMF/src; make EXE

# NOMPI library for execution without MPI

NOMPI: ENV_CHECK
	cd ${NOMPIDIR}; make LIB

#^CMP IF GM BEGIN
#	Post processing codes for BATSRUS plot files
#
PSPH:	ENV_CHECK
	cd GM/BATSRUS; make PSPH
	@echo ' '

PIDL:	ENV_CHECK
	cd GM/BATSRUS; make PIDL
	@echo ' '

#
#	Code for generating Earth trajectory for IH/BATSRUS
#
EARTH_TRAJ: ENV_CHECK
	cd GM/BATSRUS; make EARTH_TRAJ
	@echo ' '
#^CMP END GM

#^CMP IF IE BEGIN
#	Post processing code for IE/Ridley_serial plot files
#
PIONO:	ENV_CHECK
	cd IE/Ridley_serial; make PIONO
	@echo ' '

#^CMP END IE
#
#	Graphical User Interface (GUI)
#
GUI:	ENV_CHECK
	cd gui; make install
	@echo ' '

GUI_doc: 
	cd gui; make doc
	@echo ' '

GUI_start: 
	gui/start.sh
	@echo 'GUI started.'

GUI_stop:
	gui/stop.sh
	@echo 'GUI stopped.'

#					^CMP IF DOC BEGIN
#	Create the documentation files      ^CMP IF NOT REMOVEDOCTEX BEGIN
#	
PDF:	ENV_CHECK
	@cd doc/Tex; make cleanpdf; make PDF

CLEAN1 = cleanpdf #				^CMP IF NOT MAKEPDF

#	Create HTML documentation		^CMP IF DOCHTML BEGIN
HTML:	ENV_CHECK
	@cd doc/Tex; make cleanhtml; make HTML

CLEAN2 = cleanhtml #				    ^CMP IF NOT MAKEHTML
#						^CMP END DOCHTML
#					    ^CMP END REMOVEDOCTEX
#					^CMP END DOC

#			
#	General Housekeeping
#
clean: ENV_CHECK
	@echo
	rm -rf *~ doc/*~ Param/*~ TAGS
	for i in `ls -d [A-Z][A-Z]/*/ | grep -v /CVS/`; \
		do (cd $$i; make clean); done
	-(cd ESMF/ESMF_SWMF;	make clean)
	cd CON;			make clean
	cd share;		make clean
	cd util;		make clean
	@#^CMP IF DOC BEGIN
	@#^CMP IF NOT REMOVEDOCTEX BEGIN
	cd doc/Tex;             make clean
	@#^CMP END REMOVEDOCTEX
	@#^CMP END DOC
	cd bin; rm -f *.exe 
	cd lib; rm -f *.a *.so
	@echo
	@echo Clean succeeded
	@echo

distclean: ENV_CHECK rmdir
	@echo
	rm -rf *~ doc/*~ Param/*~ TAGS
	for i in `ls -d [A-Z][A-Z]/Empty/`; \
		do (cd $$i; make distclean); done
	for i in `ls -d [A-Z][A-Z]/*/ | grep -v /CVS/ | grep -v Empty`; \
		do (cd $$i; ./Config.pl -uninstall); done
	-(cd ESMF/ESMF_SWMF;	make distclean)
	cd CON;			make distclean
	@#^CMP IF DOC BEGIN
	@#^CMP IF NOT REMOVEDOCTEX BEGIN
	cd doc/Tex;             make clean ${CLEAN1} ${CLEAN2}
	@#^CMP END REMOVEDOCTEX
	@#^CMP END DOC
	@echo
	@echo Distclean succeeded
	@echo

dist:
	@echo ' '
	@echo ' NOTE: All "run" or other created directories not included!'
	@echo ' '
	./Config.pl -uninstall
	-rm -rf */*/run_test
	tar -cf tmp.tar  README
	tar -rf tmp.tar  Makefile 
	tar -rf tmp.tar  Makefile.test          #^CMP FILE TESTING
	tar -rf tmp.tar  output                 #^CMP FILE TESTING
	tar -rf tmp.tar  Copyrights
	tar -rf tmp.tar  CVS*
	tar -rf tmp.tar  .cvsignore
	tar -rf tmp.tar  Config.pl
	tar -rf tmp.tar  doc			#^CMP IF DOC
	tar -rf tmp.tar  Param
	tar -rf tmp.tar  Scripts
	tar -rf tmp.tar  share
	tar -rf tmp.tar  util
	tar -rf tmp.tar  gui
	tar -rf tmp.tar  CON
	tar -rf tmp.tar  ESMF
	for i in `ls -d [A-Z][A-Z]`; \
		do (tar -rf tmp.tar $$i); done
	@echo ' '
	gzip tmp.tar
	mv tmp.tar.gz SWMF_v${VERSION}_`date +%Y%b%d_%H%M.tgz`
	@echo ' '
	@ls -l SWMF_v*.tgz

#			
#	Create run directories
#
${RUNDIR}:
	make rundir

rundir: ENV_CHECK
	mkdir -p ${RUNDIR}/STDOUT
	cp Param/LAYOUT.DEFAULT       ${RUNDIR}/LAYOUT.in
	cp Param/PARAM.DEFAULT        ${RUNDIR}/PARAM.in
	cp share/Scripts/Restart.pl   ${RUNDIR}/Restart.pl
	cp share/Scripts/PostProc.pl  ${RUNDIR}/PostProc.pl
	touch ${RUNDIR}/core
	chmod 444 ${RUNDIR}/core
	cd ${RUNDIR}; ln -s  ${DIR}/bin/SWMF.exe . ; ln -s  ${DIR}/Param .
	cd ${GMDIR}; make rundir                 #^CMP IF GM
	cd ${IEDIR}; make rundir                 #^CMP IF IE
	cd ${IHDIR}; make rundir                 #^CMP IF IH
	cd ${IMDIR}; make rundir                 #^CMP IF IM
	cd ${LADIR}; make rundir                 #^CMP IF LA
	cd ${PSDIR}; make rundir                 #^CMP IF PS
	cd ${PWDIR}; make rundir                 #^CMP IF PW
	cd ${RBDIR}; make rundir                 #^CMP IF RB
	cd ${SCDIR}; make rundir                 #^CMP IF SC
	cd ${SPDIR}; make rundir                 #^CMP IF SP
	cd ${UADIR}; make rundir                 #^CMP IF UA
	@touch CON/Scripts/${OS}/TMP_${MACHINE}
	cp CON/Scripts/${OS}/*${MACHINE}* ${RUNDIR}/
	@rm -rf ${RUNDIR}/TMP_${MACHINE} CON/Scripts/${OS}/TMP_${MACHINE}
	@echo
	@echo Creation of ${RUNDIR} directory succeeded
	@echo

CDATE = `date +%Y%b%d`

rundir_code:
	make rundir
	./Config.pl -show > _config_show
	mkdir code_${CDATE}
	rsync -a --exclude '*.o' --exclude '*.a' --exclude '*.exe' --exclude 'run*' --exclude '*~' --exclude '*.mod' --exclude code_${CDATE} . code_${CDATE}/
	tar -czf ${RUNDIR}/code_${CDATE}.tgz code_${CDATE}
	rm -rf _config_show code_${CDATE}

#
#	Run the default code on NP processors
#
NP=2

mpirun: ENV_CHECK SWMF ${RUNDIR}
	cd ${RUNDIR}; mpirun -np ${NP} ./SWMF.exe

mprun: ENV_CHECK SWMF ${RUNDIR}
	cd ${RUNDIR}; mprun -np ${NP} ./SWMF.exe

nompirun: ENV_CHECK SWMF ${RUNDIR}
	cd ${RUNDIR}; ./SWMF.exe

ETAGS = etags

tags:	ENV_CHECK
	-$(ETAGS) ./*/*/*/*.[fF]90 ./*/*/*/*.[fF] ./*/*/*/*.for

#^CMP IF IH BEGIN
#
#	Copy and rename GM/BATSRUS/src into IH/BATSRUS/src
#

IH/BATSRUS/src/Makefile:
	rm -rf IH/BATSRUS/src IH/BATSRUS/srcUser IH/BATSRUS/srcInterface/*.f90
	mkdir -p IH/BATSRUS/src IH/BATSRUS/srcUser
	cd GM/BATSRUS/src; cp *.f90 *.h Makefile* ../../../IH/BATSRUS/src
	cd IH/BATSRUS/src; rm -f main.f90
	cd GM/BATSRUS/srcInterface/; \
		cp ModGridDescriptor.f90 ModBuffer.f90 \
		update_lagrangian_grid.f90 \
		ModRadioWaveImage.f90 ModRadioWaveRaytracing.f90 \
		ModDensityAndGradient.f90\
		../../../IH/BATSRUS/srcInterface
	cp IH/BATSRUS_share/src/IH_*.f90 IH/BATSRUS/srcInterface
	cp GM/BATSRUS/srcUser/*.f90 IH/BATSRUS/srcUser/
	cd GM/BATSRUS; \
		cp Makefile.def Makefile.conf PARAM.XML Config.pl \
			../../IH/BATSRUS/
	echo '*' > IH/BATSRUS/src/.cvsignore

# rename IH source files to avoid name conflicts
IH_SRC = src/*.f90 src/*.h srcInterface/*.f90 srcUser/*.f90

IHBATSRUS: IH/BATSRUS/src/Makefile \
		${SCRIPTDIR}/Methods.pl ${SCRIPTDIR}/Rename.pl
	cd IH/BATSRUS; \
		${SCRIPTDIR}/Methods.pl IH ${IH_SRC}; \
		${SCRIPTDIR}/Rename.pl -w -r -common=IH ${IH_SRC}
	cd IH/BATSRUS/srcInterface; \
		perl -i -pe 's?BATSRUS?IH_BATSRUS?' IH_*.f90; \
		touch Makefile.DEPEND
	cd IH/BATSRUS; \
		perl -i -pe 's/GM/IH/' Config.pl; \
		./Config.pl -u=Ih -e=Mhd

#^CMP END IH
#^CMP IF SC BEGIN
#
# configure and collect source files for SC/BATSRUS component
#
SC/BATSRUS/src/Makefile:
	rm -rf \
	SC/BATSRUS/src \
	SC/BATSRUS/srcUser \
	SC/BATSRUS/srcInterface/SC_wrapper.f90 \
	SC/BATSRUS/srcInterface/SC_get_for_sp.f90
	cd GM/BATSRUS; \
		cp -f Makefile.conf ../../SC/BATSRUS; \
		make COMP=SC DREL=TMP relax_src
	cd GM/BATSRUS/TMP; \
		mv Config.pl PARAM.XML src srcUser \
			../../../SC/BATSRUS;\
		mv srcInterface/*.f90 ../../../SC/BATSRUS/srcInterface
	rm -rf GM/BATSRUS/TMP
	cp -f IH/BATSRUS_share/src/IH_wrapper.f90 \
		SC/BATSRUS/srcInterface/SC_wrapper.f90
	cp -f IH/BATSRUS_share/src/IH_get_for_sp.f90 \
		SC/BATSRUS/srcInterface/SC_get_for_sp.f90
	cd SC/BATSRUS/srcInterface/; perl -i -pe \
	's/IH/SC/g;s/BATSRUS/SC_BATSRUS/;s/Inner/Solar/;s/Heliosphere/Corona/'\
		SC_wrapper.f90 SC_get_for_sp.f90
	cd SC/BATSRUS/src; rm -f main.f90

# rename SC source files to avoid name conflicts
SC_SRC = src/*.f90 src/*.h srcInterface/*.f90 srcUser/*.f90

SCBATSRUS: SC/BATSRUS/src/Makefile \
		${SCRIPTDIR}/Methods.pl ${SCRIPTDIR}/Rename.pl
	cd SC/BATSRUS; \
		${SCRIPTDIR}/Methods.pl SC ${SC_SRC}; \
		${SCRIPTDIR}/Rename.pl -w -r -common=SC ${SC_SRC}
	touch SC/BATSRUS/srcInterface/Makefile.DEPEND
	cd SC/BATSRUS; \
		perl -i -pe 's/GM/SC/' Config.pl; \
		./Config.pl -install=c -u=Sc -e=MhdCorona

#^CMP END SC

include Makefile.test #^CMP IF TESTING
