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
VERSION = 2.2

#
# The default target is SWMF so it is listed first
#
default : SWMF

#
# Make all common used executables.
#
ALL : SWMF PIDL PSPH

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
	@echo '    ALL         SWMF PIDL PSPH'
	@echo ' '
	@echo '    install     (to be used via SetSWMF.pl only)'
	@echo ' '
	@echo '    SWMF        (bin/SWMF.exe the main executable for SWMF)'
	@echo '    PIDL        (bin/PostIDL.exe creates .out file from *.idl files)'
	@echo '    PSPH        (bin/PostSPH.exe creates spherical tec file from sph*.tec files)'
	@echo '    PIONO       (bin/PostIONO.exe creates ionosphere tec file from idl files)'
	@echo '    LIB         (lib/libSWMF.a the SWMF library)'
	@echo '    NOMPI       (lib/libNOMPI.a for single node execution with no MPI)'
	@echo ' '
	@echo '    help        (makefile option list)'
	@echo '    clean       (remove files: *~ *.o *.kmo *.mod *.T *.lst core)'
	@echo '    distlean    (remove all files which are not in distribution)'
	@echo '    dist        (create source distribution tar file)'
	@#^CMP IF NOT REMOVEDOCTEX BEGIN
	@echo ' '
	@echo '    PDF         (make PDF  version of the documentation)' #^CMP IF DOC
	@echo '    HTML        (make HTML version of the documentation)' #^CMP IF DOCHTML
	@#^CMP END  REMOVEDOCTEX
	@echo ' '
	@echo '    rundir      (create run directory with subdirectories scripts and links)'
	@echo '    rundir MACHINE=ames (run directory with job scripts for machines at NASA Ames)'
	@echo ' '
	@echo '    mpirun      (make SWMF and mpirun SWMF.exe on 2 PEs)'
	@echo '    mpirun NP=7 (make SWMF and mpirun SWMF.exe on 7 PEs)'
	@echo '    mprun  NP=5 (make SWMF and mprun  SWMF.exe on 5 PEs)'
	@echo '    nompirun    (make SWMF and run it without MPI)'
	@echo ' '
	@echo '    tags        (create etags for emacs for easy look up in source code)'
	@#^CMP IF IH BEGIN
	@echo ' '
	@echo '    IHBATSRUS   (copy and rename GM/BATSRUS source into IH/BATSRUS)'
	@#^CMP END IH
	@#^CMP IF SC BEGIN
	@echo ' '
	@echo '    SCBATSRUS   (configure and rename GM/BATSRUS source into SC/BATSRUS)'
	@#^CMP END SC
	@echo ' '
	@echo '    ESMF_SWMF          (bin/ESMF_SWMF.exe the ESMF compatible executable)'
	@echo '    ESMF_SWMF_test     (build and run ESMF_SWMF.exe on 2 CPU-s)'
	@echo '    ESMF_SWMF_rundir   (create run directory for ESMF_SWMF.exe)'
	@echo '    ESMF_SWMF_run NP=4 (run ESMF_SWMF.exe on 4 CPU-s)'
#EOC
#
# Check the variables SWMF_ROOT=`pwd` and OS = `uname`
#
ENV_CHECK:
	@(if([ "${SWMF_ROOT}" != `/bin/pwd` ]) || ([ "${OS}" != `uname` ]); \
	then \
	  echo; \
	  echo "SWMF_ROOT='${SWMF_ROOT}' and/or OS='${OS}' are incorrect!";\
	  echo "Correcting Makefile.def";\
	  perl -pi share/Scripts/FixMakefileDef.pl Makefile.def; \
	  echo "Type 'SetSWMF.pl -s' and retry the previous make command!"; \
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
	cd share;       	make install
	cd CON;			make install
	cd ESMF/ESMF_SWMF;	make install
	for i in `ls -d [A-Z][A-Z]/*/ | grep -v /CVS/`; \
		do (cd $$i; make install); done
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
#       SWMF library
#
LIB:	ENV_CHECK
	@cd ${CONTROLDIR}; make SWMFLIB
	@echo ' '

# NOMPI library for execution without MPI

NOMPI: ENV_CHECK
	cd ${NOMPIDIR}; make LIB

#
#	Post processing codes for BATSRUS plot files
#
PSPH:	ENV_CHECK				#^CMP IF GM
	cd GM/BATSRUS; make PSPH		#^CMP IF GM
	@echo ' '				#^CMP IF GM

PIDL:	ENV_CHECK				#^CMP IF GM
	cd GM/BATSRUS; make PIDL		#^CMP IF GM
	@echo ' '				#^CMP IF GM

#
#	Post processing code for IE/Ridley_serial plot files
#
PIONO:	ENV_CHECK				#^CMP IF IE
	cd IE/Ridley_serial; make PIONO		#^CMP IF IE
	@echo ' '				#^CMP IF IE

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
	cd ESMF/ESMF_SWMF;	make clean
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

distclean_comp: ENV_CHECK rmdir
	@echo
	rm -rf *~ doc/*~ Param/*~ TAGS
	for i in `ls -d [A-Z][A-Z]/*/ | grep -v /CVS/`; \
		do (cd $$i; make distclean); done
	cd ESMF/ESMF_SWMF;	make distclean
	cd CON;			make distclean
	cd util;		make distclean
	cd share;               make distclean
	@#^CMP IF DOC BEGIN
	@#^CMP IF NOT REMOVEDOCTEX BEGIN
	cd doc/Tex;             make clean ${CLEAN1} ${CLEAN2}
	@#^CMP END REMOVEDOCTEX
	@#^CMP END DOC
	rm -f Makefile.conf Makefile.def
	@echo
	@echo Distclean succeeded
	@echo

distclean:
	rm -f */BATSRUS/src/user_routines.f90
	make distclean_comp

dist: distclean
	@echo ' '
	@echo ' NOTE: All "run" or other created directories not included!'
	@echo ' '
	tar -cf tmp.tar  README
	tar -rf tmp.tar  Makefile
	tar -rf tmp.tar  Copyrights
	tar -rf tmp.tar  CVS*
	tar -rf tmp.tar  .cvsignore
	tar -rf tmp.tar  Configure.options
	tar -rf tmp.tar  Configure.pl		#^CMP IF CONFIGURE
	tar -rf tmp.tar  doc			#^CMP IF DOC
	tar -rf tmp.tar  Param
	tar -rf tmp.tar  Scripts
	tar -rf tmp.tar  SetSWMF.pl
	tar -rf tmp.tar  share
	tar -rf tmp.tar  util
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
run:
	make rundir

rundir: ENV_CHECK
	mkdir run
	mkdir run/STDOUT
	cp Param/LAYOUT.DEFAULT       run/LAYOUT.in
	cp Param/PARAM.DEFAULT        run/PARAM.in
	cp share/Scripts/Restart.pl   run/Restart.pl
	cp share/Scripts/PostProc.pl  run/PostProc.pl
	touch run/core
	chmod 444 run/core
	cd run; ln -s  ${DIR}/bin/SWMF.exe . ; ln -s  ${DIR}/Param .
	cd ${GMDIR}; make rundir                 #^CMP IF GM
	cd ${IEDIR}; make rundir                 #^CMP IF IE
	cd ${IHDIR}; make rundir                 #^CMP IF IH
	cd ${IMDIR}; make rundir                 #^CMP IF IM
	cd ${PWDIR}; make rundir                 #^CMP IF PW
	cd ${RBDIR}; make rundir                 #^CMP IF RB
	cd ${SCDIR}; make rundir                 #^CMP IF SC
	cd ${SPDIR}; make rundir                 #^CMP IF SP
	cd ${UADIR}; make rundir                 #^CMP IF UA
	@touch CON/Scripts/${OS}/TMP_${MACHINE}
	cp CON/Scripts/${OS}/*${MACHINE}* run/
	@rm -rf run/TMP_${MACHINE} CON/Scripts/${OS}/TMP_${MACHINE}
	@echo
	@echo Creation of run directory succeeded
	@echo

#
#	Run the default code on NP processors
#
NP=2

mpirun: ENV_CHECK SWMF run
	cd run; mpirun -np ${NP} ./SWMF.exe

mprun: ENV_CHECK SWMF run
	cd run; mprun -np ${NP} ./SWMF.exe

nompirun: ENV_CHECK SWMF run
	cd run; ./SWMF.exe

ETAGS = etags

tags:	ENV_CHECK
	-$(ETAGS) ./*/*/*/*.[fF]90 ./*/*/*/*.[fF] ./*/*/*/*.for

#^CMP IF IH BEGIN
#
#	Copy and rename GM/BATSRUS/src into IH/BATSRUS/src
#

IH/BATSRUS/src/Makefile:
	mkdir -p IH/BATSRUS/src IH/BATSRUS/srcUser
	cd GM/BATSRUS/src; cp *.f90 *.f *.h Makefile* ../../../IH/BATSRUS/src
	cd IH/BATSRUS/src; rm -f main.f90 stand_alone*.f90
	cd GM/BATSRUS/srcInterface/; \
		cp ModGridDescriptor.f90 ModBuffer.f90 \
		update_lagrangian_grid.f90 \
		../../../IH/BATSRUS/srcInterface
	cp IH/BATSRUS_share/src/IH_*.f90 IH/BATSRUS/srcInterface
	cp GM/BATSRUS/srcUser/*.f90 IH/BATSRUS/srcUser/
	cd GM/BATSRUS; \
		cp Makefile.def Makefile.conf PARAM.XML PARAM.pl \
			GridSize.pl Options.pl ../../IH/BATSRUS/
	echo '*' > IH/BATSRUS/src/.cvsignore

# rename IH source files to avoid name conflicts
IH_SRC = src/*.f90 src/*.h src/*.f srcInterface/*.f90 srcUser/*.f90

IHBATSRUS: IH/BATSRUS/src/Makefile \
		${SCRIPTDIR}/Methods.pl ${SCRIPTDIR}/Rename.pl
	cd IH/BATSRUS; \
		${SCRIPTDIR}/Methods.pl IH ${IH_SRC}; \
		${SCRIPTDIR}/Rename.pl -w -r -common=IH ${IH_SRC}
	cd IH/BATSRUS/srcInterface; \
		perl -i -pe 's?BATSRUS?IH_BATSRUS?' IH_*.f90; \
		touch Makefile.DEPEND
	cd IH/BATSRUS/; ./Options.pl -u=Heliosphere

#^CMP END IH
#^CMP IF SC BEGIN
#
# configure and collect source files for SC/BATSRUS component
#
SC/BATSRUS/src/Makefile:
	cd GM/BATSRUS; \
		cp -f Makefile.conf ../../SC/BATSRUS; \
		make COMP=SC DREL=TMP relax_src
	cd GM/BATSRUS/TMP; \
		mv Makefile.def GridSize.pl Options.pl PARAM.XML src srcUser \
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
	cd SC/BATSRUS/src; rm -f main.f90 stand_alone*.f90

# rename SC source files to avoid name conflicts
SC_SRC = src/*.f90 src/*.h src/*.f srcInterface/*.f90 srcUser/*.f90

SCBATSRUS: SC/BATSRUS/src/Makefile \
		${SCRIPTDIR}/Methods.pl ${SCRIPTDIR}/Rename.pl
	cd SC/BATSRUS; \
		${SCRIPTDIR}/Methods.pl SC ${SC_SRC}; \
		${SCRIPTDIR}/Rename.pl -w -r -common=SC ${SC_SRC}
	touch SC/BATSRUS/srcInterface/Makefile.DEPEND
	cd SC/BATSRUS; ./Options.pl -e=MhdCorona -u=Heliosphere

#^CMP END SC

#
#       Targets for the ESMF compatible executable
#
ESMF_SWMF: LIB
	@cd ESMF/ESMF_SWMF/src; make EXE

ESMF_SWMF_test: ESMF_SWMF
	@make ESMF_SWMF_run

ESMF_SWMF_run: bin/ESMF_SWMF.exe ESMF_SWMF_rundir
	@cd run; rm -f PET*ESMF_LogFile; mpirun -np ${NP} ESMF_SWMF.exe; \
	perl -ne 'print if /error/i' PET*.ESMF_LogFile

ESMF_SWMF_rundir: run
	@cd ESMF/ESMF_SWMF; make rundir


test: test_comp

test_comp:
	for i in `ls -d [A-Z][A-Z]/*/ | grep -v /CVS/ | grep -v /Empty/`; \
		do ( cd $$i; make test; ); done
	ls -l [A-Z][A-Z]/*/*.diff
