#  Copyright (C) 2002 Regents of the University of Michigan,
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# A new line to test the email notification

#!QUOTE: \clearpage
#BOP
#!QUOTE: \section{Makefiles}
#!ROUTINE: Makefile - Space Weather Modeling Framework (SWMF) main Makefile
#!DESCRIPTION:
# This is the main Makefile to (re)compile executables, run the SWMF, 
# produce documentation, make distribution clean distribution, etc.
#EOP

SHELL=/bin/sh

# Allow setting CONFIG_PL="./Config.pl -verbose" if desired
CONFIG_PL = ./Config.pl

#
# The default target is SWMF so it is listed first
#
default : SWMF

# Nothing should be done parallel in this Makefile
.NOTPARALLEL:

#
# Make all common used executables.
#
ALL : SWMF PIDL EARTH_TRAJ

#
#       Definition of OS, component versions and directory structure
#
include Makefile.def
include Makefile.conf

help:
	@echo ' '
	@echo '  You can "make" the following:'
	@echo ' '
	@echo '    [default]   SWMF'
	@echo ' '
	@echo '    help        (makefile option list)'
	@echo '    info        (count lines of code)'
	@echo '    install     (to be used via Config.pl only)'
	@#^CMP IF NOT REMOVEDOCTEX BEGIN
	@echo ' '
	@echo '    MANUAL      (make all manuals)'          #^CMP IF DOC
	@echo '    PDF         (make PDF manuals for SWMF)' #^CMP IF DOC
	@#^CMP END  REMOVEDOCTEX
	@#^CFG IF TESTING BEGIN
	@echo '    test        (run all tests for the SWMF)'
	@echo '    test_help   (show all options for running the tests)'
	@echo ' '
	@#^CFG END TESTING
	@echo ' '
	@echo '    LIB         (lib/libSWMF.a the SWMF library)'
	@echo '    NOMPI       (lib/libNOMPI.a for single node execution with no MPI)'
	@echo ' '
	@echo '    SWMF        (bin/SWMF.exe the main executable for SWMF)'
	@#^CMP IF IE BEGIN
	@echo '    PIONO       (bin/PostIONO.exe creates ionosphere tec file from idl files)'
	@#^CMP END IE
	@#^CMP IF GM BEGIN
	@echo '    PIDL        (bin/PostIDL.exe post-processes *.idl files)'
	@echo '    SNAPSHOT    (SNAPSHOT.exe extract snapshots from *.outs files)'
	@echo '    INTERPOLATE (INTERPOLATE.exe interpolate from *.outs files)'
	@echo '    EARTH_TRAJ  (bin/EARTH_TRAJ.exe generates Earth trajectory)'
	@echo '    ALL         (make SWMF PIDL EARTH_TRAJ)'
	@#^CMP END GM
	@echo ' '
	@echo ' '
	@echo '    rundir      (create run directory "run/")'
	@echo '    rundir RUNDIR=`pwd`/run_test  (create run directory "run_test/")'
	@echo '    rundir MACHINE=frontera (select job scripts for frontera)'
	@echo '    rundir_code (saves code archive in run directory)'
	@echo ' '
	@echo '    parallelrun (make SWMF and run SWMF.exe on ${NP} PEs)'
	@echo '    parallelrun NP=7 (make SWMF and run SWMF.exe on 7 PEs)'
	@echo '    serialrun   (make SWMF and run it on 1 processor)'
	@echo ' '
	@echo '    clean       (remove files: *~ *.o *.kmo *.mod *.T *.lst core)'
	@echo '    cleanclones (remove src/Makefile-s from EE,IH,OH,SC/BATSRUS clones)'
	@echo '    distclean   (equivalent to: Config.pl -uninstall)'
	@echo '    dist        (create source distribution tar file)'
	@echo ' '
	@echo '    tags        (create etags for emacs for easy look up in source code)'
	@echo ' '
	@echo '    FORMATF90   (format F90 code properly and report changes)'
	@echo '    TDSETUP  FITS=magnetogram.fits (make directory to run TDSETUP.py'
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
# install for the 1st time
#

install: ENV_CHECK mkdir
	@echo VERSION ${VERSION}
	cd CON;	make install
	@if([ -n "${ESMF_DIR}" ]); then cd ESMF/ESMF_SWMF; make install; fi
	@if([ -d "GM/BATSRUS" ]); then \
		if([ -d "EE/BATSRUS" ]); \
			then cp GM/BATSRUS/Config.pl EE/BATSRUS; \
			     perl -i -pe 's/GM/EE/' EE/BATSRUS/Config.pl; \
		fi; \
		if([ -d "SC/BATSRUS" ]); \
			then cp GM/BATSRUS/Config.pl SC/BATSRUS; \
			     perl -i -pe 's/GM/SC/' SC/BATSRUS/Config.pl; \
		fi; \
		if([ -d "OH/BATSRUS" ]); \
			then cp GM/BATSRUS/Config.pl OH/BATSRUS; \
			     perl -i -pe 's/GM/OH/' OH/BATSRUS/Config.pl; \
		fi; \
		if([ -d "IH/BATSRUS" ]); \
			then cp GM/BATSRUS/Config.pl IH/BATSRUS; \
			     perl -i -pe 's/GM/IH/' IH/BATSRUS/Config.pl; \
		fi; \
	fi
	@for i in `ls -d [A-Z][A-Z]/*/ | grep -v /Empty/`; \
		do (if([ -f "$$i/Config.pl" ]); then \
			echo Installing $$i; cd $$i; ${CONFIG_PL} -install=c; \
		    fi); \
		done
	@echo
	@echo Installation succeeded
	@echo

info:
	@echo "NOTE: HYPRE, AMREX, BATSRUS clones, PT/FLEKS, PC/AMPS, PC/ALTOR/srcBATL* should not be counted!"
	@echo "Total lines of Fortran: `wc -l */*/src*/*.f* */*/src*/*/*.f* | tail -1`"
	@echo "Total lines of C++    : `wc -l */*/src*/*.c* */*/src*/*.h */*/src*/*/*.c* */*/src*/*/*.h | tail -1`"
	@echo "share/Library/        : `wc -l share/Library/src/*.f* share/Library/test/*.f* | tail -1`"
	@echo "util                  : `wc -l util/*/src*/*.f* | tail -1`"
	@echo "DATAREAD,NOMPI,TIMING : `wc -l util/DATAREAD/src*/*.f* util/NOMPI/src*/*.f* util/TIMING/src*/*.f* | tail -1`"
	@echo "CRASH                 : `wc -l util/CRASH/src*/*.f* | tail -1`"
	@echo "FISHPAK               : `wc -l util/FISHPAK/src*/*.f* | tail -1`"
	@echo "util/EMPIRICIAL/src*  : `wc -l util/EMPIRICAL/src*/*.f* | tail -1`"
	@echo "CON/*	             : `wc -l CON/*/src/*.f* | tail -1`"
	@echo "CON/Interface         : `wc -l CON/Interface/src/*.f* | tail -1`"
	@echo "srcInterface+wrappers : `wc -l ??/*/srcInterface/*.[fch]* ??/*/src/??_wrapper.f90 | tail -1`"
	@echo "Makefile-s            : `wc -l Makefile* */[mM]akefile* */*/[mM]akefile* */*/*/[mM]akefile* */*/*/*/[mM]akefile* */*/*/*/*/[mM]akefile* PT/AMPS/MakefileTest/* | tail -1`"
	@echo "Perl scripts          : `wc -l *.pl */*.pl */*/*.pl */*/*/*.pl */*/*/*/*.pl */*/*/*/*/*.pl | tail -1`"
	@echo "Python scripts        : `wc -l */*/*.py */*/*/*.py */*/*/*/*.py */*/*/*/*/*.py | tail -1`"
	@echo "Julia scripts         : `wc -l */*/*.jl */*/*/*.jl */*/*/*/*.jl */*/*/*/*/*.jl | tail -1`"
	@echo "IDL scripts           : `wc -l */*/*/*.pro */*/*/*/*.pro | tail -1`"
	@echo "Latex documentation   : `wc -l */*/*.tex */*/*/*.tex */*/*/*/*.tex | tail -1`"
	@echo "XML descriptions      : `wc -l PARAM.XML */*/*.XML | tail -1`"

mkdir: ./lib ./bin

./lib:
	-@mkdir lib

./bin:
	-@mkdir bin

rmdir:
	-@rm -rf lib
	-@rm -rf bin

GITINFO:
	@(if [ "${GITINFO}" != "NO" ]; then \
		${SCRIPTDIR}/gitall -r=f > ${CONTROLDIR}/show_git_info.h; \
	fi)

#
#	SWMF
#
SWMF:	ENV_CHECK GITINFO
	@cd ${CONTROLDIR}; $(MAKE) SWMFEXE  
	@echo ' '

#
#       SWMF library (for external code)
#
LIB:	ENV_CHECK
	@cd ${CONTROLDIR}; $(MAKE) SWMFLIB
	@echo ' '

#       ESMF driver for the SWMF and another ESMF component
#
ESMF_SWMF: LIB
	@cd ESMF/ESMF_SWMF/src; $(MAKE) EXE

# NOMPI library for execution without MPI

NOMPI: ENV_CHECK
	cd ${NOMPIDIR}; $(MAKE) LIB

#^CMP IF GM BEGIN
#	Post processing codes for BATSRUS plot files
#
PIDL:	ENV_CHECK
	cd ${SHAREDIR}; $(MAKE) PIDL
	@echo ' '

FDIPS:	ENV_CHECK
	cd ${MAGNETOGRAMDIR}; ${MAKE} FDIPS
	@echo ' '

# Convert .out data to .dat 
2TEC:	ENV_CHECK
	cd ${SHAREDIR}; $(MAKE) 2TEC
	@echo ' '

# Convert .out data to .vtk
2VTK:	ENV_CHECK
	cd ${SHAREDIR}; $(MAKE) 2VTK
	@echo ' '

SNAPSHOT: ENV_CHECK
	cd GM/BATSRUS; $(MAKE) SNAPSHOT
	@echo ' '

INTERPOLATE: ENV_CHECK
	cd GM/BATSRUS; $(MAKE) INTERPOLATE
	@echo ' '

#
#	Code for generating Earth trajectory for IH/BATSRUS
#
EARTH_TRAJ: ENV_CHECK
	cd GM/BATSRUS; $(MAKE) EARTH_TRAJ
	@echo ' '
#^CMP END GM

#^CMP IF IE BEGIN
#	Post processing code for IE/Ridley_serial plot files
#
PIONO:	ENV_CHECK
	cd IE/Ridley_serial; $(MAKE) PIONO
	@echo ' '

#^CMP END IE

#					^CMP IF DOC BEGIN
#	Create the documentation files      ^CMP IF NOT REMOVEDOCTEX BEGIN
#	
MANUAL: ENV_CHECK
	@cd doc/Tex; make cleanpdf; make PDF
	@if([ -d "util/CRASH" ]); then cd util/CRASH/doc/Tex; make PDF; fi
	@if([ -d "GM/BATSRUS" ]); then cd GM/BATSRUS; make PDF; fi #^CMP IF GM
	@if([ -d "PW/PWOM"    ]); then cd PW/PWOM;    make PDF; fi #^CMP IF PW
	@if([ -d "IM/CIMI "   ]); then cd IM/CIMI;    make PDF; fi #^CMP IF IM
	@if([ -d "PT/AMPS"    ]); then cd PT/AMPS;    make PDF; fi #^CMP IF PT

PDF:	ENV_CHECK
	@cd doc/Tex; make cleanpdf; make PDF

CLEAN1 = cleanpdf #				^CMP IF NOT MAKEPDF

#					    ^CMP END REMOVEDOCTEX
#					^CMP END DOC

#			
#	General Housekeeping
#
clean: ENV_CHECK
	@echo
	rm -rf *~ doc/*~ Param/*~ TAGS
	for i in `ls -d [A-Z][A-Z]/*/`; \
		do (echo Cleaning $$i; cd $$i; make clean); done
	-@if([ -n "${ESMF_DIR}" ]); then cd ESMF/ESMF_SWMF; make clean; fi
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

distclean: 
	${CONFIG_PL} -uninstall

allclean: ENV_CHECK rmdir
	@echo
	rm -rf *~ doc/*~ Param/*~ TAGS
	for i in `ls -d [A-Z][A-Z]/Empty/`; \
		do (echo Distcleaning $$i; cd $$i; make distclean); done
	for i in `ls -d [A-Z][A-Z]/*/ | grep -v Empty`; \
		do (echo Uninstalling $$i; cd $$i; ${CONFIG_PL} -uninstall); \
		done
	@if([ -n "${ESMF_DIR}" ]); then cd ESMF/ESMF_SWMF;make distclean; fi 
	cd CON;			make distclean
	@#^CMP IF DOC BEGIN
	@#^CMP IF NOT REMOVEDOCTEX BEGIN
	cd doc/Tex;             make clean ${CLEAN1}
	@#^CMP END REMOVEDOCTEX
	@#^CMP END DOC
	@echo
	@echo Uninstallation/distclean succeeded
	@echo

cleanclones: ENV_CHECK
	mv GM/BATSRUS/src/Makefile GM/BATSRUS/src/Makefile.safe
	rm -f ??/BATSRUS/src/Makefile
	mv GM/BATSRUS/src/Makefile.safe GM/BATSRUS/src/Makefile

TAR = tar --exclude .git -rf tmp.tar

dist:
	@echo ' '
	@echo ' NOTE: All "run" or other created directories not included!'
	@echo ' '
	${CONFIG_PL} -uninstall
	-rm -rf */*/run_test
	tar -cf tmp.tar README.md PARAM.XML LICENSE.txt Config.pl Makefile
	${TAR} Makefile.test output		#^CMP IF TESTING
	${TAR} doc				#^CMP IF DOC
	${TAR} Copyrights Param Scripts CON share util ESMF
	if([ -d gui ]); then ${TAR} gui; fi
	for i in `ls -d [A-Z][A-Z]`; \
		do (echo Tarring $$i; ${TAR} $$i); done
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
	@echo "Copy a working PARAM.in file here" > ${RUNDIR}/PARAM.in
	cp share/Scripts/Restart.pl share/Scripts/PostProc.pl share/Scripts/Resubmit.pl ${RUNDIR}/
	touch ${RUNDIR}/core
	chmod 444 ${RUNDIR}/core
	cd ${RUNDIR}; \
		ln -s ${DIR}/bin/SWMF.exe . ; \
		ln -s ${DIR}/Param . ; \
		ln -s ${DIR} SWMFDIR
	@for i in `ls -d [A-Z][A-Z]`; do \
		make rundircomp COMPDIR=$${i}DIR; \
	done
	@touch  share/JobScripts/job._TMP_${MACHINE} \
		share/JobScripts/_TMP_.${MACHINE}.pl
	@cp share/JobScripts/job.*${MACHINE}* \
	    share/JobScripts/*.${MACHINE}*.pl       ${RUNDIR}/
	@echo "cp share/JobScripts/job.*${MACHINE}* ${RUNDIR}/"
	@rm -rf ${RUNDIR}/*_TMP_* share/JobScripts/*_TMP_*
	@echo
	@echo Creation of ${RUNDIR} directory succeeded
	@echo

rundircomp:
	cd ${${COMPDIR}}; make rundir

CDATE = `date +%Y%b%d`

rundir_code:
	make rundir
	${CONFIG_PL} -show > _config_show
	mkdir code_${CDATE}
	rsync -a --exclude '*.o' --exclude '*.a' --exclude '*.exe' --exclude 'run*' --exclude '*~' --exclude '*.mod' --exclude code_${CDATE} . code_${CDATE}/
	tar -czf ${RUNDIR}/code_${CDATE}.tgz code_${CDATE}
	rm -rf _config_show code_${CDATE}

#
#	Run the default code on NP processors
#
NP=2

parallelrun: ENV_CHECK SWMF ${RUNDIR}
	cd ${RUNDIR}; ${MPIRUN} ./SWMF.exe

serialrun: ENV_CHECK SWMF ${RUNDIR}
	cd ${RUNDIR}; ${SERIAL} ./SWMF.exe

ETAGS = etags

tags:	ENV_CHECK
	-$(ETAGS) ./*/*/*/*.[fF]90 ./*/*/*/*.[fF] ./*/*/*/*.for


FORMATF90:
	-@share/Scripts/FormatFortran.pl CON/*/src/*.f90
	-@share/Scripts/FormatFortran.pl -l GM/BATSRUS/src/*.f90
	-@share/Scripts/FormatFortran.pl GM/BATSRUS/src[A-Z]*/*.f90
	-@share/Scripts/FormatFortran.pl share/Library/src/*.f90

#^CMP IF EE BEGIN
#
# collect source files for EE/BATSRUS component
#
EE/BATSRUS/src/Makefile:
	cd EE/BATSRUS; \
		rm -rf src srcBATL srcUser srcEquation; \
		mkdir src srcBATL srcUser srcEquation
	cd GM/BATSRUS/src; cp *.f90 *.h Makefile* ../../../EE/BATSRUS/src
	cd GM/BATSRUS/srcBATL; cp *.f90 Makefile* ../../../EE/BATSRUS/srcBATL
	cd GM/BATSRUS/srcInterface/; \
		cp ModGridDescriptor.f90 \
		../../../EE/BATSRUS/srcInterface
	cp GM/BATSRUS/srcUser/*.f90 EE/BATSRUS/srcUser/
	if [ -d GM/BATSRUS/srcUserExtra ]; then \
		cp GM/BATSRUS/srcUserExtra/ModUserEe*.f90 \
		   GM/BATSRUS/srcUserExtra/ModUserSwarm*.f90 \
			EE/BATSRUS/srcUser/; \
	fi
	cp GM/BATSRUS/srcEquation/*.f90 EE/BATSRUS/srcEquation/
	cd GM/BATSRUS; \
		cp Makefile.def Makefile.conf PARAM.XML Config.pl \
			../../EE/BATSRUS/
	cd EE/BATSRUS/src; rm -f main.f90

# rename EE source files to avoid name conflicts
EE_SRC = src/*.f90 src/*.h srcBATL/*.f90 srcUser/*.f90 srcEquation/*.f90 \
	srcInterface/*.f90

EEBATSRUS: EE/BATSRUS/src/Makefile \
		${SCRIPTDIR}/Methods.pl ${SCRIPTDIR}/Rename.pl
	cd EE/BATSRUS; \
		${SCRIPTDIR}/Methods.pl EE ${EE_SRC}; \
		${SCRIPTDIR}/Rename.pl -w -r -common=EE ${EE_SRC}; \
		touch srcInterface/Makefile.DEPEND; \
		perl -i -pe 's/GM/EE/' Config.pl; \
		${CONFIG_PL} -install=c -u=Ee -e=MhdEosRad

#^CMP END EE


#^CMP IF IH BEGIN
#
#	Copy and rename GM/BATSRUS/src into IH/BATSRUS/src
#

IH/BATSRUS/src/Makefile:
	cd IH/BATSRUS; \
		rm -rf src srcBATL srcUser srcEquation srcInterface/Mod*.f90; \
		mkdir src srcBATL srcUser srcEquation
	cd GM/BATSRUS/src; cp *.f90 *.h Makefile* ../../../IH/BATSRUS/src
	cd GM/BATSRUS/srcBATL; \
		cp BATL*.f90 Makefile* ../../../IH/BATSRUS/srcBATL
	cd GM/BATSRUS/srcInterface/; \
		cp ModGridDescriptor.f90 ../../../IH/BATSRUS/srcInterface
	cp GM/BATSRUS/srcUser/*.f90 IH/BATSRUS/srcUser/
	if [ -d GM/BATSRUS/srcUserExtra ]; then \
		cp GM/BATSRUS/srcUserExtra/ModUserAwsom*.f90 \
		   IH/BATSRUS/srcUser/; \
	fi
	cp GM/BATSRUS/srcEquation/*.f90 IH/BATSRUS/srcEquation/
	cd GM/BATSRUS; \
		cp Makefile.def Makefile.conf PARAM.XML Config.pl \
		../../IH/BATSRUS/
	cd IH/BATSRUS/src; rm -f main.f90

# rename IH source files to avoid name conflicts
IH_SRC = src/*.f90 src/*.h srcBATL/*.f90 srcUser/*.f90 srcEquation/*.f90 \
	srcInterface/*.f90

IHBATSRUS: IH/BATSRUS/src/Makefile \
		${SCRIPTDIR}/Methods.pl ${SCRIPTDIR}/Rename.pl
	cd IH/BATSRUS; \
		${SCRIPTDIR}/Methods.pl IH ${IH_SRC}; \
		${SCRIPTDIR}/Rename.pl -w -r -common=IH ${IH_SRC}
	cd IH/BATSRUS/srcInterface; \
		touch Makefile.DEPEND
	cd IH/BATSRUS; \
		perl -i -pe 's/GM/IH/' Config.pl; \
		${CONFIG_PL} -install=c -u=Awsom -e=Awsom

#^CMP END IH


#^CMP IF OH BEGIN
#
# configure and collect source files for OH/BATSRUS component
#
OH/BATSRUS/src/Makefile:
	cd OH/BATSRUS; \
		rm -rf src srcBATL srcUser srcEquation \
			srcInterface/OH_wrapper.f90; \
		mkdir src srcBATL srcUser srcEquation
	cd GM/BATSRUS/src; cp *.f90 *.h Makefile* ../../../OH/BATSRUS/src
	cd GM/BATSRUS/srcBATL; cp *.f90 Makefile* ../../../OH/BATSRUS/srcBATL
	cd GM/BATSRUS/srcInterface/; \
		cp ModGridDescriptor.f90 \
		../../../OH/BATSRUS/srcInterface/
	cp GM/BATSRUS/srcUser/*.f90 OH/BATSRUS/srcUser/
	cp GM/BATSRUS/srcEquation/*.f90 OH/BATSRUS/srcEquation/
	cd GM/BATSRUS; \
		cp Makefile.def Makefile.conf PARAM.XML Config.pl \
			../../OH/BATSRUS/
	cd OH/BATSRUS/src; rm -f main.f90
	cp -f IH/BATSRUS/srcInterface/IH_wrapper.f90 \
		OH/BATSRUS/srcInterface/OH_wrapper.f90

	cd OH/BATSRUS/srcInterface/; perl -i -pe \
	's/IH/OH/g;s/Ih/Oh/g;s/_sc\b/_ih/;s/SC/IH/g;s/Sc/Ih/g;s/Inner/Outer/' \
		OH_wrapper.f90

# rename OH source files to avoid name conflicts
OH_SRC = src/*.f90 src/*.h srcBATL/*.f90 srcUser/*.f90 srcEquation/*.f90 \
	srcInterface/*.f90

OHBATSRUS: OH/BATSRUS/src/Makefile \
		${SCRIPTDIR}/Methods.pl ${SCRIPTDIR}/Rename.pl
	cd OH/BATSRUS; \
		${SCRIPTDIR}/Methods.pl OH ${OH_SRC}; \
		${SCRIPTDIR}/Rename.pl -w -r -common=OH ${OH_SRC}
	touch OH/BATSRUS/srcInterface/Makefile.DEPEND
	cd OH/BATSRUS; \
		perl -i -pe 's/GM/OH/' Config.pl; \
		${CONFIG_PL} -install=c -e=OuterHelio -u=OuterHelio


#^CMP END OH


#^CMP IF SC BEGIN
#
# configure and collect source files for SC/BATSRUS component
#
SC/BATSRUS/src/Makefile:
	cd SC/BATSRUS; \
		rm -rf src srcBATL srcUser srcEquation \
		       srcInterface/SC_wrapper.f90; \
		mkdir src srcBATL srcUser srcEquation
	cd GM/BATSRUS/src; cp *.f90 *.h Makefile* ../../../SC/BATSRUS/src
	cd GM/BATSRUS/srcBATL; cp *.f90 Makefile* ../../../SC/BATSRUS/srcBATL
	cp GM/BATSRUS/srcUser/*.f90 SC/BATSRUS/srcUser/
	if [ -d GM/BATSRUS/srcUserExtra ]; then \
		cp GM/BATSRUS/srcUserExtra/ModUserAwsom*.f90 \
		   SC/BATSRUS/srcUser/; \
	fi
	cp GM/BATSRUS/srcEquation/*.f90 SC/BATSRUS/srcEquation/
	cd GM/BATSRUS; \
		cp Makefile.def Makefile.conf PARAM.XML Config.pl \
			../../SC/BATSRUS/
	cd GM/BATSRUS/srcInterface/; \
		cp ModGridDescriptor.f90 \
		../../../SC/BATSRUS/srcInterface
	cp -f IH/BATSRUS/srcInterface/IH_wrapper.f90 \
		SC/BATSRUS/srcInterface/SC_wrapper.f90
	cd SC/BATSRUS/srcInterface/; \
	perl -i -pe \
		's/IH/SC/g;s/OH/IH/;s/Inner Heliosphere/Solar Corona/' \
		SC_wrapper.f90
	cd SC/BATSRUS/src; rm -f main.f90

# rename SC source files to avoid name conflicts
SC_SRC = src/*.f90 src/*.h srcBATL/*.f90 srcUser/*.f90 srcEquation/*.f90 \
	srcInterface/*.f90

SCBATSRUS: SC/BATSRUS/src/Makefile \
		${SCRIPTDIR}/Methods.pl ${SCRIPTDIR}/Rename.pl
	cd SC/BATSRUS; \
		${SCRIPTDIR}/Methods.pl SC ${SC_SRC}; \
		${SCRIPTDIR}/Rename.pl -w -r -common=SC ${SC_SRC}
	touch SC/BATSRUS/srcInterface/Makefile.DEPEND
	cd SC/BATSRUS; \
		perl -i -pe 's/GM/SC/' Config.pl; \
		${CONFIG_PL} -install=c -u=Awsom -e=Awsom

FITS=${SCDIR}/data/input/2135_026.fits
# Creates environment to run setup for TD erruptive event generator
# It is also useful for EEGGL
# To handle the input magnetogram_file.fits in directory run_td, run:
# make TDSETUP FITS=magnetogram_file.fits RUNDIR=`pwd`/run_td
# If more than two processorrs are available, use MPIRUN='mpiexec -n 4',
# to speed up the process.
# NOTE: the run directtory will be erased and creaated from scratch. To
# avoid this (say, if only the magnetogram changes), use instead:
# make td_setup FITS=new_magnetogram_file.fits RUNDIR=`pwd`/run_td

TDSETUP:
	${MAKE} td_compile
	${MAKE} td_rundir
	${MAKE} td_setup

td_compile:
	cd ${EMPIRICALEEDIR}; \
		${MAKE} FRM;  \
		${MAKE} CMEBR
	cd ${MAGNETOGRAMDIR}; \
		${MAKE} HARMONICS; \
		${MAKE} CONVERTHARMONICS; \
		${MAKE} FDIPS

td_rundir:
	./Config.pl -v=SC/BATSRUS
	@echo "Directory ${RUNDIR} is erased!!!"
	rm -rf ${RUNDIR}
	$(MAKE) rundir
	cd ${RUNDIR}/SC; \
	perl -i -pe 's/dipole11uniform/fitsfile_01/; s/harmonics11uniform/harmonics/' \
	     HARMONICS.in; \
	perl -i -pe 's/\d+(\s+MaxOrder)/180$$1/' \
	     HARMONICS.in; \
	perl -i -pe 's/\d+(\s+MaxOrder)/180$$1/; s/\d+(\s+nR)/100$$1/' \
	     HARMONICSGRID.in;  \
	perl -i -pe 's/\d+(\s+nLon)/360$$1/; s/\d+(\s+nLat)/180$$1/' \
	     HARMONICSGRID.in
	@echo " Before stating td_setup, download the magnetogram"

td_setup:
	@echo "Copy fits file ${FITS} to  ${RUNDIR}/SC/fitsfile.fits"
	cp -f ${FITS} ${RUNDIR}/SC/fitsfile.fits
	cd ${RUNDIR}/SC; \
	python3 remap_magnetogram.py fitsfile.fits fitsfile; \
	./HARMONICS.exe |tee harmonics.log
	@echo " Reconstruction of potential field takes a while, please, wait"
	cd ${RUNDIR}/SC; \
	${MPIRUN} ./CONVERTHARMONICS.exe > convert.log
	@echo "-----------------------------------"
	@echo "To start TDSETUP do the following:"
	@echo "cd ${RUNDIR}/SC"
	@echo "python3 TDSETUP.py field_2d.out"
	@echo "-----------------------------------"
	@echo "To start GLSETUP the second command should read:"
	@echo "python3 GLSETUP.py field_2d.out -CMESpeed 600"
#^CMP END SC

include Makefile.test #^CMP IF TESTING
