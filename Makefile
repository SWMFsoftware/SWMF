#^CMP COPYRIGHT UM

#BOP
#!ROUTINE: Makefile - Space Weather Modeling Framework (SWMF) main Makefile
#!DESCRIPTION:
# This is the main Makefile to (re)compile executables, run the SWMF, 
# produce documentation, make distribution clean distribution, etc.
#EOP

SHELL=/bin/sh

#BOC
VERSION = 2.0

#
# The default target is SWMF so it is listed first
#
default : SWMF

#
#       Definition of OS, component versions and directory structure
#
include Makefile.def

help:
	@echo ' '
	@echo '  You can "make" the following:'
	@echo ' '
	@echo '    <default>   SWMF'
	@echo ' '
	@echo '    install     (to be used via SetSWMF.pl only)'
	@echo ' '
	@echo '    SWMF        (bin/SWMF.exe the main executable for SWMF)'
	@echo '    PIDL        (bin/PostIDL.exe creates .out file from *.idl files)'
	@echo '    PSPH        (bin/PostSPH.exe creates spherical tec file from sph*.tec files)'
	@echo '    PIONO       (bin/PostIONO.exe creates ionosphere tec file from idl files)'
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
	@echo '    rundir PLOT=IDL     (run directory with IDL     postprocessing only)'
	@echo '    rundir PLOT=TEC     (run directory with Tecplot postprocessing only)'
	@echo '    rundir MACHINE=ames (run directory with job scripts for machines at NASA Ames)'
	@echo ' '
	@echo '    mpirun      (make SWMF and mpirun SWMF.exe on 8 PEs)'
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
	  echo "You may now try the previous make command again ..."; \
	  exit 1; \
	fi);

#
# POSTPROCESSING SCRIPTS TO COPY TO RUN/GM OR RUN/IH
#
# this variable (PLOT) holds the names of types of scripts to copy to the run 
# directory when it is created.  The options are below. Uncomment only one
# of them.  The "???" option will copy all type of postprocessing scripts.
#
#PLOT = IDL
#PLOT = TEC
PLOT = ???

#
# RUN SCRIPTS TO COPY TO RUN
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
	cd GM/BATSRUS;		make install    #^CMP IF GM
	cd IE/Ridley_serial;	make install    #^CMP IF IE
	cd IH/BATSRUS_share;	make install    #^CMP IF IH
	cd IM/RCM;		make install    #^CMP IF IM
	cd RB/RiceV5;		make install    #^CMP IF RB
	cd SP/Kota;		make install    #^CMP IF SP
	cd UA/GITM;		make install	#^CMP IF UA
	cd UA/GITM2;		make install	#^CMP IF UA
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

# NOMPI library for execution without MPI

NOMPI: ENV_CHECK
	cd ${NOMPIDIR}; make LIB
#
# STAND-ALONE EXECUTABLES FOR COMPONENTS
#

RBM:	ENV_CHECK                               #^CMP IF RB
	@cd ${RBDIR}; make RBM                  #^CMP IF RB
	@echo ' '                               #^CMP IF RB

#
#	Post processing
#
PSPH:	ENV_CHECK				#^CMP IF GM
	cd GM/BATSRUS; make PSPH		#^CMP IF GM
	@echo ' '				#^CMP IF GM

PIDL:	ENV_CHECK				#^CMP IF GM
	cd GM/BATSRUS; make PIDL		#^CMP IF GM
	@echo ' '				#^CMP IF GM

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
	cd GM/BATSRUS;		make clean    #^CMP IF GM
	cd GM/Empty;		make clean    #^CMP IF GM
	cd IE/Ridley_serial;	make clean    #^CMP IF IE
	cd IE/Empty;		make clean    #^CMP IF IE
	cd IH/BATSRUS;		make clean    #^CMP IF IH
	cd IH/BATSRUS_share;	make clean    #^CMP IF IH
	cd IH/Empty;		make clean    #^CMP IF IH
	cd IM/RCM;		make clean    #^CMP IF IM
	cd IM/Empty;		make clean    #^CMP IF IM
	cd UA/GITM;		make clean    #^CMP IF UA
	cd UA/GITM2;		make clean    #^CMP IF UA
	cd UA/Empty;		make clean    #^CMP IF UA
	cd RB/Rice;		make clean    #^CMP IF RB
	cd RB/RiceV5;		make clean    #^CMP IF RB
	cd RB/Empty;		make clean    #^CMP IF RB
	cd SP/Kota;             make clean    #^CMP IF SP
	cd SP/FLAMPA;           make clean    #^CMP IF SP
	cd SP/Empty;            make clean    #^CMP IF SC
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
	cd GM/BATSRUS;		make distclean    #^CMP IF GM
	cd GM/Empty;		make distclean    #^CMP IF GM
	cd IE/Ridley_serial;	make distclean    #^CMP IF IE
	cd IE/Empty;		make distclean    #^CMP IF IE
	cd IH/BATSRUS;		make distclean    #^CMP IF IH
	cd IH/BATSRUS_share;	make distclean    #^CMP IF IH
	cd IH/Empty;		make distclean    #^CMP IF IH
	cd IM/RCM;		make distclean    #^CMP IF IM
	cd IM/Empty;		make distclean    #^CMP IF IM
	cd UA/GITM;		make distclean    #^CMP IF UA
	cd UA/GITM2;		make distclean    #^CMP IF UA
	cd UA/Empty;		make distclean    #^CMP IF UA
	cd RB/Rice;		make distclean    #^CMP IF RB
	cd RB/RiceV5;		make distclean    #^CMP IF RB
	cd RB/Empty;		make distclean    #^CMP IF RB
	cd SP/FLAMPA;           make distclean    #^CMP IF SP
	cd SP/Kota;             make distclean    #^CMP IF SP
	cd SP/Empty;            make distclean    #^CMP IF SP
	cd SC/BATSRUS;          make distclean    #^CMP IF SC
	cd SC/Empty;            make distclean    #^CMP IF SC
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
	@(if([ -f IH/BATSRUS/src/Makefile ]); then \
	echo 'Remove IH/BATSRUS/src'; \
	cd IH/BATSRUS; make veryclean; \
	fi)
	@(if([ -f SC/BATSRUS/src/Makefile ]); then \
	echo 'Remove SC/BATSRUS/src'; \
	cd SC/BATSRUS; make veryclean; \
	fi)
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
	tar -rf tmp.tar  Configure.options
	tar -rf tmp.tar  Configure.pl		#^CMP IF CONFIGURE
	tar -rf tmp.tar  doc			#^CMP IF DOC
	tar -rf tmp.tar  Param
	tar -rf tmp.tar  Scripts
	tar -rf tmp.tar  SetSWMF.pl
	tar -rf tmp.tar  share
	tar -rf tmp.tar  util
	tar -rf tmp.tar  CON
	tar -rf tmp.tar  GM			#^CMP IF GM
	tar -rf tmp.tar  IE			#^CMP IF IE
	tar -rf tmp.tar  IH			#^CMP IF IH
	tar -rf tmp.tar  IM			#^CMP IF IM
	tar -rf tmp.tar  UA			#^CMP IF UA
	tar -rf tmp.tar  RB                     #^CMP IF RB
	tar -rf tmp.tar  SP                     #^CMP IF SP
	tar -rf tmp.tar  SC                     #^CMP IF SC
	@echo ' '
	gzip tmp.tar
	mv tmp.tar.gz SWMF_v${VERSION}_`date +%Y%b%d_%H%M.tgz`
	@echo ' '
	@ls -l SWMF_v*.tgz

#			
#	Create run directories
#
rundir: ENV_CHECK
	mkdir run
	mkdir run/STDOUT
	cp Param/LAYOUT.DEFAULT run/LAYOUT.in
	cp Param/PARAM.DEFAULT run/PARAM.in
	cp Param/IMF.DEFAULT run/IMF.dat
	cp Scripts/Restart.pl run/Restart.pl
	touch run/core
	chmod 444 run/core
	cd run; ln -s  ${DIR}/bin/SWMF.exe . ; ln -s  ${DIR}/Param .
	cd ${GMDIR}; make rundir PLOT=${PLOT}    #^CMP IF GM
	cd ${IEDIR}; make rundir                 #^CMP IF IE
	cd ${IHDIR}; make rundir PLOT=${PLOT}    #^CMP IF IH
	cd ${IMDIR}; make rundir                 #^CMP IF IM
	cd ${UADIR}; make rundir                 #^CMP IF UA
	cd ${RBDIR}; make rundir                 #^CMP IF RB
	cd ${SPDIR}; make rundir                 #^CMP IF SP
	cd ${SCDIR}; make rundir                 #^CMP IF SC
	@touch CON/Scripts/${OS}/TMP_${MACHINE}
	cp CON/Scripts/${OS}/*${MACHINE}* run/
	@rm -rf run/TMP_${MACHINE} CON/Scripts/${OS}/TMP_${MACHINE}
	@echo
	@echo Creation of run directory succeeded
	@echo

#
#	Run the default code on NP processors
#
NP=8

mpirun: ENV_CHECK SWMF
	cd run; mpirun -np ${NP} ./SWMF.exe

mprun: ENV_CHECK SWMF
	cd run; mprun -np ${NP} ./SWMF.exe

nompirun: ENV_CHECK SWMF
	cd run; ./SWMF.exe

ETAGS = etags

tags:	ENV_CHECK
	-$(ETAGS) ./*/*/*/*.[fF]90 ./*/*/*/*.[fF] ./*/*/*/*.for

#^CMP IF IH BEGIN
#
#	Copy and rename GM/BATSRUS/src into IH/BATSRUS/src
#

IH/BATSRUS/src/Makefile:
	mkdir -p IH/BATSRUS/src
	cd GM/BATSRUS/src; cp *.f90 *.f Makefile* ../../../IH/BATSRUS/src
	cd IH/BATSRUS/src; rm -f main.f90 stand_alone*.f90
	cp GM/BATSRUS/srcInterface/ModGridDescriptor.f90 IH/BATSRUS/src
	cp GM/BATSRUS/srcInterface/update_lagrangian_grid.f90 IH/BATSRUS/src
	cp GM/BATSRUS/srcInterface/ModBuffer.f90 IH/BATSRUS/src
	cp IH/BATSRUS_share/src/IH_*.f90 IH/BATSRUS/src
	cd GM/BATSRUS; \
		cp Makefile.def PARAM.XML PARAM.pl GridSize.pl \
		../../IH/BATSRUS/
	echo '*' > IH/BATSRUS/src/.cvsignore

IHBATSRUS: IH/BATSRUS/src/Makefile \
		${SCRIPTDIR}/Methods.pl ${SCRIPTDIR}/Rename.pl
	cd IH/BATSRUS/src; \
		${SCRIPTDIR}/Methods.pl IH; \
		${SCRIPTDIR}/Rename.pl -r *.f90 *.f; \
		perl -i -pe 's[^(\s*common\s*/)(\w+/)][$$1IH_$$2]i' *.f90 *.f;\
		mv ModGridDescriptor.f90 ModBuffer.f90 \
		update_lagrangian_grid.f90 IH_*.f90 \
			../srcInterface
	cd IH/BATSRUS/srcInterface; \
		perl -i -pe 's?BATSRUS?IH_BATSRUS?' IH_*.f90; \
		touch Makefile.DEPEND

#^CMP END IH
#^CMP IF SC BEGIN
SC/BATSRUS/src/Makefile:
	cd GM/BATSRUS; \
		cp -f Makefile.conf ../../SC/BATSRUS; \
		make COMP=SC DREL=TMP relax_src
	cd GM/BATSRUS/TMP; \
		mv Makefile.def GridSize.pl PARAM.XML src ../../../SC/BATSRUS;\
		mv srcInterface/*.f90 ../../../SC/BATSRUS/src
	rm -rf GM/BATSRUS/TMP
	cp -f IH/BATSRUS_share/src/IH_wrapper.f90 \
		      SC/BATSRUS/src/SC_wrapper.f90
	cp -f IH/BATSRUS_share/src/IH_get_for_sp.f90 \
		      SC/BATSRUS/src/SC_get_for_sp.f90
	cd SC/BATSRUS/src/; perl -i -pe \
	's/IH/SC/g;s/BATSRUS/SC_BATSRUS/;s/Inner/Solar/;s/Heliosphere/Corona/'\
		SC_wrapper.f90 SC_get_for_sp.f90
	cd SC/BATSRUS/src; rm -f main.f90 stand_alone*.f90

SCBATSRUS: SC/BATSRUS/src/Makefile \
		${SCRIPTDIR}/Methods.pl ${SCRIPTDIR}/Rename.pl
	cd SC/BATSRUS/src; \
		${SCRIPTDIR}/Methods.pl SC; \
		${SCRIPTDIR}/Rename.pl -r *.f90 *.f; \
		perl -i -pe 's[^(\s*common\s*/)(\w+/)][$$1SC_$$2]i' *.f90 *.f;\
		mv ModGridDescriptor.f90 ModBuffer.f90 \
		update_lagrangian_grid.f90 SC_*.f90 \
			../srcInterface
	touch SC/BATSRUS/srcInterface/Makefile.DEPEND

#^CMP END SC

# keep this line
