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
	@echo '    PSPH    (PostSPH program creates spherical tec file from sph*.tec files)' #^CFG IF NOT SIMPLE
	@echo ' '
	@echo '    help      (makefile option list)'
	@echo '    clean     (rm -f *~ *.o *.kmo *.mod *.T *.lst core)'
	@echo '    distclean (make clean; rm -f *exe Makefile Makefile.DEPEND Makefile.OPTIONS)'
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
	@#    ^CFG IF NOT CARTESIAN BEGIN
	@echo "    spherical_src  (Make SPHERICAL directory with BATSRUS on spherical grid)"
	@echo "    spherical_conf (Make SPHERICAL directory and link it to BATSRUS_conf)"
	@echo "    covariant_src  (Make COVARIANT directory with BATSRUS on covariant grid)"
	@echo "    corelax_src    (Make CORELAX directory for the covariant version "
	@echo "                of the magnetogram-driven solar wind)"
	@#    ^CFG END CARTESIAN
	@#    ^CFG IF CARTESIAN BEGIN
	@echo "    cartesian_src  (removes source code for covariant grid)"
	@echo "    relax_src      (Make RELAX directory for the Cartesian "
	@echo "                    version of the magnetogram-driven solar wind)"
	@#    ^CFG END CARTESIAN
	@echo ' '
	@#^CFG END CONFIGURE


install: MAKEFILE_DEF
	@make install_cont;

MAKEFILE_DEF:
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		echo GMDIR=`pwd`                        >  Makefile.def; \
		echo OS=`uname`                         >> Makefile.def; \
		echo STANDALONE=${STANDALONE}           >> Makefile.def; \
		cat src/Makefile.def                    >> Makefile.def; \
	fi);

install_cont: src/Makefile.RULES src/Makefile.OPTIONS src/ModSize.f90
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cp -f share/build/Makefile.${OS}${COMPILER} Makefile.conf; \
		cp -f src/stand_alone_${STANDALONE}.f90 src/stand_alone.f90;\
		cd share; make install;\
	else \
		echo include $(DIR)/Makefile.conf > Makefile.conf; \
	fi);
	touch src/Makefile.DEPEND srcInterface/Makefile.DEPEND
	cd src; make user_routines.f90 #^CFG IF USERFILES

src/Makefile.RULES:
	cp -f src/Makefile.RULES.${OS}${COMPILER} src/Makefile.RULES

src/Makefile.OPTIONS:
	cp -f src/Makefile.OPTIONS.ORIG src/Makefile.OPTIONS

src/ModSize.f90:
	cp -f src/ModSize_${OS}.f90 src/ModSize.f90

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

#^CFG IF NOT SIMPLE BEGIN
PSPH:
	cd srcPostProc; make PSPH
	@echo ' '
	@echo Program PostSPH has been brought up to date.
	@echo ' '
#^CFG END SIMPLE

checklink:
	@(if perl -e 'exit(-l "Scripts/Run/pp_IDL/PostIDL.exe")' ; then \
	   echo "Creating link to PostIDL.exe in Scripts/Run/pp_IDL/"; \
	   cd Scripts/Run/pp_IDL/; \
	   rm -rf PostIDL.exe; \
	   ln -s ${BINDIR}/PostIDL.exe .; \
	fi)
	@#^CFG IF NOT SIMPLE BEGIN
	@(if perl -e 'exit(-l "Scripts/Run/pp_TEC/PostSPH.exe")' ; then \
	   echo "Creating link to PostSPH.exe in Scripts/Run/pp_TEC/"; \
	   cd Scripts/Run/pp_TEC/; \
	   rm -rf PostSPH.exe; \
	   ln -s ${BINDIR}/PostSPH.exe .; \
	fi)
	@#^CFG END SIMPLE

# Ths PLOT=??? will copy both IDL and TEC scripts into the run directory
PLOT = ???

# The MACHINE variable holds the machine name for which scripts should
# be copied to the run directory when it is created.  This is used mostly
# when several different machines have the same operating system,
# but they require different batch queue scripts.
# If MACHINE is empty or not defined, all scripts for the current OS will
# be copied.
#
# The default is the short name of the current machine
MACHINE = `hostname | sed -e 's/\..*//'`

rundir: checklink
	mkdir -p ${RUNDIR}/GM
	cd ${RUNDIR}/GM; mkdir restartIN restartOUT IO2
	cd ${RUNDIR}/GM; ln -s ${GMDIR}/Scripts/Run/pp_${PLOT}/[Pp]* .
	cd ${RUNDIR}/GM; ln -s ${GMDIR}/Param .
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cp -f Param/PARAM.DEFAULT run/PARAM.in; \
		touch run/core; chmod 444 run/core; \
		touch Scripts/Run/${OS}/TMP_${MACHINE}; \
		cp Scripts/Run/${OS}/*${MACHINE}* run/; \
		rm -f run/TMP_${MACHINE} Scripts/Run/${OS}/TMP_${MACHINE}; \
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
	@touch src/Makefile.DEPEND src/Makefile.RULES src/Makefile.OPTIONS
	cd src; make clean
	@touch srcInterface/Makefile.DEPEND
	cd srcInterface; make clean

distclean:
	@touch src/Makefile.DEPEND src/Makefile.RULES src/Makefile.OPTIONS
	cd src; make distclean
	@touch srcInterface/Makefile.DEPEND
	cd srcInterface; make distclean
	@				#^CFG IF DOC BEGIN
	@					#^CFG IF NOT REMOVEDOCTEX BEGIN
	cd Doc/Tex; make clean ${CLEAN1} ${CLEAN2}
	@					#^CFG END REMOVEDOCTEX
	@				#^CFG END DOC
	@(if [ -d util  ]; then cd util;  make distclean; fi);
	@(if [ -d share ]; then cd share; make distclean; fi);
	rm -f Makefile.conf *~

include Makefile_CONFIGURE #^CFG IF CONFIGURE
