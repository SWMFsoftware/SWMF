#^CMP COPYRIGHT UM

SHELL =/bin/sh

# CON needs modules from other subdirectories
SEARCH_EXTRA = -I${INTERFACEDIR} -I${LIBRARYDIR} -I${COUPLERDIR}

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND
include Makefile.RULES

#######################
# Shared libraries by #
# CON and componets:  #
#######################

LIBSHARE   = ${LIBDIR}/libSHARE.a

#######################
# Utility libraries:  #
#######################

LIBTIMING      = ${LIBDIR}/libTIMING.a
LIBINDICES     = ${LIBDIR}/libINDICES.a
LIBMAGNETOGRAM = ${LIBDIR}/libMAGNETOGRAM.a
LIBEMPIRICALEE = ${LIBDIR}/libEMPIRICALEE.a
LIBEMPIRICALGM = ${LIBDIR}/libEMPIRICALGM.a
LIBEMPIRICALIE = ${LIBDIR}/libEMPIRICALIE.a
LIBEMPIRICALUA = ${LIBDIR}/libEMPIRICALUA.a

#######################
# CON infrastructure: #
#######################

LIBLIBRARY = ${LIBDIR}/libLIBRARY.a

#######################
# Coupling toolkit:   #
#######################

LIBCOUPLER = ${LIBDIR}/libCOUPLER.a

#######################
# Components:         #
#######################
# Use alphabetic order for components!

LIBGM = ${LIBDIR}/libGM.a	#^CMP IF GM

LIBIE = ${LIBDIR}/libIE.a	#^CMP IF IE

LIBIH = ${LIBDIR}/libIH.a	#^CMP IF IH

LIBIM = ${LIBDIR}/libIM.a	#^CMP IF IM

LIBLA = ${LIBDIR}/libLA.a	#^CMP IF LA

LIBOH = ${LIBDIR}/libOH.a       #^CMP IF OH

LIBPS = ${LIBDIR}/libPS.a	#^CMP IF PS

LIBPW = ${LIBDIR}/libPW.a       #^CMP IF PW

LIBRB = ${LIBDIR}/libRB.a	#^CMP IF RB

LIBLC = ${LIBDIR}/libLC.a       #^CMP IF LC

LIBSC = ${LIBDIR}/libSC.a       #^CMP IF SC

LIBSP = ${LIBDIR}/libSP.a       #^CMP IF SP

LIBUA = ${LIBDIR}/libUA.a	#^CMP IF UA

#######################
# Interface part of   #
# super-structure:    #
#######################

LIBINTERFACE = ${LIBDIR}/libINTERFACE.a

#######################
# Control part of     #
# superstructure:     #
#######################

LIBCONTROL = ${LIBDIR}/libCONTROL.a

MY_LIB = ${LIBCONTROL}

# Note: CON_methods is actually infrastructure, but it consists of
#       external subroutines and functions only so it can be compiled here

LIBOBJECTS = \
	CON_main.o\
	CON_variables.o\
	CON_io.o\
	CON_session.o\
	CON_methods.o\
	swmf_interface.o

#
#	Main library
#

SWMFLIB = ${LIBDIR}/libSWMF.a

# SWMFLIB is a single library that contains all SWMF object files
SWMFLIB: TMP_DIRS
	rm -f ${SWMFLIB}
	${AR} ${SWMFLIB} Tmp_*/*.o
	@rm -rf Tmp_*
	@echo
	@echo SWMF library ${SWMFLIB} has been brought up to date.
	@echo

#
#	Main executable
#

SWMFEXE = ${BINDIR}/SWMF.exe

# The main program unit
EXEOBJECT = swmf.o

# Linking SWMF.exe from the libSWMF.a library fails on Mac OSX.
# As a work around we load the object files directly from Tmp_*

SWMFEXE: TMP_DIRS
	@make ${EXEOBJECT}
	${LINK.f90} -o ${SWMFEXE} ${EXEOBJECT} Tmp_*/*.o ${LBLAS} ${Lflag1}
	@rm -rf Tmp_*
	@echo
	@echo Program ${SWMFEXE} has been brought up to date.
	@echo

##############CRASH EXECUTABLE#############################
CRASHEXE = ${BINDIR}/CRASH.exe

# Linking CRASH.exe from the libSWMF.a library fails on Mac OSX.
# As a work around we load the object files directly from Tmp_*

CRASHEXE: TMP_CRASH
	@make ${EXEOBJECT}
	${LINK.f90} -o ${CRASHEXE} ${EXEOBJECT} Tmp_*/*.o ${LBLAS} ${Lflag1}
	@rm -rf Tmp_*
	@echo
	@echo Program ${CRASHEXE} has been brought up to date.
	@echo

###########################################################

# The order of compilation:
# - SHAREDIR  :     shared SWMF library                --> libSHARE.a
# - TIMINGDIR :     timing/profiling library           --> libTIMING.a
# - INDICESDIR:     space physics INDICES library      --> libINDICES.a
# - MAGNETOGRAMDIR: magnetogram reading library        --> LibMAGNETOGRAM.a
# - EMPIRICALDIR:   empirical models library           --> libEMPIRICAL.a
# - LIBRARYDIR:     CON "infrastructure"               --> libLIBRARY.a
# - COUPLERDIR:     coupling toolkit                   --> libCOUPLER.a
# - GMDIR .. UADIR: physics components                 --> libGM.a ... libUA.a
# - INTERFACEDIR:   interface part of "superstructure" --> libINTERFACE.a
# - this directory: control   part of "superstructure" --> libCONTROL.a
# - make SWMFEXE or SWMFLIB

# Create all libraries and then libCONTROL.a (MY_LIB)
# Use alphabetic order for components!
LIB:
	@cd ${SHAREDIR};             make LIB
	@cd ${TIMINGDIR};            make LIB
	@cd ${DATAREADINDICESDIR};   make LIB
	@cd ${MAGNETOGRAMDIR};       make LIB
	@cd ${EMPIRICALEEDIR};       make LIB
	@cd ${EMPIRICALGMDIR};       make LIB
	@cd ${EMPIRICALIEDIR};       make LIB
	@cd ${EMPIRICALUADIR};       make LIB
	@cd ${LIBRARYDIR};           make LIB
	@cd ${COUPLERDIR};           make LIB
	@cd ${GMDIR};                make LIB	#^CMP IF GM
	@cd ${IEDIR};                make LIB	#^CMP IF IE
	@cd ${IHDIR};                make LIB	#^CMP IF IH
	@cd ${IMDIR};                make LIB	#^CMP IF IM
	@cd ${LADIR};                make LIB	#^CMP IF LA
	@cd ${OHDIR};                make LIB   #^CMP IF OH
	@cd ${PSDIR};                make LIB	#^CMP IF PS
	@cd ${PWDIR};                make LIB	#^CMP IF PW
	@cd ${RBDIR};                make LIB	#^CMP IF RB
	@cd ${LCDIR};		     make LIB   #^CMP IF LC
	@cd ${SCDIR};                make LIB   #^CMP IF SC
	@cd ${SPDIR};                make LIB   #^CMP IF SP
	@cd ${UADIR};                make LIB	#^CMP IF UA
	@cd ${INTERFACEDIR};         make LIB
	@make DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

##########################CRASH_LIB##############
CRASH_LIB:
	@cd ${SHAREDIR};             make LIB
	@cd ${TIMINGDIR};            make LIB
	@cd ${DATAREADINDICESDIR};   make LIB
	@cd ${MAGNETOGRAMDIR};       make LIB
	@cd ${EMPIRICALEEDIR};       make LIB
	@cd ${EMPIRICALGMDIR};       make LIB
	@cd ${EMPIRICALIEDIR};       make LIB
	@cd ${EMPIRICALUADIR};       make LIB
	@cd ${CRASHDIR};             make LIB
	@cd ${LIBRARYDIR};           make LIB
	@cd ${COUPLERDIR};           make LIB
	@cd ${GMDIR};                make LIB	#^CMP IF GM
	@cd ${IEDIR};                make LIB	#^CMP IF IE
	@cd ${IHDIR};                make LIB	#^CMP IF IH
	@cd ${IMDIR};                make LIB	#^CMP IF IM
	@cd ${LADIR};                make LIB	#^CMP IF LA
	@cd ${OHDIR};                make LIB   #^CMP IF OH
	@cd ${PSDIR};                make LIB	#^CMP IF PS
	@cd ${PWDIR};                make LIB	#^CMP IF PW
	@cd ${RBDIR};                make LIB	#^CMP IF RB
	@cd ${LCDIR};		     make LIB   #^CMP IF LC
	@cd ${SCDIR};                make LIB   #^CMP IF SC
	@cd ${SPDIR};                make LIB   #^CMP IF SP
	@cd ${UADIR};                make LIB	#^CMP IF UA
	@cd ${INTERFACEDIR};         make LIB
	@make DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

###############END CRASH_LIB#####

${MY_LIB}: ${LIBOBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${LIBOBJECTS}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${LIBOBJECTS} ${EXEOBJECT}

# For now the archive extract is defined here. Should be in Makefile.conf
AR_X = ar -x

# Extract object files from libraries and collect them into SWMFLIB
# Use alphabetic order for components!
TMP_DIRS: LIB
	rm -rf Tmp_*
	mkdir Tmp_SHARE;       cd Tmp_SHARE;       ${AR_X} ${LIBSHARE}
	mkdir Tmp_TIMING;      cd Tmp_TIMING;      ${AR_X} ${LIBTIMING} 
	mkdir Tmp_INDICES;     cd Tmp_INDICES;     ${AR_X} ${LIBINDICES} 
	mkdir Tmp_MAGNETOGRAM; cd Tmp_MAGNETOGRAM; ${AR_X} ${LIBMAGNETOGRAM} 
	mkdir Tmp_EMPIRICALEE; cd Tmp_EMPIRICALEE; ${AR_X} ${LIBEMPIRICALEE} 
	mkdir Tmp_EMPIRICALGM; cd Tmp_EMPIRICALGM; ${AR_X} ${LIBEMPIRICALGM} 
	mkdir Tmp_EMPIRICALIE; cd Tmp_EMPIRICALIE; ${AR_X} ${LIBEMPIRICALIE} 
	mkdir Tmp_EMPIRICALUA; cd Tmp_EMPIRICALUA; ${AR_X} ${LIBEMPIRICALUA} 
	mkdir Tmp_LIBRARY;     cd Tmp_LIBRARY;   ${AR_X} ${LIBLIBRARY}
	mkdir Tmp_COUPLER;     cd Tmp_COUPLER;   ${AR_X} ${LIBCOUPLER}
	mkdir Tmp_INTERFACE;   cd Tmp_INTERFACE; ${AR_X} ${LIBINTERFACE}
	mkdir Tmp_CONTROL;     cd Tmp_CONTROL;   ${AR_X} ${LIBCONTROL}
	mkdir Tmp_LIBGM;       cd Tmp_LIBGM;     ${AR_X} ${LIBGM} #^CMP IF GM
	mkdir Tmp_LIBIE;       cd Tmp_LIBIE;     ${AR_X} ${LIBIE} #^CMP IF IE
	mkdir Tmp_LIBIH;       cd Tmp_LIBIH;     ${AR_X} ${LIBIH} #^CMP IF IH
	mkdir Tmp_LIBIM;       cd Tmp_LIBIM;     ${AR_X} ${LIBIM} #^CMP IF IM
	mkdir Tmp_LIBLA;       cd Tmp_LIBLA;     ${AR_X} ${LIBLA} #^CMP IF LA
	mkdir Tmp_LIBOH;       cd Tmp_LIBOH;     ${AR_X} ${LIBOH} #^CMP IF OH
	mkdir Tmp_LIBPS;       cd Tmp_LIBPS;     ${AR_X} ${LIBPS} #^CMP IF PS
	mkdir Tmp_LIBPW;       cd Tmp_LIBPW;     ${AR_X} ${LIBPW} #^CMP IF PW
	mkdir Tmp_LIBRB;       cd Tmp_LIBRB;     ${AR_X} ${LIBRB} #^CMP IF RB
	mkdir Tmp_LIBLC;       cd Tmp_LIBLC;     ${AR_X} ${LIBLC} #^CMP IF LC
	mkdir Tmp_LIBSC;       cd Tmp_LIBSC;     ${AR_X} ${LIBSC} #^CMP IF SC
	mkdir Tmp_LIBSP;       cd Tmp_LIBSP;     ${AR_X} ${LIBSP} #^CMP IF SP
	mkdir Tmp_LIBUA;       cd Tmp_LIBUA;     ${AR_X} ${LIBUA} #^CMP IF UA

#
#	cleaning
#

distclean: clean
	rm -rf Makefile.DEPEND Tmp_*
######################CRASH######################
LIBCRASH       = ${LIBDIR}/libCRASH.a

TMP_CRASH: CRASH_LIB
	rm -rf Tmp_*
	mkdir Tmp_SHARE;       cd Tmp_SHARE;       ${AR_X} ${LIBSHARE}
	mkdir Tmp_TIMING;      cd Tmp_TIMING;      ${AR_X} ${LIBTIMING} 
	mkdir Tmp_INDICES;     cd Tmp_INDICES;     ${AR_X} ${LIBINDICES} 
	mkdir Tmp_MAGNETOGRAM; cd Tmp_MAGNETOGRAM; ${AR_X} ${LIBMAGNETOGRAM} 
	mkdir Tmp_CRASH;       cd Tmp_CRASH;       ${AR_X} ${LIBCRASH} 
	mkdir Tmp_EMPIRICALEE; cd Tmp_EMPIRICALEE; ${AR_X} ${LIBEMPIRICALEE} 
	mkdir Tmp_EMPIRICALGM; cd Tmp_EMPIRICALGM; ${AR_X} ${LIBEMPIRICALGM} 
	mkdir Tmp_EMPIRICALIE; cd Tmp_EMPIRICALIE; ${AR_X} ${LIBEMPIRICALIE} 
	mkdir Tmp_EMPIRICALUA; cd Tmp_EMPIRICALUA; ${AR_X} ${LIBEMPIRICALUA} 
	mkdir Tmp_LIBRARY;     cd Tmp_LIBRARY;   ${AR_X} ${LIBLIBRARY}
	mkdir Tmp_COUPLER;     cd Tmp_COUPLER;   ${AR_X} ${LIBCOUPLER}
	mkdir Tmp_INTERFACE;   cd Tmp_INTERFACE; ${AR_X} ${LIBINTERFACE}
	mkdir Tmp_CONTROL;     cd Tmp_CONTROL;   ${AR_X} ${LIBCONTROL}
	mkdir Tmp_LIBGM;       cd Tmp_LIBGM;     ${AR_X} ${LIBGM} #^CMP IF GM
	mkdir Tmp_LIBIE;       cd Tmp_LIBIE;     ${AR_X} ${LIBIE} #^CMP IF IE
	mkdir Tmp_LIBIH;       cd Tmp_LIBIH;     ${AR_X} ${LIBIH} #^CMP IF IH
	mkdir Tmp_LIBIM;       cd Tmp_LIBIM;     ${AR_X} ${LIBIM} #^CMP IF IM
	mkdir Tmp_LIBLA;       cd Tmp_LIBLA;     ${AR_X} ${LIBLA} #^CMP IF LA
	mkdir Tmp_LIBOH;       cd Tmp_LIBOH;     ${AR_X} ${LIBOH} #^CMP IF OH
	mkdir Tmp_LIBPS;       cd Tmp_LIBPS;     ${AR_X} ${LIBPS} #^CMP IF PS
	mkdir Tmp_LIBPW;       cd Tmp_LIBPW;     ${AR_X} ${LIBPW} #^CMP IF PW
	mkdir Tmp_LIBRB;       cd Tmp_LIBRB;     ${AR_X} ${LIBRB} #^CMP IF RB
	mkdir Tmp_LIBLC;       cd Tmp_LIBLC;     ${AR_X} ${LIBLC} #^CMP IF LC
	mkdir Tmp_LIBSC;       cd Tmp_LIBSC;     ${AR_X} ${LIBSC} #^CMP IF SC
	mkdir Tmp_LIBSP;       cd Tmp_LIBSP;     ${AR_X} ${LIBSP} #^CMP IF SP
	mkdir Tmp_LIBUA;       cd Tmp_LIBUA;     ${AR_X} ${LIBUA} #^CMP IF UA
