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
LIBCRASH       = ${LIBDIR}/LibCRASH.a
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

LIBGM = ${LIBDIR}/libGM.a	#^CMP IF GM

LIBIE = ${LIBDIR}/libIE.a	#^CMP IF IE

LIBIH = ${LIBDIR}/libIH.a	#^CMP IF IH

LIBOH = ${LIBDIR}/libOH.a       #^CMP IF OH

LIBIM = ${LIBDIR}/libIM.a	#^CMP IF IM

LIBLA = ${LIBDIR}/libLA.a	#^CMP IF LA

LIBPS = ${LIBDIR}/libPS.a	#^CMP IF PS

LIBPW = ${LIBDIR}/libPW.a       #^CMP IF PW

LIBRB = ${LIBDIR}/libRB.a	#^CMP IF RB

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

#
#	Main executable
#

SWMFEXE = ${BINDIR}/SWMF.exe

SWMFEXE: SWMFLIB
	@make ${SWMFEXE}
	@echo ' '
	@echo Program ${SWMFEXE} has been brought up to date.
	@echo ' '

EXEOBJECT = swmf.o

# The order of compilation:
# - SHAREDIR  : contains shared library               --> libSHARE.a
# - TIMINGDIR : contains TIMING software library      --> libTIMING.a
# - INDICESDIR: contains INDICES software library     --> libINDICES.a
# - MAGNETOGRAMDIR: contains magnetogram software library    --> LibMAGNETOGRAM.a
# - EMPIRICALDIR: contains empirical models software library --> libEMPIRICAL.a
# - LIBRARYDIR: the CON "infrastructure"              --> libLIBRARY.a
# - COUPLERDIR: the coupling toolkit                  --> libCOUPLER.a
# - GMDIR .. UADIR: the physics components            --> libGM.a ... libUA.a
# - INTERFACEDIR: interface part of "superstructure"  --> libINTERFACE.a
# - this directory: control part of "superstructure"  --> libCONTROL.a
# - make SWMFEXE or SWMFLIB

LIB:
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
	@cd ${OHDIR};                make LIB   #^CMP IF OH
	@cd ${IMDIR};                make LIB	#^CMP IF IM
	@cd ${LADIR};                make LIB	#^CMP IF LA
	@cd ${PSDIR};                make LIB	#^CMP IF PS
	@cd ${PWDIR};                make LIB	#^CMP IF PW
	@cd ${RBDIR};                make LIB	#^CMP IF RB
	@cd ${SCDIR};                make LIB   #^CMP IF SC
	@cd ${SPDIR};                make LIB   #^CMP IF SP
	@cd ${UADIR};                make LIB	#^CMP IF UA
	@cd ${INTERFACEDIR};         make LIB
	@make DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${LIBOBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${LIBOBJECTS}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${LIBOBJECTS} ${EXEOBJECT}

# SWMFLIB is a single library that contains all object files
SWMFLIB: LIB
	make ${SWMFLIB}
	@echo
	@echo SWMF library ${SWMFLIB} has been brought up to date.
	@echo

# For now the archive extract is defined here. Should be in Makefile.conf
AR_X = ar -x

# Extract object files from libraries and collect them into SWMFLIB
${SWMFLIB}:	${LIBSHARE} \
		${LIBTIMING} \
		${LIBINDICES} \
		${LIBMAGNETOGRAM} \
		${LIBEMPIRICALEE} \
		${LIBEMPIRICALGM} \
		${LIBEMPIRICALIE} \
		${LIBEMPIRICALUA} \
		${LIBCRASH} \
		${LIBLIBRARY} \
		${LIBCOUPLER} \
		${LIBINTERFACE} \
		${LIBCONTROL} \
		${LIBGM} \
		${LIBIE} \
		${LIBIH} \
		${LIBOH} \
		${LIBIM} \
		${LIBLA} \
		${LIBPS} \
		${LIBPW} \
		${LIBRB} \
		${LIBSC} \
		${LIBSP} \
		${LIBUA}
	rm -rf Tmp_*
	mkdir Tmp_SHARE;       cd Tmp_SHARE;       ${AR_X} ${LIBSHARE}
	mkdir Tmp_TIMING;      cd Tmp_TIMING;      ${AR_X} ${LIBTIMING} 
	mkdir Tmp_INDICES;     cd Tmp_INDICES;     ${AR_X} ${LIBINDICES} 
	mkdir Tmp_MAGNETOGRAM; cd Tmp_MAGNETOGRAM; ${AR_X} ${LIBMAGNETOGRAM} 
	mkdir Tmp_EMPIRICALEE; cd Tmp_EMPIRICALEE; ${AR_X} ${LIBEMPIRICALEE} 
	mkdir Tmp_EMPIRICALGM; cd Tmp_EMPIRICALGM; ${AR_X} ${LIBEMPIRICALGM} 
	mkdir Tmp_EMPIRICALIE; cd Tmp_EMPIRICALIE; ${AR_X} ${LIBEMPIRICALIE} 
	mkdir Tmp_EMPIRICALUA; cd Tmp_EMPIRICALUA; ${AR_X} ${LIBEMPIRICALUA} 
	mkdir Tmp_CRASH;       cd Tmp_CRASH;       ${AR_X} ${LIBCRASH}
	mkdir Tmp_LIBRARY;     cd Tmp_LIBRARY;   ${AR_X} ${LIBLIBRARY}
	mkdir Tmp_COUPLER;     cd Tmp_COUPLER;   ${AR_X} ${LIBCOUPLER}
	mkdir Tmp_INTERFACE;   cd Tmp_INTERFACE; ${AR_X} ${LIBINTERFACE}
	mkdir Tmp_CONTROL;     cd Tmp_CONTROL;   ${AR_X} ${LIBCONTROL}
	mkdir Tmp_LIBGM;       cd Tmp_LIBGM;     ${AR_X} ${LIBGM} #^CMP IF GM
	mkdir Tmp_LIBIE;       cd Tmp_LIBIE;     ${AR_X} ${LIBIE} #^CMP IF IE
	mkdir Tmp_LIBIH;       cd Tmp_LIBIH;     ${AR_X} ${LIBIH} #^CMP IF IH
	mkdir Tmp_LIBOH;       cd Tmp_LIBOH;     ${AR_X} ${LIBOH} #^CMP IF OH
	mkdir Tmp_LIBIM;       cd Tmp_LIBIM;     ${AR_X} ${LIBIM} #^CMP IF IM
	mkdir Tmp_LIBLA;       cd Tmp_LIBLA;     ${AR_X} ${LIBLA} #^CMP IF LA
	mkdir Tmp_LIBPS;       cd Tmp_LIBPS;     ${AR_X} ${LIBPS} #^CMP IF PS
	mkdir Tmp_LIBPW;       cd Tmp_LIBPW;     ${AR_X} ${LIBPW} #^CMP IF PW
	mkdir Tmp_LIBRB;       cd Tmp_LIBRB;     ${AR_X} ${LIBRB} #^CMP IF RB
	mkdir Tmp_LIBSC;       cd Tmp_LIBSC;     ${AR_X} ${LIBSC} #^CMP IF SC
	mkdir Tmp_LIBSP;       cd Tmp_LIBSP;     ${AR_X} ${LIBSP} #^CMP IF SP
	mkdir Tmp_LIBUA;       cd Tmp_LIBUA;     ${AR_X} ${LIBUA} #^CMP IF UA
	rm -f ${SWMFLIB}
	@echo '${AR} ${SWMFLIB} Tmp_*/*.o'
	@${AR} ${SWMFLIB} Tmp_*/*.o
	@rm -rf Tmp_*

# Link the executable. The link procedure depends on the operating system.
${SWMFEXE}: ${EXEOBJECT} ${SWMFLIB} ${MAKEFILE_COMP_SELECT}
	@(if [ "${OS}" = "Darwin" ];         then make SWMF_for_Darwin;\
	else                                      make SWMF_general; fi)

# Link rule for all machines/compilers except for Darwin with NAG
SWMF_general:
	${LINK.f90} -o ${SWMFEXE} \
		${EXEOBJECT} -L${LIBDIR} -lSWMF ${LBLAS} ${Lflag1}

# Link rule for Darwin with the NAG compiler (and possibly other compilers)
# Note: for the xlf compiler the above general rule works too.
#
# This work around creates a temporary object file from all the F90
# libraries using the -all_load option of the linker.

SWMF_for_Darwin:
	ld -all_load -r -o LIBRARIES.o ${SWMFLIB}
	${LINK.f90} -o ${SWMFEXE} \
		${EXEOBJECT} LIBRARIES.o ${LBLAS} ${Lflag1}
	rm -f LIBRARIES.o

#
#	cleaning
#

distclean: clean
	rm -f Makefile.DEPEND

# keep this line
