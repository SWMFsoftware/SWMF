#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

# CON needs modules from other subdirectories
SEARCHDIR = -I${SHAREDIR} -I${INTERFACEDIR} -I${LIBRARYDIR} -I${COUPLERDIR}

include ../../../Makefile.def
include ../../../Makefile.conf
-include Makefile.DEPEND
-include Makefile.RULES

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
LIBDEMT        = ${LIBDIR}/libDEMT.a 
LIBEMPIRICALEE = ${LIBDIR}/libEMPIRICALEE.a
LIBEMPIRICALGM = ${LIBDIR}/libEMPIRICALGM.a
LIBEMPIRICALIE = ${LIBDIR}/libEMPIRICALIE.a
LIBEMPIRICALSC = ${LIBDIR}/libEMPIRICALSC.a
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

LIBEE = ${LIBDIR}/libEE.a	#^CMP IF EE

LIBGM = ${LIBDIR}/libGM.a	#^CMP IF GM

LIBIE = ${LIBDIR}/libIE.a	#^CMP IF IE

LIBIH = ${LIBDIR}/libIH.a	#^CMP IF IH

LIBIM = ${LIBDIR}/libIM.a	#^CMP IF IM

LIBOH = ${LIBDIR}/libOH.a       #^CMP IF OH

LIBPC = ${LIBDIR}/libPC.a       #^CMP IF PC

LIBPS = ${LIBDIR}/libPS.a	#^CMP IF PS

LIBPT = ${LIBDIR}/libPT.a	#^CMP IF PT

LIBPW = ${LIBDIR}/libPW.a       #^CMP IF PW

LIBRB = ${LIBDIR}/libRB.a	#^CMP IF RB

LIBSC = ${LIBDIR}/libSC.a       #^CMP IF SC

LIBSP = ${LIBDIR}/libSP.a       #^CMP IF SP

LIBUA = ${LIBDIR}/libUA.a	#^CMP IF UA

LIBCZ = ${LIBDIR}/libCZ.a    	#^CMP IF CZ

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

LIBOBJECTS = \
	CON_main.o\
	CON_variables.o\
	CON_io.o\
	CON_session.o\
	swmf_interface.o

#
#	Main library
#

SWMFLIB = ${LIBDIR}/libSWMF.a

# SWMFLIB is a single library that contains all SWMF object files
# Note that TMP_DIRS should not be done with parallel make !
SWMFLIB:
	make TMP_DIRS
	rm -f ${SWMFLIB}
	${AR} ${SWMFLIB} Tmp_*/*.o
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
# Note that TMP_DIRS should not be done with parallel make !
SWMFEXE:
	make TMP_DIRS
	@$(MAKE) ${EXEOBJECT}
	${LINK.f90} -o ${SWMFEXE} ${EXEOBJECT} Tmp_*/*.o ${LBLAS} ${LflagMpi} \
		${HYPRELIB} ${SPICELIB} ${FORTRANLIB}
	@echo
	@echo Program ${SWMFEXE} has been brought up to date.
	@echo

###########################################################

# The order of compilation:
# - SHAREDIR  :     shared SWMF library                --> libSHARE.a
# - TIMINGDIR :     timing/profiling library           --> libTIMING.a
# - INDICESDIR:     space physics INDICES library      --> libINDICES.a
# - MAGNETOGRAMDIR: magnetogram reading library        --> LibMAGNETOGRAM.a
# - DEMTDIR:        differential emmission measure read--> libDEMT.a
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
	@cd ${SHAREDIR};             $(MAKE) LIB
	@cd ${TIMINGDIR};            $(MAKE) LIB
	@cd ${DATAREADINDICESDIR};   $(MAKE) LIB
	@cd ${MAGNETOGRAMDIR};       $(MAKE) LIB
	@cd ${DEMTDIR}; 	     $(MAKE) LIB
	@cd ${EMPIRICALEEDIR};       $(MAKE) LIB
	@cd ${EMPIRICALGMDIR};       $(MAKE) LIB
	@cd ${EMPIRICALIEDIR};       $(MAKE) LIB
	@cd ${EMPIRICALSCDIR};       $(MAKE) LIB
	@cd ${EMPIRICALUADIR};       $(MAKE) LIB
	@cd ${LIBRARYDIR};           $(MAKE) LIB
	@cd ${COUPLERDIR};           $(MAKE) LIB
	@cd ${EEDIR};                $(MAKE) LIB	#^CMP IF EE
	@cd ${GMDIR};                $(MAKE) LIB	#^CMP IF GM
	@cd ${IEDIR};                $(MAKE) LIB	#^CMP IF IE
	@cd ${IHDIR};                $(MAKE) LIB	#^CMP IF IH
	@cd ${IMDIR};                $(MAKE) LIB	#^CMP IF IM
	@cd ${PCDIR};                $(MAKE) LIB	#^CMP IF PC
	@cd ${OHDIR};                $(MAKE) LIB	#^CMP IF OH
	@cd ${PSDIR};                $(MAKE) LIB	#^CMP IF PS
	@cd ${PTDIR};                $(MAKE) LIB	#^CMP IF PT
	@cd ${PWDIR};                $(MAKE) LIB	#^CMP IF PW
	@cd ${RBDIR};                $(MAKE) LIB	#^CMP IF RB
	@cd ${SCDIR};                $(MAKE) LIB	#^CMP IF SC
	@cd ${SPDIR};                $(MAKE) LIB	#^CMP IF SP
	@cd ${UADIR};                $(MAKE) LIB	#^CMP IF UA
	@cd ${CZDIR}; 		     $(MAKE) LIB        #^CMP IF CZ
	@cd ${INTERFACEDIR};         $(MAKE) LIB
	@make DEPEND
	$(MAKE) ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${LIBOBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${LIBOBJECTS}

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${LIBOBJECTS} ${EXEOBJECT}

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
	mkdir Tmp_DEMT;        cd Tmp_DEMT;        ${AR_X} ${LIBDEMT} 
	mkdir Tmp_EMPIRICALEE; cd Tmp_EMPIRICALEE; ${AR_X} ${LIBEMPIRICALEE} 
	mkdir Tmp_EMPIRICALGM; cd Tmp_EMPIRICALGM; ${AR_X} ${LIBEMPIRICALGM} 
	mkdir Tmp_EMPIRICALIE; cd Tmp_EMPIRICALIE; ${AR_X} ${LIBEMPIRICALIE}
	mkdir Tmp_EMPIRICALSC; cd Tmp_EMPIRICALSC; ${AR_X} ${LIBEMPIRICALSC}
	mkdir Tmp_EMPIRICALUA; cd Tmp_EMPIRICALUA; ${AR_X} ${LIBEMPIRICALUA} 
	mkdir Tmp_LIBRARY;     cd Tmp_LIBRARY;   ${AR_X} ${LIBLIBRARY}
	mkdir Tmp_COUPLER;     cd Tmp_COUPLER;   ${AR_X} ${LIBCOUPLER}
	mkdir Tmp_INTERFACE;   cd Tmp_INTERFACE; ${AR_X} ${LIBINTERFACE}
	mkdir Tmp_CONTROL;     cd Tmp_CONTROL;   ${AR_X} ${LIBCONTROL}
	mkdir Tmp_LIBEE;       cd Tmp_LIBEE;     ${AR_X} ${LIBEE} #^CMP IF EE
	mkdir Tmp_LIBGM;       cd Tmp_LIBGM;     ${AR_X} ${LIBGM} #^CMP IF GM
	mkdir Tmp_LIBIE;       cd Tmp_LIBIE;     ${AR_X} ${LIBIE} #^CMP IF IE
	mkdir Tmp_LIBIH;       cd Tmp_LIBIH;     ${AR_X} ${LIBIH} #^CMP IF IH
	mkdir Tmp_LIBIM;       cd Tmp_LIBIM;     ${AR_X} ${LIBIM} #^CMP IF IM
	mkdir Tmp_LIBOH;       cd Tmp_LIBOH;     ${AR_X} ${LIBOH} #^CMP IF OH
	mkdir Tmp_LIBPC;       cd Tmp_LIBPC;     ${AR_X} ${LIBPC} #^CMP IF PC
	mkdir Tmp_LIBPS;       cd Tmp_LIBPS;     ${AR_X} ${LIBPS} #^CMP IF PS
	mkdir Tmp_LIBPT;       cd Tmp_LIBPT;     ${AR_X} ${LIBPT} #^CMP IF PT
	mkdir Tmp_LIBPW;       cd Tmp_LIBPW;     ${AR_X} ${LIBPW} #^CMP IF PW
	mkdir Tmp_LIBRB;       cd Tmp_LIBRB;     ${AR_X} ${LIBRB} #^CMP IF RB
	mkdir Tmp_LIBSC;       cd Tmp_LIBSC;     ${AR_X} ${LIBSC} #^CMP IF SC
	mkdir Tmp_LIBSP;       cd Tmp_LIBSP;     ${AR_X} ${LIBSP} #^CMP IF SP
	mkdir Tmp_LIBUA;       cd Tmp_LIBUA;     ${AR_X} ${LIBUA} #^CMP IF UA
	mkdir Tmp_LIBCZ;       cd Tmp_LIBCZ;     ${AR_X} ${LIBCZ} #^CMP IF CZ
#
#	cleaning
#

clean: cleanfiles
	rm -rf Tmp_*

distclean: clean
