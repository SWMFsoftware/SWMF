#^CFG COPYRIGHT UM

SHELL =/bin/sh

# CON also needs modules from the Interface subdirectory
SEARCH_EXTRA = -I${INTERFACEDIR} -I${LIBRARYDIR} -I${COUPLERDIR}

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND

#######################
# Shared libraries by #
# CON and componets:  #
#######################

LIBSHARE   = ${LIBDIR}/libSHARE.a

#######################
# Software libraries: #
#######################

LIBTIMING = ${LIBDIR}/libTIMING.a

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

LIBIM = ${LIBDIR}/libIM.a	#^CMP IF IM

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

# Note: CON_methods is actually infrastructure, but it consists of
#       external subroutines and functions only so it can be compiled here

MODULES = \
	CON_variables.o\
	CON_io.o\
	CON_session.o

OBJECTS = \
	CON_main.o\
	CON_methods.o

#
#	Main executable
#

SWMFEXE = ${BINDIR}/SWMF.exe

# The order of compilation:
# - TIMINGDIR : contains TIMING software library      --> libTIMING.a
# - SHAREDIR  : contains shared library               --> libSHARE.a
# - LIBRARYDIR: the CON "infrastructure"    --> libLIBRARY.a
# - COUPLERDIR: the coupling toolkit        --> libCOUPLER.a
# - GMDIR .. UADIR: the physics components   --> libGM.a ... libUA.a
# - INTERFACEDIR: interface part of "superstructure"   --> libINTERFACE.a
# - this directory contains control part of "superstructure" --> SWMF.exe

SWMFEXE: DEPEND 
	@cd ${SHAREDIR};     make LIB
	@cd ${TIMINGDIR};    make LIB
	@cd ${LIBRARYDIR};   make LIB
	@cd ${COUPLERDIR};   make LIB
	@cd ${RBDIR};        make LIB	#^CMP IF RB
	@cd ${GMDIR};        make LIB	#^CMP IF GM
	@cd ${IEDIR};        make LIB	#^CMP IF IE
	@cd ${IHDIR};        make LIB	#^CMP IF IH
	@cd ${IMDIR};        make LIB	#^CMP IF IM
	@cd ${UADIR};        make LIB	#^CMP IF UA
	@cd ${INTERFACEDIR}; make LIB
	@make ${SWMFEXE}
	@echo ' '
	@echo Program ${SWMFEXE} has been brought up to date.
	@echo ' '

# CON uses modules from ${SHAREDIR} and the Interface subdirectory as well
DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${SEARCH_EXTRA} \
		${MODULES} ${OBJECTS}

# The order of the objects and the libraries in the link phase:
#
# - MODULES and OBJECTS come first as they use all libraries
# - LIBINTERFACE comes next as it uses all components and CON
# - component libraries come next in arbitrary order, they use CON
# - CON comes next as it uses general software libraries only
# - general software libraries come last in the order of their interdependecies

LIB1  = ${LIBIH} #	^CMP IF IH
lib1  = -lIH #		^CMP IF IH

LIB2 = ${LIBGM}	#	^CMP IF GM
lib2 = -lGM #		^CMP IF GM

LIB3 = ${LIBIM}	#	^CMP IF IM
lib3 = -lIM #		^CMP IF IM

LIB4 = ${LIBIE}	#	^CMP IF IE
lib4 = -lIE #		^CMP IF IE

LIB5 = ${LIBUA}	#	^CMP IF UA
lib5 = -lUA #		^CMP IF UA

LIB6 = ${LIBRB}	#	^CMP IF RB
lib6 = -lRB #		^CMP IF RB

${SWMFEXE}: ${MODULES} ${OBJECTS} \
		${LIBCON} ${LIBTIMING} ${LIBINTERFACE} \
		${LIB1} ${LIB2} ${LIB3} ${LIB4} ${LIB5} ${LIB6} \
		${MAKEFILE_COMP_SELECT}
	@(if [ "${COMPILE.f90}" = "xlf90" ]; then make SWMF_for_xlf; \
	elif [ "${OS}" = "Darwin" ];         then make SWMF_for_Darwin;\
	else                                      make SWMF_general; fi)

SWMF_general:
	${LINK.f90} -o ${SWMFEXE} \
		${MODULES} ${OBJECTS} \
		-L${LIBDIR} -lINTERFACE \
		${lib1} ${lib2} ${lib3} ${lib5} ${lib4} ${lib6} \
		-lCOUPLER -lLIBRARY -lTIMING -lSHARE ${LBLAS} \
		${Lflag1}

# The following target is only needed to link SWMF on Darwin 
# which has a problem with the ld linker.
# This work around creates a temporary object file from all the F90
# libraries using the -all_load option of the linker.

SWMF_for_Darwin: 
	ld -all_load -r -o LIBRARIES.o \
		${LIBSHARE} ${LIBTIMING} ${LIBLIBRARY} ${LIBCOUPLER} \
		${LIB1} ${LIB2} ${LIB3} ${LIB4} ${LIB5} ${LIB6} \
		${LIBINTERFACE}
	${LINK.f90} -o ${SWMFEXE} \
		${MODULES} ${OBJECTS} \
		LIBRARIES.o ${LBLAS} \
		${Lflag1}
	rm -f LIBRARIES.o

SWMF_for_xlf:
	rm -rf Tmp_*
	mkdir Tmp_SHARE;     cd Tmp_SHARE;     ar -x ${LIBSHARE};
	mkdir Tmp_TIMING;    cd Tmp_TIMING;    ar -x ${LIBTIMING}; 
	mkdir Tmp_LIBRARY;   cd Tmp_LIBRARY;   ar -x ${LIBLIBRARY};
	mkdir Tmp_COUPLER;   cd Tmp_COUPLER;   ar -x ${LIBCOUPLER}; 
	mkdir Tmp_INTERFACE; cd Tmp_INTERFACE; ar -x ${LIBINTERFACE}
	mkdir Tmp_LIB1;      cd Tmp_LIB1;      ar -x ${LIB1} #^CMP IF IH
	mkdir Tmp_LIB2;      cd Tmp_LIB2;      ar -x ${LIB2} #^CMP IF GM
	mkdir Tmp_LIB3;      cd Tmp_LIB3;      ar -x ${LIB3} #^CMP IF IM
	mkdir Tmp_LIB4;      cd Tmp_LIB4;      ar -x ${LIB4} #^CMP IF IE
	mkdir Tmp_LIB5;      cd Tmp_LIB5;      ar -x ${LIB5} #^CMP IF UA
	mkdir Tmp_LIB6;      cd Tmp_LIB6;      ar -x ${LIB6} #^CMP IF RB
	${LINK.f90} -o ${SWMFEXE} \
		${MODULES} ${OBJECTS} Tmp_*/*.o ${LBLAS} \
		${Lflag1}
	rm -rf Tmp_*

TESTREGISTRY: DEPEND
	cd ${SHAREDIR}; make LIB
	make test_registry.exe
	mpirun -np 4 test_registry.exe
	cat LAYOUT.in

TEST_REGISTRY_OBJ= \
	test_registry.o

test_registry.exe: DEPEND ${TEST_REGISTRY_OBJ}
	cd ${SHAREDIR}; make LIB
	cd ${LIBRARYDIR}; make LIB
	${LINK.f90} -o test_registry.exe ${SEARCH} \
		${TEST_REGISTRY_OBJ} -L${LIBDIR} -lLIBRARY -lSHARE

TEST_PHYS_OBJ = \
	test_physics.o \
	CON_methods.o \
	CON_variables.o

TESTPHYSICS: DEPEND
	cd ${SHAREDIR}; make LIB
	make test_physics.exe
	test_physics.exe

test_physics.exe: ${TEST_PHYS_OBJ}
	cd ${SHAREDIR}; make LIB
	cd ${LIBRARYDIR}; make LIB
	${LINK.f90} -o test_physics.exe ${Lflag2} ${SEARCH} \
		${TEST_PHYS_OBJ} -L${LIBDIR} -lLIBRARY -lSHARE 

#
#	cleaning
#

distclean: clean
	rm -f Makefile.DEPEND

# keep this line
