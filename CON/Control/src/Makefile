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

# CON uses modules from ${MODDIR} and the Interface subdirectory as well
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

#^CMP IF IH BEGIN
LIB1  = ${LIBIH}
lib1  = -lIH
doto1 = ${LIBDIR}/IH.o

${doto1}:	${LIBIH}
	ld -o ${doto1} ${LIBIH} -all_load -r

#^CMP END IH
#^CMP IF GM BEGIN
LIB2 = ${LIBGM}
lib2 = -lGM
doto2 = ${LIBDIR}/GM.o

${doto2}:	${LIBGM}
	ld -o ${doto2} ${LIBGM} -all_load -r

#^CMP END GM
#^CMP IF IM BEGIN
LIB3 = ${LIBIM}
lib3 = -lIM
doto3 = ${LIBDIR}/IM.o

${doto3}:	${LIBIM}
	ld -o ${doto3} ${LIBIM} -all_load -r

#^CMP END IM
#^CMP IF IE BEGIN
LIB4 = ${LIBIE}
lib4 = -lIE
doto4 = ${LIBDIR}/IE.o

${doto4}:	${LIBIE}
	ld -o ${doto4} ${LIBIE} -all_load -r

#^CMP END IE
#^CMP IF UA BEGIN
LIB5 = ${LIBUA}
lib5 = -lUA
doto5 = ${LIBDIR}/UA.o

${doto5}:	${LIBUA}
	ld -o ${doto5} ${LIBUA} -all_load -r

#^CMP END UA

dotoINTERFACE = ${LIBDIR}/INTERFACE.o
dotoTIMING = ${LIBDIR}/TIMING.o
dotoCOUPLER = ${LIBDIR}/COUPLER.o
dotoLIBRARY = ${LIBDIR}/LIBRARY.o
dotoSHARE = ${LIBDIR}/SHARE.o

${dotoINTERFACE}:	${LIBINTERFACE}
	ld -o ${dotoINTERFACE} ${LIBINTERFACE} -all_load -r

${dotoTIMING}:	${LIBTIMING}
	ld -o ${dotoTIMING} ${LIBTIMING} -all_load -r

${dotoCOUPLER}:	${LIBCOUPLER}
	ld -o ${dotoCOUPLER} ${LIBCOUPLER} -all_load -r

${dotoLIBRARY}:	${LIBLIBRARY}
	ld -o ${dotoLIBRARY} ${LIBLIBRARY} -all_load -r

${dotoSHARE}:	${LIBSHARE}
	ld -o ${dotoSHARE} ${LIBSHARE} -all_load -r

DotOs:	  ${doto1} ${doto2} ${doto3} ${doto5} ${doto4} \
	  ${dotoINTERFACE} ${dotoTIMING} \
	  ${dotoCOUPLER} ${dotoLIBRARY} ${dotoSHARE}

${SWMFEXE}: 	${MODULES} ${OBJECTS} \
		${LIBCON} ${LIBTIMING} ${LIBINTERFACE} \
		${LIB1} ${LIB2} ${LIB3} ${LIB4} ${LIB5} \
		${MAKEFILE_COMP_SELECT}
	make DotOs
	@echo ' '
	@(if [ ${OS} == "Darwin" ]; then \
	    echo \
	    ${LINK.f90} -o ${SWMFEXE} ${SEARCH} \
		${MODULES} ${OBJECTS} \
		${doto1} ${doto2} ${doto3} ${doto5} ${doto4} \
		${dotoINTERFACE} ${dotoCOUPLER} ${dotoLIBRARY} \
		${dotoTIMING} ${dotoSHARE} ${LBLAS} \
		${Lflag1}; \
	    ${LINK.f90} -o ${SWMFEXE} ${SEARCH} \
		${MODULES} ${OBJECTS} \
		${doto1} ${doto2} ${doto3} ${doto5} ${doto4} \
		${dotoINTERFACE} ${dotoCOUPLER} ${dotoLIBRARY} \
		${dotoTIMING} ${dotoSHARE} ${LBLAS} \
		${Lflag1}; \
	else \
	    echo \
	    ${LINK.f90} -o ${SWMFEXE} ${SEARCH} \
		${MODULES} ${OBJECTS} \
		-L${LIBDIR} -lINTERFACE \
		${lib1} ${lib2} ${lib3} ${lib5} ${lib4} \
		-lCOUPLER -lLIBRARY -lTIMING -lSHARE ${LBLAS} \
		${Lflag1}; \
	    ${LINK.f90} -o ${SWMFEXE} ${SEARCH} \
		${MODULES} ${OBJECTS} \
		-L${LIBDIR} -lINTERFACE \
		${lib1} ${lib2} ${lib3} ${lib5} ${lib4} \
		-lCOUPLER -lLIBRARY -lTIMING -lSHARE ${LBLAS} \
		${Lflag1}; \
	fi);


TESTREGISTRY: DEPEND
	cd ${MODDIR}; make LIB
	make test_registry.exe
	mpirun -np 4 test_registry.exe
	cat LAYOUT.in

TEST_REGISTRY_OBJ= \
	test_registry.o

test_registry.exe: DEPEND ${TEST_REGISTRY_OBJ}
	cd ${MODDIR}; make LIB
	${LINK.f90} -o test_registry.exe ${SEARCH} \
		${TEST_REGISTRY_OBJ} -L../lib -lCON

TEST_PHYS_OBJ = \
	test_physics.o \
	CON_methods.o

TESTPHYSICS: DEPEND
	cd ${MODDIR}; make LIB
	make test_physics.exe
	test_physics.exe

test_physics.exe: ${TEST_PHYS_OBJ} ${CONLIB}
	${LINK.f90} -o test_physics.exe ${Lflag2} ${SEARCH} \
		${TEST_PHYS_OBJ} -L../lib -lCON

#
#	cleaning
#

distclean: clean
	rm -f Makefile.DEPEND

# keep this line
