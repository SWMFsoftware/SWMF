#^CFG COPYRIGHT UM

SHELL =/bin/sh

# CON also needs modules from the Interface subdirectory
SEARCH_EXTRA = -I${INTERFACEDIR} -I${LIBRARYDIR} -I${COUPLERDIR}

include ../../../Makefile.def
include ../../../Makefile.conf
include Makefile.DEPEND

#######################
# Software libraries: #
#######################

LIBTIMING = ${LIBDIR}/libTIMING.a

#######################
# Shared libraries by #
# CON and componets:  #
#######################

LIBSHARE   = ${LIBDIR}/libSHARE.a

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
	@cd ${TIMINGDIR};    make LIB
	@cd ${SHAREDIR};     make LIB
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

LIB1 = ${LIBIH} #^CMP IF IH
lib1 = -lIH     #^CMP IF IH
LIB2 = ${LIBGM} #^CMP IF GM
lib2 = -lGM     #^CMP IF GM
LIB3 = ${LIBIM} #^CMP IF IM
lib3 = -lIM     #^CMP IF IM
LIB4 = ${LIBIE} #^CMP IF IE
lib4 = -lIE     #^CMP IF IE
LIB5 = ${LIBUA} #^CMP IF UA
lib5 = -lUA     #^CMP IF UA

${SWMFEXE}: 	${MODULES} ${OBJECTS} \
		${LIBCON} ${LIBTIMING} ${LIBINTERFACE} \
		${LIB1} ${LIB2} ${LIB3} ${LIB4} ${LIB5} \
		${MAKEFILE_COMP_SELECT}
	${LINK.f90} -o ${SWMFEXE} ${SEARCH} \
		${MODULES} ${OBJECTS} \
		-L${LIBDIR} -lINTERFACE \
		${lib1} ${lib2} ${lib3} ${lib5} ${lib4} \
		-lCOUPLER -lLIBRARY -lSHARE -lTIMING ${LBLAS} \
		${Lflag1}

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


