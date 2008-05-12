#^CFG COPYRIGHT UM

SHELL =/bin/sh

default : HEIDI.exe

include ../Makefile.def
include ../Makefile.conf
include Makefile.DEPEND
#
#  Makefile for 
#



MODULES = \
	ModHeidiSize.o\
	ModHeidiMain.o\
	ModHeidiCurrents.o\
	ModHeidiDGCPM.o\
	ModHeidiDrifts.o\
	ModHeidiIO.o\
	ModHeidiWaves.o\
	ModSize.o\
	ModMain.o\
	ModIonosphere.o\
	ModConductance.o\
	ModPhysics.o\
	ModIO.o\
	precision.o \
	io.o \
        CSR.o \
        MSR.o \
        LUCSRMSR.o \
        BandLU.o \
        matrix-arithmetic.o \
        flop_module.o \
        iteration-defaults.o \
        GMRES.o \
        ilu.o

OBJ_IONO = \
        iono_conductances.o\
        iono_solver.o\
        iono_coupling.o\
	read_ring_current.o\
	write_ring_current.o\
	read_amie.o\
	epencalc.o\
	convert_time_routines.o\
	iono_init.o

OBJECTS_010 = \
	heidi_setup.o\
	heidi_icbc.o\
	heidi_coefficients.o\
	heidi_operators.o\
	heidi_currents.o\
	heidi_output.o\
	heidi_waves.o\
	recipes.o\
	pbo.o\
	wgrdar.o\
	get_potential.o

MODULES_IE = \
	ModConductance.o\
	ModKind.o\
	ModAMIE_Interface.o\
	ModErrors.o\
	ModIE_Interface.o\
	ModExtras.o\
	ModFiles.o

OBJECTS_IE = \
	AMIE_Library.o\
	GM_Library.o\
	IE_End.o\
	IE_Initialize.o\
	IE_Library.o\
	IO_Library.o\
	IE_set_inputs.o\
	readAMIEoutput.o\
	wei01.o\
	wei96.o\
	weicom.o\
	mh86.o\
	hmr89.o\
	iz94.o\
	merge_str.o\
	time_routines.o\
	ihp.o

OBJECTS = ${OBJ_IONO} ${OBJECTS_IE} ${OBJECTS_010} ${MODULES_IE} ${MODULES}


DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} ${OBJECTS}

IE_LIB = libIE.a

${IE_LIB}: ${MODULES_IE} ${OBJECTS_IE}
	rm -f ${IE_LIB}
	${AR} ${IE_LIB} ${MODULES_IE} ${OBJECTS_IE}

# Object files not included into the library
#Library should be compiled first, because modules are used in main.
main.o : ${IE_LIB}

standalone : main.o
	${LINK.f90} -o test.exe main.o $ -L. -lIE ${Lflag1}

#	${LINK.f90} -o test.exe ${MODULES} ${OBJECTS} main.o 
interface: ${IE_LIB} Interface.o
	${LINK.f90} -o Interface.exe Interface.o ${MODULES} -L. -lIE ${Lflag1}

IONO_LIB = libIONO.a

${IONO_LIB}: ${MODULES} ${OBJECTS_IONO} ${IE_LIB}
	cp -f ${IE_LIB} ${IONO_LIB}
	${AR} ${IONO_LIB} ${MODULES} ${OBJECTS_IONO}

MY_LIB = libHEIDI.a

${MY_LIB}: ${OBJECTS}
	${AR} ${MY_LIB} ${OBJECTS} 

LIB : 
	make ${MY_LIB}

#IONOSPHERE

IONO:	IONO.exe

iono.o : ${IONO_LIB} 
#main.o : ${IONO_LIB}

IONO.exe :	iono.o #main.o
	${LINK.f90} -o IONO.exe iono.o -L. -lIONO ${Lflag1} 

HEIDI:  LIB
	make HEIDI.exe

SEARCH_EXTRA = -I${LIBRARYDIR} -I${COUPLERDIR}

HEIDI.exe: ${OBJECTS} heidi_main.o
	${LINK.f90} -o ${BINDIR}/HEIDI.exe heidi_main.o ${OBJECTS}

#clean:
#	/bin/rm -f *.o *.mod
#	cd ../srcInterface/srcIONO ; make clean
#	cd ../srcInterface/srcIE ; make clean

distclean: clean
	rm -f Makefile.DEPEND Makefile.RULES


# HEADERS = \
# 	numv.h
# 
# OBJECTS = \
# 	mram03.o\
# 	fct-source.o\
# 	msis86.o\
# 	iri87.o\
# 	recipes.o\
# 	plane.o\
# 	wgrdar.o
# 
# OBJECTS_4 = \
# 	mram04.o\
# 	recipes.o\
# 	pbo.o\
# 	wgrdar.o
# 
# OBJECTS_5 = \
# 	mram05.o\
# 	recipes.o\
# 	pbo.o\
# 	wgrdar.o\
# 	get_potential.o
# 
# ${OBJECTS}:	${HEADERS}
# 
# mram03.exe:	${MODULES} ${OBJECTS} ${HEADERS}
# 	$(FORTRAN) $(LFLAGS) -o mram03.exe ${MODULES} ${OBJECTS} ${IONO}
# 
# mram05:
# 	cd ../srcInterface/srcIONO ; make LIB
# 	cd ../srcInterface/srcIE ; make LIB
# 	make mram05.exe
# 
# mram05.exe:	${MODULES} ${OBJECTS_5} ${HEADERS} ../srcInterface/srcIE/libIE.a
# 	$(FORTRAN) $(LFLAGS) -o mram05.exe ${MODULES} ${OBJECTS_5} ${IONO}\
# 	-L../srcInterface/srcIE -lIE
# 
# mram04:
# 	cd ../srcInterface/srcIONO ; make
# 	make mram04.exe
# 
# mram04.exe:	${MODULES} ${OBJECTS_4} ${HEADERS}
# 	$(FORTRAN) $(LFLAGS) -o mram04.exe ${MODULES} ${OBJECTS_4} ${IONO}
# 
# mram02_epen:  mram02_epen.o fct-source.o iri87.o msis86.o recipes.o plane.o wgrdar.o numv.h
# 	$(FORTRAN) $(CFLAGS) -o mram02_epen mram02_epen.o fct-source.o msis86.o iri87.o recipes.o plane.o wgrdar.o
# 
# mram01:  mram01.o fct-source.o iri87.o msis86.o recipes.o plane.o wgrdar.o numv.h
# 	$(FORTRAN) $(CFLAGS) -o mram01 mram01.o fct-source.o msis86.o iri87.o recipes.o plane.o wgrdar.o
# 
# mram02:  mram02.o fct-source.o iri87.o msis86.o recipes.o plane.o wgrdar.o numv.h
# 	$(FORTRAN) $(CFLAGS) -o mram02 mram02.o fct-source.o msis86.o iri87.o recipes.o plane.o wgrdar.o
# 
# mram01_ios:  mram01_ios.o fct-source.o iri87.o msis86.o recipes.o plane.o wgrdar.o numv.h
# 	$(FORTRAN) $(CFLAGS) -o mram01_ios mram01_ios.o fct-source.o msis86.o iri87.o recipes.o plane.o wgrdar.o
# 
# mram01_ena:  mram01_ena.o fct-source.o iri87.o msis86.o recipes.o plane.o wgrdar.o numv.h
# 	$(FORTRAN) $(CFLAGS) -o mram01_ena mram01_ena.o fct-source.o msis86.o iri87.o recipes.o plane.o wgrdar.o
# 
# thermal:  therm1.o fct-source.o iri87.o msis86.o recipes.o plane.o numv.h
# 	$(FORTRAN) $(CCFLAGS) -o thermal therm1.o fct-source.o msis86.o iri87.o recipes.o plane.o
# 
# SOURCES = cross7.f cross7c.f diff7b.f espec7.f geomag.f getf107.f glowsub7.f iri90.f maincpl7c.f math.f msis90.f rcolum.f solzen.f ssflux.f stet7c.f stetprn7c.f stetsub7.f substet7c.f thermstet7c.f wave7.f
# 
# all:	
# 	$(FORTRAN) $(CFLAGS) -c $(SOURCES)
