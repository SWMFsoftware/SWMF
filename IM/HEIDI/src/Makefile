#
#  Makefile for 
#

FORTRAN = f95

#FORTRAN = f77
#FORTRAN = f95
#FORTRAN = mpif90

#CFLAGS = -r1 -O3 -c
#CFLAGS = -K +es -c
#CFLAGS = -r2 -N132

# for debugging:
#CFLAGS = -132 -w -dusty -gline -maxcontin=100 -C=array -c

# for runs:
#CFLAGS = -132 -maxcontin=100 -c 
CFLAGS = -132 -c -dusty -w -gline -I../srcIE -maxcontin=50 -g -O3

# for compiling smaller programs:
CCFLAGS = -132 -maxcontin=100 

# for linking:
LFORTRAN = mpif90
LFLAGS = -g

MODULES = \
	ModHeidiSize.o\
	ModHeidiMain.o\
	ModHeidiCurrents.o\
	ModHeidiDGCPM.o\
	ModHeidiDrifts.o\
	ModHeidiIO.o\
	ModHeidiWaves.o\
	../srcIONO/ModIonosphere.o\
	../srcIONO/ModConductance.o\
	../srcIONO/ModPhysics.o\
	../srcIONO/ModSize.o\
	../srcIONO/ModMain.o\
	../srcIONO/ModIO.o

OBJECTS_010 = \
	heidi_main_010.o\
	heidi_setup_011.o\
	heidi_icbc_010.o\
	heidi_coefficients_010.o\
	heidi_operators_010.o\
	heidi_currents_010.o\
	heidi_output_010.o\
	heidi_waves_010.o\
	recipes.o\
	pbo.o\
	wgrdar.o\
	get_potential.o

IONO = \
	../srcIONO/precision.o \
	../srcIONO/io.o \
	../srcIONO/CSR.o \
	../srcIONO/MSR.o \
	../srcIONO/LUCSRMSR.o \
	../srcIONO/BandLU.o \
	../srcIONO/matrix-arithmetic.o \
	../srcIONO/flop_module.o \
	../srcIONO/iteration-defaults.o \
	../srcIONO/GMRES.o \
	../srcIONO/ilu.o\
	../srcIONO/iono_conductances.o\
	../srcIONO/iono_solver.o\
	../srcIONO/iono_coupling.o\
	../srcIONO/read_ring_current.o\
	../srcIONO/write_ring_current.o\
	../srcIONO/read_amie.o\
	../srcIONO/epencalc.o\
	../srcIONO/convert_time_routines.o\
	../srcIONO/iono_init.o

heidi:
	cd ../srcIONO ; make LIB
	cd ../srcIE ; make LIB
	make heidi010.exe

heidi010:
	cd ../srcIONO ; make LIB
	cd ../srcIE ; make LIB
	make heidi010.exe

heidi010.exe:	${MODULES} ${OBJECTS_010} ${IONO} ../srcIE/libIE.a
	$(LFORTRAN) $(LFLAGS) -o heidi010.exe ${MODULES} ${OBJECTS_010} ${IONO}\
	-L../srcIE -lIE

.SUFFIXES:
.SUFFIXES: .f90 .f .o

.f.o:
	${FORTRAN} $(CFLAGS) $<

.f90.o:
	${FORTRAN} $(CFLAGS) $<

clean:
	/bin/rm -f *.o *.mod
	cd ../srcIONO ; make clean
	cd ../srcIE ; make clean

#
# Older versions :
#

HEADERS = \
	numv.h

OBJECTS = \
	mram03.o\
	fct-source.o\
	msis86.o\
	iri87.o\
	recipes.o\
	plane.o\
	wgrdar.o

OBJECTS_4 = \
	mram04.o\
	recipes.o\
	pbo.o\
	wgrdar.o

OBJECTS_5 = \
	mram05.o\
	recipes.o\
	pbo.o\
	wgrdar.o\
	get_potential.o

${OBJECTS}:	${HEADERS}

mram03.exe:	${MODULES} ${OBJECTS} ${HEADERS}
	$(FORTRAN) $(LFLAGS) -o mram03.exe ${MODULES} ${OBJECTS} ${IONO}

mram05:
	cd ../srcIONO ; make LIB
	cd ../srcIE ; make LIB
	make mram05.exe

mram05.exe:	${MODULES} ${OBJECTS_5} ${HEADERS} ../srcIE/libIE.a
	$(FORTRAN) $(LFLAGS) -o mram05.exe ${MODULES} ${OBJECTS_5} ${IONO}\
	-L../srcIE -lIE

mram04:
	cd ../srcIONO ; make
	make mram04.exe

mram04.exe:	${MODULES} ${OBJECTS_4} ${HEADERS}
	$(FORTRAN) $(LFLAGS) -o mram04.exe ${MODULES} ${OBJECTS_4} ${IONO}

mram02_epen:  mram02_epen.o fct-source.o iri87.o msis86.o recipes.o plane.o wgrdar.o numv.h
	$(FORTRAN) $(CFLAGS) -o mram02_epen mram02_epen.o fct-source.o msis86.o iri87.o recipes.o plane.o wgrdar.o

mram01:  mram01.o fct-source.o iri87.o msis86.o recipes.o plane.o wgrdar.o numv.h
	$(FORTRAN) $(CFLAGS) -o mram01 mram01.o fct-source.o msis86.o iri87.o recipes.o plane.o wgrdar.o

mram02:  mram02.o fct-source.o iri87.o msis86.o recipes.o plane.o wgrdar.o numv.h
	$(FORTRAN) $(CFLAGS) -o mram02 mram02.o fct-source.o msis86.o iri87.o recipes.o plane.o wgrdar.o

mram01_ios:  mram01_ios.o fct-source.o iri87.o msis86.o recipes.o plane.o wgrdar.o numv.h
	$(FORTRAN) $(CFLAGS) -o mram01_ios mram01_ios.o fct-source.o msis86.o iri87.o recipes.o plane.o wgrdar.o

mram01_ena:  mram01_ena.o fct-source.o iri87.o msis86.o recipes.o plane.o wgrdar.o numv.h
	$(FORTRAN) $(CFLAGS) -o mram01_ena mram01_ena.o fct-source.o msis86.o iri87.o recipes.o plane.o wgrdar.o

thermal:  therm1.o fct-source.o iri87.o msis86.o recipes.o plane.o numv.h
	$(FORTRAN) $(CCFLAGS) -o thermal therm1.o fct-source.o msis86.o iri87.o recipes.o plane.o

SOURCES = cross7.f cross7c.f diff7b.f espec7.f geomag.f getf107.f glowsub7.f iri90.f maincpl7c.f math.f msis90.f rcolum.f solzen.f ssflux.f stet7c.f stetprn7c.f stetsub7.f substet7c.f thermstet7c.f wave7.f

all:	
	$(FORTRAN) $(CFLAGS) -c $(SOURCES)
