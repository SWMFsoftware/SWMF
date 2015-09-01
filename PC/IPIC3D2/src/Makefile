include ../Makefile.def
include ./Makefile.local

OBJECTS = \
	bc/BcFields3D.o \
	bc/BcParticles.o \
	communication/ComBasic3D.o\
	communication/ComInterpNodes3D.o\
	communication/ComNodes3D.o\
	communication/ComParser3D.o\
	communication/ComParticles3D.o\
	ConfigFile/src/ConfigFile.o\
	fields/BCStructure.o\
	fields/EMfields3D.o\
	grids/Grid3DCU.o\
	inputoutput/Collective.o\
	inputoutput/phdf5.o\
	inputoutput/ParallelIO.o\
	inputoutput/Restart3D.o\
	main/iPic3Dlib.o\
	mathlib/EllipticF.o\
	particles/Particles3D.o\
	particles/Particles3Dcomm.o\
	performances/Timing.o\
	PSKOutput3D/PSKhdf5adaptor.o\
	solvers/CG.o\
	solvers/GMRES.o\
	utility/asserts.o\
	utility/debug.o

OBJECTS_EXE = \
	iPic3D.o

EXE = ${BINDIR}/IPIC3D.exe

SEARCH_C = -Iinclude

FLAGC_EXTRA = ${SEARCH_C} \
	${FLAGC_LOCAL}

include ../Makefile.conf
include Makefile.DEPEND


DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH_C} ${OBJECTS} ${OBJECTS_EXE}

#
#       Making executables and libraries
#

#
# Library for framework and also for stand alone code
#

MY_LIB  = ../libIPIC3D/libIPIC3D.a

LIB:    DEPEND
	$(MAKE) ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS} 
	${AR} ${MY_LIB} ${OBJECTS}


IPIC3D:
	$(MAKE) LIB
	$(MAKE) ${OBJECTS_EXE}
	h5pcc -o ${EXE} ${MY_LIB} ${OBJECTS_EXE} ${CPPLIB} 
	@echo ' '
	@echo Program IPIC3D has been brought up to date.
	@echo ' '       

clean:
	rm -f ${OBJECTS} ${OBJECTS_EXE} ${MY_LIB}

distclean: clean
	rm -f Makefile.DEPEND Makefile.local Makefile.RULES
