include ../Makefile.def
include ./Makefile.local

OBJECTS = \
	bc/BcFields3D.o \
	communication/Com3DNonblk.o\
	communication/VCtopology3D.o\
	ConfigFile/src/ConfigFile.o\
	fields/Moments.o\
	fields/EMfields3D.o\
	grids/Grid3DCU.o\
	inputoutput/phdf5.o\
	inputoutput/ParallelIO.o\
	inputoutput/OutputWrapperFPP.o\
	inputoutput/WriteOutputParallel.o\
	main/iPic3Dlib.o\
	main/Collective.o\
	main/Parameters.o\
	mathlib/EllipticF.o\
	particles/Particles3D.o\
	particles/Particles3Dcomm.o\
	performances/Timing.o\
	PSKOutput3D/PSKhdf5adaptor.o\
	solvers/CG.o\
	solvers/GMRES.o\
	utility/asserts.o\
	utility/debug.o\
	utility/IDgenerator.o\
	utility/Basic.o\
	utility/TimeTasks.o\
	utility/MPIdata.o\
	utility/errors.o

${OBJECTS}: Makefile.local

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
	${COMPILE.mpicxx} -o ${EXE} ${MY_LIB} ${OBJECTS_EXE} \
		-L../libIPIC3D/ -lIPIC3D ${CPPLIB} 
	@echo ' '
	@echo Program IPIC3D has been brought up to date.
	@echo ' '       

clean:
	rm -f ${OBJECTS} ${OBJECTS_EXE} ${MY_LIB}

distclean: clean
	rm -f Makefile.DEPEND Makefile.local Makefile.RULES
