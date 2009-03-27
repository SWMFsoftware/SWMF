#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libMAGNETOGRAM.a

OBJECTS = \
	ModMagnetogram.o

OBJECTS_EXE= \
	main.o \
	ModMagHarmonics.o

main.o: ModMagHarmonics.o
dipole.o: ModMagHarmonics.o

HARMONICS:
	make ${EXE}
	@echo ${EXE} has been compiled.

EXE = ${BINDIR}/HARMONICS.exe

${EXE}: ${OBJECTS_EXE}
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} -L${LIBDIR} -lSHARE ${Lflag1}
DIPOLE:
	make DIPOLE.exe

DIPOLE.exe: dipole.o ModMagHarmonics.o
	${LINK.f90} -o DIPOLE.exe dipole.o ModMagHarmonics.o \
	-L${LIBDIR} -lSHARE ${Lflag1}
LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

distclean: clean

### TESTS FOR HARMONICS ###

test: 
	@rm -f *.diff
	-@(make test_harmonics)
	@ls -l *.diff

test_harmonics:
	@echo "test_harmonics_compile..." > test_harmonics.diff
	make HARMONICS
	rm -f HARMONICS.exe
	ln -s ${BINDIR}/HARMONICS.exe HARMONICS.exe
	rm -f DIPOLE.exe
	make DIPOLE
	./DIPOLE.exe
	@echo "test_harmonics_run..." >> test_harmonics.diff
	./HARMONICS.exe > runlog
	@echo "test_harmonics_check..." >> test_harmonics.diff
	make test_harmonics_check

test_harmonics_check:
	${SCRIPTDIR}/DiffNum.pl -t -r=1e-5 -a=1e-12 \
		harmonics.dat \
		test.dat > test_harmonics.diff
	ls -l test_harmonics.diff
	rm -f fitsfile.dat harmonics.dat 
