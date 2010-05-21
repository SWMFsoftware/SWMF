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

POTENTIAL: potential_field.o
	${LINK.f90} -o POTENTIAL.exe potential_field.o \
		-L${LIBDIR} -lSHARE ${Lflag1}

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
	rm -f *.dat *.out *.log *.diff

### TESTS FOR HARMONICS ###

test: 
	@rm -f *.diff
	-@(make test_harmonics)
	-@(make test_potential)
	@ls -l *.diff

test_harmonics:
	@echo "test_harmonics_compile..." > test_harmonics.diff
	make HARMONICS
	rm -f HARMONICS.exe
	ln -s ${BINDIR}/HARMONICS.exe HARMONICS.exe
	make DIPOLE
	./DIPOLE.exe
	@echo "test_harmonics_run..." >> test_harmonics.diff
	./HARMONICS.exe > harmonics.log
	@echo "test_harmonics_check..." >> test_harmonics.diff
	make test_harmonics_check

test_harmonics_check:
	${SCRIPTDIR}/DiffNum.pl -t -r=1e-5 -a=1e-12 \
		harmonics.dat \
		harmonics.ref > test_harmonics.diff
	ls -l test_harmonics.diff

###############################################################################

test_potential:
	@echo "test_potential_compile..." > test_potential.diff
	make POTENTIAL
	make DIPOLE
	./DIPOLE.exe
	@echo "test_potential_run..." >> test_potential.diff
	./POTENTIAL.exe > potential.log
	@echo "test_potential_check..." >> test_potential.diff
	make test_potential_check

test_potential_check:
	${SCRIPTDIR}/DiffNum.pl -t -r=1e-5 -a=1e-12 \
		potential.log \
		potential.ref > test_potential.diff
	ls -l test_potential.diff
