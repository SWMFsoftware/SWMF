#^CFG COPYRIGHT UM

include ../../../Makefile.def
include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libINDICES.a

OBJECTS = \
	main.o \
	ModMagHarmonics.o

main.o: ModMagHarmonics.o

HARMONICS:
	make ${EXE}
	@echo ${EXE} has been compiled.

EXE = ${BINDIR}/HARMONICS.exe

${EXE}: ${OBJECTS}
	${LINK.f90} -o ${EXE} ${OBJECTS} ${Lflag1}

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
	cp test.in fitsfile.dat
	@echo "test_harmonics_run..." >> test_harmonics.diff
	mpirun -np 2 HARMONICS.exe > runlog
	@echo "test_harmonics_check..." >> test_harmonics.diff
	make test_harmonics_check

test_harmonics_check:
	${SCRIPTDIR}/DiffNum.pl -t -r=1e-5 -a=1e-12 \
		harmonics.dat \
		test.out > test_harmonics.diff
	ls -l test_harmonics.diff
	rm -f fitsfile.dat harmonics.dat 
