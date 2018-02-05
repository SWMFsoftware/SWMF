#  Copyright (C) 2002 Regents of the University of Michigan,
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
include Makefile.RULES

# Other required libraries
LIBSHARE  = ${LIBDIR}/libSHARE.a

#
# Post processing executables
#
${BINDIR}/SPECTRUM.exe: ${LIBDIR}/libSHARE.a spectrum.o
	${LINK.f90} -o ${BINDIR}/SPECTRUM.exe spectrum.o -L${LIBDIR} -lSHARE ${Lflag1} 

SPECTRUM:
	make ${BINDIR}/SPECTRUM.exe
	rm -f ./SPECTRUM.exe
	ln -s ${BINDIR}/SPECTRUM.exe .

${BINDIR}/SNAPSHOT.exe: ${LIBSHARE} select_snapshot.o
	${LINK.f90} -o ${BINDIR}/SNAPSHOT.exe select_snapshot.o -L${LIBDIR} -lSHARE ${Lflag2}

SNAPSHOT:
	make ${BINDIR}/SNAPSHOT.exe

${BINDIR}/ConvertRestart.exe:	ConvertRestart.o
	${LINK.f90} -o ${BINDIR}/ConvertRestart.exe ConvertRestart.o ${Lflag2}

CRST:
	make ${BINDIR}/ConvertRestart.exe

OBJECTS1 = \
	CON_geopack_internal.o\
	ModUT.o \
	${SHAREDIR}/ModNumConst.o \
	${SHAREDIR}/ModTimeConvert.o \
	${SHAREDIR}/ModConst.o

earth_traj.o: ${OBJECTS1}

EARTH_TRAJ: ${OBJECTS1} earth_traj.o
	${LINK.f90} ${SEARCH} -o ${BINDIR}/EARTH_TRAJ.exe earth_traj.o \
		${OBJECTS1} ${Lflag2}
TIME_CONV:  ${OBJECTS1} time_convert.o
	    ${LINK.f90} ${SEARCH} -o ${BINDIR}/TIME_CONV.exe time_convert.o \
		${OBJECTS1} ${Lflag2}

TIME_TOCRNO:  ${OBJECTS1} time_to_cr_no.o
	    ${LINK.f90} ${SEARCH} -o ${BINDIR}/TIME_TOCRNO.exe time_to_cr_no.o \
		${OBJECTS1} ${Lflag2}

#
# test for spectrum
#
test_spectrum:
	@echo "test_spectrum_compile..." > ../test_spectrum.diff
	$(MAKE) test_spectrum_compile
	@echo "test_spectrum_run..." >> ../test_spectrum.diff
	$(MAKE) test_spectrum_run
	@echo "test_spectrum_check..." >> ../test_spectrum.diff
	$(MAKE) test_spectrum_check

SPECTRUM_chianti_tbl.dat: ../data/SPECTRUM/SPECTRUM_chianti_tbl.dat.gz
	gunzip -c ../data/SPECTRUM/SPECTRUM_chianti_tbl.dat.gz \
	> SPECTRUM_chianti_tbl.dat

test_spectrum_compile: SPECTRUM_chianti_tbl.dat
	rm -f SPECTRUM.in ../test_spectrum.diff \
		test-spectrum.log test-spectrum.out test-spectrum.eps
	cp test-SPECTRUM.in SPECTRUM.in
	make SPECTRUM

test_spectrum_run:
	./SPECTRUM.exe > test-spectrum.log

test_spectrum_check:
	${SCRIPTDIR}/DiffNum.pl -t -r=1e-5 -a=1e-12 \
		test-spectrum.out \
		test-spectrum.ref > ../test_spectrum.diff
	@ls -l ../test_spectrum.diff

#
#	cleaning
#
clean: cleanfiles

distclean: clean
	rm -f SPECTRUM.exe SPECTRUM.in test-spectrum.log \
		test-spectrum.out test-label.out SPECTRUM_chianti_tbl.dat

