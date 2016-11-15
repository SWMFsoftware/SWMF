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
${BINDIR}/SPECTRUM.exe:	spectrum.o ${LIBDIR}/libSHARE.a
	${LINK.f90} -o ${BINDIR}/SPECTRUM.exe spectrum.o -L${LIBDIR} -lSHARE ${Lflag1} 

SPECTRUM:
	make ${BINDIR}/SPECTRUM.exe
	ln -s ${BINDIR}/SPECTRUM.exe .

${BINDIR}/PostSPH.exe:	PostSPH.o
	${LINK.f90} -o ${BINDIR}/PostSPH.exe PostSPH.o ${Lflag1}

PSPH:
	make ${BINDIR}/PostSPH.exe

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
	@echo "test_spectrum compiles..." > test_spectrum.diff
	rm -f SPECTRUM.exe SPECTRUM.in test_spectrum.diff spectrum.log spectrum.out
	make SPECTRUM
	@echo "test_spectrum copies input files..." >> test_spectrum.diff
	cp SPECTRUM.test.in SPECTRUM.in
	@echo "test_spectrum runs..." > test_spectrum.diff
	${BINDIR}/SPECTRUM.exe > spectrum.log
	@echo "test_spectrum checks result" >> test_spectrum.diff
	${SCRIPTDIR}/DiffNum.pl -t -r=1e-5 -a=1e-12 \
		spectrum.out \
		spectrum.ref > test_spectrum.diff
	@ls -l test_spectrum.diff

#
#	cleaning
#
clean: cleanfiles

distclean: clean
	rm -f SPECTRUM.exe SPECTRUM.in test_spectrum.diff spectrum.log spectrum.out

# keep this line
