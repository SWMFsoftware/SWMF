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
#	cleaning
#

clean: cleanfiles

distclean: clean

# keep this line
