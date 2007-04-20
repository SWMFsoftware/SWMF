#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf

#
# Post processing executables
#

PIDL:	PostIDL.o
	${LINK.f90} -o ${BINDIR}/PostIDL.exe PostIDL.o ${Lflag2}

PSPH:	PostSPH.o
	${LINK.f90} -o ${BINDIR}/PostSPH.exe PostSPH.o ${Lflag1}

CRST:	ConvertRestart.o
	${LINK.f90} -o ${BINDIR}/ConvertRestart.exe ConvertRestart.o ${Lflag2}

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

#
#	cleaning
#

distclean: clean

# keep this line
