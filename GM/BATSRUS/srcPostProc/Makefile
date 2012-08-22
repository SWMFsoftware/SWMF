#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
include Makefile.RULES

# Other requireed libraries
LIBSHARE  = ${LIBDIR}/libSHARE.a


#
# Post processing executables
#

${BINDIR}/PostIDL.exe: ${LIBSHARE} PostIDL.o
	${LINK.f90} -o ${BINDIR}/PostIDL.exe PostIDL.o -L${LIBDIR} -lSHARE ${Lflag2}

PIDL:
	make ${BINDIR}/PostIDL.exe

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

#
#	cleaning
#

distclean: clean

# keep this line
