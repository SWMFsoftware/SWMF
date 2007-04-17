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
	ModUT.o

earth_traj_main.o:\
	CON_geopack_internal.o \
	ModUT.o \
	${SHAREDIR}/ModNumConst.o \
	${SHAREDIR}/ModTimeConvert.o \
	${SHAREDIR}/ModConst.o

EARTH_TRAJ = ${BINDIR}/earth_traj.exe

earthtraj:
	make ${EARTH_TRAJ}
	@echo
	@echo ${EARTH_TRAJ} has been brought up to date.
	@echo

${EARTH_TRAJ}: ${OBJECTS1} earth_traj_main.o
	${LINK.f90} ${SEARCH} -o ${EARTH_TRAJ} ${OBJECTS1} earth_traj_main.o  ${LIBDIR}/libSHARE.a

#
#	cleaning
#

distclean: clean

# keep this line
