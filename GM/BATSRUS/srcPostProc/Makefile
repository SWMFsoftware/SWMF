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
#
#	cleaning
#

distclean: clean

# keep this line
