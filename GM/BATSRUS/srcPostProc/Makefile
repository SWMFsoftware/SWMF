#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf

#
# Post processing executables
#

PIDL:	PostIDL.o
	${LINK.f90} -o ${BINDIR}/PostIDL.exe PostIDL.o ${Lflag2}

#^CFG IF NOT SIMPLE BEGIN
PSPH:	PostSPH.o
	${LINK.f90} -o ${BINDIR}/PostSPH.exe PostSPH.o ${Lflag1}
#^CFG END SIMPLE

#
#	cleaning
#

distclean: clean

# keep this line
