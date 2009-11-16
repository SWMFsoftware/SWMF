#^CFG COPYRIGHT UM

SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf

# Other requireed libraries
LIBSHARE  = ${LIBDIR}/libSHARE.a


#
# Post processing executables
#

${BINDIR}/PostPwIdl.exe: ${LIBSHARE} PostPwIdl.o
	${LINK.f90} -o ${BINDIR}/PostPwIdl.exe PostPwIdl.o -L${LIBDIR} -lSHARE ${Lflag2}

PwIDL:
	make ${BINDIR}/PostPwIdl.exe

#
#	cleaning
#

distclean: clean
	rm -f *.o

# keep this line
