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

${BINDIR}/PostPwIdl_3D.exe: ${LIBSHARE} PostPwIdl_3D.o
	${LINK.f90} -o ${BINDIR}/PostPwIdl_3D.exe PostPwIdl_3D.o -L${LIBDIR} -lSHARE ${Lflag2}

PwIDL:
	make ${BINDIR}/PostPwIdl.exe
	make ${BINDIR}/PostPwIdl_3D.exe

#
#	cleaning
#

distclean: clean
	rm -f *.o

# keep this line
