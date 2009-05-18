include ../Makefile.def

#BOP
#!ROUTINE: share/Makefile - install, clean and distclean share
#!DESCRIPTION:
# This Makefile has three targets: install, clean and distclean
# The install target installs the share/Library/src,
# copies the OS and MPIVERSION dependent mpif90.h and mpif.h files into Library/src,
# copies the OS and COMPILER dependent Makefile from build to the parent directory.
# The OS, COMPILER and MPIVERSION variables should be set like in the following example:
#\begin{verbatim}
# make install OS=Linux COMPILER=ifort MPIVERSION=Altix
#\end{verbatim}
# The OS variable should always be set, the COMPILER and MPIVERSION are only needed if they
# are not the defaults.
#
# The clean and distclean targets do not need any variables.
#EOP
#BOC

Library/src/mpif.h:
	cd Library/src; cat precision.h mpif90.h > mpif.h

INSTALL_FILES = \
	Library/src/Makefile.DEPEND \
	Library/src/Makefile.RULES

install: Library/src/mpif.h
	touch ${INSTALL_FILES}
	@(if [ "${OS}" != "Darwin" ]; then \
		rm -f Library/src/ModUtilities.f90; \
	fi);

clean:
	cd Library/src; make clean
	cd Library/test;make clean
	cd Prologs;     make clean

distclean: clean
	cd Library/test;make distclean
	cd Prologs;     make distclean
	rm -f Library/src/mpif*.h *~ */*~ ${INSTALL_FILES}

#EOC
