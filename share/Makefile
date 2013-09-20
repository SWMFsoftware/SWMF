include ../Makefile.def

INSTALL_FILES = \
	Library/src/Makefile.DEPEND \
	Library/src/Makefile.RULES

install:
	touch ${INSTALL_FILES}
	@(if [ "${OS}" != "Darwin" ]; then \
		rm -f Library/src/ModUtilities.f90; \
	fi);

clean:
	cd Library/src; make clean
	cd Library/test;make clean
	cd Prologs;     make clean

distclean: clean
	cd Library/src; make distclean
	cd Library/test;make distclean
	cd Prologs;     make distclean
	rm -f Library/src/mpif*.h *~ */*~ ${INSTALL_FILES}

