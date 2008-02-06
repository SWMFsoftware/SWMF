# defines OS
include ../Makefile.def

INSTALL_FILES = \
	Control/src/Makefile.DEPEND \
	Control/src/Makefile.RULES \
	Library/src/Makefile.DEPEND \
	Coupler/src/Makefile.DEPEND \
	Interface/src/Makefile.DEPEND \
	Stubs/src/Makefile.DEPEND

install:
	touch ${INSTALL_FILES}
clean:
	cd Control/src;   make clean
	cd Library/src;   make clean
	cd Coupler/src;   make clean
	cd Interface/src; make clean
	cd Stubs/src;     make clean

distclean: clean
	rm -f ${INSTALL_FILES} *~
