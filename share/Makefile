include ../Makefile.def

INSTALL_FILES = \
	Library/src/mpif90.h \
	Library/src/mpif.h   \
	../Makefile.conf

install: ${INSTALL_FILES}
	touch Library/src/Makefile.DEPEND

Library/src/mpif90.h:
	cp include/mpif90_${OS}.h Library/src/mpif90.h

Library/src/mpif.h: Library/src/mpif90.h
	cd Library/src; cat precision.h mpif90.h > mpif.h

../Makefile.conf:
	cp build/Makefile.${OS}${COMPILER} ../Makefile.conf

distclean:
	rm -f ${INSTALL_FILES} Library/src/Makefile.DEPEND *~ */*~
