include ../Makefile.def

install: Library/src/mpif90.h Library/src/mpif.h ../Makefile.conf
	touch Library/src/Makefile.DEPEND

Library/src/mpif90.h:
	cp include/mpif90_${OS}.h Library/src/mpif90.h

Library/src/mpif.h: Library/src/mpif90.h
	cd Library/src; cat precision.h mpif90.h > mpif.h

../Makefile.conf:
	cp build/Makefile.${OS}${COMPILER} ../Makefile.conf

clean:
	cd Library/src; make clean
	cd Prologs;     make clean

distclean: clean
	cd Prologs;     make distclean
	rm -f Library/src/mpif*.h Library/src/Makefile.DEPEND *~ */*~
