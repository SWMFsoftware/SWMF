clean:
	cd NOMPI/src;       make clean
	cd TIMING/src;      make clean
	cd TIMING/srcEmpty; make clean
	cd TIMING/doc;      make clean

distclean:
	cd NOMPI/src;       make distclean
	cd TIMING/src;      make distclean
	cd TIMING/srcEmpty; make distclean
	cd TIMING/doc;      make distclean
	rm -f *~

dist: distclean
	/bin/tar -cf util.tar .
