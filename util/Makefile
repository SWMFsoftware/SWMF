install:
	touch DATAREAD/srcIndices/Makefile.DEPEND
	touch EMPIRICAL/srcIE/Makefile.DEPEND

clean:
	cd NOMPI/src;           make clean
	cd TIMING/src;          make clean
	cd TIMING/srcEmpty;     make clean
	cd TIMING/doc;          make clean
	cd DATAREAD/srcIndices; make clean
	cd EMPIRICAL/srcIE;     make clean
	cd EMPIRICAL/srcGM;     make clean
	cd EMPIRICAL/srcUA;     make clean  

distclean:
	cd NOMPI/src;           make distclean
	cd TIMING/src;          make distclean
	cd TIMING/srcEmpty;     make distclean
	cd TIMING/doc;          make distclean
	cd DATAREAD/srcIndices; make distclean
	cd EMPIRICAL/srcIE;     make distclean
	cd EMPIRICAL/srcGM;     make distclean
	cd EMPIRICAL/srcUA;     make distclean
	rm -f *~

dist: distclean
	/bin/tar -cf util.tar .
