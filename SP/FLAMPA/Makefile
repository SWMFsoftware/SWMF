
include ../../Makefile.def 

install:
	touch src/Makefile.DEPEND srcInterface/Makefile.DEPEND

LIB:    install
	cd src;          make LIB
	cd srcInterface; make LIB

rundir: 
	mkdir -p ${RUNDIR}/SP

clean:  install
	cd src; make clean
	cd srcInterface; make clean

distclean: install
	cd src; make distclean
	cd srcInterface; make distclean
	rm -f Makefile.conf
