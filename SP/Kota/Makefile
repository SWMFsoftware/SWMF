
include Makefile.def 

install: Makefile.def.orig
	touch src/Makefile.DEPEND srcInterface/Makefile.DEPEND

Makefile.def.orig:
	mv Makefile.def Makefile.def.orig
	cp Makefile.def.orig Makefile.def

LIB:
	cd src; make LIB
	cd srcInterface; make LIB

rundir: 
	mkdir -p ${RUNDIR}/SP

clean: install
	cd src; make clean
	cd srcInterface; make clean

distclean: install
	cd src; make distclean
	cd srcInterface; make distclean
	rm -f Makefile.conf Makefile.def *~
	mv Makefile.def.orig Makefile.def
