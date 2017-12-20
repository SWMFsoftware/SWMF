#  Copyright (C) 2002 Regents of the University of Michigan, 
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

DEFAULT_TARGET = MFLAMPA

MFLAMPA:
	cd src; make MFLAMPA

SHELL =/bin/sh
include Makefile.def 

install: Makefile.def.orig src/ModSize.f90
	touch src/Makefile.DEPEND srcInterface/Makefile.DEPEND

src/ModSize.f90: src/ModSize_orig.f90
	cp -f src/ModSize_orig.f90 src/ModSize.f90

Makefile.def.orig:
	mv Makefile.def Makefile.def.orig
	cp Makefile.def.orig Makefile.def

LIB:    install
	cd src;          make LIB
	cd srcInterface; make LIB

rundir: 
	mkdir -p ${RUNDIR}/SP ${RUNDIR}/SP/restartOUT ${RUNDIR}/SP/restartIN ${RUNDIR}/SP/MHDATA ${RUNDIR}/SP/IO2
	cd ${RUNDIR}; ln -s ${DIR}/bin/MFLAMPA.exe .

clean:  install
	@(if [ -r "Makefile.conf" ]; then \
		cd src; make clean; \
		cd ../srcInterface; make clean; \
	fi)

distclean: 
	./Config.pl -uninstall

allclean: install
	@(if [ -r "Makefile.conf" ]; then \
		cd src; make distclean; \
		cd ../srcInterface; make distclean; \
		rm -f ../Makefile.conf; \
	fi)
	rm -f Makefile.def *~
	mv Makefile.def.orig Makefile.def
