include Makefile.def 

install: Makefile.def.orig

Makefile.def.orig:
	mv Makefile.def Makefile.def.orig
	cp Makefile.def.orig Makefile.def

rundir: 
	mkdir -p ${RUNDIR}/SP

clean:
	cd src; make clean

distclean: clean
	rm -f Makefile.conf *~ src/io_unit_for77.f src/stdout.h
	mv Makefile.def.orig Makefile.def
