
include Makefile.def 

rundir: 
	mkdir -p ${RUNDIR}/SP

clean:
	cd src; make clean

distclean: clean
	rm -f Makefile.conf *~ src/io_unit_for77.f src/stdout.h
	echo 'STANDALONE = sp' > Makefile.def
