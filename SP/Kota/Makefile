
include Makefile.def 

rundir: 
	mkdir -p ${RUNDIR}/SP

clean:
	cd src; make clean

distclean: clean
	rm -f Makefile.conf *~
	echo 'STANDALONE = sp' > Makefile.def
