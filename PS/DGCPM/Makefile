include ../../Makefile.def

install:
	touch src/Makefile.DEPEND

LIB:
	cd src; make LIB

rundir:
	mkdir -p ${RUNDIR}/PS/Output
	mkdir -p ${RUNDIR}/PS/Input

clean:
	@touch src/Makefile.DEPEND src/Makefile.RULES
	cd src; make clean

distclean: 
	./Config.pl -uninstall

allclean:
	@touch src/Makefile.DEPEND src/Makefile.RULES
	cd src; make distclean
	rm -f *~

test:
	echo "PS/DGCPM test is incomplete" > notest.diff
