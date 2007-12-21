include Makefile.def

install:
	touch src/Makefile.DEPEND

LIB:
	cd src; make LIB

rundir:
	mkdir -p ${RUNDIR}/IE/Output

clean:
	@touch src/Makefile.DEPEND src/Makefile.RULES
	cd src; make clean

distclean:
	@touch src/Makefile.DEPEND src/Makefile.RULES
	cd src; make distclean
	rm -f *~

test:
	echo "IE/RIM test is incomplete" > notest.diff
