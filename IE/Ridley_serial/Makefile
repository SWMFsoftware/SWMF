include Makefile.def

install: src/ModSize.f90
	touch src/Makefile.DEPEND

src/ModSize.f90:
	cp -f src/ModSize_orig.f90 src/ModSize.f90

LIB:
	cd src; make LIB

PIONO:
	cd src; make PostIONO.exe

rundir:
	mkdir -p ${RUNDIR}/IE/ionosphere
	cd ${RUNDIR}/IE; cp ${IEDIR}/Scripts/pION .
	cd ${RUNDIR}/IE; ln -s ${BINDIR}/PostIONO.exe .

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
	echo "IE/Ridley_serial test is incomplete" > notest.diff
