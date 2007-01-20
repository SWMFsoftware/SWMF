include Makefile.def

install: src/ModSize.f90
	@(if [ -f src/Makefile.RULES.${OS}${COMPILER} ]; then                \
		cp -f src/Makefile.RULES.${OS}${COMPILER} src/Makefile.RULES;\
	else \
		rm -f src/Makefile.RULES; touch src/Makefile.RULES; \
	fi);
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
	@touch src/Makefile.DEPEND src/Makefile.RULES
	cd src; make distclean
	rm -f Makefile.conf Makefile.def *~

test:
	echo "IE/Ridley_serial test is incomplete" > notest.diff
