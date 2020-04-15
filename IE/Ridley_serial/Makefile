include Makefile.def

install: src/ModSize.f90

src/ModSize.f90:
	cp -f src/ModSize_orig.f90 src/ModSize.f90

LIB:
	cd src; make LIB

PIONO:
	cd src; make PostIONO.exe

rundir:
	mkdir -p ${RUNDIR}/IE/ionosphere
	cd ${RUNDIR}/IE; cp ${IEDIR}/input/cond_hal_coeffs.dat .
	cd ${RUNDIR}/IE; cp ${IEDIR}/input/cond_ped_coeffs.dat .
	cd ${RUNDIR}/IE; cp ${IEDIR}/input/cmee_hal_coeffs.dat .
	cd ${RUNDIR}/IE; cp ${IEDIR}/input/cmee_ped_coeffs.dat .
	cd ${RUNDIR}/IE; cp ${IEDIR}/Scripts/pION .
	cd ${RUNDIR}/IE; ln -s ${BINDIR}/PostIONO.exe .

clean:
	cd src; make clean

distclean: 
	./Config.pl -uninstall

allclean:
	cd src; make distclean
	rm -f *~

test:
	echo "IE/Ridley_serial test is incomplete" > notest.diff
