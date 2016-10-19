default: IPIC3D

include Makefile.conf
include Makefile.def

# Serial and parallel execution defaults:
SERIAL   =
PARALLEL = mpirun
NPFLAG   = -np
NP       = 2
MPIRUN   = ${PARALLEL} ${NPFLAG} ${NP}

help:
	@echo Makefile targets:
	@echo
	@echo 'make                        - compile IPIC3D.exe'
	@echo 'make LIB                    - compile libPW.a for SWMF'
	@echo 'make PIDL                   - compile PostIDL.exe for post processing'
	@echo 'make run                    - create run directory'
	@echo 'make test                   - test IPIC3D in stand alone mode'
	@echo 'make clean                  - remove object files'
	@echo 'make distclean              - remove all files not part of CVS'
	@echo

bin:
	mkdir bin

libIPIC3D:
	mkdir libIPIC3D

INSTALLFILES =  src/Makefile.DEPEND \
		src/Makefile.local \
		srcInterface/Makefile.DEPEND

install: bin libIPIC3D
	cd src; cp -f Makefile.local.default Makefile.local
	touch ${INSTALLFILES}

IPIC3D:
	./Config.pl -hdf5
	cd src; $(MAKE) IPIC3D

PIDL:
	cd ${SHAREDIR}; $(MAKE) PIDL

serialrun: IPIC3D
	cd ${RUNDIR}; ${SERIAL} ./IPIC3D.exe IPIC.in

LIB: 
	cd src; $(MAKE) LIB
	cd srcInterface; $(MAKE) LIB

TESTDIR = run_test

CHECKDIR  = output

test:
	-@($(MAKE) test_ipic3d)

test_ipic3d:
	@echo "test_compile..." > test_ipic3d.diff
	$(MAKE) test_ipic3d_compile
	@echo "test_rundir..." >> test_ipic3d.diff
	$(MAKE) test_ipic3d_rundir
	@echo "test_run..." >> test_ipic3d.diff
	$(MAKE) test_ipic3d_run
	@echo "test_check..." >> test_ipic3d.diff
	$(MAKE) test_ipic3d_check

test_ipic3d_compile:
	cd src; rm -f Makefile.local; \
	echo "FLAGC_LOCAL=" >Makefile.local
	./Config.pl -pseudrand=true
	-$(MAKE) IPIC3D
	cd src; cp -f Makefile.local.default Makefile.local

test_ipic3d_rundir:
	rm -rf ${TESTDIR}
	${MAKE} rundir RUNDIR=${TESTDIR} STANDALONE="YES" PWDIR=`pwd`

test_ipic3d_run:
	cd ${TESTDIR}; ${PARALLEL} ${NPFLAG} 4 ./IPIC3D.exe IPIC.in > runlog

test_ipic3d_check:
	-@(${SCRIPTDIR}/DiffNum.pl -b -a=2.e-13 -r=1e-7\
		${TESTDIR}/PC/plots/ConservedQuantities.txt \
		${CHECKDIR}/ConservedQuantities.txt\
		> test_ipic3d.diff)
	@ls -l test_ipic3d.diff

test_h5hut:
	@echo "test_compile..." > test_h5hut.diff
	$(MAKE) test_h5hut_compile
	@echo "test_rundir..." >> test_h5hut.diff
	$(MAKE) test_h5hut_rundir
	@echo "test_run..." >> test_h5hut.diff
	$(MAKE) test_h5hut_run
	@echo "test_check..." >> test_h5hut.diff
	$(MAKE) test_h5hut_check

test_h5hut_compile:
	./Config.pl -npx=2 -npy=2 -npz=1 -periodicx=true -periodicY=false -periodicZ=true -ng=3 -pseudrand=true -h5hut=true
	$(MAKE) IPIC3D

test_h5hut_rundir:
	rm -rf ${TESTDIR}
	${MAKE} rundir RUNDIR=${TESTDIR} STANDALONE="YES" PWDIR=`pwd`

test_h5hut_run:
	cd ${TESTDIR}; ${PARALLEL} ${NPFLAG} 4 ./IPIC3D.exe IPIC.in > runlog

test_h5hut_check:
	-@(${SCRIPTDIR}/DiffNum.pl -b -a=2.e-13 -r=1e-7\
		${TESTDIR}/PC/plots/ConservedQuantities.txt \
		${CHECKDIR}/ConservedQuantities.txt-GEM-PSEUDRAND\
		> test_h5hut.diff)
	@ls -l test_h5hut.diff


test_ipic3d_vtk:
	$(MAKE) clean
	rm -f ${BINDIR}/IPIC3D.exe
	@echo "test_compile..." > test_ipic3d_vtk.diff
	$(MAKE) test_ipic3d_vtk_compile
	@echo "test_rundir..." >> test_ipic3d_vtk.diff
	$(MAKE) test_ipic3d_vtk_rundir
	@echo "test_run..." >> test_ipic3d_vtk.diff
	$(MAKE) test_ipic3d_vtk_run
	@echo "test_check..." >> test_ipic3d_vtk.diff
	$(MAKE) test_ipic3d_vtk_check

test_ipic3d_vtk_compile:
	./Config.pl -npx=2 -npy=2 -npz=2 -periodicx=true -periodicy=false -periodicz=true -ng=1 -pseudrand=true -hdf5
	cd src; $(MAKE) IPIC3D

test_ipic3d_vtk_rundir:
	rm -rf ${TESTDIR}
	${MAKE} rundir RUNDIR=${TESTDIR} STANDALONE="YES" PWDIR=`pwd`

test_ipic3d_vtk_run:
	cd ${TESTDIR}; ${PARALLEL} ${NPFLAG} 8 ./IPIC3D.exe ../input/GEM3D.inp | tee runlog

test_ipic3d_vtk_check:
	-@(${SCRIPTDIR}/DiffNum.pl -b -a=2.e-13 -r=1e-7\
		${TESTDIR}/PC/plots/ConservedQuantities.txt \
		${CHECKDIR}/ConservedQuantities.txt-GEM3DVTK-PSEUDRAND\
		> test_ipic3d_vtk.diff)
	@ls -l test_ipic3d_vtk.diff



rundir:
	mkdir -p ${RUNDIR}/PC
	cd ${RUNDIR}/PC; \
		mkdir restartIN restartOUT plots;\
		ln -s ${BINDIR}/PostIDL.exe .; \
		cp    ${SCRIPTDIR}/pIDL .;
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		cd ${RUNDIR}; \
			ln -s ${BINDIR}/IPIC3D.exe .; \
			cp ${PWDIR}/input/IPIC.in .; \
			touch core ; chmod 444 core ; \
			ln -s PC/* .; \
	fi)

clean:
	touch ${INSTALLFILES}
	cd src; $(MAKE) clean
	cd srcInterface; $(MAKE) clean

distclean:
	-@(./Config.pl -uninstall)

allclean:
	-@(cd src; $(MAKE) distclean)
	-@(cd srcInterface; $(MAKE) distclean)
	-@(rm -rf Makefile.local *~ ./bin libIPIC3D ${TESTDIR})
	-@(rm test*diff )
