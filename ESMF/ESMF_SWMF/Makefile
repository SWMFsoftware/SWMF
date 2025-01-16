include ../../Makefile.def

DIFFNUM = ${SCRIPTDIR}/DiffNum.pl -BLESS=${BLESS}

# Nothing should be done parallel in this Makefile
.NOTPARALLEL:

help:
	@echo "Targets:"
	@echo "help                - this help message"
	@echo "install             - install ESMF_SWMF"
	@echo "EXE                 - compile ESMF_SWMF.exe"
	@echo "test                - run test on 2 processors"
	@echo "test NP=4           - run test on 4 processors"
	@echo "test_compile        - configure and compile for test"
	@echo "test_rundir         - create run/ directory for test"
	@echo "test_run NP=3       - create run/ directory for test"
	@echo "test_check          - write test results into test_esmf.diff"
	@echo "clean               - clean object files"
	@echo "distclean           - uninstall"

install:

ESMF_SWMF:
	cd ${DIR}; ${MAKE} test_swpc_compile
	cd ${DIR}; ${MAKE} LIB
	cd src; ${MAKE} EXE

EXE:
	cd src; ${MAKE} EXE

RUNDIR = ${DIR}/ESMF/ESMF_SWMF/run

rundir:
	rm -rf ${RUNDIR}
	cd ../../; ${MAKE} rundir RUNDIR=${RUNDIR}
	cd ${RUNDIR}; rm -f SWMF.exe EIE SWMFDIR core
	cp input/ESMF_SWMF.input ${RUNDIR}/
	cp input/fd_swmf.yaml ${RUNDIR}/ 
	cp -f input/PARAM.in ${RUNDIR}/
	cd ${RUNDIR}; ln -sf ${BINDIR}/ESMF_SWMF.exe .

test:
	@echo "test_compile..." > test_esmf.diff
	${MAKE} test_compile
	@echo "test_rundir..." >> test_esmf.diff
	${MAKE} test_rundir
	@echo "test_run..." >> test_esmf.diff
	${MAKE} test_run
	@echo "test_check..." >> test_esmf.diff
	${MAKE} test_check

test_compile:
	${MAKE} ESMF_SWMF
	cd ../..; ${MAKE} PIDL

test_rundir:
	${MAKE} rundir

test_run:
	cd ${RUNDIR}; rm -f *ESMF_LogFile; ${MPIRUN} ./ESMF_SWMF.exe | tee runlog
	cd ${RUNDIR}; ./PostProc.pl -M RESULTS

test_check:
	@(perl -ne 'print if /error/i and not /Multi_On_Error/' \
		${RUNDIR}/PET*.ESMF_LogFile > test_esmf.diff)
	@(perl -ne 'print if /ERROR/' ${RUNDIR}/RESULTS/runlog >> test_esmf.diff)
	-@(${DIFFNUM} -b -r=1e-5 -a=2e-8 \
		${RUNDIR}/RESULTS/GM/log*.log output/GM.log >> test_esmf.diff)
	-@(${DIFFNUM} -b -r=1e-5 -a=2e-8 \
		${RUNDIR}/RESULTS/IE/IE*.log output/IE.log >> test_esmf.diff)
	ls -l test_esmf.diff

clean:
	@(if([ "${ESMFMKFILE}" != "" ]); then cd src; make clean; fi)

distclean:
	@(if([ "${ESMFMKFILE}" != "" ]); then cd src; make distclean; fi)
