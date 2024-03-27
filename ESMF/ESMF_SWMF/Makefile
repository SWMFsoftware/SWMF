include ../../Makefile.def

# Nothing should be done parallel in this Makefile
.NOTPARALLEL:

install:

ESMF_SWMF:
	cd ${DIR}; ${MAKE} test_swpc_compile
	cd ${DIR}; ${MAKE} LIB
	cd src; ${MAKE} EXE

EXE:
	cd src; ${MAKE} EXE

RUNDIR = run

${RUNDIR}:
	mkdir -p ${RUNDIR}

rundir: ${RUNDIR}
	mkdir -p ${RUNDIR}/IE/ionosphere
	cp ${IEDIR}/input/cond_???_coeffs.dat run/IE/
	cp input/ESMF_SWMF.input ${RUNDIR}/
	cp input/PARAM.in ${RUNDIR}/
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

test_rundir:
	${MAKE} rundir

test_run:
	cd ${RUNDIR}; rm -f *ESMF_LogFile; ${MPIRUN} ESMF_SWMF.exe | tee runlog

test_check:
	cd ${RUNDIR}; perl -ne 'print if /error/i and not /Multi_On_Error/' \
		PET*.ESMF_LogFile > ../test_esmf.diff
	cd ${RUNDIR}; perl -ne 'print if /ERROR/' runlog >> ../test_esmf.diff
	ls -l test_esmf.diff

clean:
	@(if([ "${ESMF_DIR}" != "" ]); then cd src; make clean; fi)

distclean:
	@(if([ "${ESMF_DIR}" != "" ]); then cd src; make distclean; fi)
