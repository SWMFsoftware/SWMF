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

test: ESMF_SWMF rundir
	cd ${RUNDIR}; rm -f *.ESMF_LogFile; ${MPIRUN} ESMF_SWMF.exe
	cd ${RUNDIR}; perl -ne 'print if /error/i' PET*.ESMF_LogFile

clean:
	@(if([ "${ESMF_DIR}" != "" ]); then cd src; make clean; fi)

distclean:
	@(if([ "${ESMF_DIR}" != "" ]); then cd src; make distclean; fi)
