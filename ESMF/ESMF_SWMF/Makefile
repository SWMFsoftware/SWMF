include ../../Makefile.def

install:

EXE:
	cd src; ${MAKE} EXE

RUNDIR = run

${RUNDIR}:
	mkdir -p ${RUNDIR}


rundir: ${RUNDIR}
	cp input/ESMF_SWMF.input ${RUNDIR}/
	cp input/PARAM.in ${RUNDIR}/
	cd ${RUNDIR}; ln -sf ${BINDIR}/ESMF_SWMF.exe .

SWMFLIB = ${LIBDIR}/libSWMF.a

${SWMFLIB}:
	cd ${DIR}; ${MAKE} test_swpc_compile
	cd ${DIR}; ${MAKE} LIB

test: ${SWMFLIB} EXE rundir
	cd ${RUNDIR}; rm -f *.ESMF_LogFile; ${MPIRUN} ESMF_SWMF.exe
	cd ${RUNDIR}; perl -ne 'print if /error/i' PET*.ESMF_LogFile

clean:
	@(if([ "${ESMF_DIR}" != "" ]); then cd src; make clean; fi)

distclean:
	@(if([ "${ESMF_DIR}" != "" ]); then cd src; make distclean; fi)
