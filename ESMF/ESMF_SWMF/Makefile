#
# Default target is ESMF_SWMF so it is listed first
#
default: ESMF_SWMF

include ../../Makefile.def

DIFFNUM = ${SCRIPTDIR}/DiffNum.pl -BLESS=${BLESS}

# Nothing should be done parallel in this Makefile
.NOTPARALLEL:

help:
	@echo "Targets:"
	@echo "ESMF_SWMF (default) - compile libSWMF Geospace and ESMF_SWMF.exe"
	@echo "EXE                 - compile only ESMF_SWMF.exe"
	@echo "help                - this help message"
	@echo "install             - install ESMF_SWMF"
	@echo "test                - run test1"
	@echo "test1               - run fake IPE test on 2 processors"
	@echo "test1 NP=4          - run test on 4 processors"
	@echo "test1_compile       - configure and compile for test"
	@echo "test1_rundir        - create run/ directory for test"
	@echo "test1_run NP=3      - create run/ directory for test"
	@echo "test1_check         - write test results into test_esmf.diff"
	@echo "test2               - run real IPE two-way coupling test"
	@echo "clean               - clean object files"
	@echo "distclean           - uninstall"

install:

ESMF_SWMF:
	cd ${DIR}; ${MAKE} test_swpc_compile
	cd ${DIR}; ${MAKE} LIB
	rm -f ${BINDIR}/ESMF_SWMF.exe
	cd src; ${MAKE} EXE

EXE:
	cd src; ${MAKE} EXE

RUNDIR = ${DIR}/ESMF/ESMF_SWMF/run

# It is assumed the IPE grid and restart data are located in ${IPEDATA}
IPEDATA = ${DIR}/ESMF/ESMF_SWMF/data

rundir:
	rm -rf ${RUNDIR}
	cd ../../; ${MAKE} rundir RUNDIR=${RUNDIR}
	cd ${RUNDIR}; rm -f EIE SWMFDIR core
	cp input/ESMF_SWMF.input ${RUNDIR}/
	cp input/fd_swmf.yaml ${RUNDIR}/ 
	cp -f input/PARAM.in ${RUNDIR}/
	cd ${RUNDIR}; ln -sf ${BINDIR}/ESMF_SWMF.exe .

test:
	${MAKE} test1

test1:
	@echo "test1_compile..." > test_esmf.diff
	${MAKE} test1_compile
	@echo "test1_rundir..." >> test_esmf.diff
	${MAKE} test1_rundir
	@echo "test1_run..." >> test_esmf.diff
	${MAKE} test1_run
	@echo "test1_check..." >> test_esmf.diff
	${MAKE} test1_check

test1_compile:
	${MAKE} ESMF_SWMF EXTERNAL_IPE=no
	cd ../..; ${MAKE} PIDL

test1_rundir:
	${MAKE} rundir

test1_run:
	cd ${RUNDIR}; rm -f *ESMF_LogFile; ${MPIRUN} ./ESMF_SWMF.exe | tee runlog
	cd ${RUNDIR}; ./PostProc.pl -M RESULTS

test1_check:
	@(perl -ne 'print if /error/i and not /Multi_On_Error/' \
		${RUNDIR}/PET*.ESMF_LogFile > test_esmf.diff)
	@(perl -ne 'print if /ERROR/' ${RUNDIR}/RESULTS/runlog >> test_esmf.diff)
	-@(${DIFFNUM} -b -r=1e-5 -a=2e-8 \
		${RUNDIR}/RESULTS/GM/log*.log output/GM.log >> test_esmf.diff)
	-@(${DIFFNUM} -b -r=1e-5 -a=2e-8 \
		${RUNDIR}/RESULTS/IE/IE*.log output/IE.log >> test_esmf.diff)
	ls -l test_esmf.diff

#### Real IPE two-way coupling test##########################
test2:
	@echo "test2_compile..." > test_esmf.diff
	${MAKE} test2_compile
	@echo "test2_rundir..." >> test_esmf.diff
	${MAKE} test2_rundir
	@echo "test2_swmf_run..." >> test_esmf.diff
	${MAKE} test2_swmf_run
	@echo "test2_restart..." >> test_esmf.diff
	${MAKE} test2_restart
	@echo "test2_check..." >> test_esmf.diff
	${MAKE} test2_check

test2_compile:
	${RM} src/ESMF_grid_comp.o
	${MAKE} ESMF_SWMF EXTERNAL_IPE=yes
	cd ../..; ${MAKE} SWMF PIDL

test2_rundir:
	${MAKE} rundir		
	cp input/ESMF_SWMF.input.realIPE.twoway ${RUNDIR}/ESMF_SWMF.input
	cp input/ionprof ${RUNDIR}/
	cp input/tiros_spectra ${RUNDIR}/
	cp input/IPE.inp ${RUNDIR}/
	cp input/wam_input_f107_kp.txt ${RUNDIR}/
	ln -sf ${IPEDATA}/IPE_Grid.nc ${RUNDIR}/
	ln -sf ${IPEDATA}/IPE_State.apex.201303160000.nc ${RUNDIR}/
	ln -sf ${IPEDATA}/global_idea_coeff_W05*.dat ${RUNDIR}/

test2_swmf_run:	
	cp input/PARAM.in.swmf.only ${RUNDIR}/PARAM.in
	cd ${RUNDIR}; ${MPIRUN} ./SWMF.exe > runlog
	cd ${RUNDIR}; ./PostProc.pl RESULTS/run1

test2_restart:	
	cp input/PARAM.in.restart ${RUNDIR}/PARAM.in
	cd ${RUNDIR}; ./Restart.pl -v -i RESULTS/run1/RESTART
	cd ${RUNDIR}; ${MPIRUN} ./ESMF_SWMF.exe > runlog
	cd ${RUNDIR}; ./PostProc.pl RESULTS/run2


test2_check:
	-@(${DIFFNUM} -b -a=1e-6 -r=1e-4 \
                ${RUNDIR}/RESULTS/run2/IE/it_aur_130316_001500_000.idl \
                output/it_aur_130316_001500_000.idl.gz \
		> test2.diff)
		ls -l test2.diff

clean:
	@(if([ "${ESMFMKFILE}" != "" ]); then cd src; make clean; fi)

distclean:
	@(if([ "${ESMFMKFILE}" != "" ]); then cd src; make distclean; fi)
